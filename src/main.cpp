#include <OpenGL/gl3.h>
#include <GLFW/glfw3.h>
#include <glm/gtc/type_ptr.hpp>

#include <algorithm>
#include <atomic>
#include <chrono>
#include <exception>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "camera.h"
#include "flow.h"
#include "quantum.h"
#include "renderer.h"
#include "sampling.h"

struct Selection {
    int orbital = 0;
    unsigned seed = 42;
    OrbitalBasis basis = OrbitalBasis::Real;
    unsigned long long revision = 0;
};

static Camera g_cam;
static Selection g_requested;
static bool g_needs_generation = true, g_title_dirty = true, g_show_axes = true;
static float g_speed = 1.0f;
static bool g_mouse_down = false;
static double g_last_mx = 0, g_last_my = 0;

static void print_help() {
    std::cout << "Hydrogen Orbital Explorer (macOS)\n\n"
              << "Left/Right, Tab/Shift+Tab  Cycle states\n"
              << "1/2/3/4                   Jump to a principal shell\n"
              << "F                         Switch real orbitals / complex m states\n"
              << "+/-                       Change visual flow speed (complex mode)\n"
              << "R                         Resample with a new seed\n"
              << "A                         Show/hide coordinate axes\n"
              << "Left drag / scroll        Orbit camera / zoom\n"
              << "Esc                       Quit\n\n"
              << "Real mode: blue/orange are the signs of the full real wavefunction.\n"
              << "The real basis index labels sine/cosine combinations; it is not an m eigenvalue.\n"
              << "Complex mode: the cloud and particles share |R Y_l^m|^2; m is the magnetic\n"
              << "quantum number. The neutral cloud does not encode phase. Coloured particles\n"
              << "illustrate probability current, not electron paths. m=0 has no current.\n"
              << "Time and the plasma speed colours are scaled for viewing, not SI units.\n";
}

static std::string describe(const Selection& selection) {
    const auto& o = ORBITALS[selection.orbital];
    std::ostringstream text;
    if (selection.basis == OrbitalBasis::Real)
        text << "Real " << o.name << " [n=" << o.n << ", l=" << o.l << ", basis=" << o.m << "]";
    else
        text << "Complex [n=" << o.n << ", l=" << o.l << ", m=" << o.m << "]";
    return text.str();
}

static void request_generation() {
    ++g_requested.revision;
    g_needs_generation = true;
    g_title_dirty = true;
}

static int first_idx_for_n(int n) {
    for (int i = 0; i < NUM_ORBITALS; ++i)
        if (ORBITALS[i].n == n) return i;
    return 0;
}

static void cb_key(GLFWwindow* window, int key, int, int action, int mods) {
    if (action != GLFW_PRESS && action != GLFW_REPEAT) return;
    const bool shift = (mods & GLFW_MOD_SHIFT) != 0;
    switch (key) {
    case GLFW_KEY_ESCAPE: glfwSetWindowShouldClose(window, GLFW_TRUE); break;
    case GLFW_KEY_TAB:
    case GLFW_KEY_RIGHT:
        g_requested.orbital = (g_requested.orbital + (shift ? NUM_ORBITALS - 1 : 1)) % NUM_ORBITALS;
        request_generation(); break;
    case GLFW_KEY_LEFT:
        g_requested.orbital = (g_requested.orbital + NUM_ORBITALS - 1) % NUM_ORBITALS;
        request_generation(); break;
    case GLFW_KEY_1: case GLFW_KEY_2: case GLFW_KEY_3: case GLFW_KEY_4:
        g_requested.orbital = first_idx_for_n(key - GLFW_KEY_1 + 1);
        request_generation(); break;
    case GLFW_KEY_R: ++g_requested.seed; request_generation(); break;
    case GLFW_KEY_F:
        if (action == GLFW_PRESS) {
            g_requested.basis = g_requested.basis == OrbitalBasis::Real
                ? OrbitalBasis::Complex : OrbitalBasis::Real;
            request_generation();
        }
        break;
    case GLFW_KEY_A:
        if (action == GLFW_PRESS) { g_show_axes = !g_show_axes; g_title_dirty = true; }
        break;
    case GLFW_KEY_EQUAL: case GLFW_KEY_KP_ADD:
        g_speed = std::min(g_speed * 1.5f, 20.0f); g_title_dirty = true; break;
    case GLFW_KEY_MINUS: case GLFW_KEY_KP_SUBTRACT:
        g_speed = std::max(g_speed / 1.5f, 0.05f); g_title_dirty = true; break;
    default: break;
    }
}

static void cb_scroll(GLFWwindow*, double, double dy) { g_cam.zoom(static_cast<float>(dy)); }
static void cb_mouse_btn(GLFWwindow* window, int button, int action, int) {
    if (button != GLFW_MOUSE_BUTTON_LEFT) return;
    g_mouse_down = action == GLFW_PRESS;
    if (g_mouse_down) glfwGetCursorPos(window, &g_last_mx, &g_last_my);
}
static void cb_cursor(GLFWwindow*, double x, double y) {
    if (!g_mouse_down) return;
    g_cam.rotate(static_cast<float>(x - g_last_mx) * .005f,
                 static_cast<float>(y - g_last_my) * .005f);
    g_last_mx = x; g_last_my = y;
}

struct GeneratedData {
    Selection selection;
    std::vector<PointData> cloud;
    std::vector<FlowParticle> flow;
    std::string error;
};

// One cancellable worker owns its data until joined. No detached work can
// outlive the application, and only the main thread touches the GL context.
class Generator {
public:
    ~Generator() { cancel(); if (worker_.joinable()) worker_.join(); }
    void cancel() { cancelled_.store(true); }
    bool running() const { return worker_.joinable(); }
    void start(Selection selection) {
        cancelled_.store(false);
        done_.store(false);
        result_ = GeneratedData{};
        result_.selection = selection;
        worker_ = std::thread([this, selection] {
            const auto started = std::chrono::steady_clock::now();
            try {
                const auto& o = ORBITALS[selection.orbital];
                result_.cloud = generate_orbital(o.n, o.l, o.m, 0, selection.seed,
                                                 selection.basis, &cancelled_);
                if (!cancelled_.load() && selection.basis == OrbitalBasis::Complex)
                    result_.flow = generate_flow_particles(o.n, o.l, o.m, 0,
                                                           selection.seed + 1, &cancelled_);
                if (!cancelled_.load() && result_.cloud.empty())
                    result_.error = "Sampling returned an empty cloud";
            } catch (const std::exception& error) {
                result_.error = error.what();
            } catch (...) {
                result_.error = "Unknown sampling error";
            }
            if (!cancelled_.load())
                std::cout << "[generation] " << describe(selection) << ": "
                          << std::chrono::duration_cast<std::chrono::milliseconds>(
                                 std::chrono::steady_clock::now() - started).count() << " ms\n";
            done_.store(true, std::memory_order_release);
        });
    }
    bool take(GeneratedData& result) {
        if (!worker_.joinable() || !done_.load(std::memory_order_acquire)) return false;
        worker_.join();
        result = std::move(result_);
        return true;
    }
private:
    std::thread worker_;
    std::atomic<bool> cancelled_{false}, done_{false};
    GeneratedData result_;
};

static void run(GLFWwindow* window) {
    Renderer renderer;
    Generator generator; // joins before renderer, which dies before the GLFW context
    Selection active;
    bool has_active = false;
    std::vector<FlowParticle> flow;
    std::string generation_error;
    double previous_time = glfwGetTime();

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();
        if (glfwWindowShouldClose(window)) break;
        const double now = glfwGetTime();
        const float dt = static_cast<float>(std::min(now - previous_time, .1));
        previous_time = now;

        if (g_needs_generation) generator.cancel();
        GeneratedData completed;
        if (generator.take(completed)) {
            // Rapid navigation may finish an older request; never relabel or
            // upload its cloud as the newly selected state.
            if (completed.selection.revision == g_requested.revision) {
                generation_error = completed.error;
                if (generation_error.empty()) {
                    const auto& o = ORBITALS[completed.selection.orbital];
                    renderer.upload_cloud(completed.cloud);
                    renderer.upload_axes(static_cast<float>(orbital_radius(o.n)));
                    flow = std::move(completed.flow);
                    renderer.upload_flow(flow);
                    active = completed.selection;
                    has_active = true;
                    int width, height;
                    glfwGetFramebufferSize(window, &width, &height);
                    g_cam.fit(static_cast<float>(orbital_radius(o.n)),
                              height > 0 ? static_cast<float>(width) / height : 1.0f);
                } else std::cerr << "[generation] " << generation_error << '\n';
                g_title_dirty = true;
            }
        }
        if (g_needs_generation && !generator.running()) {
            g_needs_generation = false;
            generation_error.clear();
            generator.start(g_requested);
            g_title_dirty = true;
        }

        const bool complex_basis = has_active && active.basis == OrbitalBasis::Complex;
        const bool show_flow = complex_basis && !flow.empty();
        if (show_flow) {
            const auto& o = ORBITALS[active.orbital];
            update_flow_particles(flow, o.m, o.n, dt, g_speed);
            renderer.upload_flow(flow);
        }
        if (g_title_dirty) {
            std::ostringstream title;
            title << "Hydrogen | ";
            if (has_active) {
                title << describe(active) << " | " << renderer.point_count() << " points";
                if (complex_basis) {
                    if (flow.empty()) title << " | m=0: no current";
                    else title << " | current " << flow.size() << " particles, " << g_speed << "x";
                }
            }
            if (generator.running() || g_needs_generation)
                title << (has_active ? " | generating: " : "Generating ") << describe(g_requested);
            else if (!generation_error.empty()) title << " | generation failed (R: retry)";
            title << " | F: basis  A: axes " << (g_show_axes ? "on" : "off") << "  arrows: state";
            glfwSetWindowTitle(window, title.str().c_str());
            g_title_dirty = false;
        }

        int width, height;
        glfwGetFramebufferSize(window, &width, &height);
        if (width == 0 || height == 0) {
            glfwWaitEventsTimeout(.1);
            previous_time = glfwGetTime();
            continue;
        }
        glViewport(0, 0, width, height);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        const auto mvp = g_cam.projection(static_cast<float>(width) / height) * g_cam.view();
        if (has_active) renderer.draw(glm::value_ptr(mvp), show_flow, complex_basis, g_show_axes);
        glfwSwapBuffers(window);
    }
}

int main(int argc, char** argv) {
    if (argc > 1) {
        if (argc == 2 && (std::string(argv[1]) == "--help" || std::string(argv[1]) == "-h")) {
            print_help(); return 0;
        }
        std::cerr << "Unknown argument. Use --help for controls.\n";
        return 2;
    }
    glfwSetErrorCallback([](int code, const char* message) {
        std::cerr << "[GLFW " << code << "] " << message << '\n';
    });
    if (!glfwInit()) return 1;
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GLFW_TRUE);
    glfwWindowHint(GLFW_SAMPLES, 4);
    GLFWwindow* window = glfwCreateWindow(1280, 800, "Hydrogen Orbital Explorer", nullptr, nullptr);
    if (!window) { glfwTerminate(); return 1; }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);
    glfwSetKeyCallback(window, cb_key);
    glfwSetScrollCallback(window, cb_scroll);
    glfwSetMouseButtonCallback(window, cb_mouse_btn);
    glfwSetCursorPosCallback(window, cb_cursor);
    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_DEPTH_TEST);
    glClearColor(.04f, .04f, .09f, 1.0f);
    print_help();
    int status = 0;
    try { run(window); }
    catch (const std::exception& error) { std::cerr << "[fatal] " << error.what() << '\n'; status = 1; }
    glfwDestroyWindow(window);
    glfwTerminate();
    return status;
}
