#include <OpenGL/gl3.h>

#include "renderer.h"
#include <iostream>
#include <stdexcept>
#include <string>

// ─────────────────────────────────────────────────────────────────────────────
// Shader A – orbital cloud + coordinate axes
// value encoding:
//   ±1.0  → orbital lobe (blue / orange)
//    2.0  → X axis (red)
//    3.0  → Y axis (green)
//    4.0  → Z axis (blue-grey)
// ─────────────────────────────────────────────────────────────────────────────
static const char* ORB_VERT = R"GLSL(
#version 330 core
layout(location = 0) in vec3  aPos;
layout(location = 1) in float aVal;
out float vVal;
uniform mat4 uMVP;
uniform float uAlphaScale;   // dims the cloud when flow is shown
void main() {
    gl_Position  = uMVP * vec4(aPos, 1.0);
    gl_PointSize = (abs(aVal) < 1.5) ? 2.5 : 3.5;
    vVal = aVal;
}
)GLSL";

static const char* ORB_FRAG = R"GLSL(
#version 330 core
in  float vVal;
out vec4  FragColor;
uniform float uAlphaScale;
uniform bool uComplex;
void main() {
    // ── Axis lines (value > 1.5, drawn as GL_LINES) ─────────────────────────
    if (vVal > 1.5) {
        int axis = int(vVal + 0.5) - 2;
        vec3 col = (axis == 0) ? vec3(1.00, 0.28, 0.28)
                 : (axis == 1) ? vec3(0.28, 1.00, 0.28)
                 :                vec3(0.40, 0.60, 1.00);
        FragColor = vec4(col, 0.42);
        return;
    }
    // ── Orbital cloud (GL_POINTS, round disc) ────────────────────────────────
    vec2  pc = gl_PointCoord - vec2(0.5);
    float d  = dot(pc, pc);
    if (d > 0.25) discard;
    float alpha = uAlphaScale * 0.88 * (1.0 - smoothstep(0.16, 0.25, d));
    vec3 col = (vVal >= 0.0)
        ? vec3(0.33, 0.64, 1.00)   // positive lobe – cornflower blue
        : vec3(1.00, 0.42, 0.26);  // negative lobe – burnt orange
    if (uComplex) col = vec3(0.51, 0.72, 0.82);
    FragColor = vec4(col, alpha);
}
)GLSL";

// ─────────────────────────────────────────────────────────────────────────────
// Shader B – flow particles (plasma colormap, value = speed_norm ∈ [0,1])
// ─────────────────────────────────────────────────────────────────────────────
static const char* FLOW_VERT = R"GLSL(
#version 330 core
layout(location = 0) in vec3  aPos;
layout(location = 1) in float aSpeed;
out float vSpeed;
uniform mat4 uMVP;
void main() {
    gl_Position  = uMVP * vec4(aPos, 1.0);
    gl_PointSize = 3.2;
    vSpeed = aSpeed;
}
)GLSL";

static const char* FLOW_FRAG = R"GLSL(
#version 330 core
in  float vSpeed;
out vec4  FragColor;

// Matplotlib "plasma" colormap (5 control points)
vec3 plasma(float t) {
    t = clamp(t, 0.0, 1.0);
    vec3 c0 = vec3(0.050, 0.030, 0.528);   // dark purple
    vec3 c1 = vec3(0.494, 0.012, 0.658);   // purple
    vec3 c2 = vec3(0.800, 0.278, 0.470);   // magenta-pink
    vec3 c3 = vec3(0.973, 0.582, 0.250);   // orange
    vec3 c4 = vec3(0.940, 0.975, 0.131);   // bright yellow-green
    if (t < 0.25) return mix(c0, c1, t * 4.0);
    if (t < 0.50) return mix(c1, c2, (t - 0.25) * 4.0);
    if (t < 0.75) return mix(c2, c3, (t - 0.50) * 4.0);
    return             mix(c3, c4, (t - 0.75) * 4.0);
}

void main() {
    vec2  pc = gl_PointCoord - vec2(0.5);
    float d  = dot(pc, pc);
    if (d > 0.25) discard;
    float alpha = 0.95 * (1.0 - smoothstep(0.10, 0.25, d));
    FragColor = vec4(plasma(vSpeed), alpha);
}
)GLSL";

// ─────────────────────────────────────────────────────────────────────────────
// Internal: compile + link one shader program
// ─────────────────────────────────────────────────────────────────────────────
static unsigned int make_program(const char* vert_src, const char* frag_src) {
    auto compile = [](const char* src, GLenum type) -> unsigned int {
        unsigned int id = glCreateShader(type);
        glShaderSource(id, 1, &src, nullptr);
        glCompileShader(id);
        int ok = 0;
        glGetShaderiv(id, GL_COMPILE_STATUS, &ok);
        if (!ok) {
            char log[1024];
            glGetShaderInfoLog(id, sizeof(log), nullptr, log);
            glDeleteShader(id);
            throw std::runtime_error(std::string("Shader compilation failed: ") + log);
        }
        return id;
    };

    unsigned int v = 0, f = 0, p = 0;
    try {
        v = compile(vert_src, GL_VERTEX_SHADER);
        f = compile(frag_src, GL_FRAGMENT_SHADER);
        p = glCreateProgram();
        glAttachShader(p, v);
        glAttachShader(p, f);
        glLinkProgram(p);
        int ok = 0;
        glGetProgramiv(p, GL_LINK_STATUS, &ok);
        if (!ok) {
            char log[1024];
            glGetProgramInfoLog(p, sizeof(log), nullptr, log);
            throw std::runtime_error(std::string("Shader linking failed: ") + log);
        }
    } catch (...) {
        if (p) glDeleteProgram(p);
        if (v) glDeleteShader(v);
        if (f) glDeleteShader(f);
        throw;
    }
    glDeleteShader(v);
    glDeleteShader(f);
    return p;
}

// Internal: allocate one VAO+VBO with the shared [x y z value] layout
static void setup_vao(unsigned int vao, unsigned int vbo) {
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    constexpr int stride = 4 * sizeof(float);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride,
                          reinterpret_cast<void*>(0));
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, stride,
                          reinterpret_cast<void*>(3 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
}

// ─────────────────────────────────────────────────────────────────────────────
// Renderer lifecycle
// ─────────────────────────────────────────────────────────────────────────────
Renderer::Renderer() {
    m_shader = build_orbital_shader();
    try {
        m_flow_shader = build_flow_shader();
    } catch (...) {
        glDeleteProgram(m_shader);
        throw;
    }

    glGenVertexArrays(1, &m_cloud_vao); glGenBuffers(1, &m_cloud_vbo);
    glGenVertexArrays(1, &m_axes_vao);  glGenBuffers(1, &m_axes_vbo);
    glGenVertexArrays(1, &m_flow_vao);  glGenBuffers(1, &m_flow_vbo);

    setup_vao(m_cloud_vao, m_cloud_vbo);
    setup_vao(m_axes_vao,  m_axes_vbo);
    setup_vao(m_flow_vao,  m_flow_vbo);

    upload_axes(1.0f);
}

Renderer::~Renderer() {
    glDeleteBuffers(1, &m_cloud_vbo);   glDeleteVertexArrays(1, &m_cloud_vao);
    glDeleteBuffers(1, &m_axes_vbo);    glDeleteVertexArrays(1, &m_axes_vao);
    glDeleteBuffers(1, &m_flow_vbo);    glDeleteVertexArrays(1, &m_flow_vao);
    glDeleteProgram(m_shader);
    glDeleteProgram(m_flow_shader);
}

// ─────────────────────────────────────────────────────────────────────────────
// Data upload
// ─────────────────────────────────────────────────────────────────────────────
void Renderer::upload_cloud(const std::vector<PointData>& pts) {
    m_cloud_count = static_cast<int>(pts.size());
    if (m_cloud_count == 0) return;
    glBindBuffer(GL_ARRAY_BUFFER, m_cloud_vbo);
    glBufferData(GL_ARRAY_BUFFER,
                 m_cloud_count * (GLsizeiptr)sizeof(PointData),
                 pts.data(), GL_DYNAMIC_DRAW);
}

void Renderer::upload_axes(float half_size) {
    const float L = half_size;
    const float data[] = {
        -L, 0,  0,  2.0f,   L,  0,  0, 2.0f,
         0, -L, 0,  3.0f,   0,  L,  0, 3.0f,
         0,  0, -L, 4.0f,   0,  0,  L, 4.0f,
    };
    glBindBuffer(GL_ARRAY_BUFFER, m_axes_vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(data), data, GL_STATIC_DRAW);
}

void Renderer::upload_flow(const std::vector<FlowParticle>& pts) {
    m_flow_count = static_cast<int>(pts.size());
    if (m_flow_count == 0) return;
    glBindBuffer(GL_ARRAY_BUFFER, m_flow_vbo);
    // FlowParticle has the same binary layout as PointData: 4 × float
    glBufferData(GL_ARRAY_BUFFER,
                 m_flow_count * (GLsizeiptr)sizeof(FlowParticle),
                 pts.data(), GL_STREAM_DRAW);  // STREAM: updated every frame
}

// ─────────────────────────────────────────────────────────────────────────────
// Draw
// ─────────────────────────────────────────────────────────────────────────────
void Renderer::draw(const float* mvp, bool show_flow, bool complex_basis, bool show_axes) const {
    // When flow is shown, dim the orbital cloud so the current stands out
    float alpha_scale = show_flow ? 0.35f : 1.0f;

    // ── Orbital shader: cloud + axes ────────────────────────────────────────
    glUseProgram(m_shader);
    glUniformMatrix4fv(glGetUniformLocation(m_shader, "uMVP"),
                       1, GL_FALSE, mvp);
    glUniform1f(glGetUniformLocation(m_shader, "uAlphaScale"), alpha_scale);
    glUniform1i(glGetUniformLocation(m_shader, "uComplex"), complex_basis ? 1 : 0);

    if (show_axes) {
        glBindVertexArray(m_axes_vao);
        glDrawArrays(GL_LINES, 0, 6);
    }

    if (m_cloud_count > 0) {
        glBindVertexArray(m_cloud_vao);
        glDrawArrays(GL_POINTS, 0, m_cloud_count);
    }

    // ── Flow shader: probability-current particles ───────────────────────────
    if (show_flow && m_flow_count > 0) {
        glUseProgram(m_flow_shader);
        glUniformMatrix4fv(glGetUniformLocation(m_flow_shader, "uMVP"),
                           1, GL_FALSE, mvp);
        glBindVertexArray(m_flow_vao);
        glDrawArrays(GL_POINTS, 0, m_flow_count);
    }

    glBindVertexArray(0);
}

// ─────────────────────────────────────────────────────────────────────────────
// Shader builders
// ─────────────────────────────────────────────────────────────────────────────
unsigned int Renderer::build_orbital_shader() {
    return make_program(ORB_VERT, ORB_FRAG);
}

unsigned int Renderer::build_flow_shader() {
    return make_program(FLOW_VERT, FLOW_FRAG);
}
