#pragma once
#include <cmath>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

// ─────────────────────────────────────────────────────────────────────────────
// Spherical orbit camera.
//
// Always looks at the world origin.  The user rotates with mouse drag and
// zooms with the scroll wheel.
//
// Coordinate convention (standard OpenGL right-hand):
//   +X → right,  +Y → up,  +Z → out of screen
//
// Camera position:
//   x = radius · cos(pitch) · sin(yaw)
//   y = radius · sin(pitch)
//   z = radius · cos(pitch) · cos(yaw)
// ─────────────────────────────────────────────────────────────────────────────
class Camera {
public:
    float yaw    = 0.55f;   // horizontal angle (rad)
    float pitch  = 0.30f;   // vertical angle   (rad)
    float radius = 25.0f;   // distance from origin

    static constexpr float PITCH_LIMIT = 1.55f;   // ≈ 88.8°

    // ── Input handlers ───────────────────────────────────────────────────────
    void rotate(float dYaw, float dPitch) {
        yaw   += dYaw;
        pitch += dPitch;
        if (pitch >  PITCH_LIMIT) pitch =  PITCH_LIMIT;
        if (pitch < -PITCH_LIMIT) pitch = -PITCH_LIMIT;
    }

    // delta > 0 → zoom in, delta < 0 → zoom out
    void zoom(float delta) {
        // Multiplicative zoom feels more natural than additive
        radius *= std::pow(0.90f, delta);
        if (radius <    1.0f) radius =    1.0f;
        if (radius > 1500.0f) radius = 1500.0f;
    }

    // Set radius so the orbital (parameterised by half-size L) fills the view
    void fit(float L, float aspect = 1.0f) {
        // Fit the sampling sphere inside the narrower field of view, including
        // portrait windows. A fixed multiple clipped the sphere at 45 degrees.
        const float narrow_aspect = std::fmax(.05f, std::fmin(aspect, 1.0f));
        const float half_angle = std::atan(std::tan(glm::radians(22.5f)) * narrow_aspect);
        radius = L / std::sin(half_angle) * 1.05f;
    }

    // ── Matrices ─────────────────────────────────────────────────────────────
    glm::vec3 position() const {
        return glm::vec3(
            radius * std::cos(pitch) * std::sin(yaw),
            radius * std::sin(pitch),
            radius * std::cos(pitch) * std::cos(yaw)
        );
    }

    glm::mat4 view() const {
        return glm::lookAt(
            position(),
            glm::vec3(0.0f),          // always look at origin
            glm::vec3(0.0f, 1.0f, 0.0f)
        );
    }

    glm::mat4 projection(float aspect) const {
        return glm::perspective(
            glm::radians(45.0f),
            aspect,
            0.1f,
            4000.0f
        );
    }
};
