#pragma once
#include <vector>
#include "sampling.h"
#include "flow.h"

// ─────────────────────────────────────────────────────────────────────────────
// Renderer – OpenGL point-cloud + axis-line + flow-particle visualiser.
//
// All vertex buffers share the same layout:
//   [ x  y  z  value ]   (4 × GL_FLOAT, stride 16 bytes)
//
// Orbital cloud : value = ±1 (lobe sign)  →  blue / orange-red
// Flow particles: value = speed_norm ∈ [0,1]  →  plasma colormap
// Axis lines    : value = 2/3/4 (axis id)  →  red / green / blue
// ─────────────────────────────────────────────────────────────────────────────
class Renderer {
public:
    Renderer();
    ~Renderer();
    Renderer(const Renderer&) = delete;
    Renderer& operator=(const Renderer&) = delete;

    // ── Orbital cloud ────────────────────────────────────────────────────────
    void upload_cloud(const std::vector<PointData>& pts);

    // ── Axis lines ───────────────────────────────────────────────────────────
    void upload_axes(float half_size);

    // ── Flow particles (re-uploaded every frame because positions change) ────
    void upload_flow(const std::vector<FlowParticle>& pts);

    // ── Draw everything ──────────────────────────────────────────────────────
    // show_flow   : also draw the flow-particle layer
    // Complex states use one neutral cloud colour: they have no real lobe sign.
    void draw(const float* mvp, bool show_flow = false,
              bool complex_basis = false, bool show_axes = true) const;

    int point_count()  const { return m_cloud_count; }
    int flow_count()   const { return m_flow_count;  }

private:
    // orbital point cloud
    unsigned int m_cloud_vao   = 0;
    unsigned int m_cloud_vbo   = 0;
    int          m_cloud_count = 0;

    // coordinate axes
    unsigned int m_axes_vao = 0;
    unsigned int m_axes_vbo = 0;

    // flow particles
    unsigned int m_flow_vao   = 0;
    unsigned int m_flow_vbo   = 0;
    int          m_flow_count = 0;

    // orbital / axes shader
    unsigned int m_shader = 0;
    // flow shader (plasma colormap)
    unsigned int m_flow_shader = 0;

    unsigned int build_orbital_shader();
    unsigned int build_flow_shader();
};
