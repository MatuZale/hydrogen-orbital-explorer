#pragma once
#include <atomic>
#include <vector>

struct FlowParticle { float x, y, z, speed_norm; };

// For a complex hydrogen eigenstate, j/|psi|^2 = m/rho in the azimuthal
// direction (atomic units), where rho is cylindrical radius. This current
// does not describe the real tesseral states, whose stationary current is zero.
// Seed the matching complex density. count=0 means 25000; m=0 returns empty.
// Cancellation returns an empty vector, never a partial result.
std::vector<FlowParticle> generate_flow_particles(
    int n, int l, int m_quantum, int count = 0, unsigned seed = 99999,
    const std::atomic<bool>* cancel = nullptr);

// Exact azimuthal rotation: delta phi = m/rho^2 * visual_scale * dt.
// No per-frame angular clamp. Each point preserves radius and z, within
// float storage roundoff. dt and speed_mult must be finite and nonnegative.
// The deliberately artificial visual timescale gives one turn in 5 seconds
// at rho=n^2, speed_mult=1, for either nonzero sign of m. Rapid turns can be
// temporally aliased on screen; these tracers are not classical electron paths.
void update_flow_particles(std::vector<FlowParticle>& pts, int m_quantum,
                           int n, float dt, float speed_mult = 1.0f);
