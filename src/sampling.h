#pragma once
#include "quantum.h"
#include <atomic>
#include <vector>

// OpenGL layout: four floats. Real mode colors the full R*Y sign;
// complex mode uses neutral +1 (its continuous complex phase is not shown).
struct PointData { float x, y, z, sign; };

// Separable sampling: numerical inversion of the analytic radial CDF and
// angular rejection with the rigorous addition-theorem bound (2l+1)/(4pi).
// Samples the normalized density conditional on r <= orbital_radius(n).
// No estimated envelope or capped rejection ratio is used.
// target_points=0 selects 70000 (n<=2), 55000 (n=3), or 40000 (n=4).
// Returns exactly that count, or an empty vector if cancellation is requested.
std::vector<PointData> generate_orbital(
    int n, int l, int m, int target_points = 0, unsigned seed = 42,
    OrbitalBasis basis = OrbitalBasis::Real,
    const std::atomic<bool>* cancel = nullptr);
