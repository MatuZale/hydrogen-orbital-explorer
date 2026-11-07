#pragma once

// Hydrogen eigenstates in atomic units (a0 = hbar = electron mass = 1).
// The public orbital API supports the catalogued n = 1..4 states.
enum class OrbitalBasis { Real, Complex };
constexpr int MAX_ORBITAL_N = 4;
constexpr double ORBITAL_CONTAINED_PROBABILITY = 0.999;

struct OrbitalParams { int n, l, m; const char* name; };
extern const OrbitalParams ORBITALS[];
extern const int NUM_ORBITALS;

// Invalid quantum numbers throw std::invalid_argument.
void validate_orbital(int n, int l, int m);
double factorial(int n); // 0..20
double assoc_laguerre(int n, int alpha, double x);
// Associated Legendre WITHOUT the Condon–Shortley phase (m >= 0).
double assoc_legendre(int l, int m, double x);

// Normalized signed radial amplitude R, its square, and the probability
// density in dr (r^2 R^2). The Jacobian belongs only in radial_density.
double radial_amplitude(int n, int l, double r);
double radial_prob(int n, int l, double r);
double radial_density(int n, int l, double r);
// Analytic integral of r^2 R^2, evaluated as a finite polynomial times exp.
double radial_cdf(int n, int l, double r);
// Inverse CDF, p in [0,1); p=1 has an infinite radius.
double radial_quantile(int n, int l, double p);

// Orthonormal real tesseral harmonics: m>0 cosine, m<0 sine.
// This real basis is not an Lz eigenbasis for m != 0.
double sph_harm_real(int l, int m, double theta, double phi);
double angular_prob(int l, int m, double theta, double phi, OrbitalBasis basis);
// Normalized Cartesian probability densities in the two bases.
double psi_prob(int n, int l, int m, double x, double y, double z);
double psi_complex_prob(int n, int l, int m, double x, double y, double z);
// Full real wavefunction sign, including changes across radial nodes.
float wavefunction_sign(int n, int l, int m, double r, double theta, double phi);
// Angular sign only; retained as a math helper, not used for cloud colors.
float orbital_sign(int l, int m, double theta, double phi);

// A shared display radius for all l at fixed n. Computed from the analytic
// radial CDF at p=0.999 + 1e-10, so each truncated cloud retains at least
// 99.9% to the tested floating-point accuracy. Sampling is conditional on
// r <= this radius; the small omitted tail is not drawn.
double orbital_radius(int n);
