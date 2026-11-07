#include "quantum.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace {
constexpr double pi = 3.141592653589793238462643383279502884;
void validate_angular(int l, int m) {
    if (l < 0 || l > 10 || m < -l || m > l)
        throw std::invalid_argument("Invalid angular quantum numbers");
}
void validate_radial(int n, int l) { validate_orbital(n, l, 0); }

// In x=2r/n, the normalized radial PDF is exp(-x) P(x)/J, with
// P=x^(2l+2) [L_(n-l-1)^(2l+1)(x)]^2 and J=2n(n+l)!/(n-l-1)!.
// Integrating each monomial gives the finite upper-tail polynomial Q.
struct RadialIntegral {
    std::array<long double, 9> density{};
    std::array<long double, 9> tail{};
    int degree = 0;
};
RadialIntegral make_integral(int n, int l) {
    RadialIntegral out;
    out.degree = 2*n;
    const int k = n-l-1, alpha = 2*l+1;
    const long double norm = 2.0L*n*factorial(n+l)/factorial(k);
    std::array<long double, 4> lag{};
    for (int j=0; j<=k; ++j) {
        lag[j] = factorial(k+alpha) /
            (static_cast<long double>(factorial(k-j))*factorial(alpha+j)*factorial(j));
        if (j%2) lag[j] = -lag[j];
    }
    for (int i=0; i<=k; ++i)
        for (int j=0; j<=k; ++j)
            out.density[2*l+2+i+j] += lag[i]*lag[j]/norm;
    for (int i=0; i<=out.degree; ++i)
        for (int j=0; j<=i; ++j)
            out.tail[j] += out.density[i]*factorial(i)/factorial(j);
    return out;
}
const RadialIntegral& integral(int n, int l) {
    static const auto table = [] {
        std::array<std::array<RadialIntegral, 4>, 4> t{};
        for (int nn=1; nn<=4; ++nn)
            for (int ll=0; ll<nn; ++ll) t[nn-1][ll]=make_integral(nn,ll);
        return t;
    }();
    return table[n-1][l];
}
long double evaluate_cdf(const RadialIntegral& data, long double x) {
    if (x <= 0) return 0;
    if (x > 10000) return 1; // avoid infinity*zero in the upper-tail form
    long double value = 0;
    if (x < 1) {
        // Lower incomplete gamma series avoids 1-tail cancellation near zero.
        for (int i=0; i<=data.degree; ++i) {
            if (data.density[i] == 0) continue;
            long double term=1, series=1;
            for (int j=1; j<100; ++j) {
                term *= x/(i+1+j);
                series += term;
                if (std::abs(term) < std::abs(series)*1e-19L) break;
            }
            value += data.density[i]*std::pow(x,i+1)/(i+1)*series;
        }
        value *= std::exp(-x);
    } else {
        long double q=data.tail[data.degree];
        for (int j=data.degree-1; j>=0; --j) q=q*x+data.tail[j];
        value=1-std::exp(-x)*q;
    }
    return std::clamp(value, 0.0L, 1.0L);
}
}

const OrbitalParams ORBITALS[] = {
    // n = 1
    {1, 0,  0, "1s"},
    // n = 2
    {2, 0,  0, "2s"},
    {2, 1,  0, "2p_z"},
    {2, 1,  1, "2p_x"},
    {2, 1, -1, "2p_y"},
    // n = 3
    {3, 0,  0, "3s"},
    {3, 1,  0, "3p_z"},
    {3, 1,  1, "3p_x"},
    {3, 1, -1, "3p_y"},
    {3, 2,  0, "3d_z2"},
    {3, 2,  1, "3d_xz"},
    {3, 2, -1, "3d_yz"},
    {3, 2,  2, "3d_x2-y2"},
    {3, 2, -2, "3d_xy"},
    // n = 4
    {4, 0,  0, "4s"},
    {4, 1,  0, "4p_z"},
    {4, 1,  1, "4p_x"},
    {4, 1, -1, "4p_y"},
    {4, 2,  0, "4d_z2"},
    {4, 2,  1, "4d_xz"},
    {4, 2, -1, "4d_yz"},
    {4, 2,  2, "4d_x2-y2"},
    {4, 2, -2, "4d_xy"},
    {4, 3,  0, "4f_z3"},
    {4, 3,  1, "4f_xz2"},
    {4, 3, -1, "4f_yz2"},
    {4, 3,  2, "4f_z(x2-y2)"},
    {4, 3, -2, "4f_xyz"},
    {4, 3,  3, "4f_x(x2-3y2)"},
    {4, 3, -3, "4f_y(3x2-y2)"},
};

const int NUM_ORBITALS = static_cast<int>(sizeof(ORBITALS) / sizeof(ORBITALS[0]));

void validate_orbital(int n, int l, int m) {
    if (n < 1 || n > MAX_ORBITAL_N || l < 0 || l >= n || m < -l || m > l)
        throw std::invalid_argument("Supported orbitals require 1<=n<=4, 0<=l<n, |m|<=l");
}

double factorial(int n) {
    if (n < 0 || n > 20) throw std::invalid_argument("factorial requires 0..20");
    double result=1;
    for (int i=2; i<=n; ++i) result*=i;
    return result;
}

double assoc_laguerre(int n, int alpha, double x) {
    if (n < 0 || alpha < 0 || !std::isfinite(x))
        throw std::invalid_argument("Invalid Laguerre arguments");
    if (n == 0) return 1;
    double older=1, previous=1.0+alpha-x;
    for (int k=2; k<=n; ++k) {
        const double next=((2*k-1+alpha-x)*previous-(k-1+alpha)*older)/k;
        older=previous; previous=next;
    }
    return previous;
}

double assoc_legendre(int l, int m, double x) {
    if (l < 0 || m < 0 || m > l || !std::isfinite(x) || x < -1 || x > 1)
        throw std::invalid_argument("Invalid Legendre arguments");
    double pmm=1;
    const double s=std::sqrt(std::max(0.0,1-x*x));
    for (int i=1; i<=m; ++i) pmm *= (2*i-1)*s;
    if (l == m) return pmm;
    double previous=x*(2*m+1)*pmm;
    for (int k=m+2; k<=l; ++k) {
        const double next=((2*k-1)*x*previous-(k+m-1)*pmm)/(k-m);
        pmm=previous; previous=next;
    }
    return previous;
}

double radial_amplitude(int n, int l, double r) {
    validate_radial(n,l);
    if (std::isnan(r)) throw std::invalid_argument("Radius is NaN");
    if (r < 0 || std::isinf(r)) return 0;
    const double x=2*r/n;
    // Beyond this range the exponential is negligible and polynomial
    // intermediates should not be allowed to overflow for arbitrary inputs.
    if (x > 10000) return 0;
    const double norm=std::sqrt(4.0/(n*n*n*n)*factorial(n-l-1)/factorial(n+l));
    return norm*std::exp(-x/2)*std::pow(x,l)*assoc_laguerre(n-l-1,2*l+1,x);
}

double radial_prob(int n, int l, double r) {
    const double amplitude=radial_amplitude(n,l,r);
    return amplitude*amplitude;
}

double radial_density(int n, int l, double r) {
    const double amplitude=radial_amplitude(n,l,r);
    if (amplitude == 0) return 0;
    return (r*amplitude)*(r*amplitude);
}

double radial_cdf(int n, int l, double r) {
    validate_radial(n,l);
    if (std::isnan(r)) throw std::invalid_argument("Radius is NaN");
    return static_cast<double>(evaluate_cdf(integral(n,l),2.0L*r/n));
}

double radial_quantile(int n, int l, double p) {
    validate_radial(n,l);
    if (!std::isfinite(p) || p < 0 || p > 1)
        throw std::invalid_argument("Probability must be in [0,1]");
    if (p == 0) return 0;
    if (p == 1) return std::numeric_limits<double>::infinity();
    double lo=0, hi=n*n;
    while (radial_cdf(n,l,hi) < p) hi*=2;
    // Relative bracketing also resolves probabilities far below the RNG's
    // normal range; a fixed [0,n^2] interval loses tiny quantiles in 60 steps.
    while (radial_cdf(n,l,hi/2) >= p) hi/=2;
    for (int i=0; i<60; ++i) {
        const double mid=(lo+hi)/2;
        if (radial_cdf(n,l,mid) < p) lo=mid; else hi=mid;
    }
    return hi;
}

double sph_harm_real(int l, int m, double theta, double phi) {
    validate_angular(l,m);
    if (!std::isfinite(theta) || !std::isfinite(phi))
        throw std::invalid_argument("Angles must be finite");
    const int a=std::abs(m);
    const double p=assoc_legendre(l,a,std::cos(theta));
    const double norm=std::sqrt((2.0*l+1)/(4*pi)*factorial(l-a)/factorial(l+a));
    if (m == 0) return norm*p;
    return std::sqrt(2.0)*norm*p*(m>0 ? std::cos(m*phi) : std::sin(a*phi));
}

double angular_prob(int l, int m, double theta, double phi, OrbitalBasis basis) {
    if (basis == OrbitalBasis::Real) {
        const double y=sph_harm_real(l,m,theta,phi);
        return y*y;
    }
    validate_angular(l,m);
    if (basis != OrbitalBasis::Complex || !std::isfinite(theta) || !std::isfinite(phi))
        throw std::invalid_argument("Invalid basis or angles");
    const int a=std::abs(m);
    const double p=assoc_legendre(l,a,std::cos(theta));
    return (2.0*l+1)/(4*pi)*factorial(l-a)/factorial(l+a)*p*p;
}

namespace {
double cartesian_prob(int n,int l,int m,double x,double y,double z,OrbitalBasis basis) {
    validate_orbital(n,l,m);
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z))
        throw std::invalid_argument("Coordinates must be finite");
    const double r=std::hypot(std::hypot(x,y),z);
    const double theta=r>0 ? std::acos(std::clamp(z/r,-1.0,1.0)) : 0;
    return radial_prob(n,l,r)*angular_prob(l,m,theta,std::atan2(y,x),basis);
}
}
double psi_prob(int n,int l,int m,double x,double y,double z) {
    return cartesian_prob(n,l,m,x,y,z,OrbitalBasis::Real);
}
double psi_complex_prob(int n,int l,int m,double x,double y,double z) {
    return cartesian_prob(n,l,m,x,y,z,OrbitalBasis::Complex);
}
float wavefunction_sign(int n,int l,int m,double r,double theta,double phi) {
    validate_orbital(n,l,m);
    return radial_amplitude(n,l,r)*sph_harm_real(l,m,theta,phi) < 0 ? -1.0f : 1.0f;
}
float orbital_sign(int l,int m,double theta,double phi) {
    return sph_harm_real(l,m,theta,phi) < 0 ? -1.0f : 1.0f;
}

double orbital_radius(int n) {
    validate_radial(n,0);
    static const auto radii=[] {
        std::array<double,4> result{};
        for (int nn=1; nn<=4; ++nn)
            for (int ll=0; ll<nn; ++ll)
                result[nn-1]=std::max(result[nn-1],radial_quantile(nn,ll,
                    ORBITAL_CONTAINED_PROBABILITY+1e-10));
        return result;
    }();
    return radii[n-1];
}
