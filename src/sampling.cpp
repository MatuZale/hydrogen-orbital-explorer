#include "sampling.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <random>
#include <stdexcept>

namespace {
constexpr double pi=3.141592653589793238462643383279502884;
bool canceled(const std::atomic<bool>* flag) {
    return flag && flag->load(std::memory_order_relaxed);
}
// The table only brackets a quantile. Every sample is then refined against
// the analytic CDF, so table interpolation does not change the distribution.
class RadialSampler {
    static constexpr std::size_t intervals=512;
    int n_,l_;
    double radius_;
    std::array<double,intervals+1> cdf_{};
public:
    RadialSampler(int n,int l): n_(n),l_(l),radius_(orbital_radius(n)) {
        for (std::size_t i=0;i<=intervals;++i)
            cdf_[i]=radial_cdf(n,l,radius_*i/intervals);
    }
    double radius(double unit) const {
        const double target=unit*cdf_.back();
        const auto upper=std::lower_bound(cdf_.begin()+1,cdf_.end(),target);
        const std::size_t i=static_cast<std::size_t>(upper-cdf_.begin());
        double lo=radius_*(i-1)/intervals,hi=radius_*i/intervals;
        for (int step=0;step<32;++step) {
            const double mid=(lo+hi)/2;
            if (radial_cdf(n_,l_,mid)<target) lo=mid; else hi=mid;
        }
        return (lo+hi)/2;
    }
};
}

std::vector<PointData> generate_orbital(int n,int l,int m,int target_points,
        unsigned seed,OrbitalBasis basis,const std::atomic<bool>* cancel) {
    validate_orbital(n,l,m);
    if (target_points<0 || (basis!=OrbitalBasis::Real && basis!=OrbitalBasis::Complex))
        throw std::invalid_argument("Invalid point count or orbital basis");
    if (canceled(cancel)) return {};
    if (target_points==0) target_points=n<=2 ? 70000 : (n==3 ? 55000 : 40000);
    const RadialSampler radial(n,l);
    std::mt19937_64 rng(seed);
    auto uniform=[&] { return std::generate_canonical<double,53>(rng); };
    // Addition theorem: sum over an orthonormal harmonic basis is this
    // constant, hence each nonnegative summand is bounded by it. The same
    // bound applies to real and complex harmonics without a sampled maximum.
    const double angular_bound=(2.0*l+1)/(4*pi);
    std::vector<PointData> result;
    result.reserve(static_cast<std::size_t>(target_points));
    std::size_t attempts=0;
    while (result.size()<static_cast<std::size_t>(target_points)) {
        if ((attempts++%128)==0 && canceled(cancel)) return {};
        const double cosine=2*uniform()-1;
        const double theta=std::acos(cosine),phi=2*pi*uniform()-pi;
        if (uniform()*angular_bound>=angular_prob(l,m,theta,phi,basis)) continue;
        const double r=radial.radius(uniform());
        const double planar=r*std::sqrt(std::max(0.0,1-cosine*cosine));
        const float sign=basis==OrbitalBasis::Real
            ? wavefunction_sign(n,l,m,r,theta,phi) : 1.0f;
        result.push_back({static_cast<float>(planar*std::cos(phi)),
                          static_cast<float>(planar*std::sin(phi)),
                          static_cast<float>(r*cosine),sign});
    }
    if (canceled(cancel)) return {};
    return result;
}
