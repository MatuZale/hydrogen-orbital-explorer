#include "flow.h"
#include "sampling.h"
#include <cmath>
#include <stdexcept>

namespace {
constexpr double pi=3.141592653589793238462643383279502884;
float normalized_speed(int n,double rho) {
    // v/v_ref = (n^2/2)/rho. The visual time scale cancels in this ratio.
    return rho>0 ? static_cast<float>(-std::expm1(-0.5*n*n/rho)) : 1.0f;
}
}

std::vector<FlowParticle> generate_flow_particles(int n,int l,int m_quantum,
        int count,unsigned seed,const std::atomic<bool>* cancel) {
    validate_orbital(n,l,m_quantum);
    if (count<0) throw std::invalid_argument("Particle count must be nonnegative");
    if (m_quantum==0) return {};
    if (count==0) count=25000;
    const auto points=generate_orbital(n,l,m_quantum,count,seed,OrbitalBasis::Complex,cancel);
    std::vector<FlowParticle> result;
    result.reserve(points.size());
    for (const auto& p:points) {
        if (cancel && cancel->load(std::memory_order_relaxed)) return {};
        result.push_back({p.x,p.y,p.z,normalized_speed(n,std::hypot(p.x,p.y))});
    }
    return result;
}

void update_flow_particles(std::vector<FlowParticle>& pts,int m_quantum,int n,
        float dt,float speed_mult) {
    if (n<1 || n>MAX_ORBITAL_N || m_quantum<=-n || m_quantum>=n ||
        !std::isfinite(dt) || dt<0 || !std::isfinite(speed_mult) || speed_mult<0)
        throw std::invalid_argument("Invalid flow quantum numbers or time step");
    if (m_quantum==0 || dt==0 || speed_mult==0) return;
    const double scale=2*pi*static_cast<double>(n)*n*n*n/(std::abs(m_quantum)*5.0);
    for (auto& p:pts) {
        const double x=p.x,y=p.y,rho2=x*x+y*y;
        if (rho2==0) continue; // the axis has zero complex-state weight for m!=0
        // The ODE has constant rho and z; its exact solution is a rotation.
        // Reducing the angle improves trig evaluation, without clipping motion.
        const double angle=std::remainder(m_quantum*scale*dt*speed_mult/rho2,2*pi);
        const double c=std::cos(angle),s=std::sin(angle);
        p.x=static_cast<float>(x*c-y*s);
        p.y=static_cast<float>(x*s+y*c);
        p.speed_norm=normalized_speed(n,std::sqrt(rho2));
    }
}
