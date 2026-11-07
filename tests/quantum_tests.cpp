#include "quantum.h"
#include "sampling.h"
#include "flow.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstring>
#include <functional>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
constexpr double pi=3.141592653589793238462643383279502884;
int checks=0;
void check(bool ok,const std::string& label) {
    ++checks;
    if (!ok) throw std::runtime_error(label);
}
void near(double actual,double expected,double tolerance,const std::string& label) {
    check(std::isfinite(actual) && std::abs(actual-expected)<=tolerance,
          label+": actual="+std::to_string(actual)+", expected="+std::to_string(expected));
}
template<class Function> void invalid(Function function,const std::string& label) {
    bool caught=false;
    try { function(); } catch (const std::invalid_argument&) { caught=true; }
    check(caught,label);
}
double simpson(const std::function<double(double)>& f,double end,int intervals=16384) {
    const double h=end/intervals;
    double sum=f(0)+f(end);
    for (int i=1;i<intervals;++i) sum+=(i%2 ? 4 : 2)*f(i*h);
    return sum*h/3;
}
// Independent Gauss-Legendre quadrature on cos(theta), uniform trapezoidal phi.
// This integrates the low-degree angular densities exactly to roundoff.
double sphere_integral(const std::function<double(double,double)>& f) {
    constexpr int order=16,azimuth=32;
    double sum=0;
    for (int i=0;i<order;++i) {
        double x=std::cos(pi*(i+0.75)/(order+0.5)),derivative=0;
        for (int iteration=0;iteration<20;++iteration) {
            double older=1,previous=x;
            for (int k=2;k<=order;++k) {
                const double next=((2*k-1)*x*previous-(k-1)*older)/k;
                older=previous;previous=next;
            }
            derivative=order*(x*previous-older)/(x*x-1);
            const double correction=previous/derivative;
            x-=correction;
            if (std::abs(correction)<1e-15) break;
        }
        const double weight=2/((1-x*x)*derivative*derivative);
        for (int j=0;j<azimuth;++j)
            sum+=weight*(2*pi/azimuth)*f(std::acos(x),2*pi*(j+0.5)/azimuth);
    }
    return sum;
}
double radius(const PointData& p) { return std::hypot(std::hypot(p.x,p.y),p.z); }

void test_math() {
    near(factorial(0),1,0,"0!"); near(factorial(7),5040,0,"7!");
    near(assoc_laguerre(2,1,2),-1,1e-15,"L2^1(2)");
    near(assoc_legendre(2,0,0),-0.5,1e-15,"P2(0)");
    near(assoc_legendre(2,2,0),3,1e-15,"positive no-CS P2^2");
    near(sph_harm_real(1,1,pi/2,0),std::sqrt(3/(4*pi)),1e-15,"px convention");
    near(sph_harm_real(1,-1,pi/2,pi/2),std::sqrt(3/(4*pi)),1e-15,"py convention");
    near(psi_prob(1,0,0,0,0,0),1/pi,1e-15,"normalized 1s origin");
    near(psi_prob(2,1,1,0,0,0),0,0,"p origin node");
    near(radial_amplitude(2,0,2),0,0,"2s radial node");
    near(radial_amplitude(3,1,6),0,0,"3p radial node");
    near(radial_amplitude(3,0,1.5*(3-std::sqrt(3.0))),0,1e-15,"3s inner radial node");
    check(wavefunction_sign(2,0,0,1,1,1)==1,"2s inner phase");
    check(wavefunction_sign(2,0,0,3,1,1)==-1,"2s outer phase");
    check(wavefunction_sign(3,1,1,3,pi/2,0)==1,"3p inner phase");
    check(wavefunction_sign(3,1,1,8,pi/2,0)==-1,"3p outer phase");
    check(wavefunction_sign(3,1,1,8,pi/2,pi)==1,"radial and angular signs multiply");

    for (int l=0;l<=3;++l) {
        for (int m=-l;m<=l;++m) {
            for (auto basis:{OrbitalBasis::Real,OrbitalBasis::Complex}) {
                near(sphere_integral([&](double t,double p){return angular_prob(l,m,t,p,basis);}),
                     1,2e-13,"angular normalization");
                for (int ti=0;ti<=30;++ti) for (int pj=0;pj<32;++pj) {
                    const double prob=angular_prob(l,m,pi*ti/30,2*pi*pj/32,basis);
                    check(prob>=0 && prob<=(2*l+1)/(4*pi)+1e-14,"rigorous angular envelope");
                }
            }
            for (int other=-l;other<m;++other) {
                near(sphere_integral([&](double t,double p){return sph_harm_real(l,m,t,p)*sph_harm_real(l,other,t,p);}),
                     0,2e-13,"real harmonic orthogonality");
            }
        }
        double total=0;
        for (int m=-l;m<=l;++m) total+=angular_prob(l,m,0.73,1.14,OrbitalBasis::Real);
        near(total,(2*l+1)/(4*pi),1e-14,"addition theorem");
    }

    for (double r:{0.001,0.1,1.0,3.0,10.0}) {
        const double expected=1-std::exp(-2*r)*(1+2*r+2*r*r);
        near(radial_cdf(1,0,r),expected,3e-16,"1s exact CDF");
    }
    near(radial_quantile(1,0,1e-60),std::cbrt(0.75e-60),1e-32,
         "tiny quantile preserves relative precision");
    for (int n=1;n<=4;++n) {
        const double cutoff=orbital_radius(n);
        std::cout<<"n="<<n<<" cutoff="<<cutoff<<" a0\n";
        for (int l=0;l<n;++l) {
            check(radial_cdf(n,l,cutoff)>=ORBITAL_CONTAINED_PROBABILITY,"cutoff contains >=99.9%");
            near(simpson([&](double r){return radial_density(n,l,r);},cutoff),
                 radial_cdf(n,l,cutoff),3e-10,"radial CDF vs independent quadrature");
            for (double fraction:{0.01,0.1,0.25,0.5}) {
                const double r=cutoff*fraction;
                near(radial_cdf(n,l,r),simpson([&](double x){return radial_density(n,l,x);},r),
                     3e-10,"CDF partial integral");
            }
            double previous=0;
            for (int i=0;i<=4000;++i) {
                const double value=radial_cdf(n,l,cutoff*i/4000);
                check(value>=previous-2e-13 && value<=1,"CDF monotonicity");
                previous=value;
            }
            for (double p:{1e-7,0.001,0.1,0.5,0.9,0.999})
                near(radial_cdf(n,l,radial_quantile(n,l,p)),p,2e-13,"CDF inverse");
            near(radial_cdf(n,l,std::numeric_limits<double>::infinity()),1,0,"CDF infinity");
            near(radial_density(n,l,std::numeric_limits<double>::infinity()),0,0,"density infinity");
        }
    }
    check(orbital_radius(4)>41,"4s old radius omitted excessive tail");
}

void test_sampling() {
    const auto a=generate_orbital(2,0,0,3000,71),b=generate_orbital(2,0,0,3000,71);
    check(a.size()==3000 && b.size()==a.size(),"exact sample count");
    check(std::memcmp(a.data(),b.data(),a.size()*sizeof(PointData))==0,"repeatable seeded sampler");
    const auto c=generate_orbital(2,0,0,3000,72);
    check(std::memcmp(a.data(),c.data(),a.size()*sizeof(PointData))!=0,"distinct seed");
    int positive=0,negative=0;
    for (const auto& p:a) {
        const double r=radius(p);
        check(p.sign==(r<2 ? 1.0f : -1.0f),"sample includes full radial sign");
        if (p.sign>0) ++positive; else ++negative;
    }
    check(positive>0 && negative>0,"both 2s phases present");

    for (const auto& orbital:std::vector<OrbitalParams>{{1,0,0,"1s"},{4,0,0,"4s"},{4,3,-2,"4f"}}) {
        const int count=orbital.n==4 && orbital.l==0 ? 100000 : 20000;
        auto points=generate_orbital(orbital.n,orbital.l,orbital.m,count,1337);
        std::vector<double> cdf;
        cdf.reserve(points.size());
        int inner=0;
        const double cutoff=orbital_radius(orbital.n),mass=radial_cdf(orbital.n,orbital.l,cutoff);
        for (const auto& p:points) {
            const double r=radius(p);
            check(r<=cutoff+1e-5,"sample inside display radius");
            if (r<1) ++inner;
            cdf.push_back(radial_cdf(orbital.n,orbital.l,r)/mass);
        }
        std::sort(cdf.begin(),cdf.end());
        double ks=0;
        for (std::size_t i=0;i<cdf.size();++i)
            ks=std::max({ks,std::abs(cdf[i]-static_cast<double>(i)/count),
                        std::abs(cdf[i]-static_cast<double>(i+1)/count)});
        check(ks<0.016,"radial distribution KS bound");
        std::cout<<orbital.name<<" radial KS="<<ks<<"\n";
        if (orbital.n==4 && orbital.l==0) {
            const double expected=radial_cdf(4,0,1)/mass;
            near(static_cast<double>(inner)/count,expected,0.0007,"4s inner shell not clipped by envelope");
            std::cout<<"4s P(r<1): sample="<<static_cast<double>(inner)/count<<", analytic="<<expected<<"\n";
        }
    }
    for (auto basis:{OrbitalBasis::Real,OrbitalBasis::Complex}) {
        const auto points=generate_orbital(2,1,1,30000,927,basis);
        double xx=0,yy=0,zz=0;
        for (const auto& p:points) {
            const double r2=static_cast<double>(p.x)*p.x+static_cast<double>(p.y)*p.y+static_cast<double>(p.z)*p.z;
            xx+=p.x*p.x/r2;yy+=p.y*p.y/r2;zz+=p.z*p.z/r2;
            if (basis==OrbitalBasis::Complex) check(p.sign==1,"neutral complex cloud color");
        }
        near(xx/points.size(),basis==OrbitalBasis::Real ? 0.6 : 0.4,0.012,"basis x angular moment");
        near(yy/points.size(),basis==OrbitalBasis::Real ? 0.2 : 0.4,0.012,"basis y angular moment");
        near(zz/points.size(),0.2,0.012,"basis z angular moment");
    }
    for (int i=0;i<NUM_ORBITALS;++i) {
        const auto o=ORBITALS[i];
        check(generate_orbital(o.n,o.l,o.m,32,i).size()==32,"catalog sample valid");
    }
    const std::atomic<bool> canceled{true};
    check(generate_orbital(4,3,2,0,42,OrbitalBasis::Real,&canceled).empty(),"cloud cancellation");
    check(generate_flow_particles(4,3,2,0,42,&canceled).empty(),"flow cancellation");
}

void test_flow() {
    check(generate_flow_particles(2,1,0,0).empty(),"m0 has zero flow");
    const auto points=generate_flow_particles(2,1,1,0,123);
    check(points.size()==25000,"flow zero count selects default");
    const auto matching=generate_orbital(2,1,1,25000,123,OrbitalBasis::Complex);
    for (std::size_t i=0;i<points.size();++i) {
        check(points[i].x==matching[i].x && points[i].y==matching[i].y && points[i].z==matching[i].z,
              "flow seeds match complex density");
        check(points[i].speed_norm>=0 && points[i].speed_norm<=1,"flow colormap range");
    }
    const std::vector<FlowParticle> initial{{0.5f,0,1,0},{2,1,-3,0},{0,0,2,0},{4,0,0,0}};
    auto at30=initial,at60=initial,at120=initial;
    for (int i=0;i<30;++i) update_flow_particles(at30,1,2,1.0f/30);
    for (int i=0;i<60;++i) update_flow_particles(at60,1,2,1.0f/60);
    for (int i=0;i<120;++i) update_flow_particles(at120,1,2,1.0f/120);
    for (std::size_t i=0;i<initial.size();++i) {
        near(at30[i].x,at60[i].x,0.0002,"30/60 FPS x invariant");
        near(at30[i].y,at60[i].y,0.0002,"30/60 FPS y invariant");
        near(at30[i].x,at120[i].x,0.0002,"30/120 FPS x invariant");
        near(at30[i].y,at120[i].y,0.0002,"30/120 FPS y invariant");
        near(std::hypot(at120[i].x,at120[i].y),std::hypot(initial[i].x,initial[i].y),2e-6,"flow preserves radius");
        near(at120[i].z,initial[i].z,0,"flow preserves z");
    }
    auto single=initial;
    update_flow_particles(single,1,2,1);
    const double angle=2*pi*16/(5*0.25);
    near(single[0].x,0.5*std::cos(angle),1e-7,"unclamped exact x rotation");
    near(single[0].y,0.5*std::sin(angle),1e-7,"unclamped exact y rotation");
    auto negative=initial;
    update_flow_particles(negative,-1,2,1);
    near(negative[0].x,single[0].x,1e-7,"negative m reverses orientation x");
    near(negative[0].y,-single[0].y,1e-7,"negative m reverses orientation y");
    auto reference=initial;
    update_flow_particles(reference,1,2,5);
    near(reference[3].x,4,1e-6,"five second reference turn x");
    near(reference[3].y,0,1e-6,"five second reference turn y");
}
void test_invalid() {
    invalid([]{factorial(-1);},"invalid factorial");
    invalid([]{generate_orbital(0,0,0,1);},"n0 rejected");
    invalid([]{generate_orbital(5,0,0,1);},"outside catalog rejected");
    invalid([]{generate_orbital(2,2,0,1);},"invalid l rejected");
    invalid([]{generate_orbital(2,1,2,1);},"invalid m rejected");
    invalid([]{generate_orbital(2,1,1,-1);},"negative count rejected");
    invalid([]{generate_flow_particles(2,1,0,-1);},"negative flow count rejected");
    invalid([]{radial_quantile(2,0,1.1);},"invalid quantile rejected");
    invalid([]{radial_cdf(1,0,std::numeric_limits<double>::quiet_NaN());},"NaN radius rejected");
    invalid([]{std::vector<FlowParticle> p;update_flow_particles(p,1,2,-1);},"negative dt rejected");
}
}
int main() {
    try {
        test_math(); test_sampling(); test_flow(); test_invalid();
        std::cout<<"Passed "<<checks<<" checks.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr<<"FAIL after "<<checks<<" checks: "<<error.what()<<"\n";
        return 1;
    }
}
