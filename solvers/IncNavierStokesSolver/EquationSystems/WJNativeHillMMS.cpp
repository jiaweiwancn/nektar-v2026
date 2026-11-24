// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupledMMS.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeWallV2.h>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Nektar::WJNative
{
namespace
{
// Spatial second jets and a first time derivative for the continuous
// manufactured solution. Neither mesh metrics nor assembled actions are read.
struct Jet
{
    double v, a, r, time, aa, ar, rr;
    Jet(double v=0., double a=0., double r=0., double time=0.,
        double aa=0., double ar=0., double rr=0.)
        : v(v), a(a), r(r), time(time), aa(aa), ar(ar), rr(rr)
    {
    }
    MMSJet Export() const { return {v,a,r,time,aa,ar,rr}; }
};
Jet operator+(const Jet &x,const Jet &y)
{
    return {x.v+y.v,x.a+y.a,x.r+y.r,x.time+y.time,
            x.aa+y.aa,x.ar+y.ar,x.rr+y.rr};
}
Jet operator*(const Jet &x,const Jet &y)
{
    return {x.v*y.v,x.a*y.v+x.v*y.a,x.r*y.v+x.v*y.r,
            x.time*y.v+x.v*y.time,x.aa*y.v+2*x.a*y.a+x.v*y.aa,
            x.ar*y.v+x.a*y.r+x.r*y.a+x.v*y.ar,x.rr*y.v+2*x.r*y.r+x.v*y.rr};
}
Jet operator-(const Jet &x,const Jet &y) { return x+(-1.)*y; }
Jet Chain(const Jet &x,double f,double df,double ddf)
{
    return {f,df*x.a,df*x.r,df*x.time,df*x.aa+ddf*x.a*x.a,
            df*x.ar+ddf*x.a*x.r,df*x.rr+ddf*x.r*x.r};
}
Jet operator/(const Jet &x,const Jet &y)
{
    const double inv=1/y.v;
    return x*Chain(y,inv,-inv*inv,2*inv*inv*inv);
}
Jet Root(const Jet &x)
{
    const double y=std::sqrt(x.v);
    return Chain(x,y,.5/y,-.25/(x.v*y));
}
Jet Sin(const Jet &x) { return Chain(x,std::sin(x.v),std::cos(x.v),-std::sin(x.v)); }
Jet Cos(const Jet &x) { return Chain(x,std::cos(x.v),-std::sin(x.v),-std::cos(x.v)); }
Jet DA(const Jet &x) { return {x.a,x.aa,x.ar}; }
Jet DR(const Jet &x) { return {x.r,x.ar,x.rr}; }

struct Graph
{
    Jet height, slope, second;
};
Graph Cosine(const Geometry &g,const Jet &a)
{
    const double frequency=std::acos(-1.)/g.normalMap->lengthX;
    // This fixture is wholly inside one smooth hill patch. Its exact field
    // continues across H in a virtual tube used only by the oracle.
    if (std::abs(a.v)>=g.normalMap->lengthX)
        throw std::runtime_error("manufactured normal foot left the smooth cosine patch");
    const double amplitude=g.normalMap->amplitude/2;
    const auto phase=frequency*a;
    return {amplitude*(1+Cos(phase)),-amplitude*frequency*Sin(phase),
            -amplitude*frequency*frequency*Cos(phase)};
}

struct TubePoint
{
    double a, r;
};
TubePoint VirtualTube(const Geometry &g,double a,double rho)
{
    if (rho<=g.H+32*std::numeric_limits<double>::epsilon()*g.H)
        return {a,rho};
    const auto graph=Cosine(g,a);
    const double slope=graph.slope.v, ell=std::hypot(1.,slope), d=g.d0+g.H;
    const double roof=g.normalMap->tunnelTop;
    const double fraction=(rho-g.H)/(roof-d);
    const double x=a-(1-fraction)*d*slope/ell;
    const double z=(1-fraction)*(graph.height.v+d/ell)+fraction*roof;
    double foot=a;
    for (int iteration=0;iteration<16;++iteration)
    {
        const auto f=Cosine(g,foot);
        const double dz=f.height.v-z, p=f.slope.v;
        const double residual=foot-x+dz*p;
        if (std::abs(residual)<=2e-15)
            return {foot,((x-foot)*(-p)-dz)/std::hypot(1.,p)-g.d0};
        foot-=residual/(1+p*p+dz*f.second.v);
    }
    throw std::runtime_error("manufactured virtual tube inverse failed");
}

struct State
{
    Jet speed, angleDerivative;
    Jet u, w, Pi, k, K, epsilon; // Velocity in the local (t,n) basis.
    Vector tangent, normal;
    Tensor gradient{}, gradientA{}, gradientR{};
    double r;
};
State TubeFields(const Reference &reference,const Geometry &g,TubePoint point,double time)
{
    const Jet a(point.a,1.), r(point.r,0.,1.), t(time,0.,0.,1.);
    const double H=g.H, R=.5;
    const auto graph=Cosine(g,a);
    const auto ell=Root(1+graph.slope*graph.slope);
    State f;
    f.r=point.r;
    f.angleDerivative=graph.second/(ell*ell);
    // speed is d(physical arclength at r)/da. Its radial derivative
    // supplies the geometric terms in continuity and scalar diffusion.
    f.speed=ell-(g.d0+r)*f.angleDerivative;
    f.tangent={1/ell.v,0.,graph.slope.v/ell.v};
    f.normal={-f.tangent[2],0.,f.tangent[0]};

    const Jet S=.02*(1+.2*Sin(.9*t)), V=.07*(1+.4*Cos(.85*t));
    const Jet X=a*a*(1-a), Xprime=2*a-3*a*a;
    const Jet Y=r*r*(R-r), Yprime=2*R*r-3*r*r;
    f.u=S*r+V*X*Yprime;
    f.w=(0.-V*Xprime*Y)/f.speed;

    const auto &p=reference.Parameters();
    const Jet rootK=std::sqrt(p.KIn)*(1+.012*Sin(.6*t)+.008*a*Cos(.8*t));
    f.K=rootK*rootK;
    f.k=f.K*(1+.005*(1+.2*Sin(.7*t))*(r*r-H*H)/(H*R));
    const Jet A=rootK*rootK*rootK/p.cLength;
    const Jet q=(r-H)/R;
    f.epsilon=A/H*(1+.03*Cos(.5*t)*q+.015*q*q);
    f.Pi=.0002*(1+.3*Sin(.8*t))*((a-.5)*(a-.5)+r*r/(R*R)+.15*a*r/R);

    const std::array<Jet,4> gradient{
        (DA(f.u)-f.angleDerivative*f.w)/f.speed,DR(f.u),
        (DA(f.w)+f.angleDerivative*f.u)/f.speed,DR(f.w)};
    for (int i=0;i<2;++i) for(int j=0;j<2;++j)
    {
        const auto &entry=gradient[2*i+j];
        f.gradient[2*i][2*j]=entry.v;
        f.gradientA[2*i][2*j]=entry.a;
        f.gradientR[2*i][2*j]=entry.r;
    }
    return f;
}

Vector Rotate(const State &f,double along,double normal)
{
    return {f.tangent[0]*along+f.normal[0]*normal,0.,
            f.tangent[2]*along+f.normal[2]*normal};
}

struct Coefficients
{
    Jet E, auxiliary, Dk, De;
};
Coefficients Transport(const Reference &ref,const State &s,bool outer)
{
    const auto &p=ref.Parameters();
    Coefficients result;
    if (outer)
        result.E=s.epsilon;
    else
    {
        const auto E=ref.DependentEpsilon(s.K.v,s.r);
        const auto A=ref.Auxiliary(s.K.v,s.r);
        result.E={E.value,E.dK*s.K.a,E.dr};
        result.auxiliary={A.value,A.dK*s.K.a,A.dr};
    }
    const Jet scale=p.Cd*s.k*s.k/result.E;
    result.Dk=p.nu+scale/p.sigmaK;
    result.De=p.nu+scale/p.sigmaE;
    return result;
}

double ScalarLeft(const State &s,const Jet &scalar,const Jet &diffusion)
{
    // div(D grad q) = (d_a(D q_a/B)+d_r(B D q_r))/B,
    // with B the virtual tube's arclength per original-floor abscissa.
    const double B=s.speed.v;
    const double divDiff=(diffusion.a*scalar.a+diffusion.v*scalar.aa-
        diffusion.v*scalar.a*s.speed.a/B)/(B*B) +
        diffusion.r*scalar.r+diffusion.v*scalar.rr+
        s.speed.r*diffusion.v*scalar.r/B;
    return scalar.time+s.u.v*scalar.a/B+s.w.v*scalar.r+
        scalar.v*(s.gradient[0][0]+s.gradient[2][2])-divDiff;
}
}

MMSFields CoupledMMS::NormalFields(double a,double rho,double time) const
{
    const auto s=TubeFields(m_reference,m_geometry,
                                          VirtualTube(m_geometry,a,rho),time);
    MMSFields f;
    const auto U=Rotate(s,s.u.v,s.w.v), Ut=Rotate(s,s.u.time,s.w.time);
    f.u.value=U[0];f.w.value=U[2];f.u.dt=Ut[0];f.w.dt=Ut[2];
    f.k=s.k.Export();f.Pi=s.Pi.Export();
    const std::array<Vector,2> axes{s.tangent,s.normal};
    for (int i : {0,2}) for (int j : {0,2})
        for(int l=0;l<2;++l) for(int m=0;m<2;++m)
            f.gradient[i][j]+=axes[l][i]*s.gradient[2*l][2*m]*axes[m][j];
    return f;
}

MMSJet CoupledMMS::NormalEpsilon(double a,double rho,double time) const
{
    return TubeFields(m_reference,m_geometry,
                                    VirtualTube(m_geometry,a,rho),time).epsilon.Export();
}

Vector CoupledMMS::NormalAWF(double a,double time) const
{
    const auto s=TubeFields(m_reference,m_geometry,{a,m_geometry.H},time);
    const double force=(s.Pi.a-(2./3)*s.K.a+2*s.u.v*s.u.a)/s.speed.v;
    return Rotate(s,force,0.);
}

double CoupledMMS::NormalFlux(double a,double rho,double time,const Vector &normal,
                             bool epsilon) const
{
    const auto s=TubeFields(m_reference,m_geometry,
                                          VirtualTube(m_geometry,a,rho),time);
    const auto c=Transport(m_reference,s,rho>=m_geometry.H-1e-14);
    const auto scalar=epsilon?s.epsilon:s.k;
    const auto gradient=Rotate(s,scalar.a/s.speed.v,scalar.r);
    return -(epsilon?c.De.v:c.Dk.v)*(normal[0]*gradient[0]+normal[2]*gradient[2]);
}

MMSForcing CoupledMMS::NormalForcing(double a,double rho,double time) const
{
    const auto s=TubeFields(m_reference,m_geometry,
                                          VirtualTube(m_geometry,a,rho),time);
    MMSForcing f;
    f.outer=rho>=m_geometry.H-1e-14;
    const auto c=Transport(m_reference,s,f.outer);
    const double nu=m_reference.Parameters().nu, H=m_geometry.H, B=s.speed.v;
    const Vector localNormal{0,0,1};
    const auto wall=AdaptAmended(s.k.v,c.E.v,s.gradient,localNormal,nu,c.auxiliary.v,s.r,H);
    const auto Va=ApplyAmendedDirection(s.k.v,c.E.v,s.gradient,localNormal,
        nu,c.auxiliary.v,s.r,H,s.k.a,c.E.a,s.gradientA,c.auxiliary.a).V;
    const auto Vr=ApplyAmendedDirection(s.k.v,c.E.v,s.gradient,localNormal,
        nu,c.auxiliary.v,s.r,H,s.k.r,c.E.r,s.gradientR,c.auxiliary.r,1.).V;

    // Differentiate the physical Cartesian rows of the stress. In the
    // moving basis the axial derivative has the additional commutator.
    const auto &V=wall.V;
    const double angle=s.angleDerivative.v;
    const double divT=(Va[0][0]-angle*(V[0][2]+V[2][0]))/B+Vr[0][2];
    const double divN=(Va[2][0]+angle*(V[0][0]-V[2][2]))/B+Vr[2][2];
    const auto &G=s.gradient;
    f.momentum=Rotate(s,s.u.time+G[0][0]*s.u.v+G[0][2]*s.w.v+s.Pi.a/B-divT,
                       s.w.time+G[2][0]*s.u.v+G[2][2]*s.w.v+s.Pi.r-divN);
    if (rho<m_geometry.h)
    {
        f.awf=NormalAWF(a,time);
        for (int i : {0,2}) f.momentum[i]-=f.awf[i];
    }
    f.production=wall.production;
    f.destruction=c.E.v;
    f.diffusionK=c.Dk.v;
    f.k=ScalarLeft(s,s.k,c.Dk)-wall.production+c.E.v;
    if (f.outer)
    {
        f.diffusionE=c.De.v;
        f.epsilon=ScalarLeft(s,s.epsilon,c.De)-(1.44*wall.production-1.82*c.E.v)*c.E.v/s.k.v;
    }
    return f;
}
}
