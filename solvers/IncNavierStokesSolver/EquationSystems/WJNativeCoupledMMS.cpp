// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupledMMS.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeWallV2.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Nektar::WJNative
{
namespace
{
constexpr double Ce1=1.44, Ce2=1.82;

void Require(bool value, const char *message)
{
    if (!value)
    {
        throw std::runtime_error(message);
    }
}

void CheckTangentialTime(double x, double time)
{
    Require(std::isfinite(x) && x>=0 && x<=1 && std::isfinite(time),
            "coupled MMS requires finite time and x in [0,1]");
}

bool OuterOwner(double r, double H)
{
    // Preserve H+ ownership when subtracting the fixture's d0 perturbs an
    // exactly supplied physical interface coordinate by arithmetic roundoff.
    return r>=H-64*std::numeric_limits<double>::epsilon()*H;
}

void CheckNormal(const Vector &normal)
{
    double length=0;
    for (double value : normal)
    {
        Require(std::isfinite(value), "nonfinite coupled MMS boundary normal");
        length+=value*value;
    }
    Require(std::abs(length-1)<1e-12 && std::abs(normal[1])<1e-12,
            "coupled MMS needs a planar physical unit normal");
}

struct Transport
{
    double E=0, Ex=0, Er=0;
    double auxiliary=0, auxiliaryX=0, auxiliaryR=0;
    double Dk=0, Dkx=0, Dkr=0, De=0, Dex=0, Der=0;
};

Transport Coefficients(const Reference &reference, const MMSFields &fields,
                       const MMSJet &K, double r, const MMSJet *outerEpsilon)
{
    const auto &p=reference.Parameters();
    Transport result;
    if (outerEpsilon)
    {
        result.E=outerEpsilon->value;
        result.Ex=outerEpsilon->dx;
        result.Er=outerEpsilon->dr;
    }
    else
    {
        const auto dependent=reference.DependentEpsilon(K.value,r);
        result.E=dependent.value;
        result.Ex=dependent.dK*K.dx;
        result.Er=dependent.dr;
    }
    const auto auxiliary=reference.Auxiliary(K.value,r);
    result.auxiliary=auxiliary.value;
    result.auxiliaryX=auxiliary.dK*K.dx;
    result.auxiliaryR=auxiliary.dr;
    Require(fields.k.value>0 && result.E>0,
            "nonpositive manufactured scalar/destruction");
    // Actual transported k is used here, not the counterfactual column's k_ref.
    const double k=fields.k.value;
    const double quotient=k*k/result.E;
    const double quotientX=2*k*fields.k.dx/result.E-
                            k*k*result.Ex/(result.E*result.E);
    const double quotientR=2*k*fields.k.dr/result.E-
                            k*k*result.Er/(result.E*result.E);
    result.Dk=p.nu+p.Cd*quotient/p.sigmaK;
    result.Dkx=p.Cd*quotientX/p.sigmaK;
    result.Dkr=p.Cd*quotientR/p.sigmaK;
    if (outerEpsilon)
    {
        result.De=p.nu+p.Cd*quotient/p.sigmaE;
        result.Dex=p.Cd*quotientX/p.sigmaE;
        result.Der=p.Cd*quotientR/p.sigmaE;
    }
    Require(std::isfinite(result.Dk) && result.Dk>0 &&
                std::isfinite(result.Dkx) && std::isfinite(result.Dkr) &&
                (!outerEpsilon || (std::isfinite(result.De) && result.De>0 &&
                    std::isfinite(result.Dex) && std::isfinite(result.Der))),
            "invalid coupled MMS scalar diffusivity/direction");
    return result;
}

double ConservativeDerivative(const MMSJet &q, const MMSFields &fields)
{
    // Retain the conservative identity explicitly even though this exact
    // streamfunction velocity has u_x+w_r=0 analytically.
    return q.dt+fields.u.value*q.dx+fields.w.value*q.dr+
           q.value*(fields.u.dx+fields.w.dr);
}
}

CoupledMMS::CoupledMMS(const Reference &reference, Geometry geometry,
                       double topRadius, double bodyAmplitude, MMSVariant variant)
    : m_reference(reference), m_geometry(geometry), m_top(topRadius),
      m_bodyAmplitude(bodyAmplitude), m_variant(variant)
{
    const auto &p=reference.Parameters();
    Require(variant==MMSVariant::Historical || variant==MMSVariant::PolynomialScalarsV1 ||
            variant==MMSVariant::NormalTubeV1,
            "unsupported continuous coupled MMS variant");
    Require(bool(geometry.normalMap)==(variant==MMSVariant::NormalTubeV1),
            "normal geometry requires the independent normal-tube MMS");
    if (geometry.normalMap)
        Require(geometry.normalMap->amplitude<=.3 && bodyAmplitude==0.,
                "normal MMS fixtures require a smooth small hill patch and zero extra body force");
    Require(std::isfinite(bodyAmplitude), "nonfinite prescribed MMS body amplitude");
    Require(std::isfinite(geometry.d0) && geometry.curve==0 &&
                std::isfinite(topRadius) && std::abs(topRadius-.5)<1e-13 &&
                std::abs(geometry.H-p.H)<1e-14 &&
                std::abs(geometry.h-p.roughness)<1e-14 &&
                geometry.H>0 && geometry.H<topRadius &&
                p.KIn>0 && p.nu>0 && p.Cd>0 && p.sigmaK>0 && p.sigmaE>0,
            "coupled MMS requires the fixed flat H/h model and R=.5 fixture");
    Require(.9604*p.KIn>=p.KMinimum && 1.0404*p.KIn<=p.KMaximum,
            "manufactured K range exceeds qualified reference interval");
}

const char *CoupledMMS::VariantName() const
{
    if (m_variant==MMSVariant::NormalTubeV1)
        return "NormalTubeV1";
    return m_variant==MMSVariant::Historical?"Historical":"PolynomialScalarsV1";
}

double CoupledMMS::Radius(double x, double z, double time) const
{
    CheckTangentialTime(x,time);
    const double r=z-m_geometry.d0;
    const double tolerance=64*std::numeric_limits<double>::epsilon()*
                           std::max(m_top,std::abs(m_geometry.d0));
    Require(std::isfinite(z) && r>=0 && r<=m_top+tolerance,
            "coupled MMS point is outside r in [0,.5]");
    return r;
}

MMSJet CoupledMMS::InterfaceK(double x, double time) const
{
    CheckTangentialTime(x,time);
    const double root=std::sqrt(m_reference.Parameters().KIn);
    const double a=root*(1+.012*std::sin(.6*time));
    const double b=root*.008*std::cos(.8*time);
    const double at=root*.0072*std::cos(.6*time);
    const double bt=-root*.0064*std::sin(.8*time);
    const double s=a+b*x, st=at+bt*x;
    Require(s>0, "manufactured square-root interface amplitude is nonpositive");
    MMSJet K;
    K.value=s*s; K.dx=2*s*b; K.dt=2*s*st; K.dxx=2*b*b;
    return K;
}

MMSFields CoupledMMS::Fields(double x, double z, double time) const
{
    const double r=Radius(x,z,time), R=m_top, H=m_geometry.H;
    if (m_variant==MMSVariant::NormalTubeV1)
        return NormalFields(x,r,time);
    MMSFields fields;
    const double S=.02*(1+.1*std::cos(.9*time));
    const double St=-.0018*std::sin(.9*time);
    const double V=.05*(1+.2*std::sin(.7*time));
    const double Vt=.007*std::cos(.7*time);
    const double X=x*x*(1-x), Xx=2*x-3*x*x, Xxx=2-6*x, Xxxx=-6;
    const double Y=r*r*(R-r), Yr=2*R*r-3*r*r, Yrr=2*R-6*r, Yrrr=-6;

    // psi=S(t) r^2/2+V(t) X(x)Y(r), u=psi_r, w=-psi_x.
    fields.u.value=S*r+V*X*Yr;
    fields.u.dx=V*Xx*Yr;
    fields.u.dr=S+V*X*Yrr;
    fields.u.dt=St*r+Vt*X*Yr;
    fields.u.dxx=V*Xxx*Yr;
    fields.u.dxr=V*Xx*Yrr;
    fields.u.drr=V*X*Yrrr;
    fields.w.value=-V*Xx*Y;
    fields.w.dx=-V*Xxx*Y;
    fields.w.dr=-V*Xx*Yr;
    fields.w.dt=-Vt*Xx*Y;
    fields.w.dxx=-V*Xxxx*Y;
    fields.w.dxr=-V*Xxx*Yr;
    fields.w.drr=-V*Xx*Yrr;

    const double pressure=.0002*(1+.3*std::sin(.8*time));
    const double pressureT=.000048*std::cos(.8*time);
    const double shape=x*x+r*r/(R*R)-2./3;
    fields.Pi.value=pressure*shape;
    fields.Pi.dx=2*pressure*x;
    fields.Pi.dr=2*pressure*r/(R*R);
    fields.Pi.dt=pressureT*shape;
    fields.Pi.dxx=2*pressure;
    fields.Pi.drr=2*pressure/(R*R);

    const MMSJet K=InterfaceK(x,time);
    const double eta=.005*(1+.2*std::sin(.7*time));
    const double etaT=.0007*std::cos(.7*time);
    double F,Fr,Frr;
    if (m_variant==MMSVariant::Historical)
    {
        const double denominator=H*H*R*R;
        F=r*r*(r-H)*(r-H)/denominator;
        Fr=(4*r*r*r-6*H*r*r+2*H*H*r)/denominator;
        Frr=(12*r*r-12*H*r+2*H*H)/denominator;
    }
    else
    {
        const double denominator=H*R;
        F=r*(r-H)/denominator;
        Fr=(2*r-H)/denominator;
        Frr=2/denominator;
    }
    const double kShape=1+eta*F;
    fields.k.value=K.value*kShape;
    fields.k.dx=K.dx*kShape;
    fields.k.dr=K.value*eta*Fr;
    fields.k.dt=K.dt*kShape+K.value*etaT*F;
    fields.k.dxx=K.dxx*kShape;
    fields.k.dxr=K.dx*eta*Fr;
    fields.k.drr=K.value*eta*Frr;

    fields.gradient[0][0]=fields.u.dx;
    fields.gradient[0][2]=fields.u.dr;
    fields.gradient[2][0]=fields.w.dx;
    fields.gradient[2][2]=fields.w.dr;
    fields.gradientX[0][0]=fields.u.dxx;
    fields.gradientX[0][2]=fields.u.dxr;
    fields.gradientX[2][0]=fields.w.dxx;
    fields.gradientX[2][2]=fields.w.dxr;
    fields.gradientR[0][0]=fields.u.dxr;
    fields.gradientR[0][2]=fields.u.drr;
    fields.gradientR[2][0]=fields.w.dxr;
    fields.gradientR[2][2]=fields.w.drr;
    return fields;
}

MMSJet CoupledMMS::OuterEpsilon(double x, double z, double time) const
{
    const double r=Radius(x,z,time), H=m_geometry.H, R=m_top;
    Require(OuterOwner(r,H), "transported MMS epsilon requested below H");
    if (m_variant==MMSVariant::NormalTubeV1)
        return NormalEpsilon(x,r,time);
    const MMSJet K=InterfaceK(x,time);
    const double A=m_reference.Amplitude(K.value);
    const double Ax=m_reference.AmplitudeDerivative(K.value)*K.dx;
    const double Axx=m_reference.AmplitudeDerivative(K.value)*K.dxx+
        m_reference.AmplitudeSecondDerivative(K.value)*K.dx*K.dx;
    const double At=m_reference.AmplitudeDerivative(K.value)*K.dt;
    const double zeta=.01*std::sin(.5*time), zetaT=.005*std::cos(.5*time);
    const double g=(r-H)*(r-H)/(R*R), gr=2*(r-H)/(R*R), grr=2/(R*R);
    const double shape=1+zeta*g;
    double f,fr,frr;
    if (m_variant==MMSVariant::Historical)
    {
        f=shape/r;
        fr=zeta*gr/r-shape/(r*r);
        frr=zeta*grr/r-2*zeta*gr/(r*r)+2*shape/(r*r*r);
    }
    else
    {
        f=shape/H;
        fr=zeta*gr/H;
        frr=zeta*grr/H;
    }
    MMSJet epsilon;
    epsilon.value=A*f;
    epsilon.dx=Ax*f;
    epsilon.dr=A*fr;
    if (m_variant==MMSVariant::Historical)epsilon.dt=At*f+A*zetaT*g/r;
    else epsilon.dt=At*f+A*zetaT*g/H;
    epsilon.dxx=Axx*f;
    epsilon.dxr=Ax*fr;
    epsilon.drr=A*frr;
    Require(epsilon.value>0, "nonpositive manufactured outer epsilon");
    return epsilon;
}

Vector CoupledMMS::Body(double x, double z, double time) const
{
    const double r=Radius(x,z,time);
    // Preserve the historical zero-body arithmetic exactly and avoid generating
    // signed zeros/trigonometric terms in the default fixture.
    if (m_bodyAmplitude==0) return {};
    return {m_bodyAmplitude*2e-4*(1+x+.5*r)*std::cos(.3*time),0.,
            m_bodyAmplitude*3e-4*(1+r)*std::sin(.2*time)};
}

Vector CoupledMMS::BodyAtH(double x, double time) const
{
    return Body(x,m_geometry.d0+m_geometry.H,time);
}

Vector CoupledMMS::AWFAtH(double x, double time) const
{
    if (m_variant==MMSVariant::NormalTubeV1)
        return NormalAWF(x,time);
    const auto atH=Fields(x,m_geometry.d0+m_geometry.H,time);
    const Vector physicalBody=BodyAtH(x,time);
    Vector awf{};
    awf[0]=2*atH.u.value*atH.u.dx+atH.Pi.dx-
           (2./3)*atH.k.dx-physicalBody[0];
    return awf;
}

MMSForcing CoupledMMS::Forcing(double x, double z, double time) const
{
    const double r=Radius(x,z,time);
    if (m_variant==MMSVariant::NormalTubeV1)
        return NormalForcing(x,r,time);
    const auto fields=Fields(x,z,time);
    const auto K=InterfaceK(x,time);
    MMSForcing force;
    force.body=Body(x,z,time);
    force.outer=OuterOwner(r,m_geometry.H);
    MMSJet epsilon;
    if (force.outer)
    {
        epsilon=OuterEpsilon(x,z,time);
    }
    const auto law=Coefficients(m_reference,fields,K,r,force.outer?&epsilon:nullptr);
    const Vector normal{0,0,1};
    const double nu=m_reference.Parameters().nu;
    const auto stress=AdaptAmended(fields.k.value,law.E,fields.gradient,normal,
                                   nu,law.auxiliary,r,m_geometry.H);
    const auto stressX=ApplyAmendedDirection(fields.k.value,law.E,fields.gradient,
        normal,nu,law.auxiliary,r,m_geometry.H,fields.k.dx,law.Ex,
        fields.gradientX,law.auxiliaryX,0.);
    const auto stressR=ApplyAmendedDirection(fields.k.value,law.E,fields.gradient,
        normal,nu,law.auxiliary,r,m_geometry.H,fields.k.dr,law.Er,
        fields.gradientR,law.auxiliaryR,1.);
    if (r<=m_geometry.h)
    {
        force.awf=AWFAtH(x,time);
    }
    force.momentum[0]=fields.u.dt+fields.u.value*fields.u.dx+
        fields.w.value*fields.u.dr-stressX.V[0][0]-stressR.V[0][2]+
        fields.Pi.dx-force.awf[0];
    force.momentum[2]=fields.w.dt+fields.u.value*fields.w.dx+
        fields.w.value*fields.w.dr-stressX.V[2][0]-stressR.V[2][2]+
        fields.Pi.dr-force.awf[2];
    if (m_bodyAmplitude!=0)
    {
        force.momentum[0]-=force.body[0];
        force.momentum[2]-=force.body[2];
    }
    force.production=stress.production;
    force.destruction=law.E;
    force.diffusionK=law.Dk;
    force.k=ConservativeDerivative(fields.k,fields)-
        law.Dk*(fields.k.dxx+fields.k.drr)-
        law.Dkx*fields.k.dx-law.Dkr*fields.k.dr-
        force.production+law.E;
    if (force.outer)
    {
        force.diffusionE=law.De;
        force.epsilon=ConservativeDerivative(epsilon,fields)-
            law.De*(epsilon.dxx+epsilon.drr)-
            law.Dex*epsilon.dx-law.Der*epsilon.dr-
            (Ce1*force.production-Ce2*epsilon.value)*epsilon.value/fields.k.value;
    }
    Require(std::isfinite(force.momentum[0]) && std::isfinite(force.momentum[2]) &&
                std::isfinite(force.k) && std::isfinite(force.epsilon),
            "nonfinite continuous coupled MMS forcing");
    return force;
}

double CoupledMMS::KOutwardFlux(double x, double z, double time,
                               const Vector &normal) const
{
    CheckNormal(normal);
    const double r=Radius(x,z,time);
    if (m_variant==MMSVariant::NormalTubeV1)
        return NormalFlux(x,r,time,normal,false);
    const auto fields=Fields(x,z,time);
    const auto K=InterfaceK(x,time);
    const bool outer=OuterOwner(r,m_geometry.H);
    MMSJet epsilon;
    if (outer)
    {
        epsilon=OuterEpsilon(x,z,time);
    }
    const auto law=Coefficients(m_reference,fields,K,r,outer?&epsilon:nullptr);
    return -law.Dk*(normal[0]*fields.k.dx+normal[2]*fields.k.dr);
}

double CoupledMMS::EpsilonOutwardFlux(double x, double z, double time,
                                     const Vector &normal) const
{
    CheckNormal(normal);
    const double r=Radius(x,z,time);
    if (m_variant==MMSVariant::NormalTubeV1)
        return NormalFlux(x,r,time,normal,true);
    const auto fields=Fields(x,z,time);
    const auto K=InterfaceK(x,time);
    const auto epsilon=OuterEpsilon(x,z,time);
    const auto law=Coefficients(m_reference,fields,K,r,&epsilon);
    return -law.De*(normal[0]*epsilon.dx+normal[2]*epsilon.dr);
}
}
