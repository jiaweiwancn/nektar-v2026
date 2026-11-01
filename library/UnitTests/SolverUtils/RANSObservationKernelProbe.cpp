///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSWallGeometry.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
using namespace Nektar::SolverUtils;
namespace
{
int checks = 0;
void Check(bool condition) { ++checks;if (!condition) throw std::runtime_error("Observation kernel gate failed"); }
void Near(double a,double b,double tolerance=5e-10)
{ Check(std::isfinite(a) && std::isfinite(b) && std::abs(a-b) <= tolerance*std::max({1.,std::abs(a),std::abs(b)})); }
double Dot(const RANSVector &a,const RANSVector &b)
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
}
int main(int argc,char **argv)
{
    try
    {
        std::ofstream out(argc > 1 ? argv[1] : "rans-observation-kernels.json");
        out << std::setprecision(17) << "{\"reference_cases\":[";
        bool first = true;double maximumShearError=0.;
        for (bool cls : {false,true})
        {
            const auto reference = cls ? RANSNeutralReference::CLS() : RANSNeutralReference::Standard();
            const double star=.3,q=reference.Q(),kap=reference.Kappa(),K=q*star*star,A=star*star*star/kap;
            const auto &p=reference.Parameters();
            const double Y0=cls ? CLSYapFunction(kap*std::pow(q,1.5)/p.yapLength,p) : 0.;
            for (double nu : {0.,1.5e-5}) for (double displacement : {0.,.0035})
                for (double eta : {.008,.016,.04,.25})
                {
                    auto closure = [&](double z) {
                        RANSTensor g{};g[0][2]=star/(kap*(z-displacement));
                        return EvaluateRANSClosure(reference.Model(),K,A/(z-displacement),g,z,p,nu);
                    };
                    const auto value=closure(eta);
                    const double h=1e-5*(eta-displacement),r=eta-displacement;
                    const auto low=closure(eta-h),high=closure(eta+h);
                    const double momentum=((nu*star/(kap*(r+h))-high.reynoldsStress[0][2])-
                                           (nu*star/(kap*(r-h))-low.reynoldsStress[0][2]))/(2*h);
                    const double epsilon=(-(nu+high.eddyViscosity/p.sigmaEpsilon)*A/std::pow(r+h,2)+
                                            (nu+low.eddyViscosity/p.sigmaEpsilon)*A/std::pow(r-h,2))/(2*h)+value.epsilonSource;
                    const double Y=cls ? CLSYapFunction(kap*std::pow(q,1.5)*r/(p.yapLength*eta),p) : 0.;
                    const double expectedM=-nu*star/(kap*r*r);
                    const double expectedE=2*nu*A/std::pow(r,3)+(Y-Y0)*A*A/(K*r*r);
                    const double scaleM=star*star/r,scaleE=A*A/(K*r*r);
                    Near(momentum/scaleM,expectedM/scaleM,2e-8);
                    Near(epsilon/scaleE,expectedE/scaleE,2e-8);
                    Near(value.kSource,0);
                    if (nu==0 && (!cls || displacement==0)) Near(epsilon/scaleE,0,2e-8);
                    if (nu>0) Check(std::abs(momentum/scaleM)>1e-6);
                    if (cls && displacement>0 && nu==0) Check(std::abs(epsilon/scaleE)>1e-5);
                    if (!first) out << ',';
                    first=false;
                    out << "{\"cls\":" << cls << ",\"nu\":" << nu << ",\"d\":" << displacement
                        << ",\"eta\":" << eta << ",\"q\":" << q << ",\"kappa\":" << kap
                        << ",\"momentum_scaled\":" << momentum/scaleM << ",\"epsilon_scaled\":" << epsilon/scaleE
                        << ",\"expected_momentum_scaled\":" << expectedM/scaleM
                        << ",\"expected_epsilon_scaled\":" << expectedE/scaleE << '}';
                }
            for (double roughness : {0.,.002,.005}) for (double curvature : {-20.,0.,20.})
            {
                EquilibriumWallFunction corrected({1.5e-5,.008,roughness,.004},K,reference);
                const auto &wall=corrected.Original();
                RANSWallGeometry plane;plane.terrain=RANSTerrain::Plane;plane.slope=.4;
                auto frame=plane.At(.1);frame.curvature=curvature;
                for (const AWFVector force : {AWFVector{0.,0.},AWFVector{.03,-.01}})
                {
                    const AWFVector matching{corrected.Reference().matchingVelocity,.2};
                    const double a=.7,b=-.2,H=wall.Geometry().height;
                    for (double eta : {0.,.001,.003,.004,H})
                    {
                        const auto observed=ReconstructAWFVelocity(wall,frame,eta,matching,force,a,b);
                        const auto local=wall.Velocity(eta,matching,force);
                        Near(Dot(observed,frame.tangent),local[0]);Near(observed[1],local[1]);
                        Near(frame.Jacobian(eta)*Dot(observed,frame.normal)+eta*a+frame.Volume(eta)*b,0);
                        if (eta==H) Near(Dot(observed,frame.tangent),matching[0]);
                        if (eta==0) for (double v : observed) Near(v,0);
                    }
                    const auto a0=wall.Evaluate(matching,force);
                    const double step=1e-3*H;
                    const auto u0=wall.Velocity(H,matching,force),u1=wall.Velocity(H-step,matching,force),
                               u2=wall.Velocity(H-2*step,matching,force),u3=wall.Velocity(H-3*step,matching,force),
                               u4=wall.Velocity(H-4*step,matching,force);
                    for (int i=0;i<2;++i)
                    {
                        const double flux=wall.TotalViscosity(H)*(25*u0[i]-48*u1[i]+36*u2[i]-16*u3[i]+3*u4[i])/(12*step);
                        const double expected=a0.traction[i]+force[i]*std::max(0.,H-roughness);
                        maximumShearError=std::max(maximumShearError,std::abs(flux-expected)/std::max(1.,std::abs(expected)));
                        Near(flux,expected);
                    }
                    if (force[0]==0.) Near(a0.traction[0],star*star);
                    const auto a1=corrected.EvaluateFixedWeight(matching,force,1.);
                    Near(a0.traction[0],a1.values.traction[0],0);
                    Near(a0.meanProduction,a1.values.meanProduction,0);
                }
            }
        }
        bool rejected=false;
        try
        {
            RANSWallGeometry g;AnalyticalWallFunction wall({1e-5,.01,0,.005},.1);
            ReconstructAWFVelocity(wall,g.At(0),.02,{1.,0.},{0.,0.},0.,0.);
        }
        catch (const std::exception &) { rejected=true; }
        Check(rejected);
        out << "],\"checks\":" << checks << ",\"auxiliary_H_shear_flux_max_scaled_error\":"
            << maximumShearError << ",\"status\":\"PASS\"}\n";
        Check(out.good());std::cout << checks << " observation/reference checks passed\n";
    }
    catch (const std::exception &e) { std::cerr << e.what() << '\n';return 1; }
}
