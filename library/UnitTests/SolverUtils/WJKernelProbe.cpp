///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/WJEARSM.h>
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <SolverUtils/Turbulence/RANSRobustVariant.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
using namespace Nektar::SolverUtils;
namespace
{
void Check(bool b) { if (!b) throw std::runtime_error("WJ compiled contract check failed"); }
template<class F> void Reject(F f)
{
    bool rejected=false;
    try { f(); } catch (const std::exception &) { rejected=true; }
    Check(rejected);
}
template<class A> void Array(const A &a)
{
    std::cout << '[';
    bool first=true;
    for (const auto &x:a) { if (!first) std::cout << ','; first=false; std::cout << x; }
    std::cout << ']';
}
void Tensor(const RANSTensor &t)
{
    std::cout << '[';
    for (int i=0;i<3;++i) { if (i) std::cout << ','; Array(t[i]); }
    std::cout << ']';
}
void Values(const RANSClosureValues &v)
{
    std::cout << "{\"scalar\":";
    Array(std::array<double,19>{v.cMu,v.eddyViscosity,v.scalarViscosity,v.production,
        v.productionDev,v.kSource,v.epsilonProduction,v.epsilonDestruction,v.epsilonSource,
        v.kDiffusion,v.epsilonDiffusion,v.divergence,v.wj.n,v.wj.s,v.wj.w,v.wj.i4,
        v.wj.i5,v.wj.q,v.wj.rootDerivative});
    std::cout << ",\"R\":"; Tensor(v.reynoldsStress);
    std::cout << ",\"a\":"; Tensor(v.anisotropy);
    std::cout << ",\"dev\":"; Tensor(v.deviatoricStress);
    std::cout << ",\"remainder\":"; Tensor(v.stressRemainder);
    std::cout << ",\"beta\":"; Array(v.wj.beta);
    std::cout << ",\"root_residual\":" << v.wj.rootResidual << '}';
}
void SelfTest()
{
    const auto ref=RANSNeutralReference::WJ();
    RANSTensor g{}; g[0][2]=ref.Q();
    auto kernel=RANSClosureKernel::WJ();
    auto c=kernel.Evaluate(1,1,g);
    auto common=EvaluateRANSClosure(RANSModel::WJHighRe,1,1,g,1);
    Check(std::abs(c.production-1)<1e-13 && c.wj.n>4);
    Check(common.epsilonDestruction==1.82 && common.reynoldsStress==c.reynoldsStress);
    Check(c.quadraticStress==RANSTensor{} && c.cubicStress==RANSTensor{} && c.yap==0);
    Check(RANSParameters{}.cEpsilon2==1.92);
    Reject([&]{ EvaluateRANSClosure(RANSModel::WJHighRe,1,1,g,1,RANSParameters{}); });
    for (bool branch:{false,true}) Reject([&]{ EvaluateRANSClosureLinearization(
        RANSModel::WJHighRe,1,1,g,1,{},0,branch); });
    Reject([&]{ EvaluateRANSRobustClosure(RANSModel::WJHighRe,1,1,g,1,{},0,{}); });
    Reject([&]{ kernel.Evaluate(1,1,g,1e-5); });
    Reject([&]{ kernel.Direction(1,1,g,0,0,{},1e-5); });
    for (double bad:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
    {
        Reject([&]{ kernel.Evaluate(bad,1,g); });
        Reject([&]{ kernel.Evaluate(1,bad,g); });
    }
    auto p=WJParameters::Defaults(); p.c1=1.9;
    Reject([&]{ RANSClosureKernel::WJ(p); });
    p=WJParameters::Defaults(); p.c2=.5;
    Reject([&]{ RANSClosureKernel::WJ(p); });
    p=WJParameters::Defaults(); p.cd=0;
    Reject([&]{ RANSClosureKernel::WJ(p); });
    Reject([&]{ RANSClosureKernel::WJ({},static_cast<WJGradientPolicy>(99)); });
    RANSTensor badg=g; badg[1][1]=1;
    Reject([&]{ kernel.Evaluate(1,1,badg); });
    Reject([&]{ kernel.Direction(1,1,g,0,0,badg); });
    badg[1][1]=std::numeric_limits<double>::quiet_NaN();
    Reject([&]{ kernel.Evaluate(1,1,badg); });
    Reject([&]{ kernel.Direction(1,1,g,0,0,badg); });
    Reject([&]{ kernel.Direction(1,1,g,std::numeric_limits<double>::infinity(),0,{}); });
    AWFProfileGeometry geom{1.5e-5,1,.1,1,1};
    AWFProfileOptions opt; opt.formulation=AWFProfileFormulation::WJCalibratedBudgetAWF;
    const double K=ref.Q()*.25;
    auto wall=EquilibriumAWFProfiles::WJ(geom,K,ref,opt);
    Check(std::abs(wall.ReferenceK(0)-1.8887084080086511)<1e-11);
    Check(std::abs(wall.ReferenceKFlux(1))<1e-12 && wall.ReferenceKFlux(0)==0);
    Check(wall.RedistributionFlux(.5)==0);
    Reject([&]{ wall.Reweighted(.5); });
    opt.requestedWeight=0;
    Reject([&]{ EquilibriumAWFProfiles::WJ(geom,K,ref,opt); });
    opt.requestedWeight=1; opt.failurePolicy=AWFFailurePolicy::UseOriginal;
    Reject([&]{ EquilibriumAWFProfiles::WJ(geom,K,ref,opt); });
    opt={}; opt.formulation=AWFProfileFormulation::WJCalibratedEPAWF;
    Reject([&]{ EquilibriumAWFProfiles(geom,K,ref,{},opt); }); // no silent legacy coefficients
    Reject([&]{ EquilibriumAWFProfiles::WJ(geom,K,RANSNeutralReference::CLS(),opt); });
    geom.molecularNu=1; geom.roughness=0;
    Reject([&]{ EquilibriumAWFProfiles::WJ(geom,K,ref,opt); }); // gate failure
    opt.failurePolicy=AWFFailurePolicy::UseOriginal;
    Check(EquilibriumAWFProfiles::WJ(geom,K,ref,opt).Reference().weight==0);
    opt.failurePolicy=AWFFailurePolicy::RequirePreservation; opt.requestedWeight=0;
    auto inactive=EquilibriumAWFProfiles::WJ(geom,K,ref,opt);
    Check(inactive.Reference().cutoffIterations==0 && inactive.Reference().cutoff==0);
    Check(inactive.Evaluate(.5,{0,0},{0,0}).dissipation==inactive.Original().OriginalDissipation(.5));
    for (double bad:{0.,-1.,std::numeric_limits<double>::quiet_NaN()})
    {
        Reject([&]{ EquilibriumAWFProfiles::WJ(geom,bad,ref,opt); });
        auto invalid=geom; invalid.molecularNu=bad;
        Reject([&]{ EquilibriumAWFProfiles::WJ(invalid,K,ref,opt); });
    }
    Reject([&]{ EquilibriumAWFProfiles::WJ(geom,K,ref,opt,{-2,0,.1}); });
    std::cout << "{\"compiled_contract\":\"PASS\"}\n";
}
void Wall(std::istream &in)
{
    int formulation, fallback; double chi,K,nu,H,h,zm,k1,k2;
    AWFVector u{},f{};
    in>>formulation>>chi>>fallback>>K>>nu>>H>>h>>zm>>k1>>k2>>u[0]>>u[1]>>f[0]>>f[1];
    AWFProfileOptions options;
    options.formulation=formulation==0 ? AWFProfileFormulation::WJCalibratedEPAWF : AWFProfileFormulation::WJCalibratedBudgetAWF;
    options.requestedWeight=chi; options.quadratureOrder=24;
    options.failurePolicy=fallback ? AWFFailurePolicy::UseOriginal : AWFFailurePolicy::RequirePreservation;
    auto wall=EquilibriumAWFProfiles::WJ({nu,H,h,H,zm},K,RANSNeutralReference::WJ(),options,{k1,k2,.1});
    const auto &r=wall.Reference(); const auto &base=wall.Original();
    const auto actual=base.Evaluate(u,f);
    std::cout<<"{\"reference\":";
    Array(std::array<double,15>{r.traction,r.matchingVelocity,r.numerator,r.productionBudget,
        r.minimumBudget,r.cutoff,r.weight,base.ViscosityInterceptHeight(),base.DissipationCutoff(),
        base.Parameters().cMu,base.Parameters().cLength,base.Integrals().j0,base.Integrals().j1,
        base.Integrals().j2,actual.meanProduction});
    std::cout<<",\"traction\":"; Array(actual.traction);
    std::cout<<",\"interface_traction\":";
    Array(AWFVector{actual.traction[0]+f[0]*(H-h),actual.traction[1]+f[1]*(H-h)});
    std::vector<std::vector<double>> basis(6);
    std::cout<<",\"quadrature\":[";
    bool first=true;
    for (const auto &q:wall.Quadrature())
    {
        if (!first) std::cout<<',';
        first=false;
        const auto v=wall.Evaluate(q.r,u,f);
        const auto velocity=base.Velocity(q.r,u,f), shear=base.VelocityGradient(q.r,u,f);
        Array(std::array<double,15>{q.r,q.weight,q.jacobian,v.production,v.referenceProduction,
            v.dissipation,v.redistribution,v.kSource,wall.ReferenceKFlux(q.r),
            wall.ReferenceKGradient(q.r),base.TotalViscosity(q.r),velocity[0],velocity[1],shear[0],shear[1]});
        for (int n=0;n<6;++n) basis[n].push_back(std::pow(q.r/H,n));
    }
    std::cout<<"],\"loads\":"; Array(wall.SourceLoads(basis,u,f));
    std::cout<<",\"profile\":[";
    first=true;
    for (double z:{0.,.1*H,.5*H,H})
    {
        if (!first) std::cout<<',';
        first=false;
        const double kw=wall.ReferenceK(z), gz=9.81*z;
        // p_o(H)/rho=0 gauge; normal balance uses the actual local kw.
        const double pressure=-K/(81./20)+9.81*(H-z)+2*(K-kw)/3;
        Array(std::array<double,6>{z,kw,wall.ReferenceKFlux(z),pressure,
            pressure+gz+2*kw/3,wall.Dissipation(z)});
    }
    std::cout<<"]}";
}
}
int main(int argc,char **argv)
{
    std::cout<<std::setprecision(17);
    if (argc==1) { try { SelfTest(); } catch (const std::exception &e) { std::cerr<<e.what()<<'\n';return 1;} return 0; }
    if (std::string(argv[1])!="--stream") return 2;
    std::string line;
    while (std::getline(std::cin,line))
    {
        try
        {
            std::istringstream in(line); char command; in>>command;
            if (command=='C' || command=='D')
            {
                int completion; double k,e,nu; RANSTensor g{},dg{};
                in>>completion>>k>>e>>nu;
                for (auto &row:g) for (auto &x:row) in>>x;
                auto kernel=RANSClosureKernel::WJ({},completion ? WJGradientPolicy::DevStrainRawWorkV1 : WJGradientPolicy::Incompressible);
                if (command=='C') Values(kernel.Evaluate(k,e,g,nu));
                else { double dk,de; in>>dk>>de; for(auto &row:dg)for(auto &x:row)in>>x; Values(kernel.Direction(k,e,g,dk,de,dg,nu)); }
            }
            else if (command=='W') Wall(in);
            else if (command=='N')
            {
                double cd,z,u,H,U; in>>cd>>z>>u>>H>>U;
                WJParameters p; p.cd=cd;
                auto ref=RANSNeutralReference::WJ(p); auto v=ref.AtHeight(z,u,H,U);
                std::cout<<"{\"neutral\":";
                Array(std::array<double,11>{ref.Q(),ref.Kappa(),v.velocity,v.k,v.epsilon,
                    v.eddyViscosity,v.velocityGradient,v.epsilonGradient,v.epsilonGradientFlux,v.shearStress,cd});
                std::cout<<",\"R\":";Tensor(v.covariance); std::cout<<'}';
            }
            else throw std::invalid_argument("Unknown command");
            std::cout<<'\n';
        }
        catch (const std::exception &) { std::cout<<"{\"rejected\":true}\n"; }
    }
}
