// SPDX-License-Identifier: MIT
// Sampling adapter only. All expected values are computed in the independent
// Python oracle, which does not link to this implementation.
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <stdexcept>
using namespace Nektar::SolverUtils;
namespace
{
void Tensor(const RANSTensor &a)
{
    std::cout << '[';
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {if(i || j)std::cout << ',';std::cout << a[i][j];}
    std::cout << ']';
}
}
int main()
{
    std::cout << std::setprecision(17);
    std::string line;
    while(std::getline(std::cin,line))try
    {
        std::istringstream in(line);std::string kind,name,policy,wall;
        double nu;in>>kind>>name>>policy>>wall>>nu;
        const auto model=RANSResolvedModel::Resolve(ParseRANSModel(name),
            policy=="Natural"?RANSKappaPolicy::Natural:RANSKappaPolicy::Prescribed,
            .41,wall=="EP"?RANSWallFormulation::EPAWF:RANSWallFormulation::AcceptedBudgetAWF,nu);
        if(kind=="closure" || kind=="tangent")
        {
            double k,e;RANSTensor g{};in>>k>>e;
            for(auto &row:g)for(double &v:row)in>>v;
            if(!in)throw std::runtime_error("invalid closure probe input");
            const auto c=kind=="tangent"?EvaluateRANSClosureLinearization(model.Model(),k,e,g,1.,model.Parameters(),0.,true):EvaluateRANSClosure(model,k,e,g,1.);
            std::cout << "{\"ok\":true,\"q\":" << model.Q() << ",\"kappa\":" << model.Kappa()
                << ",\"sigma\":" << model.Parameters().sigmaEpsilon
                << ",\"cmu\":" << c.cMu << ",\"nut\":" << c.eddyViscosity
                << ",\"divergence\":" << c.divergence << ",\"production\":" << c.production
                << ",\"k_source\":" << c.kSource << ",\"eps_production\":" << c.epsilonProduction
                << ",\"eps_destruction\":" << c.epsilonDestruction << ",\"eps_source\":" << c.epsilonSource
                << ",\"yap\":" << c.yap << ",\"R\":";Tensor(c.reynoldsStress);
            std::cout << ",\"Q\":";Tensor(c.quadraticStress);
            std::cout << ",\"cubic\":";Tensor(c.cubicStress);
            std::cout << ",\"diffusion\":[" << c.momentumDiffusion << ',' << c.kDiffusion << ',' << c.epsilonDiffusion << "]}\n";
        }
        else if(kind=="wall")
        {
            double K,H,h,sample,tail,chi,ux,uy,fx,fy;
            int order,n;in>>K>>H>>h>>sample>>tail>>chi>>order>>ux>>uy>>fx>>fy>>n;
            AWFProfileOptions options;options.requestedWeight=chi;options.quadratureOrder=order;
            options.formulation=wall=="EP"?AWFProfileFormulation::LegacyEPAWF:AWFProfileFormulation::AcceptedBudgetAWF;
            AWFParameters wp=wall=="EP"?AWFParameters{}:CalibratedAWFParameters(model);
            EquilibriumAWFProfiles p({nu,H,h,sample,tail},K,model,{},options,wp);
            const auto &r=p.Reference();
            std::cout << "{\"ok\":true,\"q\":" << model.Q() << ",\"kappa\":" << model.Kappa()
                << ",\"J\":[" << p.Original().Integrals().j0 << ',' << p.Original().Integrals().j1 << ',' << p.Original().Integrals().j2
                << "],\"UH\":" << r.matchingVelocity << ",\"production_budget\":" << r.productionBudget
                << ",\"cutoff\":" << r.cutoff << ",\"budget_defect\":" << r.budgetResidual
                << ",\"chi\":" << r.weight << ",\"wall_coefficients\":[" << wp.cMu << ',' << wp.cLength << "]"
                << ",\"samples\":[";
            for(int i=0;i<n;++i)
            {
                double z;in>>z;if(!in)throw std::runtime_error("invalid wall sample input");
                auto v=p.Evaluate(z,{ux,uy},{fx,fy});
                if(i)std::cout << ',';
                const auto u=p.Original().Velocity(z,{ux,uy},{fx,fy});
                const auto g=p.Original().VelocityGradient(z,{ux,uy},{fx,fy});
                std::cout << '[' << z << ',' << v.production << ',' << v.referenceProduction << ',' << v.dissipation
                    << ',' << v.redistribution << ',' << v.kSource << ',' << p.ReferenceK(z) << ',' << p.ReferenceKFlux(z)
                    << ',' << p.RedistributionFlux(z) << ',' << u[0] << ',' << u[1] << ',' << g[0] << ',' << g[1] << ']';
            }
            std::cout << "]}\n";
        }
        else throw std::runtime_error("unknown probe operation");
    }
    catch(const std::exception &e){std::cout << "{\"ok\":false,\"error\":" << std::quoted(e.what()) << "}\n";}
}
