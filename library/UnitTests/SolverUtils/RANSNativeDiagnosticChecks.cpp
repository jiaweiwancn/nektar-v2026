// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallForwardDiagnostics.h>
#include <iostream>
#include <iomanip>
using namespace Nektar;
using namespace Nektar::SolverUtils;
namespace J=Nektar::RANSWallSolveJson;
int main()
{
    try
    {
        int checks=0;double worst=0;
        auto require=[&](bool ok){++checks;if(!ok)throw std::runtime_error("native diagnostic regression");};
        const double c=std::cos(.37),s=std::sin(.37);
        const RANSTensor rotation{{{c,0,s},{0,1,0},{-s,0,c}}};
        for(const auto spectrum: {RANSVector{.2,1.,2.},RANSVector{0.,1.,2.},RANSVector{-1./3,2./3,5./3}})
        for(double scale:{1e-280,1e-16,1.,1e16,1e280})
        for(bool rotate:{false,true})
        {
            RANSTensor a{};
            if(rotate)for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)
                a[i][j]+=scale*rotation[i][k]*spectrum[k]*rotation[j][k];
            else for(int i=0;i<3;++i)a[i][i]=scale*spectrum[i];
            const auto full=J::Eigensystem(a);const auto eigen=J::Eigenvalues(a);
            require(eigen==full.eigenvalues);
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            {
                double product=0,normal=0;for(int k=0;k<3;++k){product+=a[i][k]/scale*full.vectors[k][j];normal+=full.vectors[k][i]*full.vectors[k][j];}
                require(std::abs(product-eigen[j]/scale*full.vectors[i][j])<3e-13);
                require(std::abs(normal-(i==j?1.:0.))<3e-13);
            }
            for(int i=0;i<3;++i){const double error=std::abs(eigen[i]/scale-spectrum[i]);worst=std::max(worst,error);require(error<2e-13);}
            require(J::CovarianceAdmissible(eigen[0],scale)==(spectrum[0]>=0));
        }
        const RANSTensor zero{};require(J::Eigenvalues(zero)==RANSVector{0,0,0});
        RANSTensor tiny{{{2e-16/3,0,-1e-16},{0,2e-16/3,0},{-1e-16,0,2e-16/3}}};
        const double minimum=J::Eigenvalues(tiny)[0];require(std::abs(minimum/1e-16+1./3)<1e-14);
        require(!J::CovarianceAdmissible(minimum,1e-16));
        const auto model=RANSResolvedModel::Resolve(RANSModel::A95,RANSKappaPolicy::Natural,.41,RANSWallFormulation::EPAWF,1.5e-5);
        const double a=6.5/(std::sqrt(18.)-3)*(1-1e-14);
        RANSTensor gradient{};gradient[0][0]=-a;gradient[2][2]=a;
        const auto physical=RANSWallDiagnostics::Observe(model,1.,1.,gradient,1.);
        require(physical.closure.a95RadicandDefined && physical.closure.a95Radicand>0);
        require(physical.covarianceEigenvalues[0]>0 && !physical.tangentAvailable);
        int invalidProbes=0;for(double sign:{-1.,1.})
        {
            auto perturbed=gradient;perturbed[1][0]+=sign*1e-5;
            try{(void)EvaluateRANSClosure(model,1.,1.,perturbed,1.);}
            catch(const std::domain_error &){++invalidProbes;}
        }
        require(invalidProbes==2);
        std::cout<<std::setprecision(17)<<J::Object({{"status",J::Quote("pass")},{"checks",J::Number(checks)},
            {"maximum_scaled_eigenvalue_error",J::Number(worst)},{"tiny_tensor_minimum",J::Number(minimum)},
            {"A95_actual_radicand",J::Number(physical.closure.a95Radicand)},
            {"A95_covariance_eigenvalues",J::Numbers(physical.covarianceEigenvalues)},
            {"tangent_status",J::Quote(physical.tangentStatus)},{"old_stencil_domain_failures",J::Number(invalidProbes)}})<<std::endl;
        return 0;
    }
    catch(const std::exception &error){std::cerr<<error.what()<<std::endl;return 1;}
}
