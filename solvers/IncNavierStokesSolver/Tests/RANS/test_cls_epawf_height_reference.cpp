// SPDX-License-Identifier: MIT
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <iostream>
#include <iomanip>
#include <cmath>
using namespace Nektar::SolverUtils;
int main()
{
    const auto model=RANSNeutralReference::CLS();const double us=.4,nu=1.5e-5,h=.1,K=model.Q()*us*us;
    const EquilibriumAWFProfiles baseline({nu,2,h,1,1},K,model);
    const double zstar=2*std::exp(-model.Kappa()*baseline.Reference().matchingVelocity/us);
    std::cout<<std::setprecision(17)<<"{\"schema\":1,\"z0_baseline\":"<<zstar<<",\"rows\":[";bool first=true;
    for(double H:{2.,1.,.5})for(double fraction:{.5,1.})
    {
        AWFProfileOptions options;options.quadratureOrder=64;
        EquilibriumAWFProfiles p({nu,H,h,H/2,H*fraction},K,model,{},options);
        const auto &r=p.Reference();const auto &a=p.Original();
        const double scale=nu/(.09*2.55*std::sqrt(K))-a.ViscosityInterceptHeight();
        const double UHstar=us/model.Kappa()*std::log(H/zstar),tstar=UHstar/a.Integrals().j0;
        if(!first)std::cout<<',';first=false;
        std::cout<<"{\"H\":"<<H<<",\"h\":"<<h<<",\"Z\":4,\"z_m\":"<<H*fraction
          <<",\"H_over_h\":"<<H/h<<",\"H_over_z0\":"<<std::exp(model.Kappa()*r.matchingVelocity/us)
          <<",\"H_over_Z\":"<<H/4<<",\"velocity_scale\":"<<scale<<",\"velocity_scale_over_H\":"<<scale/H
          <<",\"Re_H\":"<<r.reynolds<<",\"U_H\":"<<r.matchingVelocity<<",\"traction\":"<<r.traction
          <<",\"z0_c\":"<<H*std::exp(-model.Kappa()*r.matchingVelocity/us)<<",\"cutoff\":"<<r.cutoff
          <<",\"cutoff_over_H\":"<<r.cutoff/H<<",\"Lambda\":"<<r.productionBudget/r.numerator
          <<",\"budget_ratio\":"<<r.productionBudget/r.minimumBudget<<",\"intrinsic_chi\":"<<r.weight
          <<",\"J0\":"<<a.Integrals().j0<<",\"K0\":"<<K<<",\"kappa\":"<<model.Kappa()
          <<",\"Zv\":"<<a.ViscosityInterceptHeight()<<",\"fixed_z0_U_H\":"<<UHstar<<",\"fixed_z0_traction\":"<<tstar
          <<",\"fixed_z0_traction_minus_ustar2\":"<<tstar-us*us<<"}";
    }
    std::cout<<"]}\n";
}
