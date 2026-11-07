// SPDX-License-Identifier: MIT
#include <SolverUtils/Turbulence/EquilibriumAWFTransition.h>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace Nektar::SolverUtils;
int main()
{
    const auto model=RANSNeutralReference::CLS();
    for(double H:{.5,2.})
    {
        AWFProfileOptions options;options.quadratureOrder=16;
        options.physicalBreakpoints={0,.01,.1,.25,H};
        const double K=model.Q()*.4*.4;
        const EquilibriumAWFProfiles active({1.5e-5,H,.1,H/2,H},K,model,{},options);
        options.requestedWeight=0;
        const EquilibriumAWFProfiles inactive({1.5e-5,H,.1,H/2,H},K,model,{},options);
        bool distinct=active.Quadrature().size()!=inactive.Quadrature().size();
        for(double chi:{1e-12,.125,.5,.9,1.})
        {
            const auto profile=active.Reweighted(chi);
            if(profile.Quadrature().size()!=active.Quadrature().size())throw std::runtime_error("positive chi changed size");
            for(int i=0;i<profile.Quadrature().size();++i)
            {
                const auto &a=profile.Quadrature()[i],&b=active.Quadrature()[i];
                if(a.r!=b.r || a.weight!=b.weight || a.jacobian!=b.jacobian)
                    throw std::runtime_error("positive chi changed a physical quadrature row");
                if(i<inactive.Quadrature().size() && a.r!=inactive.Quadrature()[i].r)distinct=true;
            }
        }
        if(!distinct)throw std::runtime_error("fixture failed to distinguish inactive/active sets");
        AWFTransitionControls c;c.start=-3;c.end=-2.25;
        AWFTransitionInput input;input.coordinate=-3.5;input.budgetMargin=active.Reference().margin;input.reynolds=active.Reference().reynolds;
        const auto on=EvaluateEquilibriumAWFTransition(input,c).value;
        input.sensors.velocity[0]=1;
        const auto off=EvaluateEquilibriumAWFTransition(input,c).value;
        if(!(on.chi>0 && off.chi==0))throw std::runtime_error("state-only gate switch fixture failed");
    }
    std::cout<<"positive chi shares Q+; state gate can switch without changing K; tail=H\n";
}
