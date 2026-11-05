// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.h>
#include <cmath>
#include <cstring>
#include <iostream>
#include <stdexcept>
using Engine=Nektar::RANSWallTimeIntegration;
using Vector=Engine::Vector;
void Check(bool value,const char *message){if(!value)throw std::runtime_error(message);}
bool Bits(const Vector &a,const Vector &b)
{return a.size()==b.size() && std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0;}
void Same(const Engine::Snapshot &a,const Engine::Snapshot &b)
{
    Check(a.order==b.order && a.steps==b.steps && a.time==b.time && a.dt==b.dt,"time/step changed");
    Check(Bits(a.state,b.state) && Bits(a.times,b.times) && a.history.size()==b.history.size(),"state/time-label bits changed");
    for(std::size_t i=0;i<a.history.size();++i)Check(Bits(a.history[i],b.history[i]),"native history slot bits changed");
}
int main()
{
    for(int order:{1,2})
    {
        bool reject=false;int callbacks=0;
        auto implicit=[&](const Vector &h,double,double lambda)
        {
            ++callbacks;
            if(reject)throw std::runtime_error("deliberate callback rejection");
            return Vector{h[0]/(1+2*lambda),h[1]/(1+3*lambda)};
        };
        auto forcing=[](const Vector&,double time){return Vector{.7+time,-.3+2*time};};
        Engine control(order,.0025,implicit,forcing),trial(order,.0025,implicit,forcing);
        control.Initialise({1.,-.5},0.);trial.Initialise({1.,-.5},0.);
        control.Advance();trial.Advance();
        const auto before=trial.Capture();
        Check(std::abs(before.history[order][0])>1e-4,"fixture has no explicit history contribution");
        reject=true;bool failed=false;
        try{trial.Advance();}catch(const std::runtime_error&){failed=true;}
        Check(failed,"rejection did not occur");Same(before,trial.Capture());
        reject=false;const auto count=callbacks;trial.Advance();control.Advance();
        Check(callbacks==count+2,"not one callback per successful physical step");
        Same(control.Capture(),trial.Capture());
        // The actual IMEX1 implicit history contains dt times the explicit
        // body contribution. This generic fixture checks its leading update.
        if(order==1)
        {
            const double expected=(1.+.0025*.7)/(1.+.0025*2.);
            Check(std::abs(before.state[0]-expected)<1e-15,"explicit forcing omitted or double-counted");
        }
    }
    std::cout<<"IMEX1/2 bitwise callback rollback and nonzero explicit history passed\n";
}
