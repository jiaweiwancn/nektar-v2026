// SPDX-License-Identifier: MIT
#include <SolverUtils/Turbulence/RANSRobustVariant.h>
#include <SolverUtils/Turbulence/WJEARSM.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace Nektar::SolverUtils;
void Near(double x,double y,double tol=3e-12)
{if(std::abs(x-y)>tol*std::max({1.,std::abs(x),std::abs(y)}))throw std::runtime_error("returned stress/action identity");}
int main()
{
    const double k=1.7,e=.83;RANSTensor G{{{.3,1.1,-.4},{-.7,-.1,.6},{.9,-.3,.2}}};
    for(int mode=0;mode<3;++mode)
    {
        RANSRobustVariant v;v.beta=mode==0?1:.63;v.gamma=mode==2?.17:0;
        v.referenceStrain={{{.1,.2,0},{.2,-.2,.1},{0,.1,.1}}};
        auto c=EvaluateRANSRobustClosure(RANSModel::CLSHighRe,k,e,G,1.2,{},0,v);
        double work=0,devWork=0;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        {
            const double dev=c.reynoldsStress[i][j]-(i==j?2*k/3:0);
            Near(c.deviatoricStress[i][j],dev);Near(c.anisotropy[i][j],dev/k);
            Near(c.reynoldsStress[i][j],(i==j?2*k/3:0)-c.eddyViscosity*(G[i][j]+G[j][i])+c.stressRemainder[i][j]);
            work-=c.reynoldsStress[i][j]*G[i][j];devWork-=dev*G[i][j];
        }
        Near(c.production,work);Near(c.productionDev,devWork);
        Near(c.productionDev,c.production+2*k/3*(G[0][0]+G[1][1]+G[2][2]));
    }
    WJParameters p;p.cd=.123;p.sigmaK=1.7;p.sigmaEpsilon=1.9;
    auto law=RANSClosureKernel::WJ(p,WJGradientPolicy::DevStrainRawWorkV1);
    const RANSTensor H{{{.4,-.2,.3},{.1,.7,-.6},{.8,.5,-.3}}};
    auto c=law.Evaluate(k,e,G),d=law.Direction(k,e,G,.23,-.17,H);
    double theta=G[0][0]+G[1][1]+G[2][2];
    Near(c.productionDev,c.production+2*k*theta/3);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {
        const double devSym=G[i][j]+G[j][i]-(i==j?2*theta/3:0);
        Near(c.deviatoricStress[i][j],-c.eddyViscosity*devSym+c.stressRemainder[i][j]);
        // A native raw-strain implicit matrix needs the isotropic correction.
        const double rawRemainder=c.stressRemainder[i][j]+(i==j?2*c.eddyViscosity*theta/3:0);
        Near(c.deviatoricStress[i][j],-c.eddyViscosity*(G[i][j]+G[j][i])+rawRemainder);
    }
    double best=1;
    for(double h:{1e-3,1e-4,1e-5,1e-6})
    {
        auto gp=G,gm=G;for(int i=0;i<3;++i)for(int j=0;j<3;++j){gp[i][j]+=h*H[i][j];gm[i][j]-=h*H[i][j];}
        auto a=law.Evaluate(k+.23*h,e-.17*h,gp),b=law.Evaluate(k-.23*h,e+.17*h,gm);
        double error=0;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)error=std::max(error,std::abs((a.reynoldsStress[i][j]-b.reynoldsStress[i][j])/(2*h)-d.reynoldsStress[i][j]));
        for(auto triple:{std::array<double,3>{a.production,b.production,d.production},{a.kSource,b.kSource,d.kSource},
            {a.epsilonSource,b.epsilonSource,d.epsilonSource},{a.kDiffusion,b.kDiffusion,d.kDiffusion},
            {a.epsilonDiffusion,b.epsilonDiffusion,d.epsilonDiffusion}})
            error=std::max(error,std::abs((triple[0]-triple[1])/(2*h)-triple[2]));
        best=std::min(best,error);
    }
    Near(best,0,1e-7);
    std::cout<<"{\"returned_fields\":\"PASS\",\"nonzero_divergence_action_error\":"<<best<<"}\n";
}
