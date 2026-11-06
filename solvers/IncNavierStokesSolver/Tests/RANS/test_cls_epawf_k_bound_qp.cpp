// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallKBoundProjection.h>
#include <iostream>
using namespace Nektar::RANSWallKBound;
void Check(bool x,const char *s){if(!x)throw std::runtime_error(s);}
int main()
{
    Projection p({{2.,.5},{.5,1.}});
    Check(Norm(p.delta)==0,"initial no-op");
    Check(p.Add({1.,0.},1.,0),"first cut");
    Check(std::abs(p.delta[0]-1.)<1e-13 && std::abs(p.delta[1]+.5)<1e-13,"consistent mass minimizer");
    Check(p.Add({0.,1.},.25,1),"second cut");
    Check(std::abs(p.delta[0]-1.)<1e-12 && std::abs(p.delta[1]-.25)<1e-12,"two active constraints");
    Check(p.Add({1.,0.},2.,2),"stronger dependent cut");
    Check(std::abs(p.delta[0]-2.)<1e-12,"dependent replacement");
    Check(p.Stationarity()<1e-12 && p.Complementarity()<1e-12,"KKT accounting");
    Projection infeasible(Matrix{{1.}});Check(infeasible.Add({1.},1.,0),"first feasible half-space");
    Check(!infeasible.Add({-1.},0.,1),"incompatible half-spaces accepted");
    Projection zero(Matrix{{1.}});Check(!zero.Add({0.},1.,0),"prescribed incompatibility accepted");
    Projection leaving(Matrix{{1.,0.},{0.,1.}});
    Check(leaving.Add({1.,0.},1.,0),"leaving first cut");
    Check(leaving.Add({1.,1.},3.,1),"leaving second cut");
    Check(std::abs(leaving.delta[0]-1.5)<1e-12 && std::abs(leaving.delta[1]-1.5)<1e-12,"leaving minimizer");
    Check(std::find(leaving.active.begin(),leaving.active.end(),0)==leaving.active.end(),"first cut did not leave");
    Check(leaving.Add({0.,1.},3.,2),"leaving third cut");
    const Matrix inventory{{1.,0.},{1.,1.},{0.,1.}};const Vector lower{1.,3.,3.};
    int rediscovered=0;
    for(int separation=0;separation<5;++separation)
    {
        int worst=-1;double violation=1e-12;
        for(int i=0;i<3;++i)if(lower[i]-Dot(inventory[i],leaving.delta)>violation)
        {worst=i;violation=lower[i]-Dot(inventory[i],leaving.delta);}
        if(worst<0)break;
        if(worst==0)++rediscovered;
        Check(leaving.Add(inventory[worst],lower[worst],worst),"rediscovered cut failed");
    }
    Check(rediscovered==1 && std::abs(leaving.delta[0]-1.)<1e-12 && std::abs(leaving.delta[1]-3.)<1e-12,"full separation did not recover (1,3)");
    for(int i=0;i<3;++i)Check(Dot(inventory[i],leaving.delta)>=lower[i]-1e-12,"inactive inventory violation");
    Check(leaving.Stationarity()<1e-12 && leaving.Complementarity()<1e-12,"rediscovered KKT diagnostics");
    std::cout<<"Mass-metric positivity QP checks passed\n";
}
