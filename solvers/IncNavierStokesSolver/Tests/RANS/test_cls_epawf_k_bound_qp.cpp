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
    Controls total;total.maximumTotal=1;
    Projection limited(Matrix{{1.,0.},{0.,1.}},total);
    Check(limited.Add({1.,0.},1.,0),"budget first solve");
    const auto saved=limited.delta;
    Check(!limited.Add({0.,1.},2.,1),"total work limit not enforced");
    Check(limited.iterations==1 && limited.totalLimitChecks==1 && limited.kktSolves==1,"completed/rejected total work counts");
    Check(limited.delta==saved && limited.solvedAdd==1 && !limited.workComplete,"failed result lost completed primal snapshot");
    Check(limited.solvedActive.size()==1 && limited.active.size()==2,"unfinished working set not separated");
    Check(limited.Stationarity()<1e-14 && limited.Complementarity()<1e-14 && limited.PoolViolation()>1.,"failed diagnostics do not describe last completed primal/dual");
    Check(std::abs(limited.Objective()-limited.DualObjective())<1e-14 && limited.DualViolation()==0,"coherent dual snapshot");
    Controls add;add.maximumPerAdd=1;add.maximumTotal=20;
    Projection perAdd(Matrix{{1.,0.},{0.,1.}},add);
    Check(perAdd.Add({1.,0.},1.,0),"per-add first solve");
    Check(!perAdd.Add({1.,1.},3.,1),"per-add limit not enforced");
    Check(perAdd.iterations==2 && perAdd.addLimitChecks==1 && perAdd.totalLimitChecks==0,"per-add/cumulative work confused");
    Check(perAdd.delta==Vector({1.,0.}) && perAdd.Stationarity()<1e-14,"per-add failure incoherent");
    Controls timed;timed.maximumSeconds=1e-15;
    Projection wall(Matrix{{1.}},timed);
    Check(!wall.Add({1.},1.,0) && wall.timeLimitChecks==1 && wall.iterations==0,"wall-time budget work accounting");
    Check(wall.delta==Vector({0.}) && wall.solvedActive.empty() && !wall.workComplete,"wall-time failed snapshot");
    Controls generous;generous.maximumPerAdd=20;generous.maximumTotal=100;
    Projection more(Matrix{{1.,0.},{0.,1.}},generous);
    Check(more.Add({1.,0.},1.,0) && more.Add({1.,1.},3.,1),"generous controls changed a finite solve");
    Check(std::abs(more.delta[0]-1.5)<1e-12 && more.maximumAddIterations==2 && more.iterations==3,"per-add completed work accounting");
    Check(more.Add({1.,1.},3.,2),"identical constraint regression");
    Check(more.repeatedRows==1 && more.repeatedConstraints==1,"duplicate telemetry");
    std::cout<<"Mass-metric positivity QP checks passed\n";
}
