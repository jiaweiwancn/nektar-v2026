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
    Projection infeasible({{1.}});Check(infeasible.Add({1.},1.,0),"first feasible half-space");
    Check(!infeasible.Add({-1.},0.,1),"incompatible half-spaces accepted");
    Projection zero({{1.}});Check(!zero.Add({0.},1.,0),"prescribed incompatibility accepted");
    std::cout<<"Mass-metric positivity QP checks passed\n";
}
