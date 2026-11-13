// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallManufacturedJet.h>
#include <algorithm>
#include <iostream>
using namespace Nektar::CT6;
int main()
{
    double error=0.;
    const auto check=[&](double actual,double expected)
    {error=std::max(error,std::abs(actual-expected)/std::max(1.,std::abs(expected)));};
    const auto x=Jet::Variable(.3,0),y=Jet::Variable(-.2,1),t=Jet::Variable(.1,2);
    const auto polynomial=x*x*y+3*y*t*t+x*y*t;
    check(polynomial.D(0).D(1).Value(),.7);
    check(polynomial.D(1).D(2).Value(),.9);
    check(polynomial.D(0).D(0).D(1).Value(),2.);
    check(polynomial.D(1).D(2).D(2).Value(),6.);
    for(double sign:{-1.,1.})
    {
        const auto f=Jet(sign*2)+x+y*t;
        const auto unit=f*Inverse(f);
        check(unit.Value(),1.);
        for(int i=1;i<35;++i)check(unit.c[i],0.);
    }
    const auto exponential=Exp(x*y+t);
    check(exponential.D(0).D(1).Value(),std::exp(.04)*(.94));
    check(exponential.D(0).D(0).D(1).D(1).Value(),std::exp(.04)*(2- .24+.0036));
    const auto logarithm=Log(2+x);
    check(logarithm.D(0).D(0).Value(),-1/(2.3*2.3));
    const auto root=Sqrt(2+x);
    check(root.D(0).D(0).D(0).Value(),.375/std::pow(2.3,2.5));
    int rejected=0;
    try{(void)Inverse(Jet(0));}catch(const std::domain_error&){++rejected;}
    try{(void)Sqrt(Jet(-1));}catch(const std::domain_error&){++rejected;}
    try{(void)Log(Jet(-1));}catch(const std::domain_error&){++rejected;}
    std::cout<<"{\"scaled_error\":"<<error<<",\"expected_domain_rejections\":"<<rejected<<"}\n";
    return error<=2e-14 && rejected==3?0:1;
}
