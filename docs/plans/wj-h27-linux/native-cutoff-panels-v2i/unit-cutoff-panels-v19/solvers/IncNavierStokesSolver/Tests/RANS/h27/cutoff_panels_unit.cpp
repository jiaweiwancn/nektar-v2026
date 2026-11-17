// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeCutoffPanels.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
using Nektar::WJNative::CubicCutoffPanels;
void Require(bool ok,const char*message){if(!ok)throw std::runtime_error(message);}
int main(int argc,char**argv)
{
    try{
        Require(argc==2,"usage: cutoff_panels_unit result.json");
        std::ifstream old(argv[1]);Require(!old.good(),"unit output is immutable");
        double maximum=0;int cases=0;
        auto roots=[&](std::array<double,4>c,double level,std::vector<double>expected){
            const auto panels=CubicCutoffPanels(c,{level});++cases;
            for(double x:panels)Require(std::isfinite(x)&&x>0&&x<1,"invalid interior partition");
            for(double x:expected){double error=1;for(double y:panels)error=std::min(error,std::abs(x-y));maximum=std::max(maximum,error);Require(error<=2e-14,"expected cubic level root missing");}
            return panels;
        };
        roots({0,1,0,0},.37,{.37});
        roots({-.08,.08,0,-.25},0,{.2,.5,.8});
        roots({.25,.25,-1,0},0,{.5});
        const double near=.25-1e-12;const double width=double(std::sqrt(.25L-(long double)near));
        roots({near,near,-1,0},0,{.5-width,.5+width});
        Require(roots({1,1,0,0},1,{}).empty(),"constant on-face trace creates artificial roots");
        Require(roots({1,1,0,0},2,{}).empty(),"constant no-crossing trace creates roots");
        Require(roots({0,1,0,0},2,{}).empty(),"monotone out-of-range level creates roots");
        Require(CubicCutoffPanels({.25,.25,-1,0},{}).empty(),"empty face set creates partitions");++cases;
        roots({0,1,0,1e-30},.31,{.31});
        bool badTrace=false,badLevel=false;
        try{CubicCutoffPanels({0,1,std::numeric_limits<double>::quiet_NaN(),0},{.5});}catch(const std::invalid_argument&){badTrace=true;}
        try{CubicCutoffPanels({0,1,0,0},{std::numeric_limits<double>::infinity()});}catch(const std::invalid_argument&){badLevel=true;}
        Require(badTrace&&badLevel,"nonfinite cutoff input accepted");
        std::ofstream out(argv[1]);out<<std::setprecision(17)<<"{\"schema\":\"WJH27CubicCutoffPanelsUnitV1\",\"status\":\"PASS_CUBIC_PANEL_UNIT_ONLY\",\"cases\":"<<cases
            <<",\"invalid_controls\":2,\"maximum_expected_root_error\":"<<maximum<<",\"full_H27_qualification\":false,\"TASK3_READY\":false}\n";
        out.close();Require(bool(out),"unit report write failed");return 0;
    }catch(const std::exception&e){std::cerr<<"Cubic cutoff panel unit: "<<e.what()<<'\n';return 1;}
}
