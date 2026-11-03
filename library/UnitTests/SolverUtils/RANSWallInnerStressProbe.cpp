///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSWallInnerStress.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Nektar::SolverUtils;
namespace
{
int checks=0;
void Check(bool condition,const std::string &name)
{++checks;if(!condition)throw std::runtime_error("Inner-stress probe: "+name);}
void Tensor(std::ostream &out,const RANSTensor &a)
{out<<'[';for(int i=0;i<3;++i){if(i)out<<',';out<<'[';for(int j=0;j<3;++j){if(j)out<<',';out<<a[i][j];}out<<']';}out<<']';}
void Values(std::ostream &out,const RANSWallInnerStressValues &v,double nu,double pi)
{
    out<<"{\"strain\":";Tensor(out,v.strain);out<<",\"turbulent\":";Tensor(out,v.turbulentStress);
    out<<",\"stress\":";Tensor(out,v.reynoldsStress);out<<",\"mean_stress\":";Tensor(out,RANSWallInnerMeanStress(v,nu,pi));
    out<<",\"s\":"<<v.strainNorm<<",\"a\":"<<v.a<<",\"c\":"<<v.c<<",\"t\":"<<v.t
       <<",\"ta\":"<<v.ta<<",\"tc\":"<<v.tc<<",\"radial\":"<<v.radialTangent
       <<",\"tangential\":"<<v.tangentialTangent<<",\"production\":"<<v.production
       <<",\"active\":"<<(v.active?"true":"false")<<",\"linear_response\":"<<(v.linearResponse?"true":"false")
       <<",\"exponential_underflow\":"<<(v.exponentialUnderflow?"true":"false")<<'}';
}
struct Case
{
    std::string name;double k=.5,nuA=.1,dk=.03,dnuA=.02;
    RANSTensor G{},H{{{.02,-.03,.01},{.04,-.01,.02},{-.02,.03,-.01}}};
    bool oneSided=false;
};
void Write(std::ostream &out,const Case &c,RANSWallInnerStressModel model)
{
    const double nu=1.5e-5,pi=.3;
    const auto base=EvaluateRANSWallInnerStress(c.k,c.nuA,c.G,model);
    const auto action=RANSWallInnerStressDirectional(c.k,c.nuA,c.G,c.dk,c.dnuA,c.H,model);
    const auto linear=EvaluateRANSWallInnerStress(c.k,c.nuA,c.G);
    if(base.linearResponse)
    {
        Check(base.reynoldsStress==linear.reynoldsStress && base.production==linear.production,
              "inactive and Linear stress/production must be exactly identical");
        Check(RANSWallInnerMeanStress(base,nu,pi)==RANSWallInnerMeanStress(linear,nu,pi),
              "inactive mean stress must retain legacy arithmetic");
    }
    if(model==RANSWallInnerStressModel::SaturatedR2 && c.k==0)
        for(const auto &row:base.reynoldsStress)for(double value:row)Check(value==0,"R2 zero-k covariance");
    double divergence=0,work=0;
    for(int i=0;i<3;++i){divergence+=c.G[i][i];for(int j=0;j<3;++j)work+=base.turbulentStress[i][j]*c.G[i][j];}
    Check(std::abs(base.production-(work-base.c*divergence))<2e-12,"full weak-divergence production");
    out<<"{\"name\":"<<std::quoted(c.name)<<",\"selector\":"<<std::quoted(RANSWallInnerStressName(model))
       <<",\"k\":"<<c.k<<",\"nuA\":"<<c.nuA<<",\"nu\":"<<nu<<",\"pi\":"<<pi<<",\"gradient\":";Tensor(out,c.G);
    out<<",\"dk\":"<<c.dk<<",\"dnuA\":"<<c.dnuA<<",\"direction\":";Tensor(out,c.H);
    out<<",\"one_sided\":"<<(c.oneSided?"true":"false")<<",\"value\":";Values(out,base,nu,pi);
    out<<",\"action\":{\"strain\":";Tensor(out,action.strain);out<<",\"turbulent\":";Tensor(out,action.turbulentStress);
    out<<",\"stress\":";Tensor(out,action.reynoldsStress);out<<",\"production\":"<<action.production
       <<",\"directional_only\":"<<(action.directionalOnly?"true":"false")<<"},\"differences\":[";
    bool first=true;for(double h:{1e-3,1e-4,1e-5,1e-6})
    {
        if(!first)out<<',';first=false;auto plus=c.G,minus=c.G;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j){plus[i][j]+=h*c.H[i][j];minus[i][j]-=h*c.H[i][j];}
        out<<"{\"step\":"<<h<<",\"plus\":";Values(out,EvaluateRANSWallInnerStress(c.k+h*c.dk,c.nuA+h*c.dnuA,plus,model),nu,pi);
        if(!c.oneSided){out<<",\"minus\":";Values(out,EvaluateRANSWallInnerStress(c.k-h*c.dk,c.nuA-h*c.dnuA,minus,model),nu,pi);}
        out<<'}';
    }out<<"]}";
}
}
int main(int argc,char **argv)
{
    try
    {
        Check(argc==3 && std::string(argv[1])=="--output","use --output REPORT.json");
        std::vector<Case> cases;
        Case active;active.name="archived_first_bad_point";active.k=.5105080130598861;
        active.nuA=.08437364924063438;active.G[0][2]=5.191590612064875;cases.push_back(active);
        Case reference=active;reference.name="neutral_reference_inactive";
        reference.G[0][2]=.16/(1.5e-5+reference.nuA);cases.push_back(reference);
        Case general;general.name="full_3d_active_weak_divergence";
        general.G={{{.3,.4,3.},{-.3,-.1,.5},{.2,-.4,-.1}}};cases.push_back(general);
        for(int column=0;column<9;++column)
        {
            auto c=general;c.name="full_gradient_column_"+std::to_string(column);c.dk=c.dnuA=0;c.H={};
            c.H[column/3][column%3]=1;cases.push_back(c);
        }
        Case zero;zero.name="zero_strain";cases.push_back(zero);
        zero.name="zero_auxiliary_viscosity";zero.nuA=0;zero.G[0][2]=2;zero.oneSided=true;cases.push_back(zero);
        zero.name="zero_k_positive_strain";zero.k=0;zero.nuA=.1;cases.push_back(zero);
        zero.name="joint_zero_k_strain";zero.G={};cases.push_back(zero);
        zero.name="joint_zero_k_viscosity";zero.nuA=0;zero.G[0][2]=2;cases.push_back(zero);
        zero.name="all_three_zero";zero.G={};cases.push_back(zero);
        Case join;join.name="saturation_join";join.G[0][2]=.8*(2*join.k/3)/(join.nuA*std::sqrt(2.));cases.push_back(join);
        std::ofstream out(argv[2]);Check(bool(out),"cannot write report");
        out<<std::setprecision(17)<<"{\"schema\":1,\"rho\":0.8,\"cases\":[";
        bool first=true;for(auto model:{RANSWallInnerStressModel::Linear,RANSWallInnerStressModel::SaturatedR2})for(const auto &c:cases)
        {if(!first)out<<',';first=false;Write(out,c,model);}
        out<<"],\"rejections\":[";first=true;
        auto reject=[&](const std::string &name,const std::function<void()> &operation)
        {
            try{operation();}catch(const std::exception &e)
            {if(!first)out<<',';first=false;out<<"{\"name\":"<<std::quoted(name)<<",\"message\":"<<std::quoted(e.what())<<'}';++checks;return;}
            Check(false,"invalid input accepted: "+name);
        };
        reject("negative_k",[]{EvaluateRANSWallInnerStress(-.1,.1,{},RANSWallInnerStressModel::SaturatedR2);});
        reject("negative_nuA",[]{EvaluateRANSWallInnerStress(.1,-.1,{});});
        reject("unknown_selector",[]{ParseRANSWallInnerStress("Clipped");});
        reject("outward_zero_k",[]{RANSWallInnerStressDirectional(0,.1,{},-.1,0,{},RANSWallInnerStressModel::SaturatedR2);});
        reject("outward_zero_nuA",[]{RANSWallInnerStressDirectional(.1,0,{},0,-.1,{},RANSWallInnerStressModel::SaturatedR2);});
        reject("nonfinite_gradient",[]{RANSTensor g{};g[0][0]=std::numeric_limits<double>::infinity();EvaluateRANSWallInnerStress(.1,.1,g);});
        out<<"],\"checks\":"<<checks<<",\"solver_or_trajectory\":false}\n";
        Check(bool(out),"failed to finish report");std::cout<<"Inner stress kernel checks: "<<checks<<'\n';return 0;
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 1;}
}
