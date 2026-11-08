///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSModels.h>
#include <SolverUtils/Turbulence/RANSRobustVariant.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace Nektar::SolverUtils;
namespace
{
std::vector<std::string> Split(const std::string &line)
{
    std::istringstream in(line);std::string value;std::vector<std::string> out;
    while (std::getline(in,value,',')) out.push_back(value);
    return out;
}
RANSTensor Stress(RANSModel model,double k,double eps,double eta,double nu,const RANSTensor &g)
{
    const auto c=EvaluateRANSClosure(model,k,eps,g,eta,{},nu);
    RANSTensor t{};
    for (int i=0;i<3;++i) for (int j=0;j<3;++j)
        t[i][j]=c.momentumDiffusion*(g[i][j]+g[j][i])-c.quadraticStress[i][j]-c.cubicStress[i][j];
    return t;
}
struct Minimum { double value=std::numeric_limits<double>::infinity(),angle=0;int side=0; };
Minimum Probe(RANSModel model,double k,double eps,double eta,double nu,double shear,double relative)
{
    RANSTensor g{};g[0][2]=shear;
    const auto base=Stress(model,k,eps,eta,nu,g);
    const double h=relative*std::max(1.,std::abs(shear));
    Minimum result;
    for (int step=0;step<=360;++step)
    {
        const double angle=std::acos(-1.)*step/360;
        const RANSVector n{{std::cos(angle),0.,std::sin(angle)}},a{{-std::sin(angle),0.,std::cos(angle)}};
        for (int side : {-1,1})
        {
            auto trial=g;
            for (int i=0;i<3;++i) for (int j=0;j<3;++j) trial[i][j]+=side*h*a[i]*n[j];
            const auto stress=Stress(model,k,eps,eta,nu,trial);
            double viscosity=0;
            for (int i=0;i<3;++i) for (int j=0;j<3;++j)
                viscosity+=a[i]*n[j]*(stress[i][j]-base[i][j])/(side*h);
            if (viscosity<result.value) result={viscosity,angle,side};
        }
    }
    return result;
}
}
// General frozen-k/epsilon physical one-sided work scan. This is a finite
// x-z incompressible rank-one scan, not a global coupled stability theorem.
int General(const char *inputPath,const char *outputPath)
{
    std::ifstream input(inputPath);std::ofstream out(outputPath);
    if(!input || !out)throw std::runtime_error("cannot open general tangent files");
    std::string line;std::getline(input,line);
    if(line!="id,distance,k,epsilon,beta,gamma,reference_shear,g00,g01,g02,g10,g11,g12,g20,g21,g22")
        throw std::runtime_error("unexpected general tangent CSV header");
    out<<std::setprecision(17)<<"id,h,min_work,angle,side,base_strain,base_rotation,trial_strain,trial_rotation,nut,c_mu";
    for(int a=0;a<3;++a)for(int b=0;b<3;++b){out<<",R"<<a<<b;}
    out<<'\n';
    while(std::getline(input,line))
    {
        if(line.empty()){continue;}
        const auto row=Split(line);if(row.size()!=16){throw std::runtime_error("general tangent row shape differs");}
        const double distance=std::stod(row[1]),k=std::stod(row[2]),epsilon=std::stod(row[3]);
        RANSRobustVariant variant;variant.beta=std::stod(row[4]);variant.gamma=std::stod(row[5]);
        variant.referenceStrain[0][2]=variant.referenceStrain[2][0]=std::stod(row[6]);
        RANSTensor g{};double scale=1;for(int a=0;a<3;++a)for(int b=0;b<3;++b){g[a][b]=std::stod(row[7+3*a+b]);scale=std::max(scale,std::abs(g[a][b]));}
        auto closure=[&](const RANSTensor &v){return EvaluateRANSRobustClosure(RANSModel::CLSHighRe,k,epsilon,v,distance,{},0,variant);};
        const auto base=closure(g);
        for(double relative:{1e-5,1e-6,1e-7})
        {
            const double h=relative*scale;Minimum minimum;RANSClosureValues critical{};
            for(int step=0;step<=720;++step)
            {
                const double angle=std::acos(-1.)*step/720;
                const RANSVector n{std::cos(angle),0,std::sin(angle)},a{-std::sin(angle),0,std::cos(angle)};
                for(int side:{-1,1})
                {
                    auto trial=g;for(int i=0;i<3;++i)for(int j=0;j<3;++j)trial[i][j]+=side*h*a[i]*n[j];
                    const auto c=closure(trial);double work=0;
                    for(int i=0;i<3;++i)for(int j=0;j<3;++j)work-=a[i]*n[j]*(c.reynoldsStress[i][j]-base.reynoldsStress[i][j])/(side*h);
                    if(work<minimum.value){minimum={work,angle,side};critical=c;}
                }
            }
            out<<row[0]<<','<<h<<','<<minimum.value<<','<<minimum.angle<<','<<minimum.side<<','
                <<base.strainInvariant<<','<<base.rotationInvariant<<','<<critical.strainInvariant<<','<<critical.rotationInvariant<<','<<base.eddyViscosity<<','<<base.cMu;
            for(const auto &v:base.reynoldsStress)for(double x:v){out<<','<<x;}
            out<<'\n';
        }
    }
    if(!out){throw std::runtime_error("general tangent write failed");}
    return 0;
}

int main(int argc,char **argv)
{
    if(argc==4 && std::string(argv[1])=="--general")
    {try{return General(argv[2],argv[3]);}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}}
    if (argc!=5)
    {
        std::cerr << "Usage: RANSStressTangentAudit standard|cls resolved_nu input.csv output.csv\n"
                     "Input: id,eta,k,epsilon,shear. Reports physical one-sided transverse stress tangents in 2D.\n";
        return 2;
    }
    try
    {
        const std::string name=argv[1];
        if (name!="standard" && name!="cls") throw std::invalid_argument("Unknown model");
        const auto model=name=="standard" ? RANSModel::StandardKEpsilon : RANSModel::CLSHighRe;
        const double nu=std::stod(argv[2]);
        std::ifstream input(argv[3]);std::ofstream output(argv[4]);
        if (!input || !output) throw std::runtime_error("Cannot open tangent audit files");
        output << std::setprecision(17) << "id,eta,k,epsilon,shear,normalized_shear,eddy_viscosity,min_transverse_tangent_h1e5,min_transverse_tangent_h1e6,angle_rad,side,R_xx,R_zz,R_xz\n";
        std::string line;std::getline(input,line);
        if (!line.empty() && line.back()=='\r') line.pop_back();
        if (line!="id,eta,k,epsilon,shear") throw std::runtime_error("Unexpected tangent input header");
        while (std::getline(input,line))
        {
            if (!line.empty() && line.back()=='\r') line.pop_back();
            const auto row=Split(line);
            if (row.size()!=5) throw std::runtime_error("Invalid tangent input row");
            const double eta=std::stod(row[1]),k=std::stod(row[2]),eps=std::stod(row[3]),shear=std::stod(row[4]);
            RANSTensor g{};g[0][2]=shear;
            const auto c=EvaluateRANSClosure(model,k,eps,g,eta,{},nu);
            const auto coarse=Probe(model,k,eps,eta,nu,shear,1e-5),fine=Probe(model,k,eps,eta,nu,shear,1e-6);
            output << row[0] << ',' << eta << ',' << k << ',' << eps << ',' << shear << ',' << k*shear/eps
                << ',' << c.eddyViscosity << ',' << coarse.value << ',' << fine.value << ',' << fine.angle << ',' << fine.side
                << ',' << c.reynoldsStress[0][0] << ',' << c.reynoldsStress[2][2] << ',' << c.reynoldsStress[0][2] << '\n';
        }
        if (!output) throw std::runtime_error("Cannot finish tangent output");
    }
    catch (const std::exception &error) { std::cerr << error.what() << '\n';return 1; }
}
