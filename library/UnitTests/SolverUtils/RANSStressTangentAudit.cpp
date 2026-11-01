///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSModels.h>
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
int main(int argc,char **argv)
{
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
