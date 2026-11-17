// SPDX-License-Identifier: MIT
#ifndef WJ_H26_COMMON_H
#define WJ_H26_COMMON_H
#include <IncNavierStokesSolver/EquationSystems/WJNativeFlow.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <boost/property_tree/json_parser.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <filesystem>
#include <chrono>
using namespace Nektar;
using namespace Nektar::WJNative;
namespace fs=std::filesystem;
namespace pt=boost::property_tree;
double Clock(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
void VectorJSON(std::ostream&o,const Vec&v){o<<'[';for(int i=0;i<v.size();++i){if(i)o<<',';o<<v[i];}o<<']';}
void MatrixJSON(std::ostream&o,const Sparse&a){o<<"{\"rows\":"<<a.rows()<<",\"cols\":"<<a.cols()<<",\"entries\":[";bool comma=false;for(int j=0;j<a.outerSize();++j)for(Sparse::InnerIterator it(a,j);it;++it){if(comma)o<<',';comma=true;o<<'['<<it.row()<<','<<it.col()<<','<<it.value()<<']';}o<<"]}";}
void InfoJSON(std::ostream&o,const SolveInfo&i){o<<"{\"iterations\":"<<i.iterations<<",\"rebuilds\":"<<i.rebuilds<<",\"momentum\":"<<i.momentum<<",\"continuity\":"<<i.continuity<<",\"gauge\":"<<i.gauge<<",\"pressure_correction\":"<<i.correction<<",\"factor_seconds\":"<<i.factorSeconds<<",\"solve_seconds\":"<<i.solveSeconds<<",\"residual_seconds\":"<<i.residualSeconds<<'}';}
std::array<Function,3>Polynomial(double t)
{
    double a=std::exp(-t),b=std::exp(-2*t);
    return {[=](double x,double z){return a*x*x*(1-x)*(2*z-3*z*z);},[=](double x,double z){return-a*(2*x-3*x*x)*z*z*(1-z);},[=](double x,double z){return b*(x*x+z*z);}};
}
Vec Exact(Flow&flow,const std::array<Function,3>&f)
{
    Vec x(2*flow.VelocitySize()+flow.PressureSize());for(int c=0;c<2;++c)x.segment(c*flow.VelocitySize(),flow.VelocitySize())=flow.Project(f[c],c);x.tail(flow.PressureSize())=flow.Project(f[2],2);
    if(flow.Gauge()){Vec one=flow.Project([](double,double){return 1.;},2);x.tail(flow.PressureSize())-=one*(flow.PressureMean().dot(x.tail(flow.PressureSize()))/flow.PressureMean().dot(one));}return x;
}
#endif
