///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_SOLVE_DIAGNOSTICS_H
#define NEKTAR_RANS_WALL_SOLVE_DIAGNOSTICS_H
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

namespace Nektar::RANSWallSolveJson
{
// Small diagnostic-only JSON writer. Nonfinite observations are explicit nulls;
// the physical residual retains its original finite-value rejection rules.
inline std::string Quote(const std::string &value)
{
    std::ostringstream out;out<<'"';
    for(unsigned char c:value)
    {
        if(c=='"' || c=='\\')out<<'\\'<<c;
        else if(c<32)out<<"\\u"<<std::hex<<std::setw(4)<<std::setfill('0')<<int(c)<<std::dec;
        else out<<c;
    }
    out<<'"';return out.str();
}
inline std::string Number(double value)
{
    if(!std::isfinite(value))return "null";
    std::ostringstream out;out<<std::setprecision(17)<<value;return out.str();
}
inline std::string Boolean(bool value){return value?"true":"false";}
inline std::string List(const std::vector<std::string> &values)
{
    std::ostringstream out;out<<'[';
    for(std::size_t i=0;i<values.size();++i){if(i)out<<',';out<<values[i];}
    out<<']';return out.str();
}
template<class Range> std::string Numbers(const Range &values)
{
    std::vector<std::string> result;for(const auto &value:values)result.push_back(Number(value));
    return List(result);
}
template<class Matrix> std::string Rows(const Matrix &values)
{
    std::vector<std::string> result;for(const auto &row:values)result.push_back(Numbers(row));
    return List(result);
}
inline std::string Object(const std::vector<std::pair<std::string,std::string>> &values)
{
    std::ostringstream out;out<<'{';
    for(std::size_t i=0;i<values.size();++i)
    {if(i)out<<',';out<<Quote(values[i].first)<<':'<<values[i].second;}
    out<<'}';return out.str();
}
inline std::string ExceptionType(const std::exception &error)
{
    if(dynamic_cast<const std::invalid_argument*>(&error))return "std::invalid_argument";
    if(dynamic_cast<const std::domain_error*>(&error))return "std::domain_error";
    if(dynamic_cast<const std::overflow_error*>(&error))return "std::overflow_error";
    if(dynamic_cast<const std::out_of_range*>(&error))return "std::out_of_range";
    if(dynamic_cast<const std::runtime_error*>(&error))return "std::runtime_error";
    return typeid(error).name();
}
template<class Tensor> std::array<double,3> Eigenvalues(Tensor a)
{
    for(int step=0;step<48;++step)
    {
        int p=0,q=1;for(int i=0;i<3;++i)for(int j=i+1;j<3;++j)
            if(std::abs(a[i][j])>std::abs(a[p][q])){p=i;q=j;}
        if(std::abs(a[p][q])<1e-15*std::max({1.,std::abs(a[0][0]),std::abs(a[1][1]),std::abs(a[2][2])}))break;
        const double angle=.5*std::atan2(2*a[p][q],a[q][q]-a[p][p]),c=std::cos(angle),s=std::sin(angle);
        const double pp=c*c*a[p][p]-2*s*c*a[p][q]+s*s*a[q][q];
        const double qq=s*s*a[p][p]+2*s*c*a[p][q]+c*c*a[q][q];
        for(int k=0;k<3;++k)if(k!=p && k!=q)
        {const double x=a[k][p],y=a[k][q];a[k][p]=a[p][k]=c*x-s*y;a[k][q]=a[q][k]=s*x+c*y;}
        a[p][p]=pp;a[q][q]=qq;a[p][q]=a[q][p]=0;
    }
    std::array<double,3> result{a[0][0],a[1][1],a[2][2]};std::sort(result.begin(),result.end());return result;
}
} // namespace Nektar::RANSWallSolveJson
#endif
