// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_CUTOFF_PANELS_H
#define NEKTAR_WJ_NATIVE_CUTOFF_PANELS_H
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Nektar::WJNative
{
/** Interior integration partitions for a Modified_A cubic edge trace.
 * K levels are inverse-cutoff values at fixed radial faces, not constraints
 * on the actual field. Every stationary point is retained when a level can
 * meet the trace; sign-changing roots are isolated on monotone intervals.
 * A constant trace identically on a face needs no interior x partition.
 */
inline std::vector<double> CubicCutoffPanels(const std::array<double,4>&c,
                                            const std::vector<double>&levels)
{
    for(double v:c)if(!std::isfinite(v))throw std::invalid_argument("nonfinite cubic cutoff trace");
    for(double v:levels)if(!std::isfinite(v))throw std::invalid_argument("nonfinite cutoff face level");
    if(levels.empty())return{};
    using Real=long double;
    const Real a=c[0],b=-Real(c[0])+c[1]+c[2]-2*Real(c[3]),
        d=-Real(c[2])+6*Real(c[3]),e=-4*Real(c[3]);
    auto value=[&](Real t){return a+t*(b+t*(d+t*e));};
    std::vector<Real>critical;
    auto addCritical=[&](Real t){if(std::isfinite(t)&&t>0&&t<1)critical.push_back(t);};
    if(e==0){if(d!=0)addCritical(-b/(2*d));}
    else{
        const Real discriminant=d*d-3*e*b;
        if(discriminant>=0){
            const Real q=-d-std::copysign(std::sqrt(discriminant),d);
            if(q!=0){addCritical(q/(3*e));addCritical(b/q);}
            else addCritical(-d/(3*e));
        }
    }
    std::sort(critical.begin(),critical.end());critical.erase(std::unique(critical.begin(),critical.end()),critical.end());
    std::vector<Real>bounds{0};bounds.insert(bounds.end(),critical.begin(),critical.end());bounds.push_back(1);
    Real lower=value(0),upper=lower;
    for(Real t:bounds){lower=std::min(lower,value(t));upper=std::max(upper,value(t));}
    const Real arithmetic=128*std::numeric_limits<Real>::epsilon()*(std::abs(a)+std::abs(b)+std::abs(d)+std::abs(e));
    bool possible=false;
    for(double level:levels)possible|=Real(level)>=lower-arithmetic&&Real(level)<=upper+arithmetic;
    if(!possible)return{};
    std::vector<double>result;
    for(Real t:critical)result.push_back(double(t));
    for(double inputLevel:levels){
        const Real level=inputLevel;
        for(Real t:critical)if(value(t)==level)result.push_back(double(t));
        for(std::size_t i=0;i+1<bounds.size();++i){
            Real left=bounds[i],right=bounds[i+1],fl=value(left)-level,fr=value(right)-level;
            if(!((fl<0&&fr>0)||(fl>0&&fr<0)))continue;
            for(int j=0;j<100;++j){
                const Real middle=(left+right)/2;
                if(middle==left||middle==right)break;
                const Real fm=value(middle)-level;
                if(fm==0){left=right=middle;break;}
                if((fm<0)==(fl<0)){left=middle;fl=fm;}else{right=middle;fr=fm;}
            }
            const Real root=(left+right)/2;
            if(root>0&&root<1)result.push_back(double(root));
        }
    }
    std::sort(result.begin(),result.end());result.erase(std::unique(result.begin(),result.end()),result.end());
    return result;
}
}
#endif
