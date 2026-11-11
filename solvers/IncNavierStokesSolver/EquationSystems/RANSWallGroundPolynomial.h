// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_GROUND_POLYNOMIAL_H
#define NEKTAR_RANS_WALL_GROUND_POLYNOMIAL_H
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <vector>
namespace Nektar::RANSWallGroundPolynomial
{
using Polynomial=std::vector<long double>;
inline long double Value(const Polynomial &a,long double x)
{long double y=0;for(auto i=a.rbegin();i!=a.rend();++i)y=y*x+*i;return y;}
inline Polynomial Derivative(const Polynomial &a)
{Polynomial d;for(int j=1;j<int(a.size());++j)d.push_back(j*a[j]);return d;}
// Recursive Rolle isolation partitions at all derivative roots. Repeated
// roots are explicitly retained using a coefficient-scaled roundoff bound.
// No fixed sampling spacing can miss a narrow negative interval.
inline std::vector<long double> Roots(Polynomial a)
{
    while(a.size()>1 && a.back()==0)a.pop_back();
    if(a.size()<2)return {};
    if(a.size()==2){const auto x=-a[0]/a[1];return std::abs(x)<1?std::vector<long double>{x}:std::vector<long double>{};}
    auto critical=Roots(Derivative(a));std::vector<long double> partition{-1};
    partition.insert(partition.end(),critical.begin(),critical.end());partition.push_back(1);
    long double norm=0;for(auto c:a)norm+=std::abs(c);
    const long double tolerance=512*std::numeric_limits<long double>::epsilon()*norm;
    std::vector<long double> out;
    for(auto x:critical)if(std::abs(Value(a,x))<=tolerance)out.push_back(x);
    for(int j=1;j<int(partition.size());++j)
    {
        long double l=partition[j-1],r=partition[j],fl=Value(a,l),fr=Value(a,r);
        if(!((fl<0 && fr>0)||(fl>0 && fr<0)))continue;
        for(int n=0;n<100;++n)
        {
            const auto mid=(l+r)/2,fm=Value(a,mid);
            if(fm==0){l=r=mid;break;}
            if((fl<0)==(fm<0)){l=mid;fl=fm;}else r=mid;
        }
        out.push_back((l+r)/2);
    }
    std::sort(out.begin(),out.end());
    out.erase(std::unique(out.begin(),out.end()),out.end());return out;
}
inline Polynomial ModifiedA(const std::array<long double,6> &c)
{
    // (1-x)/2,(1+x)/2,(1-x*x) P^(1,1)_{j-2}(x)/4.
    return {(c[0]+c[1])/2+c[2]/4-3*c[4]/16,
        (c[1]-c[0])/2+c[3]/2-3*c[5]/4,-c[2]/4+18*c[4]/16,
        -c[3]/2+10*c[5]/4,-15*c[4]/16,-7*c[5]/4};
}
inline Polynomial ModifiedA6(const std::array<long double,6> &c,long double c6)
{
    auto p=ModifiedA(c);
    // Preserve the established p3/p4/p5 evaluation exactly. The new bubble is
    // (1-x^2) P_4^(1,1)(x)/4 = (5-75x^2+175x^4-105x^6)/32.
    if(c6==0)return p;
    p.resize(7);p[0]+=5*c6/32;p[2]-=75*c6/32;
    p[4]+=175*c6/32;p[6]=-105*c6/32;return p;
}
}
#endif
