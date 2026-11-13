// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_MANUFACTURED_JET_H
#define NEKTAR_RANS_WALL_MANUFACTURED_JET_H
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
namespace Nektar::CT6
{
// Truncated Taylor algebra in (q,r,t), coefficients include 1/(i!j!k!).
// This differentiates continuous manufactured data, never native FE loads.
struct Jet
{
    std::array<double,35> c{};
    struct Powers
    {
        std::vector<std::array<int,3>> p;
        int index[5][5][5];
        std::vector<std::array<int,3>> products;
        Powers()
        {
            for(auto &a:index)for(auto &b:a)for(int &i:b)i=-1;
            for(int n=0;n<=4;++n)for(int i=0;i<=n;++i)for(int j=0;j<=n-i;++j)
            {const int k=n-i-j;index[i][j][k]=p.size();p.push_back({i,j,k});}
            for(int i=0;i<35;++i)for(int j=0;j<35;++j)
            {const auto a=p[i],b=p[j];if(a[0]+b[0]+a[1]+b[1]+a[2]+b[2]<=4)
                products.push_back({i,j,index[a[0]+b[0]][a[1]+b[1]][a[2]+b[2]]});}
        }
    };
    static const Powers &Table(){static const Powers v;return v;}
    Jet(double v=0){c[0]=v;}
    static Jet Variable(double v,int a){Jet j(v);std::array<int,3> p{};p[a]=1;j.c[Table().index[p[0]][p[1]][p[2]]]=1;return j;}
    double Value() const{return c[0];}
    Jet D(int a) const
    {
        Jet r;const auto &v=Table();
        for(int i=0;i<35;++i){auto p=v.p[i];if(p[a]){const int n=p[a]--;r.c[v.index[p[0]][p[1]][p[2]]]=n*c[i];}}
        return r;
    }
};
inline Jet operator+(Jet a,const Jet &b){for(int i=0;i<35;++i)a.c[i]+=b.c[i];return a;}
inline Jet operator-(Jet a,const Jet &b){for(int i=0;i<35;++i)a.c[i]-=b.c[i];return a;}
inline Jet operator-(Jet a){for(double &v:a.c)v=-v;return a;}
inline Jet operator*(const Jet &a,const Jet &b)
{Jet r;for(const auto &p:Jet::Table().products)r.c[p[2]]+=a.c[p[0]]*b.c[p[1]];return r;}
inline Jet operator*(Jet a,double b){for(double &v:a.c)v*=b;return a;}
inline Jet operator*(double b,Jet a){return a*b;}
inline Jet Pow(const Jet &a,double p)
{
    if(!(a.Value()>0))throw std::domain_error("CT6 Taylor power requires positive base");
    Jet d=a-a.Value(),v(1),r;double b=std::pow(a.Value(),p);
    for(int k=0;k<=4;++k){r=r+v*b;v=v*d;b*=((p-k)/(k+1))/a.Value();}return r;
}
inline Jet Inverse(const Jet &a)
{
    // A reciprocal is defined on either sign branch. Fractional powers and
    // logarithms have the separate positive-base domain used above/below.
    if(a.Value()==0 || !std::isfinite(a.Value()))
        throw std::domain_error("CT6 Taylor reciprocal requires finite nonzero base");
    const Jet d=(a-a.Value())*(-1/a.Value());Jet v(1),r;
    for(int k=0;k<=4;++k){r=r+v*(1/a.Value());v=v*d;}return r;
}
inline Jet operator/(const Jet &a,const Jet &b){return a*Inverse(b);}
inline Jet operator/(const Jet &a,double b){return a*(1/b);}
inline Jet Sqrt(const Jet &a){return Pow(a,.5);}
inline Jet Exp(const Jet &a)
{Jet d=a-a.Value(),v(1),r;double b=std::exp(a.Value());for(int k=0;k<=4;++k){r=r+v*b;v=v*d;b/=(k+1);}return r;}
inline Jet Log(const Jet &a)
{
    if(!(a.Value()>0))throw std::domain_error("CT6 Taylor log requires positive base");
    Jet d=(a-a.Value())/a.Value(),v(1),r(std::log(a.Value()));
    for(int k=1;k<=4;++k){v=v*d;r=r+v*((k%2?1.:-1.)/k);}return r;
}
using Vec=std::array<Jet,3>;
using Mat=std::array<Vec,3>;
inline Jet Dot(const Vec &a,const Vec &b){Jet x;for(int j=0;j<3;++j)x=x+a[j]*b[j];return x;}
inline Mat Product(const Mat &a,const Mat &b)
{Mat r{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)r[i][j]=r[i][j]+a[i][k]*b[k][j];return r;}
} // namespace Nektar::CT6
#endif
