// SPDX-License-Identifier: MIT
// Exact native WJ directional action versus the actual native wall adapter.
// Input rows: label k epsilon shear nu nuA r H. No PDE/production claims.
#include <IncNavierStokesSolver/EquationSystems/WJNativeFlow.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

using namespace Nektar;
using namespace Nektar::WJNative;
namespace
{
struct State
{
    std::string label;
    double k, e, shear, nu, nuA, r, H;
    Tensor G() const { Tensor g{}; g[0][2]=shear; return g; }
    WallStress Value(const Tensor& g) const
    { return Adapt(k,e,g,{0,0,1},nu,nuA,r,H); }
    Tensor Action(const Tensor& h) const
    {
        const auto dw=SolverUtils::ApplyWJClosureDirection(k,e,G(),0,0,h,
            {},SolverUtils::WJGradientPolicy::DevStrainRawWorkV1,0);
        const double chi=Value(G()).chi;
        Tensor a{};
        for(int i=0;i<3;++i) for(int j=0;j<3;++j)
        {
            // n=e_z: projector retains the tangential block and zz row.
            const bool mixed=(i==2)!=(j==2);
            const double drh=mixed ? -nuA*(h[i][j]+h[j][i])
                                   : dw.reynoldsStress[i][j];
            const double dr=(1-chi)*drh+chi*dw.reynoldsStress[i][j];
            a[i][j]=nu*(h[i][j]+h[j][i])-dr;
        }
        return a;
    }
};
std::array<Vector,2> Directions(double theta)
{
    return {Vector{std::cos(theta),0,std::sin(theta)},
            Vector{-std::sin(theta),0,std::cos(theta)}};
}
Tensor RankOne(double theta)
{
    auto v=Directions(theta);Tensor h{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)h[i][j]=v[1][i]*v[0][j];
    return h;
}
double Contract(const Tensor& a,double theta)
{
    auto v=Directions(theta);double s=0;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)s+=v[1][i]*a[i][j]*v[0][j];
    return s;
}
double Quality(const State& s,double theta)
{ return Contract(s.Action(RankOne(theta)),theta); }
void PrintTensor(const Tensor& a)
{
    std::cout<<'[';
    for(int i=0;i<3;++i){if(i)std::cout<<',';std::cout<<'[';
        for(int j=0;j<3;++j){if(j)std::cout<<',';std::cout<<a[i][j];}
        std::cout<<']';}
    std::cout<<']';
}
}
int main(int argc,char** argv)
{
    try
    {
        if(argc!=2)throw std::runtime_error("usage: WJH27PrincipalProbe states.tsv");
        std::ifstream in(argv[1]);if(!in)throw std::runtime_error("missing state input");
        std::cout<<std::setprecision(17)
            <<"{\"schema\":\"WJH27NativePrincipalV1\",\"scope\":\"constitutive feasibility only\",\"states\":[";
        bool first=true;State s;
        while(in>>s.label>>s.k>>s.e>>s.shear>>s.nu>>s.nuA>>s.r>>s.H)
        {
            if(!first)std::cout<<',';
            first=false;
            const auto value=s.Value(s.G());
            double best=std::numeric_limits<double>::infinity(),theta=0;
            constexpr int samples=2048;const double pi=std::acos(-1.);
            for(int j=0;j<samples;++j)
            {double t=pi*j/samples,a=Quality(s,t);if(a<best){best=a;theta=t;}}
            // Refine the selected local angular minimum without any stress change.
            double left=theta-pi/samples,right=theta+pi/samples;
            for(int j=0;j<55;++j)
            {double a=left+(right-left)/3,b=right-(right-left)/3;
             if(Quality(s,a)<Quality(s,b))right=b;else left=a;}
            theta=(left+right)/2;best=Quality(s,theta);
            std::cout<<"{\"label\":\""<<s.label<<"\",\"k\":"<<s.k
                <<",\"epsilon\":"<<s.e<<",\"shear\":"<<s.shear
                <<",\"nu\":"<<s.nu<<",\"nuA\":"<<s.nuA
                <<",\"r\":"<<s.r<<",\"H\":"<<s.H
                <<",\"mu\":"<<value.mu<<",\"chi\":"<<value.chi
                <<",\"production\":"<<value.production<<",\"R\":";
            PrintTensor(value.R);
            std::cout<<",\"planar_min\":"<<best
                <<",\"planar_min_over_mu\":"<<best/value.mu
                <<",\"theta\":"<<theta
                <<",\"theta30_quality\":"<<Quality(s,pi/6)
                <<",\"finite_difference\":[";
            bool fdFirst=true;
            for(double step:{2e-4,1e-4,2e-5,1e-5})
            {
                const double delta=step*s.e/s.k;auto plus=s.G(),minus=plus;
                auto h=RankOne(theta);
                for(int i=0;i<3;++i)for(int j=0;j<3;++j)
                {plus[i][j]+=delta*h[i][j];minus[i][j]-=delta*h[i][j];}
                const auto rp=s.Value(plus),rm=s.Value(minus);Tensor a{};
                for(int i=0;i<3;++i)for(int j=0;j<3;++j)
                    a[i][j]=s.nu*(h[i][j]+h[j][i])
                            -(rp.R[i][j]-rm.R[i][j])/(2*delta);
                if(!fdFirst)std::cout<<',';
                fdFirst=false;
                const double fd=Contract(a,theta);
                std::cout<<"{\"dimensionless_step\":"<<step
                    <<",\"quality\":"<<fd
                    <<",\"difference_over_mu\":"<<(fd-best)/value.mu<<'}';
            }
            std::cout<<"]}";
        }
        if(!in.eof())throw std::runtime_error("malformed state input");
        if(first)throw std::runtime_error("empty state input");
        std::cout<<"]}\n";return 0;
    }
    catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
}
