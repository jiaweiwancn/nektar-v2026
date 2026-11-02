///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSWallNormalMap.h>
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
using Vec=RANSVector;
using Tensor=RANSTensor;
const Tensor stress{{Vec{.7,.1,-.2},Vec{.1,-.3,.05},Vec{-.2,.05,1.1}}};
int checks=0;
void Check(bool good,const char *message){++checks;if(!good)throw std::runtime_error(message);}
double Dot(const Vec &a,const Vec &b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
Vec Product(const Tensor &a,const Vec &b){return {Dot(a[0],b),Dot(a[1],b),Dot(a[2],b)};}
double Maximum(const Vec &v){return std::max({std::abs(v[0]),std::abs(v[1]),std::abs(v[2])});}
void Write(std::ostream &o,const Vec &v){o<<'['<<v[0]<<','<<v[1]<<','<<v[2]<<']';}
void Write(std::ostream &o,const Tensor &v){o<<'[';for(int i=0;i<3;++i){if(i){o<<',';}Write(o,v[i]);}o<<']';}
struct Rule{std::vector<double> x,w;};
Rule Gauss(int n)
{
    Rule rule;rule.x.resize(n);rule.w.resize(n);const long double pi=std::acos(-1.L);
    for(int i=0;i<(n+1)/2;++i)
    {
        long double z=std::cos(pi*(i+.75L)/(n+.5L)),derivative=0;
        bool converged=false;
        for(int iteration=0;iteration<50;++iteration)
        {
            long double previous=1,value=z;
            for(int j=2;j<=n;++j){const long double next=((2*j-1)*z*value-(j-1)*previous)/j;previous=value;value=next;}
            derivative=n*(z*value-previous)/(z*z-1);const long double step=value/derivative;z-=step;
            if(std::abs(step)<8*std::numeric_limits<long double>::epsilon()){converged=true;break;}
        }
        Check(converged,"Gauss-Legendre root iteration failed");
        long double previous=1,value=z;
        for(int j=2;j<=n;++j){const long double next=((2*j-1)*z*value-(j-1)*previous)/j;previous=value;value=next;}
        derivative=n*(z*value-previous)/(z*z-1);
        const double weight=static_cast<double>(2/((1-z*z)*derivative*derivative));
        rule.x[i]=static_cast<double>(-z);rule.x[n-1-i]=static_cast<double>(z);rule.w[i]=rule.w[n-1-i]=weight;
    }
    return rule;
}
RANSWallNormalMap Mapping()
{
    RANSWallNormalMap m;m.kind=RANSWallNormalMap::Kind::Cylinder;m.radius=20;
    m.height=2;m.heightAmplitude=.2;m.heightWaveX=.3;return m;
}
struct Field
{
    RANSWallNormalMap::Sample metric;
    Vec velocity{},acceleration{},scalarGradient{},scalarFlux{};
    Tensor gradient{};
    double scalar=0,diffusivity=0,source=0,divergence=0,normalVelocity=0,tangentialVelocity=0;
};
Field Evaluate(const RANSWallNormalMap &m,const Vec &q)
{
    Field f;f.metric=m.At(q);const auto &g=f.metric;const double s=q[0],r=g.distance,c=1/m.radius,J=g.J;
    const double a=1+.1*std::cos(.4*s),ap=-.04*std::sin(.4*s),app=-.016*std::cos(.4*s);
    const double u=a*r,w=-ap*r*r/(2*J),us=ap*r,ur=a,ws=-app*r*r/(2*J),wr=-ap*r*(2+c*r)/(2*J*J);
    const double Hp=.06*std::cos(.3*s),rs=q[2]*Hp/m.height;
    Vec fixedS{},fixedR{};Tensor chartGradient{};
    for(int i=0;i<3;++i)
    {
        f.velocity[i]=u*g.tangent1[i]+w*g.normal[i];
        fixedS[i]=(us+c*w)*g.tangent1[i]+(ws-c*u)*g.normal[i];
        fixedR[i]=ur*g.tangent1[i]+wr*g.normal[i];
        chartGradient[i][0]=fixedS[i]+rs*fixedR[i];chartGradient[i][2]=g.scale*fixedR[i];
    }
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)f.gradient[i][j]+=chartGradient[i][k]*g.inverse[k][j];
    f.divergence=f.gradient[0][0]+f.gradient[1][1]+f.gradient[2][2];f.acceleration=Product(f.gradient,f.velocity);
    f.normalVelocity=Dot(f.velocity,g.normal);f.tangentialVelocity=Dot(f.velocity,g.tangent1);
    f.scalar=2+.1*std::sin(s)+.2*r*r;f.diffusivity=1+.1*r;
    const double ks=.1*std::cos(s),kss=-.1*std::sin(s),kr=.4*r,krr=.4;
    for(int j=0;j<3;++j)
    {
        f.scalarGradient[j]=(ks+kr*rs)*g.inverse[0][j]+kr*g.scale*g.inverse[2][j];
        f.scalarFlux[j]=f.velocity[j]*f.scalar-f.diffusivity*f.scalarGradient[j];
    }
    const double diffusion=f.diffusivity*krr+(.1+c*f.diffusivity/J)*kr+f.diffusivity*kss/(J*J);
    f.source=Dot(f.velocity,f.scalarGradient)-diffusion;return f;
}
struct Flux{double mass=0,scalar=0;Vec momentum{};};
void Add(Flux &a,const Flux &b){a.mass+=b.mass;a.scalar+=b.scalar;for(int i=0;i<3;++i)a.momentum[i]+=b.momentum[i];}
Flux Boundary(const Field &f,const Vec &area,double weight)
{
    Flux b;b.mass=weight*Dot(f.velocity,area);b.scalar=weight*Dot(f.scalarFlux,area);
    for(int i=0;i<3;++i)b.momentum[i]=f.velocity[i]*b.mass-weight*Dot(stress[i],area);
    return b;
}
void Write(std::ostream &o,const Flux &f)
{o<<"{\"mass\":"<<f.mass<<",\"momentum\":";Write(o,f.momentum);o<<",\"scalar\":"<<f.scalar<<'}';}
struct Level
{
    int order;Flux faces[4],total,volume;double maximumDivergence=0,pressure=0;Vec momentumError{};
    double massError=0,scalarError=0;
};
Level Integrate(const RANSWallNormalMap &m,int order)
{
    const auto rule=Gauss(order);Level result{};result.order=order;constexpr double left=-3,right=2.5;
    for(int i=0;i<order;++i)
    {
        const double s=(left+right)/2+(right-left)*rule.x[i]/2,sw=rule.w[i]*(right-left)/2;
        for(int top=0;top<2;++top)
        {
            const Vec q{s,0,top?m.height:0.};const auto field=Evaluate(m,q);
            Add(result.faces[top?0:1],Boundary(field,m.Cofactor(q,2,top?1:-1),sw));
        }
        for(int j=0;j<order;++j)
        {
            const double rho=m.height*(rule.x[j]+1)/2,rw=m.height*rule.w[j]/2;
            const Vec q{s,0,rho};const auto f=Evaluate(m,q);const double weight=sw*rw*f.metric.determinant;
            result.volume.mass+=weight*f.divergence;result.volume.scalar+=weight*f.source;
            for(int a=0;a<3;++a)result.volume.momentum[a]+=weight*f.acceleration[a];
            result.maximumDivergence=std::max(result.maximumDivergence,std::abs(f.divergence));
        }
        const double rho=m.height*(rule.x[i]+1)/2,rw=m.height*rule.w[i]/2;
        for(int side=0;side<2;++side)
        {
            const Vec q{side?right:left,0,rho};const auto field=Evaluate(m,q);
            Add(result.faces[side?3:2],Boundary(field,m.Cofactor(q,0,side?1:-1),rw));
        }
        const double r=rule.x[i]+1;const auto g=m.At({0,0,r});const auto t=g.tangent1;
        const double curvature=Dot(t,Product(g.Br,t));result.pressure+=rule.w[i]*r*r*curvature;
    }
    for(const auto &face:result.faces)Add(result.total,face);
    result.massError=result.total.mass-result.volume.mass;result.scalarError=result.total.scalar-result.volume.scalar;
    for(int a=0;a<3;++a)result.momentumError[a]=result.total.momentum[a]-result.volume.momentum[a];
    Check(std::abs(result.total.mass)<1e-10 && std::abs(result.massError)<1e-10,"closed mass budget failed");
    Check(Maximum(result.momentumError)<1e-10,"closed Cartesian momentum budget failed");
    Check(std::abs(result.scalarError)<1e-10,"closed scalar flux/source budget failed");
    Check(result.maximumDivergence<2e-12,"analytic mapped velocity is not incompressible");
    return result;
}
void Write(std::ostream &o,const Level &v)
{
    o<<"{\"order\":"<<v.order<<",\"faces\":{";const char *names[]{"top","ground","left","right"};
    for(int i=0;i<4;++i){if(i){o<<',';}o<<'"'<<names[i]<<"\":";Write(o,v.faces[i]);}
    o<<"},\"boundary\":";Write(o,v.total);o<<",\"volume\":";Write(o,v.volume);
    o<<",\"maximum_divergence\":"<<v.maximumDivergence<<",\"mass_error\":"<<v.massError<<",\"momentum_error\":";Write(o,v.momentumError);
    o<<",\"scalar_error\":"<<v.scalarError<<",\"pressure_increment\":"<<v.pressure<<'}';
}
void FieldRecord(std::ostream &o,const Field &f)
{
    o<<"{\"chart\":";Write(o,f.metric.chart);o<<",\"point\":";Write(o,f.metric.point);o<<",\"normal\":";Write(o,f.metric.normal);
    o<<",\"tangent\":";Write(o,f.metric.tangent1);o<<",\"U\":";Write(o,f.velocity);o<<",\"G\":";Write(o,f.gradient);
    o<<",\"r\":"<<f.metric.distance<<",\"J\":"<<f.metric.J<<",\"H\":"<<f.metric.height<<",\"determinant\":"<<f.metric.determinant
     <<",\"u_t\":"<<f.tangentialVelocity<<",\"w\":"<<f.normalVelocity<<",\"divergence\":"<<f.divergence
     <<",\"k\":"<<f.scalar<<",\"D\":"<<f.diffusivity<<",\"grad_k\":";Write(o,f.scalarGradient);
    o<<",\"scalar_flux\":";Write(o,f.scalarFlux);o<<",\"source\":"<<f.source<<'}';
}
}
int main(int argc,char **argv)
{
    if(argc!=3 || std::string(argv[1])!="--output")
    {std::cerr<<"usage: RANSCurvedNoteProbe --output REPORT.json\n";return 1;}
    try
    {
        const auto map=Mapping();map.Validate({-3,0,0},{2.5,0,2});
        std::ostringstream out;out<<std::setprecision(17)<<"{\"schema\":1,\"success\":true,\"radius\":"<<map.radius
            <<",\"curvature\":"<<1/map.radius<<",\"H0\":"<<map.height<<",\"lower\":[-3,0,0],\"upper\":[2.5,0,2],\"stress\":";Write(out,stress);
        out<<",\"levels\":[";int count=0;double pressure=0;
        for(int order:{24,48,64}){const auto level=Integrate(map,order);pressure=level.pressure;if(count++){out<<',';}Write(out,level);}out<<"],\"finite_difference_samples\":[";count=0;
        for(double s:{-2.,.7,2.1})for(double r:{.1,.8,1.7})
        {
            const double height=map.At({s,0,0}).height;const auto field=Evaluate(map,{s,0,r*map.height/height});Tensor fd{};
            const double step=2e-5;
            for(int a:{0,2})
            {
                auto plus=field.metric.point,minus=plus;plus[a]+=step;minus[a]-=step;
                const auto vp=Evaluate(map,map.Inverse(plus)),vm=Evaluate(map,map.Inverse(minus));
                for(int b=0;b<3;++b)fd[b][a]=(vp.velocity[b]-vm.velocity[b])/(2*step);
            }
            double error=0;for(int a=0;a<3;++a)for(int b=0;b<3;++b)error=std::max(error,std::abs(fd[a][b]-field.gradient[a][b]));
            Check(error<2e-10,"Cartesian inverse-map finite difference disagrees with mapped gradient");
            Check(std::abs(field.divergence)<2e-12,"sample divergence is nonzero");
            if(count++){out<<',';}out<<"{\"s\":"<<s<<",\"r\":"<<r<<",\"field\":";FieldRecord(out,field);
            out<<",\"cartesian_fd_h\":"<<step<<",\"G_fd\":";Write(out,fd);out<<",\"maximum_fd_error\":"<<error<<'}';
        }
        const Vec q{1.3,0,map.height};const auto field=Evaluate(map,q);const auto area=map.MatchingArea(q);
        const double mass=Dot(field.velocity,area)/field.metric.baseArea,naive=field.metric.J*field.normalVelocity;
        const double heightSlope=.06*std::cos(.3*q[0]),formula=field.metric.J*field.normalVelocity-heightSlope*field.tangentialVelocity;
        Check(mass<0 && naive>0,"matching-height slope omission did not reverse the mass-flux sign");
        Check(std::abs(mass-formula)<1e-12,"matching mass-flux identity failed");
        const double c=1/map.radius,closed=2-2/c+std::log1p(2*c)/(c*c);
        Check(std::abs(pressure-closed)<1e-10,"normal-pressure curvature increment failed");
        out<<"],\"matching_sample\":{\"s\":1.3,\"field\":";FieldRecord(out,field);out<<",\"area\":";Write(out,area);
        out<<",\"height_slope\":"<<heightSlope<<",\"m_H\":"<<mass<<",\"m_H_formula\":"<<formula<<",\"naive_Jw\":"<<naive
           <<"},\"pressure_closed_form\":"<<closed<<",\"checks\":"<<checks
           <<",\"scope\":\"Revised-note geometry, closed flux/source budgets and normal acceleration; no CFD timestep or inner epsilon transport equation\"}\n";
        std::ofstream file(argv[2]);file<<out.str();Check(static_cast<bool>(file),"cannot write note-probe report");
        std::cout<<"Curved note checks "<<checks<<'\n';return 0;
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 1;}
}
