// SPDX-License-Identifier: MIT
#include <SolverUtils/Turbulence/RANSWallNormalMap.h>
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <SolverUtils/Turbulence/WJEARSM.h>
#include <SolverUtils/Turbulence/RANSWallInnerStress.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace Nektar::SolverUtils;
namespace
{
void Check(bool b,const char *s){if(!b)throw std::runtime_error(s);}
template<class T> void Array(std::ostream &o,const T &v)
{o<<'[';bool first=true;for(const auto &x:v){if(!first)o<<',';first=false;o<<x;}o<<']';}
void Matrix(std::ostream &o,const RANSTensor &m)
{o<<'[';for(int i=0;i<3;++i){if(i)o<<',';Array(o,m[i]);}o<<']';}
double Dot(const RANSVector &a,const RANSVector &b){double s=0;for(int i=0;i<3;++i)s+=a[i]*b[i];return s;}
std::vector<std::array<double,2>> Gauss(int n)
{
    std::vector<std::array<double,2>> rule;
    for(int i=0;i<n;++i)
    {
        double x=std::cos(std::acos(-1.)*(i+.75)/(n+.5)),d=0;
        for(int j=0;j<100;++j)
        {double a=1,b=x;for(int k=2;k<=n;++k){double c=((2*k-1)*x*b-(k-1)*a)/k;a=b;b=c;}
         d=n*(x*b-a)/(x*x-1);const double dx=b/d;x-=dx;if(std::abs(dx)<1e-15)break;}
        rule.push_back({x,2/((1-x*x)*d*d)});
    }return rule;
}
}
int main(int argc,char **argv)
{
    try
    {
        Check(argc==2,"usage: WJTerrainGeometryProbe OUTPUT.json");
        const double nu=.051*.005/3.97,us=.051;
        const auto model=RANSNeutralReference::WJ();const double K=model.Q()*us*us;
        const auto wall=WJCalibratedAWFParameters(model);
        std::ofstream out(argv[1]);out<<std::setprecision(17)<<"{\"schema\":\"WJH1GeometryV1\",\"maps\":[";
        for(int kind=0;kind<3;++kind)
        {
            if(kind)out<<',';
            RANSWallNormalMap m;m.compositeTunnel=true;m.height=.375;m.displacement=.0875;
            m.tunnelTop=20;m.lengthX=2.5;m.amplitude=kind==0?0:1;
            m.kind=kind==2?RANSWallNormalMap::Kind::Plane:RANSWallNormalMap::Kind::CosineHill;
            m.slopeX=kind==2?.2:0;
            m.Validate({-5,0,0},{5,0,20},.1);
            out<<"{\"kind\":"<<kind<<",\"signature\":\""<<m.Signature()<<"\",\"samples\":[";
            bool first=true;double inverseError=0;
            for(double x:{-5.,-2.5,-1.9,-1.25,0.,1.25,1.9,2.5,5.})
            for(double r:{0.,.1,.375,1.,6.25,20.})for(int region:{-1,1})
            {
                if(r!=m.height && region==1)continue;
                for(int foot:{-1,1})
                {
                    const bool isFoot=m.kind==RANSWallNormalMap::Kind::CosineHill && std::abs(x)==m.lengthX;
                    if(!isFoot && foot==1)continue;
                    const RANSVector q{x,0,r};const auto g=m.At(q,isFoot?foot:0,region);
                    const auto inv=m.Inverse(g.point);
                    for(int a:{0,2})inverseError=std::max(inverseError,std::abs(inv[a]-q[a]));
                    if(!first)out<<',';first=false;
                    out<<"{\"q\":";Array(out,q);out<<",\"point\":";Array(out,g.point);
                    out<<",\"foot\":"<<(isFoot?foot:0)<<",\"region\":"<<region<<",\"jacobian\":";Matrix(out,g.jacobian);
                    out<<",\"inverse\":";Matrix(out,g.inverse);out<<",\"det\":"<<g.determinant
                       <<",\"base_area\":"<<g.baseArea<<",\"kappa_eff\":"<<g.kappa1+g.kappa2
                       <<",\"normal\":";Array(out,g.normal);out<<",\"tangent\":";Array(out,g.tangent1);
                    out<<",\"normal_region\":"<<(g.normalRegion?"true":"false")<<",\"eta\":"<<g.eta<<'}';
                }
            }
            out<<"],\"inverse_error\":"<<inverseError<<",\"integrals\":[";
            for(int order:{16,32})
            {
                if(order==32)out<<',';const auto rule=Gauss(order);
                double volume=0;RANSVector flux{},divx{},divz{};
                const std::vector<double> xs{-5,-2.5,0,2.5,5},rs{0,.375,2,6.25,20};
                for(int i=0;i<4;++i)for(int j=0;j<4;++j)
                    for(const auto &a:rule)for(const auto &b:rule)
                    {const RANSVector q{(xs[i]+xs[i+1])/2+a[0]*(xs[i+1]-xs[i])/2,0,
                        (rs[j]+rs[j+1])/2+b[0]*(rs[j+1]-rs[j])/2};
                     volume+=a[1]*b[1]*(xs[i+1]-xs[i])*(rs[j+1]-rs[j])/4*m.At(q).determinant;}
                for(int axis:{0,2})for(int side:{-1,1})for(int panel=0;panel<4;++panel)for(const auto &a:rule)
                {
                    const auto &knots=axis==0?rs:xs;
                    RANSVector q{};q[axis]=axis==0?(side<0?-5:5):(side<0?0:20);
                    q[2-axis]=(knots[panel]+knots[panel+1])/2+a[0]*(knots[panel+1]-knots[panel])/2;
                    const auto g=m.At(q);const auto c=m.Cofactor(q,axis,side);
                    const double w=a[1]*(knots[panel+1]-knots[panel])/2;
                    for(int d=0;d<3;++d){flux[d]+=w*c[d];divx[d]+=w*g.point[0]*c[d];divz[d]+=w*g.point[2]*c[d];}
                }
                out<<"{\"order\":"<<order<<",\"volume\":"<<volume<<",\"boundary_constant_flux\":";Array(out,flux);
                out<<",\"boundary_x_flux\":";Array(out,divx);out<<",\"boundary_z_flux\":";Array(out,divz);out<<'}';
            }
            out<<"]}";
        }
        out<<"],\"calibration\":[";
        for(double H:{.25,.375,.5})
        {
            if(H!=.25)out<<',';
            const double a=model.Kappa()*us,ell=H/(H/.005-1),zv=nu/a-ell;
            double lo=0,hi=H*(1-1e-12);
            for(int i=0;i<100;++i){double h=(lo+hi)/2;
                AnalyticalWallFunction awf({nu,H,h,H},K,wall);
                if(awf.ViscosityInterceptHeight()>zv)lo=h;else hi=h;}
            const double h=(lo+hi)/2;AWFProfileOptions options;options.formulation=AWFProfileFormulation::WJCalibratedBudgetAWF;
            options.quadratureOrder=32;
            auto profile=EquilibriumAWFProfiles::WJ({nu,H,h,H,H},K,model,options);
            double minK=K;for(int j=0;j<=100;++j)minK=std::min(minK,profile.ReferenceK(H*j/100));
            const auto &r=profile.Reference();const auto &awf=profile.Original();
            out<<"{\"H\":"<<H<<",\"h_s\":"<<h<<",\"z_v\":"<<awf.ViscosityInterceptHeight()
               <<",\"cutoff\":"<<r.cutoff<<",\"gate\":"<<r.productionBudget/r.minimumBudget
               <<",\"minimum_reference_k\":"<<minK<<",\"K\":"<<K<<",\"UH\":"<<r.matchingVelocity
               <<",\"nu_over_nut_H\":"<<nu/awf.EddyViscosity(H)<<'}';
        }
        int rejects=0;
        for(int kind=0;kind<5;++kind)
        {
            RANSWallNormalMap m;m.kind=RANSWallNormalMap::Kind::CosineHill;m.amplitude=1;m.lengthX=2.5;m.height=.375;
            m.compositeTunnel=kind!=0;m.displacement=.0875;m.tunnelTop=20;
            if(kind==1)m.height=2;if(kind==2)m.tunnelTop=1;
            if(kind==3)m.displacement=-1;if(kind==4)m.heightSlopeX=.1;
            try{m.Validate({-5,0,0},{5,0,20},.1);}catch(const std::invalid_argument &){++rejects;}
        }
        out<<"],\"invalid_geometry_rejections\":"<<rejects<<",\"frame_check\":{\"R_error\":";
        RANSWallNormalMap plane;plane.kind=RANSWallNormalMap::Kind::Plane;plane.slopeX=.2;
        const auto frame=plane.At({0,0,.375});RANSTensor G{},rotation{};
        const double e=us*us*us/(model.Kappa()*.375),shear=us/(model.Kappa()*.375);
        for(int i=0;i<3;++i){rotation[i][0]=frame.tangent1[i];rotation[i][1]=frame.tangent2[i];rotation[i][2]=frame.normal[i];
            for(int j=0;j<3;++j)G[i][j]=shear*frame.tangent1[i]*frame.normal[j];}
        const auto law=EvaluateWJClosure(K,e,G);RANSTensor flatG{};flatG[0][2]=shear;
        const auto flat=EvaluateWJClosure(K,e,flatG);double error=0;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j){double v=0;
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)v+=rotation[i][a]*flat.reynoldsStress[a][b]*rotation[j][b];
            error=std::max(error,std::abs(v-law.reynoldsStress[i][j]));}
        const auto tr=WJOuterTraction(law.reynoldsStress,0,frame.normal);
        out<<error<<",\"tangential_outer_traction\":"<<Dot(tr,frame.tangent1)<<",\"production\":"<<law.production<<",\"epsilon\":"<<e<<"}}\n";
        Check(out.good(),"geometry output write failed");Check(rejects==5,"invalid geometry accepted");
        return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
