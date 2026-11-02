///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSWallNormalMap.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
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
using Map=RANSWallNormalMap;
using Vec=RANSVector;
using Mat=RANSTensor;
double Dot(const Vec &a,const Vec &b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
Vec Cross(const Vec &a,const Vec &b)
{return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};}
double Norm(const Vec &v){return std::hypot(v[0],v[1],v[2]);}
Vec Column(const Mat &a,int i){return {a[0][i],a[1][i],a[2][i]};}
double Maximum(const Vec &v){return std::max({std::abs(v[0]),std::abs(v[1]),std::abs(v[2])});}
double Difference(const Vec &a,const Vec &b)
{return std::max({std::abs(a[0]-b[0]),std::abs(a[1]-b[1]),std::abs(a[2]-b[2])});}
double Difference(const Mat &a,const Mat &b)
{return std::max({Difference(a[0],b[0]),Difference(a[1],b[1]),Difference(a[2],b[2])});}
void String(std::ostream &o,const std::string &value)
{
    static const char hex[]="0123456789abcdef";o<<'"';
    for(char c:value)
    {
        const auto b=static_cast<unsigned char>(c);
        if(c=='"' || c=='\\')o<<'\\'<<c;
        else if(b<32)o<<"\\u00"<<hex[b/16]<<hex[b%16];
        else o<<c;
    }
    o<<'"';
}
void Write(std::ostream &o,const Vec &v){o<<'['<<v[0]<<','<<v[1]<<','<<v[2]<<']';}
void Write(std::ostream &o,const Mat &m)
{o<<'[';for(int i=0;i<3;++i){if(i)o<<',';Write(o,m[i]);}o<<']';}
struct Report
{
    int checks=0,rejected=0;
    double invariantError=0,mapFDError=0,inverseFDError=0,piolaFDError=0,normalFDError=0;
    double omittedHeightSlopeControl=0;
    std::vector<std::string> cases,feet,controls;
    void Check(bool good,const std::string &message)
    {++checks;if(!good)throw std::runtime_error(message);}
    void Near(double a,double b,double tolerance,const std::string &message)
    {
        const double error=std::abs(a-b)/std::max({1.,std::abs(a),std::abs(b)});
        if(std::isfinite(error))invariantError=std::max(invariantError,error);
        Check(std::isfinite(error) && error<=tolerance,message);
    }
    template<class F> void Reject(const std::string &name,F action)
    {
        std::string error;bool caught=false;
        try{action();}catch(const std::exception &e){caught=true;error=e.what();}
        Check(caught,"negative control accepted: "+name);++rejected;
        std::ostringstream o;o<<"{\"name\":";String(o,name);o<<",\"error\":";String(o,error);o<<'}';
        controls.push_back(o.str());
    }
};
void Near(Report &r,const Vec &a,const Vec &b,double tolerance,const std::string &name)
{for(int i=0;i<3;++i)r.Near(a[i],b[i],tolerance,name);}
void Near(Report &r,const Mat &a,const Mat &b,double tolerance,const std::string &name)
{for(int i=0;i<3;++i)Near(r,a[i],b[i],tolerance,name);}
const char *KindName(Map::Kind kind)
{
    switch(kind)
    {
        case Map::Kind::Flat:return "Flat";
        case Map::Kind::Plane:return "Plane";
        case Map::Kind::Cylinder:return "Cylinder";
        case Map::Kind::Gaussian:return "Gaussian";
        case Map::Kind::CosineHill:return "CosineHill";
    }
    throw std::runtime_error("unknown probe map kind");
}
void Parameters(std::ostream &o,const Map &m)
{
    o<<"{\"dimension\":"<<m.dimension<<",\"kind\":";String(o,KindName(m.kind));
    o<<",\"kind_id\":"<<static_cast<int>(m.kind)<<",\"height\":"<<m.height
     <<",\"heightSlopeX\":"<<m.heightSlopeX<<",\"heightSlopeY\":"<<m.heightSlopeY
     <<",\"heightAmplitude\":"<<m.heightAmplitude<<",\"heightWaveX\":"<<m.heightWaveX
     <<",\"heightWaveY\":"<<m.heightWaveY<<",\"slopeX\":"<<m.slopeX<<",\"slopeY\":"<<m.slopeY
     <<",\"radius\":"<<m.radius<<",\"amplitude\":"<<m.amplitude
     <<",\"lengthX\":"<<m.lengthX<<",\"lengthY\":"<<m.lengthY<<'}';
}
void Sample(std::ostream &o,const Map::Sample &s)
{
    o<<"{\"chart\":";Write(o,s.chart);o<<",\"point\":";Write(o,s.point);
    o<<",\"normal\":";Write(o,s.normal);o<<",\"tangent1\":";Write(o,s.tangent1);
    o<<",\"tangent2\":";Write(o,s.tangent2);o<<",\"B\":";Write(o,s.B);o<<",\"Br\":";Write(o,s.Br);
    o<<",\"jacobian\":";Write(o,s.jacobian);o<<",\"inverse\":";Write(o,s.inverse);
    o<<",\"determinant\":"<<s.determinant<<",\"baseArea\":"<<s.baseArea<<",\"r\":"<<s.distance
     <<",\"H\":"<<s.height<<",\"scale\":"<<s.scale<<",\"J\":"<<s.J
     <<",\"kappas\":["<<s.kappa1<<','<<s.kappa2<<"]}";
}
// Closed analytic surface derivatives, independent of the map's metric output.
struct Surface
{
    Vec point{},normal{},t1{},t2{},dn1{},dn2{};
    double H=0,H1=0,H2=0,baseArea=0,oldCurvature=0;
};
Surface Analytic(const Map &m,const Vec &q)
{
    Surface a;double f=0,fx=0,fy=0,fxx=0,fxy=0,fyy=0;
    const double wx=m.heightWaveX*q[0],wy=m.heightWaveY*q[1];
    a.H=m.height+m.heightSlopeX*q[0]+m.heightSlopeY*q[1]+m.heightAmplitude*std::sin(wx)*std::cos(wy);
    a.H1=m.heightSlopeX+m.heightAmplitude*m.heightWaveX*std::cos(wx)*std::cos(wy);
    a.H2=m.heightSlopeY-m.heightAmplitude*m.heightWaveY*std::sin(wx)*std::sin(wy);
    if(m.kind==Map::Kind::Cylinder)
    {
        const double angle=q[0]/m.radius,c=std::cos(angle),s=std::sin(angle);
        a.point={m.radius*s,q[1],m.radius*c-m.radius};a.normal={s,0,c};
        a.t1={c,0,-s};a.t2={0,1,0};a.baseArea=1;
        for(int i=0;i<3;++i)a.dn1[i]=a.t1[i]/m.radius;
        return a;
    }
    if(m.kind==Map::Kind::Flat || m.kind==Map::Kind::Plane)
    {f=m.slopeX*q[0]+m.slopeY*q[1];fx=m.slopeX;fy=m.slopeY;}
    else if(m.kind==Map::Kind::Gaussian)
    {
        const double ix=1/(m.lengthX*m.lengthX),iy=m.dimension==3?1/(m.lengthY*m.lengthY):0;
        f=m.amplitude*std::exp(-q[0]*q[0]*ix-q[1]*q[1]*iy);
        fx=-2*q[0]*ix*f;fy=-2*q[1]*iy*f;
        fxx=(4*q[0]*q[0]*ix*ix-2*ix)*f;
        fxy=4*q[0]*q[1]*ix*iy*f;fyy=(4*q[1]*q[1]*iy*iy-2*iy)*f;
    }
    else if(std::abs(q[0])<m.lengthX)
    {
        const double p=std::acos(-1.)/m.lengthX;
        f=.5*m.amplitude*(1+std::cos(p*q[0]));fx=-.5*m.amplitude*p*std::sin(p*q[0]);
        fxx=-.5*m.amplitude*p*p*std::cos(p*q[0]);
    }
    a.point={q[0],q[1],f};a.t1={1,0,fx};a.t2={0,1,fy};a.baseArea=std::sqrt(1+fx*fx+fy*fy);
    a.normal={-fx/a.baseArea,-fy/a.baseArea,1/a.baseArea};
    const Vec h1{fxx,fxy,0},h2{fxy,fyy,0};
    for(int i=0;i<3;++i)
    {
        a.dn1[i]=-h1[i]/a.baseArea-a.normal[i]*(fx*fxx+fy*fxy)/(a.baseArea*a.baseArea);
        a.dn2[i]=-h2[i]/a.baseArea-a.normal[i]*(fx*fxy+fy*fyy)/(a.baseArea*a.baseArea);
    }
    if(m.dimension==2)a.oldCurvature=fxx/std::pow(1+fx*fx,1.5);
    return a;
}
void Invariants(Report &r,const Map::Sample &s)
{
    const auto &n=s.normal;const double tolerance=3e-12;
    r.Near(Dot(n,n),1,tolerance,"unit normal");r.Near(Dot(s.tangent1,s.tangent1),1,tolerance,"unit tangent1");
    r.Near(Dot(s.tangent2,s.tangent2),1,tolerance,"unit tangent2");
    r.Near(Dot(s.tangent1,s.tangent2),0,tolerance,"orthogonal tangents");
    r.Near(Dot(s.tangent1,n),0,tolerance,"tangent1-normal orthogonality");
    Near(r,Cross(s.tangent1,s.tangent2),n,tolerance,"oriented tangent frame");
    r.Check(s.determinant>0 && s.J>0,"positive map and normal Jacobians");
    r.Near(s.determinant,Dot(Column(s.jacobian,0),Cross(Column(s.jacobian,1),Column(s.jacobian,2))),
           tolerance,"determinant from independent triple product");
    r.Near(s.determinant,s.baseArea*s.J*s.scale,tolerance,"physical volume measure");
    for(int i=0;i<3;++i)
    {
        r.Near(Dot(s.B[i],n),0,tolerance,"B tangency");r.Near(Dot(s.Br[i],n),0,tolerance,"Br tangency");
        for(int j=0;j<3;++j)
        {
            r.Near(Dot(s.inverse[i],Column(s.jacobian,j)),i==j?1.:0.,tolerance,"reciprocal Jacobian");
            r.Near(s.B[i][j],s.B[j][i],tolerance,"symmetric B");
            r.Near(s.Br[i][j],s.Br[j][i],tolerance,"symmetric Br");
            double product=0;for(int k=0;k<3;++k)product+=s.Br[i][k]*((k==j?1.:0.)+s.distance*s.B[k][j]);
            r.Near(product,s.B[i][j],tolerance,"offset curvature identity");
        }
    }
}
void InvariantErrors(std::ostream &o,const Map::Sample &s,const Vec &inversePoint)
{
    double reciprocal=0;Vec bn{},brn{};
    for(int i=0;i<3;++i)
    {
        bn[i]=Dot(s.B[i],s.normal);brn[i]=Dot(s.Br[i],s.normal);
        for(int j=0;j<3;++j)reciprocal=std::max(reciprocal,
            std::abs(Dot(s.inverse[i],Column(s.jacobian,j))-(i==j?1.:0.)));
    }
    o<<"{\"inverse_roundtrip\":"<<Difference(inversePoint,s.chart)
     <<",\"reciprocal_jacobian\":"<<reciprocal<<",\"normal_norm\":"<<std::abs(Dot(s.normal,s.normal)-1)
     <<",\"B_normal\":"<<Maximum(bn)<<",\"Br_normal\":"<<Maximum(brn)
     <<",\"determinant_triple_product\":"<<std::abs(s.determinant-
        Dot(Column(s.jacobian,0),Cross(Column(s.jacobian,1),Column(s.jacobian,2))))<<'}';
}
void ClosedChecks(Report &r,std::ostream &o,const Map &m,const Map::Sample &s)
{
    const auto a=Analytic(m,s.chart);const double distance=s.chart[2]*a.H/m.height;
    Vec point{};Mat jacobian{},B{},Br{};
    for(int i=0;i<3;++i)
    {
        point[i]=a.point[i]+distance*a.normal[i];
        jacobian[i][0]=a.t1[i]+distance*a.dn1[i]+s.chart[2]*a.H1/m.height*a.normal[i];
        jacobian[i][1]=a.t2[i]+distance*a.dn2[i]+s.chart[2]*a.H2/m.height*a.normal[i];
        jacobian[i][2]=a.H/m.height*a.normal[i];
    }
    Near(r,s.point,point,3e-12,"closed analytic position");Near(r,s.jacobian,jacobian,3e-12,"closed analytic map derivative");
    Near(r,s.normal,a.normal,3e-12,"closed analytic normal");r.Near(s.baseArea,a.baseArea,3e-12,"base surface measure");
    if(m.kind==Map::Kind::Cylinder)
    {
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        {B[i][j]=a.t1[i]*a.t1[j]/m.radius;Br[i][j]=a.t1[i]*a.t1[j]/(m.radius+distance);}
        Near(r,s.B,B,3e-12,"outward circular B sign");Near(r,s.Br,Br,3e-12,"circular offset curvature");
        r.Near(s.J,1+distance/m.radius,3e-12,"circular normal Jacobian");
    }
    else if(m.kind==Map::Kind::Flat || m.kind==Map::Kind::Plane)
    {Near(r,s.B,B,3e-12,"plane curvature zero");Near(r,s.Br,Br,3e-12,"plane offset curvature zero");r.Near(s.J,1,3e-12,"plane normal Jacobian");}
    else if(m.dimension==2)
    {
        const double norm=Norm(a.t1);
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        {
            B[i][j]=-a.oldCurvature*a.t1[i]*a.t1[j]/(norm*norm);
            Br[i][j]=B[i][j]/(1-distance*a.oldCurvature);
        }
        Near(r,s.B,B,3e-12,"B is minus old graph curvature");Near(r,s.Br,Br,3e-12,"graph offset sign");
    }
    o<<"{\"point\":";Write(o,point);o<<",\"jacobian\":";Write(o,jacobian);
    o<<",\"old_graph_curvature\":"<<a.oldCurvature<<'}';
}
void FiniteDifferences(Report &r,std::ostream &o,const Map &m,const Map::Sample &s)
{
    o<<'[';int count=0;double firstJac=0,firstInverse=0,firstPiola=0;
    for(double h:{4e-5,2e-5,1e-5})
    {
        const double cartesianH=h/10;Mat jacobian{},inverse{},normalGradient{};Vec piola{};
        if(m.dimension==2){jacobian[1][1]=1;inverse[1][1]=1;}
        for(int axis=0;axis<3;++axis)
        {
            if(m.dimension==2 && axis==1)continue;
            auto plus=s.chart,minus=s.chart;plus[axis]+=h;minus[axis]-=h;
            const auto a=m.At(plus),b=m.At(minus);
            const auto ca=m.Cofactor(plus,axis),cb=m.Cofactor(minus,axis);
            for(int i=0;i<3;++i){jacobian[i][axis]=(a.point[i]-b.point[i])/(2*h);piola[i]+=(ca[i]-cb[i])/(2*h);}
            plus=s.point;minus=s.point;plus[axis]+=cartesianH;minus[axis]-=cartesianH;
            const auto ip=m.Inverse(plus),im=m.Inverse(minus);const auto np=m.At(ip).normal,nm=m.At(im).normal;
            for(int i=0;i<3;++i)
            {inverse[i][axis]=(ip[i]-im[i])/(2*cartesianH);normalGradient[i][axis]=(np[i]-nm[i])/(2*cartesianH);}
        }
        const double je=Difference(jacobian,s.jacobian),ie=Difference(inverse,s.inverse),pe=Maximum(piola),ne=Difference(normalGradient,s.Br);
        r.mapFDError=std::max(r.mapFDError,je);r.inverseFDError=std::max(r.inverseFDError,ie);
        r.piolaFDError=std::max(r.piolaFDError,pe);r.normalFDError=std::max(r.normalFDError,ne);
        r.Check(std::isfinite(je) && je<2e-7,"finite-difference chart Jacobian");
        r.Check(std::isfinite(ie) && ie<2e-7,"finite-difference Cartesian inverse");
        r.Check(std::isfinite(pe) && pe<2e-7,"finite-difference Piola divergence");
        r.Check(std::isfinite(ne) && ne<2e-7,"finite-difference physical normal gradient");
        if(count==0){firstJac=je;firstInverse=ie;firstPiola=pe;}
        if(count==2)
        {
            // A roundoff plateau is acceptable; no spurious exact-FD or
            // asymptotic-order claim is made once the errors are this small.
            r.Check(je<=std::max(2e-8,1.25*firstJac),"chart derivative refinement sequence");
            r.Check(ie<=std::max(1e-7,1.25*firstInverse),"inverse derivative refinement sequence");
            r.Check(pe<=std::max(2e-8,1.25*firstPiola),"Piola refinement sequence");
        }
        if(count++){o<<',';}
        o<<"{\"chart_h\":"<<h<<",\"cartesian_h\":"<<cartesianH<<",\"map_jacobian\":";Write(o,jacobian);
        o<<",\"inverse_jacobian\":";Write(o,inverse);o<<",\"normal_gradient\":";Write(o,normalGradient);
        o<<",\"cofactor_divergence\":";Write(o,piola);
        o<<",\"errors\":["<<je<<','<<ie<<','<<pe<<','<<ne<<"]}";
    }
    o<<']';
}
void MatchingFlux(Report &r,std::ostream &o,const Map &m,const Vec &chart)
{
    const auto a=Analytic(m,chart);const auto matching=m.Matching(chart);const auto area=m.MatchingArea(chart);
    Vec e1{},e2{};for(int i=0;i<3;++i){e1[i]=a.t1[i]+a.H*a.dn1[i];e2[i]=a.t2[i]+a.H*a.dn2[i];}
    const double g11=Dot(e1,e1),g12=Dot(e1,e2),g22=Dot(e2,e2),det=g11*g22-g12*g12;
    const Vec U{.7,m.dimension==3?-.2:0.,.4};Vec reciprocal1{},reciprocal2{},qH{},v{},gradH{},full1{},full2{};
    const double b11=Dot(a.t1,a.t1),b12=Dot(a.t1,a.t2),b22=Dot(a.t2,a.t2),bdet=b11*b22-b12*b12;
    for(int i=0;i<3;++i)
    {
        reciprocal1[i]=(g22*e1[i]-g12*e2[i])/det;reciprocal2[i]=(g11*e2[i]-g12*e1[i])/det;
        qH[i]=a.H1*reciprocal1[i]+a.H2*reciprocal2[i];
        gradH[i]=((b22*a.H1-b12*a.H2)*a.t1[i]+(b11*a.H2-b12*a.H1)*a.t2[i])/bdet;
        full1[i]=e1[i]+a.H1*a.normal[i];full2[i]=e2[i]+a.H2*a.normal[i];
    }
    const double c1=Dot(reciprocal1,U),c2=Dot(reciprocal2,U),w=Dot(U,a.normal);
    for(int i=0;i<3;++i)v[i]=c1*a.t1[i]+c2*a.t2[i];
    const double J=Dot(Cross(e1,e2),a.normal)/a.baseArea;
    const double fromArea=Dot(U,area)/a.baseArea,fromFormula=J*(w-Dot(v,gradH)),naive=J*w;
    const auto expectedArea=Cross(full1,full2);const double beta=std::sqrt(1+Dot(qH,qH));Vec N{},expectedN{};
    for(int i=0;i<3;++i){N[i]=area[i]/Norm(area);expectedN[i]=(a.normal[i]-qH[i])/beta;}
    Near(r,area,expectedArea,4e-12,"matching area from independent tangent cross product");
    Near(r,N,expectedN,4e-12,"actual variable-height matching normal");
    r.Near(fromArea,fromFormula,4e-12,"physical matching mass flux");r.Near(J,matching.J,4e-12,"independent matching J");
    if(m.heightSlopeX!=0 || m.heightSlopeY!=0 || m.heightAmplitude!=0)
    {
        const double control=std::abs(fromFormula-naive);r.omittedHeightSlopeControl=std::max(r.omittedHeightSlopeControl,control);
        r.Check(control>1e-5,"height-slope omission must produce a nonzero control");
    }
    o<<"{\"sample\":";Sample(o,matching);o<<",\"area\":";Write(o,area);o<<",\"independent_area\":";Write(o,expectedArea);
    o<<",\"N_H\":";Write(o,N);o<<",\"naive_normal\":";Write(o,a.normal);o<<",\"q_H\":";Write(o,qH);
    o<<",\"grad_Gamma_H\":";Write(o,gradH);o<<",\"U\":";Write(o,U);o<<",\"v\":";Write(o,v);
    o<<",\"H_derivatives\":["<<a.H1<<','<<a.H2<<"],\"beta\":"<<beta<<",\"w\":"<<w
     <<",\"m_from_area\":"<<fromArea<<",\"m_from_formula\":"<<fromFormula<<",\"m_naive\":"<<naive<<'}';
}
struct Case{std::string name;Map map;Vec lower,upper;};
std::vector<Case> Cases()
{
    std::vector<Case> cases;
    auto add=[&](const std::string &name,Map m)
    {cases.push_back({name,m,{-1.3,m.dimension==3?-.9:0.,0},{1.3,m.dimension==3?.9:0.,4}});};
    for(int d:{2,3})
    {
        Map m;m.dimension=d;add("flat_"+std::to_string(d),m);
        m.kind=Map::Kind::Plane;m.slopeX=.25;m.slopeY=d==3?-.15:0;add("tilted_plane_"+std::to_string(d),m);
        m=Map{};m.dimension=d;m.kind=Map::Kind::Cylinder;add("outward_circle_"+std::to_string(d),m);
        m.heightSlopeX=.04;m.heightSlopeY=d==3?-.02:0;m.heightAmplitude=.06;m.heightWaveX=.7;m.heightWaveY=1.1;
        add("variable_height_circle_"+std::to_string(d),m);
        m=Map{};m.dimension=d;m.kind=Map::Kind::Gaussian;m.lengthX=1.2;m.lengthY=.9;add("gaussian_"+std::to_string(d),m);
        if(d==3)
        {
            m.heightSlopeX=.04;m.heightSlopeY=-.03;m.heightAmplitude=.06;m.heightWaveX=.7;m.heightWaveY=1.1;
            add("variable_height_gaussian_3",m);
        }
        m=Map{};m.dimension=d;m.kind=Map::Kind::CosineHill;add("cosine_explicit_feet_"+std::to_string(d),m);
    }
    return cases;
}
void FootChecks(Report &r,const Case &c)
{
    const auto &m=c.map;
    for(double foot:{-m.lengthX,m.lengthX})
    {
        const Vec q{foot,m.dimension==3?.17:0.,1.3};const auto left=m.At(q,-1),right=m.At(q,1);
        const auto cl=m.Cofactor(q,0,1,-1),cr=m.Cofactor(q,0,1,1);
        Near(r,left.point,right.point,3e-12,"C1 foot position");Near(r,left.normal,right.normal,3e-12,"C1 foot normal");
        Near(r,cl,cr,3e-12,"shared lateral foot cofactor");r.Check(Difference(left.B,right.B)>1e-3,"foot curvature must have distinct one-sided values");
        const double curvature=-m.amplitude*std::acos(-1.)*std::acos(-1.)/(2*m.lengthX*m.lengthX);
        Mat expectedLeft{},expectedRight{};expectedLeft[0][0]=foot>0?curvature:0;
        expectedRight[0][0]=foot<0?curvature:0;
        Near(r,left.B,expectedLeft,3e-12,"left-limit foot curvature sign");
        Near(r,right.B,expectedRight,3e-12,"right-limit foot curvature sign");
        Invariants(r,left);Invariants(r,right);
        std::ostringstream o;o<<std::setprecision(17)<<"{\"case\":";String(o,c.name);o<<",\"chart\":";Write(o,q);
        o<<",\"left\":";Sample(o,left);o<<",\"right\":";Sample(o,right);
        o<<",\"left_lateral_cofactor\":";Write(o,cl);o<<",\"right_lateral_cofactor\":";Write(o,cr);
        o<<",\"B_jump\":"<<Difference(left.B,right.B)<<'}';r.feet.push_back(o.str());
    }
}
void RunCase(Report &r,const Case &c)
{
    const auto &m=c.map;m.Validate(c.lower,c.upper);r.Check(true,"whole-domain certificate accepted");
    const bool planar=m.kind==Map::Kind::Flat || m.kind==Map::Kind::Plane;
    const bool constant=m.heightSlopeX==0 && m.heightSlopeY==0 && m.heightAmplitude==0;
    r.Check(m.PlanarConstantHeight()==(planar && constant),"planar constant-height flag");
    r.Check(m.Identity()==(m.kind==Map::Kind::Flat && constant),"identity geometry flag");
    std::ostringstream o;o<<std::setprecision(17)<<"{\"name\":";String(o,c.name);o<<",\"parameters\":";Parameters(o,m);
    o<<",\"lower\":";Write(o,c.lower);o<<",\"upper\":";Write(o,c.upper);o<<",\"minimum_stretch\":0.2,\"signature\":";String(o,m.Signature());
    o<<",\"identity\":"<<(m.Identity()?"true":"false")<<",\"planar_constant_height\":"<<(m.PlanarConstantHeight()?"true":"false");
    if(m.kind==Map::Kind::CosineHill)o<<",\"aligned_patch_edges\":["<<c.lower[0]<<','<<-m.lengthX<<','<<m.lengthX<<','<<c.upper[0]<<']';
    o<<",\"samples\":[";int count=0;
    for(const Vec q:std::vector<Vec>{{-.67,m.dimension==3?.23:0.,.35},{.17,m.dimension==3?-.31:0.,1.1},{.73,m.dimension==3?.51:0.,3.1}})
    {
        const auto s=m.At(q);Invariants(r,s);const auto inverse=m.Inverse(s.point);
        Near(r,inverse,q,3e-10,"inverse roundtrip");if(count++){o<<',';}
        o<<"{\"map\":";Sample(o,s);
        o<<",\"inverse_point\":";Write(o,inverse);o<<",\"invariant_errors\":";InvariantErrors(o,s,inverse);
        o<<",\"cofactors\":[";
        for(int axis=0;axis<3;++axis)
        {
            if(axis)o<<',';
            if(m.dimension==2 && axis==1){o<<"null";continue;}
            const auto positive=m.Cofactor(q,axis),negative=m.Cofactor(q,axis,-1);Vec sum{};
            for(int i=0;i<3;++i){sum[i]=positive[i]+negative[i];}
            Near(r,sum,Vec{},3e-12,"opposite face orientations");
            Write(o,positive);
        }
        o<<"],\"closed\":";ClosedChecks(r,o,m,s);o<<",\"finite_differences\":";FiniteDifferences(r,o,m,s);
        o<<",\"matching\":";MatchingFlux(r,o,m,q);o<<'}';
    }
    o<<"]}";r.cases.push_back(o.str());if(m.kind==Map::Kind::CosineHill)FootChecks(r,c);
}
void NegativeControls(Report &r)
{
    const Vec lo{-1,0,0},hi{1,0,4};Map m;
    m.kind=Map::Kind::Gaussian;m.amplitude=.2;
    r.Reject("graph_tube_beyond_certified_reach",[&]{m.Validate(lo,hi);});
    m=Map{};m.kind=Map::Kind::CosineHill;
    r.Reject("actual_folded_graph_offset",[&]{m.At({-m.lengthX,0,12},1);});
    m=Map{};m.kind=Map::Kind::Cylinder;
    r.Reject("wrapping_cylinder_arc",[&]{m.Validate({-3.2*m.radius,0,0},{3.2*m.radius,0,4});});
    r.Reject("cylinder_axis_inverse",[&]{m.Inverse({0,0,-m.radius});});
    m.radius=0;r.Reject("nonpositive_cylinder_radius",[&]{m.Validate(lo,hi);});
    m=Map{};m.heightAmplitude=3;m.heightWaveX=std::acos(-1.)/2;
    const double corner1=Analytic(m,{-2,0,0}).H,corner2=Analytic(m,{2,0,0}).H,inside=Analytic(m,{-1,0,0}).H;
    r.Check(corner1>0 && corner2>0 && inside<0,"between-corner height control");
    r.Reject("nonpositive_height_between_positive_corners",[&]{m.Validate({-2,0,0},{2,0,4});});
    r.Reject("nonpositive_local_height",[&]{m.At({-1,0,.3});});
    m=Map{};m.slopeY=.1;r.Reject("2D_spanwise_plane_slope",[&]{m.Validate(lo,hi);});
    m=Map{};m.heightSlopeY=.1;r.Reject("2D_spanwise_height_slope",[&]{m.Validate(lo,hi);});
    m=Map{};m.amplitude=std::numeric_limits<double>::quiet_NaN();
    r.Reject("nonfinite_parameter",[&]{m.Validate(lo,hi);});
    m=Map{};m.height=0;r.Reject("nonpositive_reference_height",[&]{m.Validate(lo,hi);});
    m=Map{};r.Reject("negative_rho_domain",[&]{m.Validate({-1,0,-.1},hi);});
    r.Reject("empty_tangential_domain",[&]{m.Validate(lo,{-1,0,4});});
    r.Reject("empty_normal_domain",[&]{m.Validate(lo,{1,0,0});});
    r.Reject("unused_2D_span_domain",[&]{m.Validate({-1,-.1,0},{1,.1,4});});
    r.Reject("negative_point_rho",[&]{m.At({0,0,-.1});});
    r.Reject("nonzero_2D_chart_span",[&]{m.At({0,.1,1});});
    r.Reject("nonfinite_chart",[&]{m.At({std::numeric_limits<double>::infinity(),0,1});});
    r.Reject("unused_2D_face_axis",[&]{m.Cofactor({0,0,1},1);});
    r.Reject("invalid_face_side",[&]{m.Cofactor({0,0,1},0,0);});
    r.Reject("below_wall_inverse",[&]{m.Inverse({.1,0,-.2});});
    r.Reject("invalid_minimum_stretch",[&]{m.Validate(lo,hi,0);});
    r.Reject("foot_side_on_plane",[&]{m.At({0,0,1},1);});
    m.kind=Map::Kind::CosineHill;r.Reject("off_foot_one_sided_metric",[&]{m.At({0,0,1},1);});
    r.Reject("invalid_foot_side",[&]{m.At({m.lengthX,0,1},2);});
    m.kind=static_cast<Map::Kind>(99);r.Reject("unknown_geometry_kind",[&]{m.Validate(lo,hi);});
}
void Array(std::ostream &o,const std::vector<std::string> &items)
{o<<'[';bool first=true;for(const auto &item:items){if(!first)o<<',';first=false;o<<item;}o<<']';}
void Save(const std::string &filename,const Report &r,const std::string &error)
{
    std::ofstream o(filename);if(!o)throw std::runtime_error("cannot open geometry report: "+filename);
    o<<std::setprecision(17)<<"{\"schema\":1,\"success\":"<<(error.empty()?"true":"false")<<",\"checks\":"<<r.checks
     <<",\"rejections\":"<<r.rejected<<",\"maximum_invariant_error\":"<<r.invariantError
     <<",\"maximum_map_fd_error\":"<<r.mapFDError<<",\"maximum_inverse_fd_error\":"<<r.inverseFDError
     <<",\"maximum_piola_fd_error\":"<<r.piolaFDError<<",\"maximum_normal_fd_error\":"<<r.normalFDError
     <<",\"height_slope_omission_control\":"<<r.omittedHeightSlopeControl<<",\"error\":";String(o,error);
    o<<",\"cases\":";Array(o,r.cases);o<<",\"foot_interfaces\":";Array(o,r.feet);
    o<<",\"negative_controls\":";Array(o,r.controls);o<<"}\n";
    if(!o)throw std::runtime_error("cannot write complete geometry report");
}
}
int main(int argc,char **argv)
{
    if(argc!=3 || std::string(argv[1])!="--output")
    {std::cerr<<"usage: RANSCurvedGeometryProbe --output REPORT.json\n";return 1;}
    Report report;std::string error;
    try
    {
        for(const auto &c:Cases())
        {
            try{RunCase(report,c);}
            catch(const std::exception &e){throw std::runtime_error(c.name+": "+e.what());}
        }
        NegativeControls(report);report.Check(report.omittedHeightSlopeControl>1e-4,"nonzero height-slope control recorded");
    }
    catch(const std::exception &e){error=e.what();}
    try{Save(argv[2],report,error);}catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
    if(!error.empty()){std::cerr<<error<<'\n';return 1;}
    std::cout<<"Curved geometry checks "<<report.checks<<", rejected controls "<<report.rejected<<'\n';
    return 0;
}
