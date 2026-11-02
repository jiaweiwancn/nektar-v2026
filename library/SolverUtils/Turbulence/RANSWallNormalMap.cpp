///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSWallNormalMap.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>

namespace Nektar::SolverUtils
{
namespace
{
void Require(bool valid, const char *message)
{
    if (!valid) throw std::invalid_argument(std::string("RANS normal map: ")+message);
}
double Dot(const RANSVector &a, const RANSVector &b)
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
RANSVector Cross(const RANSVector &a, const RANSVector &b)
{ return {a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]}; }
double Norm(const RANSVector &a)
{ return std::hypot(a[0],a[1],a[2]); }
RANSVector Column(const RANSTensor &a, int j)
{ return {a[0][j],a[1][j],a[2][j]}; }
RANSVector Multiply(const RANSTensor &a, const RANSVector &b)
{ return {Dot(a[0],b),Dot(a[1],b),Dot(a[2],b)}; }
RANSTensor Multiply(const RANSTensor &a, const RANSTensor &b)
{
    RANSTensor c{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)c[i][j]=Dot(a[i],Column(b,j));
    return c;
}
bool Finite(const RANSVector &a)
{ return std::all_of(a.begin(),a.end(),[](double x){return std::isfinite(x);}); }
RANSTensor Reciprocal(const RANSTensor &a, double determinant)
{
    RANSTensor inverse{};
    for(int i=0;i<3;++i)
    {
        const auto row=Cross(Column(a,(i+1)%3),Column(a,(i+2)%3));
        for(int j=0;j<3;++j)inverse[i][j]=row[j]/determinant;
    }
    return inverse;
}
void Parameters(const RANSWallNormalMap &m)
{
    Require(m.dimension==2 || m.dimension==3,"dimension must be two or three");
    Require(m.kind==RANSWallNormalMap::Kind::Flat || m.kind==RANSWallNormalMap::Kind::Plane ||
            m.kind==RANSWallNormalMap::Kind::Cylinder || m.kind==RANSWallNormalMap::Kind::Gaussian ||
            m.kind==RANSWallNormalMap::Kind::CosineHill,"unknown analytic geometry");
    for(double x:{m.height,m.heightSlopeX,m.heightSlopeY,m.heightAmplitude,m.heightWaveX,
                  m.heightWaveY,m.slopeX,m.slopeY,m.radius,m.amplitude,m.lengthX,m.lengthY})
        Require(std::isfinite(x),"nonfinite geometry parameter");
    Require(m.height>0,"reference height must be positive");
    Require(m.dimension==3 || (m.slopeY==0 && m.heightSlopeY==0),
            "2D geometry cannot have spanwise plane/height slopes");
    if(m.kind==RANSWallNormalMap::Kind::Cylinder)Require(m.radius>0,"cylinder radius must be positive");
    if(m.kind==RANSWallNormalMap::Kind::Gaussian || m.kind==RANSWallNormalMap::Kind::CosineHill)
        Require(m.lengthX>0,"streamwise shape length must be positive");
    if(m.kind==RANSWallNormalMap::Kind::Gaussian && m.dimension==3)
        Require(m.lengthY>0,"spanwise Gaussian length must be positive");
}
struct GraphJet
{
    double f=0, x=0, y=0, xx=0, xy=0, yy=0;
};
GraphJet Graph(const RANSWallNormalMap &m, const RANSVector &q, int footSide=0)
{
    using Kind=RANSWallNormalMap::Kind;
    GraphJet f;
    Require(footSide>=-1 && footSide<=1,"invalid foot side");
    if(m.kind==Kind::Flat || m.kind==Kind::Plane)
    {
        f.f=m.slopeX*q[0]+m.slopeY*q[1];f.x=m.slopeX;f.y=m.slopeY;
        Require(footSide==0,"foot side is only defined at cosine-hill feet");
    }
    else if(m.kind==Kind::Gaussian)
    {
        const double u=q[0]/m.lengthX,v=m.dimension==3?q[1]/m.lengthY:0.;
        f.f=m.amplitude*std::exp(-u*u-v*v);
        f.x=-2*u*f.f/m.lengthX;f.xx=(4*u*u-2)*f.f/(m.lengthX*m.lengthX);
        if(m.dimension==3)
        {
            f.y=-2*v*f.f/m.lengthY;f.yy=(4*v*v-2)*f.f/(m.lengthY*m.lengthY);
            f.xy=4*u*v*f.f/(m.lengthX*m.lengthY);
        }
        Require(footSide==0,"foot side is only defined at cosine-hill feet");
    }
    else if(m.kind==Kind::CosineHill)
    {
        const double tolerance=32*std::numeric_limits<double>::epsilon()*std::max({1.,m.lengthX,std::abs(q[0])});
        const bool foot=std::abs(std::abs(q[0])-m.lengthX)<=tolerance;
        Require(footSide==0 || foot,"one-sided foot metric requested away from a cosine foot");
        const double a=std::acos(-1.)/m.lengthX;
        if(foot)
        {
            if((q[0]<0 && footSide>0) || (q[0]>0 && footSide<0))f.xx=.5*m.amplitude*a*a;
        }
        else if(std::abs(q[0])<m.lengthX)
        {
            f.f=.5*m.amplitude*(1+std::cos(a*q[0]));
            f.x=-.5*m.amplitude*a*std::sin(a*q[0]);f.xx=-.5*m.amplitude*a*a*std::cos(a*q[0]);
        }
    }
    else Require(false,"graph data requested for a cylinder");
    for(double x:{f.f,f.x,f.y,f.xx,f.xy,f.yy})Require(std::isfinite(x),"nonfinite graph derivative");
    return f;
}
struct HeightJet { double h,x,y; };
HeightJet Height(const RANSWallNormalMap &m, const RANSVector &q)
{
    const double a=m.heightWaveX*q[0],b=m.heightWaveY*q[1];
    HeightJet h{m.height+m.heightSlopeX*q[0]+m.heightSlopeY*q[1]+m.heightAmplitude*std::sin(a)*std::cos(b),
                m.heightSlopeX+m.heightAmplitude*m.heightWaveX*std::cos(a)*std::cos(b),
                m.heightSlopeY-m.heightAmplitude*m.heightWaveY*std::sin(a)*std::sin(b)};
    Require(std::isfinite(h.h) && h.h>0 && std::isfinite(h.x) && std::isfinite(h.y),
            "nonpositive height or nonfinite height derivative");
    return h;
}
bool FlatGraph(const RANSWallNormalMap &m)
{
    using Kind=RANSWallNormalMap::Kind;
    return m.kind==Kind::Flat || m.kind==Kind::Plane ||
        ((m.kind==Kind::Gaussian || m.kind==Kind::CosineHill) && m.amplitude==0);
}
}

RANSWallNormalMap::Sample RANSWallNormalMap::At(const RANSVector &chart,int footSide) const
{
    Parameters(*this);Require(Finite(chart),"nonfinite chart coordinate");
    Require(dimension==3 || chart[1]==0,"2D chart has a nonzero unused span coordinate");
    Require(chart[2]>=0,"negative normal chart coordinate");
    const auto h=Height(*this,chart);Sample s;s.chart=chart;s.height=h.h;
    s.scale=h.h/height;s.distance=chart[2]*s.scale;
    RANSVector base{},t1{},t2{};
    if(kind==Kind::Cylinder)
    {
        Require(footSide==0,"foot side is only defined at cosine-hill feet");
        const double angle=chart[0]/radius,c=std::cos(angle),v=std::sin(angle);
        base={radius*v,chart[1],radius*(c-1)};t1={c,0,-v};t2={0,1,0};s.normal={v,0,c};
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)s.B[i][j]=t1[i]*t1[j]/radius;
        s.baseArea=1;
    }
    else
    {
        const auto f=Graph(*this,chart,footSide);base={chart[0],chart[1],f.f};
        t1={1,0,f.x};t2={0,1,f.y};s.baseArea=std::hypot(1.,f.x,f.y);
        s.normal={-f.x/s.baseArea,-f.y/s.baseArea,1/s.baseArea};
        const double g=s.baseArea*s.baseArea;
        RANSVector r1{},r2{};
        for(int i=0;i<3;++i)
        {
            r1[i]=((1+f.y*f.y)*t1[i]-f.x*f.y*t2[i])/g;
            r2[i]=((1+f.x*f.x)*t2[i]-f.x*f.y*t1[i])/g;
        }
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            s.B[i][j]=-(f.xx*r1[i]*r1[j]+f.xy*(r1[i]*r2[j]+r2[i]*r1[j])+f.yy*r2[i]*r2[j])/s.baseArea;
    }
    const double length=Norm(t1);for(int i=0;i<3;++i)s.tangent1[i]=t1[i]/length;
    s.tangent2=Cross(s.normal,s.tangent1);
    double trace=0,square=0;
    for(int i=0;i<3;++i){trace+=s.B[i][i];for(int j=0;j<3;++j)square+=s.B[i][j]*s.B[j][i];}
    const double split=std::sqrt(std::max(0.,2*square-trace*trace));
    s.kappa1=(trace+split)/2;s.kappa2=(trace-split)/2;
    const double a1=1+s.distance*s.kappa1,a2=1+s.distance*s.kappa2;
    Require(std::isfinite(a1) && std::isfinite(a2) && a1>0 && a2>0,"folded normal offset");
    s.J=a1*a2;s.determinant=s.baseArea*s.J*s.scale;
    RANSTensor a{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)a[i][j]=(i==j?1.:0.)+s.distance*s.B[i][j];
    s.Br=Multiply(s.B,Reciprocal(a,s.J));
    const auto at1=Multiply(a,t1),at2=Multiply(a,t2);
    for(int i=0;i<3;++i)
    {
        s.point[i]=base[i]+s.distance*s.normal[i];
        s.jacobian[i][0]=at1[i]+s.normal[i]*chart[2]*h.x/height;
        s.jacobian[i][1]=at2[i]+s.normal[i]*chart[2]*h.y/height;
        s.jacobian[i][2]=s.scale*s.normal[i];
    }
    const double determinant=Dot(Column(s.jacobian,0),Cross(Column(s.jacobian,1),Column(s.jacobian,2)));
    Require(std::isfinite(s.determinant) && s.determinant>0 && std::isfinite(determinant) && determinant>0 &&
            std::abs(determinant-s.determinant)<=1e-10*std::max(determinant,s.determinant),
            "invalid or inconsistent physical mapping determinant");
    s.inverse=Reciprocal(s.jacobian,s.determinant);
    Require(Finite(s.point) && Finite(s.normal) && Finite(s.tangent1) && Finite(s.tangent2),"nonfinite physical mapping");
    for(int i=0;i<3;++i)Require(Finite(s.jacobian[i]) && Finite(s.inverse[i]) && Finite(s.B[i]) && Finite(s.Br[i]),
                              "nonfinite mapping metric");
    return s;
}

RANSVector RANSWallNormalMap::Inverse(const RANSVector &point) const
{
    Parameters(*this);Require(Finite(point),"nonfinite physical inverse query");
    Require(dimension==3 || point[1]==0,"2D physical inverse query has nonzero span coordinate");
    RANSVector q{point[0],point[1],0};double r=0;
    if(kind==Kind::Cylinder)
    {
        const double radial=std::hypot(point[0],point[2]+radius);
        Require(radial>0,"cylinder axis has no unique normal foot");
        q[0]=radius*std::atan2(point[0],point[2]+radius);r=radial-radius;
    }
    else if(kind==Kind::Flat || kind==Kind::Plane)
    {
        const double length=std::hypot(1.,slopeX,slopeY);
        const RANSVector n{-slopeX/length,-slopeY/length,1/length};
        r=(point[2]-slopeX*point[0]-slopeY*point[1])/length;
        q[0]-=r*n[0];q[1]-=r*n[1];
    }
    else
    {
        const double tolerance=128*std::numeric_limits<double>::epsilon()*std::max(1.,Norm(point));
        bool converged=false;
        for(int iteration=0;iteration<100;++iteration)
        {
            const auto f=Graph(*this,q);const double z=f.f-point[2];
            const double gx=q[0]-point[0]+z*f.x,gy=dimension==3?q[1]-point[1]+z*f.y:0.;
            if(std::hypot(gx,gy)<=tolerance){converged=true;break;}
            const double a=1+f.x*f.x+z*f.xx,b=f.x*f.y+z*f.xy,c=1+f.y*f.y+z*f.yy;
            double dx=-gx,dy=-gy;
            if(dimension==2 && a>0)dx=-gx/a;
            else if(dimension==3 && a>0 && a*c-b*b>0)
            {const double determinant=a*c-b*b;dx=(-c*gx+b*gy)/determinant;dy=(b*gx-a*gy)/determinant;}
            double descent=gx*dx+gy*dy;
            if(!std::isfinite(descent) || descent>=0){dx=-gx;dy=-gy;descent=-gx*gx-gy*gy;}
            const double objective=.5*((q[0]-point[0])*(q[0]-point[0])+
                (q[1]-point[1])*(q[1]-point[1])+z*z);
            bool accepted=false;
            for(double fraction=1;fraction>=std::ldexp(1.,-40);fraction*=.5)
            {
                auto trial=q;trial[0]+=fraction*dx;trial[1]+=fraction*dy;
                const auto value=Graph(*this,trial);const double dz=value.f-point[2];
                const double next=.5*((trial[0]-point[0])*(trial[0]-point[0])+
                    (trial[1]-point[1])*(trial[1]-point[1])+dz*dz);
                const double trialGx=trial[0]-point[0]+dz*value.x,
                    trialGy=dimension==3?trial[1]-point[1]+dz*value.y:0.;
                // Near the foot, a useful Newton correction can change the
                // squared distance by less than one ulp. Still require the
                // normal-projection equation to improve in that case.
                const bool roundoffDecrease=std::abs(next-objective)<=
                    32*std::numeric_limits<double>::epsilon()*std::max(1.,objective) &&
                    std::hypot(trialGx,trialGy)<std::hypot(gx,gy);
                if(next<=objective+1e-4*fraction*descent || roundoffDecrease)
                {q=trial;accepted=true;break;}
            }
            Require(accepted,"normal-foot inverse line search failed");
        }
        Require(converged,"normal-foot inverse did not converge");
        const auto f=Graph(*this,q);const double length=std::hypot(1.,f.x,f.y);
        r=(-(point[0]-q[0])*f.x-(point[1]-q[1])*f.y+point[2]-f.f)/length;
    }
    const double tolerance=1e-10*std::max(1.,Norm(point));
    Require(std::isfinite(r) && r>=-tolerance,"physical inverse query lies below the wall");
    q[2]=std::max(0.,r)*height/Height(*this,q).h;
    const auto recovered=At(q);RANSVector error{};
    for(int i=0;i<3;++i)error[i]=recovered.point[i]-point[i];
    Require(Norm(error)<=tolerance,"physical inverse query is not on the recovered normal ray");
    return q;
}

RANSVector RANSWallNormalMap::Cofactor(const RANSVector &chart,int axis,int side,int footSide) const
{
    Require(axis>=0 && axis<3 && !(dimension==2 && axis==1),"unused or invalid chart face axis");
    Require(side==-1 || side==1,"face side must be -1 or +1");
    const auto s=At(chart,footSide);RANSVector result{};
    for(int i=0;i<3;++i)result[i]=side*s.determinant*s.inverse[axis][i];
    return result;
}
RANSWallNormalMap::Sample RANSWallNormalMap::Matching(const RANSVector &chart) const
{ auto q=chart;q[2]=height;return At(q); }
RANSVector RANSWallNormalMap::MatchingArea(const RANSVector &chart,int footSide) const
{ auto q=chart;q[2]=height;return Cofactor(q,2,1,footSide); }
bool RANSWallNormalMap::PlanarConstantHeight() const
{
    Parameters(*this);
    return FlatGraph(*this) && heightSlopeX==0 && heightSlopeY==0 &&
        (heightAmplitude==0 || heightWaveX==0);
}
bool RANSWallNormalMap::Identity() const
{
    return PlanarConstantHeight() &&
        ((kind!=Kind::Flat && kind!=Kind::Plane) || (slopeX==0 && slopeY==0));
}
void RANSWallNormalMap::Validate(const RANSVector &lower,const RANSVector &upper,double minimumStretch) const
{
    Parameters(*this);Require(Finite(lower) && Finite(upper),"nonfinite chart domain bound");
    Require(std::isfinite(minimumStretch) && minimumStretch>0 && minimumStretch<1,"invalid minimum normal stretch");
    Require(upper[0]>lower[0] && lower[2]>=0 && upper[2]>lower[2],"empty or negative-normal chart domain");
    Require(dimension==3?upper[1]>lower[1]:(lower[1]==0 && upper[1]==0),"invalid span chart domain");
    const double x0=heightSlopeX*lower[0],x1=heightSlopeX*upper[0],
        y0=heightSlopeY*lower[1],y1=heightSlopeY*upper[1];
    const double hmin=height+std::min(x0,x1)+std::min(y0,y1)-std::abs(heightAmplitude),
        hmax=height+std::max(x0,x1)+std::max(y0,y1)+std::abs(heightAmplitude);
    Require(std::isfinite(hmin) && std::isfinite(hmax) && hmin>0,"height bound does not certify a positive whole-domain height");
    const double rmax=upper[2]*hmax/height;
    Require(std::isfinite(rmax),"nonfinite whole-domain normal reach");
    if(kind==Kind::Cylinder)
    {
        Require((upper[0]-lower[0])/radius<2*std::acos(-1.),"cylinder arc wraps or repeats a normal ray");
        // Outward rays have positive radius and separate polar angles for any
        // r>=0. No unnecessary r/radius<1 restriction is imposed here.
    }
    else
    {
        double bound=0;
        if(kind==Kind::Gaussian)
        {
            const double diagonal=2*std::abs(amplitude)/(lengthX*lengthX);
            bound=diagonal;
            if(dimension==3)
            {
                const double cross=2*std::abs(amplitude)/(std::exp(1.)*lengthX*lengthY);
                bound=std::max(diagonal,2*std::abs(amplitude)/(lengthY*lengthY))+cross;
            }
        }
        else if(kind==Kind::CosineHill)
            bound=std::abs(amplitude)*std::acos(-1.)*std::acos(-1.)/(2*lengthX*lengthX);
        // n(q).(X(p)-X(q)) <= bound*|X(p)-X(q)|^2/2. Consequently a
        // point at distance r<1/bound is strictly nearer its own foot than
        // any other foot. This C1,1 proof includes both cosine-foot joins.
        Require(std::isfinite(bound) && rmax*bound<1-minimumStretch,
                "whole-domain graph bound does not certify normal-ray injectivity/minimum stretch");
    }
    for(int k=0;k<2;++k)for(int j=0;j<(dimension==3?2:1);++j)for(int i=0;i<2;++i)
        At({i?upper[0]:lower[0],j?upper[1]:lower[1],k?upper[2]:lower[2]});
}
std::string RANSWallNormalMap::Signature() const
{
    Parameters(*this);std::ostringstream out;out.imbue(std::locale::classic());out<<std::setprecision(17);
    out<<"CLS-EPAWF-NORMAL-MAP-v1:"<<dimension<<':'<<static_cast<int>(kind);
    for(double value:{height,heightSlopeX,heightSlopeY,heightAmplitude,heightWaveX,heightWaveY,
                      slopeX,slopeY,radius,amplitude,lengthX,lengthY})out<<':'<<value;
    return out.str();
}
} // namespace Nektar::SolverUtils
