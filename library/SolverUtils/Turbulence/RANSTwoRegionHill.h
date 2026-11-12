// SPDX-License-Identifier: MIT
#ifndef NEKTAR_SOLVERUTILS_RANS_TWO_REGION_HILL_H
#define NEKTAR_SOLVERUTILS_RANS_TWO_REGION_HILL_H
#include <SolverUtils/Turbulence/RANSWallNormalMap.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <locale>
#include <queue>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace Nektar::SolverUtils
{
/** Certified two-region geometry kernel. Coupled use is restricted to the
 * explicit CompositeHillV1 finite/open budget-AWF pilot; kernel verification
 * alone does not qualify mapped physical evolution.
 * s is base horizontal position (not arclength); r is physical AWF normal
 * distance. Above Gamma_H, eta is a vertical extrusion coordinate, not wall
 * distance. Incoming-frame height is a separate physical reference query.
 * Cosine feet are explicit C1/one-sided-curvature interfaces.
 */
class RANSTwoRegionHill
{
public:
    struct OuterSample
    {
        RANSVector point{};
        RANSTensor jacobian{},inverse{};
        double determinant=0;
    };
    struct WallDistance
    {
        double distance=0,errorBound=0;
        std::vector<RANSVector> footpoints;
        int boxes=0;
    };
    struct Location { bool inner; RANSVector chart; };
    RANSTwoRegionHill(double hillHeight,double halfLength,double wallHeight,
        double roof,double xmin,double xmax)
        : m_Hh(hillHeight),m_L(halfLength),m_Hw(wallHeight),m_Z(roof),m_xmin(xmin),m_xmax(xmax),
          m_inner(Map(hillHeight,halfLength,wallHeight))
    {
        for(double x:{m_Hh,m_L,m_Hw,m_Z,m_xmin,m_xmax})Check(std::isfinite(x),"nonfinite geometry");
        Check(m_Hh>=0 && m_L>0 && m_Hw>0 && m_Z>m_Hh+m_Hw && m_xmin<-m_L && m_xmax>m_L,"invalid two-region extent");
        Check(CurvatureBound()*m_Hw<1,"normal reach exceeded");
        m_inner.Validate({m_xmin,0,0},{m_xmax,0,m_Hw},1-CurvatureBound()*m_Hw-1e-12);
        // x_H'=1-Hw*f''/(1+f'^2)^(3/2) >= 1-Hw*C > 0.
        // z_H<Z; thus outer J=x_H'*(Z-z_H)>0. x_H is strictly
        // monotone, and vertical intervals are unique: global injectivity.
        // The certified inner tube lies below the same graph Gamma_H and
        // outer columns above it. They share only their common boundary.
    }
    double CurvatureBound() const { return m_Hh*std::pow(std::acos(-1.)/m_L,2)/2; }
    RANSWallNormalMap::Sample Inner(double s,double r,int footSide=0) const
    {
        Check(s>=m_xmin && s<=m_xmax && r>=0 && r<=m_Hw,"inner point outside chart");
        return m_inner.At({s,0,r},footSide);
    }
    OuterSample Outer(double s,double eta,int footSide=0) const
    {
        Check(s>=m_xmin && s<=m_xmax && eta>=0 && eta<=1,"outer point outside chart");
        const auto h=Inner(s,m_Hw,footSide);OuterSample q;
        const double roof=m_Z-h.point[2],dx=h.jacobian[0][0],dz=(1-eta)*h.jacobian[2][0];
        q.point={h.point[0],0,h.point[2]+eta*roof};
        q.jacobian[0][0]=dx;q.jacobian[1][1]=1;q.jacobian[2][0]=dz;q.jacobian[2][2]=roof;
        q.determinant=dx*roof;
        q.inverse[0][0]=1/dx;q.inverse[1][1]=1;q.inverse[2][0]=-dz/(dx*roof);q.inverse[2][2]=1/roof;
        Check(q.determinant>0,"nonpositive outer Jacobian");return q;
    }
    RANSVector OuterInverse(const RANSVector &p) const
    {
        Check(p[1]==0 && p[0]>=m_xmin && p[0]<=m_xmax,"outer physical point outside horizontal extent");
        double a=m_xmin,b=m_xmax;
        for(int j=0;j<80;++j){const double s=(a+b)/2;if(Inner(s,m_Hw).point[0]<p[0])a=s;else b=s;}
        const double s=(a+b)/2,z=Inner(s,m_Hw).point[2],eta=(p[2]-z)/(m_Z-z);
        Check(eta>=-1e-13 && eta<=1+1e-13,"physical point outside outer region");return {s,0,eta};
    }
    RANSVector InnerInverse(const RANSVector &p) const
    {
        const auto q=m_inner.Inverse(p);
        Check(q[0]>=m_xmin && q[0]<=m_xmax && q[2]>=-1e-12 && q[2]<=m_Hw+1e-12,"physical point outside inner strip");return q;
    }
    /** Observation coordinates are physical (paper vertical y -> solver z).
     * side=-1/+1 selects H-minus/H-plus at the interface; it never averages
     * a one-sided derivative. Tangential foot sides remain explicit in At.
     */
    Location Locate(const RANSVector &p,int side=0) const
    {
        Check(side>=-1 && side<=1 && p[1]==0 && p[0]>=m_xmin && p[0]<=m_xmax && p[2]<=m_Z,"invalid physical observation");
        Check(p[2]>=m_inner.At({p[0],0,0}).point[2],"physical observation lies in solid");
        double a=m_xmin,b=m_xmax;
        for(int j=0;j<80;++j){const double s=(a+b)/2;if(Inner(s,m_Hw).point[0]<p[0])a=s;else b=s;}
        const double s=(a+b)/2,z=Inner(s,m_Hw).point[2];
        if(std::abs(p[2]-z)<1e-13 && side!=0)
            return side<0?Location{true,InnerInverse(p)}:Location{false,OuterInverse(p)};
        if(p[2]<z)return {true,InnerInverse(p)};
        return {false,OuterInverse(p)};
    }
    RANSVector InterfaceArea(double s,bool inner,int footSide=0) const
    {
        const auto h=Inner(s,m_Hw,footSide);const double sign=inner?1.:-1.;
        return {-sign*h.jacobian[2][0],0,sign*h.jacobian[0][0]};
    }
    /** Global distance to the complete physical floor/hill, including the
     * flat rays outside the feet. Interval bounding boxes enclose the exact
     * cosine graph (split at both feet and crest). Branch-and-bound supplies
     * an explicit distance error, and keeps multiple minimizing witnesses.
     * This is independent of the outer eta and the incoming reference height.
     */
    WallDistance PhysicalWallDistance(const RANSVector &p,double tolerance=1e-10,bool chordEnvelope=true) const
    {
        Check(p[1]==0 && std::isfinite(p[0]) && std::isfinite(p[2]) && tolerance>0,"invalid wall-distance query");
        struct Box {double a,b,lower;bool operator<(const Box &q)const{return lower>q.lower;}};
        auto f=[&](double s){return m_inner.At({s,0,0}).point[2];};
        auto d2=[&](double s){return std::pow(s-p[0],2)+std::pow(f(s)-p[2],2);};
        double best=std::numeric_limits<double>::infinity();
        std::array<double,2> bestSide{best,best},witnesses{};
        auto sample=[&](double s)
        {
            const double d=d2(s);best=std::min(best,d);
            for(int side=0;side<2;++side)if(std::abs(s)<=m_L && (side==0?s<=0:s>=0) && d<bestSide[side])
            {bestSide[side]=d;witnesses[side]=s;}
        };
        sample(std::min(p[0],-m_L));sample(std::max(p[0],m_L));sample(0);
        auto box=[&](double a,double b)
        {
            const double fa=f(a),fb=f(b),lo=std::min(fa,fb),hi=std::max(fa,fb);
            const double dx=std::max({a-p[0],p[0]-b,0.}),dz=std::max({lo-p[2],p[2]-hi,0.});
            if(!chordEnvelope)return Box{a,b,dx*dx+dz*dz}; // retained diagnostic reference
            // A rectangle loses correlation between x and z and gives a
            // first-order distance gap on sloping segments. The cosine has
            // |f''|<=CurvatureBound(); its graph is within C*(b-a)^2/8 of
            // the chord. Distance to that segment minus the envelope is a
            // second independent, certified lower bound. No normal-footpoint
            // uniqueness or Newton convergence is assumed.
            const double vx=b-a,vz=fb-fa,length2=vx*vx+vz*vz;
            const double t=std::clamp(((p[0]-a)*vx+(p[2]-fa)*vz)/length2,0.,1.);
            const double chord=std::hypot(a+t*vx-p[0],fa+t*vz-p[2]);
            const double envelope=CurvatureBound()*(b-a)*(b-a)/8;
            const double correlated=std::max(0.,chord-envelope);
            const double floating=64*std::numeric_limits<double>::epsilon()*(1+std::abs(p[0])+std::abs(p[2])+m_Hh+m_L);
            const double lower=std::max(0.,std::sqrt(std::max(dx*dx+dz*dz,correlated*correlated))-floating);
            return Box{a,b,lower*lower};
        };
        std::priority_queue<Box> queue;queue.push(box(-m_L,0));queue.push(box(0,m_L));
        WallDistance out;
        while(!queue.empty() && std::sqrt(best)-std::sqrt(queue.top().lower)>tolerance)
        {
            Check(++out.boxes<=200000,"wall-distance interval budget exhausted");
            const auto q=queue.top();queue.pop();const double c=(q.a+q.b)/2;sample(q.a);sample(c);sample(q.b);
            for(const auto &child:{box(q.a,c),box(c,q.b)})if(child.lower<=best)queue.push(child);
        }
        out.distance=std::sqrt(best);out.errorBound=queue.empty()?0:std::max(0.,out.distance-std::sqrt(queue.top().lower));
        // Recover both sides of a symmetric tied minimum without requiring a
        // unique normal footpoint. Only distance, never an averaged normal,
        // is consumed by outer Yap. At most one witness per separated basin.
        for(double s:std::vector<double>{std::min(p[0],-m_L),std::max(p[0],m_L),witnesses[0],witnesses[1]})if(std::sqrt(d2(s))-out.distance<=2*tolerance)
            if(std::none_of(out.footpoints.begin(),out.footpoints.end(),[&](const auto &q){return std::abs(q[0]-s)<10*tolerance;}))out.footpoints.push_back({s,0,f(s)});
        return out;
    }
    void RequireCompatible(const std::string &saved) const
    {Check(saved==Signature(),"geometry identity differs from saved provenance");}
    std::string Signature() const
    {
        std::ostringstream s;s.imbue(std::locale::classic());s<<std::setprecision(17)<<"CosineNormalVerticalV1:physical-nearest-wall-distance:foot-side-explicit";
        for(double x:{m_Hh,m_L,m_Hw,m_Z,m_xmin,m_xmax})s<<':'<<x;
        return s.str();
    }
private:
    static RANSWallNormalMap Map(double H,double L,double Hw)
    {RANSWallNormalMap m;m.kind=RANSWallNormalMap::Kind::CosineHill;m.amplitude=H;m.lengthX=L;m.height=Hw;return m;}
    static void Check(bool yes,const char *s){if(!yes)throw std::invalid_argument(std::string("Two-region hill: ")+s);}
    const double m_Hh,m_L,m_Hw,m_Z,m_xmin,m_xmax;
    const RANSWallNormalMap m_inner;
};
}
#endif
