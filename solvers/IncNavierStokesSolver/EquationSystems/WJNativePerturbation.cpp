// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativePerturbation.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace Nektar::WJNative
{
namespace
{
constexpr std::array<double,5> FixedX{-2.5,-1.25,0.,1.25,2.5};
constexpr double EdgeTolerance=1e-12;

void Require(bool value,const char *message)
{if(!value)throw std::runtime_error(message);}

bool Close(double a,double b)
{return std::abs(a-b)<=EdgeTolerance*std::max({1.,std::abs(a),std::abs(b)});}

std::vector<double> Edges(const std::vector<double> &input)
{
    Require(input.size()>=2,"perturbation needs actual native edge coordinates");
    for(double value:input)Require(std::isfinite(value),"nonfinite native perturbation edge");
    std::vector<double> result=input;
    std::sort(result.begin(),result.end());
    result.erase(std::unique(result.begin(),result.end(),
        [](double a,double b){return Close(a,b);}),result.end());
    Require(result.size()>=2,"degenerate native perturbation edge inventory");
    return result;
}

int EdgeIndex(const std::vector<double> &edges,double value)
{
    int found=-1;
    for(int i=0;i<int(edges.size());++i)if(Close(edges[i],value))
    {
        Require(found<0,"ambiguous native perturbation knot coordinate");found=i;
    }
    Require(found>=0,"perturbation knot is not an actual native mesh edge");
    return found;
}

void CheckShape(PerturbationShape shape)
{
    Require(shape==PerturbationShape::PrimaryQuartile||
            shape==PerturbationShape::SecondaryLowerSkew,
            "unsupported prospective perturbation shape");
}

std::array<double,5> SelectRadial(const std::vector<double> &candidates,
    const std::array<double,3> &targets,const std::array<double,5> *exclude)
{
    Require(candidates.size()>=5,"fewer than five admissible native radial knots");
    const double lower=candidates.front(),upper=candidates.back(),span=upper-lower;
    double best=std::numeric_limits<double>::infinity();
    std::array<double,5> selected{};
    bool found=false;
    // Select by physical coordinates, not by element index on a stretched mesh.
    // Equal costs retain the lexicographically first native triple.
    for(int i=1;i<int(candidates.size())-3;++i)
        for(int j=i+1;j<int(candidates.size())-2;++j)
            for(int k=j+1;k<int(candidates.size())-1;++k)
            {
                std::array<double,5> trial{lower,candidates[i],candidates[j],candidates[k],upper};
                if(exclude&&trial==*exclude)continue;
                double cost=0;
                for(int a=0;a<3;++a){double delta=(trial[a+1]-lower)/span-targets[a];cost+=delta*delta;}
                if(cost<best){best=cost;selected=trial;found=true;}
            }
    Require(found,"requested second shape has no distinct native radial knot triple");
    return selected;
}

double PeakEdge(const std::vector<double> &edges,double lower,double upper,double fraction)
{
    const double target=lower+fraction*(upper-lower);
    double best=std::numeric_limits<double>::infinity(),selected=0;
    bool found=false;
    for(double edge:edges)if(edge>lower+EdgeTolerance&&edge<upper-EdgeTolerance)
    {
        const double error=std::abs(edge-target);
        if(error<best){best=error;selected=edge;found=true;}
    }
    Require(found,"independent epsilon bump lacks an interior native peak edge");
    return selected;
}

struct OneJet{double value=0,first=0,second=0;};

double B(const std::array<double,5> &knots,int i,int degree,double x)
{
    if(degree==0)return x>=knots[i]&&x<knots[i+1]?1.:0.;
    return (x-knots[i])/(knots[i+degree]-knots[i])*B(knots,i,degree-1,x)+
           (knots[i+degree+1]-x)/(knots[i+degree+1]-knots[i+1])*B(knots,i+1,degree-1,x);
}

double DB(const std::array<double,5> &knots,int i,int degree,double x)
{
    if(degree==0)return 0.;
    return degree/(knots[i+degree]-knots[i])*B(knots,i,degree-1,x)-
           degree/(knots[i+degree+1]-knots[i+1])*B(knots,i+1,degree-1,x);
}

OneJet Spline(const std::array<double,5> &knots,double x)
{
    // A single cubic B-spline on five distinct knots has value, first and
    // second derivative zero at both exterior support endpoints.
    if(x<=knots.front()||x>=knots.back())return {};
    return {B(knots,0,3,x),DB(knots,0,3,x),
        3/(knots[3]-knots[0])*DB(knots,0,2,x)-
        3/(knots[4]-knots[1])*DB(knots,1,2,x)};
}

OneJet HermiteBump(double r,double lower,double peak,double upper)
{
    if(r<=lower||r>=upper)return {};
    const bool rising=r<=peak;
    const double width=rising?peak-lower:upper-peak;
    const double t=rising?(r-lower)/width:(upper-r)/width;
    const double direction=rising?1.:-1.;
    return {t*t*(3-2*t),direction*6*t*(1-t)/width,(6-12*t)/(width*width)};
}

bool SameLayout(const PerturbationLayout &a,const PerturbationLayout &b)
{
    // Normalization is bound to the exact frozen mathematical shape, not a
    // tolerance-matched replacement selected on a refined mesh.
    return a.xKnots==b.xKnots&&a.rKnots==b.rKnots&&a.H==b.H&&
        a.epsilonPeakRadius==b.epsilonPeakRadius&&a.shape==b.shape;
}
}

PerturbationLayout Perturbation::SelectLayout(const std::vector<double> &nativeXEdges,
    const std::vector<double> &nativeREdges,double H,PerturbationShape shape)
{
    CheckShape(shape);
    Require(std::isfinite(H)&&H>.075&&H<.45,"interface lies outside prospective perturbation window");
    const auto x=Edges(nativeXEdges),r=Edges(nativeREdges);
    for(double knot:FixedX)EdgeIndex(x,knot);
    EdgeIndex(r,H);
    std::vector<double> candidates;
    for(double edge:r)if(edge>=.075&&edge<=.45)candidates.push_back(edge);
    Require(candidates.size()>=5&&candidates.back()-candidates.front()>=.3&&
                candidates.front()<H&&candidates.back()>H,
            "native radial edges cannot support the pinned broad cross-interface perturbation");
    const auto primary=SelectRadial(candidates,{.25,.5,.75},nullptr);
    PerturbationLayout layout;layout.xKnots=FixedX;layout.H=H;layout.shape=shape;
    layout.rKnots=shape==PerturbationShape::PrimaryQuartile?primary:
        SelectRadial(candidates,{.2,.4,.65},&primary);
    layout.epsilonPeakRadius=PeakEdge(r,H,layout.rKnots.back(),
        shape==PerturbationShape::PrimaryQuartile?.5:.65);
    return layout;
}

Perturbation::Perturbation(const Reference &reference,const PerturbationLayout &layout,
    const std::vector<double> &nativeXEdges,const std::vector<double> &nativeREdges,
    PerturbationSigns signs):m_layout(layout),m_signs(signs)
{
    CheckShape(layout.shape);
    const auto &p=reference.Parameters();m_KIn=p.KIn;m_AIn=reference.Amplitude(p.KIn);m_d0=p.d0;
    Require((signs.velocity==1||signs.velocity==-1)&&
                (signs.k==1||signs.k==-1)&&(signs.epsilon==1||signs.epsilon==-1),
            "perturbation signs must each be +1 or -1");
    Require(layout.xKnots==FixedX&&layout.H==p.H,
            "frozen perturbation streamwise knots/interface differ from the pinned shape");
    for(int i=0;i<5;++i)
    {
        Require(std::isfinite(layout.rKnots[i])&&(i==0||layout.rKnots[i]>layout.rKnots[i-1]),
                "perturbation B3 radial knots are not finite and strictly increasing");
    }
    Require(layout.rKnots.front()>=.075&&layout.rKnots.back()<=.45&&
                layout.rKnots.back()-layout.rKnots.front()>=.3&&
                layout.rKnots.front()<p.H&&p.H<layout.epsilonPeakRadius&&
                layout.epsilonPeakRadius<layout.rKnots.back(),
            "frozen perturbation layout violates support/interface/outer-peak constraints");
    const auto x=Edges(nativeXEdges),r=Edges(nativeREdges);
    for(double knot:layout.xKnots)EdgeIndex(x,knot);
    for(double knot:layout.rKnots)EdgeIndex(r,knot);
    const int interface=EdgeIndex(r,p.H);EdgeIndex(r,layout.epsilonPeakRadius);
    m_xIntervals=EdgeIndex(x,layout.xKnots.back())-EdgeIndex(x,layout.xKnots.front());
    m_rIntervals=EdgeIndex(r,layout.rKnots.back())-EdgeIndex(r,layout.rKnots.front());
    m_outerIntervals=EdgeIndex(r,layout.rKnots.back())-interface;
    Require(m_xIntervals>=4&&m_rIntervals>=4&&m_outerIntervals>=4,
            "perturbation needs at least four native intervals in x, cross-H r and outer-bump r");
    Require(x.front()<layout.xKnots.front()&&x.back()>layout.xKnots.back()&&
                r.front()<layout.rKnots.front()&&r.back()>layout.rKnots.back(),
            "perturbation support is not interior to native exterior boundaries");
    Require(Close(x.front(),-10)&&Close(x.back(),15)&&Close(r.front(),0)&&Close(r.back(),p.topRadius),
            "perturbation domain differs from retained full-height physical application");
    m_domainArea=(x.back()-x.front())*(r.back()-r.front());
    m_outerDomainArea=(x.back()-x.front())*(r.back()-p.H);
    Require(SupportArea()>=1.5,"perturbation support area in physical code units is below 1.5");
}

double Perturbation::SupportArea()const
{return (m_layout.xKnots.back()-m_layout.xKnots.front())*(m_layout.rKnots.back()-m_layout.rKnots.front());}
double Perturbation::OuterSupportArea()const
{return (m_layout.xKnots.back()-m_layout.xKnots.front())*(m_layout.rKnots.back()-m_layout.H);}

void Perturbation::SetNormalization(const PerturbationNormalization &value)
{
    Require(!m_hasNormalization,"perturbation peak normalization is immutable after first load");
    Require(SameLayout(value.layout,m_layout),"normalization record belongs to a different perturbation layout");
    Require(value.evidenceSha256.size()==64&&std::all_of(value.evidenceSha256.begin(),value.evidenceSha256.end(),
        [](char c){return (c>='0'&&c<='9')||(c>='a'&&c<='f');}),
        "perturbation normalization requires an immutable evidence SHA256");
    Require(std::isfinite(value.unitPeakLower)&&std::isfinite(value.unitPeakEstimate)&&
                std::isfinite(value.unitPeakUpper)&&value.unitPeakLower>0&&
                value.unitPeakLower<=value.unitPeakEstimate&&value.unitPeakEstimate<=value.unitPeakUpper&&
                (value.unitPeakUpper-value.unitPeakLower)/value.unitPeakEstimate<=PeakBracketRelativeWidth,
            "unit perturbation global-peak bracket is missing, invalid or insufficiently resolved");
    m_normalization=value;m_hasNormalization=true;
}

const PerturbationNormalization &Perturbation::Normalization()const
{Require(m_hasNormalization,"perturbation velocity lacks independently qualified normalization");return m_normalization;}
double Perturbation::StreamfunctionCoefficient()const
{return m_signs.velocity*VelocityPeak/Normalization().unitPeakEstimate;}

double Perturbation::Radius(double x,double z)const
{
    Require(std::isfinite(x)&&std::isfinite(z),"nonfinite perturbation evaluation coordinate");
    return z-m_d0;
}

double Perturbation::UnitStreamfunction(double x,double z)const
{
    const double r=Radius(x,z);
    return Spline(m_layout.xKnots,x).value*Spline(m_layout.rKnots,r).value;
}
Vector Perturbation::UnitVelocity(double x,double z)const
{
    const double r=Radius(x,z);const auto a=Spline(m_layout.xKnots,x),b=Spline(m_layout.rKnots,r);
    return {a.value*b.first,0.,-a.first*b.value};
}
Tensor Perturbation::UnitGradient(double x,double z)const
{
    const double r=Radius(x,z);const auto a=Spline(m_layout.xKnots,x),b=Spline(m_layout.rKnots,r);
    Tensor gradient{};
    gradient[0][0]=a.first*b.first;gradient[0][2]=a.value*b.second;
    gradient[2][0]=-a.second*b.value;gradient[2][2]=-a.first*b.first;
    return gradient;
}
Vector Perturbation::Velocity(double x,double z)const
{
    Vector value=UnitVelocity(x,z);const double scale=StreamfunctionCoefficient();
    for(double &component:value)component*=scale;return value;
}
Tensor Perturbation::VelocityGradient(double x,double z)const
{
    Tensor value=UnitGradient(x,z);const double scale=StreamfunctionCoefficient();
    for(auto &row:value)for(double &component:row)component*=scale;return value;
}

PerturbationScalarJet Perturbation::KShape(double x,double z)const
{
    const double r=Radius(x,z);const auto a=Spline(m_layout.xKnots,x);
    const auto b=HermiteBump(r,m_layout.rKnots.front(),m_layout.H,m_layout.rKnots.back());
    // Uniform fixed x knots give max B3_x=B3_x(0)=2/3 exactly.
    return {1.5*a.value*b.value,1.5*a.first*b.value,1.5*a.value*b.first};
}
PerturbationScalarJet Perturbation::IndependentEpsilonShape(double x,double z)const
{
    const double r=Radius(x,z);const auto a=Spline(m_layout.xKnots,x);
    const auto b=HermiteBump(r,m_layout.H,m_layout.epsilonPeakRadius,m_layout.rKnots.back());
    return {1.5*a.value*b.value,1.5*a.first*b.value,1.5*a.value*b.first};
}
PerturbationScalarJet Perturbation::KIncrement(double x,double z)const
{
    auto value=KShape(x,z);const double scale=m_signs.k*KRelativePeak*m_KIn;
    value.value*=scale;value.dx*=scale;value.dr*=scale;return value;
}
PerturbationScalarJet Perturbation::IndependentEpsilonIncrement(double x,double z)const
{
    const double r=Radius(x,z);const auto shape=IndependentEpsilonShape(x,z);
    if(r<=m_layout.H||r>=m_layout.rKnots.back())return {};
    const double scale=m_signs.epsilon*EpsilonRelativePeak*m_AIn;
    return {scale*shape.value/r,scale*shape.dx/r,
            scale*(shape.dr/r-shape.value/(r*r))};
}
}
