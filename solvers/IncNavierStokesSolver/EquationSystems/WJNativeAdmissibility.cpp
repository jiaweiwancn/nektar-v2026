// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeAdmissibility.h>
#include <Eigen/Eigenvalues>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Nektar::WJNative
{
namespace
{
using Real=long double;
using Polynomial=std::vector<Real>;
void Require(bool value,const char *message)
{ if(!value)throw std::invalid_argument(message); }
Real Evaluate(const Polynomial &p,Real x)
{ Real value=0;for(auto i=p.rbegin();i!=p.rend();++i)value=value*x+*i;return value; }

/** All real roots in a bounded chart. Recursive derivative roots partition
 * the interval into monotone panels, including repeated roots. Only EXACT
 * zero leading coefficients are removed. No near-degree cutoff can discard
 * a root near phi=pi: the reciprocal chart treats it near zero instead.
 */
std::vector<Real> Roots(Polynomial p)
{
    while(p.size()>1 && p.back()==0)p.pop_back();
    Real scale=0;for(Real c:p)scale=std::max(scale,std::abs(c));
    if(scale==0)return {-1,1};
    for(Real &c:p)c/=scale;
    if(p.size()==1)return {};
    if(p.size()==2)
    {
        const Real root=-p[0]/p[1];
        return std::isfinite(root)&&root>=-1&&root<=1?std::vector<Real>{root}:std::vector<Real>{};
    }
    Polynomial derivative(p.size()-1);
    for(std::size_t i=1;i<p.size();++i)derivative[i-1]=i*p[i];
    auto critical=Roots(derivative);
    std::vector<Real> boundaries{-1};
    for(Real root:critical)if(root>-1 && root<1)boundaries.push_back(root);
    boundaries.push_back(1);std::sort(boundaries.begin(),boundaries.end());
    boundaries.erase(std::unique(boundaries.begin(),boundaries.end()),boundaries.end());
    Real norm=0;for(Real c:p)norm+=std::abs(c);
    const Real tolerance=256*std::numeric_limits<Real>::epsilon()*norm;
    std::vector<Real> result;
    for(Real x:boundaries)if(std::abs(Evaluate(p,x))<=tolerance)result.push_back(x);
    for(std::size_t i=1;i<boundaries.size();++i)
    {
        Real left=boundaries[i-1],right=boundaries[i];
        Real fl=Evaluate(p,left),fr=Evaluate(p,right);
        if(fl==0 || fr==0 || std::signbit(fl)==std::signbit(fr))continue;
        while(true)
        {
            const Real middle=left+(right-left)/2;
            if(middle==left || middle==right ||
               right-left<=32*std::numeric_limits<Real>::epsilon()*(1+std::abs(middle)))break;
            const Real fm=Evaluate(p,middle);
            if(fm==0){left=right=middle;break;}
            if(std::signbit(fm)==std::signbit(fl)){left=middle;fl=fm;}
            else{right=middle;fr=fm;}
        }
        result.push_back(left+(right-left)/2);
    }
    return result;
}
Real Contract(const Tensor &a,const Tensor &b)
{
    Real value=0;for(int i=0;i<3;++i)for(int j=0;j<3;++j)value+=Real(a[i][j])*b[i][j];
    return value;
}
}

PlanarAdmissibility CheckPlanarState(double k,double epsilon,const Tensor &gradient,
    double nu,double auxiliary,double r,double H,bool principal,const WallStress*currentValue,
    ConstitutiveBackend backend,WallEvaluationStats*stats)
{
    Require(std::isfinite(k)&&k>0&&std::isfinite(epsilon)&&epsilon>0&&std::isfinite(nu)&&nu>0&&std::isfinite(auxiliary)&&auxiliary>=0&&std::isfinite(r)&&r>=0&&std::isfinite(H)&&H>0,"invalid planar current state inputs");
    for(const auto&row:gradient)for(double value:row)Require(std::isfinite(value),"nonfinite planar gradient");
    for(int i=0;i<3;++i)
        Require(gradient[1][i]==0 && gradient[i][1]==0,"planar admissibility requires inactive spanwise gradient");
    WallStress evaluated;if(!currentValue)evaluated=AdaptAmended(k,epsilon,gradient,{0,0,1},nu,auxiliary,r,H,backend,stats);
    const auto&wall=currentValue?*currentValue:evaluated;
    Require(std::isfinite(k)&&k>0&&std::isfinite(wall.mu)&&wall.mu>0,"invalid supplied current wall state");
    Eigen::Matrix3d covariance;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {
        Require(std::isfinite(wall.R[i][j]),"nonfinite current covariance");
        Require(std::abs(wall.R[i][j]-wall.R[j][i])<=128*std::numeric_limits<double>::epsilon()*
            std::max({k,std::abs(wall.R[i][j]),std::abs(wall.R[j][i])}),"nonsymmetric current covariance");
        covariance(i,j)=.5*(wall.R[i][j]+wall.R[j][i])/k;
    }
    PlanarAdmissibility result;
    if(covariance(0,1)==0&&covariance(1,2)==0){
        // This is the full 3D covariance spectrum of a planar state: its
        // spanwise eigenvalue plus both eigenvalues of the x-z block.
        result.covarianceOverK=std::min(covariance(1,1),.5*(covariance(0,0)+covariance(2,2))-
            std::hypot(.5*(covariance(0,0)-covariance(2,2)),covariance(0,2)));
    }else{
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen(covariance,Eigen::EigenvaluesOnly);
        Require(eigen.info()==Eigen::Success && eigen.eigenvalues().allFinite(),"covariance eigensolve failed");result.covarianceOverK=eigen.eigenvalues()[0];
    }
    if(!principal)return result;
    result.principalEvaluated=true;
    const Real pi=std::acos(Real(-1));
    if(wall.chi==0)
    {
        // H=l tensor n: H:nu(H+H^T)=nu and the auxiliary contribution
        // is aux*cos(2theta)^2. theta=pi/4 reaches the exact lower bound.
        result.principalOverMu=result.principalRawOverMu=nu/wall.mu;
        result.principalAngle=double(pi/4);result.stationaryCandidates=2;
        return result;
    }
    std::array<Tensor,3> basis{};
    basis[0][0][2]=-.5;basis[0][2][0]=.5;
    basis[1][0][2]=basis[1][2][0]=.5;
    basis[2][0][0]=-.5;basis[2][2][2]=.5;
    Real q[3][3]{},qnorm=0;
    for(int j=0;j<3;++j)
    {
        const auto action=ApplyAmendedDirection(k,epsilon,gradient,{0,0,1},nu,auxiliary,r,H,0,0,basis[j],0,0,backend,stats).V;
        for(int i=0;i<3;++i)
        {
            q[i][j]=Contract(basis[i],action)/wall.mu;
            Require(std::isfinite(q[i][j]),"nonfinite principal tangent action");qnorm+=std::abs(q[i][j]);
        }
    }
    const Real a=q[0][0]+(q[1][1]+q[2][2])/2,b=q[0][1]+q[1][0],c=q[0][2]+q[2][0],
        d=(q[1][1]-q[2][2])/2,e=(q[1][2]+q[2][1])/2;
    const Polynomial polynomial{c+2*e,-2*b-8*d,-12*e,-2*b+8*d,-c+2*e};
    std::vector<Real> angles{0,pi};
    for(Real t:Roots(polynomial))angles.push_back(2*std::atan(t));
    Polynomial reciprocal(polynomial.rbegin(),polynomial.rend());
    for(Real u:Roots(reciprocal))angles.push_back(2*std::atan2(Real(1),u));
    Real minimum=std::numeric_limits<Real>::infinity(),phi=0;
    for(Real angle:angles)
    {
        const Real value=a+b*std::cos(angle)+c*std::sin(angle)+d*std::cos(2*angle)+e*std::sin(2*angle);
        if(value<minimum){minimum=value;phi=angle;}
    }
    const Real residual=std::abs(-b*std::sin(phi)+c*std::cos(phi)-2*d*std::sin(2*phi)+2*e*std::cos(2*phi));
    const Real norm=1+std::abs(a)+std::abs(b)+std::abs(c)+std::abs(d)+std::abs(e)+qnorm;
    Require(std::isfinite(minimum) && residual<=1e-11L*norm,"stationary principal root resolution failed");
    // Covers floating tangent/harmonic arithmetic and angular root rounding.
    // This is a numerical reserve, not an interval-arithmetic theorem about
    // arbitrary three-dimensional modes or the full coupled PDE.
    const Real reserve=2048*std::numeric_limits<double>::epsilon()*norm;
    Real theta=phi/2;if(theta<0)theta+=pi;if(theta>=pi)theta-=pi;
    result.principalRawOverMu=double(minimum);
    result.principalReserveOverMu=double(reserve);
    result.principalOverMu=double(minimum-reserve);
    result.principalAngle=double(theta);result.stationaryResidual=double(residual);
    result.stationaryCandidates=int(angles.size());return result;
}
}
