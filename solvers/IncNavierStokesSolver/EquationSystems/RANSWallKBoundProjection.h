// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_K_BOUND_PROJECTION_H
#define NEKTAR_RANS_WALL_K_BOUND_PROJECTION_H
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace Nektar::RANSWallKBound
{
using Vector=std::vector<double>;
using Matrix=std::vector<Vector>;
inline double Dot(const Vector &a,const Vector &b)
{long double s=0;for(std::size_t i=0;i<a.size();++i)s+=(long double)a[i]*b[i];return double(s);}
inline double Norm(const Vector &a){return std::sqrt(std::max(0.,Dot(a,a)));}
inline std::string Fingerprint(const Vector &values)
{
    std::uint64_t h=1469598103934665603ULL;
    for(double x:values){std::uint64_t bits;std::memcpy(&bits,&x,8);for(int j=0;j<8;++j){h^=(bits>>(8*j))&255;h*=1099511628211ULL;}}
    std::ostringstream out;out<<std::hex<<std::setw(16)<<std::setfill('0')<<h;return out.str();
}
struct Controls
{
    int maximumPerAdd=2000,maximumTotal=2000;
    double maximumSeconds=900.;
};
struct Event
{
    std::string kind,activeFingerprint;
    int add=0,localWork=0,totalWork=0,entering=-1,leaving=-1,activeSize=0;
    double step=0,minimumQR=0,maximumQR=0,objective=0,dualViolation=0,seconds=0;
};

// Small serial consistent-mass projection. The original whitening, normalized
// rows, two-pass MGS, tolerances and working-set updates are retained.
// Optimize solves its WORKING SET. The physical separation wrapper must rescan
// all inactive cuts, including rows that previously left. A zero working-set
// primal-dual gap is not a certificate for an infeasible accumulated pool.
class Projection
{
public:
    using Clock=std::chrono::steady_clock;
    struct Cut {Vector row;double rhs;int witness;std::string rowIdentity,constraintIdentity;};
    explicit Projection(const Matrix &mass):Projection(mass,Controls{}){}
    Projection(const Matrix &mass,const Controls &limits,std::function<void(const Event&)> sink={}):
        lower(mass),n(mass.size()),y(n),delta(n),controls(limits),trace(std::move(sink)),started(Clock::now())
    {
        if(n==0 || std::any_of(mass.begin(),mass.end(),[&](const auto &r){return r.size()!=std::size_t(n);}))
            throw std::logic_error("invalid k-bound mass shape");
        if(controls.maximumPerAdd<=0 || controls.maximumTotal<=0 || !(controls.maximumSeconds>0) || !std::isfinite(controls.maximumSeconds))
            throw std::logic_error("invalid bounded projection work controls");
        for(int i=0;i<n;++i)for(int j=0;j<=i;++j)
        {
            long double s=mass[i][j];for(int k=0;k<j;++k)s-=(long double)lower[i][k]*lower[j][k];
            if(i==j){if(!(s>0))throw std::runtime_error("k-bound mass is not positive definite");lower[i][j]=std::sqrt(s);}
            else lower[i][j]=s/lower[j][j];
        }
        factorizationSeconds=Seconds();
    }
    double Seconds() const {return std::chrono::duration<double>(Clock::now()-started).count();}
    std::string ActiveFingerprint() const
    {
        std::vector<std::string> keys;for(int id:active)keys.push_back(cuts[id].constraintIdentity);
        std::sort(keys.begin(),keys.end());Vector bytes;for(const auto &s:keys)for(unsigned char c:s)bytes.push_back(c);
        return Fingerprint(bytes);
    }
    void Emit(const std::string &kind,int leaving=-1,double step=0,double qrMin=0,double qrMax=0)
    {
        if(!trace)return;
        double violation=0;for(double x:lambda)violation=std::max(violation,-x);
        trace(Event{kind,ActiveFingerprint(),adds,localIterations,iterations,entering,leaving,int(active.size()),step,qrMin,qrMax,
            .5*Dot(y,y),violation,Seconds()});
    }
    bool Add(const Vector &row,double rhs,int witness)
    {
        if(row.size()!=std::size_t(n) || !std::isfinite(rhs))throw std::logic_error("invalid projection row");
        ++adds;localIterations=0;workComplete=false;failure.clear();
        const auto rowIdentity=Fingerprint(row);auto constraint=row;constraint.push_back(rhs);const auto identity=Fingerprint(constraint);
        if(++rowVisits[rowIdentity]>1)++repeatedRows;
        if(++constraintVisits[identity]>1)
        {
            ++repeatedConstraints;
            if(std::none_of(active.begin(),active.end(),[&](int id){return cuts[id].constraintIdentity==identity;}))++rediscoveredConstraints;
        }
        Vector a=row;for(int i=0;i<n;++i){for(int j=0;j<i;++j)a[i]-=lower[i][j]*a[j];a[i]/=lower[i][i];}
        ++massFactorApplications;
        const double scale=Norm(a);
        if(!(scale>0)){failure="incompatible_prescribed_k_floor";Emit("invalid_row");return false;}
        for(auto &x:a)x/=scale;
        cuts.push_back({a,rhs/scale,witness,rowIdentity,identity});entering=cuts.size()-1;
        active.push_back(entering);lambda.push_back(0.);Emit("add");
        const auto before=Clock::now();const bool ok=Optimize();
        algebraSeconds+=std::chrono::duration<double>(Clock::now()-before).count();
        maximumAddIterations=std::max(maximumAddIterations,localIterations);return ok;
    }
    bool Optimize()
    {
        for(;;)
        {
            // Count completed loop bodies separately from rejected guards.
            ++workChecks;
            if(iterations>=controls.maximumTotal){++totalLimitChecks;failure="k_bound_total_active_set_limit";Emit("total_limit");return false;}
            if(localIterations>=controls.maximumPerAdd){++addLimitChecks;failure="k_bound_per_add_active_set_limit";Emit("add_limit");return false;}
            if(Seconds()>=controls.maximumSeconds){++timeLimitChecks;failure="k_bound_algebra_wall_time_limit";Emit("time_limit");return false;}
            ++iterations;++localIterations;
            const auto fingerprint=ActiveFingerprint();if(++workingSetVisits[fingerprint]>1)++repeatedWorkingSets;
            const int m=active.size();maximumActive=std::max(maximumActive,m);Matrix q,r(m,Vector(m));bool dependent=false;
            double qrMin=1.,qrMax=0.;
            for(int i=0;i<m;++i)
            {
                auto a=cuts[active[i]].row;
                for(int pass=0;pass<2;++pass)for(int j=0;j<q.size();++j)
                {const double c=Dot(a,q[j]);r[j][i]+=c;for(int k=0;k<n;++k)a[k]-=c*q[j][k];}
                const double length=Norm(a);qrMin=std::min(qrMin,length);qrMax=std::max(qrMax,length);
                minimumObservedQR=std::min(minimumObservedQR,length);
                if(length<2e-11)
                {
                    ++dependentEntries;
                    if(i!=m-1){failure="k_bound_old_active_rank_loss";Emit("old_rank_loss",-1,0,qrMin,qrMax);return false;}
                    Vector c(i);for(int j=0;j<i;++j)c[j]=r[j][i];
                    for(int j=i-1;j>=0;--j){for(int k=j+1;k<i;++k)c[j]-=r[j][k]*c[k];c[j]/=r[j][j];}
                    double t=std::numeric_limits<double>::infinity();int leaving=-1;
                    for(int j=0;j<i;++j)if(c[j]>1e-13 && lambda[j]/c[j]<t){t=lambda[j]/c[j];leaving=j;}
                    if(leaving<0){failure="k_bound_dependent_exchange_unresolved";Emit("dependent_no_exchange",-1,0,qrMin,qrMax);return false;}
                    for(int j=0;j<i;++j){lambda[j]-=t*c[j];}
                    lambda.back()+=t;Emit("dependent_leave",active[leaving],t,qrMin,qrMax);
                    active.erase(active.begin()+leaving);lambda.erase(lambda.begin()+leaving);
                    dependent=true;break;
                }
                r[i][i]=length;for(auto &x:a)x/=length;q.push_back(std::move(a));
            }
            if(dependent)continue;
            Vector t(m),candidate(m);
            for(int i=0;i<m;++i){t[i]=cuts[active[i]].rhs;for(int j=0;j<i;++j)t[i]-=r[j][i]*t[j];t[i]/=r[i][i];}
            candidate=t;
            for(int i=m-1;i>=0;--i){for(int j=i+1;j<m;++j)candidate[i]-=r[i][j]*candidate[j];candidate[i]/=r[i][i];}
            ++kktSolves;
            double fraction=1.;int leaving=-1;
            for(int i=0;i<m;++i)if(candidate[i]<0 && candidate[i]<lambda[i])
            {const double f=lambda[i]/(lambda[i]-candidate[i]);if(f<fraction){fraction=f;leaving=i;}}
            if(leaving>=0)
            {
                for(int i=0;i<m;++i)lambda[i]+=fraction*(candidate[i]-lambda[i]);
                Emit("dual_leave",active[leaving],fraction,qrMin,qrMax);
                active.erase(active.begin()+leaving);lambda.erase(lambda.begin()+leaving);continue;
            }
            lambda=std::move(candidate);std::fill(y.begin(),y.end(),0.);
            for(int i=0;i<m;++i)for(int j=0;j<n;++j)y[j]+=q[i][j]*t[i];
            delta=y;
            for(int i=n-1;i>=0;--i){for(int j=i+1;j<n;++j)delta[i]-=lower[j][i]*delta[j];delta[i]/=lower[i][i];}
            ++massFactorApplications;
            solvedActive=active;solvedLambda=lambda;solvedAdd=adds;workComplete=true;
            lastQRMin=qrMin;lastQRMax=qrMax;Emit("solved",-1,1.,qrMin,qrMax);return true;
        }
    }
    Vector DualAction() const
    {Vector result(n);for(int i=0;i<solvedActive.size();++i)for(int j=0;j<n;++j)result[j]+=solvedLambda[i]*cuts[solvedActive[i]].row[j];return result;}
    double Stationarity() const
    {auto r=y;for(int i=0;i<solvedActive.size();++i)for(int j=0;j<n;++j)r[j]-=solvedLambda[i]*cuts[solvedActive[i]].row[j];return Norm(r);}
    double Complementarity() const
    {double r=0;for(int i=0;i<solvedActive.size();++i)r=std::max(r,std::abs(solvedLambda[i]*(Dot(cuts[solvedActive[i]].row,y)-cuts[solvedActive[i]].rhs)));return r;}
    double DualViolation() const {double r=0;for(double x:solvedLambda)r=std::max(r,-x);return r;}
    double PoolViolation() const {double r=0;for(const auto &c:cuts)r=std::max(r,c.rhs-Dot(c.row,y));return r;}
    double Objective() const {return .5*Dot(y,y);}
    double DualObjective() const
    {double b=0;for(int i=0;i<solvedActive.size();++i)b+=solvedLambda[i]*cuts[solvedActive[i]].rhs;const auto action=DualAction();return b-.5*Dot(action,action);}
    Matrix lower;int n;
    Vector y,delta,lambda,solvedLambda;
    std::vector<Cut> cuts;
    std::vector<int> active,solvedActive;
    int iterations=0,localIterations=0,maximumAddIterations=0,workChecks=0,totalLimitChecks=0,addLimitChecks=0,timeLimitChecks=0;
    int adds=0,solvedAdd=0,entering=-1,massFactorApplications=0,kktSolves=0,dependentEntries=0,maximumActive=0;
    int repeatedRows=0,repeatedConstraints=0,rediscoveredConstraints=0,repeatedWorkingSets=0;
    bool workComplete=true;
    double minimumObservedQR=1.,lastQRMin=0,lastQRMax=0,algebraSeconds=0,factorizationSeconds=0;
    Controls controls;std::string failure;
    std::map<std::string,int> rowVisits,constraintVisits,workingSetVisits;
private:
    std::function<void(const Event&)> trace;Clock::time_point started;
};
}
#endif
