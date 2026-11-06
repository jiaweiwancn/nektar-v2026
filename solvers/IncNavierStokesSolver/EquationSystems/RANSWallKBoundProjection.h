// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_K_BOUND_PROJECTION_H
#define NEKTAR_RANS_WALL_K_BOUND_PROJECTION_H
#include <algorithm>
#include <cmath>
#include <limits>
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

// Small serial consistent-mass projection. Whitening, QR and dual active-set
// work are algebraic operations; this class has no transport inverse or fields.
// Add/Optimize maintain the working active set. Inactive physical constraints
// can become violated after another cut is added: the separation wrapper must
// rescan its complete physical inventory, including previously dropped cuts.
class Projection
{
public:
    struct Cut {Vector row;double rhs;int witness;};
    explicit Projection(const Matrix &mass):lower(mass),n(mass.size()),y(n),delta(n)
    {
        if(n==0 || std::any_of(mass.begin(),mass.end(),[&](const auto &r){return r.size()!=std::size_t(n);}))
            throw std::logic_error("invalid k-bound mass shape");
        for(int i=0;i<n;++i)for(int j=0;j<=i;++j)
        {
            long double s=mass[i][j];for(int k=0;k<j;++k)s-=(long double)lower[i][k]*lower[j][k];
            if(i==j){if(!(s>0))throw std::runtime_error("k-bound mass is not positive definite");lower[i][j]=std::sqrt(s);}
            else lower[i][j]=s/lower[j][j];
        }
    }
    bool Add(const Vector &row,double rhs,int witness)
    {
        Vector a=row;for(int i=0;i<n;++i){for(int j=0;j<i;++j)a[i]-=lower[i][j]*a[j];a[i]/=lower[i][i];}
        ++massFactorApplications;
        const double scale=Norm(a);
        if(!(scale>0)){failure="incompatible_prescribed_k_floor";return false;}
        for(auto &x:a)x/=scale;
        cuts.push_back({a,rhs/scale,witness});active.push_back(cuts.size()-1);lambda.push_back(0.);
        return Optimize();
    }
    bool Optimize()
    {
        while(++iterations<=2000)
        {
            const int m=active.size();Matrix q,r(m,Vector(m));bool dependent=false;
            for(int i=0;i<m;++i)
            {
                auto a=cuts[active[i]].row;
                for(int pass=0;pass<2;++pass)for(int j=0;j<q.size();++j)
                {const double c=Dot(a,q[j]);r[j][i]+=c;for(int k=0;k<n;++k)a[k]-=c*q[j][k];}
                const double length=Norm(a);
                if(length<2e-11)
                {
                    // A dependent entering row can replace a weaker face.
                    // Increase its dual weight while preserving y, until an
                    // old positive dual coefficient reaches zero.
                    if(i!=m-1){failure="k_bound_old_active_rank_loss";return false;}
                    Vector c(i);for(int j=0;j<i;++j)c[j]=r[j][i];
                    for(int j=i-1;j>=0;--j){for(int k=j+1;k<i;++k)c[j]-=r[j][k]*c[k];c[j]/=r[j][j];}
                    double t=std::numeric_limits<double>::infinity();int leaving=-1;
                    for(int j=0;j<i;++j)if(c[j]>1e-13 && lambda[j]/c[j]<t){t=lambda[j]/c[j];leaving=j;}
                    if(leaving<0){failure="k_bound_incompatible_dependent_constraints";return false;}
                    for(int j=0;j<i;++j){lambda[j]-=t*c[j];}
                    lambda.back()+=t;
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
                active.erase(active.begin()+leaving);lambda.erase(lambda.begin()+leaving);continue;
            }
            lambda=std::move(candidate);std::fill(y.begin(),y.end(),0.);
            for(int i=0;i<m;++i)for(int j=0;j<n;++j)y[j]+=q[i][j]*t[i];
            delta=y;
            for(int i=n-1;i>=0;--i){for(int j=i+1;j<n;++j)delta[i]-=lower[j][i]*delta[j];delta[i]/=lower[i][i];}
            ++massFactorApplications;
            return true;
        }
        failure="k_bound_active_set_iteration_limit";return false;
    }
    double Stationarity() const
    {auto r=y;for(int i=0;i<active.size();++i)for(int j=0;j<n;++j)r[j]-=lambda[i]*cuts[active[i]].row[j];return Norm(r);}
    double Complementarity() const
    {double r=0;for(int i=0;i<active.size();++i)r=std::max(r,std::abs(lambda[i]*(Dot(cuts[active[i]].row,y)-cuts[active[i]].rhs)));return r;}
    Matrix lower;int n;
    Vector y,delta,lambda;
    std::vector<Cut> cuts;
    std::vector<int> active;
    int iterations=0,massFactorApplications=0,kktSolves=0;
    std::string failure;
};
}
#endif
