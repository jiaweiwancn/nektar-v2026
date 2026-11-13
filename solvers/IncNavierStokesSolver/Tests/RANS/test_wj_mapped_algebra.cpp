// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallMappedAlgebra.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallKBoundProjection.h>
#include <fstream>
#include <iostream>
#include <cmath>

using namespace Nektar;
using namespace Nektar::RANSWallAlgebra;
namespace
{
void close(const Vector &a,const Vector &b,double tolerance=1e-11)
{Require(a.size()==b.size(),"fixture shape differs");Vector d=a;for(int i=0;i<d.size();++i)d[i]-=b[i];Require(Maximum(d)<=tolerance*std::max(1.,Maximum(b)),"mapped/dense fixture parity failed");}
template<class F> void rejects(F f)
{bool caught=false;try{f();}catch(const std::exception&){caught=true;}Require(caught,"unsupported mapped operation did not reject");}
double dot(const Vector &a,const Vector &b)
{long double s=0;for(int i=0;i<a.size();++i)s+=(long double)a[i]*b[i];return s;}
MappedMatrix sparse(const DenseMatrix &a,bool mass=false)
{MappedMatrix s(a.size(),a[0].size(),true);for(int i=0;i<a.size();++i)for(int j=0;j<a[i].size();++j)s.Set(i,j,a[i][j]);if(mass)s.ConsistentMass();return s;}
}

int main(int argc,char **argv)
{
    try
    {
        if(argc==2){std::ifstream maps("/proc/self/maps");std::ofstream out(argv[1]);out<<maps.rdbuf();Require(out.good(),"fixture native-exec maps write failed");}
        const int n=35;DenseMatrix mass(n,Vector(n));
        // Scrambled native IDs; a path graph with nontrivial signs/couplings.
        for(int i=0;i<n;++i)mass[i][i]=2.+.03*i;
        for(int i=0;i<n-1;++i){int a=(i*12)%n,b=((i+1)*12)%n;mass[a][b]=mass[b][a]=-.2;}
        auto sm=sparse(mass,true);MappedBandFactor mf(sm,true);Vector x(n),r(n),b(n);
        for(int i=0;i<n;++i){x[i]=std::sin(.7*i)+.1;r[i]=std::cos(.2*i);b[i]=.3-.01*i;}
        close(Multiply(sm,x),Multiply(mass,x));close(mf.Apply(b),Solve(mass,b));
        close(mf.Recover(mf.Whiten(r)),Solve(mass,r));
        Require(std::abs(dot(mf.Whiten(r),mf.Whiten(b))-dot(r,Solve(mass,b)))<2e-12,"band whitening is not the consistent mass metric");
        Require(mf.LowerBandwidth()<=2 && mf.StorageBytes()<n*n*8,"RCM did not expose narrow actual graph");
        rejects([&]{MappedBandFactor bad(sm,true,1);});
        rejects([&]{(void)static_cast<const DenseMatrix&>(sm);});
        rejects([&]{sm.PressureDense(4);});
        // Nonsymmetric, pivot-requiring LU and transpose recovery.
        DenseMatrix a=mass;for(int i=0;i<n-1;++i)a[i][i+1]+=.13;
        a[0][0]=0.;a[0][1]=2.;a[1][0]=5.;auto sa=sparse(a);MappedBandFactor af(sa);
        close(Multiply(a,af.Apply(b)),b);DenseMatrix at(n,Vector(n));
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)at[i][j]=a[j][i];
        close(Multiply(at,af.Apply(b,true)),b);
        // Restriction preserves explicit native order and dependent columns.
        std::vector<int> free{8,2,5,1};auto sub=sm.Restrict(free,free);Vector v{.2,-.3,.7,.8};
        Vector full(n);for(int i=0;i<free.size();++i)full[free[i]]=v[i];auto action=Multiply(mass,full);Vector expected;
        for(int row:free)expected.push_back(action[row]);close(sub.Apply(v),expected);
        auto epsilon=sparse(DenseMatrix{{3.,-.7,.2},{-.4,4.,-.3},{.2,-.6,2.}});
        epsilon.ClearRow(0);epsilon.Set(0,0,1.);
        Require(epsilon.Get(1,0)==-.4 && epsilon.Get(2,0)==.2,"dependent trace transport columns were erased");
        close(epsilon.Apply(MappedBandFactor(epsilon).Apply(Vector{.1,.2,.3})),Vector{.1,.2,.3});
        // Activated retained-pool and near-dependent cut: same finite QP,
        // different exact Cholesky coordinates. Compare physical/KKT results.
        DenseMatrix qm{{2.,.2,0.,0.},{.2,3.,.1,0.},{0.,.1,4.,.2},{0.,0.,.2,2.}};
        auto factor=std::make_shared<MappedBandFactor>(sparse(qm,true),true);
        RANSWallKBound::Controls controls;controls.retainedPool=true;
        RANSWallKBound::Projection dense(qm,controls);
        RANSWallKBound::Projection::Metric metric{4,[factor](const Vector &v){return factor->Whiten(v);},
            [factor](const Vector &v){return factor->Recover(v);},RANSWallKBound::Projection::Clock::now()};
        RANSWallKBound::Projection band(metric,controls);
        std::vector<Vector> cuts{{1,0,0,0},{-1,1,0,0},{1,0,1e-13,0},{0,0,0,1}};
        Vector rhs{.4,.1,.45,.2};
        for(int i=0;i<cuts.size();++i)
        {Require(dense.Add(cuts[i],rhs[i],i) && band.Add(cuts[i],rhs[i],i),"activated bound fixture failed");close(band.delta,dense.delta,1e-9);}
        Require(band.dependentEntries>0 && dense.dependentEntries>0,"near-dependent metric fixture not exercised");
        Require(band.poolComplete && band.PoolViolation()<1e-12 && band.Stationarity()<1e-12 && band.Complementarity()<1e-12 && band.DualViolation()<1e-12,"mapped bound KKT/pool certificate failed");
        Require(std::abs(.5*dot(band.delta,Multiply(qm,band.delta))-band.Objective())<2e-12,"bound objective uses wrong metric/permutation");
        // Streamed nonsymmetric pressure Schur with actual MASS inverse.
        const DenseMatrix M{{2.,.2,0.},{.2,3.,.1},{0.,.1,2.}},A{{5.,.7,0.},{-.2,4.,.2},{.1,0.,3.}},
            G{{1.,.2},{-.3,1.},{.4,.5}},D{{.2,1.,.3},{1.,-.4,.1}};
        auto bm=sparse(M,true),ba=sparse(A);MappedBandFactor fm(bm,true),fa(ba);
        Vector ru{.1,-.2,.3},rp{.002,-.001};auto predictor=fa.Apply(ru);
        DenseMatrix schur(2,Vector(2)),reference(2,Vector(2));double dt=.01;
        for(int j=0;j<2;++j)
        {
            Vector col;for(int i=0;i<3;++i)col.push_back(G[i][j]);auto response=fm.Apply(col),ref=Solve(M,col);
            for(auto &v:response)v*=dt;for(auto &v:ref)v*=dt;
            auto s=Multiply(D,response),t=Multiply(D,ref);for(int i=0;i<2;++i){schur[i][j]=s[i];reference[i][j]=t[i];}
        }
        Require(std::abs(schur[0][1]-schur[1][0])>1e-5,"Schur fixture accidentally symmetric");
        auto prhs=Multiply(D,predictor);for(int i=0;i<2;++i)prhs[i]+=rp[i];auto dp=Solve(schur,prhs);
        auto gdp=Multiply(G,dp),correction=fm.Apply(gdp);for(auto &v:correction)v*=dt;
        close(Multiply(D,correction),prhs);Vector du=predictor;for(int i=0;i<3;++i)du[i]-=correction[i];
        auto continuity=Multiply(D,du);for(int i=0;i<2;++i)continuity[i]+=rp[i];close(continuity,Vector(2));
        auto wrong=fa.Apply(gdp);for(auto &v:wrong)v*=dt;
        Vector different=wrong;for(int i=0;i<3;++i)different[i]-=correction[i];Require(Maximum(different)>1e-5,"fixture cannot detect forbidden full-momentum inverse");
        std::cout<<"PASS mapped mass/action/restriction, pivoted band LU, exact whitening/recovery, activated/near-dependent retained bounds, streamed nonsymmetric consistent-mass Schur; no governing step\n";
        return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
