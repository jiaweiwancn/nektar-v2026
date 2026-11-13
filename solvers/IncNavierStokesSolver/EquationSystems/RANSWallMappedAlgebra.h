// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_MAPPED_ALGEBRA_H
#define NEKTAR_RANS_WALL_MAPPED_ALGEBRA_H

#include <IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.h>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <deque>
#include <limits>
#include <memory>
#include <numeric>
#include <set>

namespace Nektar::RANSWallAlgebra
{
// A narrow assembly sink for the custom wall path. Legacy callers retain their
// dense rows; mapped callers explicitly select sparse rows before allocation.
// A sparse matrix is never implicitly expanded into a global dense square.
class MappedMatrix
{
public:
    MappedMatrix() = default;
    MappedMatrix(const DenseMatrix &a): dense(a), n(a.size()), m(a.empty()?0:a[0].size()) {}
    MappedMatrix(int rows,const Vector &row):dense(rows,row),n(rows),m(row.size()) {}
    MappedMatrix(int rows,int cols,bool sparse):n(rows),m(cols),compressed(sparse)
    {
        Require(rows>=0 && cols>=0,"negative wall matrix shape");
        if(sparse)entries.resize(n);else dense.assign(n,Vector(m));
    }
    void assign(int rows,const Vector &row) {*this=MappedMatrix(rows,row);}
    std::size_t size() const {return n;}
    int Columns() const {return m;}
    bool empty() const {return n==0;}
    bool Sparse() const {return compressed;}
    bool MassMetric() const {return massMetric;}
    MappedMatrix &ConsistentMass() {massMetric=true;return *this;}
    Vector &operator[](int i) {Require(!compressed,"dense row requested from mapped matrix");return dense.at(i);}
    const Vector &operator[](int i) const {Require(!compressed,"dense row requested from mapped matrix");return dense.at(i);}
    auto begin() const {Require(!compressed,"dense iteration requested from mapped matrix");return dense.begin();}
    auto end() const {Require(!compressed,"dense iteration requested from mapped matrix");return dense.end();}
    operator const DenseMatrix &() const {Require(!compressed,"implicit dense mapped matrix conversion is forbidden");return dense;}
    double Get(int i,int j) const
    {
        Require(i>=0 && i<n && j>=0 && j<m,"wall matrix index outside shape");
        if(!compressed)return dense[i][j];
        const auto it=entries[i].find(j);return it==entries[i].end()?0.:it->second;
    }
    void Add(int i,int j,double x)
    {
        Require(i>=0 && i<n && j>=0 && j<m,"wall matrix insertion outside shape");
        if(compressed){if(x!=0.)entries[i][j]+=x;}else dense[i][j]+=x;
    }
    void Set(int i,int j,double x)
    {
        Require(i>=0 && i<n && j>=0 && j<m,"wall matrix insertion outside shape");
        if(compressed){if(x==0.)entries[i].erase(j);else entries[i][j]=x;}else dense[i][j]=x;
    }
    void ClearRow(int i)
    {if(compressed)entries.at(i).clear();else std::fill(dense.at(i).begin(),dense.at(i).end(),0.);}
    template<class F> void Each(int i,F f) const
    {
        if(compressed){for(auto [j,v]:entries.at(i))if(v!=0.)f(j,v);}
        else for(int j=0;j<m;++j)if(dense[i][j]!=0.)f(j,dense[i][j]);
    }
    Vector Apply(const Vector &x) const
    {
        Require(x.size()==std::size_t(m),"mapped action vector shape mismatch");
        Vector y(n);for(int i=0;i<n;++i){long double sum=0.;Each(i,[&](int j,double v){sum+=(long double)v*x[j];});y[i]=sum;}
        Maximum(y);return y;
    }
    MappedMatrix Restrict(const std::vector<int> &rows,const std::vector<int> &cols) const
    {
        MappedMatrix a(rows.size(),cols.size(),compressed);std::vector<int> index(m,-1);
        for(int j=0;j<cols.size();++j){Require(cols[j]>=0 && cols[j]<m && index[cols[j]]<0,"invalid mapped restriction");index[cols[j]]=j;}
        for(int i=0;i<rows.size();++i){Require(rows[i]>=0 && rows[i]<n,"invalid mapped restriction row");Each(rows[i],[&](int j,double v){if(index[j]>=0)a.Add(i,index[j],v);});}
        if(massMetric && rows==cols)a.ConsistentMass();return a;
    }
    std::uint64_t Nonzeros() const
    {std::uint64_t result=0;for(int i=0;i<n;++i)Each(i,[&](int,double){++result;});return result;}
    DenseMatrix PressureDense(int maximum=1024) const
    {
        Require(n<=maximum && m<=maximum,"only bounded pressure matrices may be expanded");
        DenseMatrix result(n,Vector(m));for(int i=0;i<n;++i)Each(i,[&](int j,double v){result[i][j]=v;});return result;
    }
private:
    DenseMatrix dense;
    std::vector<std::map<int,double>> entries;
    int n=0,m=0;bool compressed=false,massMetric=false;
};

// Explicit native->band permutation: perm[band_index]=native_index.
// Reverse Cuthill-McKee uses the union graph, even for nonsymmetric operators.
inline std::vector<int> BandPermutation(const MappedMatrix &a)
{
    const int n=a.size();Require(a.Columns()==n && n>0,"band factor needs a square matrix");
    std::vector<std::vector<int>> graph(n);
    for(int i=0;i<n;++i)a.Each(i,[&](int j,double){if(i!=j){graph[i].push_back(j);graph[j].push_back(i);}});
    for(auto &row:graph){std::sort(row.begin(),row.end());row.erase(std::unique(row.begin(),row.end()),row.end());}
    auto less=[&](int i,int j){return graph[i].size()!=graph[j].size()?graph[i].size()<graph[j].size():i<j;};
    std::vector<int> seeds(n),order;std::iota(seeds.begin(),seeds.end(),0);std::stable_sort(seeds.begin(),seeds.end(),less);
    std::vector<bool> seen(n,false);
    for(int seed:seeds)if(!seen[seed])
    {
        std::deque<int> queue{seed};seen[seed]=true;
        while(!queue.empty())
        {
            int i=queue.front();queue.pop_front();order.push_back(i);auto neighbors=graph[i];std::stable_sort(neighbors.begin(),neighbors.end(),less);
            for(int j:neighbors)if(!seen[j]){seen[j]=true;queue.push_back(j);}
        }
    }
    std::reverse(order.begin(),order.end());return order;
}

class MappedBandFactor
{
public:
    // Factor-storage ceiling is checked with 64-bit arithmetic BEFORE allocation.
    explicit MappedBandFactor(const MappedMatrix &a,bool positiveMass=false,
        std::uint64_t maximumBytes=512ULL*1024*1024):n(a.size()),spd(positiveMass),perm(BandPermutation(a)),inverse(n)
    {
        for(int i=0;i<n;++i)inverse[perm[i]]=i;
        for(int i=0;i<n;++i)a.Each(i,[&](int j,double v){
            kl=std::max(kl,inverse[i]-inverse[j]);ku=std::max(ku,inverse[j]-inverse[i]);
            if(spd)Require(std::abs(v-a.Get(j,i))<=1e-12*std::max(1.,std::abs(v)),"consistent mass is not symmetric");});
        kd=std::max(kl,ku);ld=spd?kd+1:2*kl+ku+1;
        const std::uint64_t bytes=std::uint64_t(n)*ld*sizeof(double);
        Require(bytes<=maximumBytes,"mapped band fill exceeds factor storage reservation");
        values.assign(std::uint64_t(n)*ld,0.);pivots.resize(spd?0:n);
        for(int i=0;i<n;++i)a.Each(i,[&](int j,double v){
            int r=inverse[i],c=inverse[j];
            if(spd){if(r>=c)values[r-c+ld*c]=v;}
            else values[kl+ku+r-c+ld*c]=v;});
        int info=0;
        if(spd)Lapack::Dpbtrf('L',n,kd,values.data(),ld,info);
        else Lapack::Dgbtrf(n,n,kl,ku,values.data(),ld,pivots.data(),info);
        Require(info==0,spd?"mapped consistent mass Cholesky failed":"mapped nonsymmetric band LU failed");
    }
    Vector Apply(const Vector &rhs,bool transpose=false) const
    {
        auto x=Permute(rhs);int info=0;
        if(spd)Lapack::Dpbtrs('L',n,kd,1,values.data(),ld,x.data(),n,info);
        else Lapack::Dgbtrs(transpose?'T':'N',n,kl,ku,1,values.data(),ld,pivots.data(),x.data(),n,info);
        Require(info==0,"mapped band triangular solve failed");return RecoverPermutation(x);
    }
    Vector Whiten(const Vector &nativeRow) const
    {
        Require(spd,"bound whitening requires consistent mass Cholesky");auto x=Permute(nativeRow);
        for(int i=0;i<n;++i)
        {long double s=x[i];for(int j=std::max(0,i-kd);j<i;++j)s-=(long double)L(i,j)*x[j];x[i]=s/L(i,i);}
        Maximum(x);return x; // y lives in band-factor coordinates.
    }
    Vector Recover(const Vector &y) const
    {
        Require(spd && y.size()==std::size_t(n),"invalid bound recovery");auto x=y;
        for(int i=n-1;i>=0;--i)
        {long double s=x[i];for(int j=i+1;j<=std::min(n-1,i+kd);++j)s-=(long double)L(j,i)*x[j];x[i]=s/L(i,i);}
        return RecoverPermutation(x);
    }
    int LowerBandwidth() const {return kl;}
    int UpperBandwidth() const {return ku;}
    std::uint64_t StorageBytes() const {return values.size()*sizeof(double)+pivots.size()*sizeof(int);}
    const std::vector<int> &Permutation() const {return perm;}
private:
    Vector Permute(const Vector &x) const
    {Require(x.size()==std::size_t(n),"mapped factor vector shape mismatch");Vector y(n);for(int i=0;i<n;++i)y[i]=x[perm[i]];return y;}
    Vector RecoverPermutation(const Vector &y) const
    {Vector x(n);for(int i=0;i<n;++i)x[perm[i]]=y[i];Maximum(x);return x;}
    double L(int i,int j) const {return values[i-j+ld*j];}
    int n,kl=0,ku=0,kd=0,ld=0;bool spd;
    std::vector<int> perm,inverse,pivots;Vector values;
};

inline Vector Multiply(const MappedMatrix &a,const Vector &x){return a.Apply(x);}
inline Vector Solve(const MappedMatrix &a,const Vector &b)
{return a.Sparse()?MappedBandFactor(a,a.MassMetric()).Apply(b):Solve(static_cast<const DenseMatrix&>(a),b);}
}
#endif
