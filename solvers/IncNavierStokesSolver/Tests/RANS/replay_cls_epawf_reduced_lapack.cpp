// SPDX-License-Identifier: MIT
// Standalone Linux affine-model replay using the producing executable's
// LAPACK/BLAS libraries. This does not assemble governing-equation fields.
// Use a zero-multiplier packet and independently verified orthogonal field
// columns. A positive diagnostic penalty changes only this saved small model.
#include <IncNavierStokesSolver/EquationSystems/RANSWallReducedModel.h>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>

namespace R=Nektar::RANSWallReduced;
extern "C" {
void dgetrf_(const int *,const int *,double *,const int *,int *,int *);
void dgetrs_(const char *,const int *,const int *,const double *,const int *,
             const int *,double *,const int *,int *,std::size_t);
}
class Factor
{
    int n;
    R::Vector lu;
    std::vector<int> pivots;
public:
    explicit Factor(const R::Matrix &a):n(a.size()),lu(n*n),pivots(n)
    {
        if(n<=0 || std::any_of(a.begin(),a.end(),[&](const auto &r){return r.size()!=a.size();}))
            throw std::logic_error("factor shape");
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)lu[i+n*j]=a[i][j];
        int info;dgetrf_(&n,&n,lu.data(),&n,pivots.data(),&info);
        if(info<0)throw std::logic_error("factor argument");
        if(info)throw std::runtime_error("singular model");
    }
    R::Vector Apply(R::Vector rhs) const
    {
        if(rhs.size()!=static_cast<std::size_t>(n))throw std::logic_error("rhs shape");
        auto pivot=pivots;
        int info,one=1;char trans='N';
        dgetrs_(&trans,&n,&one,lu.data(),&n,pivot.data(),rhs.data(),&n,&info,1);
        if(info<0)throw std::logic_error("solve argument");
        if(info || std::any_of(rhs.begin(),rhs.end(),[](double v){return !std::isfinite(v);}))
            throw std::runtime_error("model solve failed");
        return rhs;
    }
};
struct Cut {R::Vector a;double b;};
int main(int argc,char **argv)
{
    if(argc!=5)throw std::logic_error("zero_multiplier_packet radius step_budget diagnostic_penalty");
    std::ifstream stream(argv[1],std::ios::binary);
    std::uint64_t shape[2]{};
    stream.read(reinterpret_cast<char *>(shape),sizeof(shape));
    if(!stream.good() || !shape[0] || shape[0]>512 || shape[1]>1000000)
        throw std::logic_error("packet shape");
    const int n=shape[0];
    auto read=[&](R::Vector &v){stream.read(reinterpret_cast<char *>(v.data()),8*v.size());};
    R::Matrix h(n,R::Vector(n)),columns(n,R::Vector(n));
    for(auto &row:h)read(row);
    R::Vector g(n),initial(n);read(g);read(initial);
    std::vector<Cut> cuts;
    for(std::uint64_t i=0;i<shape[1];++i){Cut c{R::Vector(n),0.};read(c.a);stream.read(reinterpret_cast<char *>(&c.b),8);cuts.push_back(std::move(c));}
    if(!stream.good() || stream.peek()!=std::char_traits<char>::eof() || R::Length(initial)!=0.)
        throw std::logic_error("packet payload");
    const double radius=std::stod(argv[2]),penalty=std::stod(argv[4]);
    const int budget=std::stoi(argv[3]);
    for(int i=0;i<n;++i){columns[i][i]=1.;h[i][i]+=penalty;}
    const auto result=R::SolveRadiusQuadratic<Factor>(h,g,cuts,initial,columns,radius,budget);
    const auto &q=result.quadratic;
    std::cout<<std::setprecision(17)<<"{\"solved\":"<<(q.solved?"true":"false")
        <<",\"failure\":\""<<q.failure<<"\",\"dimension\":"<<n
        <<",\"cuts\":"<<cuts.size()<<",\"pruned\":"<<result.ballRedundantCuts
        <<",\"norm\":"<<result.norm<<",\"raw_norm\":"<<result.rawNorm
        <<",\"radius_multiplier\":"<<result.multiplier
        <<",\"unconstrained_available\":"<<(result.unconstrainedAvailable?"true":"false")
        <<",\"unconstrained_failure\":\""<<result.unconstrainedFailure<<"\""
        <<",\"failed_positive_queries\":"<<result.failedPositiveQueries
        <<",\"diagnostic_penalty\":"<<penalty<<",\"steps\":"<<result.steps
        <<",\"solves\":"<<result.solves<<",\"primal\":"<<q.primal
        <<",\"stationarity\":"<<q.stationarity<<",\"complementarity\":"<<q.complementarity
        <<",\"radius_complementarity\":"<<result.complementarity
        <<",\"scope\":\"LAPACK affine model replay in the verified orthogonal metric; not original nonlinear convergence\",\"solution\":[";
    for(int i=0;i<n;++i)std::cout<<(i?",":"")<<q.solution[i];
    std::cout<<"]}\n";
    return q.solved?0:1;
}
