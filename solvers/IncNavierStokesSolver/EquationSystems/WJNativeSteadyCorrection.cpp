// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeSteadyCorrection.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <numeric>

namespace Nektar::WJNative
{
namespace
{
using Clock=std::chrono::steady_clock;
using Vec=SteadyVector;
using Triplet=Eigen::Triplet<double>;

struct Timer
{
    explicit Timer(double &seconds):elapsed(seconds),start(Clock::now()){}
    ~Timer(){elapsed+=std::chrono::duration<double>(Clock::now()-start).count();}
    double &elapsed;
    Clock::time_point start;
};
}

SteadySparseInverse::SteadySparseInverse(std::vector<std::vector<int>> rowsByColumn,
    SparseMatrix exact,Vec unknownScale,Vec rowScale)
    :m_support(std::move(rowsByColumn)),m_exact(std::move(exact)),
     m_unknownScale(std::move(unknownScale)),m_rowScale(std::move(rowScale))
{
    const int n=m_unknownScale.size();
    if(m_support.size()!=std::size_t(n)||m_rowScale.size()!=n||m_exact.rows()!=n||m_exact.cols()!=n)
        throw std::invalid_argument("sparse steady correction inventory differs");
    std::vector<int> order(n);
    std::iota(order.begin(),order.end(),0);
    for(auto &rows:m_support)
    {
        std::sort(rows.begin(),rows.end());
        rows.erase(std::unique(rows.begin(),rows.end()),rows.end());
        if(!rows.empty()&&(rows.front()<0||rows.back()>=n))
            throw std::invalid_argument("sparse steady correction support contains an invalid row");
    }
    std::stable_sort(order.begin(),order.end(),[&](int lhs,int rhs){
        return m_support[lhs].size()>m_support[rhs].size();
    });
    std::vector<std::vector<int>> colorsByRow(n);
    std::vector<int> unavailable(n,-1);
    for(int col:order)
    {
        if(m_support[col].empty())continue;
        for(int row:m_support[col])
            for(int color:colorsByRow[row])unavailable[color]=col;
        int color=0;
        while(unavailable[color]==col)++color;
        if(color==int(m_colors.size()))m_colors.emplace_back();
        m_colors[color].push_back(col);
        for(int row:m_support[col])colorsByRow[row].push_back(color);
    }
}

void SteadySparseInverse::Prepare(const Vec &base,const Action &residual,const Action &correction,
    const SteadyKrylovOptions &options,const std::function<void()> &deadline,SteadySparseWork &work)
{
    const int n=base.size();
    work.colors=m_colors.size();
    std::vector<Triplet> values,shiftValues;
    {
        Timer timer(work.assemblySeconds);
        std::size_t entries=m_exact.nonZeros();
        for(const auto &rows:m_support)entries+=rows.size();
        values.reserve(entries);
        for(int col=0;col<m_exact.outerSize();++col)
            for(SparseMatrix::InnerIterator it(m_exact,col);it;++it)
                if(it.value()!=0.)
                    values.emplace_back(it.row(),it.col(),m_rowScale[it.row()]*it.value()*m_unknownScale[col]);
        const double scale=std::max(1.,base.cwiseQuotient(m_unknownScale).cwiseAbs().maxCoeff());
        for(const auto &columns:m_colors)
        {
            Vec direction=Vec::Zero(n);
            std::vector<int> owner(n,-1);
            for(int col:columns)
            {
                direction[col]=m_unknownScale[col];
                for(int row:m_support[col])
                {
                    if(owner[row]!=-1)throw std::logic_error("sparse correction color has conflicting rows");
                    owner[row]=col;
                }
            }
            Vec differential;
            double relative=options.fdRelativeStep;
            for(int reduction=0;reduction<=options.maximumFDReductions;++reduction)
            {
                deadline();
                const double step=relative*scale;
                const Vec plus=base+step*direction,minus=base-step*direction;
                if(!plus.allFinite()||!minus.allFinite())
                    throw SteadyCorrectionFailure("SPARSE_FD_ARITHMETIC: nonfinite colored displacement");
                try
                {
                    ++work.residualEvaluations;
                    const Vec high=residual(plus);
                    deadline();
                    ++work.residualEvaluations;
                    const Vec low=residual(minus);
                    differential=(high-low)/(2.*step);
                    break;
                }
                catch(const SteadyInadmissibleTrial &error)
                {
                    ++work.rejectedSamples;
                    work.lastRejection=error.what();
                    relative*=options.fdReduction;
                    if(relative<options.fdMinimumRelativeStep)break;
                }
            }
            if(differential.size()==0)
                throw SteadyCorrectionFailure("SPARSE_FD_INFEASIBLE: "+work.lastRejection);
            if(!differential.allFinite())
                throw SteadyCorrectionFailure("SPARSE_FD_ARITHMETIC: nonfinite colored action");

            differential-=m_exact*direction;
            const Vec scaled=differential.cwiseProduct(m_rowScale);
            const Vec shift=correction?correction(direction).cwiseProduct(m_rowScale).eval():Vec();
            const double tolerance=2e-8*std::max(1.,scaled.cwiseAbs().maxCoeff());
            for(int row=0;row<n;++row)
            {
                const int col=owner[row];
                const double extra=correction?shift[row]:0.;
                if(col<0)
                {
                    const double observed=std::max(std::abs(scaled[row]),std::abs(extra));
                    work.maximumUnassignedAction=std::max(work.maximumUnassignedAction,observed);
                    if(observed>tolerance)
                        throw SteadyCorrectionFailure("SPARSE_SUPPORT: nonzero action outside the declared row support at row "+
                                                      std::to_string(row));
                    continue;
                }
                values.emplace_back(row,col,scaled[row]);
                if(extra!=0.)shiftValues.emplace_back(row,col,extra);
            }
        }
        m_original.resize(n,n);
        m_original.setFromTriplets(values.begin(),values.end());
        m_original.makeCompressed();
        work.nonzeros=m_original.nonZeros();
    }
    deadline();
    {
        Timer timer(work.factorSeconds);
        if(shiftValues.empty())
            m_factor.compute(m_original);
        else
        {
            SparseMatrix shift(n,n);
            shift.setFromTriplets(shiftValues.begin(),shiftValues.end());
            SparseMatrix total=m_original+shift;
            total.makeCompressed();
            m_factor.compute(total);
        }
        if(m_factor.info()!=Eigen::Success)
            throw SteadyCorrectionFailure("SPARSE_FACTOR: "+m_factor.lastErrorMessage());
        work.factorNonzeros=m_factor.nnzL()+m_factor.nnzU();
    }
    deadline();
}

Vec SteadySparseInverse::OriginalAction(const Vec &direction)const
{
    return (m_original*direction.cwiseQuotient(m_unknownScale)).cwiseQuotient(m_rowScale);
}

Vec SteadySparseInverse::Solve(const Vec &rhs)const
{
    Vec correction=m_factor.solve(rhs.cwiseProduct(m_rowScale));
    if(m_factor.info()!=Eigen::Success||!correction.allFinite())
        throw SteadyCorrectionFailure("SPARSE_BACKSOLVE: invalid correction");
    return correction.cwiseProduct(m_unknownScale);
}
}
