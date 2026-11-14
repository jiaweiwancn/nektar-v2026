///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallSerialSparse.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <LibUtilities/LinearAlgebra/NekLinSysIter.h>
#include <fstream>
#include <numeric>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#ifdef NEKTAR_RANS_HAVE_EIGEN_SPARSE
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#endif

namespace Nektar::RANSWallAlgebra
{
double &SparseMatrix::Entry(int row,int column)
{
    Require(!m_data && row>=0 && row<m_rows && column>=0 && column<m_columns,
            "invalid write to finalized sparse matrix");
    return m_assembly[row][column];
}
void SparseMatrix::ClearRow(int row)
{ Require(!m_data && row>=0 && row<m_rows,"invalid sparse row replacement");m_assembly[row].clear(); }
double SparseMatrix::Get(int row,int column) const
{
    Require(row>=0 && row<m_rows && column>=0 && column<m_columns,"sparse coefficient index");
    if(!m_data){const auto it=m_assembly[row].find(column);return it==m_assembly[row].end()?0.:it->second;}
    const auto begin=m_data->ids.begin()+m_data->ptr[row],end=m_data->ids.begin()+m_data->ptr[row+1];
    const auto it=std::lower_bound(begin,end,column);
    return it==end || *it!=column?0.:m_data->values[it-m_data->ids.begin()];
}
void SparseMatrix::Finalize() const
{
    if(m_data)return;
    auto data=std::make_shared<Data>();data->ptr.push_back(0);
    for(const auto &row:m_assembly)
    {
        for(auto [j,value]:row)if(value!=0.)
        { Require(std::isfinite(value),"nonfinite sparse entry");data->ids.push_back(j);data->values.push_back(value); }
        Require(data->values.size()<=20000000,"sparse matrix exceeds declared20M nonzero budget");
        data->ptr.push_back(data->ids.size());
    }
    m_data=data;m_assembly.clear();m_assembly.shrink_to_fit();
}
Vector SparseMatrix::Apply(const Vector &x) const
{
    Require(x.size()==static_cast<std::size_t>(m_columns),"sparse action shape");Vector y(m_rows);
    for(int i=0;i<m_rows;++i){long double sum=0;Each(i,[&](int j,double v){sum+=(long double)v*x[j];});y[i]=sum;}
    return y;
}
SparseMatrix SparseMatrix::Block(int r0,int r1,int c0,int c1) const
{
    Require(r0>=0 && r1<=m_rows && r0<=r1 && c0>=0 && c1<=m_columns && c0<=c1,"sparse block shape");
    SparseMatrix out(r1-r0,c1-c0);
    for(int i=r0;i<r1;++i)Each(i,[&](int j,double v){if(j>=c0 && j<c1)out.Entry(i-r0,j-c0)=v;});
    out.Finalize();return out;
}
SparseMatrix SparseMatrix::Restrict(const std::vector<int> &ids) const
{
    Require(m_rows==m_columns,"nonsquare restricted mass");std::vector<int> inverse(m_rows,-1);
    for(int i=0;i<ids.size();++i){Require(ids[i]>=0 && ids[i]<m_rows && inverse[ids[i]]<0,"invalid sparse free IDs");inverse[ids[i]]=i;}
    SparseMatrix out(ids.size(),ids.size());
    for(int i=0;i<ids.size();++i)Each(ids[i],[&](int j,double v){if(inverse[j]>=0)out.Entry(i,inverse[j])=v;});
    out.Finalize();return out;
}
DenseMatrix SparseMatrix::DenseForVerification(int limit) const
{
    Require(m_rows<=limit && m_columns<=limit,"optional dense sparse-export limit exceeded");
    DenseMatrix out(m_rows,Vector(m_columns));for(int i=0;i<m_rows;++i)Each(i,[&](int j,double v){out[i][j]=v;});return out;
}
void SparseMatrix::Write(const std::string &file) const
{
    Finalize();std::ofstream out(file,std::ios::binary);out.write("RANSCSR1",8);
    for(std::uint64_t n:{std::uint64_t(m_rows),std::uint64_t(m_columns),std::uint64_t(NonZeros())})out.write(reinterpret_cast<const char*>(&n),8);
    for(std::uint64_t n:m_data->ptr)out.write(reinterpret_cast<const char*>(&n),8);
    for(std::uint32_t n:m_data->ids)out.write(reinterpret_cast<const char*>(&n),4);
    out.write(reinterpret_cast<const char*>(m_data->values.data()),m_data->values.size()*8);
    out.close();Require(out.good(),"could not save CSR snapshot");
}

struct SparseFactor::Impl
{
    SparseMatrix matrix;bool mass;
#ifdef NEKTAR_RANS_HAVE_EIGEN_SPARSE
    using Matrix=Eigen::SparseMatrix<double,Eigen::ColMajor,int>;
    Matrix coefficients;
    Eigen::SparseLU<Matrix,Eigen::COLAMDOrdering<int>> lu;
    Eigen::SimplicialLDLT<Matrix> ldlt;
#endif
    std::size_t nonzeros=0;
    Impl(const SparseMatrix &m,bool positive):matrix(m),mass(positive){}
};
SparseFactor::SparseFactor(const SparseMatrix &m,bool positiveMass):m_impl(std::make_shared<Impl>(m,positiveMass))
{
    RANSWallTiming::Scope timing("sparse.factorize");
    Require(m.Rows()>0 && m.Rows()==m.Columns(),"sparse factor shape");m.Finalize();
#ifdef NEKTAR_RANS_HAVE_EIGEN_SPARSE
    using Row=Eigen::SparseMatrix<double,Eigen::RowMajor,int>;
    Eigen::Map<const Row> map(m.Rows(),m.Columns(),m.NonZeros(),m.Pointers().data(),m.Indices().data(),m.Values().data());
    m_impl->coefficients=map;
    if(positiveMass)
    {
        m_impl->ldlt.compute(m_impl->coefficients);
        Require(m_impl->ldlt.info()==Eigen::Success && m_impl->ldlt.vectorD().minCoeff()>0,"nonpositive/singular consistent sparse mass");
        m_impl->nonzeros=m_impl->ldlt.matrixL().nestedExpression().nonZeros();
    }
    else
    {
        m_impl->lu.analyzePattern(m_impl->coefficients);m_impl->lu.factorize(m_impl->coefficients);
        Require(m_impl->lu.info()==Eigen::Success,"singular sparse stage/preconditioner");
        m_impl->nonzeros=m_impl->lu.nnzL()+m_impl->lu.nnzU();
    }
    Require(m_impl->nonzeros<=30000000,"sparse factor exceeds declared30M fill budget: rows="+
        std::to_string(m.Rows())+", matrix_nnz="+std::to_string(m.NonZeros())+", factor_nnz="+
        std::to_string(m_impl->nonzeros)+", mass="+std::to_string(positiveMass));
#else
    Require(false,"AcceptedSerialSparseV1 requires an explicitly Eigen-enabled isolated build");
#endif
}
Vector SparseFactor::Apply(const Vector &rhs) const
{
    RANSWallTiming::Scope timing("sparse.solve_validate");
    Require(rhs.size()==static_cast<std::size_t>(m_impl->matrix.Rows()),"sparse solve shape");
#ifdef NEKTAR_RANS_HAVE_EIGEN_SPARSE
    const auto solve=[&](const Vector &r)
    {
        Eigen::Map<const Eigen::VectorXd> b(r.data(),r.size());Eigen::VectorXd x;
        if(m_impl->mass)x=m_impl->ldlt.solve(b);else x=m_impl->lu.solve(b);
        Require(x.allFinite(),"nonfinite sparse solution");return Vector(x.data(),x.data()+x.size());
    };
    Vector result=solve(rhs);const double bound=2e-12*Maximum(rhs)+1e-25;double initial=0.;
    for(int refinement=0;refinement<=3;++refinement)
    {
        // Extended accumulation checks the ORIGINAL unscaled sparse system.
        // Refinement reuses its factor and never changes a coefficient, mask,
        // physical source, governing stage or acceptance tolerance.
        auto residual=m_impl->matrix.Apply(result);for(int i=0;i<rhs.size();++i)residual[i]-=rhs[i];
        const double defect=Maximum(residual);if(refinement==0)initial=defect;
        if(defect<=bound)
        {
            if(refinement)std::cerr<<std::setprecision(17)<<"RANS sparse refinement {\"rows\":"<<rhs.size()
                <<",\"mass\":"<<(m_impl->mass?"true":"false")<<",\"steps\":"<<refinement
                <<",\"initial_residual\":"<<initial<<",\"final_residual\":"<<defect<<",\"bound\":"<<bound<<"}\n";
            return result;
        }
        if(refinement==3)
        {
            std::ostringstream context;context<<std::setprecision(17)<<"sparse direct true residual exceeds contract: rows="<<rhs.size()
                <<", mass="<<m_impl->mass<<", nnz="<<m_impl->matrix.NonZeros()<<", rhs="<<Maximum(rhs)
                <<", solution="<<Maximum(result)<<", initial="<<initial<<", final="<<defect<<", bound="<<bound<<", refinements=3";
            Require(false,context.str());
        }
        RANSWallTiming::Scope refinementTiming("sparse.residual_refinement");
        for(auto &v:residual)v=-v;const auto correction=solve(residual);
        for(int i=0;i<result.size();++i)result[i]+=correction[i];
    }
    Require(false,"unreachable sparse refinement exit");return {};
#else
    Require(false,"sparse factor not compiled");return {};
#endif
}
std::size_t SparseFactor::FactorNonZeros() const { return m_impl->nonzeros; }

SparseKrylovResult SparsePressureGMRES(const LibUtilities::SessionReaderSharedPtr &session,
    const LibUtilities::CommSharedPtr &comm,const Vector &rhs,
    const std::function<Vector(const Vector&)> &action,
    const std::function<Vector(const Vector&)> &preconditioner)
{
    RANSWallTiming::Scope timing("sparse.pressure_gmres");
    Require(comm->GetSize()==1 && !rhs.empty(),"sparse pressure is serial only");
    struct Adapter
    {
        std::function<Vector(const Vector&)> lhs,pre;
        using In=LibUtilities::NekSysOperators::InArrayType;
        using Out=LibUtilities::NekSysOperators::OutArrayType;
        void Lhs(In &x,Out &y,const bool&) { auto v=lhs(Vector(x.data(),x.data()+x.size()));std::copy(v.begin(),v.end(),y.data()); }
        void Pre(In &x,Out &y,const bool&) { auto v=pre(Vector(x.data(),x.data()+x.size()));std::copy(v.begin(),v.end(),y.data()); }
    } adapter{action,preconditioner};
    SparseKrylovResult out;out.solution.assign(rhs.size(),0.);const double scale=Maximum(rhs);if(scale==0)return out;
    // Nektar's small-RHS magnitude convention otherwise creates an absolute
    // floor. Normalize the linear system, then verify its UNSCALED residual.
    Array<OneD,NekDouble> b(rhs.size()),x(rhs.size(),0.);Array<OneD,int> unique(rhs.size(),1);
    for(int i=0;i<rhs.size();++i)b[i]=rhs[i]/scale;
    LibUtilities::NekSysKey key;key.m_NekLinSysMaxIterations=1000;key.m_LinSysMaxStorage=80;
    key.m_KrylovMaxHessMatBand=0;key.m_NekLinSysTolerance=2e-13/std::sqrt(double(rhs.size()));
    auto solver=LibUtilities::GetNekLinSysIterFactory().CreateInstance("GMRES",session,comm,rhs.size(),key);
    LibUtilities::NekSysOperators ops;ops.DefineNekSysLhsEval(&Adapter::Lhs,&adapter);ops.DefineNekSysPrecon(&Adapter::Pre,&adapter);
    solver->SetSysOperators(ops);solver->SetUniversalUniqueMap(unique);
    out.iterations=solver->SolveSystem(rhs.size(),b,x,0);
    for(int i=0;i<rhs.size();++i)out.solution[i]=scale*x[i];
    auto residual=action(out.solution);for(int i=0;i<rhs.size();++i)residual[i]-=rhs[i];out.trueResidual=Maximum(residual);
    Require(out.trueResidual<=2e-12*scale+1e-25,"matrix-free pressure true residual exceeds contract");return out;
}
}
