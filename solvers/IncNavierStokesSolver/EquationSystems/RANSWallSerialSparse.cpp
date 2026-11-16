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
#include <atomic>
#include <cstdlib>
#ifdef NEKTAR_RANS_HAVE_SCOTCH_ORDERING
#include <scotch.h>
#endif
#ifdef NEKTAR_RANS_HAVE_EIGEN_SPARSE
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#endif

namespace Nektar::RANSWallAlgebra
{
namespace {
std::atomic_size_t factorizationCount{0};
// Optional bounded audit output; no coefficient, solve or acceptance change.
void SparseReceipt(const std::string &record)
{
    const char *path=std::getenv("NEKTAR_RANS_SPARSE_RECEIPTS");if(!path || !*path)return;
    std::ofstream out(path,std::ios::app);out<<record<<'\n';out.close();
    Require(out.good(),"could not save sparse numerical receipt");
}
}
std::size_t SparseFactorizationCount() { return factorizationCount.load(std::memory_order_relaxed); }
double &SparseMatrix::Entry(int row,int column)
{
    Require(!m_data && row>=0 && row<m_rows && column>=0 && column<m_columns,
            "invalid write to finalized sparse matrix");
    return m_hash?m_hashAssembly[row][column]:m_assembly[row][column];
}
void SparseMatrix::ClearRow(int row)
{ Require(!m_data && row>=0 && row<m_rows,"invalid sparse row replacement");if(m_hash)m_hashAssembly[row].clear();else m_assembly[row].clear(); }
double SparseMatrix::Get(int row,int column) const
{
    Require(row>=0 && row<m_rows && column>=0 && column<m_columns,"sparse coefficient index");
    if(!m_data)
    {
        if(m_hash){const auto it=m_hashAssembly[row].find(column);return it==m_hashAssembly[row].end()?0.:it->second;}
        const auto it=m_assembly[row].find(column);return it==m_assembly[row].end()?0.:it->second;
    }
    const auto begin=m_data->ids.begin()+m_data->ptr[row],end=m_data->ids.begin()+m_data->ptr[row+1];
    const auto it=std::lower_bound(begin,end,column);
    return it==end || *it!=column?0.:m_data->values[it-m_data->ids.begin()];
}
void SparseMatrix::Finalize() const
{
    if(m_data)return;
    auto data=std::make_shared<Data>();data->ptr.push_back(0);
    const auto append=[&](const auto &row)
    {
        for(auto [j,value]:row)if(value!=0.)
        { Require(std::isfinite(value),"nonfinite sparse entry");data->ids.push_back(j);data->values.push_back(value); }
        Require(data->values.size()<=20000000,"sparse matrix exceeds declared20M nonzero budget");
        data->ptr.push_back(data->ids.size());
    };
    if(m_hash)
    {
        // Each Entry += contribution still occurs in the original order on
        // the same binary64 value. Only lookup/storage differs. Sort final
        // column IDs so the canonical CSR and every subsequent action retain
        // exactly the original order, including explicit-zero elimination.
        std::vector<std::pair<int,double>> ordered;
        for(const auto &row:m_hashAssembly)
        {
            ordered.assign(row.begin(),row.end());
            std::sort(ordered.begin(),ordered.end(),[](const auto &a,const auto &b){return a.first<b.first;});
            append(ordered);
        }
    }
    else for(const auto &row:m_assembly)append(row);
    m_data=data;m_assembly.clear();m_assembly.shrink_to_fit();m_hashAssembly.clear();m_hashAssembly.shrink_to_fit();
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

std::vector<int> SparseScotchPermutation(const SparseMatrix &matrix,const std::vector<int> &groups)
{
#ifdef NEKTAR_RANS_HAVE_SCOTCH_ORDERING
    Require(matrix.Rows()>0 && matrix.Rows()==matrix.Columns(),"Scotch graph shape");
    const int n=matrix.Rows();Require(groups.empty() || groups.size()==n,"Scotch vector grouping shape");
    std::vector<int> labels=groups,node(n);if(labels.empty()){labels.resize(n);std::iota(labels.begin(),labels.end(),0);}
    std::sort(labels.begin(),labels.end());labels.erase(std::unique(labels.begin(),labels.end()),labels.end());
    const int ng=labels.size();
    for(int i=0;i<n;++i)node[i]=groups.empty()?i:std::lower_bound(labels.begin(),labels.end(),groups[i])-labels.begin();
    std::vector<std::vector<SCOTCH_Num>> neighbors(ng);
    for(int i=0;i<n;++i)matrix.Each(i,[&](int j,double)
    {if(node[i]!=node[j]){neighbors[node[i]].push_back(node[j]);neighbors[node[j]].push_back(node[i]);}});
    std::vector<SCOTCH_Num> pointers(ng+1),edges;
    for(int i=0;i<ng;++i)
    {
        auto &row=neighbors[i];std::sort(row.begin(),row.end());
        row.erase(std::unique(row.begin(),row.end()),row.end());
        edges.insert(edges.end(),row.begin(),row.end());pointers[i+1]=edges.size();
    }
    SCOTCH_Graph graph;Require(SCOTCH_graphInit(&graph)==0,"Scotch graph init");
    struct GraphExit { SCOTCH_Graph *g;~GraphExit(){SCOTCH_graphExit(g);} } graphExit{&graph};
    Require(SCOTCH_graphBuild(&graph,0,ng,pointers.data(),pointers.data()+1,
        nullptr,nullptr,edges.size(),edges.data(),nullptr)==0 && SCOTCH_graphCheck(&graph)==0,
        "invalid canonical union Scotch graph");
    SCOTCH_Strat strategy;Require(SCOTCH_stratInit(&strategy)==0,"Scotch strategy init");
    struct StrategyExit { SCOTCH_Strat *s;~StrategyExit(){SCOTCH_stratExit(s);} } strategyExit{&strategy};
    // Fixed Scotch 6.1.3 strategy builder: quality + safety, level zero,
    // imbalance .1. Its resolved strategy is recorded by the CT14 probe.
    Require(SCOTCH_stratGraphOrderBuild(&strategy,SCOTCH_STRATQUALITY|SCOTCH_STRATSAFETY,0,.1)==0,
        "Scotch ND0 strategy unavailable");
    std::cerr<<(groups.empty()?"ScotchND0V1":"ScotchVectorND0V1")<<" version="<<SCOTCH_VERSION<<'.'<<SCOTCH_RELEASE<<'.'<<SCOTCH_PATCHLEVEL
        <<" num_bytes="<<sizeof(SCOTCH_Num)<<" scalar_rows="<<n<<" graph_vertices="<<ng<<" directed_edges="<<edges.size()
        <<" seed=0 flags=QUALITY|SAFETY level=0 balance=0.1 strategy=";
    Require(SCOTCH_stratSave(&strategy,stderr)==0,"could not record Scotch strategy");
    std::cerr<<'\n';
    SCOTCH_randomSeed(0);SCOTCH_randomReset();
    std::vector<SCOTCH_Num> forward(ng),inverse(ng);
    Require(SCOTCH_graphOrder(&graph,&strategy,forward.data(),inverse.data(),nullptr,nullptr,nullptr)==0,
        "Scotch ND0 ordering failed");
    std::vector<bool> seen(ng,false);std::vector<std::vector<int>> members(ng);
    for(int i=0;i<n;++i)members[node[i]].push_back(i);
    for(int i=0;i<ng;++i)
    {
        Require(forward[i]>=0 && forward[i]<ng && !seen[forward[i]],"invalid Scotch forward permutation");
        seen[forward[i]]=true;
        Require(inverse[forward[i]]==i,"Scotch inverse permutation mismatch");
    }
    std::vector<int> result(n);int next=0;
    for(int j=0;j<ng;++j)for(int old:members[inverse[j]])result[old]=next++;
    Require(next==n,"Scotch expanded permutation coverage");
    return result;
#else
    Require(false,"ScotchND0V1 requires explicitly linked Scotch support");return {};
#endif
}

std::vector<int> SparseVectorCOLAMDPermutation(const SparseMatrix &matrix,const std::vector<int> &groups)
{
#ifdef NEKTAR_RANS_HAVE_EIGEN_SPARSE
    const int n=matrix.Rows();Require(n>0 && matrix.Columns()==n && groups.size()==n,"vector COLAMD group shape");
    auto labels=groups;std::sort(labels.begin(),labels.end());labels.erase(std::unique(labels.begin(),labels.end()),labels.end());
    std::vector<int> node(n);const int ng=labels.size();std::vector<std::vector<int>> members(ng);
    for(int i=0;i<n;++i){node[i]=std::lower_bound(labels.begin(),labels.end(),groups[i])-labels.begin();members[node[i]].push_back(i);}
    using Matrix=Eigen::SparseMatrix<double,Eigen::ColMajor,int>;
    std::vector<Eigen::Triplet<double,int>> entries;entries.reserve(matrix.NonZeros());
    for(int i=0;i<n;++i)matrix.Each(i,[&](int j,double){entries.emplace_back(node[i],node[j],1.);});
    Matrix quotient(ng,ng);quotient.setFromTriplets(entries.begin(),entries.end(),[](double,double){return 1.;});quotient.makeCompressed();
    Eigen::PermutationMatrix<Eigen::Dynamic,Eigen::Dynamic,int> perm;
    Eigen::COLAMDOrdering<int>{}(quotient,perm);
    // Eigen Ordering.h stores old_to_new in indices(), unlike a list of
    // original columns in their new order. Expand each native vector group.
    std::vector<int> inverse(ng,-1),result(n);int next=0;
    for(int i=0;i<ng;++i){const int j=perm.indices()[i];Require(j>=0 && j<ng && inverse[j]<0,"invalid vector COLAMD permutation");inverse[j]=i;}
    for(int j=0;j<ng;++j)for(int i:members[inverse[j]])result[i]=next++;
    Require(next==n,"vector COLAMD expansion coverage");return result;
#else
    Require(false,"vector COLAMD needs explicitly enabled Eigen");return {};
#endif
}

struct SparseFactor::Impl
{
    SparseMatrix matrix;bool mass;SparsePredictorOrdering ordering;
    std::vector<int> oldToNew;
    std::size_t id=0;
#ifdef NEKTAR_RANS_HAVE_EIGEN_SPARSE
    using Matrix=Eigen::SparseMatrix<double,Eigen::ColMajor,int>;
    Matrix coefficients;
    Eigen::SparseLU<Matrix,Eigen::COLAMDOrdering<int>> lu;
    Eigen::SparseLU<Matrix,Eigen::AMDOrdering<int>> amd;
    Eigen::SparseLU<Matrix,Eigen::NaturalOrdering<int>> natural;
    Eigen::SimplicialLDLT<Matrix> ldlt;
#endif
    std::size_t nonzeros=0;
    Impl(const SparseMatrix &m,bool positive,SparsePredictorOrdering order):matrix(m),mass(positive),ordering(order){}
};
SparseFactor::SparseFactor(const SparseMatrix &m,bool positiveMass,SparsePredictorOrdering ordering,const std::vector<int> &groups):m_impl(std::make_shared<Impl>(m,positiveMass,ordering))
{
    m_impl->id=1+factorizationCount.fetch_add(1,std::memory_order_relaxed);
    RANSWallTiming::Scope timing("sparse.factorize");
    Require(m.Rows()>0 && m.Rows()==m.Columns(),"sparse factor shape");m.Finalize();
    Require(!positiveMass || ordering==SparsePredictorOrdering::COLAMD,"predictor ordering cannot change mass LDLT");
    Require((ordering==SparsePredictorOrdering::ScotchVectorND0V1 || ordering==SparsePredictorOrdering::ScotchVectorNDKeepV1 || ordering==SparsePredictorOrdering::VectorCOLAMDV1)==!groups.empty(),"explicit native grouping required only for vector orderings");
#ifdef NEKTAR_RANS_HAVE_EIGEN_SPARSE
    using Row=Eigen::SparseMatrix<double,Eigen::RowMajor,int>;
    Eigen::Map<const Row> map(m.Rows(),m.Columns(),m.NonZeros(),m.Pointers().data(),m.Indices().data(),m.Values().data());
    if(ordering==SparsePredictorOrdering::ScotchND0V1 || ordering==SparsePredictorOrdering::ScotchVectorND0V1 || ordering==SparsePredictorOrdering::ScotchVectorNDKeepV1 || ordering==SparsePredictorOrdering::VectorCOLAMDV1)
    {
        m_impl->oldToNew=ordering==SparsePredictorOrdering::VectorCOLAMDV1?SparseVectorCOLAMDPermutation(m,groups):SparseScotchPermutation(m,groups);
        std::vector<Eigen::Triplet<double,int>> entries;entries.reserve(m.NonZeros());
        for(int i=0;i<m.Rows();++i)m.Each(i,[&](int j,double v)
        {entries.emplace_back(m_impl->oldToNew[i],m_impl->oldToNew[j],v);});
        m_impl->coefficients.resize(m.Rows(),m.Columns());
        m_impl->coefficients.setFromTriplets(entries.begin(),entries.end());
        Require(m_impl->coefficients.nonZeros()==m.NonZeros(),"Scotch changed canonical nonzeros");
        for(int i=0;i<m.Rows();++i)m.Each(i,[&](int j,double v)
        {Require(m_impl->coefficients.coeff(m_impl->oldToNew[i],m_impl->oldToNew[j])==v,
            "Scotch permuted matrix differs from canonical coefficient");});
    }
    else m_impl->coefficients=map;
    if(positiveMass)
    {
        m_impl->ldlt.compute(m_impl->coefficients);
        Require(m_impl->ldlt.info()==Eigen::Success && m_impl->ldlt.vectorD().minCoeff()>0,"nonpositive/singular consistent sparse mass");
        m_impl->nonzeros=m_impl->ldlt.matrixL().nestedExpression().nonZeros();
    }
    else
    {
        const auto factor=[&](auto &lu)
        {
            lu.analyzePattern(m_impl->coefficients);lu.factorize(m_impl->coefficients);
            Require(lu.info()==Eigen::Success,"singular sparse stage/preconditioner");
            m_impl->nonzeros=lu.nnzL()+lu.nnzU();
        };
        if(ordering==SparsePredictorOrdering::COLAMD)factor(m_impl->lu);
        else if(ordering==SparsePredictorOrdering::AMD)factor(m_impl->amd);
        else
        {
            // Retain the supplied ND order. In installed Eigen this flag
            // disables column-etree postordering and selects heap-relaxed
            // supernodes only; generic LU arithmetic and pivot threshold stay.
            m_impl->natural.isSymmetric(ordering==SparsePredictorOrdering::ScotchVectorNDKeepV1);
            factor(m_impl->natural);
        }
    }
    Require(m_impl->nonzeros<=40000000,"sparse factor exceeds declared40M fill budget: rows="+
        std::to_string(m.Rows())+", matrix_nnz="+std::to_string(m.NonZeros())+", factor_nnz="+
        std::to_string(m_impl->nonzeros)+", mass="+std::to_string(positiveMass));
    SparseReceipt("{\"event\":\"factor\",\"id\":"+std::to_string(m_impl->id)+
        ",\"rows\":"+std::to_string(m.Rows())+",\"matrix_nnz\":"+std::to_string(m.NonZeros())+
        ",\"factor_nnz\":"+std::to_string(m_impl->nonzeros)+",\"mass\":"+(positiveMass?"true":"false")+
        ",\"ordering\":\""+(positiveMass?std::string("LDLT"):ordering==SparsePredictorOrdering::VectorCOLAMDV1?std::string("VectorCOLAMDV1"):ordering==SparsePredictorOrdering::ScotchVectorNDKeepV1?std::string("ScotchVectorNDKeepV1"):
            ordering==SparsePredictorOrdering::ScotchVectorND0V1?
            std::string("ScotchVectorND0V1"):ordering==SparsePredictorOrdering::ScotchND0V1?
            std::string("ScotchND0V1"):ordering==SparsePredictorOrdering::COLAMD?std::string("COLAMD"):
            ordering==SparsePredictorOrdering::AMD?std::string("AMD"):std::string("Natural"))+"\"}");
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
        Eigen::Map<const Eigen::VectorXd> original(r.data(),r.size());Eigen::VectorXd b=original,x;
        if(!m_impl->oldToNew.empty())
            for(int i=0;i<r.size();++i)b[m_impl->oldToNew[i]]=r[i];
        if(m_impl->mass)x=m_impl->ldlt.solve(b);
        else if(m_impl->ordering==SparsePredictorOrdering::COLAMD)x=m_impl->lu.solve(b);
        else if(m_impl->ordering==SparsePredictorOrdering::AMD)x=m_impl->amd.solve(b);
        else x=m_impl->natural.solve(b);
        Require(x.allFinite(),"nonfinite sparse solution");Vector result(x.data(),x.data()+x.size());
        if(!m_impl->oldToNew.empty())
            for(int i=0;i<r.size();++i)result[i]=x[m_impl->oldToNew[i]];
        return result;
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
            if(std::getenv("NEKTAR_RANS_SPARSE_RECEIPTS"))
            {
                std::ostringstream receipt;receipt<<std::setprecision(17)<<"{\"event\":\"solve\",\"id\":"<<m_impl->id
                    <<",\"rhs_maximum\":"<<Maximum(rhs)<<",\"true_residual\":"<<defect<<",\"bound\":"<<bound
                    <<",\"refinements\":"<<refinement<<'}';SparseReceipt(receipt.str());
            }
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
        for(auto &v:residual)
        {
            v=-v;
        }
        const auto correction=solve(residual);
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
