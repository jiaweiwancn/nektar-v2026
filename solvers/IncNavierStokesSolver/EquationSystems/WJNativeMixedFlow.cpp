// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeMixedFlow.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeFailures.h>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/QR>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace Nektar::WJNative
{
namespace
{
using T=Eigen::Triplet<double>;
void Check(bool ok,const char*text){if(!ok)throw std::runtime_error(text);}
double Clock(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}

class WorkTimer
{
public:
    explicit WorkTimer(double *total) : m_total(total), m_start(total ? Clock() : 0.) {}
    ~WorkTimer()
    {
        if (m_total)
        {
            *m_total += Clock() - m_start;
        }
    }

private:
    double *m_total;
    double m_start;
};

template <typename Apply, typename Precondition, typename CountIteration>
Vec KrylovCorrection(const Vec &residual, int dimension, double projectedTolerance,
                     Apply apply, Precondition precondition, CountIteration countIteration)
{
    const double beta = residual.norm();
    RequireLinear(std::isfinite(beta) && beta > 0, "invalid Krylov residual");
    Mat v(residual.size(), dimension + 1), z(residual.size(), dimension);
    Mat h = Mat::Zero(dimension + 1, dimension);
    Vec rhs = Vec::Zero(dimension + 1), best;
    v.col(0) = residual / beta;
    rhs[0] = beta;
    int used = 0;
    for (int j = 0; j < dimension; ++j)
    {
        z.col(j) = precondition(v.col(j));
        Vec action = apply(z.col(j));
        for (int pass = 0; pass < 2; ++pass)
        {
            for (int i = 0; i <= j; ++i)
            {
                const double value = v.col(i).dot(action);
                h(i, j) += value;
                action.noalias() -= value * v.col(i);
            }
        }
        h(j + 1, j) = action.norm();
        if (h(j + 1, j) > 0)
        {
            v.col(j + 1) = action / h(j + 1, j);
        }
        const auto column = h.topLeftCorner(j + 2, j + 1);
        best = column.colPivHouseholderQr().solve(rhs.head(j + 2));
        RequireLinear(best.allFinite(), "nonfinite Krylov least-squares solve");
        countIteration();
        used = j + 1;
        if ((rhs.head(j + 2) - column * best).norm() < projectedTolerance ||
            h(j + 1, j) < 1e-28)
        {
            break;
        }
    }
    return z.leftCols(used) * best;
}

bool SameCompressedMatrix(const Sparse &a, const Sparse &b)
{
    return a.rows() == b.rows() && a.cols() == b.cols() && a.nonZeros() == b.nonZeros() &&
           a.isCompressed() && b.isCompressed() &&
           std::equal(a.outerIndexPtr(), a.outerIndexPtr() + a.outerSize() + 1, b.outerIndexPtr()) &&
           std::equal(a.innerIndexPtr(), a.innerIndexPtr() + a.nonZeros(), b.innerIndexPtr()) &&
           std::equal(a.valuePtr(), a.valuePtr() + a.nonZeros(), b.valuePtr());
}

Sparse Gather(const Sparse&a,const std::vector<int>&ids)
{
    std::vector<int>inverse(a.rows(),-1);for(int i=0;i<int(ids.size());++i)inverse[ids[i]]=i;
    std::vector<T>t;
    for(int c=0;c<a.outerSize();++c)for(Sparse::InnerIterator it(a,c);it;++it)
        if(inverse[it.row()]>=0&&inverse[it.col()]>=0)t.emplace_back(inverse[it.row()],inverse[it.col()],it.value());
    Sparse result(ids.size(),ids.size());result.setFromTriplets(t.begin(),t.end());return result;
}
Mat NativeBasis(const LocalRegions::ExpansionSharedPtr&e,int derivative)
{
    int nq=e->GetTotPoints(),nc=e->GetNcoeffs();Mat b(nq,nc);
    Array<OneD,NekDouble>coef(nc,0.),p(nq),dx(nq),dz(nq);
    for(int j=0;j<nc;++j){Vmath::Zero(nc,coef,1);coef[j]=1;e->BwdTrans(coef,p);
        if(derivative>=0)e->PhysDeriv(p,dx,dz);
        for(int q=0;q<nq;++q)b(q,j)=derivative<0?p[q]:derivative==0?dx[q]:dz[q];}
    return b;
}
Vec At(const LocalRegions::ExpansionSharedPtr&e,const Mat&basis,double xi,double eta)
{
    Array<OneD,NekDouble>point(2),p(e->GetTotPoints());point[0]=xi;point[1]=eta;Vec b(basis.cols());
    for(int j=0;j<basis.cols();++j){for(int q=0;q<basis.rows();++q)p[q]=basis(q,j);b[j]=e->StdPhysEvaluate(point,p);}return b;
}
std::pair<Vec,Vec>Gauss(int n)
{
    Vec x(n),w(n);
    for(int i=0;i<n;++i){double z=std::cos(std::acos(-1.)*(i+.75)/(n+.5)),d=0;
        for(int j=0;j<40;++j){double a=1,b=z;for(int k=2;k<=n;++k){double c=((2*k-1)*z*b-(k-1)*a)/k;a=b;b=c;}d=n*(z*b-a)/(z*z-1);double step=b/d;z-=step;if(std::abs(step)<2e-16)break;}
        x[i]=z;w[i]=2/((1-z*z)*d*d);}
    return{x,w};
}
}

MixedFlow::MixedFlow(const std::vector<MultiRegions::ContFieldSharedPtr>&fields,
                     Geometry geometry, int pressureCorrectors)
    :m_fields(fields),m_core({fields.at(0),fields.at(0),fields.at(2)},geometry),
     m_pressureCorrectors(pressureCorrectors)
{
    Check(fields.size()==3,"H27 MixedFlow requires serial u,v,p fields");
    Check(pressureCorrectors >= 0, "negative number of mixed pressure correctors");
    const int n=VelocitySize(),np=PressureSize();auto vm=fields[1]->GetLocalToGlobalMap();
    Check(vm->GetNumGlobalCoeffs()==n,"component global Q3 spaces have different dimensions");
    Check(fields[1]->GetExpSize()==int(Elements().size()),"component element inventories differ");
    m_vToU.assign(n,-1);m_vSign.assign(n,0.);std::vector<int>reverse(n,-1);
    for(int e=0;e<int(Elements().size());++e)
    {
        const auto&a=Elements()[e];auto v=fields[1]->GetExp(e);
        Check(v->GetGeom()->GetGlobalID()==a.geometry&&v->GetNcoeffs()==16,"component native geometry/basis ordering differs");
        Check(v->GetNumBases()==2,"component field is not two dimensional");
        for(int d=0;d<2;++d)Check(v->GetBasisNumModes(d)==4&&v->GetBasis(d)->GetBasisType()==a.velocity->GetBasis(d)->GetBasisType(),"component Q3 basis definitions differ");
        for(int j=0;j<16;++j){int l=fields[1]->GetCoeff_Offset(e)+j,id=vm->GetLocalToGlobalMap(l),u=a.u[j];double sign=vm->GetLocalToGlobalSign(l)/a.us[j];
            Check(std::abs(std::abs(sign)-1)<1e-14,"non-unit native component sign");
            Check(m_vToU[id]<0||(m_vToU[id]==u&&m_vSign[id]==sign),"inconsistent native v signed permutation");
            Check(reverse[u]<0||reverse[u]==id,"native component permutation is not bijective");
            m_vToU[id]=u;m_vSign[id]=sign;reverse[u]=id;}
        m_volume+=a.w.sum();
    }
    Check(std::all_of(m_vToU.begin(),m_vToU.end(),[](int v){return v>=0;}),"incomplete native v permutation");
    std::array<std::vector<bool>,2>dir{std::vector<bool>(n,false),std::vector<bool>(n,false)};
    for(int i=0;i<m_core.DirichletSize();++i)dir[0][i]=true;
    for(int i=0;i<vm->GetNumGlobalDirBndCoeffs();++i)dir[1][m_vToU[i]]=true;
    m_fullToFree.assign(2*n,-1);
    for(int c=0;c<2;++c)for(int i=0;i<n;++i)if(!dir[c][i]){m_free[c].push_back(i);m_fullToFree[c*n+i]=int(m_freeFull.size());m_freeFull.push_back(c*n+i);}
    Check(!m_free[0].empty()&&!m_free[1].empty()&&m_volume>0,"empty component increment space/volume");
    m_sharedRiesz=m_free[0]==m_free[1];
    for(int c=0;c<2;++c)if(c==0||!m_sharedRiesz){m_velocityRiesz[c]=std::make_unique<Eigen::SimplicialLDLT<Sparse>>();m_velocityRiesz[c]->compute(Gather(m_core.Mass(),m_free[c]));Check(m_velocityRiesz[c]->info()==Eigen::Success,"component mass Riesz factor failed");}
    m_pressureRiesz.compute(m_core.PressureMass());Check(m_pressureRiesz.info()==Eigen::Success,"pressure mass Riesz factor failed");
    const Vec one=PressureOne();Sparse effective=m_core.Gradient()-m_core.AWF();double norm2=0;
    for(int c=0;c<effective.outerSize();++c)for(Sparse::InnerIterator it(effective,c);it;++it)if(m_fullToFree[it.row()]>=0)norm2+=it.value()*it.value();
    Vec fullAction=effective*one;double action2=0;for(int id:m_freeFull)action2+=fullAction[id]*fullAction[id];
    Check(norm2>0&&one.norm()>0,"degenerate restricted pressure action");
    m_constantPressureAction=std::sqrt(action2/norm2)/one.norm();
    // A constant pressure is a gauge only when its ACTUAL restricted action
    // vanishes. Tangential Neumann data alone need not remove that nullspace.
    m_gauge=m_constantPressureAction<=1024*std::numeric_limits<double>::epsilon();
    const double cscale=m_core.AWF().norm()*one.norm();
    Check(cscale==0||(m_core.AWF()*one).norm()/cscale<1e-12,"AWF sampler changes constant pressure");
    Eigen::Index pivot=0;one.cwiseAbs().maxCoeff(&pivot);m_pressurePivot=int(pivot);
    m_scale=Vec::Ones(Size());
    // Unit U/L and acceleration reference scales; alpha/time mass is absent.
    for(int i=0;i<PressureOffset();++i)m_scale[i]=1/std::sqrt(m_core.Mass().coeff(m_freeFull[i]%n,m_freeFull[i]%n));
    for(int i=0;i<np;++i)m_scale[PressureOffset()+i]=1/std::sqrt(m_core.PressureMass().coeff(i,i));
    if(m_gauge)m_scale[Size()-1]=1/(.01*m_volume);
}

Sparse MixedFlow::System(double alpha,const Sparse&k,double awf)const
{
    const int n=VelocitySize(),offset=PressureOffset();std::vector<T>t;
    for(int c=0;c<k.outerSize();++c)for(Sparse::InnerIterator it(k,c);it;++it){int i=m_fullToFree[it.row()],j=m_fullToFree[it.col()];if(i>=0&&j>=0)t.emplace_back(i,j,it.value());}
    for(int component=0;component<2;++component)for(int c=0;c<m_core.Mass().outerSize();++c)for(Sparse::InnerIterator it(m_core.Mass(),c);it;++it){int i=m_fullToFree[component*n+it.row()],j=m_fullToFree[component*n+it.col()];if(i>=0&&j>=0)t.emplace_back(i,j,alpha*it.value());}
    Sparse g=m_core.Gradient()-awf*m_core.AWF();
    for(int c=0;c<g.outerSize();++c)for(Sparse::InnerIterator it(g,c);it;++it){int i=m_fullToFree[it.row()];if(i>=0)t.emplace_back(i,offset+it.col(),it.value());}
    for(int c=0;c<m_core.Divergence().outerSize();++c)for(Sparse::InnerIterator it(m_core.Divergence(),c);it;++it){int j=m_fullToFree[it.col()];if(j>=0)t.emplace_back(offset+it.row(),j,it.value());}
    if(m_gauge)for(int i=0;i<PressureSize();++i){t.emplace_back(offset+i,Size()-1,PressureMean()[i]);t.emplace_back(Size()-1,offset+i,PressureMean()[i]);}
    Sparse a(Size(),Size());a.setFromTriplets(t.begin(),t.end());return a;
}

Vec MixedFlow::VelocityDiagonal(double alpha, const Sparse &diffusion) const
{
    Vec diagonal(PressureOffset());
    const int n = VelocitySize();
    for (int i = 0; i < diagonal.size(); ++i)
    {
        const int id = m_freeFull[i];
        diagonal[i] = alpha * m_core.Mass().coeff(id % n, id % n) +
                      diffusion.coeff(id, id);
    }
    return diagonal;
}

Vec MixedFlow::ProjectNative(const Function&f,const MultiRegions::ContFieldSharedPtr&field)const
{
    Array<OneD,NekDouble>x(field->GetTotPoints()),z(x.size()),dummy(x.size()),v(x.size()),coef(field->GetNcoeffs());
    field->GetCoords(x,z,dummy);for(int i=0;i<int(x.size());++i)v[i]=f(x[i],z[i]);field->FwdTransBndConstrained(v,coef);
    for(int b=0;b<int(field->GetBndConditions().size());++b)if(field->GetBndConditions()[b]->GetBoundaryConditionType()==SpatialDomains::eDirichlet){auto edge=field->GetBndCondExpansions()[b];Array<OneD,NekDouble>ex(edge->GetTotPoints()),ez(ex.size()),ed(ex.size()),ev(ex.size());edge->GetCoords(ex,ez,ed);for(int i=0;i<int(ex.size());++i)ev[i]=f(ex[i],ez[i]);edge->FwdTransBndConstrained(ev,edge->UpdateCoeffs());}
    field->ImposeDirichletConditions(coef);field->UpdateCoeffs()=coef;field->AvgAssemble();field->GlobalToLocal();auto map=field->GetLocalToGlobalMap();Array<OneD,NekDouble>g(map->GetNumGlobalCoeffs());map->LocalToGlobal(field->GetCoeffs(),g,true);Vec result(g.size());for(int i=0;i<int(g.size());++i)result[i]=g[i];return result;
}
Vec MixedFlow::Project(const Function&f,int component)const
{
    Check(component>=0&&component<3,"invalid mixed component");Vec native=ProjectNative(f,m_fields[component]);if(component!=1)return native;
    Vec canonical(native.size());for(int j=0;j<int(m_vToU.size());++j)canonical[m_vToU[j]]=m_vSign[j]*native[j];return canonical;
}
Vec MixedFlow::NativeVelocity(const Vec&full,int component)const
{
    Check((component==0||component==1)&&full.size()==2*VelocitySize()+PressureSize(),"invalid velocity component/full state");Vec result(VelocitySize());
    if(component==0)return full.head(VelocitySize());
    for(int j=0;j<int(m_vToU.size());++j)result[j]=m_vSign[j]*full[VelocitySize()+m_vToU[j]];return result;
}
Vec MixedFlow::Lift(const Vec&full)const
{
    Check(full.size()==2*VelocitySize()+PressureSize(),"mixed full state size");Vec result=Vec::Zero(full.size());for(int i=0;i<2*VelocitySize();++i)if(m_fullToFree[i]<0)result[i]=full[i];return result;
}
Vec MixedFlow::Restrict(const Vec&full)const
{
    Check(full.size()==2*VelocitySize()+PressureSize(),"mixed full load/state size");
    Vec result=Vec::Zero(Size());for(int i=0;i<PressureOffset();++i)result[i]=full[m_freeFull[i]];result.segment(PressureOffset(),PressureSize())=full.tail(PressureSize());return result;
}
Vec MixedFlow::Expand(const Vec&free,const Vec&lift)const
{
    Vec result=lift;for(int i=0;i<PressureOffset();++i)result[m_freeFull[i]]=free[i];result.tail(PressureSize())=free.segment(PressureOffset(),PressureSize());return result;
}
Vec MixedFlow::LiftedRhs(double alpha,const Sparse&k,double awf,const Vec&load,const Vec&lift)const
{return Restrict(load-m_core.ApplyFull(alpha,k,awf,lift));}
double MixedFlow::VelocityDual(const Vec&load)const
{
    Check(load.size()>=2*VelocitySize(),"velocity dual load size");double norm=0;
    for(int c=0;c<2;++c){Vec v(m_free[c].size());for(int i=0;i<int(m_free[c].size());++i)v[i]=load[c*VelocitySize()+m_free[c][i]];Vec z=m_velocityRiesz[m_sharedRiesz?0:c]->solve(v);norm+=v.dot(z);}
    return std::sqrt(std::max(0.,norm)/m_volume);
}
SolveInfo MixedFlow::Residual(const Sparse&a,const Vec&b,const Vec&x)const
{
    Vec residual=a*x-b;
    if(m_gauge){residual.segment(PressureOffset(),PressureSize())-=PressureMean()*x[Size()-1];residual[Size()-1]=PressureMean().dot(x.segment(PressureOffset(),PressureSize()));}
    return RawResidualNorm(residual);
}
SolveInfo MixedFlow::RawResidualNorm(const Vec&residual)const
{
    Check(residual.size()==Size()&&residual.allFinite(),"invalid raw mixed residual");
    SolveInfo info;Vec full=Vec::Zero(2*VelocitySize()+PressureSize());for(int i=0;i<PressureOffset();++i)full[m_freeFull[i]]=residual[i];info.momentum=VelocityDual(full);
    info.continuity=ContinuityNorm(residual.segment(PressureOffset(),PressureSize()));
    info.gauge=m_gauge?std::abs(residual[Size()-1])/(.01*m_volume):0;return info;
}

double MixedFlow::ContinuityNorm(const Vec &residual) const
{
    const Vec riesz = m_pressureRiesz.solve(residual);
    Check(m_pressureRiesz.info() == Eigen::Success && riesz.allFinite(),
          "raw continuity Riesz solve failed");
    return std::sqrt(std::max(0., residual.dot(riesz)) / m_volume);
}

bool MixedFlow::NeedsFactor(double change) const
{
    return (!m_factor && !m_pressureCorrection) || change > .25 ||
           m_lastIterations > 30 ||
           (m_initialIterations > 0 && m_lastIterations > 2 * m_initialIterations);
}

double MixedFlow::Prepare(const Sparse &a, bool amdOrdering, MixedFlowWork *work,
                          const Vec &positiveDiagonal)
{
    const double began = Clock();
    WorkTimer setup(work ? &work->setup : nullptr);
    if (m_pressureCorrectors > 0)
    {
        ResetPreconditioner();
        PreparePressureCorrection(a, positiveDiagonal, work);
        return Clock() - began;
    }
    Sparse factor=a;
    if(m_gauge){int pivot=PressureOffset()+m_pressurePivot,last=Size()-1;std::vector<T>t;for(int c=0;c<a.outerSize();++c)for(Sparse::InnerIterator it(a,c);it;++it)if(it.row()!=pivot&&it.col()!=pivot&&it.row()!=last&&it.col()!=last)t.emplace_back(it.row(),it.col(),it.value());t.emplace_back(pivot,pivot,1);t.emplace_back(last,last,1);factor.setZero();factor.setFromTriplets(t.begin(),t.end());}
    ResetPreconditioner();
    if(amdOrdering)
    {
        m_amdFactor=std::make_unique<Eigen::SparseLU<Sparse,Eigen::AMDOrdering<int>>>();
        m_amdFactor->compute(factor);
        RequireLinear(m_amdFactor->info()==Eigen::Success,"AMD mixed component preconditioner failed");
    }
    else
    {
        m_factor=std::make_unique<Eigen::SparseLU<Sparse>>();
        m_factor->compute(factor);
        RequireLinear(m_factor->info()==Eigen::Success,"mixed component preconditioner failed");
    }
    return Clock()-began;
}

void MixedFlow::PreparePressureCorrection(const Sparse &a,
                                           const Vec &positiveDiagonal,
                                           MixedFlowWork *work)
{
    const int nv = PressureOffset(), np = PressureSize();
    Check(a.rows() == Size() && a.cols() == Size() &&
          (positiveDiagonal.size() == 0 || positiveDiagonal.size() == nv),
          "pressure correction block inventory differs");
    auto data = std::make_unique<PressureCorrection>();
    const Sparse momentum = a.topLeftCorner(nv, nv);
    {
        WorkTimer timer(work ? &work->velocityFactor : nullptr);
        data->velocityFactor.compute(momentum);
        RequireLinear(data->velocityFactor.info() == Eigen::Success,
                      "pressure correction velocity factor failed");
    }

    Sparse pressure;
    {
        WorkTimer timer(work ? &work->pressureAssembly : nullptr);
        const Vec diagonal = positiveDiagonal.size() ? positiveDiagonal : Vec(momentum.diagonal());
        RequireLinear(diagonal.allFinite() && (diagonal.array() > 0).all(),
                      "pressure correction needs a positive momentum diagonal");
        data->correctionGradient = a.block(0, nv, nv, np);
        for (int col = 0; col < np; ++col)
        {
            for (Sparse::InnerIterator entry(data->correctionGradient, col); entry; ++entry)
            {
                entry.valueRef() /= diagonal[entry.row()];
            }
        }
        data->divergence = a.block(nv, 0, np, nv);
        pressure = PressureSystem(data->divergence * data->correctionGradient);
        if (m_pressureCorrectors > 1)
        {
            data->system = a;
        }
    }
    {
        WorkTimer timer(work ? &work->pressureFactor : nullptr);
        data->pressureFactor.compute(pressure);
        RequireLinear(data->pressureFactor.info() == Eigen::Success,
                      "pressure correction Schur factor failed");
    }
    m_pressureCorrection = std::move(data);
}

Sparse MixedFlow::PressureSystem(Sparse pressure) const
{
    if (m_gauge)
    {
        // Eliminating velocity gives S deltaPi - mean deltaLambda = D predictor - rhs_c.
        const int np = PressureSize();
        std::vector<T> entries;
        entries.reserve(pressure.nonZeros() + 2 * np);
        for (int col = 0; col < np; ++col)
        {
            for (Sparse::InnerIterator entry(pressure, col); entry; ++entry)
            {
                entries.emplace_back(entry.row(), col, entry.value());
            }
            entries.emplace_back(col, np, -PressureMean()[col]);
            entries.emplace_back(np, col, PressureMean()[col]);
        }
        pressure.resize(np + 1, np + 1);
        pressure.setFromTriplets(entries.begin(), entries.end());
    }
    return pressure;
}

Vec MixedFlow::PressureUpdate(const Vec &rhs, MixedFlowWork *work) const
{
    Check(bool(m_pressureCorrection) && rhs.size() == Size(),
          "pressure correction requires prepared flow blocks");
    const auto &data = *m_pressureCorrection;
    const int nv = PressureOffset(), np = PressureSize();
    Vec predictor;
    {
        WorkTimer timer(work ? &work->predictor : nullptr);
        predictor = data.velocityFactor.solve(rhs.head(nv));
        RequireLinear(data.velocityFactor.info() == Eigen::Success && predictor.allFinite(),
                      "pressure correction momentum predictor failed");
    }
    Vec pressureIncrement;
    {
        WorkTimer timer(work ? &work->pressureSolve : nullptr);
        Vec pressureRhs(Size() - nv);
        pressureRhs.head(np) = data.divergence * predictor - rhs.segment(nv, np);
        if (m_gauge)
        {
            pressureRhs[np] = rhs[Size() - 1];
        }
        pressureIncrement = data.pressureFactor.solve(pressureRhs);
        RequireLinear(data.pressureFactor.info() == Eigen::Success &&
                      pressureIncrement.allFinite(), "pressure correction solve failed");
    }
    WorkTimer timer(work ? &work->correction : nullptr);
    Vec update(Size());
    update.head(nv) = predictor - data.correctionGradient * pressureIncrement.head(np);
    update.tail(Size() - nv) = pressureIncrement;
    RequireLinear(update.allFinite(), "pressure correction produced a nonfinite update");
    if (work)
    {
        ++work->pressureCorrections;
    }
    return update;
}

void MixedFlow::ApplyPressureCorrections(const Sparse &system, const Vec &rhs, Vec &x,
                                          MixedFlowWork *work) const
{
    for (int count = 0; count < m_pressureCorrectors; ++count)
    {
        x += PressureUpdate(rhs - system * x, work);
    }
}

Vec MixedFlow::Precondition(const Vec &rhs, MixedFlowWork *work) const
{
    if (m_pressureCorrection)
    {
        if (m_pressureCorrectors == 1)
        {
            return PressureUpdate(rhs, work);
        }
        Vec result = Vec::Zero(Size());
        ApplyPressureCorrections(m_pressureCorrection->system, rhs, result, work);
        return result;
    }
    Vec b=rhs;double mean=0,lambda=0;
    if(m_gauge){mean=b[Size()-1];lambda=PressureOne().dot(b.segment(PressureOffset(),PressureSize()))/m_volume;b.segment(PressureOffset(),PressureSize())-=PressureMean()*lambda;b[PressureOffset()+m_pressurePivot]=0;b[Size()-1]=0;}
    Vec x;
    if(m_amdFactor)
    {
        x=m_amdFactor->solve(b);
        RequireLinear(m_amdFactor->info()==Eigen::Success&&x.allFinite(),
                      "AMD mixed component preconditioner application failed");
    }
    else
    {
        RequireLinear(bool(m_factor),"mixed component preconditioner is not prepared");
        x=m_factor->solve(b);
        RequireLinear(m_factor->info()==Eigen::Success&&x.allFinite(),
                      "mixed component preconditioner application failed");
    }
    if(m_gauge){x.segment(PressureOffset(),PressureSize())-=PressureOne()*((PressureMean().dot(x.segment(PressureOffset(),PressureSize()))-mean)/m_volume);x[Size()-1]=lambda;}return x;
}
SolveInfo MixedFlow::Solve(const Sparse &a, const Vec &b, Vec &x, double change,
                           double tolerance, MixedFlowWork *work)
{
    Check(x.size()==Size()&&b.size()==Size(),"mixed component system size");SolveInfo info;double began=Clock();
    if (NeedsFactor(change))
    {
        info.factorSeconds = Prepare(a, false, work);
        ++info.rebuilds;
    }
    int used=0;auto actual=Residual(a,b,x);
    while(std::max({actual.momentum,actual.continuity,actual.gauge})>tolerance&&used<400)
    {
        const Vec delta = KrylovCorrection(
            (b - a * x).cwiseProduct(m_scale), std::min(30, 400 - used),
            tolerance * std::sqrt(m_volume),
            [&](const Vec &z) -> Vec { return (a * z).cwiseProduct(m_scale); },
            [&](const Vec &v) { return Precondition(v.cwiseQuotient(m_scale), work); },
            [&] { ++used; if (work) ++work->krylovIterations; });
        x += delta;
        info.correction += delta.segment(PressureOffset(), PressureSize()).norm();
        actual = Residual(a, b, x);
    }
    double check=Clock();auto final=Residual(a,b,x);info.residualSeconds=Clock()-check;info.momentum=final.momentum;info.continuity=final.continuity;info.gauge=final.gauge;info.iterations=used;info.solveSeconds=Clock()-began-info.factorSeconds-info.residualSeconds;
    if(std::max({info.momentum,info.continuity,info.gauge})>tolerance){std::ostringstream error;error<<std::setprecision(17)<<"mixed component CURRENT true residual failed: momentum="<<info.momentum<<", continuity="<<info.continuity<<", gauge="<<info.gauge<<", tolerance="<<tolerance<<", iterations="<<used;throw NativeLinearFailure(error.str());}
    if(!m_initialIterations&&used)m_initialIterations=used;m_lastIterations=used;return info;
}

SolveInfo MixedFlow::Correct(const Sparse &system, const Vec &rhs, Vec &x, MixedFlowWork &work,
                             double continuityTolerance, bool symmetricMomentum)
{
    Check(m_pressureCorrectors > 0 && x.size() == Size() && rhs.size() == Size() &&
          system.rows() == Size() && system.cols() == Size() &&
          std::isfinite(continuityTolerance) && continuityTolerance > 0,
          "invalid fixed pressure correction system or tolerance");
    const int nv = PressureOffset(), np = PressureSize(), nq = Size() - nv;
    // Both inner targets leave room for the final direct continuity check.
    const double pressureTolerance = (symmetricMomentum ? .1 : .75) * continuityTolerance;
    const Vec scale = m_scale.tail(nq);
    SolveInfo info;
    const double setupStart = Clock();
    Sparse momentum, gradient, divergence, massGradient, scaledMomentum;
    Vec inverseMassDiagonal, velocityScale;
    Eigen::SimplicialLDLT<Sparse> velocityFactor;
    Eigen::BiCGSTAB<Sparse, Eigen::IncompleteLUT<double>> velocityIteration;
    {
        WorkTimer setup(&work.setup);
        momentum = system.topLeftCorner(nv, nv);
        {
            WorkTimer timer(&work.velocityFactor);
            if (symmetricMomentum)
            {
                velocityFactor.compute(momentum);
                RequireLinear(velocityFactor.info() == Eigen::Success &&
                              (velocityFactor.vectorD().array() > 0).all(),
                              "fixed Picard momentum factor is not positive definite");
            }
            else
            {
                // Diagonal equilibration and incomplete LU bound the setup
                // cost even with the sampled wall term's nonlocal couplings.
                velocityScale = momentum.diagonal().cwiseAbs().array().rsqrt().matrix();
                RequireLinear(velocityScale.allFinite() && (velocityScale.array() > 0).all(),
                              "invalid diagonal scaling for implicit Picard momentum");
                scaledMomentum = velocityScale.asDiagonal() * momentum * velocityScale.asDiagonal();
                scaledMomentum.prune([](int, int, double value) { return value != 0.; });
                velocityIteration.preconditioner().setDroptol(1e-4);
                velocityIteration.preconditioner().setFillfactor(4);
                velocityIteration.setTolerance(2e-13);
                velocityIteration.setMaxIterations(150);
                velocityIteration.compute(scaledMomentum);
                RequireLinear(velocityIteration.info() == Eigen::Success,
                              "implicit Picard momentum ILUT failed");
            }
        }
        Sparse pressure;
        {
            WorkTimer timer(&work.pressureAssembly);
            gradient = system.block(0, nv, nv, np);
            divergence = system.block(nv, 0, np, nv);
            inverseMassDiagonal = m_scale.head(nv).array().square().matrix();
            massGradient = gradient;
            for (int col = 0; col < np; ++col)
            {
                for (Sparse::InnerIterator entry(massGradient, col); entry; ++entry)
                {
                    entry.valueRef() *= inverseMassDiagonal[entry.row()];
                }
            }
            pressure = PressureSystem(divergence * massGradient);
            pressure.makeCompressed();
        }
        if (!m_fixedPressureFactor || !SameCompressedMatrix(pressure, m_fixedPressureMatrix))
        {
            WorkTimer timer(&work.pressureFactor);
            auto factor = std::make_unique<Eigen::SparseLU<Sparse>>();
            factor->compute(pressure);
            RequireLinear(factor->info() == Eigen::Success,
                          "fixed Picard pressure preconditioner factor failed");
            m_fixedPressureMatrix = std::move(pressure);
            m_fixedPressureFactor = std::move(factor);
        }
    }
    info.factorSeconds = Clock() - setupStart;
    ++info.rebuilds;

    const auto velocitySolve = [&](const Vec &load) -> Vec
    {
        if (symmetricMomentum)
        {
            Vec velocity = velocityFactor.solve(load);
            RequireLinear(velocityFactor.info() == Eigen::Success && velocity.allFinite(),
                          "fixed Picard momentum solve failed");
            return velocity;
        }
        const Vec scaledLoad = load.cwiseProduct(velocityScale);
        const Vec velocity = velocityIteration.solve(scaledLoad);
        const double defect = (scaledMomentum * velocity - scaledLoad).norm();
        if (!(velocity.allFinite() && defect <= 2e-12 * scaledLoad.norm()))
        {
            std::ostringstream error;
            error << "implicit Picard momentum did not converge in "
                  << velocityIteration.iterations() << " linear iterations";
            throw NativeLinearFailure(error.str());
        }
        return velocity.cwiseProduct(velocityScale);
    };
    const auto schur = [&](const Vec &pressure) -> Vec
    {
        Vec action(nq);
        action.head(np) = divergence * velocitySolve(gradient * pressure.head(np));
        if (m_gauge)
        {
            action.head(np) -= PressureMean() * pressure[np];
            action[np] = PressureMean().dot(pressure.head(np));
        }
        return action;
    };
    const auto precondition = [&](const Vec &scaledRhs) -> Vec
    {
        // Least-squares commutator using L = D Mdiag^-1 G and the frozen A.
        // It accelerates the pressure-only iteration. The actual action uses A^-1.
        const Vec load = scaledRhs.cwiseQuotient(scale);
        const auto &pressureFactor = *m_fixedPressureFactor;
        const Vec first = pressureFactor.solve(load);
        RequireLinear(pressureFactor.info() == Eigen::Success && first.allFinite(),
                      "fixed Picard pressure preconditioner solve failed");
        const Vec weighted = (momentum * (massGradient * first.head(np))).cwiseProduct(inverseMassDiagonal);
        Vec middle(nq);
        middle.head(np) = divergence * weighted;
        if (m_gauge)
        {
            middle.head(np) -= PressureMean() * first[np];
            middle[np] = load[np];
        }
        Vec result = pressureFactor.solve(middle);
        RequireLinear(pressureFactor.info() == Eigen::Success && result.allFinite(),
                      "fixed Picard pressure preconditioner solve failed");
        return result;
    };

    const double began = Clock();
    for (int pass = 0; pass < m_pressureCorrectors; ++pass)
    {
        const Vec residual = rhs - system * x;
        Vec predictor;
        {
            WorkTimer timer(&work.predictor);
            predictor = velocitySolve(residual.head(nv));
        }
        Vec pressure = Vec::Zero(nq);
        {
            WorkTimer timer(&work.pressureSolve);
            Vec pressureRhs(nq);
            pressureRhs.head(np) = divergence * predictor - residual.segment(nv, np);
            if (m_gauge)
            {
                pressureRhs[np] = residual[Size() - 1];
            }
            Vec defect = pressureRhs;
            const auto norm = [&]()
            {
                const double divergenceNorm = ContinuityNorm(defect.head(np));
                return m_gauge ? std::max(divergenceNorm, std::abs(defect[np] * scale[np]))
                               : divergenceNorm;
            };
            const auto scaledAction = [&](const Vec &p) -> Vec { return schur(p).cwiseProduct(scale); };
            const int maximumIterations = symmetricMomentum ? 240 : 400;
            const int restart = symmetricMomentum ? 60 : 100;
            int iterations = 0;
            double current = norm();
            while (current > pressureTolerance && iterations < maximumIterations)
            {
                const Vec scaledDefect = defect.cwiseProduct(scale);
                // Relate the diagonal-scaled Krylov norm to the actual
                // continuity norm at each restart of the nonsymmetric solve.
                const double projectedTolerance = symmetricMomentum
                    ? pressureTolerance * std::sqrt(m_volume)
                    : .8 * pressureTolerance * (scaledDefect.norm() / current);
                pressure += KrylovCorrection(scaledDefect,
                    std::min(restart, maximumIterations - iterations),
                    projectedTolerance, scaledAction, precondition,
                    [&] { ++iterations; ++work.pressureKrylovIterations; });
                defect = pressureRhs - schur(pressure);
                current = norm();
            }
            if (current > pressureTolerance)
            {
                std::ostringstream error;
                error << std::setprecision(17) << "fixed Picard pressure residual failed: residual="
                      << current << ", tolerance=" << pressureTolerance << ", iterations=" << iterations;
                throw NativeLinearFailure(error.str());
            }
        }
        WorkTimer timer(&work.correction);
        x.head(nv) += predictor - velocitySolve(gradient * pressure.head(np));
        x.tail(nq) += pressure;
        info.correction += pressure.head(np).norm();
        ++work.pressureCorrections;
    }
    const double solveEnd = Clock();
    const auto final = Residual(system, rhs, x);
    info.momentum = final.momentum;
    info.continuity = final.continuity;
    info.gauge = final.gauge;
    work.fixedMomentum = std::max(work.fixedMomentum, final.momentum);
    work.fixedContinuity = std::max(work.fixedContinuity, final.continuity);
    work.fixedGauge = std::max(work.fixedGauge, final.gauge);
    if (!symmetricMomentum)
    {
        RequireLinear(std::max({final.momentum, final.continuity, final.gauge}) <=
                      continuityTolerance,
                      "implicit Picard correction failed the direct frozen-system residual check");
    }
    info.solveSeconds = solveEnd - began;
    info.residualSeconds = Clock() - solveEnd;
    return info;
}

void MixedFlow::SetTransportedKSpace(const MultiRegions::ContFieldSharedPtr&field)
{
    Check(bool(field)&&field->GetComm()->GetSize()==1&&field->GetExpSize()==int(Elements().size()),"native Q3 k must be serial over the full domain");m_kField=field;
    if (m_core.Mapping().normalMap)
    {
        m_scalarAction = m_core.NormalSampledAction(field);
        return;
    }
    auto map=field->GetLocalToGlobalMap();const int nk=map->GetNumGlobalCoeffs();std::vector<T>triplets;auto rule=Gauss(16);
    struct Trace{LocalRegions::ExpansionSharedPtr exp;Mat dx,dz;std::vector<int>ids;};std::vector<Trace>traces;
    for(int e=0;e<int(Elements().size());++e){auto exp=field->GetExp(e);Check(exp->GetGeom()->GetGlobalID()==Elements()[e].geometry&&exp->GetNcoeffs()==16,"transported k geometry/Q3 basis differs");for(int d=0;d<2;++d)Check(exp->GetBasisNumModes(d)==4,"transported k is not tensor Q3");Trace t;t.exp=exp;t.dx=NativeBasis(exp,0);t.dz=NativeBasis(exp,1);for(int j=0;j<16;++j){int local=field->GetCoeff_Offset(e)+j;t.ids.push_back(map->GetLocalToGlobalMap(local));double sign=map->GetLocalToGlobalSign(local);t.dx.col(j)*=sign;t.dz.col(j)*=sign;}traces.push_back(std::move(t));}
    const auto geo=m_core.Mapping();
    for(const auto&e:Elements())
    {
        double low=std::max(0.,e.rl),high=std::min(e.rr,geo.h);if(high<=low)continue;int owner=-1;
        for(int j=0;j<int(Elements().size());++j){const auto&t=Elements()[j];if(std::abs(t.xl-e.xl)<1e-10&&std::abs(t.xr-e.xr)<1e-10&&t.rl<=geo.H+1e-12&&t.rr>geo.H+1e-12){owner=j;break;}}
        Check(owner>=0,"missing outer-owned Q3 k trace at H");const auto&native=traces[owner];const auto&element=Elements()[owner];Mat lx=Mat::Zero(16,16),lz=lx;
        for(int i=0;i<16;++i){double xi=rule.first[i],x=(e.xl+e.xr)/2+(e.xr-e.xl)*xi/2,slope=geo.Slope(x),etaH=2*(geo.H-element.rl)/(element.rr-element.rl)-1;Vec tangent=(At(native.exp,native.dx,xi,etaH)+slope*At(native.exp,native.dz,xi,etaH))/(1+slope*slope);
            for(int j=0;j<16;++j){double r=(low+high)/2+(high-low)*rule.first[j]/2,eta=2*(r-e.rl)/(e.rr-e.rl)-1;Vec b=At(e.velocity,e.B,xi,eta);double weight=rule.second[i]*rule.second[j]*(e.xr-e.xl)*(high-low)/4;lx.noalias()+=weight*b*tangent.transpose();lz.noalias()+=weight*slope*b*tangent.transpose();}}
        for(int i=0;i<16;++i)for(int j=0;j<16;++j){triplets.emplace_back(e.u[i],native.ids[j],lx(i,j));triplets.emplace_back(VelocitySize()+e.u[i],native.ids[j],lz(i,j));}
    }
    m_scalarAction.resize(2*VelocitySize(),nk);m_scalarAction.setFromTriplets(triplets.begin(),triplets.end());
}
Vec MixedFlow::ProjectTransportedK(const Function&f)const
{Check(bool(m_kField),"Q3 k space has not been registered");return ProjectNative(f,m_kField);}
Vec MixedFlow::KnownScalarKLoad(const Vec&k)const
{Check(bool(m_kField)&&k.size()==m_scalarAction.cols(),"actual transported Q3 k coefficient size differs");return(-2./3)*(m_scalarAction*k);}
}
