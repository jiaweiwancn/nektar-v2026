// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeMixedFlow.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeFailures.h>
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

MixedFlow::MixedFlow(const std::vector<MultiRegions::ContFieldSharedPtr>&fields,Geometry geometry)
    :m_fields(fields),m_core({fields.at(0),fields.at(0),fields.at(2)},geometry)
{
    Check(fields.size()==3,"H27 MixedFlow requires serial u,v,p fields");
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
    Vec v=residual.segment(PressureOffset(),PressureSize()),z=m_pressureRiesz.solve(v);
    Check(m_pressureRiesz.info()==Eigen::Success&&z.allFinite(),"raw continuity Riesz solve failed");
    info.continuity=std::sqrt(std::max(0.,v.dot(z))/m_volume);
    info.gauge=m_gauge?std::abs(residual[Size()-1])/(.01*m_volume):0;return info;
}
double MixedFlow::Prepare(const Sparse&a,bool amdOrdering)
{
    double began=Clock();Sparse factor=a;
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
Vec MixedFlow::Precondition(const Vec&rhs)const
{
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
SolveInfo MixedFlow::Solve(const Sparse&a,const Vec&b,Vec&x,double change,double tolerance)
{
    Check(x.size()==Size()&&b.size()==Size(),"mixed component system size");SolveInfo info;double began=Clock();
    if(!m_factor||change>.25||m_lastIterations>30||(m_initialIterations>0&&m_lastIterations>2*m_initialIterations)){info.factorSeconds=Prepare(a);++info.rebuilds;}
    int used=0;auto actual=Residual(a,b,x);
    while(std::max({actual.momentum,actual.continuity,actual.gauge})>tolerance&&used<400)
    {
        Vec residual=(b-a*x).cwiseProduct(m_scale);double beta=residual.norm();RequireLinear(std::isfinite(beta)&&beta>0,"invalid mixed Krylov residual");
        int dimension=std::min(30,400-used);Mat v=Mat::Zero(Size(),dimension+1),z(Size(),dimension),h=Mat::Zero(dimension+1,dimension);v.col(0)=residual/beta;Vec rhs=Vec::Zero(dimension+1);rhs[0]=beta;Vec best;int done=0;
        for(int j=0;j<dimension;++j){z.col(j)=Precondition(v.col(j).cwiseQuotient(m_scale));Vec action=(a*z.col(j)).cwiseProduct(m_scale);for(int pass=0;pass<2;++pass)for(int i=0;i<=j;++i){double value=v.col(i).dot(action);h(i,j)+=value;action.noalias()-=value*v.col(i);}h(j+1,j)=action.norm();if(h(j+1,j)>0)v.col(j+1)=action/h(j+1,j);best=h.block(0,0,j+2,j+1).colPivHouseholderQr().solve(rhs.head(j+2));++used;done=j+1;if((rhs.head(j+2)-h.block(0,0,j+2,j+1)*best).norm()<tolerance*std::sqrt(m_volume)||h(j+1,j)<1e-28)break;}
        Vec delta=z.leftCols(done)*best;x+=delta;info.correction+=delta.segment(PressureOffset(),PressureSize()).norm();actual=Residual(a,b,x);
    }
    double check=Clock();auto final=Residual(a,b,x);info.residualSeconds=Clock()-check;info.momentum=final.momentum;info.continuity=final.continuity;info.gauge=final.gauge;info.iterations=used;info.solveSeconds=Clock()-began-info.factorSeconds-info.residualSeconds;
    if(std::max({info.momentum,info.continuity,info.gauge})>tolerance){std::ostringstream error;error<<std::setprecision(17)<<"mixed component CURRENT true residual failed: momentum="<<info.momentum<<", continuity="<<info.continuity<<", gauge="<<info.gauge<<", tolerance="<<tolerance<<", iterations="<<used;throw NativeLinearFailure(error.str());}
    if(!m_initialIterations&&used)m_initialIterations=used;m_lastIterations=used;return info;
}

void MixedFlow::SetTransportedKSpace(const MultiRegions::ContFieldSharedPtr&field)
{
    Check(bool(field)&&field->GetComm()->GetSize()==1&&field->GetExpSize()==int(Elements().size()),"native Q3 k must be serial over the full domain");m_kField=field;auto map=field->GetLocalToGlobalMap();const int nk=map->GetNumGlobalCoeffs();std::vector<T>triplets;auto rule=Gauss(16);
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
