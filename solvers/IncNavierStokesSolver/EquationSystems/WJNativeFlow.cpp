// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeFlow.h>
#include <Eigen/QR>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace Nektar::WJNative
{
namespace
{
using Triplet=Eigen::Triplet<double>;
void Require(bool v,const char*m){if(!v)throw std::runtime_error(m);}
double Now(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
std::pair<Vec,Vec>Gauss(int n)
{
    Vec x(n),w(n);
    for(int i=0;i<n;++i)
    {
        double z=std::cos(M_PI*(i+.75)/(n+.5)),d=0;
        for(int k=0;k<30;++k){double p0=1,p1=z;for(int j=2;j<=n;++j){double p=((2*j-1)*z*p1-(j-1)*p0)/j;p0=p1;p1=p;}d=n*(z*p1-p0)/(z*z-1);double dz=p1/d;z-=dz;if(std::abs(dz)<2e-16)break;}
        x[i]=z;w[i]=2/((1-z*z)*d*d);
    }
    return {x,w};
}
Mat Basis(const LocalRegions::ExpansionSharedPtr&e,int derivative)
{
    int nq=e->GetTotPoints(),nc=e->GetNcoeffs();Mat b(nq,nc);
    Array<OneD,NekDouble>c(nc,0.),p(nq),dx(nq),dz(nq);
    for(int j=0;j<nc;++j){Vmath::Zero(nc,c,1);c[j]=1;e->BwdTrans(c,p);if(derivative>=0)e->PhysDeriv(p,dx,dz);for(int q=0;q<nq;++q)b(q,j)=derivative<0?p[q]:derivative==0?dx[q]:dz[q];}
    return b;
}
Vec Evaluate(const LocalRegions::ExpansionSharedPtr&e,const Mat&basis,double xi,double eta)
{
    Array<OneD,NekDouble>loc(2);loc[0]=xi;loc[1]=eta;Array<OneD,NekDouble>v(e->GetTotPoints());Vec b(basis.cols());
    for(int j=0;j<basis.cols();++j){for(int q=0;q<basis.rows();++q)v[q]=basis(q,j);b[j]=e->StdPhysEvaluate(loc,v);}return b;
}
}

WallStress Adapt(double k,double epsilon,const Tensor&G,const Vector&n,double nu,double nuA,double r,double H)
{
    Require(k>0&&epsilon>0&&nu>0&&nuA>=0&&H>0,"invalid prescribed wall stress inputs");
    double nn=0;for(double v:n)nn+=v*v;Require(std::abs(nn-1)<1e-12,"wall normal is not unit");
    const auto w=SolverUtils::EvaluateWJClosure(k,epsilon,G,{},SolverUtils::WJGradientPolicy::DevStrainRawWorkV1,0);
    Tensor P{},M{},Rh{};Vector Mn{},b{};double normal=0;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){P[i][j]=(i==j?1.:0.)-n[i]*n[j];M[i][j]=G[i][j]+G[j][i];Mn[i]+=M[i][j]*n[j];normal+=n[i]*w.reynoldsStress[i][j]*n[j];}
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)b[i]-=nuA*P[i][j]*Mn[j];
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){Rh[i][j]=normal*n[i]*n[j]+b[i]*n[j]+n[i]*b[j];for(int a=0;a<3;++a)for(int c=0;c<3;++c)Rh[i][j]+=P[i][a]*w.reynoldsStress[a][c]*P[c][j];}
    const double s=std::clamp((r/H-.4)/.35,0.,1.);WallStress out;out.chi=s*s*s*(10+s*(-15+6*s));out.mu=nu+(1-out.chi)*nuA+out.chi*w.eddyViscosity;
    Require(out.mu>0,"nonpositive implicit diffusion");
    for(int i=0;i<3;++i)for(int j=0;j<3;++j){out.R[i][j]=(1-out.chi)*Rh[i][j]+out.chi*w.reynoldsStress[i][j];out.V[i][j]=nu*M[i][j]-out.R[i][j]+(i==j?2*k/3:0);out.remainder[i][j]=out.V[i][j]-out.mu*M[i][j];out.production-=out.R[i][j]*G[i][j];}
    return out;
}

Flow::Flow(const std::vector<MultiRegions::ContFieldSharedPtr>&fields,Geometry geometry):m_fields(fields),m_geometry(geometry)
{
    Require(fields.size()==3&&fields[0]->GetComm()->GetSize()==1,"H26 starts with serial u,v,p ContFields");
    auto um=fields[0]->GetLocalToGlobalMap(),pm=fields[2]->GetLocalToGlobalMap(),vm=fields[1]->GetLocalToGlobalMap();
    m_nu=um->GetNumGlobalCoeffs();m_np=pm->GetNumGlobalCoeffs();m_nd=um->GetNumGlobalDirBndCoeffs();m_nf=m_nu-m_nd;
    Require(vm->GetNumGlobalCoeffs()==m_nu&&vm->GetNumGlobalDirBndCoeffs()==m_nd,"velocity homogeneous increment spaces differ");
    for(auto bc:fields[0]->GetBndConditions())if(bc->GetBoundaryConditionType()==SpatialDomains::eNeumann||bc->GetBoundaryConditionType()==SpatialDomains::eRobin)m_gauge=false;
    for(int e=0;e<fields[0]->GetExpSize();++e)
    {
        Element a;a.velocity=fields[0]->GetExp(e);a.pressure=fields[2]->GetExp(e);a.geometry=a.velocity->GetGeom()->GetGlobalID();
        Require(a.velocity->GetNcoeffs()==16&&a.pressure->GetNcoeffs()==9,"H26 requires actual tensor Q3/Q2 spaces");
        int nq=a.velocity->GetTotPoints();Require(nq==a.pressure->GetTotPoints(),"Q3/Q2 quadrature grids differ");
        Array<OneD,NekDouble>x(nq),z(nq),dummy(nq),one(nq,1.),weight(nq);a.velocity->GetCoords(x,z,dummy);a.velocity->MultiplyByQuadratureMetric(one,weight);
        a.x=Vec(nq);a.z=Vec(nq);a.w=Vec(nq);a.xl=a.rl=1e100;a.xr=a.rr=-1e100;
        for(int q=0;q<nq;++q){a.x[q]=x[q];a.z[q]=z[q];a.w[q]=weight[q];a.xl=std::min(a.xl,x[q]);a.xr=std::max(a.xr,x[q]);double r=z[q]-geometry.Ground(x[q]);a.rl=std::min(a.rl,r);a.rr=std::max(a.rr,r);m_volume+=weight[q];}
        a.B=Basis(a.velocity,-1);a.Dx=Basis(a.velocity,0);a.Dz=Basis(a.velocity,1);a.P=Basis(a.pressure,-1);a.Px=Basis(a.pressure,0);a.Pz=Basis(a.pressure,1);
        for(int j=0;j<16;++j){int l=fields[0]->GetCoeff_Offset(e)+j;a.u.push_back(um->GetLocalToGlobalMap(l));double s=um->GetLocalToGlobalSign(l);a.us.push_back(s);a.B.col(j)*=s;a.Dx.col(j)*=s;a.Dz.col(j)*=s;Require(a.u.back()==vm->GetLocalToGlobalMap(l)&&s==vm->GetLocalToGlobalSign(l),"velocity signed map differs");}
        for(int j=0;j<9;++j){int l=fields[2]->GetCoeff_Offset(e)+j;a.p.push_back(pm->GetLocalToGlobalMap(l));double s=pm->GetLocalToGlobalSign(l);a.ps.push_back(s);a.P.col(j)*=s;a.Px.col(j)*=s;a.Pz.col(j)*=s;}
        m_elements.push_back(std::move(a));
    }
    Require(m_volume>0,"nonpositive physical volume");m_length=std::sqrt(m_volume);m_forceScale=1+1/m_length+(6.423173803526447e-5+.3819178635467331*.051*m_elements.back().rr/2)/(m_length*m_length)+.01/m_length;
    BuildInvariant();BuildAWF();m_scale=RowScale();
}

void Flow::BuildInvariant()
{
    std::vector<Triplet>m,q,d;m_mean=Vec::Zero(m_np);
    for(const auto&e:m_elements)
    {
        Mat lm=e.B.transpose()*e.w.asDiagonal()*e.B,lq=e.P.transpose()*e.w.asDiagonal()*e.P;
        Mat dx=e.P.transpose()*e.w.asDiagonal()*e.Dx,dz=e.P.transpose()*e.w.asDiagonal()*e.Dz;
        Vec mean=e.P.transpose()*e.w;
        for(int i=0;i<16;++i)for(int j=0;j<16;++j)m.emplace_back(e.u[i],e.u[j],lm(i,j));
        for(int i=0;i<9;++i){m_mean[e.p[i]]+=mean[i];for(int j=0;j<9;++j)q.emplace_back(e.p[i],e.p[j],lq(i,j));for(int j=0;j<16;++j){d.emplace_back(e.p[i],e.u[j],dx(i,j));d.emplace_back(e.p[i],m_nu+e.u[j],dz(i,j));}}
    }
    m_mass.resize(m_nu,m_nu);m_mass.setFromTriplets(m.begin(),m.end());m_qmass.resize(m_np,m_np);m_qmass.setFromTriplets(q.begin(),q.end());m_div.resize(m_np,2*m_nu);m_div.setFromTriplets(d.begin(),d.end());m_grad=-Sparse(m_div.transpose());
    m_freeMass=m_mass.block(m_nd,m_nd,m_nf,m_nf);m_velocityRiesz.compute(m_freeMass);m_pressureRiesz.compute(m_qmass);
    Require(m_velocityRiesz.info()==Eigen::Success&&m_pressureRiesz.info()==Eigen::Success,"native sparse mass/Riesz factor failed");
}

void Flow::BuildAWF()
{
    std::vector<Triplet>t;auto quadrature=Gauss(16);
    for(const auto&e:m_elements)
    {
        double upper=std::min(e.rr,m_geometry.h),lower=std::max(e.rl,0.);if(upper<=lower)continue;
        const Element*trace=nullptr;
        for(const auto&p:m_elements)if(std::abs(p.xl-e.xl)<1e-10&&std::abs(p.xr-e.xr)<1e-10&&p.rl<=m_geometry.H+1e-12&&p.rr>m_geometry.H+1e-12){trace=&p;break;}
        Require(trace,"missing native pressure element at H");
        Mat localX=Mat::Zero(16,9),localZ=localX;
        for(int a=0;a<16;++a)
        {
            double xi=quadrature.first[a],x=(e.xl+e.xr)/2+(e.xr-e.xl)*xi/2,hp=m_geometry.Slope(x),etaH=2*(m_geometry.H-trace->rl)/(trace->rr-trace->rl)-1;
            Vec px=Evaluate(trace->pressure,trace->Px,xi,etaH),pz=Evaluate(trace->pressure,trace->Pz,xi,etaH);
            // P_H grad p = (1,h') (p_x+h'p_z)/(1+h'^2).
            Vec tangent=(px+hp*pz)/(1+hp*hp);
            for(int b=0;b<16;++b){double r=(lower+upper)/2+(upper-lower)*quadrature.first[b]/2,eta=2*(r-e.rl)/(e.rr-e.rl)-1;Vec v=Evaluate(e.velocity,e.B,xi,eta);double w=quadrature.second[a]*quadrature.second[b]*(e.xr-e.xl)*(upper-lower)/4;localX.noalias()+=w*v*tangent.transpose();localZ.noalias()+=w*hp*v*tangent.transpose();}
        }
        for(int i=0;i<16;++i)for(int j=0;j<9;++j){t.emplace_back(e.u[i],trace->p[j],localX(i,j));t.emplace_back(m_nu+e.u[i],trace->p[j],localZ(i,j));}
    }
    m_awf.resize(2*m_nu,m_np);m_awf.setFromTriplets(t.begin(),t.end());
}

Sparse Flow::Strain(const Function&mu,bool transpose)const
{
    std::vector<Triplet>t;t.reserve(m_elements.size()*1024);
    for(const auto&e:m_elements)
    {
        Vec w=e.w;for(int q=0;q<w.size();++q){double v=mu(e.x[q],e.z[q]);Require(std::isfinite(v)&&v>=0,"invalid physical implicit viscosity");w[q]*=v;}
        std::array<Mat,2>grad{e.Dx,e.Dz};Mat lap=e.Dx.transpose()*w.asDiagonal()*e.Dx+e.Dz.transpose()*w.asDiagonal()*e.Dz;
        for(int i=0;i<2;++i)for(int j=0;j<2;++j){Mat k=Mat::Zero(16,16);if(i==j)k=lap;if(transpose)k.noalias()+=grad[j].transpose()*w.asDiagonal()*grad[i];for(int a=0;a<16;++a)for(int b=0;b<16;++b)t.emplace_back(i*m_nu+e.u[a],j*m_nu+e.u[b],k(a,b));}
    }
    Sparse k(2*m_nu,2*m_nu);k.setFromTriplets(t.begin(),t.end());return k;
}

Sparse Flow::System(double alpha,const Sparse&K,double awf)const
{
    std::vector<Triplet>t;auto id=[&](int v){int c=v/m_nu,g=v%m_nu;return g<m_nd?-1:c*m_nf+g-m_nd;};
    for(int k=0;k<K.outerSize();++k)for(Sparse::InnerIterator it(K,k);it;++it){int i=id(it.row()),j=id(it.col());if(i>=0&&j>=0)t.emplace_back(i,j,it.value());}
    for(int c=0;c<2;++c)for(int k=0;k<m_mass.outerSize();++k)for(Sparse::InnerIterator it(m_mass,k);it;++it)if(it.row()>=m_nd&&it.col()>=m_nd)t.emplace_back(c*m_nf+it.row()-m_nd,c*m_nf+it.col()-m_nd,alpha*it.value());
    Sparse g=m_grad-awf*m_awf;
    for(int k=0;k<g.outerSize();++k)for(Sparse::InnerIterator it(g,k);it;++it){int i=id(it.row());if(i>=0)t.emplace_back(i,2*m_nf+it.col(),it.value());}
    for(int k=0;k<m_div.outerSize();++k)for(Sparse::InnerIterator it(m_div,k);it;++it){int j=id(it.col());if(j>=0)t.emplace_back(2*m_nf+it.row(),j,it.value());}
    if(m_gauge)for(int i=0;i<m_np;++i){t.emplace_back(2*m_nf+i,Size()-1,m_mean[i]);t.emplace_back(Size()-1,2*m_nf+i,m_mean[i]);}
    Sparse a(Size(),Size());a.setFromTriplets(t.begin(),t.end());return a;
}

Vec Flow::ApplyFull(double alpha,const Sparse&K,double awf,const Vec&x)const
{
    Require(x.size()==2*m_nu+m_np,"full flow state size");Vec r=Vec::Zero(x.size());r.head(2*m_nu)=K*x.head(2*m_nu)+(m_grad-awf*m_awf)*x.tail(m_np);
    for(int c=0;c<2;++c)r.segment(c*m_nu,m_nu)+=alpha*(m_mass*x.segment(c*m_nu,m_nu));r.tail(m_np)=m_div*x.head(2*m_nu);return r;
}
Vec Flow::Load(const std::array<Function,2>&f)const
{
    Vec r=Vec::Zero(2*m_nu+m_np);for(const auto&e:m_elements)for(int c=0;c<2;++c){Vec v(e.w.size());for(int q=0;q<v.size();++q)v[q]=e.w[q]*f[c](e.x[q],e.z[q]);Vec b=e.B.transpose()*v;for(int j=0;j<16;++j)r[c*m_nu+e.u[j]]+=b[j];}return r;
}
Vec Flow::StressLoad(const std::function<Tensor(double,double)>&T)const
{
    Vec r=Vec::Zero(2*m_nu+m_np);for(const auto&e:m_elements){std::array<Vec,4>s;for(auto&v:s)v=Vec(e.w.size());for(int q=0;q<e.w.size();++q){auto a=T(e.x[q],e.z[q]);for(int i=0;i<2;++i)for(int j=0;j<2;++j)s[2*i+j][q]=e.w[q]*a[i?2:0][j?2:0];}for(int c=0;c<2;++c){Vec b=-e.Dx.transpose()*s[2*c]-e.Dz.transpose()*s[2*c+1];for(int j=0;j<16;++j)r[c*m_nu+e.u[j]]+=b[j];}}return r;
}
Vec Flow::AWFLoad(const std::function<Vector(double,double)>&sample)const
{
    Vec result=Vec::Zero(2*m_nu+m_np);auto rule=Gauss(16);
    for(const auto&e:m_elements)
    {
        double lower=std::max(0.,e.rl),upper=std::min(e.rr,m_geometry.h);if(upper<=lower)continue;
        for(int a=0;a<16;++a){double xi=rule.first[a],x=(e.xl+e.xr)/2+(e.xr-e.xl)*xi/2;auto f=sample(x,m_geometry.Ground(x)+m_geometry.H);for(int b=0;b<16;++b){double r=(lower+upper)/2+(upper-lower)*rule.first[b]/2,eta=2*(r-e.rl)/(e.rr-e.rl)-1;Vec v=Evaluate(e.velocity,e.B,xi,eta);double w=rule.second[a]*rule.second[b]*(e.xr-e.xl)*(upper-lower)/4;for(int j=0;j<16;++j){result[e.u[j]]+=w*v[j]*f[0];result[m_nu+e.u[j]]+=w*v[j]*f[2];}}}
    }return result;
}
Vec Flow::KnownKLoad(const Vec&k)const
{
    Require(k.size()==m_np,"known mechanical-pressure k trace size");Vec r=Vec::Zero(2*m_nu+m_np);r.head(2*m_nu)=(-2./3)*(m_awf*k);return r;
}
Vec Flow::MassAction(const Vec&x)const
{
    Vec y=Vec::Zero(2*m_nu+m_np);for(int c=0;c<2;++c)y.segment(c*m_nu,m_nu)=m_mass*x.segment(c*m_nu,m_nu);return y;
}
Vec Flow::Convection(const Vec&x)const
{
    Vec out=Vec::Zero(2*m_nu+m_np);
    for(const auto&e:m_elements){Vec u(16),v(16);for(int j=0;j<16;++j){u[j]=x[e.u[j]];v[j]=x[m_nu+e.u[j]];}Vec up=e.B*u,vp=e.B*v;for(int c=0;c<2;++c){const auto&coeff=c==0?u:v;Vec conv=up.cwiseProduct(e.Dx*coeff)+vp.cwiseProduct(e.Dz*coeff);Vec rhs=e.B.transpose()*e.w.asDiagonal()*conv;for(int j=0;j<16;++j)out[c*m_nu+e.u[j]]+=rhs[j];}}
    return out;
}
double Flow::PreparePreconditioner(const Sparse&a)
{
    double t=Now();m_factor=std::make_unique<Eigen::SparseLU<Sparse>>();m_factor->analyzePattern(a);m_factor->factorize(a);Require(m_factor->info()==Eigen::Success,"sparse native-block preconditioner factorization failed");m_initialIterations=m_lastIterations=0;return Now()-t;
}
Vec Flow::TractionLoad(int boundary,const std::function<Vector(double,double,const Vector&)>&fun)const
{
    Vec result=Vec::Zero(2*m_nu+m_np);auto rule=Gauss(20);
    for(const auto&e:m_elements)
    {
        bool face=boundary==0?std::abs(e.xl-m_elements.front().xl)<1e-10:boundary==1?std::abs(e.xr-m_elements.back().xr)<1e-10:boundary==2?std::abs(e.rl-m_elements.front().rl)<1e-10:std::abs(e.rr-m_elements.back().rr)<1e-10;if(!face)continue;
        for(int q=0;q<20;++q){double xi=boundary<2?(boundary?1.:-1.):rule.first[q],eta=boundary<2?rule.first[q]:(boundary==2?-1.:1.),x=(e.xl+e.xr)/2+(e.xr-e.xl)*xi/2,r=(e.rl+e.rr)/2+(e.rr-e.rl)*eta/2,hp=m_geometry.Slope(x);Vector n{};double jac=0;if(boundary<2){n={boundary?1.:-1.,0,0};jac=(e.rr-e.rl)/2;}else{double s=boundary==2?-1.:1.,norm=std::sqrt(1+hp*hp);n={-s*hp/norm,0,s/norm};jac=(e.xr-e.xl)*norm/2;}Vec b=Evaluate(e.velocity,e.B,xi,eta);auto t=fun(x,r+m_geometry.Ground(x),n);for(int i=0;i<16;++i){result[e.u[i]]+=b[i]*t[0]*jac*rule.second[q];result[m_nu+e.u[i]]+=b[i]*t[2]*jac*rule.second[q];}}
    }return result;
}
Vec Flow::Project(const Function&f,int c)const
{
    auto field=m_fields[c];Array<OneD,NekDouble>x(field->GetTotPoints()),z(x.size()),d(x.size()),v(x.size()),coef(field->GetNcoeffs());field->GetCoords(x,z,d);for(int q=0;q<x.size();++q)v[q]=f(x[q],z[q]);field->FwdTransLocalElmt(v,coef);
    for(int b=0;b<field->GetBndConditions().size();++b)if(field->GetBndConditions()[b]->GetBoundaryConditionType()==SpatialDomains::eDirichlet){auto boundary=field->GetBndCondExpansions()[b];int nq=boundary->GetTotPoints();Array<OneD,NekDouble>bx(nq),bz(nq),bd(nq),values(nq);boundary->GetCoords(bx,bz,bd);for(int q=0;q<nq;++q)values[q]=f(bx[q],bz[q]);boundary->FwdTrans(values,boundary->UpdateCoeffs());}
    field->ImposeDirichletConditions(coef);field->UpdateCoeffs()=coef;field->AvgAssemble();field->GlobalToLocal();auto map=field->GetLocalToGlobalMap();Array<OneD,NekDouble>g(map->GetNumGlobalCoeffs());map->LocalToGlobal(field->GetCoeffs(),g,true);Vec out(g.size());for(int i=0;i<g.size();++i)out[i]=g[i];return out;
}
Vec Flow::Lift(const Vec&x)const{Vec l=Vec::Zero(x.size());for(int c=0;c<2;++c)l.segment(c*m_nu,m_nd)=x.segment(c*m_nu,m_nd);return l;}
Vec Flow::Restrict(const Vec&x)const{Vec y=Vec::Zero(Size());for(int c=0;c<2;++c)y.segment(c*m_nf,m_nf)=x.segment(c*m_nu+m_nd,m_nf);y.segment(2*m_nf,m_np)=x.tail(m_np);return y;}
Vec Flow::Expand(const Vec&x,const Vec&lift)const{Vec y=lift;for(int c=0;c<2;++c)y.segment(c*m_nu+m_nd,m_nf)=x.segment(c*m_nf,m_nf);y.tail(m_np)=x.segment(2*m_nf,m_np);return y;}
Vec Flow::LiftedRhs(double alpha,const Sparse&K,double awf,const Vec&load,const Vec&lift)const{return Restrict(load-ApplyFull(alpha,K,awf,lift));}
Vec Flow::RowScale()const{Vec s=Vec::Ones(Size());for(int c=0;c<2;++c)for(int i=0;i<m_nf;++i)s[c*m_nf+i]=1/(std::sqrt(m_mass.coeff(i+m_nd,i+m_nd))*m_forceScale);for(int i=0;i<m_np;++i)s[2*m_nf+i]=m_length/std::sqrt(m_qmass.coeff(i,i));if(m_gauge)s[Size()-1]=1/(.01*m_volume);return s;}

SolveInfo Flow::Residual(const Sparse&a,const Vec&b,const Vec&x)const
{
    Vec r=a*x-b;SolveInfo info;double norm=0;for(int c=0;c<2;++c){Vec v=r.segment(c*m_nf,m_nf),z=m_velocityRiesz.solve(v);norm+=v.dot(z);}info.momentum=std::sqrt(std::max(0.,norm)/m_volume)/m_forceScale;
    Vec v=r.segment(2*m_nf,m_np);if(m_gauge)v-=m_mean*x[Size()-1];Vec z=m_pressureRiesz.solve(v);info.continuity=std::sqrt(std::max(0.,v.dot(z))/m_volume)*m_length;info.gauge=m_gauge?std::abs(m_mean.dot(x.segment(2*m_nf,m_np)))/(.01*m_volume):0;return info;
}
SolveInfo Flow::Solve(const Sparse&a,const Vec&b,Vec&x,double relativeChange,double tolerance)
{
    Require(x.size()==Size()&&b.size()==Size(),"mixed solve size mismatch");SolveInfo info;double begin=Now();
    if(!m_factor||relativeChange>.25||m_lastIterations>30||(m_initialIterations>0&&m_lastIterations>2*m_initialIterations)){info.factorSeconds+=PreparePreconditioner(a);info.rebuilds++;}
    int used=0;double correction=0;
    for(int restart=0;restart<14;++restart)
    {
        Vec r=(b-a*x).cwiseProduct(m_scale);double beta=r.norm();if(beta<tolerance*std::sqrt(m_volume))break;
        int dimension=std::min(30,400-used);Require(dimension>0,"mixed Krylov action cap exceeded");Mat V=Mat::Zero(Size(),dimension+1),Z(Size(),dimension),H=Mat::Zero(dimension+1,dimension);V.col(0)=r/beta;Vec g=Vec::Zero(dimension+1);g[0]=beta;Vec best;int done=0;
        for(int j=0;j<dimension;++j)
        {
            Z.col(j)=m_factor->solve(V.col(j).cwiseQuotient(m_scale));Require(m_factor->info()==Eigen::Success,"preconditioner application failed");Vec v=(a*Z.col(j)).cwiseProduct(m_scale);
            for(int pass=0;pass<2;++pass)for(int i=0;i<=j;++i){double h=V.col(i).dot(v);H(i,j)+=h;v.noalias()-=h*V.col(i);}H(j+1,j)=v.norm();if(H(j+1,j)>0)V.col(j+1)=v/H(j+1,j);
            best=H.block(0,0,j+2,j+1).colPivHouseholderQr().solve(g.head(j+2));double err=(g.head(j+2)-H.block(0,0,j+2,j+1)*best).norm();done=j+1;++used;
            if(err<tolerance*std::sqrt(m_volume)||H(j+1,j)<1e-28)break;
        }
        Vec delta=Z.leftCols(done)*best;correction+=delta.segment(2*m_nf,m_np).norm();x+=delta;
        auto trueR=Residual(a,b,x);if(std::max({trueR.momentum,trueR.continuity,trueR.gauge})<tolerance*20)break;
        if(used>=400)throw std::runtime_error("current native mixed system did not converge");
    }
    double t=Now();auto exact=Residual(a,b,x);info.residualSeconds=Now()-t;info.momentum=exact.momentum;info.continuity=exact.continuity;info.gauge=exact.gauge;info.iterations=used;info.correction=correction;info.solveSeconds=Now()-begin-info.factorSeconds-info.residualSeconds;
    if(!m_initialIterations&&used)m_initialIterations=used;m_lastIterations=used;
    Require(std::max({info.momentum,info.continuity,info.gauge})<=1e-9,"true current momentum/continuity/gauge residual failed");return info;
}
Vec Flow::Physical(const Vec&x,int c)const
{
    Vec out(m_fields[c]->GetTotPoints());for(const auto&e:m_elements){Vec coeff(c<2?16:9);for(int j=0;j<coeff.size();++j)coeff[j]=c<2?x[c*m_nu+e.u[j]]:x[2*m_nu+e.p[j]];Vec v=(c<2?e.B:e.P)*coeff;int offset=m_fields[c]->GetPhys_Offset(&e-&m_elements[0]);out.segment(offset,v.size())=v;}return out;
}
std::array<double,3>Flow::Error(const Vec&x,const std::array<Function,3>&f)const
{
    std::array<double,3>norm{};for(const auto&e:m_elements)for(int c=0;c<3;++c){Vec coeff(c<2?16:9);for(int j=0;j<coeff.size();++j)coeff[j]=c<2?x[c*m_nu+e.u[j]]:x[2*m_nu+e.p[j]];Vec v=(c<2?e.B:e.P)*coeff;for(int q=0;q<v.size();++q)norm[c]+=e.w[q]*std::pow(v[q]-f[c](e.x[q],e.z[q]),2);}for(double&v:norm)v=std::sqrt(v/m_volume);return norm;
}
double Flow::StrongDivergence(const Vec&x)const
{
    double norm=0;for(const auto&e:m_elements){Vec u(16),v(16);for(int j=0;j<16;++j){u[j]=x[e.u[j]];v[j]=x[m_nu+e.u[j]];}Vec div=e.Dx*u+e.Dz*v;norm+=div.cwiseProduct(div).dot(e.w);}return std::sqrt(norm/m_volume);
}
double Flow::VelocityDual(const Vec&r)const
{
    double value=0;for(int c=0;c<2;++c){Vec f=r.segment(c*m_nu+m_nd,m_nf),z=m_velocityRiesz.solve(f);value+=f.dot(z);}return std::sqrt(std::max(0.,value)/m_volume);
}
double Flow::PressureNorm(const Vec&p)const{return std::sqrt(std::max(0.,p.dot(m_qmass*p))/m_volume);}
}
