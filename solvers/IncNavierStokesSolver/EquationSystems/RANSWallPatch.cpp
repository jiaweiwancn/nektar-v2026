///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallPatch.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>

namespace Nektar
{
using namespace RANSWallAlgebra;
using namespace SolverUtils;
using WallMatrix=RANSWallAlgebra::DenseMatrix;
namespace
{
bool Equal(double a,double b){return std::abs(a-b)<1e-11;}
void Add(Vector &r,int offset,const RANSWallPoint &q,double value,double fx,double fr,double weight)
{
    for(int j=0;j<q.ids.size();++j)
        r[offset+q.ids[j]] += weight*(q.phi[j]*value+q.dx[j]*fx+q.dr[j]*fr);
}
}
RANSWallPatchData::RANSWallPatchData()
{
    RANSWallTimeFunction zero=[](double,double,double){return 0.;};
    innerBody.fill(zero);outerBody.fill(zero);boundary.fill(zero);
    innerKSource=outerKSource=outerEpsilonSource=zero;
}
RANSWallPatch::RANSWallPatch(std::shared_ptr<Layout> layout,const RANSWallPatchParameters &parameters)
    : m_layout(std::move(layout)),m_parameters(parameters)
{
    parameters.geometry.Validate();
    Require(Equal(parameters.geometry.height,m_layout->Height()),"layout/profile H mismatch");
    Require(parameters.tangentialOrder>=4 && parameters.normalOrder>=8,"underresolved patch quadrature request");
    m_options.requestedWeight=parameters.chi;m_options.quadratureOrder=parameters.normalOrder;
    std::set<std::pair<double,double>> columns;
    for(const auto &e:m_layout->innerVelocity->Elements())
    {
        columns.emplace(e.left,e.right);m_options.physicalBreakpoints.push_back(e.bottom);
        m_options.physicalBreakpoints.push_back(e.top);
    }
    auto rule=LibUtilities::PointsManager()[LibUtilities::PointsKey(parameters.tangentialOrder,LibUtilities::eGaussGaussLegendre)];
    const auto x=rule->GetZ(),w=rule->GetW();
    for(const auto &column:columns)for(int i=0;i<x.size();++i)
        m_rays.push_back({column.first+(x[i]+1)*(column.second-column.first)/2,
            w[i]*(column.second-column.first)/2,-1,nullptr,{}});
    for(const auto &q:m_layout->outerVelocity->Quadrature(parameters.tangentialOrder))
        m_outer.push_back({m_layout->velocity->Point(q.x,q.r),m_layout->k->Point(q.x,q.r),
                          m_layout->outerPressure->Point(q.x,q.r),m_layout->epsilon->Point(q.x,q.r),q.weight});
    const double volume=(m_layout->Right()-m_layout->Left())*m_layout->Top();
    for(int b=0;b<6;++b)m_mass[b]=m_layout->Field(static_cast<Layout::Block>(b)).Mass();
    m_pressureMean[0]=m_layout->innerPressure->NativeLoad([](double,double){return 1.;});
    m_pressureMean[1]=m_layout->outerPressure->NativeLoad([](double,double){return 1.;});
    for(auto &load:m_pressureMean)for(double &value:load)value/=volume;
}
std::shared_ptr<EquilibriumAWFProfiles> RANSWallPatch::Profile(double k) const
{
    return std::make_shared<EquilibriumAWFProfiles>(m_parameters.geometry,k,m_parameters.model,AWFNormalMetric{},m_options);
}
double RANSWallPatch::SampleK(const Vector &s,double x) const
{
    const auto &g=m_parameters.geometry;
    return g.SampleK([&](double r){return m_layout->k->Point(x,r,Equal(r,g.height)?-1:0).Value(s,m_layout->offset[Layout::K]);});
}
double RANSWallPatch::Forcing(const Vector &s,double x,double time,const RANSWallPatchData &data) const
{
    const auto &l=*m_layout;const double H=l.Height();
    const auto v=l.velocity->Point(x,H,1),p=l.outerPressure->Point(x,H,1),k=l.k->Point(x,H,1);
    const double u=v.Value(s,l.offset[Layout::U]);
    return 2*u*v.Gradient(s,l.offset[Layout::U])[0]+p.Gradient(s,l.offset[Layout::OuterPressure])[0]
           -2./3*k.Gradient(s,l.offset[Layout::K])[0]-data.outerBody[0](x,H,time);
}
Vector RANSWallPatch::Project(const std::array<RANSWallTimeFunction,6> &f,double time) const
{
    Vector result(m_layout->Size());
    for(int b=0;b<6;++b)
    {
        auto values=m_layout->Field(static_cast<Layout::Block>(b)).Project([&](double x,double r){return f[b](x,r,time);});
        std::copy(values.begin(),values.end(),result.begin()+m_layout->offset[b]);
    }
    double mean=0;
    for(int b=0;b<2;++b)for(int j=0;j<m_pressureMean[b].size();++j)
        mean+=m_pressureMean[b][j]*result[m_layout->offset[Layout::InnerPressure+b]+j];
    for(int b=Layout::InnerPressure;b<=Layout::OuterPressure;++b)
    {
        auto constant=m_layout->Field(static_cast<Layout::Block>(b)).Project([](double,double){return 1.;});
        for(int j=0;j<constant.size();++j)result[m_layout->offset[b]+j]-=mean*constant[j];
    }
    return result;
}
std::map<int,double> RANSWallPatch::FixedConstraints(double time,const RANSWallPatchData &data) const
{
    const auto &l=*m_layout;std::map<int,double> result;
    for(int b=0;b<4;++b)
    {
        const auto boundary=l.Field(static_cast<Layout::Block>(b)).Boundary(
            [&](double x,double r){return b<2 && Equal(r,0)?0.:data.boundary[b](x,r,time);},
            [&](int edge,const RANSWallNativeField::Element &e){
                return (edge==0 && b<2 && Equal(e.bottom,0)) || (edge==1 && Equal(e.right,l.Right())) ||
                    (edge==2 && Equal(e.top,l.Top())) || (edge==3 && Equal(e.left,l.Left())); });
        for(auto [id,value]:boundary)result[l.offset[b]+id]=value;
    }
    return result;
}
Vector RANSWallPatch::Lift(const Vector &state,double time,const RANSWallPatchData &data,
    std::map<int,double> &constraints) const
{
    Require(state.size()==static_cast<std::size_t>(m_layout->Size()),"state/layout mismatch");
    auto s=state;constraints=FixedConstraints(time,data);
    for(auto [row,value]:constraints)s[row]=value;
    const auto eps=m_layout->epsilon->Boundary([&](double x,double){return Profile(SampleK(s,x))->Dissipation(m_layout->Height());},
        [&](int edge,const RANSWallNativeField::Element &e){return edge==0 && Equal(e.bottom,m_layout->Height());});
    for(auto [id,value]:eps)
    {
        const int row=m_layout->offset[Layout::Epsilon]+id;
        if(constraints.count(row))Require(std::abs(constraints[row]-value)<1e-8*std::max(1.,std::abs(value)),"incompatible epsilon interface/external corner data");
        constraints[row]=value;s[row]=value;
    }
    return s;
}
std::map<int,double> RANSWallPatch::Constraints(const Vector &state,double time,const RANSWallPatchData &data) const
{std::map<int,double> constraints;Lift(state,time,data,constraints);return constraints;}

Vector RANSWallPatch::Residual(const Vector &state,const Vector &history,double alpha,double time,const RANSWallPatchData &data) const
{
    Require(alpha>=0 && (alpha==0 || history.size()==state.size()),"stage history/storage mismatch");
    const auto &l=*m_layout;std::map<int,double> constraints;
    const auto s=Lift(state,time,data,constraints);Vector residual(l.Size());
    for(auto &ray:m_rays)
    {
        const double k=SampleK(s,ray.x);
        if(k!=ray.k)
        {
            auto profile=Profile(k);std::vector<Point> points;
            for(const auto &q:profile->Quadrature())
                points.push_back({l.velocity->Point(ray.x,q.r),l.k->Point(ray.x,q.r),
                    l.innerPressure->Point(ray.x,q.r),{},ray.weight*q.weight*q.jacobian});
            ray.k=k;ray.profile=std::move(profile);ray.points=std::move(points);
        }
    }
    auto evaluate=[&](const Point &q,bool inner,const EquilibriumAWFProfiles *profile,double forcing)
    {
        const int ou=l.offset[Layout::U],ow=l.offset[Layout::W],ok=l.offset[Layout::K],oe=l.offset[Layout::Epsilon],
            op=l.offset[inner?Layout::InnerPressure:Layout::OuterPressure];
        const double u=q.v.Value(s,ou),w=q.v.Value(s,ow),k=q.k.Value(s,ok),pi=q.p.Value(s,op);
        const auto gu=q.v.Gradient(s,ou),gw=q.v.Gradient(s,ow),gk=q.k.Gradient(s,ok);
        Require(k>=0 && std::isfinite(k),"invalid transported k="+std::to_string(k)+
            " at x="+std::to_string(q.v.x)+", r="+std::to_string(q.v.r));
        RANSTensor gradient{};gradient[0][0]=gu[0];gradient[0][2]=gu[1];gradient[2][0]=gw[0];gradient[2][2]=gw[1];
        RANSTensor stress{};double dk,sk,de=0,se=0,ep=0;
        if(inner)
        {
            const auto &original=profile->Original();const double nt=original.EddyViscosity(q.v.r),d=original.TotalViscosity(q.v.r);
            double production=0;
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            {
                const double strain=gradient[i][j]+gradient[j][i];
                const double reynolds=(i==j?2*k/3:0)-nt*strain;
                production-=reynolds*gradient[i][j];stress[i][j]=d*strain-(i==j?pi:0);
            }
            dk=m_parameters.geometry.molecularNu+nt/m_parameters.model.Parameters().sigmaK;
            sk=profile->AtProduction(q.v.r,production).kSource+data.innerKSource(q.v.x,q.v.r,time);
        }
        else
        {
            ep=q.e.Value(s,oe);
            const auto c=EvaluateRANSClosure(m_parameters.model.Model(),k,ep,gradient,q.v.r,m_parameters.model.Parameters(),0);
            dk=c.kDiffusion;de=c.epsilonDiffusion;
            sk=c.kSource+data.outerKSource(q.v.x,q.v.r,time);
            se=c.epsilonSource+data.outerEpsilonSource(q.v.x,q.v.r,time);
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)
                stress[i][j]=-c.reynoldsStress[i][j]+(i==j?2*k/3-pi:0);
        }
        const auto &body=inner?data.innerBody:data.outerBody;
        const double fx=body[0](q.v.x,q.v.r,time)+(inner && q.v.r<=m_parameters.geometry.roughness?forcing:0),
            fr=body[1](q.v.x,q.v.r,time);
        auto storage=[&](const RANSWallPoint &p,int o,double value){return alpha==0?0:alpha*(value-p.Value(history,o));};
        Add(residual,ou,q.v,storage(q.v,ou,u)-fx,stress[0][0]-u*u,stress[0][2]-u*w,q.weight);
        Add(residual,ow,q.v,storage(q.v,ow,w)-fr,stress[2][0]-w*u,stress[2][2]-w*w,q.weight);
        Add(residual,ok,q.k,storage(q.k,ok,k)-sk,dk*gk[0]-u*k,dk*gk[1]-w*k,q.weight);
        Add(residual,op,q.p,gu[0]+gw[1],0,0,q.weight);
        if(!inner)
        {
            const auto ge=q.e.Gradient(s,oe);
            Add(residual,oe,q.e,storage(q.e,oe,ep)-se,de*ge[0]-u*ep,de*ge[1]-w*ep,q.weight);
        }
    };
    for(const auto &ray:m_rays)
    {
        const double forcing=Forcing(s,ray.x,time,data);
        for(const auto &q:ray.points)evaluate(q,true,ray.profile.get(),forcing);
    }
    for(const auto &q:m_outer)evaluate(q,false,nullptr,0);
    for(int b=0;b<2;++b)for(int j=0;j<m_pressureMean[b].size();++j)
    {
        const int row=l.offset[Layout::InnerPressure+b]+j;const double c=m_pressureMean[b][j];
        residual[row]+=s[l.offset[Layout::Gauge]]*c;
        residual[l.offset[Layout::Gauge]]+=s[row]*c;
    }
    for(auto [row,value]:constraints)residual[row]=state[row]-value;
    Maximum(residual);return residual;
}
RANSWallPatch::WallMatrix RANSWallPatch::Jacobian(const Vector &s,const Vector &history,double alpha,double time,
    const RANSWallPatchData &data,double step) const
{
    WallMatrix result(s.size(),Vector(s.size()));
    for(int j=0;j<s.size();++j)
    {
        auto plus=s,minus=s;const double h=step*std::max(1.,std::abs(s[j]));plus[j]+=h;minus[j]-=h;
        const auto a=Residual(plus,history,alpha,time,data),b=Residual(minus,history,alpha,time,data);
        for(int i=0;i<s.size();++i)result[i][j]=(a[i]-b[i])/(2*h);
    }
    return result;
}
std::vector<int> RANSWallPatch::CondensableVelocity(const std::map<int,double> &constraints) const
{
    const auto &l=*m_layout;std::set<int> outer;
    for(const auto &q:l.velocity->Quadrature())if(q.r>l.Height())outer.insert(q.ids.begin(),q.ids.end());
    std::vector<int> result;
    for(int b=0;b<2;++b)for(int j=0;j<l.size[b];++j)
        if(!outer.count(j) && !constraints.count(l.offset[b]+j))result.push_back(l.offset[b]+j);
    Require(!result.empty(),"no free inner velocity modes to condense");return result;
}
Vector RANSWallPatch::Solve(Vector s,const Vector &history,double alpha,double time,
    const RANSWallPatchData &data,bool condensed,int &iterations,double tolerance,int maximumIterations) const
{
    for(iterations=0;iterations<=maximumIterations;++iterations)
    {
        const auto r=Residual(s,history,alpha,time,data);const double error=Maximum(r);
        if(error<tolerance)return s;
        Require(iterations<maximumIterations,"nonlinear stage iteration limit");
        const auto jac=Jacobian(s,history,alpha,time,data);Vector step;
        if(condensed)step=RANSWallCondensation(jac,r,CondensableVelocity(Constraints(s,time,data))).Solve();
        else {auto rhs=r;for(double &x:rhs)x=-x;step=RANSWallAlgebra::Solve(jac,rhs);}
        bool accepted=false;
        for(double fraction=1;fraction>=1./1048576;fraction/=2)
        {
            auto trial=s;for(int j=0;j<s.size();++j)trial[j]+=fraction*step[j];
            try {if(Maximum(Residual(trial,history,alpha,time,data))<error){s=std::move(trial);accepted=true;break;}}
            catch(const std::exception &){ }
        }
        Require(accepted,"nonlinear stage line search rejected every admissible trial");
    }
    throw std::runtime_error("unreachable stage solve");
}
Vector RANSWallPatch::ResidualNorms(const Vector &r,const std::map<int,double> &constraints) const
{
    Vector norms;
    for(int b=0;b<6;++b)
    {
        Vector load;WallMatrix mass;
        for(int i=0;i<m_layout->size[b];++i)if(!constraints.count(m_layout->offset[b]+i))
        {
            load.push_back(r[m_layout->offset[b]+i]);mass.emplace_back();
            for(int j=0;j<m_layout->size[b];++j)if(!constraints.count(m_layout->offset[b]+j))mass.back().push_back(m_mass[b][i][j]);
        }
        if(load.empty()){norms.push_back(0);continue;}
        const auto inverse=RANSWallAlgebra::Solve(mass,load);
        const double square=std::inner_product(load.begin(),load.end(),inverse.begin(),0.);
        Require(square>=-1e-22,"invalid mass dual norm");norms.push_back(std::sqrt(std::max(0.,square)));
    }
    norms.push_back(std::abs(r[m_layout->offset[Layout::Gauge]]));
    double bc=0;for(auto [row,value]:constraints){(void)value;bc=std::max(bc,std::abs(r[row]));}norms.push_back(bc);
    return norms;
}
void RANSWallPatch::Store(const Vector &s) const
{
    for(int b=0;b<6;++b)m_layout->Field(static_cast<Layout::Block>(b)).Store(s,m_layout->offset[b]);
    const auto &l=*m_layout;
    for(const auto &sub:{l.innerVelocity,l.outerVelocity})
        sub->Store(sub->Project([&](double x,double r){return l.velocity->Point(x,r).Value(s,l.offset[Layout::U]);}));
}
RANSWallCondensation::RANSWallCondensation(const WallMatrix &j,const Vector &r,const std::vector<int> &inner):m_inner(inner)
{
    std::set<int> set(inner.begin(),inner.end());Require(set.size()==inner.size(),"duplicate condensed index");
    for(int i=0;i<r.size();++i)if(!set.count(i))m_retained.push_back(i);
    m_ii=WallMatrix(inner.size(),Vector(inner.size()));m_ir=WallMatrix(inner.size(),Vector(m_retained.size()));
    for(int i=0;i<inner.size();++i)
    {
        m_innerResidual.push_back(r[inner[i]]);
        for(int k=0;k<inner.size();++k)m_ii[i][k]=j[inner[i]][inner[k]];
        for(int k=0;k<m_retained.size();++k)m_ir[i][k]=j[inner[i]][m_retained[k]];
    }
    const auto inverse=RANSWallAlgebra::Solve(m_ii,m_innerResidual);
    m_schur=WallMatrix(m_retained.size(),Vector(m_retained.size()));m_residual=Vector(m_retained.size());
    for(int i=0;i<m_retained.size();++i)
    {
        m_residual[i]=r[m_retained[i]];
        for(int a=0;a<inner.size();++a)m_residual[i]-=j[m_retained[i]][inner[a]]*inverse[a];
    }
    for(int k=0;k<m_retained.size();++k)
    {
        Vector rhs;for(const auto &row:m_ir)rhs.push_back(row[k]);const auto solution=RANSWallAlgebra::Solve(m_ii,rhs);
        for(int i=0;i<m_retained.size();++i)
        {
            m_schur[i][k]=j[m_retained[i]][m_retained[k]];
            for(int a=0;a<inner.size();++a)m_schur[i][k]-=j[m_retained[i]][inner[a]]*solution[a];
        }
    }
}
Vector RANSWallCondensation::Recover(const Vector &retainedStep) const
{
    Require(retainedStep.size()==m_retained.size(),"invalid Schur vector");
    auto rhs=Multiply(m_ir,retainedStep);for(int i=0;i<rhs.size();++i)rhs[i]=-rhs[i]-m_innerResidual[i];
    const auto inner=RANSWallAlgebra::Solve(m_ii,rhs);Vector full(m_inner.size()+m_retained.size());
    for(int i=0;i<m_inner.size();++i)full[m_inner[i]]=inner[i];
    for(int i=0;i<m_retained.size();++i)full[m_retained[i]]=retainedStep[i];
    return full;
}
Vector RANSWallCondensation::Solve() const
{
    auto rhs=m_residual;for(double &x:rhs)x=-x;return Recover(RANSWallAlgebra::Solve(m_schur,rhs));
}
}
