///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <set>
#include <sstream>
#include <tuple>

namespace Nektar
{
using namespace RANSWallAlgebra;
using namespace SolverUtils;
namespace
{
using WallMatrix=RANSWallAlgebra::DenseMatrix;
bool Same(double a,double b){return std::abs(a-b)<1e-11;}
class Factor
{
public:
    explicit Factor(const WallMatrix &a):n(a.size()),lu(n*n),pivots(n)
    {
        Require(n>0,"empty factorisation");
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)lu[i+n*j]=a[i][j];
        int info;Lapack::Dgetrf(n,n,lu.data(),n,pivots.data(),info);
        Require(info==0,"singular evolution stage/mass matrix");
    }
    Vector Apply(Vector rhs) const
    {
        auto pivot=pivots;
        int info;Lapack::Dgetrs('N',n,1,lu.data(),n,pivot.data(),rhs.data(),n,info);
        Require(info==0,"failed evolution factor solve");Maximum(rhs);return rhs;
    }
private:
    int n;Vector lu;std::vector<int> pivots;
};
void Add(Vector &r,int offset,const RANSWallGlobalPoint &q,double value,
         const std::array<double,3> &flux,double weight)
{
    for(int j=0;j<q.ids.size();++j)
    {
        long double v=q.phi[j]*value;
        for(int a=0;a<3;++a)v+=(long double)q.derivative[a][j]*flux[a];
        r[offset+q.ids[j]]+=weight*static_cast<double>(v);
    }
}
double Volume(const RANSWallGlobalField &field)
{
    double volume=0;
    for(const auto &e:field.Elements())
    {
        double v=(e.upper[0]-e.lower[0])*(e.upper[2]-e.lower[2]);
        if(field.Dimension()==3)v*=e.upper[1]-e.lower[1];
        volume+=v;
    }
    return volume;
}
}
struct RANSWallEvolutionPatch::SolverCache
{
    WallMatrix jacobian;std::shared_ptr<Factor> factor;
    double alpha=-1;
    std::vector<std::vector<int>> free;
    std::vector<std::shared_ptr<Factor>> massFactors;
};
RANSWallEvolutionData::RANSWallEvolutionData(int dimension)
{
    RANSWallEvolutionFunction zero=[](const RANSWallCoordinate&,double){return 0.;};
    innerBody.assign(dimension,zero);outerBody.assign(dimension,zero);boundary.assign(dimension+2,zero);
    innerKSource=outerKSource=outerEpsilonSource=zero;
}
RANSWallEvolutionPatch::RANSWallEvolutionPatch(const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph,const LibUtilities::CommSharedPtr &comm,
    const RANSWallPatchParameters &parameters):m_dimension(graph->GetMeshDimension()),m_comm(comm),m_parameters(parameters)
{
    Require((m_dimension==2 || m_dimension==3) && session->GetComm()->GetSize()==1,
            "evolution requires a full serial native mesh on each collective rank");
    parameters.geometry.Validate();
    auto full=MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(session,graph,true,"u");
    auto pressure=MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(session,graph,true,"p");
    std::vector<unsigned int> all,inner,outer;
    lower.fill(1e100);upper.fill(-1e100);
    const int normal=m_dimension-1;
    for(int e=0;e<full->GetExpSize();++e)
    {
        const auto geom=full->GetExp(e)->GetGeom();double bottom=1e100,top=-1e100;
        for(int j=0;j<geom->GetNumVerts();++j)
        {
            double x,y,z;geom->GetVertex(j)->GetCoords(x,y,z);
            RANSWallCoordinate point{x,m_dimension==3?y:0.,m_dimension==3?z:y};
            bottom=std::min(bottom,point[2]);top=std::max(top,point[2]);
            for(int a=0;a<3;++a){lower[a]=std::min(lower[a],point[a]);upper[a]=std::max(upper[a],point[a]);}
        }
        Require(bottom>=-1e-12 && !(bottom<parameters.geometry.height-1e-12 && top>parameters.geometry.height+1e-12),
                "evolution element crosses ground or H");
        all.push_back(e);if(top<=parameters.geometry.height+1e-12)inner.push_back(e);else outer.push_back(e);
        Require(full->GetExp(e)->GetBasisNumModes(normal)>=4,"evolution velocity/scalar degree must be at least three");
    }
    Require(!inner.empty() && !outer.empty() && Same(lower[2],0),"incomplete wall/outer physical domains");
    auto velocity=std::make_shared<RANSWallGlobalField>(session,full,all,m_dimension);
    innerVelocity=std::make_shared<RANSWallGlobalField>(session,full,inner,m_dimension);
    outerVelocity=std::make_shared<RANSWallGlobalField>(session,full,outer,m_dimension);
    fields.assign(m_dimension,velocity);
    fields.push_back(std::make_shared<RANSWallGlobalField>(session,full,all,m_dimension));
    fields.push_back(std::make_shared<RANSWallGlobalField>(session,full,outer,m_dimension));
    fields.push_back(std::make_shared<RANSWallGlobalField>(session,pressure,inner,m_dimension));
    fields.push_back(std::make_shared<RANSWallGlobalField>(session,pressure,outer,m_dimension));
    offset.resize(Gauge()+1);size.resize(Gauge()+1);
    for(int b=0;b<Gauge();++b){size[b]=fields[b]->Size();offset[b+1]=offset[b]+size[b];}
    size[Gauge()]=1;
    Require(fields[K()]->Size()==velocity->Size(),"different velocity/k polynomial spaces are not supported");
    for(int b=0;b<Gauge();++b)
    {
        m_mass.push_back(b>0 && b<m_dimension?m_mass[0]:fields[b]->Mass());
        m_constants.push_back(fields[b]->Project([](const RANSWallCoordinate&){return 1.;}));
    }
    const double volume=Volume(*velocity);
    for(int side=0;side<2;++side)
    {
        const int b=InnerPressure()+side;m_pressureMean[side]=Multiply(m_mass[b],m_constants[b]);
        for(double &v:m_pressureMean[side])v/=volume;
    }
    m_options.requestedWeight=parameters.chi;m_options.quadratureOrder=parameters.normalOrder;
    std::set<std::array<double,4>> columns;
    for(const auto &e:innerVelocity->Elements())
    {
        columns.insert({e.lower[0],e.upper[0],e.lower[1],e.upper[1]});
        m_options.physicalBreakpoints.push_back(e.lower[2]);m_options.physicalBreakpoints.push_back(e.upper[2]);
    }
    const int order=parameters.tangentialOrder;
    Require(order>=4 && parameters.normalOrder>=8,"underresolved evolution quadrature");
    auto rule=LibUtilities::PointsManager()[LibUtilities::PointsKey(order,LibUtilities::eGaussGaussLegendre)];
    const auto z=rule->GetZ(),w=rule->GetW();
    for(const auto &column:columns)for(int j=0;j<(m_dimension==3?order:1);++j)for(int i=0;i<order;++i)
    {
        Ray ray;ray.coordinate={column[0]+(z[i]+1)*(column[1]-column[0])/2,
            m_dimension==3?column[2]+(z[j]+1)*(column[3]-column[2])/2:0.,parameters.geometry.sampleHeight};
        ray.weight=w[i]*(column[1]-column[0])/2;
        if(m_dimension==3)ray.weight*=w[j]*(column[3]-column[2])/2;
        ray.sample=fields[K()]->Point(ray.coordinate,Same(ray.coordinate[2],parameters.geometry.height)?-1:0);
        auto H=ray.coordinate;H[2]=parameters.geometry.height;
        ray.outerV=velocity->Point(H,1);ray.outerP=fields[OuterPressure()]->Point(H,1);ray.outerK=fields[K()]->Point(H,1);
        m_rays.push_back(std::move(ray));
    }
    for(const auto &q:outerVelocity->Quadrature(order))if(q.geometry%m_comm->GetSize()==m_comm->GetRank())
        m_outer.push_back({velocity->Point(q.coordinate),fields[OuterPressure()]->Point(q.coordinate),
            fields[Epsilon()]->Point(q.coordinate),q.weight});
    m_solver=std::make_shared<SolverCache>();m_solver->free.resize(Gauge());m_solver->massFactors.resize(Gauge());
    Vector seed(Size());for(int i=0;i<size[K()];++i)seed[offset[K()]+i]=parameters.model.Q()*.16*m_constants[K()][i];
    Refresh(seed);for(const auto &ray:m_rays)for(const auto &q:ray.points)m_bodyInner.push_back(q);
}
std::shared_ptr<EquilibriumAWFProfiles> RANSWallEvolutionPatch::Profile(double k) const
{return std::make_shared<EquilibriumAWFProfiles>(m_parameters.geometry,k,m_parameters.model,AWFNormalMetric{},m_options);}
void RANSWallEvolutionPatch::Broadcast(Vector &value) const {m_comm->Bcast(value,0);}
double RANSWallEvolutionPatch::SampleK(const Vector &s,const RANSWallCoordinate &point) const
{
    auto p=point;p[2]=m_parameters.geometry.sampleHeight;
    const double k=fields[K()]->Point(p,Same(p[2],m_parameters.geometry.height)?-1:0).Value(s,offset[K()]);
    Require(k>0 && std::isfinite(k),"invalid physical sampled K");return k;
}
std::array<double,2> RANSWallEvolutionPatch::Forcing(const Vector &s,const RANSWallCoordinate &point,double time,
    const RANSWallEvolutionData &data) const
{
    auto p=point;p[2]=m_parameters.geometry.height;
    const auto v=fields[0]->Point(p,1),pressure=fields[OuterPressure()]->Point(p,1),k=fields[K()]->Point(p,1);
    std::array<double,2> force{};
    const auto gp=pressure.Gradient(s,offset[OuterPressure()]),gk=k.Gradient(s,offset[K()]);
    for(int a=0;a<m_dimension-1;++a)
    {
        const int aa=PhysicalAxis(a);const double ua=v.Value(s,offset[a]);const auto ga=v.Gradient(s,offset[a]);
        force[a]=gp[aa]-2*gk[aa]/3-data.outerBody[a](p,time);
        for(int b=0;b<m_dimension-1;++b)
        {
            const int bb=PhysicalAxis(b);const auto gb=v.Gradient(s,offset[b]);
            force[a]+=v.Value(s,offset[b])*ga[bb]+ua*gb[bb];
        }
    }
    return force;
}
void RANSWallEvolutionPatch::Refresh(const Vector &s) const
{
    for(auto &ray:m_rays)
    {
        const double k=ray.sample.Value(s,offset[K()]);
        if(k==ray.cachedK)continue;
        auto profile=Profile(k);std::vector<VolumePoint> points;
        for(const auto &q:profile->Quadrature())
        {
            auto p=ray.coordinate;p[2]=q.r;auto v=fields[0]->Point(p);
            if(v.geometry%m_comm->GetSize()!=m_comm->GetRank())continue;
            const auto kp=fields[K()]->Point(p);Require(kp.ids==v.ids,"native k/velocity maps differ");
            points.push_back({std::move(v),fields[InnerPressure()]->Point(p),{},ray.weight*q.weight*q.jacobian});
        }
        ray.cachedK=k;ray.profile=std::move(profile);ray.points=std::move(points);
    }
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::Project(const std::vector<RANSWallEvolutionFunction> &f,double time) const
{
    Require(f.size()==fields.size(),"incomplete evolution initial/reference functions");Vector s(Size());
    for(int b=0;b<Gauge();++b)
    {
        const auto values=fields[b]->Project([&](const RANSWallCoordinate &p){return f[b](p,time);});
        std::copy(values.begin(),values.end(),s.begin()+offset[b]);
    }
    double mean=0;for(int side=0;side<2;++side)for(int i=0;i<m_pressureMean[side].size();++i)
        mean+=m_pressureMean[side][i]*s[offset[InnerPressure()+side]+i];
    for(int b=InnerPressure();b<=OuterPressure();++b)for(int i=0;i<size[b];++i)s[offset[b]+i]-=mean*m_constants[b][i];
    Broadcast(s);return s;
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::NeutralReference(double us,const std::array<double,2> &direction) const
{
    Require(us>0 && std::abs(direction[0]*direction[0]+direction[1]*direction[1]-1)<1e-12,"invalid reference friction velocity/direction");
    Require(m_dimension==3 || std::abs(direction[1])<1e-15,"2D reference direction must be streamwise");
    const double K0=m_parameters.model.Q()*us*us,H=m_parameters.geometry.height,Z=upper[2];
    const auto profile=Profile(K0);const double UH=profile->Reference().matchingVelocity;
    const auto speed=[&](double r){return r<=H?profile->Original().Velocity(r,{UH,0},{0,0})[0]:UH+us/m_parameters.model.Kappa()*std::log(r/H);};
    const auto velocity=fields[0]->RadialProject(speed,{{0,0,0},{H,0,UH},
        {H,-1,us*us/profile->Original().TotalViscosity(H)},{H,1,m_parameters.model.Q()*profile->Dissipation(H)/K0},
        {Z,0,speed(Z)},{Z,-1,us/(m_parameters.model.Kappa()*Z)}});
    Vector s(Size());
    for(int a=0;a<m_dimension-1;++a)for(int i=0;i<size[a];++i)s[offset[a]+i]=direction[a]*velocity[i];
    for(int i=0;i<size[K()];++i)s[offset[K()]+i]=K0*m_constants[K()][i];
    const double A=us*us*us/m_parameters.model.Kappa();
    const auto ep=fields[Epsilon()]->RadialProject([&](double r){return A/r;},
        {{H,0,profile->Dissipation(H)},{Z,0,A/Z},{Z,-1,-A/(Z*Z)}});
    std::copy(ep.begin(),ep.end(),s.begin()+offset[Epsilon()]);
    RANSTensor G{};for(int a=0;a<m_dimension-1;++a)G[PhysicalAxis(a)][2]=direction[a]*us/(m_parameters.model.Kappa()*H);
    const auto closure=EvaluateRANSClosure(m_parameters.model.Model(),K0,A/H,G,H,m_parameters.model.Parameters(),0);
    const double jump=closure.reynoldsStress[2][2]-2*K0/3;
    for(int b=InnerPressure();b<=OuterPressure();++b)
    {
        const double pressure=b==InnerPressure()?jump*(Z-H)/Z:-jump*H/Z;
        for(int i=0;i<size[b];++i)s[offset[b]+i]=pressure*m_constants[b][i];
    }
    Broadcast(s);return s;
}
std::map<int,double> RANSWallEvolutionPatch::ExternalConstraints(double time,const RANSWallEvolutionData &data) const
{
    if(!data.fixedBoundary.empty())return data.fixedBoundary;
    std::map<int,double> result;
    for(int b=0;b<InnerPressure();++b)
    {
        const auto values=fields[b]->Boundary([&](const RANSWallCoordinate &p){return b<m_dimension && Same(p[2],0)?0.:data.boundary[b](p,time);},
            [&](int axis,int side,const RANSWallGlobalField::Element &e){
                if(axis!=2)return Same(side?e.upper[axis]:e.lower[axis],side?upper[axis]:lower[axis]);
                return side?Same(e.upper[2],upper[2]):(b<m_dimension && Same(e.lower[2],0));});
        for(auto [id,value]:values)result[offset[b]+id]=value;
    }
    return result;
}
std::map<int,double> RANSWallEvolutionPatch::FixedReferenceBoundary(const Vector &ref) const
{
    RANSWallEvolutionData zero(m_dimension);auto result=ExternalConstraints(0,zero);
    for(auto &[row,value]:result)value=ref[row];
    return result;
}
std::map<int,double> RANSWallEvolutionPatch::Constraints(const Vector &state,double time,const RANSWallEvolutionData &data) const
{
    auto fixed=ExternalConstraints(time,data);auto s=state;for(auto [row,value]:fixed)s[row]=value;
    const auto ep=fields[Epsilon()]->Boundary([&](const RANSWallCoordinate &p){return Profile(SampleK(s,p))->Dissipation(m_parameters.geometry.height);},
        [&](int axis,int side,const RANSWallGlobalField::Element &e){return axis==2 && side==0 && Same(e.lower[2],m_parameters.geometry.height);});
    for(auto [id,value]:ep)
    {
        const int row=offset[Epsilon()]+id;
        if(fixed.count(row))Require(std::abs(fixed[row]-value)<1e-9*std::max(1.,std::abs(value)),"incompatible epsilon corner data");
        fixed[row]=value;
    }
    return fixed;
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::Lift(const Vector &state,double time,const RANSWallEvolutionData &data) const
{auto s=state;for(auto [row,value]:Constraints(state,time,data))s[row]=value;return s;}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::Residual(const Vector &state,const Vector &history,double alpha,double time,
    const RANSWallEvolutionData &data,bool anchored,bool imposeConstraints) const
{
    Require(state.size()==static_cast<size_t>(Size()) && history.size()==state.size() && alpha>=0,"invalid evolution stage state/history");
    const auto constraints=Constraints(state,time,data);auto s=state;for(auto [row,value]:constraints)s[row]=value;
    Vector r(Size());int bad=0;
    try
    {
        Refresh(s);
        auto evaluate=[&](const VolumePoint &q,bool inner,const EquilibriumAWFProfiles *profile,const std::array<double,2> &force)
        {
            RANSTensor G{},stress{};std::array<double,3> u{};
            for(int a=0;a<m_dimension;++a){u[PhysicalAxis(a)]=q.v.Value(s,offset[a]);G[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);}
            const double k=q.v.Value(s,offset[K()]),pi=q.p.Value(s,offset[inner?InnerPressure():OuterPressure()]);
            const auto gk=q.v.Gradient(s,offset[K()]);Require(k>=0 && std::isfinite(k),"negative/nonfinite transported k");
            double dk,sk,de=0,se=0,epsilon=0;
            if(inner)
            {
                const double nt=profile->Original().EddyViscosity(q.v.coordinate[2]);double production=0;
                for(int a=0;a<3;++a)for(int b=0;b<3;++b)
                {
                    const double strain=G[a][b]+G[b][a];production-=((a==b?2*k/3:0)-nt*strain)*G[a][b];
                    stress[a][b]=(m_parameters.geometry.molecularNu+nt)*strain-(a==b?pi:0);
                }
                dk=m_parameters.geometry.molecularNu+nt/m_parameters.model.Parameters().sigmaK;
                sk=profile->AtProduction(q.v.coordinate[2],production).kSource+data.innerKSource(q.v.coordinate,time);
            }
            else
            {
                epsilon=q.e.Value(s,offset[Epsilon()]);
                const auto c=EvaluateRANSClosure(m_parameters.model.Model(),k,epsilon,G,q.v.coordinate[2],m_parameters.model.Parameters(),0);
                dk=c.kDiffusion;de=c.epsilonDiffusion;sk=c.kSource+data.outerKSource(q.v.coordinate,time);
                se=c.epsilonSource+data.outerEpsilonSource(q.v.coordinate,time);
                for(int a=0;a<3;++a)for(int b=0;b<3;++b)stress[a][b]=-c.reynoldsStress[a][b]+(a==b?2*k/3-pi:0);
            }
            const auto &body=inner?data.innerBody:data.outerBody;
            for(int a=0;a<m_dimension;++a)
            {
                const int aa=PhysicalAxis(a);std::array<double,3> flux{};
                for(int b=0;b<3;++b)flux[b]=stress[aa][b]-u[aa]*u[b];
                const double rough=inner && a<m_dimension-1 && q.v.coordinate[2]<=m_parameters.geometry.roughness?force[a]:0;
                Add(r,offset[a],q.v,-body[a](q.v.coordinate,time)-rough,flux,q.weight);
            }
            std::array<double,3> kflux{};for(int b=0;b<3;++b)kflux[b]=dk*gk[b]-u[b]*k;
            Add(r,offset[K()],q.v,-sk,kflux,q.weight);
            Add(r,offset[inner?InnerPressure():OuterPressure()],q.p,G[0][0]+G[1][1]+G[2][2],{},q.weight);
            if(!inner)
            {
                const auto ge=q.e.Gradient(s,offset[Epsilon()]);std::array<double,3> flux{};
                for(int b=0;b<3;++b)flux[b]=de*ge[b]-u[b]*epsilon;
                Add(r,offset[Epsilon()],q.e,-se,flux,q.weight);
            }
        };
        for(const auto &ray:m_rays)
        {
            const auto force=Forcing(s,ray.coordinate,time,data);
            for(const auto &q:ray.points)evaluate(q,true,ray.profile.get(),force);
        }
        for(const auto &q:m_outer)evaluate(q,false,nullptr,{});
    }
    catch(const std::exception &){bad=1;}
    m_comm->AllReduce(bad,LibUtilities::ReduceMax);Require(!bad,"inadmissible evolution volume/sampler state on a collective rank");
    m_comm->AllReduce(r,LibUtilities::ReduceSum);
    if(alpha)for(int b=0;b<InnerPressure();++b)
    {
        Vector delta(size[b]);for(int i=0;i<size[b];++i)delta[i]=s[offset[b]+i]-history[offset[b]+i];
        const auto load=Multiply(m_mass[b],delta);for(int i=0;i<size[b];++i)r[offset[b]+i]+=alpha*load[i];
    }
    for(int side=0;side<2;++side)for(int i=0;i<m_pressureMean[side].size();++i)
    {
        const int row=offset[InnerPressure()+side]+i;const double w=m_pressureMean[side][i];
        r[row]+=s[offset[Gauge()]]*w;r[offset[Gauge()]]+=s[row]*w;
    }
    if(anchored && m_eta)for(int i=0;i<DifferentialSize();++i)r[i]-=m_eta*m_defect[i];
    if(imposeConstraints)for(auto [row,value]:constraints)r[row]=state[row]-value;
    Maximum(r);return r;
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::Norms(const Vector &r,const std::map<int,double> &constraints,double us) const
{
    Vector norms;const double K0=m_parameters.model.Q()*us*us,E0=us*us*us/(m_parameters.model.Kappa()*m_parameters.geometry.height);
    for(int b=0;b<Gauge();++b)
    {
        std::vector<int> free;Vector load;
        for(int i=0;i<size[b];++i)if(!constraints.count(offset[b]+i)){free.push_back(i);load.push_back(r[offset[b]+i]);}
        if(free.empty()){norms.push_back(0);continue;}
        if(free!=m_solver->free[b] || !m_solver->massFactors[b])
        {
            WallMatrix mass(free.size(),Vector(free.size()));
            for(int i=0;i<free.size();++i)for(int j=0;j<free.size();++j)mass[i][j]=m_mass[b][free[i]][free[j]];
            m_solver->free[b]=free;m_solver->massFactors[b]=std::make_shared<Factor>(mass);
        }
        const auto inverse=m_solver->massFactors[b]->Apply(load);long double square=0;
        for(int i=0;i<load.size();++i)square+=(long double)load[i]*inverse[i];
        Require(square>=-1e-21,"invalid evolution mass dual norm");
        const double scale=b<m_dimension?us:(b==K()?K0:(b==Epsilon()?E0:1.));
        const double rate=b>=InnerPressure()?us/m_parameters.geometry.height:scale*us/m_parameters.geometry.height;
        norms.push_back(std::sqrt(std::max(0.L,square))/(rate*std::sqrt(Volume(*fields[b]))));
    }
    norms.push_back(std::abs(r[offset[Gauge()]])/(us*us));double bc=0;
    for(auto [row,value]:constraints)
    {
        (void)value;int b=0;while(b+1<offset.size() && row>=offset[b+1])++b;
        const double scale=b<m_dimension?us:(b==K()?K0:E0);bc=std::max(bc,std::abs(r[row])/scale);
    }
    norms.push_back(bc);return norms;
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::StateNorms(const Vector &s,const Vector &reference,double us) const
{
    Vector norms;
    for(int b=0;b<Gauge();++b)
    {
        Vector delta(size[b]);for(int i=0;i<size[b];++i)delta[i]=s[offset[b]+i]-reference[offset[b]+i];
        const auto load=Multiply(m_mass[b],delta);long double sum=0;for(int i=0;i<size[b];++i)sum+=(long double)delta[i]*load[i];
        const double scale=b<m_dimension?us:(b==K()?m_parameters.model.Q()*us*us:
            (b==Epsilon()?us*us*us/(m_parameters.model.Kappa()*m_parameters.geometry.height):us*us));
        norms.push_back(std::sqrt(std::max(0.L,sum))/(scale*std::sqrt(Volume(*fields[b]))));
    }
    return norms;
}
std::vector<int> RANSWallEvolutionPatch::CondensableVelocity(const std::map<int,double> &constraints) const
{
    std::set<int> outerIds;for(const auto &q:outerVelocity->Quadrature(4))
    {const auto v=fields[0]->Point(q.coordinate);outerIds.insert(v.ids.begin(),v.ids.end());}
    std::vector<int> result;for(int b=0;b<m_dimension;++b)for(int i=0;i<size[b];++i)
        if(!outerIds.count(i) && !constraints.count(offset[b]+i))result.push_back(offset[b]+i);
    Require(!result.empty(),"no free inner velocity modes to condense");return result;
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::Solve(Vector s,const Vector &history,double alpha,double time,
    const RANSWallEvolutionData &data,bool condensed,int &iterations,double tolerance,int maximumIterations) const
{
    Require(Size()<=3000,"replicated dense evolution solve is limited to 3000 unknowns; larger meshes permit reference audits only");
    bool refreshed=false;double previous=1e100;
    for(iterations=0;iterations<=maximumIterations;++iterations)
    {
        const auto r=Residual(s,history,alpha,time,data);const auto constraints=Constraints(s,time,data);
        const double error=Maximum(Norms(r,constraints));if(error<tolerance)return s;
        Require(iterations<maximumIterations,"evolution stage nonlinear iteration limit, scaled residual="+std::to_string(error));
        if(!m_solver->factor || !Same(alpha,m_solver->alpha) || (error>previous*.8 && !refreshed))
        {
            if(m_comm->GetRank()==0)std::cout<<"EP-AWF stage Jacobian at time "<<time<<", alpha "<<alpha<<", scaled residual "<<error<<std::endl;
            WallMatrix jac(Size(),Vector(Size()));
            for(int j=0;j<Size();++j)
            {
                // Lifted essential/dependent columns have no physical-load
                // variation. Their only nonzero entry is their constraint row.
                // Dependence of the epsilon trace on free k columns is still
                // differentiated below through the full coupled residual.
                if(constraints.count(j)){jac[j][j]=1.;continue;}
                auto plus=s,minus=s;const double h=2e-6*std::max(1.,std::abs(s[j]));plus[j]+=h;minus[j]-=h;
                const auto a=Residual(plus,history,alpha,time,data),b=Residual(minus,history,alpha,time,data);
                for(int i=0;i<Size();++i)jac[i][j]=(a[i]-b[i])/(2*h);
            }
            m_solver->jacobian=std::move(jac);m_solver->alpha=alpha;
            int bad=0;
            if(m_comm->GetRank()==0)
            {
                try{m_solver->factor=std::make_shared<Factor>(m_solver->jacobian);}
                catch(const std::exception &){bad=1;}
            }
            m_comm->Bcast(bad,0);Require(!bad,"collective stage factorisation failed");
            if(m_comm->GetRank()!=0)m_solver->factor=std::make_shared<Factor>(WallMatrix{{1.}});
            refreshed=true;
        }
        Vector step(Size());
        int bad=0;
        if(m_comm->GetRank()==0)
        {
            try
            {
                if(condensed)step=RANSWallCondensation(m_solver->jacobian,r,CondensableVelocity(constraints)).Solve();
                else{auto rhs=r;for(double &v:rhs)v=-v;step=m_solver->factor->Apply(std::move(rhs));}
            }
            catch(const std::exception &){bad=1;}
        }
        m_comm->Bcast(bad,0);Require(!bad,"collective Newton/Schur solve failed");
        Broadcast(step);bool accepted=false;
        for(double fraction=1;fraction>=1./1048576;fraction/=2)
        {
            auto trial=s;for(int i=0;i<Size();++i)trial[i]+=fraction*step[i];
            try
            {
                const auto next=Residual(trial,history,alpha,time,data);
                if(Maximum(Norms(next,Constraints(trial,time,data)))<error){s=std::move(trial);accepted=true;break;}
            }
            catch(const std::exception &){}
        }
        if(!accepted && !refreshed){m_solver->factor.reset();previous=1e100;continue;}
        Require(accepted,"evolution stage rejected every admissible Newton trial");previous=error;
    }
    throw std::runtime_error("unreachable evolution solve");
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::ExplicitBody(const Vector &s,double time,const RANSWallEvolutionData &data) const
{
    Vector load(Size());
    auto add=[&](const VolumePoint &q,bool inner)
    {
        auto H=q.v.coordinate;H[2]=m_parameters.geometry.height;
        for(int a=0;a<m_dimension;++a)
        {
            double body=(inner?data.innerBody:data.outerBody)[a](q.v.coordinate,time);
            if(inner && a<m_dimension-1 && q.v.coordinate[2]<=m_parameters.geometry.roughness)body-=data.outerBody[a](H,time);
            Add(load,offset[a],q.v,body,{},q.weight);
        }
    };
    for(const auto &q:m_bodyInner)add(q,true);
    for(const auto &q:m_outer)add(q,false);
    m_comm->AllReduce(load,LibUtilities::ReduceSum);Vector result(DifferentialSize());
    const auto constraints=Constraints(s,time,data);Norms(load,constraints);
    for(int b=0;b<m_dimension;++b)
    {
        const auto &free=m_solver->free[b];Vector rhs;for(int i:free)rhs.push_back(load[offset[b]+i]);
        if(!rhs.empty())
        {const auto values=m_solver->massFactors[b]->Apply(rhs);for(int j=0;j<free.size();++j)result[offset[b]+free[j]]=values[j];}
    }
    return result;
}
void RANSWallEvolutionPatch::SetReference(const Vector &ref,const RANSWallEvolutionData &data,double time,double eta)
{
    Require(eta==0 || eta==1,"only fixed global eta=0 or eta=1 is supported");
    const auto raw=Residual(ref,ref,0,time,data,false);const auto norms=Norms(raw,Constraints(ref,time,data));
    for(int b=InnerPressure();b<norms.size();++b)Require(norms[b]<1e-9,"reference violates an algebraic constraint");
    m_reference=ref;m_referenceLoads=Residual(ref,ref,0,time,data,false,false);m_defect=m_referenceLoads;
    m_referenceConstraints=Constraints(ref,time,data);
    for(int i=DifferentialSize();i<Size();++i)m_defect[i]=0;
    for(auto [row,value]:Constraints(ref,time,data)){(void)value;m_defect[row]=0;}
    m_eta=eta;
}
void RANSWallEvolutionPatch::RestoreReference(const Vector &ref,const Vector &defect,double eta,
    const RANSWallEvolutionData &data,double time)
{
    SetReference(ref,data,time,eta);Require(defect.size()==m_defect.size(),"checkpoint defect shape differs");
    for(int i=0;i<Size();++i)Require(std::abs(defect[i]-m_defect[i])<1e-10*std::max(1.,std::abs(m_defect[i])),"checkpoint fixed defect differs");
    m_defect=defect;
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::DefectIntegrals() const
{
    Vector integrals(InnerPressure());if(m_defect.empty())return integrals;
    for(int b=0;b<InnerPressure();++b)for(int i=0;i<size[b];++i)integrals[b]+=m_defect[offset[b]+i]*m_constants[b][i];
    return integrals;
}
double RANSWallEvolutionPatch::AnchorEnergyWork(const Vector &s) const
{
    double work=0;if(m_defect.empty())return work;
    for(int b=0;b<m_dimension;++b)for(int i=0;i<size[b];++i)work+=s[offset[b]+i]*m_defect[offset[b]+i];
    for(int i=0;i<size[K()];++i)work+=m_constants[K()][i]*m_defect[offset[K()]+i];
    return m_eta*work;
}
std::string RANSWallEvolutionPatch::Signature() const
{
    std::ostringstream text;text<<std::setprecision(17)<<"CLS-EPAWF-M4:replicated-CG:"<<m_dimension<<':'<<m_parameters.geometry.height
        <<':'<<m_parameters.geometry.molecularNu<<':'<<m_parameters.geometry.roughness<<':'<<m_parameters.geometry.sampleHeight
        <<':'<<m_parameters.geometry.tailStart<<':'<<m_parameters.chi<<':'<<m_parameters.model.Q()<<':'<<m_parameters.tangentialOrder<<':'<<m_parameters.normalOrder;
    for(int b=0;b<Gauge();++b)
    {
        text<<"|"<<size[b];for(const auto &e:fields[b]->Elements())
        {text<<':'<<e.geometry<<':'<<e.degree;for(int a=0;a<3;++a)text<<':'<<e.lower[a]<<':'<<e.upper[a];}
        for(int e=0;e<fields[b]->Field()->GetExpSize();++e)
        {
            const auto exp=fields[b]->Field()->GetExp(e);
            for(int d=0;d<m_dimension;++d)text<<':'<<exp->GetBasis(d)->GetBasisType()<<':'<<exp->GetBasis(d)->GetNumModes()
                <<':'<<exp->GetBasis(d)->GetPointsType()<<':'<<exp->GetBasis(d)->GetNumPoints();
        }
    }
    return text.str();
}
}
