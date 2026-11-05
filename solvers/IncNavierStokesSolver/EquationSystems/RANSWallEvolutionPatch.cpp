///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numeric>
#include <limits>
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

}
struct RANSWallEvolutionPatch::SolverCache
{
    WallMatrix jacobian;std::shared_ptr<Factor> factor;
    double alpha=-1;
    double buildTime=-1;
    int builds=0,uses=0,buildIteration=-1;
    Vector buildState,buildHistory;
    std::vector<std::vector<int>> free;
    std::vector<std::shared_ptr<Factor>> massFactors;
};
RANSWallEvolutionData::RANSWallEvolutionData(int dimension)
{
    RANSWallEvolutionFunction zero=[](const RANSWallCoordinate&,double){return 0.;};
    innerBody.assign(dimension,zero);outerBody.assign(dimension,zero);boundary.assign(dimension+2,zero);
    innerKSource=outerKSource=outerEpsilonSource=zero;
}
RANSWallNormalMap LoadRANSWallNormalMap(const LibUtilities::SessionReaderSharedPtr &session,int dimension,double height)
{
    RANSWallNormalMap m;m.dimension=dimension;m.height=height;std::string name;
    session->LoadSolverInfo("RANSWallMap",name,"Flat");
    const std::map<std::string,RANSWallNormalMap::Kind> kinds{{"Flat",RANSWallNormalMap::Kind::Flat},{"Plane",RANSWallNormalMap::Kind::Plane},
        {"Cylinder",RANSWallNormalMap::Kind::Cylinder},{"Gaussian",RANSWallNormalMap::Kind::Gaussian},{"CosineHill",RANSWallNormalMap::Kind::CosineHill}};
    Require(kinds.count(name),"unknown wall map");m.kind=kinds.at(name);
    for(auto entry:std::vector<std::pair<std::string,double*>>{{"SlopeX",&m.slopeX},{"SlopeY",&m.slopeY},{"Radius",&m.radius},
        {"Amplitude",&m.amplitude},{"LengthX",&m.lengthX},{"LengthY",&m.lengthY},{"HeightSlopeX",&m.heightSlopeX},
        {"HeightSlopeY",&m.heightSlopeY},{"HeightAmplitude",&m.heightAmplitude},{"HeightWaveX",&m.heightWaveX},{"HeightWaveY",&m.heightWaveY}})
        session->LoadParameter("RANSWallMap"+entry.first,*entry.second,*entry.second);
    return m;
}
RANSWallCoordinate RANSWallEvolutionPatch::Physical(const RANSWallCoordinate &p) const
{return Mapped()?m_mapping.At(p).point:p;}
double RANSWallEvolutionPatch::Distance(const RANSWallCoordinate &p) const
{return Mapped()?p[2]*(m_mapping.height+m_mapping.heightSlopeX*p[0]+m_mapping.heightSlopeY*p[1]+
    m_mapping.heightAmplitude*std::sin(m_mapping.heightWaveX*p[0])*std::cos(m_mapping.heightWaveY*p[1]))/m_mapping.height:p[2];}
RANSWallGlobalPoint RANSWallEvolutionPatch::Transform(RANSWallGlobalPoint q) const
{
    if(!Mapped())return q;
    int side=0;
    if(m_mapping.kind==RANSWallNormalMap::Kind::CosineHill)
        for(double foot:{-m_mapping.lengthX,m_mapping.lengthX})if(Same(q.coordinate[0],foot))
            for(const auto &e:fields[0]->Elements())if(e.geometry==q.geometry)side=(e.lower[0]+e.upper[0])/2<foot?-1:1;
    const auto g=m_mapping.At(q.coordinate,side);const auto old=q.derivative;
    for(int a=0;a<3;++a)for(int i=0;i<q.ids.size();++i)
    {q.derivative[a][i]=0;for(int b=0;b<3;++b)q.derivative[a][i]+=old[b][i]*g.inverse[b][a];}
    q.weight*=g.determinant;return q;
}
RANSWallGlobalPoint RANSWallEvolutionPatch::Point(int block,const RANSWallCoordinate &p,int side) const
{return Transform(fields[block]->Point(p,side));}
std::vector<RANSWallGlobalPoint> RANSWallEvolutionPatch::Quadrature(int block,int order) const
{
    if(!order && Mapped())order=std::max(m_parameters.normalOrder,m_parameters.tangentialOrder);
    auto q=fields[block]->Quadrature(order);for(auto &p:q)p=Transform(std::move(p));return q;
}
std::vector<RANSWallEvolutionPatch::FacePoint> RANSWallEvolutionPatch::Surface(int block,int axis,double position,int side,int order) const
{
    Require((axis==0 || axis==2 || (m_dimension==3 && axis==1)) && (side==-1 || side==1),"invalid chart face");
    if(!order)order=std::max(m_parameters.normalOrder,m_parameters.tangentialOrder);
    const auto rule=LibUtilities::PointsManager()[LibUtilities::PointsKey(order,LibUtilities::eGaussGaussLegendre)];
    const auto z=rule->GetZ(),w=rule->GetW();std::vector<int> free;
    for(int a=0;a<3;++a)if(a!=axis && (m_dimension==3 || a!=1))free.push_back(a);
    std::vector<FacePoint> result;
    for(int e=0;e<fields[block]->Elements().size();++e)
    {
        const auto &element=fields[block]->Elements()[e];if(!Same(side<0?element.lower[axis]:element.upper[axis],position))continue;
        for(int j=0;j<(m_dimension==3?order:1);++j)for(int i=0;i<order;++i)
        {
            RANSWallCoordinate p{};p[axis]=position;double weight=1;const std::array<int,2> index{i,j};
            for(int a=0;a<free.size();++a){const int d=free[a];p[d]=element.lower[d]+(z[index[a]]+1)*(element.upper[d]-element.lower[d])/2;weight*=w[index[a]]*(element.upper[d]-element.lower[d])/2;}
            int footSide=0;if(m_mapping.kind==RANSWallNormalMap::Kind::CosineHill)
                for(double foot:{-m_mapping.lengthX,m_mapping.lengthX})if(Same(p[0],foot))footSide=(element.lower[0]+element.upper[0])/2<foot?-1:1;
            result.push_back({Transform(fields[block]->ElementPoint(e,p)),m_mapping.Cofactor(p,axis,side,footSide),weight});
        }
    }
    return result;
}
double RANSWallEvolutionPatch::DependentEpsilon(const Vector &s,const RANSWallCoordinate &p) const
{
    Require(!m_transitionEnabled,"transition dissipation evaluation requires its current time and forcing data");
    return WallDissipation(s,p,*Profile(SampleK(s,p),p));
}
RANSWallEvolutionPatch::WallMatrix RANSWallEvolutionPatch::PhysicalMass(int block,int region) const
{
    if(!Mapped() && region<0)return fields[block]->Mass();
    const bool regional=region>=0 && (block<m_dimension || block==K());
    if(regional && !m_regionalMass[region].empty())return m_regionalMass[region];
    const int n=fields[block]->Size();WallMatrix mass(n,Vector(n));
    for(const auto &q:Quadrature(block))
    {
        if(region==0 && q.coordinate[2]>m_parameters.geometry.height)continue;
        if(region==1 && q.coordinate[2]<m_parameters.geometry.height)continue;
        for(int i=0;i<q.ids.size();++i)for(int j=0;j<q.ids.size();++j)
            mass[q.ids[i]][q.ids[j]]+=q.weight*q.phi[i]*q.phi[j];
    }
    if(regional)m_regionalMass[region]=mass;
    return mass;
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::PhysicalProject(int block,const RANSWallEvolutionFunction &f,double time) const
{
    const auto function=[&](const RANSWallCoordinate &p){return f(p,time);};
    auto values=fields[block]->Project(function);
    if(!Mapped() || m_mapping.PlanarConstantHeight())return values;
    const auto quadrature=Quadrature(block);
    bool constant=true;const double first=function(quadrature.front().coordinate);
    for(const auto &q:quadrature)if(function(q.coordinate)!=first){constant=false;break;}
    if(constant)return values;
    for(int pass=0;pass<3;++pass)
    {
        Vector load(values.size());for(const auto &q:quadrature)
        {const double defect=function(q.coordinate)-q.Value(values);for(int j=0;j<q.ids.size();++j)load[q.ids[j]]+=q.weight*q.phi[j]*defect;}
        const auto correction=RANSWallAlgebra::Solve(m_mass[block],load);
        for(int j=0;j<values.size();++j)values[j]+=correction[j];
    }
    return values;
}
RANSWallEvolutionPatch::RANSWallEvolutionPatch(const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph,const LibUtilities::CommSharedPtr &comm,
    const RANSWallPatchParameters &parameters,const RANSWallNormalMap &mapping):
    m_dimension(graph->GetMeshDimension()),m_comm(comm),m_parameters(parameters),m_mapping(mapping)
{
    Require((m_dimension==2 || m_dimension==3) && session->GetComm()->GetSize()==1,
            "evolution requires a full serial native mesh on each collective rank");
    parameters.geometry.Validate();m_mapping.dimension=m_dimension;m_mapping.height=parameters.geometry.height;
    std::string solveDiagnostics,solveReplay,solveIterations;
    session->LoadSolverInfo("RANSWallSolveDiagnostics",solveDiagnostics,"Off");
    session->LoadSolverInfo("RANSWallSolveDiagnosticReplay",solveReplay,"Off");
    session->LoadSolverInfo("RANSWallSolveDiagnosticIterations",solveIterations,"12");
    session->LoadSolverInfo("RANSWallSolveDiagnosticFile",m_solveDiagnosticFile,"case.nonlinear.json");
    Require(solveDiagnostics=="Off" || solveDiagnostics=="Failure","unknown nonlinear diagnostic mode");
    Require(solveReplay=="Off" || solveReplay=="Fresh" || solveReplay=="AdmissibleFresh","unknown same-stage diagnostic replay mode");
    m_solveDiagnostics=solveDiagnostics=="Failure";m_solveDiagnosticReplay=solveReplay!="Off";
    m_solveDiagnosticAdmissible=solveReplay=="AdmissibleFresh";
    Require(!m_solveDiagnosticReplay || m_solveDiagnostics,"same-stage replay requires first-failure diagnostics");
    m_solveDiagnosticIterations=std::stoi(solveIterations);
    Require(m_solveDiagnosticIterations>0 && m_solveDiagnosticIterations<=20,"diagnostic same-stage retry is limited to 20 iterations");
    Require(!m_solveDiagnostics || m_comm->GetSize()==1,"detailed first-failure diagnostics currently require one MPI rank");
    session->LoadSolverInfo("RANSWallNonlinearPolicy",m_nonlinearPolicy,"Legacy");
    session->LoadSolverInfo("RANSWallFeasibleFile",m_feasibleTraceFile,"case.feasible.jsonl");
    session->LoadSolverInfo("RANSWallR3ProgressFile",m_r3ProgressFile,"case.r3-progress.jsonl");
    session->LoadSolverInfo("RANSWallActionPolicy",m_actionPolicy,"Legacy");
    session->LoadSolverInfo("RANSWallGlobalization",m_globalization,"FractionToBoundary");
    Require(m_globalization=="FractionToBoundary" ||
        (GroundSafeguard() && m_nonlinearPolicy=="ScalarFeasibleR3"),
        "CoupledR4 and JointRestorationR5 require ScalarFeasibleR3 and VerifiedR3Action");
    session->LoadSolverInfo("RANSWallDirectionAccuracy",m_directionAccuracy,"Legacy");
    Require(m_directionAccuracy=="Legacy" || (JointRestorationR5() && (DirectionReplayR6() || RefinedDirectionR6())),
        "R6 direction accuracy requires JointRestorationR5");
    session->LoadSolverInfo("RANSWallR7Policy",m_r7Policy,"Off");
    Require(m_r7Policy=="Off" || ((CompareR7() || InteriorProbeR7() || AdaptiveR7()) && RefinedDirectionR6()),
        "R7 policies require RefinedR6");
    if(CompareR7())
    {
        std::string displacement;
        session->LoadSolverInfo("RANSWallR7PreviousDisplacement",displacement,"0");
        std::size_t end=0;
        m_r7PreviousDisplacement=std::stod(displacement,&end);
        Require(end==displacement.size() && m_r7PreviousDisplacement>0. && m_r7PreviousDisplacement<=2.,
            "R7 comparison requires a measured previous displacement");
    }
    if(InteriorProbeR7() || AdaptiveR7())
    {
        std::string radius,interior;
        session->LoadSolverInfo("RANSWallR7Radius",radius,"0");
        session->LoadSolverInfo("RANSWallR7InteriorTarget",interior,"Off");
        std::size_t end=0;
        m_r7Radius=std::stod(radius,&end);
        Require(end==radius.size() && m_r7Radius>=.125 && m_r7Radius<=2. &&
            (interior=="On" || (AdaptiveR7() && interior=="Off")), "Invalid R7 radius or interior target");
        m_r7Interior=interior=="On";
    }
    std::string r8;
    session->LoadSolverInfo("RANSWallR8Experiment",r8,"Off");
    m_r8Experiment=r8=="On";m_r8Pilot=r8=="Pilot";
    Require(r8=="Off" || ((m_r8Experiment || m_r8Pilot) && RefinedDirectionR6() && m_r7Policy=="Off"),
        "R8 fixed-state experiment requires RefinedR6 with the R7 policy off");
    std::string r9;
    session->LoadSolverInfo("RANSWallR9Experiment",r9,"Off");
    m_r9Experiment=r9=="On";
    Require(r9=="Off" || (m_r9Experiment && RefinedDirectionR6() && m_r7Policy=="Off" && r8=="Off"),
        "R9 comparison requires only RefinedR6");
    std::string r10;
    session->LoadSolverInfo("RANSWallR10Experiment",r10,"Off");
    m_r10Experiment=r10=="On";
    Require(r10=="Off" || (m_r10Experiment && RefinedDirectionR6() && m_r7Policy=="Off" &&
        r8=="Off" && r9=="Off"),"R10 requires its isolated RefinedR6 comparison");
    std::string r11;
    session->LoadSolverInfo("RANSWallR11Experiment",r11,"Off");
    m_r11Experiment=r11=="On";
    Require(r11=="Off" || (m_r11Experiment && RefinedDirectionR6() && m_r7Policy=="Off" &&
        r8=="Off" && r9=="Off" && r10=="Off"),"R11 requires its isolated RefinedR6 comparison");
    session->LoadSolverInfo("RANSWallR12Experiment",m_r12Experiment,"Off");
    session->LoadParameter("RANSWallR12CutLimit",m_r12CutLimit,96);
    Require(m_r12Experiment=="Off" || ((m_r12Experiment=="Replay" || m_r12Experiment=="Complete" ||
        m_r12Experiment=="Verify" || m_r12Experiment=="Compare") && RefinedDirectionR6() &&
        m_r7Policy=="Off" && r8=="Off" && r9=="Off" && r10=="Off" && r11=="Off" &&
        m_r12CutLimit>=1 && m_r12CutLimit<=512),"R12 requires its isolated bounded RefinedR6 policy");
    Require(m_actionPolicy=="Legacy" || m_actionPolicy=="VerifiedR3Action","unknown wall action policy");
    Require((m_nonlinearPolicy=="ScalarFeasibleR3")== (m_actionPolicy=="VerifiedR3Action"),
        "VerifiedR3Action and ScalarFeasibleR3 must be selected together");
    Require(m_nonlinearPolicy=="Legacy" || m_nonlinearPolicy=="ScalarFeasibleR2" ||
        m_nonlinearPolicy=="ScalarFeasibleR3","unknown wall nonlinear policy");
    Require(m_nonlinearPolicy=="Legacy" || m_comm->GetSize()==1,
        "ScalarFeasibleR2 currently requires one MPI rank");
    session->LoadParameter("RANSWallOwnershipShift",m_ownershipShift,0);Require(m_ownershipShift>=0,"negative ownership shift");
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
    m_mapping.Validate(lower,upper,.1);
    std::string transition;session->LoadSolverInfo("RANSWallTransition",transition,"Fixed");
    Require(transition=="Fixed" || transition=="Upstream","RANSWallTransition must be Fixed or Upstream");
    m_transitionEnabled=transition=="Upstream";
    Require(m_mapping.PlanarConstantHeight() || parameters.chi==0 || m_transitionEnabled,
        "curved corrections require the prescribed upstream transition");
    for(const auto &point:std::vector<RANSWallCoordinate>{lower,upper,{lower[0],upper[1],0},{upper[0],lower[1],0}})
    {
        const double minH=m_mapping.At(point).height-std::abs(m_mapping.heightAmplitude)*2;
        Require(parameters.geometry.roughness<minH && parameters.geometry.sampleHeight<=minH && parameters.geometry.tailStart<=minH,
            "physical h/rK/rm do not fit every mapped ray");
    }
    if(m_mapping.kind==RANSWallNormalMap::Kind::CosineHill)
        for(int e=0;e<full->GetExpSize();++e)
        {
            double a=1e100,b=-1e100;const auto g=full->GetExp(e)->GetGeom();
            for(int j=0;j<g->GetNumVerts();++j){double x,y,z;g->GetVertex(j)->GetCoords(x,y,z);a=std::min(a,x);b=std::max(b,x);}
            for(double foot:{-m_mapping.lengthX,m_mapping.lengthX})Require(!(a<foot-1e-12 && b>foot+1e-12),"C1 hill foot requires an explicit native patch face");
        }
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
        m_mass.push_back(b>0 && b<m_dimension?m_mass[0]:PhysicalMass(b));
        m_constants.push_back(fields[b]->Project([](const RANSWallCoordinate&){return 1.;}));
    }
    m_volumes.resize(Gauge());
    for(int b=0;b<Gauge();++b)m_volumes[b]=std::inner_product(m_constants[b].begin(),m_constants[b].end(),Multiply(m_mass[b],m_constants[b]).begin(),0.);
    const double volume=m_volumes[0];
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
        ray.coordinate[2]=parameters.geometry.sampleHeight/m_mapping.At(ray.coordinate).scale;
        ray.weight*=m_mapping.At(ray.coordinate).baseArea;
        ray.sample=Point(K(),ray.coordinate,Same(ray.coordinate[2],parameters.geometry.height)?-1:0);
        auto H=ray.coordinate;H[2]=parameters.geometry.height;
        ray.outerV=Point(0,H,1);ray.outerP=Point(OuterPressure(),H,1);ray.outerK=Point(K(),H,1);
        m_rays.push_back(std::move(ray));
    }
    for(const auto &q:outerVelocity->Quadrature(order))if(Owner(q.geometry)==m_comm->GetRank())
        m_outer.push_back({Point(0,q.coordinate),Point(OuterPressure(),q.coordinate),
            Point(Epsilon(),q.coordinate),q.weight*(Mapped()?m_mapping.At(q.coordinate).determinant:1.)});
    m_solver=std::make_shared<SolverCache>();m_solver->free.resize(Gauge());m_solver->massFactors.resize(Gauge());
    InitialiseRobust(session);
    Vector seed(Size());for(int i=0;i<size[K()];++i)seed[offset[K()]+i]=parameters.model.Q()*.16*m_constants[K()][i];
    Refresh(seed);for(const auto &ray:m_rays)for(const auto &q:ray.points)m_bodyInner.push_back(q);
    if(m_transitionEnabled)InitialiseTransition(session);
    if(m_nonlinearPolicy=="ScalarFeasibleR2" || m_nonlinearPolicy=="ScalarFeasibleR3")
        Require(m_dimension==2 && !Mapped() && m_transitionEnabled && m_cappedDestruction,
            "ScalarFeasibleR2 prototype requires a flat 2D upstream patch and capped destruction");
}
std::shared_ptr<EquilibriumAWFProfiles> RANSWallEvolutionPatch::Profile(double k,const RANSWallCoordinate &chart,double weight) const
{
    if(weight<0)weight=m_transitionEnabled?0:m_parameters.chi;
    Require(k>0 && std::isfinite(k),"invalid physical sampled K");
    const auto metric=m_mapping.At(chart);auto geometry=m_parameters.geometry;geometry.height=metric.height;
    auto options=m_options;options.requestedWeight=weight;for(double &r:options.physicalBreakpoints)r*=metric.scale;
    const AWFNormalMetric normal{metric.kappa1,metric.kappa2};
    if(!m_transitionEnabled)return std::make_shared<EquilibriumAWFProfiles>(geometry,k,m_parameters.model,normal,options);
    // K and geometry determine the cutoff; velocity, pressure and epsilon
    // change only the analytic gate. Never reuse a previously active cutoff
    // on an inactive branch or solve one before the admissibility gates pass.
    if(m_profileCache.size()>4096)m_profileCache.clear();
    auto &entry=m_profileCache[{chart[0],chart[1],k}];
    const int active=weight>0?1:0;
    if(!entry[active]){options.requestedWeight=active;entry[active]=std::make_shared<EquilibriumAWFProfiles>(geometry,k,m_parameters.model,normal,options);}
    if(weight==0 || weight==1)return entry[active];
    return std::make_shared<EquilibriumAWFProfiles>(entry[1]->Reweighted(weight));
}
void RANSWallEvolutionPatch::Broadcast(Vector &value) const {m_comm->Bcast(value,0);}
double RANSWallEvolutionPatch::SampleK(const Vector &s,const RANSWallCoordinate &point) const
{
    auto p=point;p[2]=m_parameters.geometry.sampleHeight/m_mapping.At(point).scale;
    const double k=Point(K(),p,Same(p[2],m_parameters.geometry.height)?-1:0).Value(s,offset[K()]);
    Require(k>0 && std::isfinite(k),"invalid physical sampled K");return k;
}
std::array<double,3> RANSWallEvolutionPatch::CartesianForcing(const Vector &s,const RANSWallCoordinate &point,double time,
    const RANSWallEvolutionData &data,bool includeBody) const
{
    auto p=point;p[2]=m_parameters.geometry.height;const auto metric=m_mapping.At(p);
    const auto v=Point(0,p,1),pressure=Point(OuterPressure(),p,1),k=Point(K(),p,1);
    RANSTensor G{},P{};RANSVector u{},ut{},base{};
    for(int a=0;a<m_dimension;++a){const int aa=PhysicalAxis(a);u[aa]=v.Value(s,offset[a]);G[aa]=v.Gradient(s,offset[a]);}
    const auto &n=metric.normal;double w=0;for(int a=0;a<3;++a)w+=u[a]*n[a];
    for(int a=0;a<3;++a){ut[a]=u[a]-w*n[a];for(int b=0;b<3;++b)P[a][b]=(a==b?1.:0.)-n[a]*n[b];}
    RANSTensor gradient{};
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)
    {
        for(int c=0;c<3;++c)gradient[a][b]+=P[a][c]*G[c][b]-n[a]*metric.Br[c][b]*u[c];
        gradient[a][b]-=w*metric.Br[a][b];
    }
    double divergence=0;for(int a=0;a<3;++a)divergence+=gradient[a][a];
    const auto gp=pressure.Gradient(s,offset[OuterPressure()]),gk=k.Gradient(s,offset[K()]);
    for(int a=0;a<3;++a){base[a]=gp[a]-2*gk[a]/3+ut[a]*divergence;for(int b=0;b<3;++b)base[a]+=gradient[a][b]*ut[b];}
    if(includeBody)for(int a=0;a<m_dimension;++a)base[PhysicalAxis(a)]-=data.outerBody[a](p,time);
    RANSVector force{};for(int a=0;a<3;++a)for(int b=0;b<3;++b)force[a]+=P[a][b]*base[b];return force;
}
std::array<double,2> RANSWallEvolutionPatch::Forcing(const Vector &s,const RANSWallCoordinate &point,double time,
    const RANSWallEvolutionData &data,bool includeBody) const
{
    const auto f=CartesianForcing(s,point,time,data,includeBody);const auto g=m_mapping.At(point);
    std::array<double,2> result{};for(int a=0;a<3;++a){result[0]+=g.tangent1[a]*f[a];result[1]+=g.tangent2[a]*f[a];}return result;
}
void RANSWallEvolutionPatch::Refresh(const Vector &s,double time,const RANSWallEvolutionData *data) const
{
    for(auto &ray:m_rays)
    {
        const double k=ray.sample.Value(s,offset[K()]);
        auto profile=m_transitionReady && data?StateProfile(s,ray.coordinate,time,*data):Profile(k,ray.coordinate);
        const double weight=profile->Reference().weight;
        if(m_cappedDestruction)ray.destruction=WallNormalization(s,ray.coordinate,*profile);
        if(k==ray.cachedK && weight==ray.cachedWeight)continue;
        if(k==ray.cachedK && (weight>0)==(ray.cachedWeight>0))
        {ray.profile=std::move(profile);ray.cachedWeight=weight;continue;}
        std::vector<VolumePoint> points;
        for(const auto &q:profile->Quadrature())
        {
            auto p=ray.coordinate;p[2]=q.r/m_mapping.At(ray.coordinate).scale;auto v=Point(0,p);
            if(Owner(v.geometry)!=m_comm->GetRank())continue;
            const auto kp=Point(K(),p);Require(kp.ids==v.ids,"native k/velocity maps differ");
            points.push_back({std::move(v),Point(InnerPressure(),p),{},ray.weight*q.weight*q.jacobian});
        }
        ray.cachedK=k;ray.cachedWeight=weight;ray.profile=std::move(profile);ray.points=std::move(points);
    }
}
RANSWallEvolutionPatch::InnerDestruction RANSWallEvolutionPatch::DestructionValueR9(
    const ScalarDomainResult &domain,double time,const RANSWallEvolutionData &data,
    const std::map<int,double> &constraints) const
{
    Require(domain.admissible && m_cappedDestruction && m_comm->GetSize()==1 && !Mapped(),
        "R9 destruction requires the qualified serial flat scalar domain");
    const auto &s=domain.effectiveState;
    Refresh(s,time,&data);
    InnerDestruction out;
    out.value.resize(Size());
    // b, m, lambda and profile weight bounds. Count each cap and
    // limiter branch on the actual ordered moving quadrature.
    out.ranges={1.,0.,1.,0.,1.,0.,1.,0.,0.,0.};
    for(const auto &ray:m_rays)
    {
        const auto &n=ray.destruction;
        const auto &profile=*ray.profile;
        const std::array<double,4> ranges{n.b,n.m,n.lambda,profile.Reference().weight};
        for(int i=0;i<ranges.size();++i)
        {
            out.ranges[2*i]=std::min(out.ranges[2*i],ranges[i]);
            out.ranges[2*i+1]=std::max(out.ranges[2*i+1],ranges[i]);
        }
        out.ranges[8]+=n.capped;
        out.cap+=n.capped?'C':'P';
        const double weight=profile.Reference().weight;
        out.profile+=weight==0?'B':weight==1?'E':'M';
        for(const auto &q:ray.points)
        {
            const double k=Point(K(),q.v.coordinate).Value(s,offset[K()]);
            const bool limited=k/profile.K()<m_destruction.theta;
            out.limiter+=limited?'L':'S';
            out.ranges[9]+=limited;
            out.panels+=std::to_string(q.v.geometry)+",";
            // At zero production the native k source is -D.
            const double destruction=-WallSource(s,q.v.coordinate,0.,profile,&n).kSource;
            Add(out.value,offset[K()],q.v,destruction,{},q.weight);
        }
        out.limiter+='|';out.panels+='|';
    }
    for(const auto &[row,value]:constraints)out.value[row]=0.;
    return out;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::Project(const std::vector<RANSWallEvolutionFunction> &f,double time) const
{
    Require(f.size()==fields.size(),"incomplete evolution initial/reference functions");Vector s(Size());
    for(int b=0;b<Gauge();++b)
    {
        const auto values=PhysicalProject(b,f[b],time);
        std::copy(values.begin(),values.end(),s.begin()+offset[b]);
    }
    double mean=0;for(int side=0;side<2;++side)for(int i=0;i<m_pressureMean[side].size();++i)
        mean+=m_pressureMean[side][i]*s[offset[InnerPressure()+side]+i];
    for(int b=InnerPressure();b<=OuterPressure();++b)for(int i=0;i<size[b];++i)s[offset[b]+i]-=mean*m_constants[b][i];
    Broadcast(s);return s;
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::NeutralReference(double us,const std::array<double,2> &direction) const
{
    Require(m_mapping.PlanarConstantHeight(),"neutral reference requires a planar constant-height map");
    return FlatReference(us,direction);
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::FlatReference(double us,const std::array<double,2> &direction) const
{
    Require(us>0 && std::abs(direction[0]*direction[0]+direction[1]*direction[1]-1)<1e-12,"invalid reference friction velocity/direction");
    Require(m_dimension==3 || std::abs(direction[1])<1e-15,"2D reference direction must be streamwise");
    const double K0=m_parameters.model.Q()*us*us,H=m_parameters.geometry.height,Z=upper[2];
    auto options=m_options;if(m_transitionEnabled)options.requestedWeight=1;
    const auto profile=m_transitionEnabled?std::make_shared<EquilibriumAWFProfiles>(m_parameters.geometry,K0,m_parameters.model,AWFNormalMetric{},options):Profile(K0);
    const double UH=profile->Reference().matchingVelocity;
    const auto speed=[&](double r){return r<=H?profile->Original().Velocity(r,{UH,0},{0,0})[0]:UH+us/m_parameters.model.Kappa()*std::log(r/H);};
    const auto velocity=fields[0]->RadialProject(speed,{{0,0,0},{H,0,UH},
        {H,-1,us*us/profile->Original().TotalViscosity(H)},{H,1,m_parameters.model.Q()*profile->Dissipation(H)/K0},
        {Z,0,speed(Z)},{Z,-1,us/(m_parameters.model.Kappa()*Z)}});
    Vector s(Size());
    auto flat=m_mapping;if(!m_mapping.PlanarConstantHeight()){flat.kind=RANSWallNormalMap::Kind::Flat;flat.slopeX=flat.slopeY=0;}
    const auto frame=flat.At({});RANSVector flow{};
    for(int a=0;a<3;++a)flow[a]=direction[0]*frame.tangent1[a]+direction[1]*frame.tangent2[a];
    for(int a=0;a<m_dimension;++a)for(int i=0;i<size[a];++i)s[offset[a]+i]=flow[PhysicalAxis(a)]*velocity[i];
    for(int i=0;i<size[K()];++i)s[offset[K()]+i]=K0*m_constants[K()][i];
    const double A=us*us*us/m_parameters.model.Kappa();
    const auto ep=fields[Epsilon()]->RadialProject([&](double r){return A/r;},
        {{H,0,profile->Dissipation(H)},{Z,0,A/Z},{Z,-1,-A/(Z*Z)}});
    std::copy(ep.begin(),ep.end(),s.begin()+offset[Epsilon()]);
    RANSTensor G{};for(int a=0;a<3;++a)for(int b=0;b<3;++b)G[a][b]=flow[a]*frame.normal[b]*us/(m_parameters.model.Kappa()*H);
    const auto closure=EvaluateRANSClosure(m_parameters.model.Model(),K0,A/H,G,H,m_parameters.model.Parameters(),0);
    double jump=-2*K0/3;for(int a=0;a<3;++a)for(int b=0;b<3;++b)jump+=frame.normal[a]*closure.reynoldsStress[a][b]*frame.normal[b];
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
    if(Mapped())
    {
        if(m_boundaryFlux.empty())
        {
            m_boundaryFlux.assign(DifferentialSize(),0.);
            for(int axis=0;axis<3;++axis)if(m_dimension==3 || axis!=1)for(int side:{-1,1})
                for(const auto &face:Surface(0,axis,side<0?lower[axis]:upper[axis],side))
                    for(int a=0;a<m_dimension;++a)for(int j=0;j<face.test.ids.size();++j)
                        m_boundaryFlux[offset[a]+face.test.ids[j]]+=face.weight*face.area[PhysicalAxis(a)]*face.test.phi[j];
        }
        if(m_boundaryDirection.empty())
        {
            // Physical surface-L2 projection subject to the one exact net-flux
            // constraint. Fix edge/corner and side traces; distribute the
            // correction over all Cartesian top-face interior modes.
            const auto excluded=fields[0]->Boundary([](const RANSWallCoordinate&){return 0.;},
                [&](int axis,int side,const RANSWallGlobalField::Element &e){return axis!=2?Same(side?e.upper[axis]:e.lower[axis],side?upper[axis]:lower[axis]):!side && Same(e.lower[2],0);});
            const auto top=fields[0]->Boundary([](const RANSWallCoordinate&){return 0.;},
                [&](int axis,int side,const RANSWallGlobalField::Element &e){return axis==2 && side && Same(e.upper[2],upper[2]);});
            std::vector<int> ids;std::map<int,int> index;
            for(auto [id,value]:top)if(!excluded.count(id)){(void)value;index[id]=ids.size();ids.push_back(id);}
            Require(!ids.empty(),"no top interior modes for conservative boundary projection");
            WallMatrix mass(ids.size(),Vector(ids.size()));double area=0;
            for(const auto &face:Surface(0,2,upper[2],1))
            {
                const double measure=face.weight*std::hypot(face.area[0],face.area[1],face.area[2]);area+=measure;
                for(int i=0;i<face.test.ids.size();++i)if(index.count(face.test.ids[i]))
                    for(int j=0;j<face.test.ids.size();++j)if(index.count(face.test.ids[j]))
                        mass[index.at(face.test.ids[i])][index.at(face.test.ids[j])]+=measure*face.test.phi[i]*face.test.phi[j];
            }
            Factor factor(mass);m_boundaryDirection.assign(DifferentialSize(),0.);double denominator=0;
            for(int a=0;a<m_dimension;++a)
            {
                Vector load;for(int id:ids)load.push_back(m_boundaryFlux[offset[a]+id]);
                const auto direction=factor.Apply(load);
                for(int j=0;j<ids.size();++j){m_boundaryDirection[offset[a]+ids[j]]=direction[j];denominator+=load[j]*direction[j];}
            }
            Require(denominator>0 && area>0,"singular conservative boundary flux constraint");
            for(double &v:m_boundaryDirection)v/=denominator;
            m_boundaryCost=1/std::sqrt(denominator*area);
        }
        double flux=0;for(auto [row,value]:result)if(row<offset[K()])flux+=m_boundaryFlux[row]*value;
        m_boundaryCorrection=std::abs(flux)*m_boundaryCost;
        for(auto &[row,value]:result)if(row<offset[K()])value-=flux*m_boundaryDirection[row];
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
    if(m_transitionEnabled)
    {
        auto fixed=ExternalConstraints(time,data);auto s=ExternalLift(state,time,data);
        const auto trace=DissipationTrace(s,time,data);
        for(auto [id,value]:trace[0])fixed[offset[Epsilon()]+id]=value+trace[1].at(id);
        return fixed;
    }
    auto fixed=ExternalConstraints(time,data);auto s=state;for(auto [row,value]:fixed)s[row]=value;
    const auto ep=fields[Epsilon()]->Boundary([&](const RANSWallCoordinate &p){const auto profile=Profile(SampleK(s,p),p);return WallDissipation(s,p,*profile);},
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
{
    auto s=state;
    // The independent outer epsilon sample contains dependent trace modes.
    // Solve that trace constraint at initialization; stage Newton retains its
    // columns explicitly and never substitutes a lagged sensor value.
    for(int iteration=0;iteration<(m_transitionEnabled?100:1);++iteration)
    {
        double change=0;for(auto [row,value]:Constraints(s,time,data)){change=std::max(change,std::abs(s[row]-value));s[row]=value;}
        if(!m_transitionEnabled || change<2e-14)return s;
    }
    throw std::runtime_error("dependent transition trace lifting did not converge");
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::Residual(const Vector &state,const Vector &history,double alpha,double time,
    const RANSWallEvolutionData &data,bool anchored,bool imposeConstraints,int region) const
{
    if(m_solveDiagnostics){m_solveResidualException.clear();m_solveInvalidContext.clear();}
    Require(state.size()==static_cast<size_t>(Size()) && history.size()==state.size() && alpha>=0,"invalid evolution stage state/history");
    const auto constraints=Constraints(state,time,data);auto s=m_transitionEnabled?ExternalLift(state,time,data):state;
    if(!m_transitionEnabled)for(auto [row,value]:constraints)s[row]=value;
    Vector r(Size());int bad=0;
    try
    {
        Refresh(s,time,&data);
        auto evaluate=[&](const VolumePoint &q,bool inner,const EquilibriumAWFProfiles *profile,const std::array<double,3> &force,
            const RANSCappedDestructionValues *normalization)
        {
            RANSTensor G{},stress{};std::array<double,3> u{};
            for(int a=0;a<m_dimension;++a){u[PhysicalAxis(a)]=q.v.Value(s,offset[a]);G[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);}
            const double k=q.v.Value(s,offset[K()]),pi=q.p.Value(s,offset[inner?InnerPressure():OuterPressure()]);
            const auto gk=q.v.Gradient(s,offset[K()]);Require(k>=0 && std::isfinite(k),"negative/nonfinite transported k");
            double dk,sk,de=0,se=0,epsilon=0;
            if(inner)
            {
                const double nt=profile->Original().EddyViscosity(Distance(q.v.coordinate));
                const auto innerStress=InnerStress(k,nt,G);
                stress=RANSWallInnerMeanStress(innerStress,m_parameters.geometry.molecularNu,pi);
                dk=m_parameters.geometry.molecularNu+nt/m_parameters.model.Parameters().sigmaK;
                sk=WallSource(s,q.v.coordinate,innerStress.production,*profile,normalization).kSource+data.innerKSource(q.v.coordinate,time);
            }
            else
            {
                epsilon=q.e.Value(s,offset[Epsilon()]);
                const auto c=OuterClosure(k,epsilon,G,q.v.coordinate);
                dk=c.kDiffusion;de=c.epsilonDiffusion;sk=c.kSource+data.outerKSource(q.v.coordinate,time);
                se=c.epsilonSource+data.outerEpsilonSource(q.v.coordinate,time);
                for(int a=0;a<3;++a)for(int b=0;b<3;++b)stress[a][b]=-c.reynoldsStress[a][b]+(a==b?2*k/3-pi:0);
            }
            for(int a=0;a<m_dimension;++a)
            {
                const int aa=PhysicalAxis(a);std::array<double,3> flux{};
                for(int b=0;b<3;++b)flux[b]=stress[aa][b]-u[aa]*u[b];
                const double rough=inner && Distance(q.v.coordinate)<=m_parameters.geometry.roughness?force[PhysicalAxis(a)]:0;
                Add(r,offset[a],q.v,-rough,flux,q.weight);
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
        if(region!=1)for(const auto &ray:m_rays)
        {
            const auto force=CartesianForcing(s,ray.coordinate,time,data,false);
            for(const auto &q:ray.points)evaluate(q,true,ray.profile.get(),force,m_cappedDestruction?&ray.destruction:nullptr);
        }
        if(region!=0)for(const auto &q:m_outer)evaluate(q,false,nullptr,{},nullptr);
        const auto body=OwnedBodyLoad(time,data,region);
        for(int i=0;i<Size();++i)r[i]-=body[i];
    }
    catch(const std::exception &error)
    {
        bad=1;
        if(m_solveDiagnostics)m_solveResidualException=RANSWallSolveJson::Object({
            {"type",RANSWallSolveJson::Quote(RANSWallSolveJson::ExceptionType(error))},
            {"message",RANSWallSolveJson::Quote(error.what())}});
    }
    m_comm->AllReduce(bad,LibUtilities::ReduceMax);Require(!bad,"inadmissible evolution volume/sampler state on a collective rank");
    m_comm->AllReduce(r,LibUtilities::ReduceSum);
    if(alpha)for(int b=0;b<InnerPressure();++b)
    {
        Vector delta(size[b]);for(int i=0;i<size[b];++i)delta[i]=s[offset[b]+i]-history[offset[b]+i];
        Vector load;
        if(region<0)load=Multiply(m_mass[b],delta);
        else
        {
            const auto mass=b==Epsilon()?(region==1?m_mass[b]:WallMatrix(size[b],Vector(size[b]))):PhysicalMass(b,region);
            load=Multiply(mass,delta);
        }for(int i=0;i<size[b];++i)r[offset[b]+i]+=alpha*load[i];
    }
    for(int side=0;side<2;++side)if(region<0 || side==region)for(int i=0;i<m_pressureMean[side].size();++i)
    {
        const int row=offset[InnerPressure()+side]+i;const double w=m_pressureMean[side][i];
        r[row]+=s[offset[Gauge()]]*w;r[offset[Gauge()]]+=s[row]*w;
    }
    if(anchored && m_eta && region<0)for(int i=0;i<DifferentialSize();++i)r[i]-=m_eta*m_defect[i];
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
        norms.push_back(std::sqrt(std::max(0.L,square))/(rate*std::sqrt(m_volumes[b])));
    }
    norms.push_back(std::abs(r[offset[Gauge()]])/(us*us));double bc=0;
    for(auto [row,value]:constraints)
    {
        (void)value;int b=0;while(b+1<offset.size() && row>=offset[b+1])++b;
        const double scale=b<m_dimension?us:(b==K()?K0:E0);bc=std::max(bc,std::abs(r[row])/scale);
    }
    norms.push_back(bc);return norms;
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::DualProductsR10(
    const Vector &a,const Vector &b,const std::map<int,double> &constraints) const
{
    Norms(a,constraints,m_robustUstar); // Populate the identical free consistent-mass factors.
    Vector products;
    const double us=m_robustUstar,H=m_parameters.geometry.height;
    const double K0=m_parameters.model.Q()*us*us,E0=us*us*us/(m_parameters.model.Kappa()*H);
    for(int block=0;block<Gauge();++block)
    {
        const auto &free=m_solver->free[block];
        if(free.empty()){products.push_back(0.);continue;}
        Vector load;
        for(int id:free)load.push_back(b[offset[block]+id]);
        const auto inverse=m_solver->massFactors[block]->Apply(load);
        long double sum=0;
        for(int i=0;i<free.size();++i)sum+=(long double)a[offset[block]+free[i]]*inverse[i];
        const double scale=block<m_dimension?us:block==K()?K0:block==Epsilon()?E0:1.;
        const double rate=block>=InnerPressure()?us/H:scale*us/H;
        products.push_back(sum/(rate*rate*m_volumes[block]));
    }
    return products;
}

std::vector<RANSWallEvolutionPatch::Vector> RANSWallEvolutionPatch::MeritRowsR10(
    const Vector &residual,const WallMatrix &jac,const std::map<int,double> &constraints) const
{
    const double us=m_robustUstar,H=m_parameters.geometry.height;
    const double K0=m_parameters.model.Q()*us*us,E0=us*us*us/(m_parameters.model.Kappa()*H);
    const auto norms=Norms(residual,constraints,us);
    const double maximum=Maximum(norms);
    Require(maximum>0.,"R10 slope rows require a nonzero native merit");
    const double tie=64*std::numeric_limits<double>::epsilon()*maximum;
    std::vector<Vector> rows;
    for(int block=0;block<Gauge();++block)if(norms[block]>=maximum-tie)
    {
        Vector load;
        for(int id:m_solver->free[block])load.push_back(residual[offset[block]+id]);
        const auto dual=m_solver->massFactors[block]->Apply(load);
        const double scale=block<m_dimension?us:block==K()?K0:block==Epsilon()?E0:1.;
        const double rate=block>=InnerPressure()?us/H:scale*us/H;
        const double denominator=rate*rate*m_volumes[block]*norms[block];
        Vector row(Size());
        for(int col=0;col<Size();++col)
        {
            long double value=0.;
            for(int i=0;i<dual.size();++i)value+=(long double)dual[i]*jac[offset[block]+m_solver->free[block][i]][col];
            row[col]=value/denominator;
        }
        rows.push_back(std::move(row));
    }
    if(norms[Gauge()]>=maximum-tie)
    {
        auto row=jac[offset[Gauge()]];
        for(double &v:row)v*=std::copysign(1.,residual[offset[Gauge()]])/(us*us);
        rows.push_back(std::move(row));
    }
    if(norms.back()>=maximum-tie)for(const auto &[id,value]:constraints)
    {
        (void)value;
        int block=0;while(block+1<offset.size() && id>=offset[block+1])++block;
        const double scale=block<m_dimension?us:block==K()?K0:E0;
        if(std::abs(residual[id])/scale<maximum-tie)continue;
        auto row=jac[id];
        for(double &v:row)v*=std::copysign(1.,residual[id])/scale;
        rows.push_back(std::move(row));
    }
    Require(!rows.empty(),"R10 native merit has no active block");
    return rows;
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
        norms.push_back(std::sqrt(std::max(0.L,sum))/(scale*std::sqrt(m_volumes[b])));
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
    if(m_nonlinearPolicy=="ScalarFeasibleR3")
    {
        Require(!condensed,"ScalarFeasibleR3 retains the full coupled dense backend");
        return SolveScalarFeasibleR3(std::move(s),history,alpha,time,data,iterations,tolerance,maximumIterations);
    }
    if(m_nonlinearPolicy=="ScalarFeasibleR2")
        return SolveScalarFeasible(std::move(s),history,alpha,time,data,condensed,iterations,tolerance,maximumIterations);
    namespace J=RANSWallSolveJson;
    const bool recording=m_solveDiagnostics && !m_solveDiagnosticWritten;
    std::vector<std::string> records;
    std::string assemblyFailure="null";
    const int entryBuilds=m_solver->builds,entryAge=m_solver->uses;
    Vector scales(Size(),.16);
    const double us=.4,K0=m_parameters.model.Q()*us*us;
    const double E0=us*us*us/(m_parameters.model.Kappa()*m_parameters.geometry.height);
    for(int b=0;b<Gauge();++b)for(int i=0;i<size[b];++i)
        scales[offset[b]+i]=b<m_dimension?us:b==K()?K0:b==Epsilon()?E0:us*us;
    // The gauge multiplier enters continuity and has inverse-time units.
    scales[offset[Gauge()]]=us/m_parameters.geometry.height;
    auto constraintJSON=[&](const std::map<int,double> &constraints)
    {
        std::vector<std::string> entries;
        for(auto [row,value]:constraints)entries.push_back(J::Numbers(Vector{double(row),value}));
        return J::List(entries);
    };
    auto physical=[&](const Vector &state)
    {
        try
        {
            auto actual=m_transitionEnabled?ExternalLift(state,time,data):state;
            if(!m_transitionEnabled)for(auto [row,value]:Constraints(state,time,data))actual[row]=value;
            Refresh(actual,time,&data);
            double ki=1e300,ko=1e300,eo=1e300,ci=1e300,co=1e300,kr=1e300;
            double bmin=1,lmin=1,active=0,capped=0,limiterGap=1e300,capGap=1e300;
            double gradientError=0,kRatioMin=1e300,kRatioMax=-1e300,epsRatioMin=1e300,epsRatioMax=-1e300;
            double invariantGap=1e300;
            Vector worst;RANSTensor worstG{},worstS{};
            std::string branches,limiterBranches,capBranches;
            for(const auto &ray:m_rays)
            {
                if(m_cappedDestruction)
                {
                    bmin=std::min(bmin,ray.destruction.b);lmin=std::min(lmin,ray.destruction.lambda);
                    capGap=std::min(capGap,std::abs(ray.destruction.b-m_destruction.delta));
                    capBranches.push_back(ray.destruction.b<m_destruction.delta?'C':'P');
                }
                for(const auto &q:ray.points)
                {
                    const double k=q.v.Value(actual,offset[K()]),K=ray.profile->K();
                    const double nt=ray.profile->Original().EddyViscosity(Distance(q.v.coordinate));
                    RANSTensor G{},S{};for(int a=0;a<m_dimension;++a)G[PhysicalAxis(a)]=q.v.Gradient(actual,offset[a]);
                    for(int a=0;a<3;++a)for(int b=0;b<3;++b)S[a][b]=G[a][b]+G[b][a];
                    const auto inner=InnerStress(k,nt,G);
                    const double largest=J::Eigenvalues(S)[2],cov=J::Eigenvalues(inner.reynoldsStress)[0];
                    if(cov<ci)
                    {
                        ci=cov;worstG=G;worstS=S;
                        worst={q.v.coordinate[0],q.v.coordinate[1],q.v.coordinate[2],k,K,k/K,nt,largest,cov,k>0?3*nt*largest/(2*k):1e300};
                    }
                    ki=std::min(ki,k);kr=std::min(kr,k/K);
                    const double gap=k/K-m_destruction.theta;
                    limiterGap=std::min(limiterGap,std::abs(gap));limiterBranches.push_back(gap<0?'L':'P');
                    if(m_cappedDestruction && gap<0)active+=q.weight;
                    if(m_cappedDestruction && ray.destruction.capped)capped+=q.weight;
                }
            }
            for(const auto &q:m_outer)
            {
                const double k=q.v.Value(actual,offset[K()]),eps=q.e.Value(actual,offset[Epsilon()]);
                RANSTensor G{};for(int a=0;a<m_dimension;++a)G[PhysicalAxis(a)]=q.v.Gradient(actual,offset[a]);
                const auto c=OuterClosure(k,eps,G,q.v.coordinate);
                ko=std::min(ko,k);eo=std::min(eo,eps);co=std::min(co,J::Eigenvalues(c.reynoldsStress)[0]);
                const double gap=c.strainInvariant-c.rotationInvariant;
                invariantGap=std::min(invariantGap,std::abs(gap));branches.push_back(gap>0?'S':gap<0?'R':'T');
                const double g0=us/(m_parameters.model.Kappa()*Distance(q.v.coordinate));
                for(int a=0;a<3;++a)for(int b=0;b<3;++b)
                    gradientError=std::max(gradientError,std::abs(G[a][b]/g0-m_referenceFlow[a]*m_referenceNormal[b]));
                const double e0=us*us*us/(m_parameters.model.Kappa()*Distance(q.v.coordinate));
                kRatioMin=std::min(kRatioMin,k/K0);kRatioMax=std::max(kRatioMax,k/K0);
                epsRatioMin=std::min(epsRatioMin,eps/e0);epsRatioMax=std::max(epsRatioMax,eps/e0);
            }
            double nativeKi=1e300,nativeKo=1e300,nativeE=1e300;Vector nativeWorst;
            for(int e=0;e<fields[K()]->Elements().size();++e)
            {
                const auto &el=fields[K()]->Elements()[e];const auto exp=fields[K()]->Field()->GetExp(e);
                const int nx=exp->GetNumPoints(0),ny=m_dimension==3?exp->GetNumPoints(1):1;
                const int nz=exp->GetNumPoints(m_dimension-1);const bool inner=el.upper[2]<=m_parameters.geometry.height;
                for(int iz=0;iz<nz;++iz)for(int iy=0;iy<ny;++iy)for(int ix=0;ix<nx;++ix)
                {
                    RANSWallCoordinate p{};const std::array<int,3> indices{ix,iy,iz};
                    for(int a=0;a<3;++a)if(m_dimension==3 || a!=1)
                    {
                        const int axis=m_dimension==2 && a==2?1:a;
                        const double z=exp->GetBasis(axis)->GetZ()[indices[a]];
                        p[a]=el.lower[a]+(z+1)*(el.upper[a]-el.lower[a])/2;
                    }
                    const auto q=fields[K()]->ElementPoint(e,p);const double k=q.Value(actual,offset[K()]);
                    if(inner)
                    {if(k<nativeKi){nativeKi=k;nativeWorst={p[0],p[1],p[2],k};}}
                    else
                    {
                        nativeKo=std::min(nativeKo,k);
                        nativeE=std::min(nativeE,Point(Epsilon(),p,std::abs(p[2]-m_parameters.geometry.height)<1e-12?1:0).Value(actual,offset[Epsilon()]));
                    }
                }
            }
            return J::Object({{"status",J::Quote("evaluated_source_quadrature")},
                {"inner_stress_model",J::Quote(InnerStressName())},
                {"worst_point_last_ratio_scope",J::Quote("original nuA*S ratio, not the saturated stress ratio")},
                {"minimum_inner_k",J::Number(ki)},{"minimum_outer_k",J::Number(ko)},
                {"minimum_outer_epsilon",J::Number(eo)},{"minimum_inner_covariance",J::Number(ci)},
                {"minimum_native_inner_k",J::Number(nativeKi)},{"minimum_native_outer_k",J::Number(nativeKo)},
                {"minimum_native_outer_epsilon",J::Number(nativeE)},{"native_inner_k_point",J::Numbers(nativeWorst)},
                {"native_finite_scalar_admissible",J::Boolean(nativeKi>=0 && nativeKo>0 && nativeE>0)},
                {"minimum_outer_covariance",J::Number(co)},{"minimum_inner_k_over_sampled_K",J::Number(kr)},
                {"worst_inner_point",J::Numbers(worst)},{"worst_inner_gradient",J::Rows(worstG)},
                {"worst_inner_strain",J::Rows(worstS)},{"minimum_b",J::Number(bmin)},
                {"minimum_lambda",J::Number(lmin)},{"limiter_active_volume",J::Number(active)},
                {"normalization_capped_volume",J::Number(capped)},
                {"outer_normalized_gradient_entry_error",J::Number(gradientError)},
                {"outer_k_over_K0",J::Numbers(Vector{kRatioMin,kRatioMax})},
                {"outer_epsilon_over_epsilon0",J::Numbers(Vector{epsRatioMin,epsRatioMax})},
                {"minimum_invariant_branch_gap",J::Number(invariantGap)},
                {"minimum_limiter_knot_gap",J::Number(limiterGap)},
                {"minimum_normalization_knot_gap",J::Number(capGap)},
                {"outer_branch_pattern",J::Quote(branches)},
                {"limiter_branch_pattern",J::Quote(limiterBranches)},
                {"cap_branch_pattern",J::Quote(capBranches)}});
        }
        catch(const std::exception &error)
        {return J::Object({{"status",J::Quote("diagnostic_evaluation_failed")},
            {"type",J::Quote(J::ExceptionType(error))},{"message",J::Quote(error.what())}});}
    };
    auto assemble=[&](const Vector &state,const std::map<int,double> &constraints)
    {
        WallMatrix jac(Size(),Vector(Size()));
        for(int j=0;j<Size();++j)
        {
            if(constraints.count(j) && !(m_transitionEnabled && j>=offset[Epsilon()] && j<offset[InnerPressure()] && m_traceIds.count(j-offset[Epsilon()])))
            {jac[j][j]=1.;continue;}
            auto plus=state,minus=state;const double h=2e-6*std::max(1.,std::abs(state[j]));plus[j]+=h;minus[j]-=h;
            std::string side="plus";
            try
            {
                const auto a=Residual(plus,history,alpha,time,data);side="minus";
                const auto b=Residual(minus,history,alpha,time,data);
                for(int i=0;i<Size();++i)jac[i][j]=(a[i]-b[i])/(2*h);
            }
            catch(const std::exception &error)
            {
                if(recording)assemblyFailure=J::Object({{"column",J::Number(j)},
                    {"side",J::Quote(side)},{"step",J::Number(h)},{"type",J::Quote(J::ExceptionType(error))},
                    {"message",J::Quote(error.what())},
                    {"invalid_context",m_solveInvalidContext.empty()?"null":m_solveInvalidContext},
                    {"volume_exception",m_solveResidualException.empty()?"null":m_solveResidualException}});
                throw;
            }
        }
        return jac;
    };
    std::vector<std::string> diagnosticColumns;
    auto admissibleAssemble=[&](const Vector &state,const std::map<int,double> &constraints)
    {
        // Diagnostic only: same residual, dimensionless column steps and
        // admissible central/one-sided differences. Never used by Solve's
        // production refresh policy or accepted time integration.
        const auto base=Residual(state,history,alpha,time,data);
        WallMatrix jac(Size(),Vector(Size()));diagnosticColumns.clear();
        for(int j=0;j<Size();++j)
        {
            if(constraints.count(j) && !(m_transitionEnabled && j>=offset[Epsilon()] && j<offset[InnerPressure()] && m_traceIds.count(j-offset[Epsilon()])))
            {jac[j][j]=1.;diagnosticColumns.push_back(J::Object({{"column",J::Number(j)},{"stencil",J::Quote("lifted_constraint_identity")}}));continue;}
            Vector a,b,one;double hUsed=0;int side=0;bool central=false;
            std::vector<std::string> attempts;
            const double first=2e-6*std::max(scales[j],std::abs(state[j]));
            for(double shrink:{1.,.25,.0625,.015625,.00390625})
            {
                const double h=first*shrink;Vector plus=state,minus=state;plus[j]+=h;minus[j]-=h;
                Vector p,m;bool pValid=false,mValid=false;std::string pe="null",me="null";
                try{p=Residual(plus,history,alpha,time,data);pValid=true;}
                catch(const std::exception &error)
                {pe=J::Object({{"type",J::Quote(J::ExceptionType(error))},{"message",J::Quote(error.what())},
                    {"invalid_context",m_solveInvalidContext.empty()?"null":m_solveInvalidContext}});}
                try{m=Residual(minus,history,alpha,time,data);mValid=true;}
                catch(const std::exception &error)
                {me=J::Object({{"type",J::Quote(J::ExceptionType(error))},{"message",J::Quote(error.what())},
                    {"invalid_context",m_solveInvalidContext.empty()?"null":m_solveInvalidContext}});}
                attempts.push_back(J::Object({{"step",J::Number(h)},{"plus_valid",J::Boolean(pValid)},
                    {"minus_valid",J::Boolean(mValid)},{"plus_exception",pe},{"minus_exception",me}}));
                if(pValid && mValid){a=std::move(p);b=std::move(m);hUsed=h;central=true;break;}
                if(pValid || mValid){one=pValid?std::move(p):std::move(m);side=pValid?1:-1;hUsed=h;}
            }
            if(!central && !side)
            {
                diagnosticColumns.push_back(J::Object({{"column",J::Number(j)},{"stencil",J::Quote("no_admissible_coordinate_probe")},{"attempts",J::List(attempts)}}));
                throw std::runtime_error("no admissible diagnostic finite-difference side for column "+std::to_string(j));
            }
            std::string stencil=central?"central":"one_sided_first_order";
            if(central)for(int i=0;i<Size();++i)jac[i][j]=(a[i]-b[i])/(2*hUsed);
            else
            {
                bool secondValid=false;Vector secondState=state,second;
                secondState[j]+=2*side*hUsed;
                try{second=Residual(secondState,history,alpha,time,data);secondValid=true;}
                catch(const std::exception &){}
                for(int i=0;i<Size();++i)jac[i][j]=secondValid?
                    side*(-3*base[i]+4*one[i]-second[i])/(2*hUsed):side*(one[i]-base[i])/hUsed;
                if(secondValid)stencil="one_sided_second_order";
            }
            diagnosticColumns.push_back(J::Object({{"column",J::Number(j)},{"step",J::Number(hUsed)},
                {"side",J::Number(central?0:side)},{"stencil",J::Quote(stencil)},{"attempts",J::List(attempts)}}));
        }
        return jac;
    };
    struct Search { bool accepted=false;Vector state;std::string observations; };
    auto search=[&](const Vector &state,const Vector &step,double error,bool detail)
    {
        Search result;std::vector<std::string> trials;
        for(double fraction=1;fraction>=1./1048576;fraction/=2)
        {
            auto trial=state;for(int i=0;i<Size();++i)trial[i]+=fraction*step[i];
            try
            {
                const auto next=Residual(trial,history,alpha,time,data);
                const auto constraints=Constraints(trial,time,data);const auto norms=Norms(next,constraints);
                const double nextError=Maximum(norms);const bool accepted=nextError<error;
                if(detail)trials.push_back(J::Object({{"fraction",J::Number(fraction)},
                    {"evaluation_succeeded","true"},{"classification",J::Quote(accepted?"accepted_descent":"valid_non_descent")},
                    {"residual",J::Numbers(next)},{"scaled_block_norms",J::Numbers(norms)},
                    {"maximum_scaled_residual",J::Number(nextError)},
                    {"constraints",constraintJSON(constraints)},{"physical",physical(trial)}}));
                if(accepted){result.accepted=true;result.state=std::move(trial);break;}
            }
            catch(const std::exception &error)
            {
                if(detail)trials.push_back(J::Object({{"fraction",J::Number(fraction)},
                    {"evaluation_succeeded","false"},{"classification",J::Quote("invalid_residual_or_norm_evaluation")},
                    {"type",J::Quote(J::ExceptionType(error))},{"message",J::Quote(error.what())},
                    {"invalid_context",m_solveInvalidContext.empty()?"null":m_solveInvalidContext},
                    {"volume_exception",m_solveResidualException.empty()?"null":m_solveResidualException}}));
            }
        }
        result.observations=J::List(trials);return result;
    };
    auto stepRecord=[&](const Vector &state,const Vector &r,const WallMatrix &jac,const Vector &step,
                        const std::map<int,double> &constraints)
    {
        auto linear=Multiply(jac,step);double scaledStep=0;
        for(int i=0;i<Size();++i){linear[i]+=r[i];scaledStep=std::max(scaledStep,std::abs(step[i])/scales[i]);}
        return J::Object({{"linear_residual",J::Numbers(linear)},
            {"linear_scaled_block_norms",J::Numbers(Norms(linear,constraints))},
            {"maximum_dimensionless_step",J::Number(scaledStep)},
            {"state",J::Numbers(state)},{"step",J::Numbers(step)}});
    };
    bool refreshed=false;double previous=1e100;
    try
    {
    for(iterations=0;iterations<=maximumIterations;++iterations)
    {
        const auto r=Residual(s,history,alpha,time,data);const auto constraints=Constraints(s,time,data);
        const auto norms=Norms(r,constraints);const double error=Maximum(norms);
        if(error<tolerance)
        {
            if(Mapped())Require(std::abs(s[offset[Gauge()]])/(.4/m_parameters.geometry.height*m_volumes[0])<1e-9,
                "incompatible boundary flux: gauge reaction would create a mass source");
            return s;
        }
        Require(iterations<maximumIterations,"evolution stage nonlinear iteration limit, scaled residual="+std::to_string(error));
        std::string reason="cached";
        if(!m_solver->factor || !Same(alpha,m_solver->alpha) || (error>previous*.8 && !refreshed))
        {
            reason=!m_solver->factor?"no_factor":!Same(alpha,m_solver->alpha)?"changed_alpha":"slow_progress_first_refresh";
            if(m_comm->GetRank()==0)std::cout<<"EP-AWF stage Jacobian at time "<<time<<", alpha "<<alpha<<", scaled residual "<<error<<std::endl;
            auto jac=assemble(s,constraints);
            m_solver->jacobian=std::move(jac);m_solver->alpha=alpha;
            ++m_solver->builds;m_solver->uses=0;m_solver->buildTime=time;m_solver->buildIteration=iterations;
            if(recording){m_solver->buildState=s;m_solver->buildHistory=history;}
            int bad=0;
            if(m_comm->GetRank()==0)
            {
                try{m_solver->factor=std::make_shared<Factor>(m_solver->jacobian);}
                catch(const std::exception &error)
                {
                    bad=1;
                    if(recording)assemblyFailure=J::Object({{"classification",J::Quote("linear_factorization_failure")},
                        {"type",J::Quote(J::ExceptionType(error))},{"message",J::Quote(error.what())}});
                }
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
            catch(const std::exception &error)
            {
                bad=1;
                if(recording)assemblyFailure=J::Object({{"classification",J::Quote("linear_solve_failure")},
                    {"type",J::Quote(J::ExceptionType(error))},{"message",J::Quote(error.what())}});
            }
        }
        m_comm->Bcast(bad,0);Require(!bad,"collective Newton/Schur solve failed");
        Broadcast(step);
        auto searched=search(s,step,error,recording);const bool accepted=searched.accepted;
        if(recording)records.push_back(J::Object({{"iteration",J::Number(iterations)},
            {"jacobian_build",J::Number(m_solver->builds)},{"jacobian_age",J::Number(m_solver->uses)},
            {"jacobian_build_time",J::Number(m_solver->buildTime)},{"rebuild_reason",J::Quote(reason)},
            {"refreshed_flag",J::Boolean(refreshed)},{"state",J::Numbers(s)},
            {"residual",J::Numbers(r)},{"scaled_block_norms",J::Numbers(norms)},
            {"constraints",constraintJSON(constraints)},
            {"linear_step",stepRecord(s,r,m_solver->jacobian,step,constraints)},
            {"physical",physical(s)},{"trials",searched.observations},{"accepted",J::Boolean(accepted)}}));
        ++m_solver->uses;
        if(accepted)s=std::move(searched.state);
        if(!accepted && !refreshed){m_solver->factor.reset();previous=1e100;continue;}
        Require(accepted,"evolution stage rejected every admissible Newton trial");previous=error;
    }
    throw std::runtime_error("unreachable evolution solve");
    }
    catch(const std::exception &failure)
    {
        if(recording)
        {
            try
            {
                std::map<int,double> constraints;Vector base,baseNorms;
                std::string baseException="null";bool baseSucceeded=false;
                try
                {
                    constraints=Constraints(s,time,data);base=Residual(s,history,alpha,time,data);
                    baseNorms=Norms(base,constraints);baseSucceeded=true;
                }
                catch(const std::exception &error)
                {baseException=J::Object({{"type",J::Quote(J::ExceptionType(error))},
                    {"message",J::Quote(error.what())},
                    {"volume_exception",m_solveResidualException.empty()?"null":m_solveResidualException}});}
                std::string freshJSON="null",actionsJSON="null",retryJSON="null",freshColumnsJSON="[]";
                if(m_solveDiagnosticReplay && baseSucceeded)
                {
                    auto candidate=s;std::vector<std::string> retry;
                    try
                    {
                        WallMatrix fresh;std::string freshFailure="null";
                        try
                        {
                            fresh=m_solveDiagnosticAdmissible?admissibleAssemble(s,constraints):assemble(s,constraints);
                            freshJSON=J::Rows(fresh);
                        }
                        catch(const std::exception &error)
                        {freshFailure=J::Object({{"type",J::Quote(J::ExceptionType(error))},
                            {"message",J::Quote(error.what())},
                            {"invalid_context",m_solveInvalidContext.empty()?"null":m_solveInvalidContext}});}
                        freshColumnsJSON=J::List(diagnosticColumns);
                        std::vector<std::string> directions;
                        for(const std::string kind:{"velocity","free_k","epsilon_with_trace","pressure_and_gauge","joint"})
                        {
                            Vector direction(Size());
                            for(int i=0;i<Size();++i)
                            {
                                int block=0;while(block<Gauge() && i>=offset[block+1])++block;
                                const bool selected=kind=="joint" || (kind=="velocity" && block<m_dimension) ||
                                    (kind=="free_k" && block==K()) || (kind=="epsilon_with_trace" && block==Epsilon()) ||
                                    (kind=="pressure_and_gauge" && block>=InnerPressure());
                                if(selected && (!constraints.count(i) || (m_transitionEnabled && block==Epsilon() && m_traceIds.count(i-offset[Epsilon()]))))
                                    direction[i]=scales[i]*std::sin(1.317*(i+1));
                            }
                            if(m_solveDiagnosticAdmissible && (kind=="free_k" || kind=="joint"))
                            {
                                const auto bubble=fields[K()]->Project([&](const RANSWallCoordinate &p)
                                {
                                    if(p[2]>m_parameters.geometry.height)return 0.;
                                    const double x=(p[0]-lower[0])/(upper[0]-lower[0]);
                                    const double r=p[2]/m_parameters.geometry.height;
                                    return K0*16*x*(1-x)*r*(1-r);
                                });
                                for(int i=0;i<size[K()];++i)direction[offset[K()]+i]=constraints.count(offset[K()]+i)?0:bubble[i];
                            }
                            std::vector<std::string> levels;
                            for(double h:{2e-4,5e-5,1e-5,2e-6,5e-7})
                            {
                                std::vector<std::string> values;
                                for(int sign:{1,-1})
                                {
                                    auto state=s;for(int i=0;i<Size();++i)state[i]+=sign*h*direction[i];
                                    try
                                    {
                                        const auto r=Residual(state,history,alpha,time,data);
                                        values.push_back(J::Object({{"evaluation_succeeded","true"},
                                            {"residual",J::Numbers(r)},{"constraints",constraintJSON(Constraints(state,time,data))},
                                            {"physical",physical(state)}}));
                                    }
                                    catch(const std::exception &error)
                                    {values.push_back(J::Object({{"evaluation_succeeded","false"},
                                        {"type",J::Quote(J::ExceptionType(error))},{"message",J::Quote(error.what())},
                                        {"invalid_context",m_solveInvalidContext.empty()?"null":m_solveInvalidContext},
                                        {"volume_exception",m_solveResidualException.empty()?"null":m_solveResidualException}}));}
                                }
                                levels.push_back(J::Object({{"step",J::Number(h)},{"plus",values[0]},{"minus",values[1]}}));
                            }
                            directions.push_back(J::Object({{"kind",J::Quote(kind)},{"direction",J::Numbers(direction)},
                                {"cached_action",J::Numbers(Multiply(m_solver->jacobian,direction))},
                                {"fresh_action",fresh.empty()?"null":J::Numbers(Multiply(fresh,direction))},{"levels",J::List(levels)}}));
                        }
                        actionsJSON=J::List(directions);
                        if(fresh.empty())retryJSON=J::Object({{"diagnostic_unqualified","true"},
                            {"history_advanced","false"},{"terminal",J::Quote("fresh_jacobian_assembly_unavailable")},
                            {"state",J::Numbers(candidate)},{"assembly_exception",freshFailure}});
                        else
                        {
                        bool converged=false;std::string terminal="iteration_limit";
                        for(int iteration=0;iteration<m_solveDiagnosticIterations;++iteration)
                        {
                            const auto r=Residual(candidate,history,alpha,time,data),cNorms=Norms(r,Constraints(candidate,time,data));
                            const double error=Maximum(cNorms);
                            if(error<tolerance){converged=true;terminal="converged_same_stage";break;}
                            const auto c=Constraints(candidate,time,data);
                            const auto jac=iteration==0?fresh:(m_solveDiagnosticAdmissible?admissibleAssemble(candidate,c):assemble(candidate,c));
                            auto rhs=r;for(double &v:rhs)v=-v;const auto step=Factor(jac).Apply(std::move(rhs));
                            auto searched=search(candidate,step,error,true);
                            retry.push_back(J::Object({{"iteration",J::Number(iteration)},{"state",J::Numbers(candidate)},
                                {"residual",J::Numbers(r)},{"scaled_block_norms",J::Numbers(cNorms)},
                                {"constraints",constraintJSON(c)},{"linear_step",stepRecord(candidate,r,jac,step,c)},
                                {"physical",physical(candidate)},{"trials",searched.observations},
                                {"accepted",J::Boolean(searched.accepted)}}));
                            if(!searched.accepted){terminal="fresh_jacobian_valid_or_invalid_non_descent";break;}
                            candidate=std::move(searched.state);
                        }
                        const auto final=Residual(candidate,history,alpha,time,data);
                        const auto finalNorms=Norms(final,Constraints(candidate,time,data));
                        if(Maximum(finalNorms)<tolerance){converged=true;terminal="converged_same_stage";}
                        retryJSON=J::Object({{"diagnostic_unqualified","true"},{"history_advanced","false"},
                            {"maximum_iterations",J::Number(m_solveDiagnosticIterations)},
                            {"converged",J::Boolean(converged)},{"terminal",J::Quote(terminal)},
                            {"iterations",J::List(retry)},{"state",J::Numbers(candidate)},
                            {"residual",J::Numbers(final)},{"scaled_block_norms",J::Numbers(finalNorms)},
                            {"constraints",constraintJSON(Constraints(candidate,time,data))},{"physical",physical(candidate)}});
                        }
                    }
                    catch(const std::exception &error)
                    {retryJSON=J::Object({{"diagnostic_unqualified","true"},{"history_advanced","false"},
                        {"terminal",J::Quote("diagnostic_evaluation_or_linear_failure")},
                        {"state",J::Numbers(candidate)},{"iterations",J::List(retry)},
                        {"type",J::Quote(J::ExceptionType(error))},{"message",J::Quote(error.what())}});}
                }
                const auto record=J::Object({{"schema","1"},{"kind",J::Quote("first_failed_nonlinear_stage")},
                    {"diagnostic_unqualified","true"},{"original_policy_unchanged","true"},
                    {"admissible_dimensionless_diagnostic_columns",J::Boolean(m_solveDiagnosticAdmissible)},
                    {"diagnostic_column_admissibility_scope",J::Quote("residual evaluation domain only; complete physical admissibility is checked separately")},
                    {"failure_type",J::Quote(J::ExceptionType(failure))},{"failure",J::Quote(failure.what())},
                    {"time",J::Number(time)},{"alpha",J::Number(alpha)},{"tolerance",J::Number(tolerance)},
                    {"maximum_iterations",J::Number(maximumIterations)},{"condensed",J::Boolean(condensed)},
                    {"state",J::Numbers(s)},{"history",J::Numbers(history)},
                    {"constraints",constraintJSON(constraints)},{"residual",J::Numbers(base)},
                    {"base_evaluation_succeeded",J::Boolean(baseSucceeded)},{"base_exception",baseException},
                    {"scaled_block_norms",J::Numbers(baseNorms)},
                    {"dimensionless_variable_scales",J::Numbers(scales)},{"offsets",J::Numbers(offset)},
                    {"sizes",J::Numbers(size)},{"stage_entry_builds",J::Number(entryBuilds)},
                    {"stage_entry_jacobian_age",J::Number(entryAge)},
                    {"cached_jacobian",J::Rows(m_solver->jacobian)},
                    {"cached_jacobian_state",J::Numbers(m_solver->buildState)},
                    {"cached_jacobian_history",J::Numbers(m_solver->buildHistory)},
                    {"cached_jacobian_time",J::Number(m_solver->buildTime)},
                    {"cached_jacobian_iteration",J::Number(m_solver->buildIteration)},
                    {"assembly_failure",assemblyFailure},{"iterations",J::List(records)},
                    {"physical",physical(s)},{"fresh_jacobian",freshJSON},
                    {"diagnostic_jacobian_columns",freshColumnsJSON},
                    {"last_retry_jacobian_columns",J::List(diagnosticColumns)},
                    {"directional_checks",actionsJSON},{"fresh_jacobian_retry",retryJSON}});
                Require(!std::ifstream(m_solveDiagnosticFile).good(),"refusing to overwrite an existing nonlinear diagnostic");
                std::ofstream out(m_solveDiagnosticFile);Require(out.good(),"cannot open nonlinear diagnostic output");
                out<<record<<'\n';out.close();Require(out.good(),"cannot finish nonlinear diagnostic output");
                m_solveDiagnosticWritten=true;
            }
            catch(const std::exception &error)
            {std::cerr<<"EP-AWF nonlinear diagnostic export: "<<error.what()<<std::endl;}
        }
        throw;
    }
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::OwnedBodyLoad(double time,const RANSWallEvolutionData &data,int region) const
{
    Vector load(Size());
    auto add=[&](const VolumePoint &q,bool inner)
    {
        auto H=q.v.coordinate;H[2]=m_parameters.geometry.height;
        for(int a=0;a<m_dimension;++a)
        {
            double body=(inner?data.innerBody:data.outerBody)[a](q.v.coordinate,time);
            if(inner && Distance(q.v.coordinate)<=m_parameters.geometry.roughness)
            {
                const auto n=m_mapping.At(H).normal;double bn=0;
                for(int b=0;b<m_dimension;++b)bn+=data.outerBody[b](H,time)*n[PhysicalAxis(b)];
                body-=data.outerBody[a](H,time)-bn*n[PhysicalAxis(a)];
            }
            Require(std::isfinite(body),"nonfinite prescribed body force");
            if(body!=0.)Add(load,offset[a],q.v,body,{},q.weight);
        }
    };
    // Prescribed IMEX body loads use one fixed physical quadrature, including
    // with transition sensors. Otherwise K-dependent panels would change the
    // explicit load between the accepted state and the implicit stage.
    if(region!=1)for(const auto &q:m_bodyInner)add(q,true);
    if(region!=0)for(const auto &q:m_outer)add(q,false);
    return load;
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::ExplicitBody(const Vector &s,double time,const RANSWallEvolutionData &data) const
{
    Vector load(Size());int bad=0;
    try{load=OwnedBodyLoad(time,data);}catch(const std::exception &){bad=1;}
    m_comm->AllReduce(bad,LibUtilities::ReduceMax);Require(!bad,"inadmissible prescribed body load on a collective rank");
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
    Require(!m_transitionEnabled || eta==0,"every upstream transition/terrain run requires eta=0");
    const auto raw=Residual(ref,ref,0,time,data,false);const auto norms=Norms(raw,Constraints(ref,time,data));
    for(int b=InnerPressure();b<norms.size();++b)
        if(!((Mapped() || m_transitionEnabled) && eta==0 && b<Gauge()))Require(norms[b]<1e-9,"reference violates an algebraic constraint");
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
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::RoughnessBudget(const Vector &state,double time,const RANSWallEvolutionData &data) const
{
    const auto s=Lift(state,time,data);Refresh(s,time,&data);Vector sum(5);double normal=0;
    for(const auto &ray:m_rays)
    {
        const auto force=CartesianForcing(s,ray.coordinate,time,data,true);const auto n=m_mapping.At(ray.coordinate).normal;
        double dot=0;for(int a=0;a<3;++a)dot+=force[a]*n[a];normal=std::max(normal,std::abs(dot));
        for(const auto &q:ray.points)if(Distance(q.v.coordinate)<=m_parameters.geometry.roughness)
        {
            for(int a=0;a<3;++a)sum[a]+=q.weight*force[a];
            for(int a=0;a<m_dimension;++a)sum[3]+=q.weight*q.v.Value(s,offset[a])*force[PhysicalAxis(a)];
        }
    }
    m_comm->AllReduce(sum,LibUtilities::ReduceSum);m_comm->AllReduce(normal,LibUtilities::ReduceMax);sum[4]=normal;return sum;
}
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::EnergyBudget(const Vector &state,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data) const
{
    // R4 reports direct and weak work on the same retained epsilon-H state.
    const auto s=GroundSafeguard()?ExternalLift(state,time,data):Lift(state,time,data);
    const auto raw=Residual(state,history,alpha,time,data,false,false);Vector result(9);
    for(int a=0;a<m_dimension;++a)for(int i=0;i<size[a];++i)result[0]+=s[offset[a]+i]*raw[offset[a]+i];
    for(int i=0;i<size[K()];++i)result[0]+=m_constants[K()][i]*raw[offset[K()]+i];
    auto integrate=[&](const VolumePoint &q,bool inner,const EquilibriumAWFProfiles *profile,const RANSVector &force,
        const RANSCappedDestructionValues *normalization)
    {
        RANSTensor G{},stress{};RANSVector u{};
        for(int a=0;a<m_dimension;++a){const int aa=PhysicalAxis(a);u[aa]=q.v.Value(s,offset[a]);G[aa]=q.v.Gradient(s,offset[a]);}
        const double k=q.v.Value(s,offset[K()]),pi=q.p.Value(s,offset[inner?InnerPressure():OuterPressure()]);double sk;
        if(inner)
        {
            const double nt=profile->Original().EddyViscosity(Distance(q.v.coordinate));
            const auto innerStress=InnerStress(k,nt,G);
            stress=RANSWallInnerMeanStress(innerStress,m_parameters.geometry.molecularNu,pi);
            sk=WallSource(s,q.v.coordinate,innerStress.production,*profile,normalization).kSource+data.innerKSource(q.v.coordinate,time);
            if(Distance(q.v.coordinate)<=m_parameters.geometry.roughness)for(int a=0;a<3;++a)result[4]-=q.weight*u[a]*force[a];
        }
        else
        {
            const auto c=OuterClosure(k,q.e.Value(s,offset[Epsilon()]),G,q.v.coordinate);
            sk=c.kSource+data.outerKSource(q.v.coordinate,time);
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)stress[a][b]=-c.reynoldsStress[a][b]+(a==b?2*k/3-pi:0);
        }
        result[2]-=q.weight*sk;
        for(int a=0;a<3;++a)for(int b=0;b<3;++b){result[2]+=q.weight*G[a][b]*stress[a][b];result[3]-=q.weight*G[a][b]*u[a]*u[b];}
    };
    for(const auto &ray:m_rays){const auto f=CartesianForcing(s,ray.coordinate,time,data,false);for(const auto &q:ray.points)integrate(q,true,ray.profile.get(),f,m_cappedDestruction?&ray.destruction:nullptr);}
    for(const auto &q:m_outer)integrate(q,false,nullptr,{},nullptr);
    const auto body=OwnedBodyLoad(time,data);for(int a=0;a<m_dimension;++a)for(int j=0;j<size[a];++j)result[4]-=s[offset[a]+j]*body[offset[a]+j];
    for(int i=1;i<result.size();++i){double v=result[i];m_comm->AllReduce(v,LibUtilities::ReduceSum);result[i]=v;}
    // Storage is the actual physical mass operator used by the stage, not
    // the possibly lower-order constitutive quadrature. These replicated
    // global matrices are contracted after the owned-volume reductions.
    for(int b=0;b<=K();++b)
    {
        Vector current(size[b]),delta(size[b]);
        for(int i=0;i<size[b];++i){current[i]=s[offset[b]+i];delta[i]=current[i]-history[offset[b]+i];}
        const auto change=Multiply(m_mass[b],delta),stored=Multiply(m_mass[b],current);
        for(int i=0;i<size[b];++i)
        {
            const double test=b==K()?m_constants[b][i]:current[i];
            result[1]+=alpha*test*change[i];result[8]+=(b==K()?1.:.5)*test*stored[i];
        }
    }
    result[6]=result[1]+result[2]+result[3]+result[4];result[7]=result[0]-result[6];return result;
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
    if(Mapped())text<<':'<<m_mapping.Signature();
    if(m_innerStress==RANSWallInnerStressModel::SaturatedR2)
        text<<":M8-R2:InnerStress:SaturatedR2:v1:rho=0.8";
    if(m_variantName!="Original" || m_cappedDestruction)
    {
        text<<":M8:"<<m_variantName;
        for(double v:RobustParameters())text<<':'<<v;
        for(double v:m_referenceFlow)text<<':'<<v;
        for(double v:m_referenceNormal)text<<':'<<v;
        text<<':'<<m_robustUstar<<":incoming-flat-strain-extension";
    }
    if(m_transitionEnabled)
    {
        text<<":upstream:"<<m_transition.start<<':'<<m_transition.end<<':'<<m_sensorA<<':'<<m_sensorB<<':'<<m_baselineStart;
        const auto &c=m_transition.state;text<<':'<<c.on<<':'<<c.off<<':'<<c.marginOff<<':'<<c.marginOn<<':'<<c.reynoldsOff<<':'<<c.reynoldsOn;
        for(double w:c.weights)text<<':'<<w;
        for(double v:m_transitionReference)text<<':'<<v;
    }
    return text.str();
}
}
