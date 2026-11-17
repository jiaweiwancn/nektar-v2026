// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupling.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeWallV2.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeFailures.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeAdmissibility.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeInitialFields.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeFlowTangent.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace Nektar::WJNative
{
namespace
{
namespace D=RANSWallR3;
using Triplet=Eigen::Triplet<double>;
void Require(bool ok,const std::string&message){if(!ok)throw std::runtime_error(message);}
double Clock(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
Vec Local(const Vec&state,const std::vector<int>&ids,int offset=0)
{Vec result(ids.size());for(int i=0;i<result.size();++i)result[i]=state[offset+ids[i]];return result;}
void Add(Vec&state,const std::vector<int>&ids,const Vec&local,int offset=0)
{for(int i=0;i<local.size();++i)state[offset+ids[i]]+=local[i];}
std::array<double,2>CubicBounds(const std::array<double,4>&c)
{
    // Modified_A on [-1,1] written in t=(xi+1)/2 in [0,1].
    double a=c[0],b=-c[0]+c[1]+c[2]-2*c[3],d=-c[2]+6*c[3],e=-4*c[3];
    auto value=[&](double t){return a+t*(b+t*(d+t*e));};
    std::array<double,2>bounds{std::min(value(0),value(1)),std::max(value(0),value(1))};
    auto include=[&](double t){if(t>0&&t<1){double v=value(t);bounds[0]=std::min(bounds[0],v);bounds[1]=std::max(bounds[1],v);}};
    if(e==0){if(d!=0)include(-b/(2*d));}
    else{double disc=4*d*d-12*e*b;if(disc>=0){double root=std::sqrt(disc),q=-d-std::copysign(root/2,d);if(q==0)include(-d/(3*e));else{include(q/(3*e));include(b/q);}}}
    return bounds;
}
std::array<double,2>CubicTrace(const std::array<double,4>&c,double t,double width)
{
    // Exact native Modified_A Q3 edge polynomial and physical x derivative.
    const double bubble=t*(1-t),linear=c[2]+(4*t-2)*c[3];
    return {(1-t)*c[0]+t*c[1]+bubble*linear,
        (-c[0]+c[1]+(1-2*t)*linear+4*bubble*c[3])/width};
}
void VectorJSON(std::ostream&out,const Vec&v)
{out<<'[';for(int i=0;i<v.size();++i){if(i)out<<',';out<<v[i];}out<<']';}
std::string StateHash(const CoupledState&s)
{D::Sha256 hash;for(const auto*v:{&s.flow,&s.k,&s.epsilon})hash.Update(v->data(),v->size()*sizeof(double));return hash.Finish();}
CoupledState Combination(const CoupledState&a,double wa,const CoupledState&b,double wb)
{return {wa*a.flow+wb*b.flow,wa*a.k+wb*b.k,wa*a.epsilon+wb*b.epsilon};}
}

Coupling::Coupling(std::vector<MultiRegions::ContFieldSharedPtr>fields,Geometry geometry,
    std::shared_ptr<const Reference>reference,CoupledOptions options)
    :m_fields(std::move(fields)),m_geometry(geometry),m_reference(std::move(reference)),m_options(std::move(options)),
     m_flow({m_fields.at(0),m_fields.at(1),m_fields.at(4)},geometry),
     m_rules(m_flow.Core(),m_options.quadratureX,m_options.quadratureR),
     m_k(m_fields.at(2),m_flow.Core(),[this](double,double){return m_reference->Parameters().KIn;},{0}),
     m_e(m_fields.at(3),m_flow.Core(),[this](double,double z){return m_reference->Amplitude(m_reference->Parameters().KIn)/(z-m_geometry.d0);},{0,2})
{
    Require(m_fields.size()==5&&bool(m_reference),"coupling requires native u,v,k,outer-epsilon,p and a qualified reference");
    Require(m_options.dt>0&&m_options.linearTolerance>0&&m_options.linearTolerance<=1e-10&&
        m_options.nonlinearTolerance>0&&m_options.nonlinearTolerance<=1e-8&&m_options.maximumIterations>=2,
        "invalid coupled timestep or convergence policy");
    Require(m_options.sourceIdentity.size()==64,"coupled scientific source identity is missing");
    Require(m_options.andersonDepth>=0&&m_options.andersonDepth<=8,"Anderson depth must lie in [0,8]");
    Require(m_options.stressSplit=="PositiveShearV1"||m_options.stressSplit=="IsotropicMuV1","unknown numerical stress split");
    Require(m_options.steadyVelocityPC=="FullGradientV1"||m_options.steadyVelocityPC=="PositiveFrozenV1","unknown steady velocity preconditioner");
    Require(m_options.advancement=="PicardV1"||m_options.advancement=="CoupledJFNKV1","unknown implicit coupled solver");
    Require(m_options.stepPreconditioner=="RebuildEveryNewtonV1"||m_options.stepPreconditioner=="OperatorDriftV1"||m_options.stepPreconditioner=="CurrentDataDriftV1",
        "unknown implicit Newton preconditioner policy");
    Require(std::isfinite(m_options.stepPCRelativeDrift)&&m_options.stepPCRelativeDrift>=0&&m_options.stepPCRelativeDrift<=.25&&
        m_options.stepPCMaximumUses>=1&&m_options.stepPCMaximumUses<=64,"invalid Newton preconditioner refresh policy");
    Require(m_options.advancement=="CoupledJFNKV1"||m_options.stepPreconditioner=="RebuildEveryNewtonV1",
        "Newton preconditioner reuse is not a Picard option");
    if(m_options.advancement=="CoupledJFNKV1"){
        Require(m_options.andersonDepth==0,"Anderson acceleration is a Picard-only option");
        Require(std::isfinite(m_options.stepNewton.maximumSeconds)&&m_options.stepNewton.maximumSeconds>0&&
            m_options.stepNewton.acceptanceNorm==m_options.nonlinearTolerance,"implicit Newton requires explicit time cap and current coupled tolerance");
    }
    Require(m_geometry.curve==0&&m_geometry.H==m_reference->Parameters().H&&m_geometry.h==m_reference->Parameters().roughness,
        "coupled v2 supports only the pinned flat geometry/model heights");
    for(const auto&field:m_fields)for(int e=0;e<field->GetExpSize();++e)for(int d=0;d<2;++d)
        Require(field->GetExp(e)->GetBasis(d)->GetBasisType()==LibUtilities::eModified_A,"coupled fields require the declared native Modified_A basis");
    m_flow.SetTransportedKSpace(m_fields[2]);BuildTopology();
    Require(std::isfinite(m_options.mmsBodyAmplitude),"nonfinite MMS body amplitude");
    Require(m_options.manufactured||m_options.mmsBodyAmplitude==0,"the prescribed MMS body is verification-only");
    Require(m_options.mmsVariant=="Historical"||m_options.mmsVariant=="PolynomialScalarsV1","unknown continuous MMS variant");
    Require(m_options.manufactured||m_options.mmsVariant=="Historical","a continuous MMS variant is verification-only");
    if(m_options.manufactured)m_mms=std::make_unique<CoupledMMS>(*m_reference,m_geometry,m_top,m_options.mmsBodyAmplitude,
        m_options.mmsVariant=="Historical"?MMSVariant::Historical:MMSVariant::PolynomialScalarsV1);
    else Require(std::abs(m_geometry.d0-m_reference->Parameters().d0)<1e-14&&
        std::abs(m_top-m_reference->Parameters().topRadius)<1e-12&&std::abs(m_left+10)<1e-12&&std::abs(m_right-15)<1e-12,
        "application geometry differs from the retained physical domain");
    std::ostringstream signature;signature<<std::setprecision(17)<<"WJ-H27-coupled-v2-BDF2-EX2|"<<WallModelV2<<'|'
        <<m_options.sourceIdentity<<'|'<<m_reference->Identity()<<'|'<<m_options.dt<<'|'<<m_options.manufactured<<'|'
        <<m_options.quadratureX<<'|'<<m_options.quadratureR<<'|'<<m_geometry.d0<<'|'<<m_geometry.H<<'|'<<m_geometry.h<<'|'
        <<m_options.linearTolerance<<'|'<<m_options.nonlinearTolerance<<'|'<<m_options.maximumIterations<<'|'<<m_options.mmsBodyAmplitude<<'|'<<int(m_options.backend)<<'|'<<m_options.mmsVariant<<'|'<<m_options.andersonDepth<<'|'<<m_options.stressSplit<<'|'<<m_options.steadyVelocityPC;
    const auto&newton=m_options.stepNewton;
    signature<<'|'<<m_options.advancement<<'|'<<newton.maximumSeconds<<'|'<<newton.acceptanceNorm<<'|'
        <<newton.maximumNewtonIterations<<'|'<<newton.restartDimension<<'|'<<newton.maximumLinearCycles<<'|'<<newton.fdRelativeStep;
    signature<<'|'<<m_options.stepPreconditioner<<'|'<<m_options.stepPCRelativeDrift<<'|'<<m_options.stepPCMaximumUses;
    for(const auto&e:m_flow.Elements()){
        signature<<'|'<<e.geometry<<','<<e.xl<<','<<e.xr<<','<<e.rl<<','<<e.rr;
        for(int a:e.u)signature<<','<<a;for(double a:e.us)signature<<','<<a;
        for(int a:e.p)signature<<','<<a;for(double a:e.ps)signature<<','<<a;
    }
    for(const auto&ids:m_flow.FreeVelocityIds()){signature<<"|component-mask:"<<ids.size()<<':';for(int a:ids)signature<<','<<a;}
    for(const auto*scalar:{&m_k,&m_e}){
        signature<<'|'<<scalar->Size();for(int a:scalar->FreeIds())signature<<','<<a;
        const auto&field=scalar->Field();auto map=field->GetLocalToGlobalMap();
        for(int j=0;j<field->GetNcoeffs();++j)signature<<','<<map->GetLocalToGlobalMap(j)<<','<<map->GetLocalToGlobalSign(j);
    }
    m_signature=D::Hash(signature.str());
    CheckTraceEvaluation();
}

void Coupling::BuildTopology()
{
    const auto&elements=m_flow.Elements();Require(!elements.empty(),"empty coupled mesh");
    m_left=1e100;m_right=-1e100;m_top=0;
    std::map<int,std::vector<std::pair<int,int>>>edges;
    for(int i=0;i<int(elements.size());++i){const auto&e=elements[i];Require(m_geometryIndex.emplace(e.geometry,i).second,"duplicate coupled geometry");
        Require(e.rl>=-1e-12&&!(e.rl<m_geometry.H-1e-12&&e.rr>m_geometry.H+1e-12),"H must be a complete native mesh-edge layer");
        Require(m_e.Owns(e.geometry)==(e.rl>=m_geometry.H-1e-12),"outer epsilon inventory differs from native H layer");
        m_left=std::min(m_left,e.xl);m_right=std::max(m_right,e.xr);m_top=std::max(m_top,e.rr);m_volume+=(e.xr-e.xl)*(e.rr-e.rl);
        for(int side=0;side<4;++side)edges[e.velocity->GetGeom()->GetEid(side)].push_back({i,side});
    }
    Require(std::abs(m_volume-(m_right-m_left)*m_top)<1e-10*m_volume,"flow mesh is not the complete supported rectangle");
    constexpr int sideFromNative[4]={2,1,3,0};
    int interfaces=0;
    for(const auto&entry:edges){const auto&owners=entry.second;Require(owners.size()==1||owners.size()==2,"native edge is nonmanifold");
        if(owners.size()==1){auto[e,side]=owners[0];m_boundaryFaces.emplace_back(e,sideFromNative[side]);}
        else{
            const auto&a=elements[owners[0].first];const auto&b=elements[owners[1].first];
            if(m_e.Owns(a.geometry)!=m_e.Owns(b.geometry)){
                const auto&outer=m_e.Owns(a.geometry)?a:b;const auto&inner=m_e.Owns(a.geometry)?b:a;
                Require(std::abs(outer.rl-m_geometry.H)<1e-12&&std::abs(inner.rr-m_geometry.H)<1e-12&&
                    std::abs(outer.xl-inner.xl)<1e-12&&std::abs(outer.xr-inner.xr)<1e-12,"inconsistent native H edge ownership");++interfaces;
            }
        }
    }
    Require(interfaces>0,"no native inner/outer H edges");
    // Boundary-region numbers in a retained XML are arbitrary. Check the
    // component masks against the actual exterior topology, using the two
    // interior edge modes so a neighbouring corner law cannot mask an error.
    constexpr int edgeModes[4][2]={{8,12},{9,13},{2,3},{6,7}};
    std::array<std::set<int>,2> free;
    for(int component=0;component<2;++component)
        free[component].insert(m_flow.FreeVelocityIds()[component].begin(),m_flow.FreeVelocityIds()[component].end());
    for(auto [index,side]:m_boundaryFaces){const auto&e=elements[index];
        const double position=side==0?e.xl:side==1?e.xr:side==2?e.rl:e.rr;
        const double expected=side==0?m_left:side==1?m_right:side==2?0.:m_top;
        Require(std::abs(position-expected)<1e-11,"unsupported nonrectangular exterior boundary");
        for(int component=0;component<2;++component){
            const bool essential=m_options.manufactured||side==0||side==2||(side==3&&component==1);
            for(int mode:edgeModes[side])Require((free[component].count(e.u[mode])==0)==essential,
                "native component boundary mask disagrees with the declared physical side law");
        }
    }
    Require(m_flow.Gauge()==m_options.manufactured,
        "actual pressure nullspace disagrees with the supported closed-MMS/open-outlet boundary law");
    for(int index=0;index<int(elements.size());++index){int owner=index;std::set<int>visited;
        while(std::abs(elements[owner].rl-m_geometry.H)>1e-12){
            Require(visited.insert(owner).second,"cyclic native H owner traversal");
            const int side=elements[owner].rl<m_geometry.H?2:0;
            const auto&adjacent=edges.at(elements[owner].velocity->GetGeom()->GetEid(side));
            Require(adjacent.size()==2,"native column does not reach H through internal edges");
            owner=adjacent[0].first==owner?adjacent[1].first:adjacent[0].first;
        }
        const auto&e=elements[index];const auto&a=elements[owner];
        Require(std::abs(a.xl-e.xl)<1e-11&&std::abs(a.xr-e.xr)<1e-11&&m_e.Owns(a.geometry),"unsupported noncolumnar H trace ownership");
        m_hOwner[e.geometry]=owner;
    }
    auto field=m_fields[2];auto map=field->GetLocalToGlobalMap();
    for(int i=0;i<field->GetExpSize();++i){auto exp=field->GetExp(i);const auto&e=Owner(exp->GetGeom()->GetGlobalID());if(std::abs(e.rl-m_geometry.H)>1e-12)continue;
        Require(exp->GetBasis(0)->GetBasisType()==LibUtilities::eModified_A&&exp->GetBasis(1)->GetBasisType()==LibUtilities::eModified_A,
            "direct H trace requires the declared native Modified_A basis");
        Require(e.velocity->GetBasis(0)->GetBasisType()==LibUtilities::eModified_A&&e.velocity->GetBasis(1)->GetBasisType()==LibUtilities::eModified_A,
            "direct velocity H trace requires native Modified_A");
        KTrace trace;for(int j=0;j<4;++j){const int local=field->GetCoeff_Offset(i)+j;trace.ids[j]=map->GetLocalToGlobalMap(local);trace.signs[j]=map->GetLocalToGlobalSign(local);}m_hKTrace[e.geometry]=trace;
    }
}
const Element&Coupling::Owner(int geometry)const{return m_flow.Elements().at(m_geometryIndex.at(geometry));}
const Element&Coupling::HOwner(const Element&e)const{return m_flow.Elements().at(m_hOwner.at(e.geometry));}
NativeRule Coupling::HRule(const NativeRule&rule)const
{return m_rules.Points(HOwner(Owner(rule.geometry)),rule.x,Vec::Constant(rule.x.size(),m_geometry.d0+m_geometry.H));}
double Coupling::InterfaceK(const Element&e,double x,const CoupledState&s)const
{
    const auto&owner=HOwner(e);const auto&trace=m_hKTrace.at(owner.geometry);std::array<double,4>c{};
    for(int j=0;j<4;++j)c[j]=s.k[trace.ids[j]]*trace.signs[j];
    const double K=CubicTrace(c,(x-owner.xl)/(owner.xr-owner.xl),owner.xr-owner.xl)[0];
    RequireState(std::isfinite(K)&&K>=m_reference->Parameters().KMinimum&&K<=m_reference->Parameters().KMaximum,"current sampled interface K left the qualified reference interval");m_reference->Cutoff(K);return K;
}
std::array<double,2>Coupling::InterfaceVelocity(const Element&e,double x,const CoupledState&s)const
{
    const auto&owner=HOwner(e);std::array<double,4>c{};for(int j=0;j<4;++j)c[j]=s.flow[owner.u[j]]*owner.us[j];
    return CubicTrace(c,(x-owner.xl)/(owner.xr-owner.xl),owner.xr-owner.xl);
}
void Coupling::CheckTraceEvaluation()const
{
    CoupledState synthetic;synthetic.k=m_k.Project([&](double,double){return m_reference->Parameters().KIn;});
    synthetic.flow=Vec(2*m_flow.VelocitySize()+m_flow.PressureSize());
    for(int j=0;j<synthetic.k.size();++j)synthetic.k[j]+=1e-4*m_reference->Parameters().KIn*std::sin(double(j+1));
    for(int j=0;j<synthetic.flow.size();++j)synthetic.flow[j]=std::cos(.7*(j+1));
    for(const auto&entry:m_hKTrace){const auto&e=Owner(entry.first);for(double t:{0.,.135,.555,.8,1.}){
        double x=e.xl+t*(e.xr-e.xl);auto rule=m_rules.Points(e,Vec::Constant(1,x),Vec::Constant(1,m_geometry.d0+m_geometry.H));
        const double K=m_k.Evaluate(rule,synthetic.k)[0],fast=InterfaceK(e,x,synthetic);
        Require(std::abs(K-fast)<=2e-13*m_reference->Parameters().KIn,"direct native H scalar trace disagrees with owned basis evaluation");
        auto velocity=InterfaceVelocity(e,x,synthetic);Vec local=Local(synthetic.flow,rule.u);
        Require(std::abs((rule.B*local)[0]-velocity[0])<=2e-13&&std::abs((rule.Dx*local)[0]-velocity[1])*(e.xr-e.xl)<=2e-13,
            "direct native H velocity trace/gradient disagrees with owned basis evaluation");
    }}
}
double Coupling::NativeRadius(double physicalZ)const
{
    const double r=physicalZ-m_geometry.d0,top=m_options.manufactured?.5:m_reference->Parameters().topRadius;
    // Trusted native geometry coordinates only. Scientific Reference inputs
    // remain strict: no state, parameter, K interval or cutoff is clamped.
    for(double edge:{0.,m_geometry.H,top})if(std::abs(r-edge)<=64*std::numeric_limits<double>::epsilon()*std::max({m_geometry.H,std::abs(r),std::abs(edge)}))return edge;
    return r;
}
double Coupling::NativeX(double x)const
{
    const double left=m_options.manufactured?0.:-10.,right=m_options.manufactured?1.:15.;
    for(double edge:{left,right})if(std::abs(x-edge)<=64*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(x),std::abs(edge)}))return edge;
    return x;
}

Vec Coupling::FlowLift(double time)const
{
    if(!m_mms&&m_physicalFlowLift.size())return m_physicalFlowLift;
    Vec full=Vec::Zero(2*m_flow.VelocitySize()+m_flow.PressureSize());
    for(int component=0;component<2;++component){auto f=[&,component](double x,double z){
        if(m_mms){auto a=m_mms->Fields(NativeX(x),NativeRadius(z)+m_geometry.d0,time);return component==0?a.u.value:a.w.value;}
        return component==0?m_reference->InletU(NativeRadius(z)):0.;};
        full.segment(component*m_flow.VelocitySize(),m_flow.VelocitySize())=m_flow.Project(f,component);}
    Vec lift=m_flow.Lift(full);if(!m_mms)m_physicalFlowLift=lift;return lift;
}
CoupledState Coupling::DirichletState(const CoupledState&s,double time)const
{
    CoupledState result=s;result.flow=m_flow.Expand(m_flow.Restrict(s.flow),FlowLift(time));
    Vec kLift;
    if(!m_mms&&m_physicalKLift.size())kLift=m_physicalKLift;
    else{kLift=m_k.EssentialLift([&](double x,double z,const Vector&){return m_mms?m_mms->Fields(NativeX(x),NativeRadius(z)+m_geometry.d0,time).k.value:m_reference->InletK(NativeRadius(z));},m_rules);if(!m_mms)m_physicalKLift=kLift;}
    result.k=m_k.Expand(m_k.Restrict(s.k),kLift);
    auto eLift=m_e.EssentialLift([&](double x,double z,const Vector&n){
        if(n[2]<-.5){
            const Element*owner=nullptr;for(const auto&e:m_flow.Elements())if(std::abs(e.rl-m_geometry.H)<1e-12&&x>=e.xl-1e-12&&x<=e.xr+1e-12){owner=&e;break;}
            Require(owner!=nullptr,"dependent epsilon trace has no outer owner");
            return m_reference->Amplitude(InterfaceK(*owner,x,result))/m_geometry.H;
        }
        return m_mms?m_mms->OuterEpsilon(NativeX(x),NativeRadius(z)+m_geometry.d0,time).value:m_reference->InletE(NativeRadius(z));
    },m_rules);
    result.epsilon=m_e.Expand(m_e.Restrict(s.epsilon),eLift);return result;
}

CoupledHistory Coupling::Initialise()const
{
    CoupledHistory h;h.dt=m_options.dt;
    h.current.flow=Vec::Zero(2*m_flow.VelocitySize()+m_flow.PressureSize());
    for(int component=0;component<3;++component){auto f=[&,component](double x,double z){
            if(m_mms){auto a=m_mms->Fields(NativeX(x),NativeRadius(z)+m_geometry.d0,0);return component==0?a.u.value:component==1?a.w.value:a.Pi.value;}
            if(component==0)return m_reference->InletU(NativeRadius(z));
            if(component==1)return 0.;
            // A fresh normal-traction-compatible seed, with no imposed old
            // normal-stress jump. This is not a stationary-flow assertion.
            double r=NativeRadius(z),K=m_reference->Parameters().KIn,k=m_reference->InletK(r);
            double E=r<m_geometry.H?m_reference->Ew(K,r):m_reference->InletE(r);Tensor G{};G[0][2]=m_reference->InletUGradient(r);
            auto a=AdaptAmended(k,E,G,{0,0,1},m_reference->Parameters().nu,r<m_geometry.H?m_reference->Aux(K,r):0.,r,m_geometry.H,m_options.backend,&m_constitutiveStats);
            return a.V[2][2]-K/m_reference->Parameters().Ne;
        };
        Vec value=m_flow.Project(f,component);
        if(component<2)h.current.flow.segment(component*m_flow.VelocitySize(),m_flow.VelocitySize())=value;
        else h.current.flow.tail(m_flow.PressureSize())=value;
    }
    h.current.k=m_k.Project([&](double x,double z){return m_mms?m_mms->Fields(NativeX(x),NativeRadius(z)+m_geometry.d0,0).k.value:m_reference->InletK(NativeRadius(z));});
    h.current.epsilon=m_e.Project([&](double x,double z){return m_mms?m_mms->OuterEpsilon(NativeX(x),NativeRadius(z)+m_geometry.d0,0).value:m_reference->InletE(NativeRadius(z));});
    if(!m_mms)m_physicalFlowLift=m_flow.Lift(h.current.flow);
    return FreshHistory(h.current);
}

CoupledHistory Coupling::FreshHistory(const CoupledState&state)const
{
    CoupledHistory h;h.dt=m_options.dt;h.current=DirichletState(state,0);
    if(m_flow.Gauge())h.current.flow.tail(m_flow.PressureSize())-=m_flow.PressureOne()*
        (m_flow.PressureMean().dot(h.current.flow.tail(m_flow.PressureSize()))/m_volume);
    h.previous=h.current;CoupledTiming timing;auto data=Evaluate(h.current,0,true,timing);
    h.advection=Advection(h.current,data.kRules);h.previousAdvection=h.advection;CheckHistory(h);return h;
}

CoupledHistory Coupling::ReadInitialFields(const std::string&path,bool allowCandidate)
{
    InitialFields initial(path);InitialFieldGeometry geometry;
    geometry.d0=m_geometry.d0;geometry.H=m_geometry.H;geometry.h=m_geometry.h;geometry.left=m_left;geometry.right=m_right;geometry.top=m_top;geometry.curve=m_geometry.curve;
    const std::string expectedSource=m_options.initialSourceIdentity.empty()?m_options.sourceIdentity:m_options.initialSourceIdentity;
    Require(expectedSource.size()==64,"initial-field source binding is missing");
    Require(expectedSource==m_options.sourceIdentity||allowCandidate,"a different declared initial source is allowed only as a steady-preparation guess");
    initial.ValidateBindings(expectedSource,m_reference->Identity(),WallModelV2,m_options.manufactured,m_options.mmsVariant,geometry);
    Require(initial.Metadata().gauge==m_flow.Gauge(),"initial-field pressure-gauge binding differs");
    Require(allowCandidate||initial.Metadata().qualifiedReference,"advancement requires a qualified reference initial packet");
    std::vector<InitialFieldElement>inventory;
    for(const auto&e:m_flow.Elements()){
        InitialFieldElement cell;cell.geometry=e.geometry;cell.xl=e.xl;cell.xr=e.xr;cell.rl=e.rl;cell.rr=e.rr;cell.u=e.u;cell.us=e.us;cell.p=e.p;cell.ps=e.ps;
        for(int field:{2,3}){auto f=m_fields[field];auto map=f->GetLocalToGlobalMap();auto&ids=field==2?cell.k:cell.epsilon;auto&signs=field==2?cell.ks:cell.epsilonSigns;
            for(int i=0;i<f->GetExpSize();++i)if(f->GetExp(i)->GetGeom()->GetGlobalID()==e.geometry){for(int j=0;j<16;++j){int local=f->GetCoeff_Offset(i)+j;ids.push_back(map->GetLocalToGlobalMap(local));signs.push_back(map->GetLocalToGlobalSign(local));}break;}}
        inventory.push_back(std::move(cell));
    }
    const bool same=initial.SameInventory(inventory,{m_flow.VelocitySize(),m_flow.PressureSize(),m_k.Size(),m_e.Size()});
    CoupledState state;
    if(same)state={initial.Flow(),initial.K(),initial.Epsilon()};
    else{
        state.flow=Vec::Zero(2*m_flow.VelocitySize()+m_flow.PressureSize());
        for(int component=0;component<3;++component){auto function=[&](double x,double z){return initial.Evaluate(component<2?component:4,x,z);};auto value=m_flow.Project(function,component);
            if(component<2)state.flow.segment(component*m_flow.VelocitySize(),m_flow.VelocitySize())=value;else state.flow.tail(m_flow.PressureSize())=value;}
        state.k=m_k.Project([&](double x,double z){return initial.Evaluate(2,x,z);});
        state.epsilon=m_e.Project([&](double x,double z){return initial.Evaluate(3,x,z);});
    }
    auto history=FreshHistory(state);
    if(!allowCandidate)Require(SteadyResidual(history.current).Maximum()<=m_options.nonlinearTolerance,
        "imported reference is not stationary for this native recipe; prepare a current reference first");
    m_initialPacketSha=initial.Metadata().packetSha256;m_initialSourceStateSha=initial.Metadata().stateSha256;
    m_initialProjected=!same;m_initialCandidate=!initial.Metadata().qualifiedReference;return history;
}

Vec Coupling::PackFree(const CoupledState&s)const
{
    const int nf=m_flow.Size(),nk=m_k.FreeSize(),ne=m_e.FreeSize();Vec result(nf+nk+ne);
    result.head(nf)=m_flow.Restrict(s.flow);result.segment(nf,nk)=m_k.Restrict(s.k);result.tail(ne)=m_e.Restrict(s.epsilon);return result;
}
CoupledState Coupling::FromFree(const Vec&x,double time)const
{
    const int nf=m_flow.Size(),nk=m_k.FreeSize(),ne=m_e.FreeSize();Require(x.size()==nf+nk+ne&&x.allFinite(),"invalid coupled free unknown inventory");
    return DirichletState({m_flow.Expand(x.head(nf),FlowLift(time)),m_k.Expand(x.segment(nf,nk),Vec::Zero(m_k.Size())),m_e.Expand(x.tail(ne),Vec::Zero(m_e.Size()))},time);
}
Vec Coupling::UnknownScales()const
{
    const int nf=m_flow.Size(),nk=m_k.FreeSize(),ne=m_e.FreeSize();Vec scale=Vec::Ones(nf+nk+ne);
    scale.segment(m_flow.PressureOffset(),m_flow.PressureSize()).setConstant(m_reference->Parameters().tau0);
    for(int i=0;i<nk;++i){int id=m_k.FreeIds()[i];scale[nf+i]=1/(m_k.ResidualRowScale()[i]*std::sqrt(m_k.Mass().coeff(id,id)));}
    for(int i=0;i<ne;++i){int id=m_e.FreeIds()[i];scale[nf+nk+i]=1/(m_e.ResidualRowScale()[i]*std::sqrt(m_e.Mass().coeff(id,id)));}
    return scale;
}
CoupledHistory Coupling::SteadyCandidateHistory(const Vec&unknowns)const{return FreshHistory(FromFree(unknowns,0));}

Coupling::Data Coupling::Evaluate(const CoupledState&s,double time,bool principal,CoupledTiming&timing)const
{
    const double evaluateStart=Clock(),attributedBefore=timing.Sum();Data d;double start=Clock();const auto&p=m_reference->Parameters();
    auto kb=m_k.PolynomialBounds(s.k),eb=m_e.PolynomialBounds(s.epsilon);
    if(kb[0]<=0||eb[0]<=0){std::ostringstream error;error<<std::setprecision(17)<<"transported scalar positivity is not certified between quadrature points: k=["<<kb[0]<<','<<kb[1]<<"], epsilon=["<<eb[0]<<','<<eb[1]<<']';throw NativeInadmissibleState(error.str());}
    d.minK=kb[0];d.minE=std::min(eb[0],m_reference->Amplitude(p.KMinimum)/m_geometry.H);
    // The interface polynomial must remain in the qualified reference table
    // over its complete edge, rather than only at a chosen node set.
    auto km=m_fields[2]->GetLocalToGlobalMap();
    for(int ei=0;ei<m_fields[2]->GetExpSize();++ei){const auto&e=Owner(m_fields[2]->GetExp(ei)->GetGeom()->GetGlobalID());if(std::abs(e.rl-m_geometry.H)>1e-12)continue;
        std::array<double,4>c{};for(int j=0;j<4;++j){int l=m_fields[2]->GetCoeff_Offset(ei)+j;c[j]=s.k[km->GetLocalToGlobalMap(l)]*km->GetLocalToGlobalSign(l);}
        auto bound=CubicBounds(c);RequireState(bound[0]>=p.KMinimum&&bound[1]<=p.KMaximum,"current interface K polynomial left the qualified reference interval");
    }
    timing.guards+=Clock()-start;start=Clock();
    d.kRules=m_rules.Cutoff([&](const Element&e,double x){return m_reference->Cutoff(InterfaceK(e,x,s));});
    timing.wall+=Clock()-start;
    std::vector<double>coefficient,kCoefficient,eCoefficient;const int n=m_flow.VelocitySize();
    auto samples=[&](const NativeRule&r,bool flowPart,Vec&mu,Vec&Dk,Vec&De,Vec&sourceK,Vec&sourceE,Vec&reactionE,std::vector<Tensor>&rem,Vec&bulk,Vec&shear,std::vector<Tensor>&viscous){
        Vec u=Local(s.flow,r.u),v=Local(s.flow,r.u,n),kx=m_k.Evaluate(r,s.k),eps;
        Vec ux=r.Dx*u,uz=r.Dz*u,vx=r.Dx*v,vz=r.Dz*v;
        const bool outer=m_e.Owns(r.geometry);
        if(outer)eps=m_e.Evaluate(r,s.epsilon);
        Vec interface;
        if(!outer){interface=Vec(r.x.size());for(int q=0;q<interface.size();++q)interface[q]=InterfaceK(Owner(r.geometry),r.x[q],s);}
        const int nq=r.w.size();mu=Vec(nq);Dk=Vec(nq);De=Vec(nq);sourceK=Vec(nq);sourceE=Vec(nq);reactionE=Vec(nq);rem.resize(nq);bulk=Vec(nq);shear=Vec(nq);viscous.resize(nq);
        if(flowPart)d.velocityInputs.emplace_back(nq);
        for(int q=0;q<nq;++q){const double radius=r.z[q]-m_geometry.d0,k=kx[q];
            RequireState(std::isfinite(k)&&k>0,"current transported k is not strictly positive");
            double E=outer?eps[q]:m_reference->Ew(interface[q],radius),aux=outer?0:m_reference->Aux(interface[q],radius);
            RequireState(std::isfinite(E)&&E>0,"current transported/dependent epsilon is not strictly positive");
            Tensor G{};G[0][0]=ux[q];G[0][2]=uz[q];G[2][0]=vx[q];G[2][2]=vz[q];
            if(flowPart)d.velocityInputs.back()[q]={k,E,aux,radius,G};
            double began=Clock();auto a=AdaptAmended(k,E,G,{0,0,1},p.nu,aux,radius,m_geometry.H,m_options.backend,&m_constitutiveStats);timing.constitutive+=Clock()-began;
            Require(a.mu>0&&std::isfinite(a.production),"invalid current complete constitutive data");
            mu[q]=a.mu;rem[q]=a.remainder;viscous[q]=a.V;
            shear[q]=m_options.stressSplit=="PositiveShearV1"?(1-a.chi)*aux:0.;
            bulk[q]=a.mu-shear[q];
            // B=b M+a[(PMn) tensor n+n tensor (PMn)]. On the flat
            // supported geometry only the xz/zx shear receives a. Keep the
            // COMPLETE V-B remainder, including weak-divergence terms.
            if(shear[q]!=0)for(int i:{0,2})for(int j:{0,2})
                if(i==j)rem[q][i][j]+=shear[q]*(G[i][j]+G[j][i]);
            Dk[q]=p.nu+p.Cd*k*k/(p.sigmaK*E);
            De[q]=outer?p.nu+p.Cd*k*k/(p.sigmaE*E):0.;
            sourceK[q]=a.production-E;sourceE[q]=outer?1.44*a.production*E/k:0.;reactionE[q]=outer?1.82*E/k:0.;
            d.minK=std::min(d.minK,k);d.minE=std::min(d.minE,E);
            if(flowPart){coefficient.push_back(bulk[q]);coefficient.push_back(shear[q]);kCoefficient.push_back(Dk[q]);if(outer){eCoefficient.push_back(De[q]);eCoefficient.push_back(reactionE[q]);}}
            // Covariance is checked at every current assembly point. Principal
            // checks include the continuum of supported planar directions.
            began=Clock();auto guard=CheckPlanarState(k,E,G,p.nu,aux,radius,m_geometry.H,principal&&flowPart,&a,m_options.backend,&m_constitutiveStats);
            d.minCovariance=std::min(d.minCovariance,guard.covarianceOverK);
            if(principal&&flowPart)d.minPrincipal=std::min(d.minPrincipal,guard.principalOverMu);
            if(guard.covarianceOverK<=1e-10||(principal&&flowPart&&guard.principalOverMu<=1e-8)){
                std::ostringstream error;error<<std::setprecision(17)<<"current constitutive admissibility failed: geometry="<<r.geometry<<", x="<<r.x[q]<<", r="<<radius<<", k="<<k<<", E="<<E<<", auxiliary="<<aux<<", chi="<<a.chi<<", covariance/k="<<guard.covarianceOverK<<", principal/mu="<<guard.principalOverMu<<", Gxx="<<G[0][0]<<", Gxz="<<G[0][2]<<", Gzx="<<G[2][0]<<", Gzz="<<G[2][2];throw NativeInadmissibleState(error.str());
            }
            timing.guards+=Clock()-began;
        }
    };
    start=Clock();const double pointTimeBefore=timing.constitutive+timing.guards;
    auto accountSamples=[&]{timing.assembly+=Clock()-start-
        (timing.constitutive+timing.guards-pointTimeBefore);};
    try{
    std::vector<Vec>staticDk,staticSourceK;
    for(const auto&rp:m_rules.Static()){
        Vec mu,Dk,De,sk,se,re,bulk,shear;std::vector<Tensor>rem,viscous;samples(*rp,true,mu,Dk,De,sk,se,re,rem,bulk,shear,viscous);d.mu.push_back(std::move(mu));d.remainder.push_back(std::move(rem));d.bulk.push_back(std::move(bulk));d.shear.push_back(std::move(shear));d.viscousStress.push_back(std::move(viscous));
        // epsilon is supported only on its original outer element inventory.
        d.De.push_back(std::move(De));d.sourceE.push_back(std::move(se));d.reactionE.push_back(std::move(re));
        staticDk.push_back(std::move(Dk));staticSourceK.push_back(std::move(sk));
    }
    for(int i=0;i<int(d.kRules.size());++i){const auto&rp=d.kRules[i];
        if(rp.get()==m_rules.Static()[i].get()){d.Dk.push_back(std::move(staticDk[i]));d.sourceK.push_back(std::move(staticSourceK[i]));}
        else{Vec mu,Dk,De,sk,se,re,bulk,shear;std::vector<Tensor>rem,viscous;samples(*rp,false,mu,Dk,De,sk,se,re,rem,bulk,shear,viscous);d.Dk.push_back(std::move(Dk));d.sourceK.push_back(std::move(sk));}
    }
    d.coefficients=Vec(coefficient.size());for(int i=0;i<d.coefficients.size();++i)d.coefficients[i]=coefficient[i];
    d.kCoefficients=Vec(kCoefficient.size());for(int i=0;i<d.kCoefficients.size();++i)d.kCoefficients[i]=kCoefficient[i];
    d.eCoefficients=Vec(eCoefficient.size());for(int i=0;i<d.eCoefficients.size();++i)d.eCoefficients[i]=eCoefficient[i];
    }catch(...){accountSamples();throw;}
    // Native interpolation, gradients and complete source/coefficient data
    // assembly, excluding the separately timed WJ kernels and point guards.
    accountSamples();
    // Compare the actual essential coefficients with the freshly rebuilt
    // nonlinear trace projection. Pointwise projection error is a separate
    // spatial error, never an algebraic boundary condition to overconstrain.
    start=Clock();auto refreshed=DirichletState(s,time);
    const double scale=m_reference->Amplitude(p.KIn)/m_geometry.H;
    d.traceError=std::max({(m_e.Lift(s.epsilon)-m_e.Lift(refreshed.epsilon)).cwiseAbs().maxCoeff()/scale,
        (m_k.Lift(s.k)-m_k.Lift(refreshed.k)).cwiseAbs().maxCoeff()/p.KIn,
        (m_flow.Lift(s.flow)-m_flow.Lift(refreshed.flow)).cwiseAbs().maxCoeff()});
    timing.wall+=Clock()-start;
    if(!m_mms){
        start=Clock();for(const auto&[index,side]:m_boundaryFaces)if(side==1){const auto&e=m_flow.Elements()[index];std::array<double,4>c{};for(int j=0;j<4;++j)c[j]=s.flow[e.u[1+4*j]]*e.us[1+4*j];RequireState(CubicBounds(c)[0]>=-1e-12,"physical scalar outlet backflow is unsupported");}timing.guards+=Clock()-start;
    }
    // Native basis actions, field/gradient sampling and coefficient-vector
    // construction are equation assembly work. Subtract the disjoint point,
    // wall and guard timers already recorded above before assigning it.
    timing.assembly+=(Clock()-evaluateStart)-(timing.Sum()-attributedBefore);
    return d;
}

Sparse Coupling::VelocityTangent(const Data&d,CoupledTiming&timing)const
{
    Require(d.velocityInputs.size()==m_rules.Static().size(),"steady tangent current input inventory differs");
    const double start=Clock();std::vector<OwnedPlanarVelocityTangents>samples; samples.reserve(m_rules.Static().size());
    for(int i=0;i<int(m_rules.Static().size());++i){
        const auto&r=*m_rules.Static()[i];Require(d.velocityInputs[i].size()==std::size_t(r.w.size()),"steady tangent current point inventory differs");
        OwnedPlanarVelocityTangents cell;cell.geometry=r.geometry;cell.points.reserve(r.w.size());
        for(const auto&input:d.velocityInputs[i])cell.points.push_back(EvaluatePlanarVelocityTangent(input.k,input.E,input.gradient,
            m_reference->Parameters().nu,input.aux,input.r,m_geometry.H,m_options.backend,&m_constitutiveStats));
        samples.push_back(std::move(cell));
    }
    timing.constitutive+=Clock()-start;const double assembly=Clock();
    auto matrix=AssemblePlanarVelocityTangent(m_rules.Static(),m_flow.VelocitySize(),samples);
    timing.assembly+=Clock()-assembly;return matrix;
}

Sparse Coupling::Strain(const Data&d)const
{
    std::vector<Triplet>t;t.reserve(m_rules.Static().size()*1024);const int n=m_flow.VelocitySize();
    for(int e=0;e<int(m_rules.Static().size());++e){const auto&r=*m_rules.Static()[e];Vec w=r.w.cwiseProduct(d.bulk[e]),wa=r.w.cwiseProduct(d.shear[e]);std::array<Mat,2>g{r.Dx,r.Dz};Mat lap=r.Dx.transpose()*w.asDiagonal()*r.Dx+r.Dz.transpose()*w.asDiagonal()*r.Dz;
        for(int i=0;i<2;++i)for(int j=0;j<2;++j){Mat k=g[j].transpose()*w.asDiagonal()*g[i];if(i==j)k+=lap;
            k.noalias()+=g[1-i].transpose()*wa.asDiagonal()*g[1-j];
            for(int a=0;a<16;++a)for(int b=0;b<16;++b)t.emplace_back(i*n+r.u[a],j*n+r.u[b],k(a,b));}}
    Sparse result(2*n,2*n);result.setFromTriplets(t.begin(),t.end());return result;
}
Vec Coupling::RemainderLoad(const Data&d)const
{
    const int n=m_flow.VelocitySize();Vec load=Vec::Zero(2*n+m_flow.PressureSize());
    for(int e=0;e<int(m_rules.Static().size());++e){const auto&r=*m_rules.Static()[e];std::array<Vec,4>stress;for(auto&v:stress)v=Vec(r.w.size());
        for(int q=0;q<r.w.size();++q)for(int i=0;i<2;++i)for(int j=0;j<2;++j)stress[2*i+j][q]=r.w[q]*d.remainder[e][q][i?2:0][j?2:0];
        for(int i=0;i<2;++i)Add(load,r.u,-r.Dx.transpose()*stress[2*i]-r.Dz.transpose()*stress[2*i+1],i*n);}
    return load;
}

Vec Coupling::MomentumAdvection(const CoupledState&s)const
{
    const int n=m_flow.VelocitySize();Vec load=Vec::Zero(2*n+m_flow.PressureSize());
    for(const auto&rp:m_rules.Static()){
        const auto&r=*rp;Vec u=Local(s.flow,r.u),v=Local(s.flow,r.u,n),up=r.B*u,vp=r.B*v;
        for(int component=0;component<2;++component){const Vec&c=component==0?u:v;Vec a=up.cwiseProduct(r.Dx*c)+vp.cwiseProduct(r.Dz*c);
            // EX2 applies to the complete declared tangential convective
            // roughness acceleration, as well as the domain convection.
            if(component==0&&Owner(r.geometry).rl<m_geometry.h){const auto&e=Owner(r.geometry);
                for(int q=0;q<a.size();++q)if(r.z[q]-m_geometry.d0<m_geometry.h){auto h=InterfaceVelocity(e,r.x[q],s);a[q]-=2*h[0]*h[1];}}
            Add(load,r.u,r.B.transpose()*r.w.asDiagonal()*a,component*n);
        }
    }
    return load;
}
CoupledState Coupling::Advection(const CoupledState&s,const Rules&kRules)const
{return {MomentumAdvection(s),m_k.Advection(kRules,s.flow,s.k),m_e.Advection(m_rules.Static(),s.flow,s.epsilon)};}

Sparse Coupling::MomentumTransportJacobian(const CoupledState&s)const
{
    const int n=m_flow.VelocitySize();std::vector<Triplet>triplets;
    for(const auto&rp:m_rules.Static()){
        const auto&r=*rp;const Vec u=Local(s.flow,r.u),v=Local(s.flow,r.u,n),up=r.B*u,vp=r.B*v;
        const std::array<Vec,2>gx{Vec(r.Dx*u),Vec(r.Dx*v)},gz{Vec(r.Dz*u),Vec(r.Dz*v)};
        const Mat adv=up.asDiagonal()*r.Dx+vp.asDiagonal()*r.Dz;
        for(int i=0;i<2;++i)for(int j=0;j<2;++j){Mat local=r.B.transpose()*r.w.asDiagonal()*
                ((j==0?gx[i]:gz[i]).asDiagonal()*r.B);
            if(i==j)local+=r.B.transpose()*r.w.asDiagonal()*adv;
            for(int a=0;a<16;++a)for(int b=0;b<16;++b)triplets.emplace_back(i*n+r.u[a],j*n+r.u[b],local(a,b));
        }
        if(Owner(r.geometry).rl<m_geometry.h){auto h=HRule(r);const Vec hu=Local(s.flow,h.u),uh=h.B*hu,uxh=h.Dx*hu;
            Vec weights=r.w;for(int q=0;q<weights.size();++q)if(r.z[q]-m_geometry.d0>=m_geometry.h)weights[q]=0;
            Mat sampled=-2*r.B.transpose()*weights.asDiagonal()*(uxh.asDiagonal()*h.B+uh.asDiagonal()*h.Dx);
            for(int a=0;a<16;++a)for(int b=0;b<16;++b)triplets.emplace_back(r.u[a],h.u[b],sampled(a,b));
        }
    }
    Sparse result(2*n,2*n);result.setFromTriplets(triplets.begin(),triplets.end());return result;
}

Coupling::Loads Coupling::Forcing(double time,CoupledTiming&timing,bool freezeManufactured)const
{
    if(!m_mms&&m_hasPhysicalForcing)return m_physicalForcing;
    double start=Clock();const int n=m_flow.VelocitySize();Loads load{Vec::Zero(2*n+m_flow.PressureSize()),Vec::Zero(m_k.Size()),Vec::Zero(m_e.Size())};
    m_bulkBodyAction=Vec::Zero(2*n);m_sampledBodyAction=Vec::Zero(2*n);
    if(m_mms){
        std::vector<Vec>esource;
        for(const auto&rp:m_rules.Static()){
            const auto&r=*rp;Vec fx(r.w.size()),fz(fx.size()),fe(fx.size()),bx(fx.size()),bz(fx.size()),sampled=Vec::Zero(fx.size());
            for(int q=0;q<fx.size();++q){auto f=m_mms->Forcing(r.x[q],r.z[q],time);bx[q]=f.body[0];bz[q]=f.body[2];
                if(r.z[q]-m_geometry.d0<m_geometry.h)sampled[q]=-m_mms->BodyAtH(r.x[q],time)[0];
                fx[q]=f.momentum[0]+bx[q]+sampled[q];fz[q]=f.momentum[2]+bz[q];fe[q]=f.outer?f.epsilon:0.;
                if(freezeManufactured){auto fields=m_mms->Fields(r.x[q],r.z[q],time);fx[q]-=fields.u.dt;fz[q]-=fields.w.dt;if(f.outer)fe[q]-=m_mms->OuterEpsilon(r.x[q],r.z[q],time).dt;}}
            Add(load.flow,r.u,r.B.transpose()*r.w.asDiagonal()*fx);Add(load.flow,r.u,r.B.transpose()*r.w.asDiagonal()*fz,n);esource.push_back(std::move(fe));
            Add(m_bulkBodyAction,r.u,r.B.transpose()*r.w.asDiagonal()*bx);Add(m_bulkBodyAction,r.u,r.B.transpose()*r.w.asDiagonal()*bz,n);
            Add(m_sampledBodyAction,r.u,r.B.transpose()*r.w.asDiagonal()*sampled);
        }
        load.epsilon=m_e.Source(m_rules.Static(),esource);
        // The continuous forcing has its OWN exact cutoff surface. It is
        // integrated independently of the current numerical cutoff surface.
        auto exactRules=m_rules.Cutoff([&](const Element&,double x){return m_reference->Cutoff(m_mms->InterfaceK(x,time).value);});
        std::vector<Vec>ksource;
        for(const auto&rp:exactRules){Vec f(rp->w.size());for(int q=0;q<f.size();++q){f[q]=m_mms->Forcing(rp->x[q],rp->z[q],time).k;if(freezeManufactured)f[q]-=m_mms->Fields(rp->x[q],rp->z[q],time).k.dt;}ksource.push_back(std::move(f));}
        load.k=m_k.Source(exactRules,ksource);
        load.k+=m_k.NaturalFlux([&](double x,double z,const Vector&normal){return m_mms->KOutwardFlux(x,z,time,normal);},m_rules,
            [&](const Element&,double x){return m_reference->Cutoff(m_mms->InterfaceK(NativeX(x),time).value);});
        load.epsilon+=m_e.NaturalFlux([&](double x,double z,const Vector&normal){
            // H is Dirichlet only: its reaction is diagnosed from the solve.
            return normal[2]<-.5?0.:m_mms->EpsilonOutwardFlux(x,z,time,normal);
        },m_rules);
    }else{
        // The application has no manufactured/balance forcing. Body force
        // is exactly zero for this pinned neutral, flat physical case.
        const auto&p=m_reference->Parameters();
        for(const auto&[index,side]:m_boundaryFaces)if(side==1||side==3){auto r=m_rules.Face(m_flow.Elements()[index],side,20);Vec tx(r.w.size()),tz(tx.size());
            for(int q=0;q<tx.size();++q){if(side==1){auto t=m_reference->OutletTraction(r.z[q]-m_geometry.d0);tx[q]=t[0];tz[q]=t[1];}else{tx[q]=p.tau0;tz[q]=0;}}
            Add(load.flow,r.u,r.B.transpose()*r.w.asDiagonal()*tx);Add(load.flow,r.u,r.B.transpose()*r.w.asDiagonal()*tz,n);}
        load.epsilon+=m_e.NaturalFlux([&](double,double,const Vector&normal){return normal[2]>.5?m_reference->TopEpsilonOutwardFlux():0.;},m_rules);
    }
    timing.forcing+=Clock()-start;if(!m_mms){m_physicalForcing=load;m_hasPhysicalForcing=true;}return load;
}
Vec Coupling::KVolumeLoad(const CoupledState&,const Data&d)const{return m_k.Source(d.kRules,d.sourceK);}
Vec Coupling::EVolumeLoad(const CoupledState&,const Data&d)const{return m_e.Source(m_rules.Static(),d.sourceE);}
double Coupling::Drift(const Vec&current,const Vec&previous)const
{
    if(previous.size()!=current.size())return 1.;
    // The shear coefficient is exactly zero in the outer region. A physical
    // molecular-viscosity reference makes this refresh-only metric defined;
    // no coefficient in any PDE/operator is floored or changed.
    return ((current-previous).cwiseAbs().array()/previous.cwiseAbs().array().max(m_reference->Parameters().nu)).maxCoeff();
}

Coupling::Loads Coupling::CurrentActions(const CoupledState&s,const Data&d)const
{
    const int n=m_flow.VelocitySize();Loads action{Vec::Zero(2*n+m_flow.PressureSize()),
        m_k.DiffusionAction(d.kRules,d.Dk,s.k)-KVolumeLoad(s,d),
        m_e.DiffusionAction(m_rules.Static(),d.De,s.epsilon,d.reactionE)-EVolumeLoad(s,d)};
    for(int e=0;e<int(m_rules.Static().size());++e){const auto&r=*m_rules.Static()[e];std::array<Vec,4>stress;for(auto&v:stress)v=Vec(r.w.size());
        for(int q=0;q<r.w.size();++q)for(int i=0;i<2;++i)for(int j=0;j<2;++j)
            stress[2*i+j][q]=r.w[q]*d.viscousStress[e][q][i?2:0][j?2:0];
        for(int i=0;i<2;++i)Add(action.flow,r.u,r.Dx.transpose()*stress[2*i]+r.Dz.transpose()*stress[2*i+1],i*n);
    }
    action.flow.head(2*n)+=(m_flow.Core().Gradient()-m_flow.Core().AWF())*s.flow.tail(m_flow.PressureSize())-m_flow.KnownScalarKLoad(s.k);
    action.flow.tail(m_flow.PressureSize())=m_flow.Core().Divergence()*s.flow.head(2*n);
    if(m_options.checkCurrentActions){
        const Sparse K=Strain(d),Kk=m_k.Diffusion(d.kRules,d.Dk),Ke=m_e.Diffusion(m_rules.Static(),d.De,d.reactionE);
        Loads assembled{m_flow.Core().ApplyFull(0,K,1,s.flow)-RemainderLoad(d),Kk*s.k-KVolumeLoad(s,d),Ke*s.epsilon-EVolumeLoad(s,d)};
        assembled.flow.head(2*n)-=m_flow.KnownScalarKLoad(s.k);
        const auto&p=m_reference->Parameters();
        int component=0;++m_actionChecks;
        for(const auto&pair:{std::make_pair(&action.flow,&assembled.flow),std::make_pair(&action.k,&assembled.k),std::make_pair(&action.epsilon,&assembled.epsilon)}){
            const double difference=(*pair.first-*pair.second).cwiseAbs().maxCoeff();
            // Fixed physical scales, independent of alpha or current RHS.
            const double scale=pair.first==&action.flow?1.:pair.first==&action.k?p.KIn:m_reference->Amplitude(p.KIn)/m_geometry.H;
            m_actionDifference[component]=std::max(m_actionDifference[component],difference/scale);++component;
            Require(difference/scale<=1e-10,"direct/assembled complete current all-row action disagreement");
        }
    }
    return action;
}

CoupledResidual Coupling::Residual(const CoupledState&s,const Data&d,const Loads&external,double alpha,double time,CoupledTiming&timing,const CoupledHistory*history)const
{
    const double start=Clock();(void)time;auto action=CurrentActions(s,d);
    if(alpha>0){
        Require(history!=nullptr,"transient residual requires complete accepted history");
        // Form BDF differences BEFORE mass action, including essential and
        // dependent-Gamma history. This avoids cancellation of large mass RHSs.
        CoupledState delta{alpha*(s.flow-history->current.flow),alpha*(s.k-history->current.k),alpha*(s.epsilon-history->current.epsilon)};
        if(history->step>0){const double a=.5/history->dt;
            delta.flow+=a*(history->previous.flow-history->current.flow);
            delta.k+=a*(history->previous.k-history->current.k);
            delta.epsilon+=a*(history->previous.epsilon-history->current.epsilon);}
        action.flow+=m_flow.Core().MassAction(delta.flow);action.k+=m_k.Mass()*delta.k;action.epsilon+=m_e.Mass()*delta.epsilon;
    }
    action.flow-=external.flow;action.k-=external.k;action.epsilon-=external.epsilon;
    auto raw=m_flow.Restrict(action.flow);if(m_flow.Gauge())raw[raw.size()-1]=m_flow.PressureMean().dot(s.flow.tail(m_flow.PressureSize()));
    const auto flow=m_flow.RawResidualNorm(raw);
    const double kr=m_k.FreeResidualNorm(m_k.Restrict(action.k)),er=m_e.FreeResidualNorm(m_e.Restrict(action.epsilon));
    timing.residual+=Clock()-start;return {flow.momentum,flow.continuity,flow.gauge,kr,er,d.traceError};
}

CoupledStepInfo Coupling::Advance(CoupledHistory&h,bool injectRetry)
{
    if(m_options.advancement=="CoupledJFNKV1")return AdvanceNewton(h,injectRetry);
    const double begin=Clock();CheckHistory(h);m_constitutiveStats={};m_lastInputHistory=h;m_hasInputHistory=true;CoupledStepInfo info;const double time=(h.step+1)*h.dt,alpha=(h.step==0?1.:1.5)/h.dt;
    m_lastAttempt={};auto saveAttempt=[&]{info.seconds=Clock()-begin;info.constitutiveStats=m_constitutiveStats;info.timing.other=info.seconds-info.timing.Sum();m_lastAttempt=info;};
    auto history=Combination(h.current,h.step==0?1/h.dt:2/h.dt,h.previous,h.step==0?0:-.5/h.dt);
    auto adv=Combination(h.advection,h.step==0?1:2,h.previousAdvection,h.step==0?0:-1);
    Loads external=Forcing(time,info.timing);external.flow-=adv.flow;external.k-=adv.k;external.epsilon-=adv.epsilon;
    Loads base=external;base.flow+=m_flow.Core().MassAction(history.flow);
    base.k+=m_k.Mass()*history.k;base.epsilon+=m_e.Mass()*history.epsilon;
    double boundaryStart=Clock();CoupledState s=DirichletState(h.current,time);info.timing.wall+=Clock()-boundaryStart;
    const Vec accelerationScale=m_options.andersonDepth?UnknownScales():Vec();
    auto alphaDrift=[&](double saved){return saved>0?std::abs(alpha/saved-1):1.;};
    // Every retry starts from the accepted history at the SAME dt. No failed
    // provisional state, EX2 history or scientific time can be committed.
    for(int attempt=0;attempt<2;++attempt){
        if(attempt){ResetPreconditioners();boundaryStart=Clock();s=DirichletState(h.current,time);info.timing.wall+=Clock()-boundaryStart;++info.retries;}
        bool complete=false,recoverable=false;
        std::vector<Vec>andersonF,andersonG;
        auto data=Evaluate(s,time,false,info.timing);
        try{for(int iteration=0;iteration<m_options.maximumIterations;++iteration){
            ++info.nonlinearIterations;
            const double iterationStart=Clock();const Vec input=m_options.andersonDepth?PackFree(s):Vec();
            double start=Clock();
            const double drift=std::max(Drift(data.coefficients,m_preconditionerCoefficients),alphaDrift(m_flowFactorAlpha));info.coefficientDrift=std::max(info.coefficientDrift,drift);
            Sparse K=Strain(data),A=m_flow.System(alpha,K,1.);Vec flowLoad=base.flow+RemainderLoad(data);flowLoad.head(2*m_flow.VelocitySize())+=m_flow.KnownScalarKLoad(s.k);
            Vec lift=m_flow.Lift(s.flow),rhs=m_flow.LiftedRhs(alpha,K,1.,flowLoad,lift),free=m_flow.Restrict(s.flow);info.timing.assembly+=Clock()-start;
            if(injectRetry&&attempt==0&&iteration==0)throw NativeLinearFailure("qualification injected recoverable Krylov failure");
            auto fi=m_flow.Solve(A,rhs,free,drift,m_options.linearTolerance);s.flow=m_flow.Expand(free,lift);info.flowIterations+=fi.iterations;info.rebuilds+=fi.rebuilds;info.timing.flow+=fi.solveSeconds+fi.residualSeconds;info.timing.rebuild+=fi.factorSeconds;
            if(fi.rebuilds){m_preconditionerCoefficients=data.coefficients;m_flowFactorAlpha=alpha;}
            // Refresh P, E_w, current positive scalar diffusivities and all
            // stresses after the flow update, BEFORE the k solve.
            data=Evaluate(s,time,false,info.timing);start=Clock();Sparse Kk=m_k.Diffusion(data.kRules,data.Dk);Vec klift=m_k.Lift(s.k),krhs=m_k.LiftedRhs(alpha,Kk,base.k+KVolumeLoad(s,data),klift),kf=m_k.Restrict(s.k);auto Ak=m_k.System(alpha,Kk);info.timing.assembly+=Clock()-start;
            const double kd=std::max(Drift(data.kCoefficients,m_kPreconditionerCoefficients),alphaDrift(m_kFactorAlpha));info.coefficientDrift=std::max(info.coefficientDrift,kd);
            auto ki=m_k.Solve(Ak,krhs,kf,kd,m_options.linearTolerance);s.k=m_k.Expand(kf,klift);info.kIterations+=ki.iterations;info.rebuilds+=ki.rebuilds;info.timing.scalar+=ki.solveSeconds+ki.residualSeconds;info.timing.rebuild+=ki.factorSeconds;
            if(ki.rebuilds){m_kPreconditionerCoefficients=data.kCoefficients;m_kFactorAlpha=alpha;}
            // K(H) has changed: replace the dependent outer epsilon trace and
            // recompute the entire closure before solving outer epsilon.
            start=Clock();s=DirichletState(s,time);info.timing.wall+=Clock()-start;data=Evaluate(s,time,false,info.timing);start=Clock();Sparse Ke=m_e.Diffusion(m_rules.Static(),data.De,data.reactionE);Vec elift=m_e.Lift(s.epsilon),erhs=m_e.LiftedRhs(alpha,Ke,base.epsilon+EVolumeLoad(s,data),elift),ef=m_e.Restrict(s.epsilon);auto Ae=m_e.System(alpha,Ke);info.timing.assembly+=Clock()-start;
            const double ed=std::max(Drift(data.eCoefficients,m_ePreconditionerCoefficients),alphaDrift(m_eFactorAlpha));info.coefficientDrift=std::max(info.coefficientDrift,ed);
            auto ei=m_e.Solve(Ae,erhs,ef,ed,m_options.linearTolerance);s.epsilon=m_e.Expand(ef,elift);info.epsilonIterations+=ei.iterations;info.rebuilds+=ei.rebuilds;info.timing.scalar+=ei.solveSeconds+ei.residualSeconds;info.timing.rebuild+=ei.factorSeconds;
            if(ei.rebuilds){m_ePreconditionerCoefficients=data.eCoefficients;m_eFactorAlpha=alpha;}
            // Acceptance always tests the complete CURRENT coupled equations.
            data=Evaluate(s,time,false,info.timing);info.residual=Residual(s,data,external,alpha,time,info.timing,&h);
            bool accelerated=false;
            if(m_options.andersonDepth){
                start=Clock();const Vec fixed=PackFree(s),g=(fixed-input).cwiseQuotient(accelerationScale);
                andersonF.push_back(fixed);andersonG.push_back(g);
                while(int(andersonF.size())>m_options.andersonDepth+1){andersonF.erase(andersonF.begin());andersonG.erase(andersonG.begin());}
                const int columns=int(andersonF.size())-1;Vec proposal;
                if(columns>0&&info.residual.Maximum()>m_options.nonlinearTolerance){
                    Mat dg(g.size(),columns),df(g.size(),columns);
                    for(int j=0;j<columns;++j){dg.col(j)=andersonG[j+1]-andersonG[j];df.col(j)=andersonF[j+1]-andersonF[j];}
                    Eigen::ColPivHouseholderQR<Mat>qr(dg);qr.setThreshold(1e-10);
                    const Vec coefficients=qr.solve(g);
                    // A bounded numerical extrapolation, never a physical
                    // floor. Every usable proposal is checked afresh below.
                    if(qr.rank()>0&&coefficients.allFinite()&&coefficients.norm()<=1e3){proposal=fixed-df*coefficients;if(!proposal.allFinite())proposal.resize(0);}
                }
                info.timing.nonlinear+=Clock()-start;
                if(proposal.size()){
                    ++info.accelerationProposals;
                    for(double fraction:{1.,.5,.25,.125}){
                        try{
                            start=Clock();auto candidate=FromFree(fixed+fraction*(proposal-fixed),time);info.timing.wall+=Clock()-start;
                            auto fresh=Evaluate(candidate,time,false,info.timing);auto residual=Residual(candidate,fresh,external,alpha,time,info.timing,&h);
                            if(residual.Maximum()<info.residual.Maximum()){
                                s=std::move(candidate);data=std::move(fresh);info.residual=residual;accelerated=true;++info.accelerationAccepted;break;
                            }
                        }catch(const NativeInadmissibleState&){++info.accelerationInadmissible;}
                    }
                }
            }
            if(m_options.recordNonlinearTrace)info.nonlinearTrace.push_back({info.residual,accelerated,Clock()-iterationStart});
            if(info.residual.Maximum()<=m_options.nonlinearTolerance){
                auto final=Evaluate(s,time,true,info.timing);info.residual=Residual(s,final,external,alpha,time,info.timing,&h);
                Require(info.residual.Maximum()<=m_options.nonlinearTolerance,"final refreshed coupled residual changed acceptance");
                info.minimumK=final.minK;info.minimumEpsilon=final.minE;info.minimumCovarianceOverK=final.minCovariance;info.minimumPrincipalOverMu=final.minPrincipal;
                const double advStart=Clock();auto nextAdvection=Advection(s,final.kRules);info.timing.assembly+=Clock()-advStart;
                h.previous=h.current;h.current=s;h.previousAdvection=h.advection;h.advection=std::move(nextAdvection);++h.step;h.time=time;
                complete=true;break;
            }
        }}catch(const NativeLinearFailure&e){info.retryReason=e.what();recoverable=true;if(attempt==1){saveAttempt();throw;}}
        catch(const NativeInadmissibleState&){saveAttempt();throw;}
        if(complete){saveAttempt();CheckHistory(h);return info;}
        if(!recoverable)break;
    }
    saveAttempt();std::ostringstream error;error<<std::setprecision(17)<<"current coupled iteration did not converge; accepted history preserved, residuals=["<<info.residual.momentum<<','<<info.residual.continuity<<','<<info.residual.gauge<<','<<info.residual.k<<','<<info.residual.epsilon<<','<<info.residual.trace<<"], nonlinear iterations="<<info.nonlinearIterations<<", linear iterations=["<<info.flowIterations<<','<<info.kIterations<<','<<info.epsilonIterations<<']';throw std::runtime_error(error.str());
}

CoupledResidual Coupling::SteadyResidual(const CoupledState&s,double time)
{
    return AuditSteadyState(s,time).residual;
}
CoupledStepInfo Coupling::AuditSteadyState(const CoupledState&s,double time)
{
    const double began=Clock();CoupledStepInfo info;auto d=Evaluate(s,time,true,info.timing);auto load=Forcing(time,info.timing,true);auto adv=Advection(s,d.kRules);
    load.flow-=adv.flow;load.k-=adv.k;load.epsilon-=adv.epsilon;info.residual=Residual(s,d,load,0,time,info.timing);
    info.minimumK=d.minK;info.minimumEpsilon=d.minE;info.minimumCovarianceOverK=d.minCovariance;info.minimumPrincipalOverMu=d.minPrincipal;
    info.seconds=Clock()-began;info.timing.other=info.seconds-info.timing.Sum();return info;
}

SteadyKrylovResult Coupling::SolveSteady(CoupledHistory&history,const SteadyKrylovOptions&options,
    CoupledStepInfo&final,double preconditionerMass)
{
    const double began=Clock();CheckHistory(history);m_lastInputHistory=history;m_hasInputHistory=true;
    Require(history.step==0&&history.time==0,"steady preparation requires a fresh step-zero state");
    Require(std::isfinite(preconditionerMass)&&preconditionerMass>=0,"invalid steady preconditioner mass");
    Require(options.acceptanceNorm>0&&options.acceptanceNorm<=m_options.nonlinearTolerance,
        "steady acceptance must meet the declared current coupled tolerance");
    ResetPreconditioners();struct ColdExit{Coupling&owner;~ColdExit(){owner.ResetPreconditioners();}}cold{*this};
    m_constitutiveStats={};final={};CoupledTiming&timing=final.timing;
    const int nf=m_flow.Size(),nk=m_k.FreeSize(),ne=m_e.FreeSize(),size=nf+nk+ne;
    const auto forcing=Forcing(0,timing,true);
    auto pack=[&](const CoupledState&s){return PackFree(s);};
    auto unpack=[&](const Vec&x){return FromFree(x,0);};
    auto norms=[&](const Vec&raw){Require(raw.size()==size&&raw.allFinite(),"invalid native steady residual inventory");
        auto f=m_flow.RawResidualNorm(raw.head(nf));return CoupledResidual{f.momentum,f.continuity,f.gauge,
            m_k.FreeResidualNorm(raw.segment(nf,nk)),m_e.FreeResidualNorm(raw.tail(ne)),0};};
    Vec validatedUnknowns;std::shared_ptr<const Data>validatedData;
    auto currentData=[&](const Vec&x,const CoupledState&s){
        if(validatedData&&x.size()==validatedUnknowns.size()&&(x.array()==validatedUnknowns.array()).all())return validatedData;
        return std::shared_ptr<const Data>(std::make_shared<Data>(Evaluate(s,0,false,timing)));
    };
    auto equations=[&](const Vec&x){
        auto s=unpack(x);auto owned=currentData(x,s);const auto&d=*owned;const double start=Clock();auto adv=Advection(s,d.kRules);auto action=CurrentActions(s,d);
        action.flow+=adv.flow-forcing.flow;action.k+=adv.k-forcing.k;action.epsilon+=adv.epsilon-forcing.epsilon;
        Vec raw(size);raw.head(nf)=m_flow.Restrict(action.flow);
        if(m_flow.Gauge()){
            raw.segment(m_flow.PressureOffset(),m_flow.PressureSize())+=m_flow.PressureMean()*x[nf-1];
            raw[nf-1]=m_flow.PressureMean().dot(s.flow.tail(m_flow.PressureSize()));
        }
        raw.segment(nf,nk)=m_k.Restrict(action.k);raw.tail(ne)=m_e.Restrict(action.epsilon);
        if(m_options.checkCurrentActions){
            const Sparse K=Strain(d),Kk=m_k.Diffusion(d.kRules,d.Dk),Ke=m_e.Diffusion(m_rules.Static(),d.De,d.reactionE);
            Vec flowLoad=forcing.flow-adv.flow+RemainderLoad(d);flowLoad.head(2*m_flow.VelocitySize())+=m_flow.KnownScalarKLoad(s.k);
            Vec assembled(size);assembled.head(nf)=m_flow.System(0,K,1)*x.head(nf)-m_flow.LiftedRhs(0,K,1,flowLoad,m_flow.Lift(s.flow));
            assembled.segment(nf,nk)=m_k.System(0,Kk)*x.segment(nf,nk)-m_k.LiftedRhs(0,Kk,forcing.k-adv.k+KVolumeLoad(s,d),m_k.Lift(s.k));
            assembled.tail(ne)=m_e.System(0,Ke)*x.tail(ne)-m_e.LiftedRhs(0,Ke,forcing.epsilon-adv.epsilon+EVolumeLoad(s,d),m_e.Lift(s.epsilon));
            Vec scale=Vec::Ones(size);scale.segment(nf,nk).setConstant(m_reference->Parameters().KIn);scale.tail(ne).setConstant(m_reference->Amplitude(m_reference->Parameters().KIn)/m_geometry.H);
            const double difference=(raw-assembled).cwiseQuotient(scale).cwiseAbs().maxCoeff();
            ++m_augmentedChecks;m_augmentedDifference=std::max(m_augmentedDifference,difference);
            m_checkedWeakDivergence=std::max(m_checkedWeakDivergence,(m_flow.Core().Divergence()*s.flow.head(2*m_flow.VelocitySize())).norm());
            if(m_flow.Gauge())m_checkedGaugeMultiplier=std::max(m_checkedGaugeMultiplier,std::abs(x[nf-1]));
            Require(difference<=1e-10,"direct/assembled augmented steady residual disagreement");
        }
        timing.residual+=Clock()-start;return raw;
    };
    Vec sx=UnknownScales(),sr(size);
    sr.head(nf)=m_flow.ResidualRowScale();sr.segment(nf,nk)=m_k.ResidualRowScale();sr.tail(ne)=m_e.ResidualRowScale();
    CoupledState pcState;Sparse epsilonPCFull;
    SteadyKrylovCallbacks callbacks;
    callbacks.ValidTrial=[&](const Vec&x){
        try{auto s=unpack(x);auto k=m_k.PolynomialBounds(s.k),e=m_e.PolynomialBounds(s.epsilon);
            if(k[0]<=0||e[0]<=0)return SteadyTrialValidity{false,"transported scalar polynomial is not strictly positive"};
            return SteadyTrialValidity{};
        }catch(const NativeInadmissibleState&error){return SteadyTrialValidity{false,error.what()};}
    };
    callbacks.ValidateNonlinearState=[&](const Vec&x){
        try{
            auto s=unpack(x);auto full=std::make_shared<Data>(Evaluate(s,0,true,timing));
            validatedUnknowns=x;validatedData=std::move(full);
        }catch(const NativeInadmissibleState&error){throw SteadyInadmissibleTrial(error.what());}
    };
    callbacks.Residual=[&](const Vec&x)->Vec{try{return equations(x);}catch(const NativeInadmissibleState&error){throw SteadyInadmissibleTrial(error.what());}};
    callbacks.TrueNorm=[&](const Vec&raw){return norms(raw).Maximum();};
    callbacks.PreparePreconditioner=[&](const Vec&x,const Vec&){
        pcState=unpack(x);auto owned=currentData(x,pcState);const auto&d=*owned;
        // Only the inverse approximation changes. Actual steady equations
        // and every FD action still use the complete current nonlinear V.
        Sparse flowPC=(m_options.steadyVelocityPC=="FullGradientV1"?VelocityTangent(d,timing):Strain(d))+MomentumTransportJacobian(pcState);
        Sparse kPC=m_k.Diffusion(d.kRules,d.Dk)+m_k.TransportMatrix(d.kRules,pcState.flow);
        Sparse ePC=m_e.Diffusion(m_rules.Static(),d.De,d.reactionE)+m_e.TransportMatrix(m_rules.Static(),pcState.flow);
        // This optional mass belongs solely to the inverse approximation.
        // The actual residual and every finite-difference J action use alpha=0.
        epsilonPCFull=ePC+preconditionerMass*m_e.Mass();
        m_flow.PrepareLinearPreconditioner(m_flow.System(preconditionerMass,flowPC,1));
        m_k.PrepareLinearPreconditioner(m_k.System(preconditionerMass,kPC));
        m_e.PrepareLinearPreconditioner(m_e.System(preconditionerMass,ePC));
    };
    callbacks.Precondition=[&](const Vec&rhs){
        Vec result(size);result.segment(nf,nk)=m_k.ApplyLinearPreconditioner(rhs.segment(nf,nk));
        Vec dk=m_k.Expand(result.segment(nf,nk),Vec::Zero(m_k.Size()));
        Vec cross=Vec::Zero(2*m_flow.VelocitySize()+m_flow.PressureSize());cross.head(2*m_flow.VelocitySize())=m_flow.KnownScalarKLoad(dk);
        result.head(nf)=m_flow.ApplyLinearPreconditioner(rhs.head(nf)+m_flow.Restrict(cross));
        auto trace=m_e.EssentialLift([&](double x,double,const Vector&normal){
            if(normal[2]>=-.5)return 0.;const Element*owner=nullptr;
            for(const auto&e:m_flow.Elements())if(std::abs(e.rl-m_geometry.H)<1e-12&&x>=e.xl-1e-12&&x<=e.xr+1e-12){owner=&e;break;}
            Require(owner!=nullptr,"steady Gamma derivative has no outer owner");
            auto rule=m_rules.Points(*owner,Vec::Constant(1,x),Vec::Constant(1,m_geometry.d0+m_geometry.H));
            return m_reference->AmplitudeDerivative(InterfaceK(*owner,x,pcState))*m_k.Evaluate(rule,dk)[0]/m_geometry.H;
        },m_rules);
        result.tail(ne)=m_e.ApplyLinearPreconditioner(rhs.tail(ne)-m_e.Restrict(epsilonPCFull*trace));return result;
    };
    if(m_options.checkCurrentActions){
        // Uncommitted audit only: explicitly exercise a weakly divergent
        // free-velocity state and the otherwise often-zero gauge multiplier.
        Vec probe=pack(history.current);probe[0]+=1e-7;
        if(m_flow.Gauge())probe[nf-1]=.003;
        (void)equations(probe);
    }
    SteadyKrylov solver(std::move(callbacks),sx,sr,options);auto result=solver.Solve(pack(history.current));
    if(result.converged){
        auto state=unpack(result.x);auto d=Evaluate(state,0,true,timing);auto load=forcing;auto adv=Advection(state,d.kRules);load.flow-=adv.flow;load.k-=adv.k;load.epsilon-=adv.epsilon;
        final.residual=Residual(state,d,load,0,0,timing);
        if(final.residual.Maximum()>options.acceptanceNorm){result.converged=false;result.status="NOT_CONVERGED";result.reason="fresh physical steady residual failed after eliminating the gauge multiplier";}
        else{
            CoupledHistory fresh;fresh.dt=m_options.dt;fresh.current=state;fresh.previous=state;fresh.advection=adv;fresh.previousAdvection=adv;CheckHistory(fresh);history=std::move(fresh);
            final.minimumK=d.minK;final.minimumEpsilon=d.minE;final.minimumCovarianceOverK=d.minCovariance;final.minimumPrincipalOverMu=d.minPrincipal;
        }
    }
    final.nonlinearIterations=int(result.statistics.acceptedNewtonSteps);final.rebuilds=int(3*result.statistics.preconditionerPreparations);
    final.seconds=Clock()-began;final.constitutiveStats=m_constitutiveStats;final.timing.other=final.seconds-final.timing.Sum();
    return result;
}
void Coupling::ResetPreconditioners()
{
    m_flow.ResetPreconditioner();m_k.ResetPreconditioner();m_e.ResetPreconditioner();
    m_preconditionerCoefficients.resize(0);m_kPreconditionerCoefficients.resize(0);m_ePreconditionerCoefficients.resize(0);
    m_flowFactorAlpha=m_kFactorAlpha=m_eFactorAlpha=-1;
    m_timePCReady=false;m_timePCAlpha=-1;m_timePCUses=0;
    for(auto&matrix:m_timePCOperators)matrix.resize(0,0);
    m_timePCEpsilonFull.resize(0,0);
    m_timePCData.resize(0);
}

void Coupling::CheckHistory(const CoupledHistory&h)const
{
    Require(h.dt==m_options.dt&&std::isfinite(h.time)&&std::abs(h.time-h.step*h.dt)<=2e-13*std::max(1.,h.time),"coupled scientific time/history differs");
    for(const auto*s:{&h.current,&h.previous,&h.advection,&h.previousAdvection})Require(s->flow.size()==2*m_flow.VelocitySize()+m_flow.PressureSize()&&s->k.size()==m_k.Size()&&s->epsilon.size()==m_e.Size()&&s->flow.allFinite()&&s->k.allFinite()&&s->epsilon.allFinite(),"coupled history inventory/nonfinite data");
}
void Coupling::Store(const CoupledState&s)const
{
    for(int component=0;component<5;++component){Vec global=component<2?m_flow.NativeVelocity(s.flow,component):component==2?s.k:component==3?s.epsilon:Vec(s.flow.tail(m_flow.PressureSize()));
        Array<OneD,NekDouble>view(global.size(),global.data());auto field=m_fields[component];field->GetLocalToGlobalMap()->GlobalToLocal(view,field->UpdateCoeffs());field->BwdTrans(field->GetCoeffs(),field->UpdatePhys());field->SetPhysState(true);}
}
void Coupling::WriteHistory(const std::filesystem::path&path,const CoupledHistory&h)const
{
    CheckHistory(h);Require(!std::filesystem::exists(path),"coupled checkpoint is immutable");std::ostringstream stream(std::ios::binary);stream<<"WJ_H27_NATIVE_COUPLED_V2\n"<<m_signature<<'\n';
    stream.write(reinterpret_cast<const char*>(&h.step),sizeof(h.step));stream.write(reinterpret_cast<const char*>(&h.time),sizeof(double));stream.write(reinterpret_cast<const char*>(&h.dt),sizeof(double));
    for(const auto*s:{&h.current,&h.previous,&h.advection,&h.previousAdvection})for(const auto*v:{&s->flow,&s->k,&s->epsilon}){std::uint64_t n=v->size();stream.write(reinterpret_cast<const char*>(&n),sizeof(n));stream.write(reinterpret_cast<const char*>(v->data()),n*sizeof(double));}
    auto payload=stream.str();payload+=D::Hash(payload);auto tmp=path.string()+".tmp";int fd=::open(tmp.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);D::Check(fd>=0,"coupled checkpoint create");D::WriteAll(fd,payload.data(),payload.size());D::Check(::fsync(fd)==0,"coupled checkpoint sync");::close(fd);std::filesystem::rename(tmp,path);D::SyncDirectory(path);
}
CoupledHistory Coupling::ReadHistory(const std::filesystem::path&path)const
{
    std::ifstream file(path,std::ios::binary);std::string bytes((std::istreambuf_iterator<char>(file)),{});Require(bytes.size()>64&&bytes.size()<128*1024*1024,"invalid coupled checkpoint size");auto payload=bytes.substr(0,bytes.size()-64);Require(D::Hash(payload)==bytes.substr(bytes.size()-64),"coupled checkpoint checksum mismatch");std::istringstream stream(payload,std::ios::binary);std::string magic,identity;std::getline(stream,magic);std::getline(stream,identity);Require(magic=="WJ_H27_NATIVE_COUPLED_V2"&&identity==m_signature,"coupled checkpoint scientific identity differs");
    CoupledHistory h;stream.read(reinterpret_cast<char*>(&h.step),sizeof(h.step));stream.read(reinterpret_cast<char*>(&h.time),sizeof(double));stream.read(reinterpret_cast<char*>(&h.dt),sizeof(double));
    for(auto*s:{&h.current,&h.previous,&h.advection,&h.previousAdvection}){int which=0;for(auto*v:{&s->flow,&s->k,&s->epsilon}){std::uint64_t n=0;stream.read(reinterpret_cast<char*>(&n),sizeof(n));const int expected=which==0?2*m_flow.VelocitySize()+m_flow.PressureSize():which==1?m_k.Size():m_e.Size();++which;Require(n==std::uint64_t(expected),"coupled checkpoint vector inventory differs");v->resize(n);stream.read(reinterpret_cast<char*>(v->data()),n*sizeof(double));}}
    Require(stream.good()&&stream.peek()==std::char_traits<char>::eof(),"coupled checkpoint has truncated/trailing bytes");CheckHistory(h);
    // Reconstruct scientific derived data and compare the explicit histories.
    // Normal factory cold restarts use a new Coupling instance. External
    // same-instance reloads call ResetPreconditioners for cold inverse setup.
    CoupledTiming timing;auto data=Evaluate(h.current,h.time,true,timing);Require(data.traceError<1e-10,"restored dependent trace differs");auto adv=Advection(h.current,data.kRules);
    for(const auto pair:{std::make_pair(&adv.flow,&h.advection.flow),std::make_pair(&adv.k,&h.advection.k),std::make_pair(&adv.epsilon,&h.advection.epsilon)})Require((*pair.first-*pair.second).norm()<=1e-12*std::max(pair.first->norm(),1e-30),"restored explicit history differs from current state");
    auto previous=Evaluate(h.previous,h.step>0?(h.step-1)*h.dt:0,true,timing);auto oldAdv=Advection(h.previous,previous.kRules);
    for(const auto pair:{std::make_pair(&oldAdv.flow,&h.previousAdvection.flow),std::make_pair(&oldAdv.k,&h.previousAdvection.k),std::make_pair(&oldAdv.epsilon,&h.previousAdvection.epsilon)})Require((*pair.first-*pair.second).norm()<=1e-12*std::max(pair.first->norm(),1e-30),"restored previous explicit history differs");
    return h;
}

std::array<double,5>Coupling::MMSErrors(const CoupledState&s,double time)const
{
    Require(bool(m_mms),"MMS errors requested for application");std::array<double,5>error{};double outerVolume=0;const int n=m_flow.VelocitySize();
    for(const auto&rp:m_rules.Static()){
        const auto&r=*rp;Vec u=r.B*Local(s.flow,r.u),v=r.B*Local(s.flow,r.u,n),p=r.P*Local(s.flow,r.p,2*n),k=m_k.Evaluate(r,s.k),e;
        bool outer=m_e.Owns(r.geometry);if(outer)e=m_e.Evaluate(r,s.epsilon);
        for(int q=0;q<r.w.size();++q){auto exact=m_mms->Fields(r.x[q],r.z[q],time);error[0]+=r.w[q]*(std::pow(u[q]-exact.u.value,2)+std::pow(v[q]-exact.w.value,2));error[1]+=r.w[q]*std::pow(p[q]-exact.Pi.value,2);error[2]+=r.w[q]*std::pow((k[q]-exact.k.value)/m_reference->Parameters().KIn,2);
            if(outer){const double ex=m_mms->OuterEpsilon(r.x[q],r.z[q],time).value;error[3]+=r.w[q]*std::pow((e[q]-ex)/ex,2);outerVolume+=r.w[q];}
            Vec ul=Local(s.flow,r.u),vl=Local(s.flow,r.u,n);Tensor G{};G[0][0]=r.Dx.row(q).dot(ul);G[0][2]=r.Dz.row(q).dot(ul);G[2][0]=r.Dx.row(q).dot(vl);G[2][2]=r.Dz.row(q).dot(vl);
            for(int i:{0,2})for(int j:{0,2})error[4]+=r.w[q]*std::pow(G[i][j]-exact.gradient[i][j],2);
        }
    }
    for(int i=0;i<5;++i)error[i]=std::sqrt(error[i]/(i==3?outerVolume:m_volume));return error;
}
void Coupling::Report(std::ostream&out,const CoupledHistory&h,const CoupledStepInfo&info,bool coefficients)const
{
    out<<std::setprecision(17)<<"{\"schema\":\"WJH27NativeCoupledStateV2\",\"model\":\""<<WallModelV2<<"\",\"signature\":\""<<m_signature<<"\",\"source_identity\":\""<<m_options.sourceIdentity<<"\",\"reference_identity\":\""<<m_reference->Identity()<<"\",\"manufactured\":"<<(m_mms?"true":"false")<<",\"step\":"<<h.step<<",\"time\":"<<h.time<<",\"dt\":"<<h.dt<<",\"state_sha256\":\""<<StateHash(h.current)<<"\",\"previous_sha256\":\""<<StateHash(h.previous)<<"\",\"advection_sha256\":\""<<StateHash(h.advection)<<"\",\"previous_advection_sha256\":\""<<StateHash(h.previousAdvection)<<"\",\"sizes\":["<<m_flow.VelocitySize()<<','<<m_flow.PressureSize()<<','<<m_k.Size()<<','<<m_e.Size()<<"],\"gauge\":"<<(m_flow.Gauge()?"true":"false")<<",\"nonlinear_iterations\":"<<info.nonlinearIterations<<",\"linear_iterations\":["<<info.flowIterations<<','<<info.kIterations<<','<<info.epsilonIterations<<"],\"rebuilds\":"<<info.rebuilds<<",\"retries\":"<<info.retries<<",\"seconds\":"<<info.seconds<<",\"coefficient_drift\":"<<info.coefficientDrift<<",\"residuals\":["<<info.residual.momentum<<','<<info.residual.continuity<<','<<info.residual.gauge<<','<<info.residual.k<<','<<info.residual.epsilon<<','<<info.residual.trace<<"],\"minimum_k\":"<<info.minimumK<<",\"minimum_epsilon\":"<<info.minimumEpsilon<<",\"minimum_covariance_over_k\":"<<info.minimumCovarianceOverK<<",\"minimum_principal_over_mu\":"<<info.minimumPrincipalOverMu<<",\"disjoint_seconds\":{\"constitutive\":"<<info.timing.constitutive<<",\"wall\":"<<info.timing.wall<<",\"assembly\":"<<info.timing.assembly<<",\"flow\":"<<info.timing.flow<<",\"scalar\":"<<info.timing.scalar<<",\"rebuild\":"<<info.timing.rebuild<<",\"guards\":"<<info.timing.guards<<",\"forcing\":"<<info.timing.forcing<<",\"residual\":"<<info.timing.residual<<",\"nonlinear\":"<<info.timing.nonlinear<<",\"io\":"<<info.timing.io<<",\"other\":"<<info.timing.other<<'}';
    out<<",\"quadratureX\":"<<m_options.quadratureX<<",\"quadratureR\":"<<m_options.quadratureR<<",\"linear_tolerance\":"<<m_options.linearTolerance<<",\"nonlinear_tolerance\":"<<m_options.nonlinearTolerance<<",\"maximum_iterations\":"<<m_options.maximumIterations<<",\"mms_body_amplitude\":"<<m_options.mmsBodyAmplitude<<",\"retry_reason\":"<<RANSWallSolveJson::Quote(info.retryReason);
    out<<",\"mms_variant\":"<<RANSWallSolveJson::Quote(m_options.mmsVariant);
    out<<",\"implicit_solver\":"<<RANSWallSolveJson::Quote(m_options.advancement)
        <<",\"coupled_newton\":{\"krylov_iterations\":"<<info.coupledKrylovIterations
        <<",\"jacobian_actions\":"<<info.coupledJacobianActions
        <<",\"residual_evaluations\":"<<info.coupledResidualEvaluations
        <<",\"preconditioner_applications\":"<<info.coupledPreconditionerApplications
        <<",\"rejected_trials\":"<<info.coupledRejectedTrials
        <<",\"counts_complete\":"<<(info.coupledStatisticsComplete?"true":"false")
        <<",\"step_seconds_cap\":"<<m_options.stepNewton.maximumSeconds
        <<",\"uncommitted_guard\":\"positivity/root/covariance\",\"physical_acceptance_guard\":\"fresh full planar principal and complete current residual\"";
    WriteStepPreconditionerFields(out,m_options,info);out<<'}';
    out<<",\"nonlinear_trace_scope\":"<<RANSWallSolveJson::Quote(m_options.advancement=="CoupledJFNKV1"?
        "Uncommitted Newton bases and Armijo candidates; FD excluded; records are not accepted updates; times are inclusive observations":
        "Completed Picard iterations and guarded acceleration proposals; times are inclusive observations");
    out<<",\"steady_velocity_preconditioner\":"<<RANSWallSolveJson::Quote(m_options.steadyVelocityPC);
    out<<",\"stress_split\":"<<RANSWallSolveJson::Quote(m_options.stressSplit)<<",\"all_row_action_crosschecks\":"<<(m_options.checkCurrentActions?"true":"false");
    out<<",\"current_action_checks\":{\"all_row_calls\":"<<m_actionChecks<<",\"maximum_scaled_differences_flow_k_epsilon\":["<<m_actionDifference[0]<<','<<m_actionDifference[1]<<','<<m_actionDifference[2]<<"],\"augmented_calls\":"<<m_augmentedChecks<<",\"maximum_augmented_scaled_difference\":"<<m_augmentedDifference<<",\"maximum_weak_divergence_load_norm\":"<<m_checkedWeakDivergence<<",\"maximum_gauge_multiplier\":"<<m_checkedGaugeMultiplier<<'}';
    out<<",\"basis\":\"Modified_A\",\"anderson_depth\":"<<m_options.andersonDepth<<",\"acceleration\":{\"proposals\":"<<info.accelerationProposals<<",\"accepted\":"<<info.accelerationAccepted<<",\"inadmissible_trials\":"<<info.accelerationInadmissible<<'}';
    out<<",\"initial_field_transfer\":{\"packet_sha256\":"<<RANSWallSolveJson::Quote(m_initialPacketSha)<<",\"source_state_sha256\":"<<RANSWallSolveJson::Quote(m_initialSourceStateSha)<<",\"projected\":"<<(m_initialProjected?"true":"false")<<",\"source_was_candidate\":"<<(m_initialCandidate?"true":"false")<<'}';
    if(m_options.recordNonlinearTrace){out<<",\"nonlinear_trace\":[";for(int i=0;i<int(info.nonlinearTrace.size());++i){const auto&t=info.nonlinearTrace[i];if(i)out<<',';out<<"{\"residuals\":["<<t.residual.momentum<<','<<t.residual.continuity<<','<<t.residual.gauge<<','<<t.residual.k<<','<<t.residual.epsilon<<','<<t.residual.trace<<"],\"accelerated\":"<<(t.accelerated?"true":"false")<<",\"seconds\":"<<t.seconds<<'}';}out<<']';}

    out<<",\"constitutive_backend\":\""<<(m_options.backend==ConstitutiveBackend::Historical?"Historical":"FastFull3DV1")<<"\",\"constitutive_calls\":{\"values\":"<<info.constitutiveStats.values<<",\"directions\":"<<info.constitutiveStats.directions<<",\"fallbacks\":"<<info.constitutiveStats.fallbacks<<",\"root_iterations\":"<<info.constitutiveStats.rootIterations<<",\"maximum_normalized_root_residual\":"<<info.constitutiveStats.maximumRootResidual<<'}';
    if(m_mms){auto errors=MMSErrors(h.current,h.time);out<<",\"mms_rms_errors_U_Pi_k_epsilon_gradU\":[";for(int i=0;i<5;++i){if(i)out<<',';out<<errors[i];}out<<']';}
    if(coefficients){
        auto state=[&](const CoupledState&s){out<<"{\"flow\":";VectorJSON(out,s.flow);out<<",\"k\":";VectorJSON(out,s.k);out<<",\"epsilon\":";VectorJSON(out,s.epsilon);out<<'}';};
        auto integers=[&](const std::vector<int>&v){out<<'[';for(int i=0;i<int(v.size());++i){if(i)out<<',';out<<v[i];}out<<']';};
        auto numbers=[&](const std::vector<double>&v){out<<'[';for(int i=0;i<int(v.size());++i){if(i)out<<',';out<<v[i];}out<<']';};
        out<<",\"coefficients\":";state(h.current);
        auto essential=DirichletState(h.current,h.time);
        out<<",\"essential_flow_lift\":";VectorJSON(out,m_flow.Lift(essential.flow));
        out<<",\"essential_k_lift\":";VectorJSON(out,m_k.Lift(essential.k));
        out<<",\"essential_epsilon_lift\":";VectorJSON(out,m_e.Lift(essential.epsilon));
        out<<",\"mms_bulk_body_rhs\":";VectorJSON(out,m_bulkBodyAction);out<<",\"mms_sampled_body_rhs\":";VectorJSON(out,m_sampledBodyAction);
        out<<",\"current_advection\":";state(h.advection);out<<",\"previous_advection\":";state(h.previousAdvection);
        out<<",\"geometry\":{\"d0\":"<<m_geometry.d0<<",\"H\":"<<m_geometry.H<<",\"h\":"<<m_geometry.h<<",\"left\":"<<m_left<<",\"right\":"<<m_right<<",\"top\":"<<m_top<<"},\"free_velocity_ids\":[";integers(m_flow.FreeVelocityIds()[0]);out<<',';integers(m_flow.FreeVelocityIds()[1]);out<<"],\"k_free_ids\":";integers(m_k.FreeIds());out<<",\"epsilon_free_ids\":";integers(m_e.FreeIds());
        if(m_hasInputHistory){const auto&a=m_lastInputHistory;out<<",\"input_history\":{\"step\":"<<a.step<<",\"time\":"<<a.time<<",\"dt\":"<<a.dt<<",\"current\":";state(a.current);out<<",\"previous\":";state(a.previous);out<<",\"advection\":";state(a.advection);out<<",\"previous_advection\":";state(a.previousAdvection);out<<'}';}
        out<<",\"elements\":[";int count=0;
        for(const auto&e:m_flow.Elements()){
            if(count++)out<<',';out<<"{\"geometry\":"<<e.geometry<<",\"bounds\":["<<e.xl<<','<<e.xr<<','<<e.rl<<','<<e.rr<<"],\"u\":";integers(e.u);out<<",\"us\":";numbers(e.us);out<<",\"p\":";integers(e.p);out<<",\"ps\":";numbers(e.ps);
            for(int field:{2,3}){std::vector<int>ids;std::vector<double>signs;auto f=m_fields[field];auto map=f->GetLocalToGlobalMap();
                for(int j=0;j<f->GetExpSize();++j)if(f->GetExp(j)->GetGeom()->GetGlobalID()==e.geometry){for(int k=0;k<16;++k){int index=f->GetCoeff_Offset(j)+k;ids.push_back(map->GetLocalToGlobalMap(index));signs.push_back(map->GetLocalToGlobalSign(index));}break;}
                out<<",\""<<(field==2?"k":"epsilon")<<"\":";integers(ids);out<<",\""<<(field==2?"ks":"epsilon_signs")<<"\":";numbers(signs);
            }
            out<<'}';
        }
        out<<']';
    }
    out<<'}';
}
}
