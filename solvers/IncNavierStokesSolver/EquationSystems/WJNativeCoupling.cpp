// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupling.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeCutoffPanels.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeWallV2.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeFailures.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeAdmissibility.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeInitialFields.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeFlowTangent.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativePerturbation.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <Eigen/Eigenvalues>
#include <Eigen/QR>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
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
template<class Write> void PublishDiagnostic(const std::filesystem::path &path,Write write)
{
    const auto temporary=path.string()+".tmp."+std::to_string(::getpid());
    const int fd=::open(temporary.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);
    D::Check(fd>=0,"cannot reserve diagnostic packet");
    try
    {
        write(fd);
        D::Check(::fsync(fd)==0,"cannot sync diagnostic packet");
    }
    catch(...)
    {
        ::close(fd);
        ::unlink(temporary.c_str());
        throw;
    }
    const int closed=::close(fd);
    try
    {
        D::Check(closed==0,"cannot close diagnostic packet");
        std::filesystem::create_hard_link(temporary,path);
    }
    catch(...)
    {
        ::unlink(temporary.c_str());
        throw;
    }
    D::Check(::unlink(temporary.c_str())==0,"cannot unlink diagnostic temporary");
    D::SyncDirectory(std::filesystem::absolute(path));
}
CoupledState Combination(const CoupledState&a,double wa,const CoupledState&b,double wb)
{return {wa*a.flow+wb*b.flow,wa*a.k+wb*b.k,wa*a.epsilon+wb*b.epsilon};}
std::array<double,2> HillWindow(double a,double left,double right,double band)
{
    const double width=a<0 ? -band-left : right-band;
    const double xi=std::clamp(a<0 ? (a-left)/width : (right-a)/width,0.,1.);
    return {xi*xi*xi*(10+xi*(-15+6*xi)),
            (a<0 ? 1. : -1.)*30*xi*xi*(1-xi)*(1-xi)/width};
}
std::array<double,2> SmoothTransition(double value,double lower,double upper)
{
    const double t=std::clamp((value-lower)/(upper-lower),0.,1.);
    return {t*t*t*(10+t*(-15+6*t)),30*t*t*(1-t)*(1-t)/(upper-lower)};
}
}

Coupling::Coupling(std::vector<MultiRegions::ContFieldSharedPtr>fields,Geometry geometry,
    std::shared_ptr<const Reference>reference,CoupledOptions options)
    :m_fields(std::move(fields)),m_geometry(geometry),m_reference(std::move(reference)),m_options(std::move(options)),
     m_flow({m_fields.at(0),m_fields.at(1),m_fields.at(4)},geometry,
            m_options.flowPreconditioner == "PressureCorrection" ? m_options.pressureCorrectors : 0),
     m_rules(m_flow.Core(),m_options.quadratureX,m_options.quadratureR),
     m_k(m_fields.at(2),m_flow.Core(),[this](double,double){return m_reference->Parameters().KIn;},{0}),
     m_e(m_fields.at(3),m_flow.Core(),[this](double,double z){return m_reference->Amplitude(m_reference->Parameters().KIn)/(z-m_geometry.d0);},{0,2})
{
    Require(m_fields.size()==5&&bool(m_reference),"coupling requires native u,v,k,outer-epsilon,p and a qualified reference");
    Require(m_options.dt>0&&m_options.linearTolerance>0&&m_options.linearTolerance<=1e-10&&
        m_options.nonlinearTolerance>0&&m_options.nonlinearTolerance<=1e-8&&m_options.maximumIterations>=2,
        "invalid coupled timestep or convergence policy");
    Require(m_options.sourceIdentity.size()==64,"coupled scientific source identity is missing");
    m_zeroTransverseOutlet=m_options.outletPolicy=="ZeroTransverseV1";
    Require(m_zeroTransverseOutlet||m_options.outletPolicy=="FullTractionV1","unknown WJ outlet policy");
    Require(!m_zeroTransverseOutlet||(!m_options.manufactured&&
        m_reference->Policy()==ReferencePolicy::ExactInnerFiniteNuV1&&m_reference->Transport()==InnerKTransport::ScalarV2),
        "zero-transverse outlet requires the exact-inner ScalarCdV2 physical model");
    Require(m_options.diagnosticOutletX==15.||
        (m_zeroTransverseOutlet&&std::isfinite(m_options.diagnosticOutletX)&&m_options.diagnosticOutletX>15.),
        "outlet location diagnostic requires a downstream zero-transverse outlet in the supported range");
    m_outletIdentity=D::Hash(m_reference->BoundaryIdentity()+"|"+m_options.outletPolicy);
    Require(m_options.andersonDepth>=0&&m_options.andersonDepth<=8,"Anderson depth must lie in [0,8]");
    Require(m_options.stressSplit=="PositiveShearV1"||m_options.stressSplit=="IsotropicMuV1","unknown numerical stress split");
    Require(m_options.steadyVelocityPC=="FullGradientV1"||m_options.steadyVelocityPC=="PositiveFrozenV1","unknown steady velocity preconditioner");
    Require(m_options.advancement=="PicardV1"||m_options.advancement=="CoupledJFNKV1","unknown implicit coupled solver");
    Require(m_options.newtonAdvection == "Explicit" ||
            (m_options.newtonAdvection == "Implicit" && m_options.advancement == "CoupledJFNKV1"),
            "WJNewtonAdvection requires Explicit or CoupledJFNKV1 with Implicit");
    Require(m_options.newtonAdvection != "Implicit" || m_options.nonlinearTolerance <= 1e-10,
            "Implicit Newton requires a refreshed coupled tolerance of at most 1e-10");
    Require(m_options.picardAdvection == "Explicit" ||
            (m_options.picardAdvection == "Implicit" && m_options.advancement == "PicardV1"),
            "WJPicardAdvection requires Explicit or PicardV1 with Implicit");
    Require(m_options.picardStressLinearization == "FrozenStrain" ||
            (m_options.picardStressLinearization == "VelocityTangent" &&
             m_options.picardAdvection == "Implicit"),
            "WJPicardStressLinearization requires FrozenStrain or Implicit with VelocityTangent");
    const bool pressureCorrection = m_options.flowPreconditioner == "PressureCorrection";
    const bool fixedBudget = m_options.picardStop == "FixedCorrectors";
    Require(pressureCorrection || m_options.flowPreconditioner == "SparseLU",
            "unknown WJ flow preconditioner");
    Require(m_options.pressureCorrectors > 0, "WJPressureCorrectors must be positive");
    Require(fixedBudget || m_options.picardStop == "CoupledResidual", "unknown WJ Picard stopping rule");
    Require(fixedBudget ? (m_options.advancement == "PicardV1" && pressureCorrection &&
                          m_options.picardOuterCorrectors > 0 && m_options.andersonDepth == 0)
                       : m_options.picardOuterCorrectors == 0,
            "fixed budgets require PicardV1, PressureCorrection, positive outer count and zero Anderson depth");
    Require(!pressureCorrection || m_options.stepPreconditioner != "CurrentDataAMDOrderingV1",
            "WJ PressureCorrection uses the COLAMD ordering for both sparse factors");
    Require(m_options.stepPreconditioner=="RebuildEveryNewtonV1"||m_options.stepPreconditioner=="OperatorDriftV1"||
            m_options.stepPreconditioner=="CurrentDataDriftV1"||m_options.stepPreconditioner=="CurrentDataAMDOrderingV1",
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
        "coupled v2 requires the pinned model heights on the supported chart mesh");
    Require(!m_geometry.normalMap || (m_reference->Policy()==ReferencePolicy::ExactInnerFiniteNuV1 &&
            m_reference->Transport()==InnerKTransport::ScalarV2),
        "the normal layer requires the exact finite-viscosity ScalarCdV2 reference");
    Require(m_options.hillInitialSlowdown==0. ||
            (m_geometry.normalMap && !m_options.manufactured &&
             m_options.hillInitialSlowdown>0. && m_options.hillInitialSlowdown<=.9),
            "initial slowdown requires a physical hill and an amplitude in (0,.9]");
    Require(m_options.hillInitialEpsilonBump==0. ||
            (m_geometry.normalMap && !m_options.manufactured &&
             m_options.hillInitialEpsilonBump>0. && m_options.hillInitialEpsilonBump<=6.),
            "initial epsilon bump requires a physical hill and an amplitude in (0,6]");
    for(const auto&field:m_fields)for(int e=0;e<field->GetExpSize();++e)for(int d=0;d<2;++d)
        Require(field->GetExp(e)->GetBasis(d)->GetBasisType()==LibUtilities::eModified_A,"coupled fields require the declared native Modified_A basis");
    m_flow.SetTransportedKSpace(m_fields[2]);BuildTopology();
    m_cutoffRange=m_reference->CutoffRadiusRange();
    for(const auto&e:m_flow.Elements())for(double radius:{e.rl,e.rr})
        if(radius>=m_cutoffRange[0]&&radius<=m_cutoffRange[1]&&!m_cutoffFaceK.count(radius))
            m_cutoffFaceK.emplace(radius,m_reference->KAtCutoff(radius));
    Require(std::isfinite(m_options.mmsBodyAmplitude),"nonfinite MMS body amplitude");
    Require(m_options.manufactured||m_options.mmsBodyAmplitude==0,"the prescribed MMS body is verification-only");
    Require(m_options.mmsVariant=="Historical"||m_options.mmsVariant=="PolynomialScalarsV1"||
            m_options.mmsVariant=="NormalTubeV1","unknown continuous MMS variant");
    Require(m_options.manufactured||m_options.mmsVariant=="Historical","a continuous MMS variant is verification-only");
    Require(m_reference->Transport()==InnerKTransport::ScalarV2||!m_options.manufactured,
            "inner-k diagnostic transport is restricted to the unforced physical steady fixture");
    if(m_options.manufactured)m_mms=std::make_unique<CoupledMMS>(*m_reference,m_geometry,m_top,m_options.mmsBodyAmplitude,
        m_options.mmsVariant=="NormalTubeV1"?MMSVariant::NormalTubeV1:
        m_options.mmsVariant=="Historical"?MMSVariant::Historical:MMSVariant::PolynomialScalarsV1);
    else Require(std::abs(m_geometry.d0-m_reference->Parameters().d0)<1e-14&&
        std::abs(m_top-m_reference->Parameters().topRadius)<1e-12&&std::abs(m_left+10)<1e-12&&
        std::abs(m_right-m_options.diagnosticOutletX)<1e-12,
        "application geometry differs from the retained physical domain");
    const char *integration = m_options.UsesImplicitAdvection() ?
        "BDF2-ImplicitAdvection" : "BDF2-EX2";
    std::ostringstream signature;signature<<std::setprecision(17)<<"WJ-H27-coupled-v2-"<<integration<<'|'<<WallModelV2<<"|CubicCutoffFacePanelsV1|"
        <<m_options.sourceIdentity<<'|'<<m_reference->Identity()<<'|'<<m_options.dt<<'|'<<m_options.manufactured<<'|'
        <<m_options.quadratureX<<'|'<<m_options.quadratureR<<'|'<<m_geometry.d0<<'|'<<m_geometry.H<<'|'<<m_geometry.h<<'|'
        <<m_options.linearTolerance<<'|'<<m_options.nonlinearTolerance<<'|'<<m_options.maximumIterations<<'|'<<m_options.mmsBodyAmplitude<<'|'<<int(m_options.backend)<<'|'<<m_options.mmsVariant<<'|'<<m_options.andersonDepth<<'|'<<m_options.stressSplit<<'|'<<m_options.steadyVelocityPC;
    if (m_options.picardStressLinearization == "VelocityTangent")
    {
        signature << "|PicardVelocityTangent";
    }
    const auto&newton=m_options.stepNewton;
    signature<<'|'<<m_options.advancement<<'|'<<newton.maximumSeconds<<'|'<<newton.acceptanceNorm<<'|'
        <<newton.maximumNewtonIterations<<'|'<<newton.restartDimension<<'|'<<newton.maximumLinearCycles<<'|'<<newton.fdRelativeStep;
    if (m_options.advancement == "CoupledJFNKV1")
    {
        signature << "|NewtonAdvection:" << m_options.newtonAdvection
                  << "|forcing:" << newton.forcingInitial << ',' << newton.forcingMinimum
                  << ',' << newton.forcingMaximum << ',' << newton.forcingGamma
                  << ',' << newton.forcingExponent;
    }
    signature<<'|'<<m_options.stepPreconditioner<<'|'<<m_options.stepPCRelativeDrift<<'|'<<m_options.stepPCMaximumUses;
    if (pressureCorrection)
    {
        signature << '|' << m_options.flowPreconditioner << '|' << m_options.pressureCorrectors
                  << '|' << m_options.picardStop << '|' << m_options.picardOuterCorrectors;
        if (fixedBudget)
        {
            signature << "|MomentumSchurV1";
        }
    }
    if(m_reference->Policy()!=ReferencePolicy::Historical)
        signature<<'|'<<m_reference->PolicyName()<<'|'<<m_reference->TransportName();
    if(m_reference->Transport()!=InnerKTransport::ScalarV2)
        signature<<'|'<<m_reference->TransportName()<<'|'<<m_reference->BoundaryIdentity();
    signature<<"|outlet:"<<m_outletIdentity;
    if (m_geometry.normalMap)
        signature << "|NormalCosineHillV1|" << m_geometry.normalMap->Signature()
                  << "|initial_slowdown:" << m_options.hillInitialSlowdown
                  << "|initial_epsilon_bump:" << m_options.hillInitialEpsilonBump;
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
    if (m_geometry.normalMap)
    {
        m_volume=0.;
        for (const auto &e : elements)
            m_volume += e.w.sum();
    }
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
            const bool essential=m_options.manufactured||side==0||side==2||
                (component==1&&(side==3||(side==1&&m_zeroTransverseOutlet)));
            for(int mode:edgeModes[side])Require((free[component].count(e.u[mode])==0)==essential,
                "native component boundary mask disagrees with the declared physical side law");
        }
        if(side==1&&!m_options.manufactured){
            const bool ground=std::abs(e.rl)<1e-12;
            const bool top=std::abs(e.rr-m_top)<1e-12;
            for(int mode:{1,5,9,13}){
                const bool sharedGround=mode==1&&ground;
                const bool sharedTop=mode==5&&top;
                const bool essentialW=m_zeroTransverseOutlet||sharedGround||sharedTop;
                const int id=e.u[mode];
                Require((free[0].count(id)==0)==sharedGround&&(free[1].count(id)==0)==essentialW,
                        "native outlet vertex or edge component does not match the declared law");
                if(essentialW)m_outletReactionRows.emplace(id,sharedGround||sharedTop);
            }
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
const Element& Coupling::MatchingOwner(double a) const
{
    for (const auto &entry : m_hKTrace)
    {
        const auto &cell = Owner(entry.first);
        if (a >= cell.xl-1e-12 && a <= cell.xr+1e-12)
            return cell;
    }
    throw std::runtime_error("dependent epsilon trace has no outer matching owner");
}
NativeRule Coupling::HRule(const NativeRule&rule)const
{return m_rules.HorizontalTrace(HOwner(Owner(rule.geometry)),rule.x,2);}
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
    const double width=owner.xr-owner.xl, xi=(x-owner.xl)/width;
    const auto u=CubicTrace(c,xi,width);
    if (!m_geometry.normalMap)
        return u;
    for (int j=0;j<4;++j)
        c[j]=s.flow[m_flow.VelocitySize()+owner.u[j]]*owner.us[j];
    const auto w=CubicTrace(c,xi,width);
    const auto metric=m_geometry.At(x,m_geometry.H,owner.xl,owner.xr,false);
    const auto &t=metric.tangent1;
    const auto &n=metric.normal;
    const double curvature=metric.kappa1+metric.kappa2;
    return {t[0]*u[0]+t[2]*w[0],
            (t[0]*u[1]+t[2]*w[1])/metric.determinant -
                curvature/(1+m_geometry.H*curvature)*(n[0]*u[0]+n[2]*w[0])};
}
std::vector<double>Coupling::CutoffXBreaks(const Element&e,const std::array<double,4>&c,double left,double right)const
{
    std::vector<double>levels;
    for(double radius:{e.rl,e.rr}){auto found=m_cutoffFaceK.find(radius);if(found!=m_cutoffFaceK.end())levels.push_back(found->second);}
    auto parameters=CubicCutoffPanels(c,levels);std::vector<double>points;
    const double roundoff=64*std::numeric_limits<double>::epsilon()*std::max({std::abs(e.xl),std::abs(e.xr),e.xr-e.xl});
    for(double t:parameters){const double x=left+(right-left)*t;if(x>e.xl+roundoff&&x<e.xr-roundoff)points.push_back(x);}
    return points;
}
std::vector<double>Coupling::CurrentCutoffXBreaks(const Element&e,const CoupledState&s)const
{
    const auto&owner=HOwner(e);const auto&trace=m_hKTrace.at(owner.geometry);std::array<double,4>c{};
    for(int j=0;j<4;++j)c[j]=s.k[trace.ids[j]]*trace.signs[j];
    return CutoffXBreaks(e,c,owner.xl,owner.xr);
}
std::array<double,2> Coupling::InterfaceBounds(const CoupledState &state,bool requireQualified) const
{
    const auto &p=m_reference->Parameters();
    std::array<double,2> result{std::numeric_limits<double>::infinity(),
                                -std::numeric_limits<double>::infinity()};
    for (const auto &[id,trace] : m_hKTrace)
    {
        std::array<double,4> coefficients;
        for (int j=0; j<4; ++j)
            coefficients[j]=state.k[trace.ids[j]]*trace.signs[j];
        const auto bounds=CubicBounds(coefficients);
        if (requireQualified && !(bounds[0]>=p.KMinimum && bounds[1]<=p.KMaximum))
        {
            const auto &owner=Owner(id);
            std::ostringstream error;
            error << std::setprecision(17) << "interface K polynomial left the qualified reference interval"
                  << " in geometry " << id << " (chart a " << owner.xl << " to " << owner.xr
                  << "), K/Kin=[" << bounds[0]/p.KIn << ',' << bounds[1]/p.KIn << ']';
            throw NativeInadmissibleState(error.str());
        }
        result[0]=std::min(result[0],bounds[0]);
        result[1]=std::max(result[1],bounds[1]);
    }
    return result;
}
std::vector<double>Coupling::MMSCutoffXBreaks(const Element&e,double time)const
{
    Require(bool(m_mms),"exact cutoff panels require the declared continuous MMS");
    // Both declared MMS variants have the same exact quadratic interface K.
    // Its own surface is independent of the current numerical trace.
    const double width=e.xr-e.xl;const auto middle=m_mms->InterfaceK((e.xl+e.xr)/2,time);
    const std::array<double,4>c{middle.value-width*middle.dx/2+width*width*middle.dxx/8,
        middle.value+width*middle.dx/2+width*width*middle.dxx/8,-width*width*middle.dxx/2,0.};
    for(double t:{.211,.533,.877}){
        const double exact=m_mms->InterfaceK(e.xl+width*t,time).value;
        Require(std::abs(CubicTrace(c,t,width)[0]-exact)<=64*std::numeric_limits<double>::epsilon()*m_reference->Parameters().KIn,
            "continuous MMS interface K is not the declared quadratic trace");
    }
    return CutoffXBreaks(e,c,e.xl,e.xr);
}
void Coupling::CheckTraceEvaluation()const
{
    CoupledState synthetic;synthetic.k=m_k.Project([&](double,double){return m_reference->Parameters().KIn;});
    synthetic.flow=Vec(2*m_flow.VelocitySize()+m_flow.PressureSize());
    for(int j=0;j<synthetic.k.size();++j)synthetic.k[j]+=1e-4*m_reference->Parameters().KIn*std::sin(double(j+1));
    for(int j=0;j<synthetic.flow.size();++j)synthetic.flow[j]=std::cos(.7*(j+1));
    for(const auto&entry:m_hKTrace){const auto&e=Owner(entry.first);for(double t:{0.,.135,.555,.8,1.}){
        double x=e.xl+t*(e.xr-e.xl);auto rule=m_rules.HorizontalTrace(e,Vec::Constant(1,x),2);
        const double K=m_k.Evaluate(rule,synthetic.k)[0],fast=InterfaceK(e,x,synthetic);
        Require(std::abs(K-fast)<=2e-13*m_reference->Parameters().KIn,"direct native H scalar trace disagrees with owned basis evaluation");
        auto velocity=InterfaceVelocity(e,x,synthetic);
        const auto normal=rule.WallNormal(0);
        const std::array<double,2> tangent{normal[2],-normal[0]};
        const Mat ds=tangent[0]*rule.Dx+tangent[1]*rule.Dz;
        double traceValue=0., derivative=0.;
        for (int i=0;i<2;++i)
        {
            const Vec local=Local(synthetic.flow,rule.u,i*m_flow.VelocitySize());
            traceValue+=tangent[i]*rule.B.row(0).dot(local);
            derivative+=tangent[i]*ds.row(0).dot(local)
                -rule.MatchingCurvature(0)*normal[2*i]*rule.B.row(0).dot(local);
        }
        Require(std::abs(traceValue-velocity[0])<=2e-13&&std::abs(derivative-velocity[1])*(e.xr-e.xl)<=2e-13,
            "direct native H velocity trace/gradient disagrees with owned basis evaluation");
    }}
    if (m_geometry.normalMap && m_options.checkCurrentActions)
    {
        synthetic.epsilon=Vec::Zero(m_e.Size());
        synthetic=DirichletState(synthetic,0.);
        Vec dk=Vec::Zero(m_k.Size());
        for (int id : m_k.FreeIds())
            dk[id]=.01*m_reference->Parameters().KIn*std::sin(.3*(id+1));
        auto trace=[&](double displacement)
        {
            CoupledState state=synthetic;
            state.k+=displacement*dk;
            return m_e.Lift(DirichletState(state,0.).epsilon);
        };
        const double step=.001;
        const Vec finiteDifference=(trace(-2*step)-8*trace(-step)+
            8*trace(step)-trace(2*step))/(12*step);
        const auto derivative=DependentTraceDirection(synthetic,dk);
        m_dependentTraceDerivativeError=(finiteDifference-derivative).norm()/derivative.norm();
        Require(m_dependentTraceDerivativeError<=1e-8,"dependent epsilon lift derivative mismatch");
    }
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
    const double left=m_options.manufactured?0.:-10.,right=m_options.manufactured?1.:m_options.diagnosticOutletX;
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
    Vec lift=m_flow.Lift(full);
    if (m_mms && m_geometry.normalMap)
    {
        // The exact closed fixture has zero net mass flux. Preserve that
        // compatibility when projecting its nonpolynomial boundary values.
        // Only right-side interior trace modes take the correction.
        const int n=m_flow.VelocitySize();
        const Vec flux=m_flow.Core().Divergence().transpose()*m_flow.PressureOne();
        Vec direction=Vec::Zero(2*n);
        for (const auto &[index,side] : m_boundaryFaces)
            if (side==1)
                for (int mode : {9,13})
                {
                    const int id=m_flow.Elements()[index].u[mode];
                    direction[id]=flux[id];
                    direction[n+id]=flux[n+id];
                }
        Require(direction.squaredNorm()>0,"closed lift has no right-side flux mode");
        lift.head(2*n)-=direction*(flux.dot(lift.head(2*n))/direction.squaredNorm());
    }
    if(!m_mms)m_physicalFlowLift=lift;
    return lift;
}
CoupledState Coupling::DirichletState(const CoupledState&s,double time)const
{
    CoupledState result=s;result.flow=m_flow.Expand(m_flow.Restrict(s.flow),FlowLift(time));
    Vec kLift;
    if(!m_mms&&m_physicalKLift.size())kLift=m_physicalKLift;
    else{kLift=m_k.EssentialLift([&](double x,double z,const Vector&){return m_mms?m_mms->Fields(NativeX(x),NativeRadius(z)+m_geometry.d0,time).k.value:m_reference->InletK(NativeRadius(z));},m_rules);if(!m_mms)m_physicalKLift=kLift;}
    result.k=m_k.Expand(m_k.Restrict(s.k),kLift);
    auto eLift=m_e.EssentialLift([&](double x,double z,const Vector&){
        if(NativeRadius(z)==m_geometry.H)
            return m_reference->Amplitude(InterfaceK(MatchingOwner(x),x,result))/m_geometry.H;
        return m_mms?m_mms->OuterEpsilon(NativeX(x),NativeRadius(z)+m_geometry.d0,time).value:m_reference->InletE(NativeRadius(z));
    },m_rules);
    result.epsilon=m_e.Expand(m_e.Restrict(s.epsilon),eLift);return result;
}

Vec Coupling::DependentTraceDirection(const CoupledState &state, const Vec &dk) const
{
    return m_e.EssentialLift([&](double a,double z,const Vector&){
        if (NativeRadius(z)!=m_geometry.H)
            return 0.;
        const auto &owner=MatchingOwner(a);
        const auto rule=m_rules.HorizontalTrace(owner,Vec::Constant(1,a),2);
        return m_reference->AmplitudeDerivative(InterfaceK(owner,a,state))*
            m_k.Evaluate(rule,dk)[0]/m_geometry.H;
    },m_rules);
}
Vector Coupling::InitialVelocity(double a, double r) const
{
    if (!m_geometry.normalMap)
        return {m_reference->InletU(r),0.,0.};

    const auto &map=*m_geometry.normalMap;
    const bool outer=r>m_geometry.H;
    const auto chart=m_geometry.At(a,r,a,a,outer);
    double level=r;
    Vector gradient{chart.inverse[2][0],0.,chart.inverse[2][2]};
    // The streamfunction is the integral of the retained approach U in
    // level. Match the physical distance and gradient at H using the
    // linear normal projection of the ruled displacement there. Recover
    // horizontal approach streamlines above the transition in physical space.
    if (outer && map.amplitude!=0.)
    {
        const auto matching=m_geometry.At(a,m_geometry.H,a,a,false);
        const auto &normal=matching.normal;
        const auto &tangent=matching.tangent1;
        const double curvature=matching.kappa1+matching.kappa2;
        const double radialX=chart.jacobian[0][2],radialZ=chart.jacobian[2][2];
        const double slope=normal[0]*radialX+normal[2]*radialZ;
        const double slopeA=curvature*matching.baseArea*(tangent[0]*radialX+tangent[2]*radialZ)
            +normal[0]/(map.tunnelTop-m_geometry.d0-m_geometry.H);
        const auto [blend,blendR]=SmoothTransition(r,.35,4.);
        const double distance=r-m_geometry.H;
        const double matched=m_geometry.H+distance*slope;
        const double vertical=chart.point[2]-m_geometry.d0;
        level=(1-blend)*matched+blend*vertical;
        const double radialDerivative=(1-blend)*slope+blendR*(vertical-matched);
        for (int i : {0,2})
            gradient[i]=(i==2?blend:0.)+radialDerivative*gradient[i]+
                (1-blend)*distance*slopeA*chart.inverse[0][i];
    }
    if (m_options.hillInitialSlowdown!=0.)
    {
        // Change the initial streamfunction only. Blend back to the full
        // approach profile at the inlet, outlet and roof.
        const auto [window,windowA]=HillWindow(a,m_left,m_right,map.lengthX+.5);
        const auto [turned,turnedR]=SmoothTransition(level,2.5,9.);
        const double fade=1-turned, fadeR=-turnedR;
        const double alpha=m_options.hillInitialSlowdown;
        for(int i : {0,2})
            gradient[i]=(1-alpha*window*(fade+level*fadeR))*gradient[i]-
                alpha*windowA*level*fade*chart.inverse[0][i];
        level*=1-alpha*window*fade;
    }
    const double velocity=m_reference->InletU(level);
    return {velocity*gradient[2],0.,-velocity*gradient[0]};
}

double Coupling::SeedGradientFit::Error(const Vec &velocity) const
{
    const double squared=velocity.dot(gram*velocity)-2*velocity.dot(load)+targetSquared;
    return std::sqrt(std::max(0.,squared)/weight);
}

Coupling::SeedGradientFit Coupling::GeometrySeedGradientFit(const Vec &oldFlow,
                                                            const Geometry &source) const
{
    const int count=m_flow.VelocitySize();
    const double H=m_geometry.H, foot=m_geometry.normalMap->lengthX;
    const double footLength=m_options.geometrySeedFootLength;
    const double collarLength=m_options.geometrySeedCollarLength;
    const double collarReach=m_options.geometrySeedCollarReach;
    SeedGradientFit fit;
    fit.load=Vec::Zero(2*count);
    // Evaluate the reflected quintic directly to retain nonnegative
    // weights even within roundoff of the upper end of the transition.
    const auto fade=[](double value,double lower,double upper)
    {
        return SmoothTransition(upper-value,0.,upper-lower)[0];
    };
    std::vector<Vec> weights;
    weights.reserve(m_flow.Elements().size());
    for (const auto &cell : m_flow.Elements())
    {
        const bool outer=cell.rl>=H-1e-12;
        const bool hasFoot=outer && footLength>0.;
        const bool hasCollar=collarLength>0. && cell.rl<collarReach;
        Vec &weightsHere=weights.emplace_back(Vec::Zero(cell.w.size()));
        if (!hasFoot && !hasCollar)
            continue;
        const auto u=Local(oldFlow,cell.u), w=Local(oldFlow,cell.u,count);
        const Vec ux=cell.Dx*u, uz=cell.Dz*u, wx=cell.Dx*w, wz=cell.Dz*w;
        Mat gradX=Mat::Zero(cell.w.size(),2), gradZ=gradX;
        for (int q=0; q<cell.w.size(); ++q)
        {
            const double a=cell.x[q], r=cell.z[q]-m_geometry.d0;
            const double footWeight=hasFoot ? footLength*footLength*fade(std::abs(a-foot),.25,.80)
                *SmoothTransition(r,4*H,5.5*H)[0]*fade(r,9*H,11*H) : 0.;
            const double collarWeight=hasCollar ? collarLength*collarLength*fade(std::abs(a),.2*foot,.6*foot)
                *SmoothTransition(r,.15*H,.4*H)[0]*fade(r,.8*H,collarReach) : 0.;
            weightsHere[q]=footWeight+collarWeight;
            if (weightsHere[q]==0.)
                continue;
            const auto old=source.At(a,r,cell.xl,cell.xr,outer);
            const auto next=m_geometry.At(a,r,cell.xl,cell.xr,outer);
            Eigen::Matrix2d gradient, J, inverse;
            gradient << ux[q],uz[q],wx[q],wz[q];
            J << next.jacobian[0][0],next.jacobian[0][2],next.jacobian[2][0],next.jacobian[2][2];
            inverse << old.inverse[0][0],old.inverse[0][2],old.inverse[2][0],old.inverse[2][2];
            // Recover the old Cartesian gradient of the source polynomial
            // at the same chart label, from its target-metric basis action.
            gradient=(gradient*J*inverse).eval();
            Eigen::Matrix2d weightedGradient=footWeight*gradient;
            double targetEnergy=footWeight*gradient.squaredNorm();
            if (collarWeight>0.)
            {
                // Compare crown gradients in corresponding wall frames.
                const auto &s=old.tangent1, &t=next.tangent1;
                const double c=s[0]*t[0]+s[2]*t[2], p=s[0]*t[2]-s[2]*t[0];
                Eigen::Matrix2d rotation;
                rotation << c,-p,p,c;
                const Eigen::Matrix2d rotated=rotation*gradient*rotation.transpose();
                weightedGradient+=collarWeight*rotated;
                targetEnergy+=collarWeight*rotated.squaredNorm();
            }
            const double weight=cell.w[q];
            gradX.row(q)=weight*weightedGradient.col(0).transpose();
            gradZ.row(q)=weight*weightedGradient.col(1).transpose();
            fit.targetSquared+=weight*targetEnergy;
            fit.weight+=weight*weightsHere[q];
        }
        for (int component=0; component<2; ++component)
        {
            const Vec load=cell.Dx.transpose()*gradX.col(component)+cell.Dz.transpose()*gradZ.col(component);
            Add(fit.load,cell.u,load,component*count);
        }
    }
    Require(fit.weight>0.,"geometry seed gradient fit has no quadrature support");
    fit.gram=m_flow.Core().StrainSamples(weights,false);
    return fit;
}

Vec Coupling::ProjectContinuity(const Vec &flow, const SeedGradientFit *fit) const
{
    // Constrain continuity while minimizing an H1 velocity correction.
    // A pure mass projection can give steep shear corrections in the very
    // thin ground cells. The regularizer acts only in this initial lift.
    // Discard the projection multiplier. Advancement uses the full G_eff.
    const int nv=2*m_flow.VelocitySize();
    Sparse regularizer(nv,nv);
    if (!m_options.manufactured)
    {
        const double length=.4*m_geometry.H;
        regularizer=m_flow.Core().Strain([length](double,double){return length*length;},false);
    }
    Vec target=m_flow.Core().MassAction(flow);
    target.head(nv)+=regularizer*flow.head(nv);
    if (fit)
    {
        regularizer+=fit->gram;
        target.head(nv)+=fit->load;
    }
    const auto system=m_flow.System(1.,regularizer,0.);
    Eigen::SparseLU<Sparse> factor;
    factor.compute(system);
    Require(factor.info()==Eigen::Success,"initial velocity mass projection factorization failed");
    const Vec lift=m_flow.Lift(flow);
    const Vec rhs=m_flow.LiftedRhs(1.,regularizer,0.,target,lift);
    const Vec projected=factor.solve(rhs);
    Require(factor.info()==Eigen::Success&&projected.allFinite(),
            "initial velocity mass projection failed");
    Vec result=m_flow.Expand(projected,lift);
    result.tail(m_flow.PressureSize())=flow.tail(m_flow.PressureSize());
    return result;
}

double Coupling::ContinuityNorm(const Vec &flow) const
{
    Vec residual=Vec::Zero(flow.size());
    residual.tail(m_flow.PressureSize())=
        m_flow.Core().Divergence()*flow.head(2*m_flow.VelocitySize());
    return m_flow.RawResidualNorm(m_flow.Restrict(residual)).continuity;
}

Vec Coupling::RemapVelocity(const Vec &flow, const Geometry &source)
{
    const int count=m_flow.VelocitySize();
    Vec result=flow;
    if (source.normalMap->amplitude==m_geometry.normalMap->amplitude)
        return result;

    // Preserve the contravariant volume flux on each chart patch. Its
    // tangential velocity can jump where the mapping derivatives change,
    // so project the difference into the continuous Cartesian space with
    // the physical mass rule. The boundary and continuity lifts follow.
    Vec correction=Vec::Zero(2*count);
    m_initialFluxRoundtripError=0.;
    m_initialTangentialTraceCorrection=0.;
    using Sample=SolverUtils::RANSWallNormalMap::Sample;
    auto piola=[&](const Sample &old,const Sample &next,double u,double w)
    {
        const double qa=old.determinant*(old.inverse[0][0]*u+old.inverse[0][2]*w);
        const double qr=old.determinant*(old.inverse[2][0]*u+old.inverse[2][2]*w);
        const double x=(next.jacobian[0][0]*qa+next.jacobian[0][2]*qr)/next.determinant;
        const double z=(next.jacobian[2][0]*qa+next.jacobian[2][2]*qr)/next.determinant;
        const double recoveredA=next.determinant*(next.inverse[0][0]*x+next.inverse[0][2]*z);
        const double recoveredR=next.determinant*(next.inverse[2][0]*x+next.inverse[2][2]*z);
        m_initialFluxRoundtripError=std::max({m_initialFluxRoundtripError,
            std::abs(recoveredA-qa),std::abs(recoveredR-qr)});
        return std::array<double,2>{x,z};
    };
    const double foot=m_geometry.normalMap->lengthX;
    for (const auto &cell : m_flow.Elements())
    {
        if (cell.xr<=-foot || cell.xl>=foot)
            continue;
        const bool outer=cell.rl>=m_geometry.H-1e-12;
        const Vec u=cell.B*Local(flow,cell.u);
        const Vec w=cell.B*Local(flow,cell.u,count);
        const auto &matching=HOwner(cell);
        std::array<double,4> traceU{},traceW{};
        if (outer)
            for (int j=0;j<4;++j)
            {
                traceU[j]=matching.us[j]*flow[matching.u[j]];
                traceW[j]=matching.us[j]*flow[count+matching.u[j]];
            }
        Vec du(u.size()),dw(w.size());
        for (int i=0;i<u.size();++i)
        {
            const double a=cell.x[i], radius=cell.z[i]-m_geometry.d0;
            const auto old=source.At(a,radius,cell.xl,cell.xr,outer);
            const auto next=m_geometry.At(a,radius,cell.xl,cell.xr,outer);
            auto velocity=piola(old,next,u[i],w[i]);
            // Both chart owners give the same Piola normal flux through H.
            // Their transverse chart derivatives generally differ. Remove
            // only the resulting tangential jump from the outer seed.
            const double taperEnd=m_geometry.H+m_options.geometrySeedTaper;
            const double blend=outer ? 1-SmoothTransition(radius,m_geometry.H,taperEnd)[0] : 1.;
            if (outer && blend>0.)
            {
                const double H=m_geometry.H;
                const double width=matching.xr-matching.xl;
                const double coordinate=(a-matching.xl)/width;
                const double uH=CubicTrace(traceU,coordinate,width)[0];
                const double wH=CubicTrace(traceW,coordinate,width)[0];
                const auto newInner=m_geometry.At(a,H,cell.xl,cell.xr,false);
                const auto innerVelocity=piola(source.At(a,H,cell.xl,cell.xr,false),newInner,uH,wH);
                const auto outerVelocity=piola(source.At(a,H,cell.xl,cell.xr,true),
                                       m_geometry.At(a,H,cell.xl,cell.xr,true),uH,wH);
                const auto &t=newInner.tangent1;
                const double jump=t[0]*(outerVelocity[0]-innerVelocity[0])+t[2]*(outerVelocity[1]-innerVelocity[1]);
                m_initialTangentialTraceCorrection=std::max(m_initialTangentialTraceCorrection,std::abs(jump));
                velocity[0]-=jump*t[0];
                velocity[1]-=jump*t[2];
            }
            // Keep the previous Cartesian field in the outer bulk, then
            // project the represented profile with the target divergence.
            du[i]=blend*(velocity[0]-u[i]);
            dw[i]=blend*(velocity[1]-w[i]);
        }
        Add(correction,cell.u,cell.B.transpose()*cell.w.cwiseProduct(du));
        Add(correction,cell.u,cell.B.transpose()*cell.w.cwiseProduct(dw),count);
    }
    Eigen::SimplicialLDLT<Sparse> mass(m_flow.Core().Mass());
    Require(mass.info()==Eigen::Success,"geometry seed velocity mass factorization failed");
    result.head(count)+=mass.solve(correction.head(count));
    result.segment(count,count)+=mass.solve(correction.tail(count));
    Require(mass.info()==Eigen::Success && result.allFinite(),
            "geometry seed velocity mass projection failed");
    return result;
}

CoupledState Coupling::InitialState()const
{
    CoupledHistory h;h.dt=m_options.dt;
    h.current.flow=Vec::Zero(2*m_flow.VelocitySize()+m_flow.PressureSize());
    for(int component=0;component<3;++component){auto f=[&,component](double x,double z){
            if(m_mms){auto a=m_mms->Fields(NativeX(x),NativeRadius(z)+m_geometry.d0,0);return component==0?a.u.value:component==1?a.w.value:a.Pi.value;}
            if(component<2)return InitialVelocity(x,NativeRadius(z))[2*component];
            if(m_reference->Policy()==ReferencePolicy::ExactInnerFiniteNuV1)
                return m_reference->InletModifiedPressure(NativeRadius(z));
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
    h.current.epsilon=m_e.Project([&](double x,double z){
        const double r=NativeRadius(z);
        if (m_mms)
            return m_mms->OuterEpsilon(NativeX(x),r+m_geometry.d0,0).value;
        double level=r;
        if (m_geometry.normalMap && m_geometry.normalMap->amplitude!=0.)
        {
            const auto metric=m_geometry.At(x,r,x,x,true);
            const double blend=SmoothTransition(r,.35,4.)[0];
            level=(1-blend)*r+blend*(metric.point[2]-m_geometry.d0);
        }
        double E=m_reference->InletE(level);
        if (m_options.hillInitialEpsilonBump!=0.)
        {
            // This initial outer profile leaves the fixed inlet and the
            // dependent H trace intact. It supplies no equation source.
            // Use a separate entrance length for epsilon. The onset of the
            // slowdown's strain involves two horizontal derivatives of its
            // streamfunction, so the velocity window itself rises too slowly
            // to support that strain near the entrance.
            const double edge=std::min(x-m_left,m_right-x);
            const double fraction=std::min(1.,edge/2.);
            const double window=fraction*(2-fraction);
            const double y=r-m_geometry.H, rise=-std::expm1(-y/.035);
            const double fade=1-SmoothTransition(r,3.5,5.5)[0];
            E*=1+m_options.hillInitialEpsilonBump*window*rise*rise*fade;
            // The broad velocity transition can have appreciable strain
            // where the logarithmic inlet epsilon is already small. Give
            // that outer seed region its own smooth perturbation scale.
            const double midLayer=SmoothTransition(r,.8,2.)[0]
                *(1-SmoothTransition(r,7.,14.)[0]);
            E+=.1*m_reference->InletE(m_geometry.H)
                *m_options.hillInitialEpsilonBump*window*midLayer;
        }
        return E;
    });
    if(!m_mms && !m_geometry.normalMap)m_physicalFlowLift=m_flow.Lift(h.current.flow);
    if(m_geometry.normalMap)
    {
        h.current=DirichletState(h.current,0.);
        m_initialContinuityBefore=ContinuityNorm(h.current.flow);
        h.current.flow=ProjectContinuity(h.current.flow);
        m_initialContinuityAfter=ContinuityNorm(h.current.flow);
        Require(m_initialContinuityAfter<=1e-9,
                "normal-chart initialization does not satisfy represented continuity");
    }
    if (m_geometry.normalMap && !m_mms)
    {
        const auto &p=m_reference->Parameters();
        const auto K=InterfaceBounds(h.current);
        std::cout << std::setprecision(17) << "WJ-H27 initial normal-chart continuity "
                  << m_initialContinuityBefore << " -> " << m_initialContinuityAfter
                  << ", whole matching K/Kin [" << K[0]/p.KIn << ',' << K[1]/p.KIn
                  << ']' << std::endl;
    }
    return h.current;
}

CoupledHistory Coupling::Initialise()const
{
    return FreshHistory(InitialState());
}

CoupledState Coupling::LiftInitialState(const CoupledState&state)const
{
    try
    {
        auto lifted=DirichletState(state,0.);
        if(m_flow.Gauge())lifted.flow.tail(m_flow.PressureSize())-=m_flow.PressureOne()*
            (m_flow.PressureMean().dot(lifted.flow.tail(m_flow.PressureSize()))/m_volume);
        return lifted;
    }
    catch(const NativeInadmissibleState &error)
    {
        if(m_exportRejections)
            m_rejectedCandidate=DiagnosticCandidate{DiagnosticRole::RejectedInitialization,
                state,0,0.,error.what(),std::nullopt};
        throw;
    }
}

CoupledHistory Coupling::FreshHistory(const CoupledState&state)const
{
    CoupledHistory h;h.dt=m_options.dt;h.current=LiftInitialState(state);
    h.previous=h.current;CoupledTiming timing;
    Data data;
    try{data=Evaluate(h.current,0,true,timing);}
    catch(const NativeInadmissibleState&error){
        if(m_exportRejections)
            m_rejectedCandidate=DiagnosticCandidate{DiagnosticRole::RejectedInitialization,
                std::move(h.current),0,0.,error.what(),std::nullopt};
        throw;
    }
    h.advection=Advection(h.current,data.kRules);h.previousAdvection=h.advection;CheckHistory(h);return h;
}

CoupledHistory Coupling::ReadInitialFields(const std::string &path, bool allowCandidate,
                                            std::optional<double> geometrySeedHeight)
{
    auto history=FreshHistory(ImportInitialState(path,allowCandidate,geometrySeedHeight,nullptr));
    if(!allowCandidate && !geometrySeedHeight)
        Require(SteadyResidual(history.current).Maximum()<=m_options.nonlinearTolerance,
        "imported reference is not stationary for this native recipe; prepare a current reference first");
    return history;
}

CoupledState Coupling::PrepareState(const std::string &path,const std::string &seedSha256,
                                   std::optional<double> geometrySeedHeight,bool allowUncommittedGuess)
{
    Require(!allowUncommittedGuess||(!seedSha256.empty()&&!geometrySeedHeight),
            "uncommitted guess requires a digest-bound fixed-height diagnostic");
    return LiftInitialState(path.empty()?InitialState():
        ImportInitialState(path,true,geometrySeedHeight,&seedSha256,allowUncommittedGuess));
}

Vec Coupling::DownstreamVelocityPulse(double radius,double axialPeak) const
{
    const double foot=m_geometry.normalMap->lengthX;
    const int n=m_flow.VelocitySize();
    double lower=radius,upper=radius;
    for(const auto &cell:m_flow.Elements())
        if(cell.xl>=foot-1e-12)
        {
            if(std::abs(cell.rr-radius)<1e-12)lower=std::min(lower,cell.rl);
            if(std::abs(cell.rl-radius)<1e-12)upper=std::max(upper,cell.rr);
        }
    Require(lower>m_geometry.H&&lower<radius&&radius<upper&&upper<m_top&&m_right>foot,
            "velocity pulse requires a downstream outer radial face with two interior owners");
    const double span=m_right-foot;
    // The downstream flat chart is Cartesian. A C1 piecewise cubic
    // streamfunction gives an exactly representable curl on its Q3 cells.
    // The outlet derivative in x vanishes, preserving its transverse trace.
    const double scale=axialPeak*std::min(radius-lower,upper-radius)/1.5;
    auto sample=[&](double a,double z,int component) {
        const double r=z-m_geometry.d0;
        if(a<=foot||r<=lower||r>=upper)return 0.;
        const double x=(a-foot)/span;
        const bool rising=r<=radius;
        const double width=rising?radius-lower:upper-radius;
        const double y=rising?(r-lower)/width:(upper-r)/width;
        if(component==0)
            return scale*x*x*(3.-2.*x)*6.*y*(1.-y)/width*(rising?1.:-1.);
        return -scale*6.*x*(1.-x)/span*y*y*(3.-2.*y);
    };
    Vec delta=Vec::Zero(2*n+m_flow.PressureSize());
    for(int component=0;component<2;++component)
        delta.segment(component*n,n)=m_flow.Project(
            [&](double a,double z){return sample(a,z,component);},component);
    // Remove roundoff in the homogeneous essential and matching traces.
    delta-=m_flow.Lift(delta);
    for(const auto &cell:m_flow.Elements())
        if(std::abs(cell.rl-m_geometry.H)<1e-12)
            for(int j=0;j<4;++j)
                delta[cell.u[j]]=delta[n+cell.u[j]]=0.;
    const double continuity=ContinuityNorm(delta);
    Require(continuity<=1e-12,"represented velocity pulse does not satisfy discrete incompressibility");
    double peak=0.,gradient=0.;
    for(const auto &cell:m_flow.Elements())
        for(int component=0;component<2;++component)
        {
            const auto local=Local(delta,cell.u,component*n);
            peak=std::max(peak,(cell.B*local).cwiseAbs().maxCoeff());
            gradient=std::max({gradient,(cell.Dx*local).cwiseAbs().maxCoeff(),
                               (cell.Dz*local).cwiseAbs().maxCoeff()});
        }
    Require(peak>0.&&gradient>0.,"velocity pulse rounded to zero");
    std::cout<<std::setprecision(17)<<"WJ-H27 original velocity pulse signed axial peak "<<axialPeak
             <<" support ["<<foot<<','<<m_right<<"] x ["<<lower<<','<<upper
             <<"], sampled speed "<<peak<<", sampled gradient "<<gradient
             <<", scaled weak continuity "<<continuity<<std::endl;
    return delta;
}

CoupledHistory Coupling::ReleaseOriginalSteady(const std::string &path,const std::string &receipt,
    const std::string &receiptSha,int pulseCell,double relativeKPeak,
    double velocityPulseRadius,double velocityPulsePeak)
{
    Require(!m_mms&&m_geometry.normalMap&&m_options.advancement=="CoupledJFNKV1"&&
            m_options.UsesImplicitAdvection()&&std::isfinite(relativeKPeak)&&
            (relativeKPeak==0.?pulseCell<0:pulseCell>=0)&&std::isfinite(velocityPulsePeak)&&
            std::isfinite(velocityPulseRadius)&&(velocityPulsePeak==0.||velocityPulseRadius>m_geometry.H),
            "original steady release requires a fresh implicit physical hill trajectory");
    const InitialFields initial(path,false,true);
    initial.ValidateOriginalSteadyRelease(receipt,receiptSha,m_options.sourceIdentity,
        m_reference->Identity(),m_signature,m_options.quadratureX,m_options.quadratureR,
        m_geometry.normalMap->lengthX);
    const auto &metadata=initial.Metadata();
    auto state=PrepareState(path,metadata.packetSha256,std::nullopt);
    Require(!m_initialProjected&&StateHash(state)==metadata.stateSha256,
            "released coefficient inventory or rebuilt dependent lift differs");
    const auto result=AuditSteadyState(state);
    Require(result.residual.Maximum()<=std::min(1e-10,m_options.nonlinearTolerance),
            "fresh original steady residual failed before physical release");
    const auto released=FreshHistory(state);
    if(relativeKPeak==0.&&velocityPulsePeak==0.)return released;
    if(relativeKPeak!=0.)
    {
        // The Modified_A tensor interior mode (2,2) peaks at 1/16.
        Require(m_geometryIndex.count(pulseCell)>0,"unknown original steady k pulse cell");
        const int id=m_k.CellIds(pulseCell)[10];
        const auto &element=*std::find_if(initial.Elements().begin(),initial.Elements().end(),
            [pulseCell](const auto &cell){return cell.geometry==pulseCell;});
        Require(element.k[10]==id&&std::find(m_k.FreeIds().begin(),m_k.FreeIds().end(),id)!=m_k.FreeIds().end(),
                "pulse must use a free k interior mode");
        state.k[id]+=16.*element.ks[10]*m_reference->Parameters().KIn*relativeKPeak;
        Require(state.k[id]!=released.current.k[id],"original steady k pulse rounded to zero");
    }
    if(velocityPulsePeak!=0.)
        state.flow+=DownstreamVelocityPulse(velocityPulseRadius,velocityPulsePeak);
    return FreshHistory(state);
}

CoupledState Coupling::ImportInitialState(const std::string &path,bool allowCandidate,
    std::optional<double> geometrySeedHeight,const std::string *diagnosticSeedSha256,
    bool allowUncommittedGuess,const TransferObserver &observe)
{
    InitialFields initial(path,geometrySeedHeight.has_value(),diagnosticSeedSha256!=nullptr,allowUncommittedGuess);
    const bool acceptedTimeSeed=initial.Metadata().role==InitialFieldRole::AcceptedTimeState;
    const bool diagnosticSeed=initial.Metadata().role==InitialFieldRole::Diagnostic;
    InitialFieldGeometry geometry;
    geometry.d0=m_geometry.d0;geometry.H=m_geometry.H;geometry.h=m_geometry.h;geometry.left=m_left;geometry.right=m_right;geometry.top=m_top;geometry.curve=m_geometry.curve;
    if (m_geometry.normalMap)
        geometry.normalMap=m_geometry.normalMap->Signature();
    const std::string expectedSource=m_options.initialSourceIdentity.empty()?m_options.sourceIdentity:m_options.initialSourceIdentity;
    Require(expectedSource.size()==64,"initial-field source binding is missing");
    const auto &packet=m_options.geometrySeedPacketSha256;
    const bool pinnedGeometrySeed=geometrySeedHeight && !packet.empty();
    if (pinnedGeometrySeed)
        Require(packet==initial.Metadata().packetSha256,"geometry seed packet SHA256 differs");
    if (diagnosticSeedSha256 && !diagnosticSeedSha256->empty())
        Require(*diagnosticSeedSha256==initial.Metadata().packetSha256,
                "diagnostic seed packet SHA256 differs");
    Require(!diagnosticSeed || pinnedGeometrySeed || (diagnosticSeedSha256&&!diagnosticSeedSha256->empty()),
            "diagnostic seeds must be explicitly pinned by packet SHA256");
    Require(expectedSource==m_options.sourceIdentity || allowCandidate || pinnedGeometrySeed,
            "a different initial source requires a steady-preparation guess or an explicitly pinned geometry seed packet");
    const auto &seedReference=m_options.initialReferenceIdentity;
    Require(seedReference.empty()||(allowCandidate&&m_reference->Transport()==InnerKTransport::AuxiliaryHistoricalFixtureV1&&
            seedReference==m_reference->BoundaryIdentity()),"a different reference is allowed only as an explicitly bound v2 steady fixture guess");
    Geometry seedGeometry=m_geometry;
    if (geometrySeedHeight)
    {
        Require(m_geometry.normalMap && !m_options.manufactured &&
                std::isfinite(*geometrySeedHeight) && *geometrySeedHeight>=0. &&
                *geometrySeedHeight<=1. &&
                (initial.Metadata().qualifiedReference || acceptedTimeSeed || diagnosticSeed) &&
                initial.Metadata().outletPolicy==m_options.outletPolicy &&
                m_options.hillInitialSlowdown==0. && m_options.hillInitialEpsilonBump==0.,
                "geometry continuation requires a compatible accepted physical seed");
        seedGeometry.normalMap->amplitude=*geometrySeedHeight;
        seedGeometry.ValidateNormalMesh(m_left,m_right,m_top);
        geometry.normalMap=seedGeometry.normalMap->Signature();
    }
    const double seedRight=initial.GeometryData().right;
    const bool extendSeed=!geometrySeedHeight&&allowCandidate&&seedRight>=15.&&m_right>seedRight+1e-12;
    if(extendSeed){
        Require(m_zeroTransverseOutlet&&initial.Metadata().qualifiedReference&&
            initial.Metadata().outletPolicy==m_options.outletPolicy,
            "outlet extension requires a converged short-domain zero-transverse guess");
        geometry.right=seedRight;
    }
    initial.ValidateBindings(expectedSource,seedReference.empty()?m_reference->Identity():seedReference,
                             WallModelV2,m_options.manufactured,m_options.mmsVariant,geometry);
    Require(initial.Metadata().gauge==m_flow.Gauge(),"initial-field pressure-gauge binding differs");
    Require(!diagnosticSeed || initial.Metadata().outletPolicy==m_options.outletPolicy,
            "diagnostic seed outlet policy differs");
    Require(allowCandidate || initial.Metadata().qualifiedReference || acceptedTimeSeed,
            "advancement requires a qualified reference or geometry seed");
    std::vector<InitialFieldElement>inventory;
    for(const auto&e:m_flow.Elements()){
        InitialFieldElement cell;cell.geometry=e.geometry;cell.xl=e.xl;cell.xr=e.xr;cell.rl=e.rl;cell.rr=e.rr;cell.u=e.u;cell.us=e.us;cell.p=e.p;cell.ps=e.ps;
        for(int field:{2,3}){auto f=m_fields[field];auto map=f->GetLocalToGlobalMap();auto&ids=field==2?cell.k:cell.epsilon;auto&signs=field==2?cell.ks:cell.epsilonSigns;
            for(int i=0;i<f->GetExpSize();++i)if(f->GetExp(i)->GetGeom()->GetGlobalID()==e.geometry){for(int j=0;j<16;++j){int local=f->GetCoeff_Offset(i)+j;ids.push_back(map->GetLocalToGlobalMap(local));signs.push_back(map->GetLocalToGlobalSign(local));}break;}}
        inventory.push_back(std::move(cell));
    }
    const bool same=initial.SameInventory(inventory,{m_flow.VelocitySize(),m_flow.PressureSize(),m_k.Size(),m_e.Size()});
    const bool refined=geometrySeedHeight && !same && initial.CoversRefinement(inventory);
    Require(!geometrySeedHeight || same || refined,
            "geometry continuation requires the same signed chart inventory or a strict nested refinement");
    Require(!acceptedTimeSeed || refined || *geometrySeedHeight!=m_geometry.normalMap->amplitude,
            "same-height transient continuation requires a cold checkpoint");
    CoupledState state;
    if(same)state={initial.Flow(),initial.K(),initial.Epsilon()};
    else{
        // Extrapolate the accepted original exit section solely for the
        // initial projection. Current operators use the actual extended mesh.
        auto evaluate=[&](unsigned field,double x,double z){
            return initial.Evaluate(field,extendSeed?std::min(x,geometry.right):x,z);
        };
        state.flow=Vec::Zero(2*m_flow.VelocitySize()+m_flow.PressureSize());
        for(int component=0;component<3;++component){auto function=[&](double x,double z){return evaluate(component<2?component:4,x,z);};auto value=m_flow.Project(function,component);
            if(component<2)state.flow.segment(component*m_flow.VelocitySize(),m_flow.VelocitySize())=value;else state.flow.tail(m_flow.PressureSize())=value;}
        state.k=m_k.Project([&](double x,double z){return evaluate(2,x,z);});
        state.epsilon=m_e.Project([&](double x,double z){return evaluate(3,x,z);});
    }
    if (geometrySeedHeight)
    {
        if (observe) observe("nested_polynomial_transfer",state);
        std::optional<SeedGradientFit> fit;
        if (m_options.geometrySeedFootLength>0. || m_options.geometrySeedCollarLength>0.)
            fit=GeometrySeedGradientFit(state.flow,seedGeometry);
        state.flow=RemapVelocity(state.flow,seedGeometry);
        state=DirichletState(state,0.);
        if (observe) observe("essential_and_dependent_lifts",state);
        m_initialContinuityBefore=ContinuityNorm(state.flow);
        const int nv=2*m_flow.VelocitySize();
        const double gradientBefore=fit ? fit->Error(state.flow.head(nv)) : 0.;
        state.flow=ProjectContinuity(state.flow,fit ? &*fit : nullptr);
        if (observe) observe("constrained_continuity_projection",state);
        if (fit)
        {
            m_initialGradientFit=std::array<double,2>{gradientBefore,fit->Error(state.flow.head(nv))};
            std::cout<<std::setprecision(17)<<"WJ-H27 geometry seed weighted gradient RMS "
                     <<m_initialGradientFit->front()<<" -> "<<m_initialGradientFit->back()
                     <<std::endl;
        }
        m_initialContinuityAfter=ContinuityNorm(state.flow);
        Require(m_initialContinuityAfter<=1e-9,"geometry seed continuity projection failed");
        m_geometrySeedHeight=geometrySeedHeight;
        m_geometrySeedRefined=refined;
        m_geometrySeedStep=initial.Metadata().step;
        m_geometrySeedTime=initial.Metadata().time;
        std::cout<<std::setprecision(17)<<"WJ-H27 geometry seed "<<*geometrySeedHeight
                 <<" -> "<<m_geometry.normalMap->amplitude<<", weak continuity "
                 <<m_initialContinuityBefore<<" -> "<<m_initialContinuityAfter
                 <<", mapped flux roundtrip "<<m_initialFluxRoundtripError<<std::endl;
    }
    m_initialPacketSha=initial.Metadata().packetSha256;m_initialSourceStateSha=initial.Metadata().stateSha256;
    if(diagnosticSeedSha256)
    {
        const auto &metadata=initial.Metadata();
        m_importedMerit=metadata.meritController;
        m_importedMeritPolicy=metadata.meritPolicy;
        m_importedMeritCompatible=same&&!geometrySeedHeight&&!allowUncommittedGuess&&
            metadata.sourceIdentity==m_options.sourceIdentity&&metadata.signature==m_signature&&
            metadata.equationIdentity==RANSWallR3::Hash(m_signature+"|steady");
    }
    m_initialProjected=!same||geometrySeedHeight.has_value();
    m_initialCandidate=initial.Metadata().role==InitialFieldRole::Candidate || diagnosticSeed;
    return state;
}

CoupledHistory Coupling::PerturbInitialHistory(const CoupledHistory&base,const std::string&certificate,
                                               const std::string&sha256)
{
    Require(!m_mms&&!m_geometry.normalMap&&m_zeroTransverseOutlet&&
            m_reference->Policy()==ReferencePolicy::ExactInnerFiniteNuV1&&
            !m_initialPacketSha.empty()&&!m_initialProjected&&!m_initialCandidate&&
            base.step==0&&base.time==0,
            "Bf perturbation requires a qualified same-mesh physical initial state");
    std::set<double> xEdges,rEdges;
    for(const auto &cell:m_flow.Elements())
    {
        xEdges.insert(cell.xl);xEdges.insert(cell.xr);
        rEdges.insert(cell.rl);rEdges.insert(cell.rr);
    }
    const auto shape=Perturbation::ReadBfCertificate(*m_reference,{xEdges.begin(),xEdges.end()},
        {rEdges.begin(),rEdges.end()},certificate,sha256,m_options.sourceIdentity,StateHash(base.current));
    CoupledState state=base.current;
    for(int component=0;component<2;++component)
    {
        const auto delta=m_flow.Project([&](double x,double z){
            return shape.Velocity(x,z)[component==0?0:2];
        },component);
        state.flow.segment(component*m_flow.VelocitySize(),m_flow.VelocitySize())+=delta;
    }
    state.k+=m_k.Project([&](double x,double z){return shape.KIncrement(x,z).value;});
    state.epsilon+=m_e.Project([&](double x,double z){return shape.IndependentEpsilonIncrement(x,z).value;});
    // Rebuild essential traces from current K before creating either explicit history.
    return FreshHistory(state);
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

std::unique_ptr<SteadySparseInverse> Coupling::MakeSteadySparseInverse(
    const Vec &unknownScale,const Vec &rowScale)const
{
    const int n=m_flow.VelocitySize(),nf=m_flow.Size(),nk=m_k.FreeSize();
    const int size=unknownScale.size();
    std::array<std::vector<int>,2> velocity{std::vector<int>(n,-1),std::vector<int>(n,-1)};
    int row=0;
    for(int component=0;component<2;++component)
        for(int id:m_flow.FreeVelocityIds()[component])velocity[component][id]=row++;
    auto idsToFree=[](const Scalar &scalar,int offset) {
        std::vector<int> map(scalar.Size(),-1);
        for(int id:scalar.FreeIds())map[id]=offset++;
        return map;
    };
    const auto k=idsToFree(m_k,nf),epsilon=idsToFree(m_e,nf+nk);
    std::vector<std::vector<int>> columns(size);
    auto add=[](std::vector<int> &inventory,const auto &ids,const std::vector<int> &free) {
        for(int id:ids)if(free[id]>=0)inventory.push_back(free[id]);
    };
    auto connect=[&](const std::vector<int> &rows,const std::vector<int> &sources){
        for(int col:sources)columns[col].insert(columns[col].end(),rows.begin(),rows.end());
    };
    for(const auto &cell:m_flow.Elements())
    {
        std::vector<int> momentum;
        for(const auto &free:velocity)add(momentum,cell.u,free);
        auto rows=momentum;
        add(rows,m_k.CellIds(cell.geometry),k);
        const bool outer=m_e.Owns(cell.geometry);
        if(outer)add(rows,m_e.CellIds(cell.geometry),epsilon);
        auto sources=rows;
        // WJ production and stress use cell fields. The inner wall data and
        // the moving k quadrature use the matching K from the same column.
        // The only outer lift depending on K is the projected bottom trace.
        const auto &owner=HOwner(cell);
        if(!outer||std::abs(cell.rl-m_geometry.H)<1e-12)
            add(sources,m_hKTrace.at(owner.geometry).ids,k);
        connect(rows,sources);
        if(cell.rl<m_geometry.h)
        {
            // Sampled convection uses the four velocity trace modes. Use the
            // actual transported-k owner inventory for the sampled k action,
            // including rounded basis entries from the native trace mapping.
            std::vector<int> trace;
            const std::array<int,4> modes{owner.u[0],owner.u[1],owner.u[2],owner.u[3]};
            for(const auto &free:velocity)add(trace,modes,free);
            add(trace,m_k.CellIds(owner.geometry),k);
            connect(momentum,trace);
        }
    }
    const Sparse mixed=m_flow.System(0.,Sparse(2*n,2*n),1.);
    std::vector<Triplet> entries;
    for(int col=0;col<mixed.outerSize();++col)
        for(Sparse::InnerIterator it(mixed,col);it;++it)
            if(it.value()!=0.)entries.emplace_back(it.row(),col,it.value());
    Sparse algebraic(size,size);
    algebraic.setFromTriplets(entries.begin(),entries.end());
    return std::make_unique<SteadySparseInverse>(std::move(columns),std::move(algebraic),unknownScale,rowScale);
}

CoupledHistory Coupling::SteadyCandidateHistory(const Vec&unknowns)const{return FreshHistory(FromFree(unknowns,0));}

double Coupling::RegularizationShape(const NativeRule &rule,int q,const DiagnosticSteadyOptions &options)const
{
    // A polynomial of physical position, common on both sides of matching.
    // Added traction vanishes at the inlet, outlet and top.
    const auto &cell=Owner(rule.geometry);
    const double radius=rule.z[q]-m_geometry.d0;
    const auto point=m_geometry.normalMap?m_geometry.At(rule.x[q],radius,
        cell.xl,cell.xr,m_e.Owns(rule.geometry)).point:Vector{rule.x[q],0.,rule.z[q]};
    const double x=(point[0]-m_left)/(m_right-m_left);
    const double z=(point[2]-m_geometry.d0)/m_top;
    const double span=4.*x*(1.-x),height=1.-std::pow(z,8);
    const double ramp=options.regularizationFull>0.?
        SmoothTransition(radius,options.regularizationStart,options.regularizationFull)[0]:1.;
    return span*span*height*height*ramp;
}

Coupling::Data Coupling::Evaluate(const CoupledState&s,double time,PrincipalPolicy principal,
    CoupledTiming&timing,const DiagnosticSteadyOptions *diagnostic)const
{
    const double evaluateStart=Clock(),attributedBefore=timing.Sum();Data d;double start=Clock();const auto&p=m_reference->Parameters();
    const int subdivision=diagnostic?diagnostic->scalarSubdivisionDepth:0;
    auto kb=m_k.PolynomialBounds(s.k,subdivision),eb=m_e.PolynomialBounds(s.epsilon,subdivision);
    if(kb[0]<=0||eb[0]<=0){std::ostringstream error;error<<std::setprecision(17)<<"transported scalar positivity is not certified between quadrature points: k=["<<kb[0]<<','<<kb[1]<<"], epsilon=["<<eb[0]<<','<<eb[1]<<']';throw NativeInadmissibleState(error.str());}
    d.minK=kb[0];d.minE=std::min(eb[0],m_reference->Amplitude(p.KMinimum)/m_geometry.H);
    // The interface polynomial must remain in the qualified reference table
    // over its complete edge, rather than only at a chosen node set.
    const auto KBounds=InterfaceBounds(s);
    d.minE=std::min(eb[0],m_reference->Amplitude(KBounds[0])/m_geometry.H);
    timing.guards+=Clock()-start;start=Clock();
    const CutoffPartition partition{m_cutoffRange,[&](const Element&e){return CurrentCutoffXBreaks(e,s);}};
    d.kRules=m_rules.Cutoff([&](const Element&e,double x){return m_reference->Cutoff(InterfaceK(e,x,s));},&partition);
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
            const auto normal=r.WallNormal(q);
            if(flowPart)d.velocityInputs.back()[q]={k,E,aux,radius,G,normal};
            double began=Clock();auto a=AdaptAmended(k,E,G,normal,p.nu,aux,radius,m_geometry.H,m_options.backend,&m_constitutiveStats);timing.constitutive+=Clock()-began;
            Require(a.mu>0&&std::isfinite(a.production),"invalid current complete constitutive data");
            mu[q]=a.mu;rem[q]=a.remainder;viscous[q]=a.V;
            shear[q]=m_options.stressSplit=="PositiveShearV1"?(1-a.chi)*aux:0.;
            bulk[q]=a.mu-shear[q];
            // Positive planar split, with the full V minus this split kept
            // in the remainder. S is t tensor n + n tensor t.
            if(shear[q]!=0)
            {
                const Vector t{normal[2],0.,-normal[0]};
                Tensor S{};
                double projected=0.;
                for (int i : {0,2}) for (int j : {0,2})
                {
                    S[i][j]=t[i]*normal[j]+normal[i]*t[j];
                    projected+=S[i][j]*G[i][j];
                }
                for (int i : {0,2}) for (int j : {0,2})
                    rem[q][i][j]+=shear[q]*(G[i][j]+G[j][i]-S[i][j]*projected);
            }
            Dk[q]=m_reference->KDiffusivity(k,E,aux,outer);
            De[q]=outer?p.nu+p.Cd*k*k/(p.sigmaE*E):0.;
            sourceK[q]=a.production-E;sourceE[q]=outer?1.44*a.production*E/k:0.;reactionE[q]=outer?1.82*E/k:0.;
            d.minK=std::min(d.minK,k);d.minE=std::min(d.minE,E);
            if(flowPart){coefficient.push_back(bulk[q]);coefficient.push_back(shear[q]);kCoefficient.push_back(Dk[q]);if(outer){eCoefficient.push_back(De[q]);eCoefficient.push_back(reactionE[q]);}}
            // Covariance is checked at every current assembly point. Principal
            // checks include the continuum of supported planar directions.
            const bool checkPrincipal=principal!=PrincipalPolicy::Skip&&flowPart;
            began=Clock();auto guard=CheckPlanarState(k,E,G,p.nu,aux,radius,m_geometry.H,checkPrincipal,&a,m_options.backend,&m_constitutiveStats,normal);
            d.minCovariance=std::min(d.minCovariance,guard.covarianceOverK);
            if(checkPrincipal)
            {
                d.minPrincipal=std::min(d.minPrincipal,guard.principalOverMu);
                if(diagnostic)
                {
                    const double added=diagnostic->eta>0.?diagnostic->eta*diagnostic->regularizationViscosity*
                        RegularizationShape(r,q,*diagnostic):0.;
                    d.minModifiedPrincipal=std::min(d.minModifiedPrincipal,guard.principalOverMu+added/a.mu);
                }
            }
            if(guard.covarianceOverK<=1e-10||(principal==PrincipalPolicy::Enforce&&flowPart&&guard.principalOverMu<=1e-8)){
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
            m_reference->Parameters().nu,input.aux,input.r,m_geometry.H,m_options.backend,&m_constitutiveStats,input.normal));
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
        std::array<Mat,2> projected{g[1],g[0]};
        if(m_geometry.normalMap && wa.squaredNorm()!=0)
        {
            for(int q=0;q<wa.size();++q)
            {
                const auto n=r.WallNormal(q);
                const std::array<double,2> t{n[2],-n[0]};
                for(int i=0;i<2;++i)
                    projected[i].row(q)=(t[i]*n[0]+n[2*i]*t[0])*r.Dx.row(q)
                        +(t[i]*n[2]+n[2*i]*t[1])*r.Dz.row(q);
            }
        }
        for(int i=0;i<2;++i)for(int j=0;j<2;++j){Mat k=g[j].transpose()*w.asDiagonal()*g[i];if(i==j)k+=lap;
            k.noalias()+=projected[i].transpose()*wa.asDiagonal()*projected[j];
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

Vec Coupling::MomentumAdvection(const CoupledState&s,Vec *sampledAction)const
{
    const int n=m_flow.VelocitySize();Vec load=Vec::Zero(2*n+m_flow.PressureSize());
    if(sampledAction)*sampledAction=Vec::Zero(load.size());
    for(const auto&rp:m_rules.Static()){
        const auto&r=*rp;Vec u=Local(s.flow,r.u),v=Local(s.flow,r.u,n),up=r.B*u,vp=r.B*v;
        std::array<Vec,2> a{up.cwiseProduct(r.Dx*u)+vp.cwiseProduct(r.Dz*u),
                            up.cwiseProduct(r.Dx*v)+vp.cwiseProduct(r.Dz*v)};
        const auto &cell=Owner(r.geometry);
        if(cell.rl<m_geometry.h)
        {
            std::array<Vec,2> sampled;
            if(sampledAction)for(auto &v:sampled)v=Vec::Zero(r.w.size());
            for(int q=0;q<r.w.size();++q) if(r.z[q]-m_geometry.d0<m_geometry.h)
            {
                const auto h=InterfaceVelocity(cell,r.x[q],s);
                const auto normal=r.WallNormal(q);
                const double amplitude=2*h[0]*h[1];
                a[0][q]-=normal[2]*amplitude;
                a[1][q]+=normal[0]*amplitude;
                if(sampledAction){
                    sampled[0][q]=-normal[2]*amplitude;
                    sampled[1][q]=normal[0]*amplitude;
                }
            }
            if(sampledAction)
                for(int i=0;i<2;++i)
                    Add(*sampledAction,r.u,r.B.transpose()*r.w.asDiagonal()*sampled[i],i*n);
        }
        for(int i=0;i<2;++i)
            Add(load,r.u,r.B.transpose()*r.w.asDiagonal()*a[i],i*n);
    }
    return load;
}
CoupledState Coupling::Advection(const CoupledState &state, const Rules &kRules,
                                CoupledStepTimingDetail *timing,Vec *sampledMomentum) const
{
    const double began = Clock();
    Vec momentum = MomentumAdvection(state,sampledMomentum);
    const double scalarBegan = Clock();
    Vec k = m_k.Advection(kRules, state.flow, state.k);
    Vec epsilon = m_e.Advection(m_rules.Static(), state.flow, state.epsilon);
    if (timing)
    {
        timing->momentumAdvection += scalarBegan - began;
        timing->scalarAdvection += Clock() - scalarBegan;
    }
    return {std::move(momentum), std::move(k), std::move(epsilon)};
}

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
        if(Owner(r.geometry).rl<m_geometry.h){auto h=HRule(r);
            Vec weights=r.w;for(int q=0;q<weights.size();++q)if(r.z[q]-m_geometry.d0>=m_geometry.h)weights[q]=0;
            std::array<Mat,2> values{h.B,h.B}, derivatives{h.Dx,h.Dx};
            std::array<Vec,2> tangents{Vec::Ones(weights.size()),Vec::Zero(weights.size())};
            for(int q=0;q<weights.size();++q)
            {
                const auto normal=h.WallNormal(q);
                tangents[0][q]=normal[2];
                tangents[1][q]=-normal[0];
                const auto ds=(normal[2]*h.Dx.row(q)-normal[0]*h.Dz.row(q)).eval();
                for (int j=0;j<2;++j)
                {
                    values[j].row(q)=tangents[j][q]*h.B.row(q);
                    derivatives[j].row(q)=tangents[j][q]*ds-
                        h.MatchingCurvature(q)*normal[2*j]*h.B.row(q);
                }
            }
            const std::array<Vec,2> local{Local(s.flow,h.u),Local(s.flow,h.u,n)};
            const Vec ut=values[0]*local[0]+values[1]*local[1];
            const Vec dut=derivatives[0]*local[0]+derivatives[1]*local[1];
            for(int i=0;i<2;++i) for(int j=0;j<2;++j)
            {
                const Mat sampled=-2*r.B.transpose()*weights.cwiseProduct(tangents[i]).asDiagonal()*
                    (dut.asDiagonal()*values[j]+ut.asDiagonal()*derivatives[j]);
                for(int a=0;a<16;++a)for(int b=0;b<16;++b)
                    triplets.emplace_back(i*n+r.u[a],j*n+h.u[b],sampled(a,b));
            }
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
        const CutoffPartition exactPartition{m_cutoffRange,[&](const Element&e){return MMSCutoffXBreaks(e,time);}};
        auto exactRules=m_rules.Cutoff([&](const Element&,double x){return m_reference->Cutoff(m_mms->InterfaceK(x,time).value);},&exactPartition);
        std::vector<Vec>ksource;
        for(const auto&rp:exactRules){Vec f(rp->w.size());for(int q=0;q<f.size();++q){f[q]=m_mms->Forcing(rp->x[q],rp->z[q],time).k;if(freezeManufactured)f[q]-=m_mms->Fields(rp->x[q],rp->z[q],time).k.dt;}ksource.push_back(std::move(f));}
        load.k=m_k.Source(exactRules,ksource);
        load.k+=m_k.NaturalFlux([&](double x,double z,const Vector&normal){return m_mms->KOutwardFlux(x,z,time,normal);},m_rules,
            [&](const Element&,double x){return m_reference->Cutoff(m_mms->InterfaceK(NativeX(x),time).value);});
        load.epsilon+=m_e.NaturalFlux([&](double x,double z,const Vector&normal){
            // H is Dirichlet only: its reaction is diagnosed from the solve.
            return NativeRadius(z)==m_geometry.H?0.:m_mms->EpsilonOutwardFlux(x,z,time,normal);
        },m_rules);
    }else{
        // The application has no manufactured/balance forcing. Body force
        // is exactly zero for this pinned neutral physical case.
        const double topTraction=m_reference->TopTangentialTraction();
        if(m_options.checkOutletActions)m_outletLoad=Vec::Zero(2*n);
        for(const auto&[index,side]:m_boundaryFaces)if(side==1||side==3){auto r=m_rules.Face(m_flow.Elements()[index],side,20);Vec tx(r.w.size()),tz(tx.size());
            for(int q=0;q<tx.size();++q){if(side==1){auto t=m_reference->OutletTraction(r.z[q]-m_geometry.d0);tx[q]=t[0];tz[q]=m_zeroTransverseOutlet?0.:t[1];}else{tx[q]=topTraction;tz[q]=0;}}
            const Vec normal=r.B.transpose()*r.w.asDiagonal()*tx;
            const Vec transverse=r.B.transpose()*r.w.asDiagonal()*tz;
            Add(load.flow,r.u,normal);Add(load.flow,r.u,transverse,n);
            if(side==1&&m_options.checkOutletActions){
                Add(m_outletLoad,r.u,normal);Add(m_outletLoad,r.u,transverse,n);
            }
        }
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

Vec Coupling::FixedKDiffusionPart(const Rules &rules,const CoupledState &state,const Vec *direction)const
{
    Require(std::string(m_reference->TransportName())=="ScalarCdV2",
            "fixed k diffusion tangent uses the ScalarCdV2 law");
    const auto &p=m_reference->Parameters();
    Vec dk,de;
    if(direction)
    {
        dk=m_k.Expand(direction->segment(m_flow.Size(),m_k.FreeSize()),Vec::Zero(m_k.Size()));
        de=m_e.Expand(direction->tail(m_e.FreeSize()),DependentTraceDirection(state,dk));
    }
    std::vector<Vec> diffusion,tangentD;
    for(const auto &entry:rules)
    {
        const auto &q=*entry;
        const Vec k=m_k.Evaluate(q,state.k);
        const bool outer=m_e.Owns(q.geometry);
        Vec epsilon,deltaEpsilon;
        if(outer)
        {
            epsilon=m_e.Evaluate(q,state.epsilon);
            if(direction)deltaEpsilon=m_e.Evaluate(q,de);
        }
        else
        {
            const auto &owner=HOwner(Owner(q.geometry));
            const auto &trace=m_hKTrace.at(owner.geometry);
            std::array<double,4> dkTrace{};
            if(direction)
                for(int i=0;i<4;++i)dkTrace[i]=trace.signs[i]*dk[trace.ids[i]];
            const double width=owner.xr-owner.xl;
            epsilon.resize(q.w.size());
            if(direction)deltaEpsilon.resize(q.w.size());
            for(int i=0;i<q.w.size();++i)
            {
                const double K=InterfaceK(owner,q.x[i],state);
                const auto e=m_reference->DependentEpsilon(K,q.z[i]-m_geometry.d0);
                epsilon[i]=e.value;
                if(direction)
                    deltaEpsilon[i]=e.dK*CubicTrace(dkTrace,(q.x[i]-owner.xl)/width,width)[0];
            }
        }
        auto &D=diffusion.emplace_back(k.size());
        for(int i=0;i<k.size();++i)D[i]=m_reference->KDiffusivity(k[i],epsilon[i],0.,outer);
        if(direction)
        {
            const Vec deltaK=m_k.Evaluate(q,dk);
            tangentD.emplace_back(p.Cd/p.sigmaK*(2.*k.array()*deltaK.array()/epsilon.array()-
                k.array().square()*deltaEpsilon.array()/epsilon.array().square()));
        }
    }
    if(direction)
        return m_k.Restrict(m_k.FixedRuleDiffusionTangent(rules,diffusion,tangentD,state.k,dk));
    return m_k.Restrict(m_k.DiffusionAction(rules,diffusion,state.k));
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
        if (m_options.picardAdvection == "Implicit")
        {
            const auto advection = Advection(s, d.kRules);
            action.flow += advection.flow;
            action.k += advection.k;
            action.epsilon += advection.epsilon;
        }
    }
    action.flow-=external.flow;action.k-=external.k;action.epsilon-=external.epsilon;
    auto raw=m_flow.Restrict(action.flow);if(m_flow.Gauge())raw[raw.size()-1]=m_flow.PressureMean().dot(s.flow.tail(m_flow.PressureSize()));
    const auto flow=m_flow.RawResidualNorm(raw);
    const double kr=m_k.FreeResidualNorm(m_k.Restrict(action.k)),er=m_e.FreeResidualNorm(m_e.Restrict(action.epsilon));
    timing.residual+=Clock()-start;return {flow.momentum,flow.continuity,flow.gauge,kr,er,d.traceError};
}

CoupledStepInfo Coupling::Advance(CoupledHistory&h,bool injectRetry)
{
    Require(m_reference->Transport()==InnerKTransport::ScalarV2,"inner-k diagnostic transport cannot advance a trajectory");
    if(m_options.advancement=="CoupledJFNKV1")return AdvanceNewton(h,injectRetry);
    const bool fixedBudget = m_options.picardStop == "FixedCorrectors";
    const bool implicitAdvection = m_options.picardAdvection == "Implicit";
    const bool velocityTangent = m_options.picardStressLinearization == "VelocityTangent";
    const int maximumSweeps = fixedBudget ? m_options.picardOuterCorrectors : m_options.maximumIterations;
    const double begin=Clock();CheckHistory(h);m_constitutiveStats={};m_lastInputHistory=h;m_hasInputHistory=true;CoupledStepInfo info;const double time=(h.step+1)*h.dt,alpha=(h.step==0?1.:1.5)/h.dt;
    m_lastAttempt={};auto saveAttempt=[&]{info.seconds=Clock()-begin;info.constitutiveStats=m_constitutiveStats;info.timing.other=info.seconds-info.timing.Sum();m_lastAttempt=info;};
    if (fixedBudget)
    {
        ResetPreconditioners();
    }
    auto history=Combination(h.current,h.step==0?1/h.dt:2/h.dt,h.previous,h.step==0?0:-.5/h.dt);
    Loads external=Forcing(time,info.timing);
    if (!implicitAdvection)
    {
        const auto adv=Combination(h.advection,h.step==0?1:2,h.previousAdvection,h.step==0?0:-1);
        external.flow-=adv.flow;external.k-=adv.k;external.epsilon-=adv.epsilon;
    }
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
        try{for(int iteration=0;iteration<maximumSweeps;++iteration){
            ++info.nonlinearIterations;
            const double iterationStart=Clock();const Vec input=m_options.andersonDepth?PackFree(s):Vec();
            double start=Clock();
            const double attributedBefore = info.timing.Sum();
            const double drift=std::max(Drift(data.coefficients,m_preconditionerCoefficients),alphaDrift(m_flowFactorAlpha));info.coefficientDrift=std::max(info.coefficientDrift,drift);
            Sparse K=Strain(data);
            Vec flowLoad=base.flow+RemainderLoad(data);
            const int velocitySize=2*m_flow.VelocitySize();
            flowLoad.head(velocitySize)+=m_flow.KnownScalarKLoad(s.k);
            if (velocityTangent)
            {
                Sparse tangent = VelocityTangent(data, info.timing);
                const auto velocity = s.flow.head(velocitySize);
                // The full current stress carries the same affine remainder
                // when its frozen strain is replaced by its velocity tangent.
                flowLoad.head(velocitySize) += tangent*velocity - K*velocity;
                K = std::move(tangent);
            }
            if (implicitAdvection)
            {
                const Sparse convectionJacobian = MomentumTransportJacobian(s);
                const Vec convection = MomentumAdvection(s);
                // The complete convection A(u) is quadratic, including the H trace.
                // Its affine linearization is A'(u) u_next - A(u).
                if (m_options.checkCurrentActions)
                {
                    const Vec defect = convectionJacobian * s.flow.head(velocitySize) -
                                       2 * convection.head(velocitySize);
                    Require(defect.cwiseAbs().maxCoeff() <= 1e-10,
                            "assembled momentum convection Jacobian disagrees with its quadratic action");
                }
                K += convectionJacobian;
                flowLoad += convection;
            }
            const Sparse A=m_flow.System(alpha,K,1.);
            Vec lift=m_flow.Lift(s.flow),rhs=m_flow.LiftedRhs(alpha,K,1.,flowLoad,lift),free=m_flow.Restrict(s.flow);
            info.timing.assembly += (Clock()-start) - (info.timing.Sum()-attributedBefore);
            if(injectRetry&&attempt==0&&iteration==0)throw NativeLinearFailure("qualification injected recoverable Krylov failure");
            const double flowStart = Clock(), setupBefore = info.flowWork.setup;
            SolveInfo fi;
            try
            {
                fi = fixedBudget ? m_flow.Correct(A, rhs, free, info.flowWork,
                                                m_options.nonlinearTolerance, !implicitAdvection)
                                 : m_flow.Solve(A, rhs, free, drift, m_options.linearTolerance,
                                                &info.flowWork);
            }
            catch (...)
            {
                const double setup = info.flowWork.setup - setupBefore;
                info.timing.rebuild += setup;
                info.timing.flow += Clock() - flowStart - setup;
                info.flowIterations = info.flowWork.krylovIterations;
                throw;
            }
            s.flow=m_flow.Expand(free,lift);info.flowIterations=info.flowWork.krylovIterations;info.rebuilds+=fi.rebuilds;info.timing.flow+=fi.solveSeconds+fi.residualSeconds;info.timing.rebuild+=fi.factorSeconds;
            if(fi.rebuilds){m_preconditionerCoefficients=data.coefficients;m_flowFactorAlpha=alpha;}
            // Refresh P, E_w, current positive scalar diffusivities and all
            // stresses after the flow update, BEFORE the k solve.
            data=Evaluate(s,time,false,info.timing);start=Clock();Sparse Kk=m_k.Diffusion(data.kRules,data.Dk);
            if (implicitAdvection)
            {
                Kk += m_k.TransportMatrix(data.kRules, s.flow);
            }
            Vec klift=m_k.Lift(s.k),krhs=m_k.LiftedRhs(alpha,Kk,base.k+KVolumeLoad(s,data),klift),kf=m_k.Restrict(s.k);auto Ak=m_k.System(alpha,Kk);info.timing.assembly+=Clock()-start;
            const double kd=std::max(Drift(data.kCoefficients,m_kPreconditionerCoefficients),alphaDrift(m_kFactorAlpha));info.coefficientDrift=std::max(info.coefficientDrift,kd);
            auto ki=m_k.Solve(Ak,krhs,kf,kd,m_options.linearTolerance);s.k=m_k.Expand(kf,klift);info.kIterations+=ki.iterations;info.rebuilds+=ki.rebuilds;info.timing.scalar+=ki.solveSeconds+ki.residualSeconds;info.timing.rebuild+=ki.factorSeconds;
            info.timingDetail.kFactor += ki.factorSeconds;
            if(ki.rebuilds){m_kPreconditionerCoefficients=data.kCoefficients;m_kFactorAlpha=alpha;}
            // K(H) has changed: replace the dependent outer epsilon trace and
            // recompute the entire closure before solving outer epsilon.
            start=Clock();s=DirichletState(s,time);info.timing.wall+=Clock()-start;data=Evaluate(s,time,false,info.timing);start=Clock();Sparse Ke=m_e.Diffusion(m_rules.Static(),data.De,data.reactionE);
            if (implicitAdvection)
            {
                Ke += m_e.TransportMatrix(m_rules.Static(), s.flow);
            }
            Vec elift=m_e.Lift(s.epsilon),erhs=m_e.LiftedRhs(alpha,Ke,base.epsilon+EVolumeLoad(s,data),elift),ef=m_e.Restrict(s.epsilon);auto Ae=m_e.System(alpha,Ke);info.timing.assembly+=Clock()-start;
            const double ed=std::max(Drift(data.eCoefficients,m_ePreconditionerCoefficients),alphaDrift(m_eFactorAlpha));info.coefficientDrift=std::max(info.coefficientDrift,ed);
            auto ei=m_e.Solve(Ae,erhs,ef,ed,m_options.linearTolerance);s.epsilon=m_e.Expand(ef,elift);info.epsilonIterations+=ei.iterations;info.rebuilds+=ei.rebuilds;info.timing.scalar+=ei.solveSeconds+ei.residualSeconds;info.timing.rebuild+=ei.factorSeconds;
            info.timingDetail.epsilonFactor += ei.factorSeconds;
            if(ei.rebuilds){m_ePreconditionerCoefficients=data.eCoefficients;m_eFactorAlpha=alpha;}
            // Refresh all current equation residuals for diagnostics and acceptance.
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
            if (fixedBudget ? iteration + 1 == maximumSweeps
                            : info.residual.Maximum() <= m_options.nonlinearTolerance)
            {
                auto final=Evaluate(s,time,true,info.timing);info.residual=Residual(s,final,external,alpha,time,info.timing,&h);
                if (fixedBudget)
                {
                    const auto &res = info.residual;
                    const double tolerance = m_options.nonlinearTolerance;
                    const std::array<double, 6> rows{res.momentum, res.continuity, res.gauge,
                                                     res.k, res.epsilon, res.trace};
                    RequireState(std::all_of(rows.begin(), rows.end(),
                                             [](double value) { return std::isfinite(value); }) &&
                                 res.continuity <= tolerance && res.gauge <= tolerance &&
                                 res.trace <= tolerance,
                                 "fixed Picard budget failed final continuity, gauge or trace checks");
                }
                else
                {
                    Require(info.residual.Maximum() <= m_options.nonlinearTolerance,
                            "final refreshed coupled residual changed acceptance");
                }
                info.minimumK=final.minK;info.minimumEpsilon=final.minE;info.minimumCovarianceOverK=final.minCovariance;info.minimumPrincipalOverMu=final.minPrincipal;
                const double advStart=Clock();auto nextAdvection=Advection(s,final.kRules);info.timing.assembly+=Clock()-advStart;
                CoupledHistory next;
                next.step = h.step + 1;
                next.time = time;
                next.dt = h.dt;
                next.current = std::move(s);
                next.previous = h.current;
                next.advection = std::move(nextAdvection);
                next.previousAdvection = h.advection;
                CheckHistory(next);
                h = std::move(next);
                info.acceptedByFixedBudget = fixedBudget;
                complete=true;break;
            }
        }}catch(const NativeLinearFailure&e){info.retryReason=e.what();recoverable=true;if(fixedBudget||attempt==1){saveAttempt();throw;}}
        catch(...){saveAttempt();throw;}
        if(complete){saveAttempt();return info;}
        if(!recoverable)break;
    }
    saveAttempt();std::ostringstream error;error<<std::setprecision(17)<<"current coupled iteration did not converge; accepted history preserved, residuals=["<<info.residual.momentum<<','<<info.residual.continuity<<','<<info.residual.gauge<<','<<info.residual.k<<','<<info.residual.epsilon<<','<<info.residual.trace<<"], nonlinear iterations="<<info.nonlinearIterations<<", linear iterations=["<<info.flowIterations<<','<<info.kIterations<<','<<info.epsilonIterations<<']';throw std::runtime_error(error.str());
}

CoupledResidual Coupling::SteadyResidual(const CoupledState&s,double time)
{
    return AuditSteadyState(s,time).residual;
}

Vec Coupling::SteadyRows(const Loads &action,const CoupledState &state,double gaugeMultiplier)const
{
    const int nf=m_flow.Size(),nk=m_k.FreeSize(),ne=m_e.FreeSize();
    Vec raw(nf+nk+ne);
    raw.head(nf)=m_flow.Restrict(action.flow);
    if(m_flow.Gauge())
    {
        raw.segment(m_flow.PressureOffset(),m_flow.PressureSize())+=m_flow.PressureMean()*gaugeMultiplier;
        raw[nf-1]=m_flow.PressureMean().dot(state.flow.tail(m_flow.PressureSize()));
    }
    raw.segment(nf,nk)=m_k.Restrict(action.k);
    raw.tail(ne)=m_e.Restrict(action.epsilon);
    return raw;
}

CoupledResidual Coupling::SteadyRowNorms(const Vec &raw)const
{
    const int nf=m_flow.Size(),nk=m_k.FreeSize(),ne=m_e.FreeSize();
    Require(raw.size()==nf+nk+ne&&raw.allFinite(),"invalid native steady residual inventory");
    const auto flow=m_flow.RawResidualNorm(raw.head(nf));
    return {flow.momentum,flow.continuity,flow.gauge,m_k.FreeResidualNorm(raw.segment(nf,nk)),
            m_e.FreeResidualNorm(raw.tail(ne)),0.};
}

void Coupling::WriteOriginalTransferAudit(const std::filesystem::path &path,
    const CoupledState &root,const std::string &rootPacketSha,
    const std::string &parent,const std::string &parentPacketSha,
    const SteadyKrylovOptions &options)
{
    Require(m_geometry.normalMap&&!m_mms&&!m_flow.Gauge()&&options.evaluateInitialOnly&&
            rootPacketSha.size()==64&&parentPacketSha.size()==64,
            "original transfer audit requires two bound steady fields on a physical traction-outlet chart");
    Require(!std::filesystem::exists(path),"original transfer audit is immutable");
    const InitialFields parentField(parent,true,true,false);
    Require(parentField.Metadata().originalDiagnosticSteady&&
            parentField.Metadata().packetSha256==parentPacketSha,
            "transfer audit requires a bound original steady parent");
    const double began=Clock();
    std::vector<std::pair<std::string,CoupledState>> stages;
    auto observe=[&](const char *name,const CoupledState &state){stages.emplace_back(name,state);};
    const auto transferred=LiftInitialState(ImportInitialState(parent,true,
        m_geometry.normalMap->amplitude,&parentPacketSha,false,observe));
    Require(m_geometrySeedRefined,"transfer audit requires a strict nested child");
    observe("prepared",transferred);
    observe("solved",root);
    const Vec x=PackFree(root),p=PackFree(transferred),delta=x-p;
    Require(StateHash(FromFree(x,0.))==StateHash(root)&&StateHash(FromFree(p,0.))==StateHash(transferred),
            "transfer audit changed essential or dependent lift coefficients");
    const int nf=m_flow.Size(),nk=m_k.FreeSize(),ne=m_e.FreeSize(),nv=2*m_flow.VelocitySize();
    Vec rowScale(nf+nk+ne);
    rowScale.head(nf)=m_flow.ResidualRowScale();
    rowScale.segment(nf,nk)=m_k.ResidualRowScale();
    rowScale.tail(ne)=m_e.ResidualRowScale();
    CoupledTiming timing;
    const auto forcing=Forcing(0.,timing,true);
    using TermRows=std::map<std::string,Vec>;
    auto evaluate=[&](const Vec &unknown,TermRows *terms=nullptr)->Vec {
        Require(Clock()-began<options.maximumSeconds,"original transfer audit exceeded its wall budget");
        try
        {
            const auto s=FromFree(unknown,0.);
            const auto data=Evaluate(s,0.,PrincipalPolicy::Skip,timing);
            const auto current=CurrentActions(s,data);
            Vec sampledConvection;
            const auto adv=Advection(s,data.kRules,nullptr,terms?&sampledConvection:nullptr);
            const Loads total{current.flow+adv.flow-forcing.flow,
                current.k+adv.k-forcing.k,current.epsilon+adv.epsilon-forcing.epsilon};
            const Vec rows=SteadyRows(total,s,0.);
            if(!terms)return rows;

            auto flow=[&](const char *name,const Vec &v,int offset=0) {
                Vec local=Vec::Zero(current.flow.size()),raw=Vec::Zero(rows.size());
                local.segment(offset,v.size())=v;
                raw.head(nf)=m_flow.Restrict(local);
                terms->emplace(name,std::move(raw));
            };
            const auto pressure=s.flow.tail(m_flow.PressureSize());
            const Vec ordinary=m_flow.Core().Gradient()*pressure;
            const Vec sampledPressure=-m_flow.Core().AWF()*pressure;
            const Vec sampledK=-m_flow.KnownScalarKLoad(s.k);
            flow("ordinary_pressure",ordinary);
            flow("sampled_pressure",sampledPressure);
            flow("sampled_k",sampledK);
            flow("sampled_convection",sampledConvection);
            flow("constitutive_momentum",current.flow.head(nv)-ordinary-sampledPressure-sampledK);
            flow("local_momentum_advection",adv.flow-sampledConvection);
            flow("continuity",current.flow.tail(m_flow.PressureSize()),nv);
            flow("exterior_traction",-forcing.flow);
            Vec scalar=Vec::Zero(rows.size());
            scalar.segment(nf,nk)=rows.segment(nf,nk);
            terms->emplace("k_equation",scalar);
            scalar.segment(nf,nk).setZero();
            scalar.tail(ne)=rows.tail(ne);
            terms->emplace("epsilon_equation",std::move(scalar));
            return rows;
        }
        catch(const NativeInadmissibleState &error){throw SteadyInadmissibleTrial(error.what());}
    };
    auto feasible=[&](const Vec &unknown){
        try
        {
            const auto s=FromFree(unknown,0.);
            if(m_k.PolynomialBounds(s.k)[0]<=0.||m_e.PolynomialBounds(s.epsilon)[0]<=0.)
                return SteadyTrialValidity{false,"transported scalar polynomial is not strictly positive"};
            return SteadyTrialValidity{};
        }
        catch(const NativeInadmissibleState &error){return SteadyTrialValidity{false,error.what()};}
    };
    TermRows preparedTerms,rootTerms;
    const Vec fp=evaluate(p,&preparedTerms),fx=evaluate(x,&rootTerms);
    Require(SteadyRowNorms(fx).Maximum()<=options.acceptanceNorm,
            "transfer audit root does not satisfy the current original equations");
    SteadyKrylovCallbacks callbacks;
    callbacks.ValidTrial=feasible;
    callbacks.Residual=[&](const Vec &v){return evaluate(v);};
    callbacks.TrueNorm=[&](const Vec &v){return SteadyRowNorms(v).Maximum();};
    callbacks.Precondition=[](const Vec &v){return v;}; // A directional probe never applies an inverse.
    const Vec scales=UnknownScales();
    const SteadyKrylov probe(std::move(callbacks),scales,rowScale,options);
    const auto derivative=probe.ProbeDirection(x,delta,{2.,1.,.5,.25});
    Require(derivative.completed&&derivative.actions.size()==4&&
            derivative.finiteDifferences.back().method=="central",
            "cannot resolve current transfer direction: "+derivative.reason);
    const auto &j=derivative.actions.back();
    const double h=derivative.finiteDifferences.back().parameterStep;
    TermRows plus,minus;
    (void)evaluate(x+h*delta,&plus);
    (void)evaluate(x-h*delta,&minus);
    const Vec remainder=fp-fx+j;
    const Vec spacingChange=j-derivative.actions[derivative.actions.size()-2];
    Vec exactContinuity=Vec::Zero(x.size());
    exactContinuity.segment(m_flow.PressureOffset(),m_flow.PressureSize())=
        m_flow.Core().Divergence()*(root.flow.head(nv)-transferred.flow.head(nv));

    std::ofstream out(path);
    namespace J=RANSWallSolveJson;
    auto writeNorms=[&](const Vec &raw){
        const auto n=SteadyRowNorms(raw);
        out<<'['<<n.momentum<<','<<n.continuity<<','<<n.k<<','<<n.epsilon<<']';
    };
    out<<std::setprecision(17)<<"{\"schema\":\"WJH27NativeOriginalTransferAuditV1\",\"physical_history\":false"
       <<",\"source_identity\":"<<J::Quote(m_options.sourceIdentity)<<",\"signature\":"<<J::Quote(m_signature)
       <<",\"reference_identity\":"<<J::Quote(m_reference->Identity())<<",\"eta\":0"
       <<",\"root_packet_sha256\":"<<J::Quote(rootPacketSha)
       <<",\"parent_packet_sha256\":"<<J::Quote(parentPacketSha)
       <<",\"parent_coefficient_sha256\":"<<J::Quote(parentField.Metadata().stateSha256)
       <<",\"row_order\":[\"momentum\",\"continuity\",\"k\",\"epsilon\"]"
       <<",\"free_rows\":["<<m_flow.PressureOffset()<<','<<m_flow.PressureSize()<<','<<nk<<','<<ne<<']'
       <<",\"scaled_delta_infinity\":"<<delta.cwiseQuotient(scales).lpNorm<Eigen::Infinity>()
       <<",\"F_prepared\":";writeNorms(fp);
    out<<",\"F_root\":";writeNorms(fx);
    out<<",\"J_root_delta\":";writeNorms(j);
    out<<",\"remainder\":";writeNorms(remainder);
    out<<",\"halving_change\":";writeNorms(spacingChange);
    out<<",\"exact_linear_continuity_action\":";writeNorms(exactContinuity);
    Vec continuityClosure=exactContinuity;
    continuityClosure.segment(m_flow.PressureOffset(),m_flow.PressureSize())+=
        (fp-fx).segment(m_flow.PressureOffset(),m_flow.PressureSize());
    out<<",\"exact_linear_continuity_closure\":";writeNorms(continuityClosure);
    out<<",\"spacing_records\":[";
    bool comma=false;
    for(const auto &record:derivative.finiteDifferences)
    {
        if(comma)out<<',';
        comma=true;
        out<<"{\"method\":"<<J::Quote(record.method)
           <<",\"parameter_step\":"<<record.parameterStep<<",\"relative_step\":"<<record.relativeStep<<'}';
    }
    out<<"],\"successive_scaled_relative_change\":"<<J::Numbers(derivative.successiveRelativeScaledChanges)
       <<",\"terms\":{";
    comma=false;
    Vec fpSum=Vec::Zero(x.size()),fxSum=fpSum,jSum=fpSum;
    for(const auto &[name,term]:preparedTerms)
    {
        const Vec &other=rootTerms.at(name);
        const Vec action=(plus.at(name)-minus.at(name))/(2.*h);
        fpSum+=term;fxSum+=other;jSum+=action;
        if(comma)out<<',';
        comma=true;
        out<<J::Quote(name)<<":{\"prepared_minus_root\":";writeNorms(term-other);
        out<<",\"action\":";writeNorms(action);
        out<<",\"remainder\":";writeNorms(term-other+action);out<<'}';
    }
    out<<"},\"partition_closure\":{\"F_prepared\":";writeNorms(fpSum-fp);
    out<<",\"F_root\":";writeNorms(fxSum-fx);
    out<<",\"J_root_delta\":";writeNorms(jSum-j);
    out<<"},\"root_sampled_residual_rows\":{";
    comma=false;
    for(const char *name:{"pressure","k","convection"})
    {
        if(comma)out<<',';
        comma=true;
        out<<J::Quote(name)<<':';
        VectorJSON(out,rootTerms.at(std::string("sampled_")+name).head(m_flow.PressureOffset()));
    }
    out<<"},\"stages\":{";
    comma=false;
    for(const auto &[name,state]:stages)
    {
        if(comma)out<<',';
        comma=true;
        out<<J::Quote(name)<<":{\"state_sha256\":"<<J::Quote(StateHash(state))
           <<",\"weak_continuity\":"<<ContinuityNorm(state.flow)<<",\"coefficients\":{\"flow\":";
        VectorJSON(out,state.flow);out<<",\"k\":";VectorJSON(out,state.k);
        out<<",\"epsilon\":";VectorJSON(out,state.epsilon);out<<"}}";
    }
    out<<"},\"audit_seconds\":"<<Clock()-began<<"}\n"<<std::flush;
    Require(bool(out),"original transfer audit write failed");
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
    CheckHistory(history);m_lastInputHistory=history;m_hasInputHistory=true;
    Require(history.step==0&&history.time==0,"steady preparation requires a fresh step-zero state");
    return SearchSteady(history.current,options,final,preconditionerMass,&history,nullptr);
}

SteadyKrylovResult Coupling::SolveDiagnosticSteady(const CoupledState &seed,
    const SteadyKrylovOptions &options,CoupledStepInfo &final,double preconditionerMass,
    const DiagnosticSteadyOptions &diagnostic)
{
    for(double value:{diagnostic.pseudoTime,diagnostic.correctionDiffusion,diagnostic.eta,
                      diagnostic.regularizationViscosity,diagnostic.regularizationStart,
                      diagnostic.regularizationFull})
        Require(std::isfinite(value)&&value>=0.,"invalid diagnostic steady parameter");
    Require(diagnostic.pseudoTime==0.||std::isfinite(1./diagnostic.pseudoTime),
            "diagnostic reciprocal pseudo-time must be representable");
    Require(diagnostic.eta==0.||(diagnostic.regularizationViscosity>0.&&!m_mms&&
            std::isfinite(diagnostic.eta*diagnostic.regularizationViscosity)),
            "diagnostic diffusion requires a finite positive prescribed physical-case coefficient");
    Require((diagnostic.regularizationStart==0.&&diagnostic.regularizationFull==0.)||
            (m_geometry.normalMap&&diagnostic.regularizationStart<diagnostic.regularizationFull&&
             diagnostic.regularizationFull<=m_geometry.H),
            "diagnostic radial ramp must lie inside the normal wall layer");
    Require(diagnostic.scalarSubdivisionDepth>=0&&diagnostic.scalarSubdivisionDepth<=6,
            "invalid diagnostic scalar Bernstein depth");
    m_diagnosticOptions=diagnostic;
    m_diagnosticCorrections.clear();
    m_sparseSteadyWork.clear();
    m_steadyBaseObservations.clear();
    m_steadyActionAudits.clear();
    m_steadyEpsilonRemainders.clear();
    m_steadyPatternSeconds=m_steadyAuditSeconds=0.;
    m_steadyEpsilonSeconds=0.;
    m_steadyEpsilonEvaluations=0;
    m_steadyAuditResiduals=0;
    m_steadyDiagnosticCandidate.reset();
    m_steadyUncommitted.reset();
    auto controlled=options;
    Require(options.merit||options.evaluateInitialOnly||!m_importedMerit,
            "an original-merit checkpoint needs explicit continuation or re-anchoring");
    m_currentMerit.reset();
    m_meritPolicy.clear();
    if(controlled.merit)
    {
        Require(!options.evaluateInitialOnly&&diagnostic.eta==0.&&diagnostic.correctionDiffusion==0.&&
                diagnostic.pseudoTime==0.,"original-merit policy requires an original correction and explicit mass");
        auto &policy=*controlled.merit;
        std::ostringstream binding;
        binding<<std::setprecision(17)<<SteadyMeritState::Identity;
        for(double value:{policy.initialMass,policy.reanchorMass,options.forcingInitial,options.forcingMinimum,
                          options.forcingMaximum,options.forcingGamma,options.forcingExponent,
                          options.fdRelativeStep,options.fdMinimumRelativeStep,options.armijo})
            binding<<'|'<<value;
        if(policy.initialTrialFraction>0.)binding<<"|trial-fraction:"<<policy.initialTrialFraction;
        if(policy.curvature)
            binding<<"|coupled-trial-v1:"<<policy.curvature->maximumProducts
                   <<'|'<<policy.curvature->maximumScaledSize<<'|'<<policy.curvature->budgetFraction;
        if(diagnostic.fixedKTangent)binding<<"|fixed-rule-k-tangent-v1";
        m_meritPolicy=RANSWallR3::Hash(binding.str());
        if(diagnostic.meritEntry=="Resume")
        {
            Require(m_importedMerit&&m_importedMeritCompatible&&m_importedMeritPolicy==m_meritPolicy&&
                    StateHash(FromFree(PackFree(seed),0.))==m_initialSourceStateSha,
                    "original-merit resume requires the exact accepted coefficients, policy and equation/source bindings");
            policy.restored=m_importedMerit;
            Require(diagnostic.actionAuditPrefix.empty()||
                    (policy.restored->usedMass==0.&&policy.restored->nextMass==0.),
                    "unshifted action audit cannot restore a shifted controller");
        }
        else
            Require(diagnostic.meritEntry=="Reanchor"||!m_importedMerit,
                    "original-merit restart requires explicit Resume or Reanchor");
    }
    return SearchSteady(seed,controlled,final,preconditionerMass,nullptr,&diagnostic);
}

SteadyKrylovResult Coupling::SearchSteady(const CoupledState &seed,const SteadyKrylovOptions &options,
    CoupledStepInfo &final,double preconditionerMass,CoupledHistory *physicalHistory,
    const DiagnosticSteadyOptions *diagnostic)
{
    const double began=Clock();
    Require(std::isfinite(preconditionerMass)&&preconditionerMass>=0,"invalid steady preconditioner mass");
    Require(options.acceptanceNorm>0&&options.acceptanceNorm<=m_options.nonlinearTolerance,
        "steady acceptance must meet the declared current coupled tolerance");
    ResetPreconditioners();struct ColdExit{Coupling&owner;~ColdExit(){owner.ResetPreconditioners();}}cold{*this};
    m_constitutiveStats={};final={};CoupledTiming&timing=final.timing;
    const int nf=m_flow.Size(),nk=m_k.FreeSize(),ne=m_e.FreeSize(),size=nf+nk+ne;
    const int nv=2*m_flow.VelocitySize();
    const auto forcing=Forcing(0,timing,true);
    Sparse positiveDiffusion,regularization;
    const bool shifted=diagnostic&&(options.merit||diagnostic->pseudoTime>0.||diagnostic->correctionDiffusion>0.);
    const bool modified=diagnostic&&diagnostic->eta>0.;
    if(shifted||modified)
    {
        auto positiveForm=[&](double coefficient,bool useProfile) {
            Data weights;
            for(const auto&cell:m_rules.Static())
            {
                weights.bulk.emplace_back(Vec::Constant(cell->w.size(),coefficient));
                weights.shear.emplace_back(Vec::Zero(cell->w.size()));
                if(useProfile)
                    for(int q=0;q<cell->w.size();++q)
                        weights.bulk.back()[q]*=RegularizationShape(*cell,q,*diagnostic);
            }
            return Strain(weights);
        };
        if(diagnostic->correctionDiffusion>0.)positiveDiffusion=positiveForm(1.,false);
        if(modified)regularization=positiveForm(diagnostic->eta*diagnostic->regularizationViscosity,true);
    }
    auto pack=[&](const CoupledState&s){return PackFree(s);};
    auto unpack=[&](const Vec&x){return FromFree(x,0);};
    auto norms=[&](const Vec&raw){return SteadyRowNorms(raw);};
    Vec validatedUnknowns;std::shared_ptr<const Data>validatedData;
    auto currentData=[&](const Vec&x,const CoupledState&s){
        if(validatedData&&x.size()==validatedUnknowns.size()&&(x.array()==validatedUnknowns.array()).all())return validatedData;
        return std::shared_ptr<const Data>(
            std::make_shared<Data>(Evaluate(s,0,PrincipalPolicy::Skip,timing,diagnostic)));
    };
    auto equations=[&](const Vec&x){
        auto s=unpack(x);auto owned=currentData(x,s);const auto&d=*owned;const double start=Clock();auto adv=Advection(s,d.kRules);auto action=CurrentActions(s,d);
        action.flow+=adv.flow-forcing.flow;action.k+=adv.k-forcing.k;action.epsilon+=adv.epsilon-forcing.epsilon;
        if(modified)action.flow.head(nv).noalias()+=regularization*s.flow.head(nv);
        const Vec raw=SteadyRows(action,s,m_flow.Gauge()?x[nf-1]:0.);
        if(m_options.checkCurrentActions){
            Sparse K=Strain(d);
            if(modified)K+=regularization;
            const Sparse Kk=m_k.Diffusion(d.kRules,d.Dk),Ke=m_e.Diffusion(m_rules.Static(),d.De,d.reactionE);
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
    const int npOffset=m_flow.PressureOffset();
    const std::array<int,5> offsets{0,npOffset,nf,nf+nk,size};
    if(diagnostic)
        for(int block=0;block<4;++block)
        {
            const int offset=offsets[block],n=offsets[block+1]-offset;
            const auto x=sx.segment(offset,n),r=sr.segment(offset,n);
            m_steadyScaleExtrema[block]={x.minCoeff(),x.maxCoeff(),r.minCoeff(),r.maxCoeff()};
        }
    const bool sparse=diagnostic&&diagnostic->sparseInverse;
    const double patternStart=Clock();
    auto coupledInverse=sparse?MakeSteadySparseInverse(sx,sr):nullptr;
    if(diagnostic)
    {
        m_steadyPatternSeconds=sparse?Clock()-patternStart:0.;
        m_steadyBaseSha256=StateHash(unpack(pack(seed)));
    }
    CoupledState pcState, auditedFactorState;Sparse epsilonPCFull;
    double initialBaseNorm=0.,correctionMass=0.,correctionDiffusion=0.,auditedFactorMass=0.;
    auto blockNorms=[&](const Vec &raw) {
        const auto r=norms(raw);
        return std::array<double,6>{r.momentum,r.continuity,r.gauge,r.k,r.epsilon,r.trace};
    };
    auto deadline=[&]{
        if(Clock()-began>=options.maximumSeconds)
            throw SteadyCorrectionFailure("DEADLINE: bounded coupled inverse preparation exhausted");
    };
    SteadyKrylovCallbacks callbacks;
    callbacks.ValidTrial=[&](const Vec&x){
        try{auto s=unpack(x);
            const int depth=diagnostic?diagnostic->scalarSubdivisionDepth:0;
            const auto k=m_k.PolynomialBounds(s.k,depth),e=m_e.PolynomialBounds(s.epsilon,depth);
            if(k[0]<=0||e[0]<=0)return SteadyTrialValidity{false,"transported scalar polynomial is not strictly positive"};
            return SteadyTrialValidity{};
        }catch(const NativeInadmissibleState&error){return SteadyTrialValidity{false,error.what()};}
    };
    callbacks.ValidateNonlinearState=[&](const Vec&x){
        try{
            auto s=unpack(x);
            auto full=std::make_shared<Data>(Evaluate(s,0,
                diagnostic?PrincipalPolicy::Skip:PrincipalPolicy::Enforce,timing,diagnostic));
            validatedUnknowns=x;validatedData=std::move(full);
        }catch(const NativeInadmissibleState&error){throw SteadyInadmissibleTrial(error.what());}
    };
    callbacks.Residual=[&](const Vec&x)->Vec{try{return equations(x);}catch(const NativeInadmissibleState&error){throw SteadyInadmissibleTrial(error.what());}};
    callbacks.TrueNorm=[&](const Vec&raw){return norms(raw).Maximum();};
    double meritBaseNorm=0.;
    Rules fixedKRules;
    if(diagnostic&&diagnostic->fixedKTangent)
    {
        // The original cutoff factory marks every cell intersecting the
        // certified range as moving, for any admissible current trace.
        for(const auto &rule:m_rules.Static())
        {
            const auto &cell=Owner(rule->geometry);
            if(cell.rl>=.4*m_geometry.H||cell.rr<=0.||
               cell.rr<m_cutoffRange[0]||cell.rl>m_cutoffRange[1])
                fixedKRules.push_back(rule);
        }
        std::cout<<"WJ-H27 exact fixed-rule k diffusion tangent "<<fixedKRules.size()
                 <<" cells, "<<m_rules.Static().size()-fixedKRules.size()<<" moving cells kept in FD"<<std::endl;
        callbacks.ExactPartResidual=[&](const Vec &x) {
            Vec raw=Vec::Zero(size);
            raw.segment(nf,nk)=FixedKDiffusionPart(fixedKRules,unpack(x));
            return raw;
        };
        callbacks.ExactPartAction=[&](const Vec &x,const Vec &direction) {
            Vec raw=Vec::Zero(size);
            raw.segment(nf,nk)=FixedKDiffusionPart(fixedKRules,unpack(x),&direction);
            return raw;
        };
    }
    callbacks.PrepareCorrection=[&](const Vec&x,const Vec&raw){
        pcState=unpack(x);
        if(options.merit)meritBaseNorm=norms(raw).Maximum();
        else if(shifted)
        {
            const double norm=norms(raw).Maximum();
            if(initialBaseNorm==0.)initialBaseNorm=norm;
            // Grow pseudo-time and relax correction diffusion from measured
            // reduction in the unshifted equation norm. Freeze this base.
            const double ratio=std::clamp(norm/initialBaseNorm,1e-8,1.);
            correctionMass=diagnostic->pseudoTime>0.?ratio/diagnostic->pseudoTime:0.;
            correctionDiffusion=std::sqrt(ratio)*diagnostic->correctionDiffusion;
            m_diagnosticCorrections.push_back({norm,correctionMass,correctionDiffusion});
        }
    };
    if(options.merit)
    {
        const std::array<int,7> freeOffsets{0,int(m_flow.FreeVelocityIds()[0].size()),npOffset,
            nf-(m_flow.Gauge()?1:0),nf,nf+nk,size};
        callbacks.ScaledUnknownBlocks=[freeOffsets](const Vec &x) {
            std::array<double,6> result{};
            for(int i=0;i<6;++i)
                if(const int count=freeOffsets[i+1]-freeOffsets[i];count>0)
                    result[i]=x.segment(freeOffsets[i],count).cwiseAbs().maxCoeff();
            return result;
        };
        callbacks.SetCorrectionMass=[&](double mass) {
            correctionMass=mass;
            m_diagnosticCorrections.push_back({meritBaseNorm,mass,0.});
        };
        callbacks.ObserveMeritState=[&](const SteadyMeritState &state) {
            m_currentMerit=state;
        };
        callbacks.ObserveMeritDirection=[](const SteadyMeritRecord &r) {
            std::cout<<std::setprecision(12)<<"WJ-H27 original merit base "<<r.base
                     <<" proposal "<<r.proposal<<" mass "<<r.mass<<" forcing "<<r.forcing
                     <<" shifted error "<<r.linearNorm<<" action spacing uncertainty "<<r.actionUncertainty
                     <<" verified fraction "<<r.step<<" "<<r.decision<<std::endl;
        };
    }
    auto differentialShift=[&](const CoupledState&delta,double mass,double diffusion) {
        Vec action=Vec::Zero(size);
        Vec momentum=mass*m_flow.Core().MassAction(delta.flow);
        if(diffusion>0.)
            momentum.head(nv).noalias()+=diffusion*(positiveDiffusion*delta.flow.head(nv));
        action.head(nf)=m_flow.Restrict(momentum);
        if(mass>0.)
        {
            action.segment(nf,nk)=mass*m_k.Restrict(m_k.Mass()*delta.k);
            action.tail(ne)=mass*m_e.Restrict(m_e.Mass()*delta.epsilon);
        }
        return action;
    };
    auto shiftAt=[&](const CoupledState &base,const Vec&dx,double mass,double diffusion) {
        CoupledState delta;
        delta.flow=m_flow.Expand(dx.head(nf),Vec::Zero(nv+m_flow.PressureSize()));
        if(mass>0.)
        {
            delta.k=m_k.Expand(dx.segment(nf,nk),Vec::Zero(m_k.Size()));
            delta.epsilon=m_e.Expand(dx.tail(ne),DependentTraceDirection(base,delta.k));
        }
        return differentialShift(delta,mass,diffusion);
    };
    auto applyShift=[&](const Vec&dx,double mass,double diffusion) {
        return shiftAt(pcState,dx,mass,diffusion);
    };
    if(shifted)callbacks.CorrectionAction=[&](const Vec&,const Vec&dx) {
        return applyShift(dx,correctionMass,correctionDiffusion);
    };
    if(diagnostic&&!diagnostic->actionAuditPrefix.empty())
    {
        std::set<int> matching;
        for(const auto &[id,trace]:m_hKTrace)
            matching.insert(trace.ids.begin(),trace.ids.end());
        std::vector<int> traceRows;
        for(int i=0;i<nk;++i)
            if(matching.count(m_k.FreeIds()[i]))traceRows.push_back(nf+i);
        callbacks.TraceDirection=[traceRows](const Vec &d) {
            Vec trace=Vec::Zero(d.size());
            for(int row:traceRows)trace[row]=d[row];
            return trace;
        };
        D::Sha256 ordering;
        for(const auto *ids:{&m_flow.FreeVelocityIds()[0],&m_flow.FreeVelocityIds()[1],&m_k.FreeIds(),&m_e.FreeIds()})
            ordering.Update(ids->data(),ids->size()*sizeof(int));
        const auto orderingSha=ordering.Finish();
        callbacks.ObserveActionAudit=[&,orderingSha,sequence=std::size_t{0}](const SteadyKrylovActionAudit &observed) mutable {
            const double start=Clock();
            Require(correctionMass==0.&&correctionDiffusion==0.,"action audit requires an unshifted current correction");
            namespace J=RANSWallSolveJson;
            std::optional<SteadyKrylovActionAudit> localized;
            auto raw=[&](const char *name)->Vec {
                for(const auto &[key,array]:observed.arrays)
                    if(key==name)return array.col(0);
                throw std::runtime_error(std::string("missing action audit array ")+name);
            };
            auto spacing=[&](const char *purpose) {
                for(const auto &s:observed.stencils)
                    if(s.purpose==purpose&&s.method=="central")return s.parameterStep;
                throw std::runtime_error(std::string("term audit requires central stencil ")+purpose);
            };
            if(diagnostic->auditScalarAction&&observed.stage=="seed-direction")
            {
                localized=observed;
                auto &audit=*localized;
                const Vec x=raw("base"),direction=raw("direction");
                struct ScalarActionSample
                {
                    Vec k;
                    Rules rules;
                    std::vector<Vec> diffusion;
                    std::array<Vec,7> residual;
                };
                const auto terms=[&](const Vec &at) {
                    const auto s=unpack(at);
                    const auto d=currentData(at,s); // Rebuild trace and moving rules on every side.
                    const auto &fixed=m_rules.Static();
                    const Vec ed=m_e.DiffusionAction(fixed,d->De,s.epsilon);
                    return ScalarActionSample{s.k,d->kRules,d->Dk,{
                        m_k.Restrict(m_k.DiffusionAction(d->kRules,d->Dk,s.k)),
                        m_k.Restrict(-KVolumeLoad(s,*d)),
                        m_k.Restrict(m_k.Advection(d->kRules,s.flow,s.k)),
                        m_e.Restrict(ed),
                        m_e.Restrict(m_e.DiffusionAction(fixed,d->De,s.epsilon,d->reactionE)-ed),
                        m_e.Restrict(-EVolumeLoad(s,*d)),
                        m_e.Restrict(m_e.Advection(fixed,s.flow,s.epsilon))}};
                };
                const double h=spacing("spacing_0"),half=spacing("spacing_1");
                const auto plus=terms(x+h*direction),minus=terms(x-h*direction);
                const auto finePlus=terms(x+half*direction),fineMinus=terms(x-half*direction);
                const std::array<const char *,7> labels{"k_diffusion","k_source","k_advection",
                    "epsilon_diffusion","epsilon_destruction","epsilon_source","epsilon_advection"};
                Vec difference=Vec::Zero(size);
                for(std::size_t i=0;i<labels.size();++i)
                {
                    const bool k=i<3;
                    Vec term=Vec::Zero(size);
                    term.segment(k?nf:nf+nk,k?nk:ne)=
                        (finePlus.residual[i]-fineMinus.residual[i])/(2.*half)-
                        (plus.residual[i]-minus.residual[i])/(2.*h);
                    difference+=term;
                    audit.arrays.emplace_back(std::string("raw_spacing_term_")+labels[i],std::move(term));
                }
                const Vec measured=raw("raw_spacing_change_1");
                Vec closure=Vec::Zero(size);
                closure.tail(nk+ne)=(difference-measured).tail(nk+ne);
                audit.arrays.emplace_back("raw_scalar_term_closure",std::move(closure));
                auto stableDifference=[&](const ScalarActionSample &a,const ScalarActionSample &b,double step)->Vec {
                    return m_k.Restrict(m_k.FixedRuleDiffusionDifferenceCorrection(
                        a.rules,a.diffusion,a.k,b.rules,b.diffusion,b.k))/(2.*step);
                };
                const Vec coarseCorrection=stableDifference(plus,minus,h);
                const Vec fineCorrection=stableDifference(finePlus,fineMinus,half);
                Vec stableCoarse=raw("raw_action_0"),stableFine=raw("raw_action_1");
                stableCoarse.segment(nf,nk)+=coarseCorrection;
                stableFine.segment(nf,nk)+=fineCorrection;
                audit.arrays.emplace_back("raw_stable_fixed_k_spacing_change",stableFine-stableCoarse);
                audit.arrays.emplace_back("raw_stable_fixed_k_action_0",stableCoarse);
                audit.arrays.emplace_back("raw_stable_fixed_k_action_1",stableFine);
                audit.arrays.emplace_back("raw_stable_fixed_k_model_0",raw("raw_base_residual")+stableCoarse);
                audit.arrays.emplace_back("raw_stable_fixed_k_taylor_0",raw("raw_difference_0")-h*stableFine);
                audit.arrays.emplace_back("raw_stable_fixed_k_taylor_1",raw("raw_difference_1")-half*stableFine);
                const auto base=unpack(x);
                const auto d=currentData(x,base);
                Rules fixed;
                std::array<std::vector<Vec>,5> diffusion;
                for(std::size_t i=0;i<m_rules.Static().size();++i)
                {
                    const auto &rule=m_rules.Static()[i];
                    if(d->kRules[i]!=rule||plus.rules[i]!=rule||minus.rules[i]!=rule||
                       finePlus.rules[i]!=rule||fineMinus.rules[i]!=rule)continue;
                    fixed.push_back(rule);
                    diffusion[0].push_back(d->Dk[i]);diffusion[1].push_back(plus.diffusion[i]);
                    diffusion[2].push_back(minus.diffusion[i]);diffusion[3].push_back(finePlus.diffusion[i]);
                    diffusion[4].push_back(fineMinus.diffusion[i]);
                }
                const Vec tangent=FixedKDiffusionPart(fixed,base,&direction);
                const auto fixedAction=[&](const Vec &k,std::size_t i) {
                    return m_k.Restrict(m_k.DiffusionAction(fixed,diffusion[i],k));
                };
                const Vec fixedBase=fixedAction(base.k,0),fixedPlus=fixedAction(plus.k,1),
                    fixedMinus=fixedAction(minus.k,2),fixedFinePlus=fixedAction(finePlus.k,3),
                    fixedFineMinus=fixedAction(fineMinus.k,4);
                Vec tangentCoarse=raw("raw_action_0"),tangentFine=raw("raw_action_1"),tangentRows=Vec::Zero(size);
                tangentCoarse.segment(nf,nk)+=tangent-(fixedPlus-fixedMinus)/(2.*h);
                tangentFine.segment(nf,nk)+=tangent-(fixedFinePlus-fixedFineMinus)/(2.*half);
                tangentRows.segment(nf,nk)=tangent;
                audit.arrays.emplace_back("raw_exact_fixed_k_tangent",tangentRows);
                audit.arrays.emplace_back("raw_tangent_replaced_spacing_change",tangentFine-tangentCoarse);
                audit.arrays.emplace_back("raw_tangent_replaced_action_0",tangentCoarse);
                audit.arrays.emplace_back("raw_tangent_replaced_model_0",raw("raw_base_residual")+tangentCoarse);
                tangentRows.segment(nf,nk)=fixedPlus-fixedBase-h*tangent;
                audit.arrays.emplace_back("raw_exact_fixed_k_taylor_0",tangentRows);
                tangentRows.segment(nf,nk)=fixedFinePlus-fixedBase-half*tangent;
                audit.arrays.emplace_back("raw_exact_fixed_k_taylor_1",tangentRows);
                audit.measurements.emplace_back("tangent_fixed_rules",fixed.size());
                audit.measurements.emplace_back("unchanged_moving_rules",m_rules.Static().size()-fixed.size());
                for(bool k:{true,false})
                {
                    const int offset=k?nf:nf+nk,count=k?nk:ne;
                    Eigen::Index row;
                    measured.segment(offset,count).cwiseProduct(sr.segment(offset,count)).
                        cwiseAbs().maxCoeff(&row);
                    const auto &scalar=k?m_k:m_e;
                    const int id=scalar.FreeIds()[row];
                    const auto label=k?"k":"epsilon";
                    audit.measurements.emplace_back(std::string("peak_")+label+"_free_row",offset+row);
                    audit.measurements.emplace_back(std::string("peak_")+label+"_scalar_id",id);
                    std::ostringstream owners;
                    owners<<std::setprecision(17);
                    const auto &elements=m_flow.Elements();
                    for(std::size_t i=0;i<elements.size();++i)
                    {
                        const auto &cell=elements[i];
                        if(!scalar.Owns(cell.geometry))continue;
                        const auto &ids=scalar.CellIds(cell.geometry);
                        if(std::find(ids.begin(),ids.end(),id)!=ids.end())
                            owners<<cell.geometry<<" ["<<cell.xl<<','<<cell.xr<<','<<cell.rl<<','<<cell.rr
                                  <<"] moving rule "<<(d->kRules[i]!=m_rules.Static()[i])<<' ';
                    }
                    audit.notes.emplace_back(std::string("peak_")+label+"_owners",owners.str());
                }
                audit.measurements.emplace_back("scalar_attribution_seconds",Clock()-start);
            }
            if(diagnostic->auditMomentumAction&&observed.stage=="seed-direction")
            {
                localized=observed;
                auto &audit=*localized;
                const Vec x=raw("base"),direction=raw("direction");
                const auto terms=[&](const Vec &at) {
                    const auto s=unpack(at);
                    const auto d=currentData(at,s);
                    const Vec pressure=s.flow.tail(m_flow.PressureSize());
                    const Vec gradient=m_flow.Core().Gradient()*pressure;
                    const Vec sampledPressure=-m_flow.Core().AWF()*pressure;
                    const Vec sampledK=-m_flow.KnownScalarKLoad(s.k);
                    const auto current=CurrentActions(s,*d);
                    Vec sampledConvection;
                    const Vec transport=MomentumAdvection(s,&sampledConvection);
                    return std::array<Vec,6>{gradient,sampledPressure,sampledK,
                        current.flow.head(nv)-gradient-sampledPressure-sampledK,
                        transport.head(nv)-sampledConvection.head(nv),sampledConvection.head(nv)};
                };
                const double h=spacing("spacing_0"),half=spacing("spacing_1");
                const auto plus=terms(x+h*direction),minus=terms(x-h*direction),
                    finePlus=terms(x+half*direction),fineMinus=terms(x-half*direction);
                const std::array<const char *,6> labels{"pressure","sampled_pressure","sampled_k",
                    "constitutive","local_advection","sampled_advection"};
                Vec total=Vec::Zero(size),linear=Vec::Zero(size);
                for(std::size_t i=0;i<labels.size();++i)
                {
                    Vec rows=Vec::Zero(nv+m_flow.PressureSize());
                    rows.head(nv)=(finePlus[i]-fineMinus[i])/(2.*half)-(plus[i]-minus[i])/(2.*h);
                    Vec term=Vec::Zero(size);
                    term.head(nf)=m_flow.Restrict(rows);
                    total+=term;
                    if(i<3)linear+=term;
                    audit.arrays.emplace_back(std::string("raw_spacing_momentum_")+labels[i],std::move(term));
                }
                const Vec measured=raw("raw_spacing_change_1");
                Vec closure=Vec::Zero(size),remain=Vec::Zero(size);
                closure.head(npOffset)=(total-measured).head(npOffset);
                remain.head(npOffset)=(measured-linear).head(npOffset);
                audit.arrays.emplace_back("raw_momentum_term_closure",std::move(closure));
                audit.arrays.emplace_back("raw_momentum_spacing_without_linear_terms",std::move(remain));
                Eigen::Index row;
                measured.head(npOffset).cwiseProduct(sr.head(npOffset)).cwiseAbs().maxCoeff(&row);
                const auto &free=m_flow.FreeVelocityIds();
                const int component=row<static_cast<Eigen::Index>(free[0].size())?0:1;
                const int id=free[component][row-(component?free[0].size():0)];
                std::ostringstream owners;
                owners<<std::setprecision(17);
                for(const auto &cell:m_flow.Elements())
                    if(std::find(cell.u.begin(),cell.u.end(),id)!=cell.u.end())
                        owners<<cell.geometry<<" ["<<cell.xl<<','<<cell.xr<<','<<cell.rl<<','<<cell.rr<<"] ";
                audit.measurements.emplace_back("peak_momentum_free_row",row);
                audit.measurements.emplace_back("peak_momentum_component",component);
                audit.notes.emplace_back("peak_momentum_owners",owners.str());
                audit.measurements.emplace_back("momentum_attribution_seconds",Clock()-start);
            }
            if(observed.stage=="frozen-size-finite-trial"&&
                std::any_of(observed.arrays.begin(),observed.arrays.end(),[](const auto &entry) {
                    return entry.first=="raw_actual_residual";
                }))
            {
                localized=observed;
                const auto field=unpack(raw("trial_free"));
                const auto k=InterfaceBounds(field);
                const double kin=m_reference->Parameters().KIn;
                const int depth=diagnostic->scalarSubdivisionDepth;
                localized->notes.emplace_back("observed_trial_coefficient_sha256",StateHash(field));
                localized->measurements.emplace_back("whole_k_lower",m_k.PolynomialBounds(field.k,depth)[0]);
                localized->measurements.emplace_back("whole_epsilon_lower",m_e.PolynomialBounds(field.epsilon,depth)[0]);
                localized->measurements.emplace_back("whole_K_over_Kin_lower",k[0]/kin);
                localized->measurements.emplace_back("whole_K_over_Kin_upper",k[1]/kin);
                try
                {
                    const auto scientific=Evaluate(field,0,PrincipalPolicy::Observe,timing,diagnostic);
                    localized->measurements.emplace_back("sampled_principal_over_mu",scientific.minPrincipal);
                    localized->measurements.emplace_back("sampled_covariance_over_k",scientific.minCovariance);
                    localized->measurements.emplace_back("projected_trace_error",scientific.traceError);
                    localized->notes.emplace_back("principal_observation","OBSERVED");
                }
                catch(const NativeInadmissibleState &error)
                {
                    localized->notes.emplace_back("principal_observation","UNRESOLVED");
                    localized->notes.emplace_back("principal_observation_reason",error.what());
                }
            }
            const auto &audit=localized?*localized:observed;
            const double writing=Clock();
            const std::uint16_t endian=1;
            Require(sizeof(double)==8&&std::numeric_limits<double>::is_iec559&&
                    *reinterpret_cast<const unsigned char*>(&endian)==1,"action audit requires little-endian binary64");
            const auto prefix=diagnostic->actionAuditPrefix+audit.FileSuffix(++sequence);
            const auto binary=prefix+".bin";
            PublishDiagnostic(binary,[&](int fd) {
                for(const auto &[name,values]:audit.arrays)
                    D::WriteAll(fd,values.data(),values.size()*sizeof(double));
            });
            std::ostringstream out;
            out<<std::setprecision(17)<<"{\"schema\":\"WJH27SteadyActionAuditV1\",\"role\":\"NONACCEPTED_ACTION_VECTOR\""
                <<",\"physical_history\":false,\"stage\":"<<J::Quote(audit.stage)
                <<",\"source_identity\":"<<J::Quote(m_options.sourceIdentity)
                <<",\"reference_identity\":"<<J::Quote(m_reference->Identity())
                <<",\"signature\":"<<J::Quote(m_signature)
                <<",\"equation_identity\":"<<J::Quote(D::Hash(m_signature+"|steady"))
                <<",\"initial_packet_sha256\":"<<J::Quote(m_initialPacketSha)
                <<",\"base_coefficient_sha256\":"<<J::Quote(StateHash(unpack(audit.arrays.front().second.col(0))))
                <<",\"session_sha256\":[";
            for(std::size_t i=0;i<diagnostic->checkpointSessionDigests.size();++i)
            {
                if(i)out<<',';
                out<<J::Quote(diagnostic->checkpointSessionDigests[i]);
            }
            out<<"],\"probe_input_sha256\":"<<J::Quote(diagnostic->actionProbeSha)
                <<",\"free_order\":\"free u, free v, all Pi, optional gauge, free k, free outer epsilon\""
                <<",\"free_id_sha256\":"<<J::Quote(orderingSha)
                <<",\"free_counts\":["<<m_flow.FreeVelocityIds()[0].size()<<','<<m_flow.FreeVelocityIds()[1].size()
                <<','<<m_flow.PressureSize()<<','<<(m_flow.Gauge()?1:0)<<','<<nk<<','<<ne<<']'
                <<",\"quadrature\":["<<m_options.quadratureX<<','<<m_options.quadratureR<<']'
                <<",\"current_fd\":"<<options.fdRelativeStep<<",\"correction_mass\":"<<correctionMass
                <<",\"fixed_k_tangent\":"<<J::Boolean(diagnostic->fixedKTangent)
                <<",\"inverse_only_mass\":"<<preconditionerMass
                <<",\"newton\":"<<audit.newton<<",\"cycle\":"<<audit.cycle<<",\"iteration\":"<<audit.iteration
                <<",\"proposal\":"<<audit.proposal<<",\"capture_sequence\":"<<sequence
                <<",\"factor_serial\":"<<audit.factorSerial<<",\"factor\":";
            if(audit.factorSerial)
            {
                const auto &work=m_sparseSteadyWork.at(audit.factorSerial-1);
                out<<"{\"assembly_base_sha256\":"<<J::Quote(work.baseSha256)
                   <<",\"total_mass\":"<<work.inverseMass<<",\"fd_relative_step\":"<<work.fdRelativeStep
                   <<",\"colors\":"<<work.colors<<",\"matrix_nonzeros\":"<<work.nonzeros
                   <<",\"factor_nonzeros\":"<<work.factorNonzeros<<'}';
            }
            else out<<"null";
            std::vector<std::pair<std::string,std::string>> measures;
            for(const auto &[name,value]:audit.measurements)measures.emplace_back(name,J::Number(value));
            std::vector<std::pair<std::string,std::string>> notes;
            for(const auto &[name,value]:audit.notes)notes.emplace_back(name,J::Quote(value));
            out<<",\"measurements\":"<<J::Object(measures)<<",\"notes\":"<<J::Object(notes)
               <<",\"raw_payload\":"<<D::Pin(binary)
               <<",\"encoding\":\"little-endian binary64, column major\",\"arrays\":[";
            std::size_t offset=0;
            for(const auto &[name,values]:audit.arrays)
            {
                if(offset)out<<',';
                const auto bytes=values.size()*sizeof(double);
                out<<"{\"name\":"<<J::Quote(name)<<",\"offset_bytes\":"<<offset
                   <<",\"rows\":"<<values.rows()<<",\"columns\":"<<values.cols()
                   <<",\"sha256\":"<<J::Quote(D::Hash(values.data(),bytes));
                if(name.rfind("raw_",0)==0&&values.cols()==1)
                {
                    const Vec raw=values.col(0);
                    out<<",\"scaled_l2\":"<<sr.cwiseProduct(raw).stableNorm()
                       <<",\"physical_blocks\":"<<J::Numbers(blockNorms(raw));
                }
                out<<'}';
                offset+=bytes;
            }
            out<<"],\"stencils\":[";
            for(std::size_t i=0;i<audit.stencils.size();++i)
            {
                const auto &fd=audit.stencils[i];
                if(i)out<<',';
                out<<"{\"purpose\":"<<J::Quote(fd.purpose)<<",\"method\":"<<J::Quote(fd.method)
                   <<",\"reduction\":"<<fd.reduction<<",\"relative\":"<<fd.relativeStep
                   <<",\"parameter\":"<<fd.parameterStep<<",\"coordinate\":"<<fd.scaledCoordinateStep
                   <<",\"plus_displacement\":"<<fd.plusDisplacement<<",\"minus_displacement\":"<<fd.minusDisplacement
                   <<",\"half_displacement\":"<<(fd.halfDisplacement>0.?J::Number(fd.halfDisplacement):"null")
                   <<",\"scaled_direction_infinity\":"<<fd.scaledDirectionInfinityNorm
                   <<",\"plus\":"<<J::Boolean(fd.plusAdmissible)<<",\"minus\":"<<J::Boolean(fd.minusAdmissible)
                   <<",\"half\":"<<J::Boolean(fd.halfAdmissible)<<",\"plus_reason\":"<<J::Quote(fd.plusReason)
                   <<",\"minus_reason\":"<<J::Quote(fd.minusReason)<<",\"half_reason\":"<<J::Quote(fd.halfReason)<<'}';
            }
            out<<"]}\n";
            const auto contents=out.str();
            PublishDiagnostic(prefix+".json",[&](int fd){D::WriteAll(fd,contents.data(),contents.size());});
            timing.io+=Clock()-writing;
            std::cout<<"WJ-H27 action audit "<<audit.stage<<" published in "<<Clock()-start<<" s"<<std::endl;
        };
        if(sparse)
        {
            callbacks.ColoredAction=[&](const Vec &d){return coupledInverse->OriginalAction(d);};
            callbacks.FactoredAction=[&](const Vec &d)->Vec {
                return coupledInverse->OriginalAction(d)+shiftAt(auditedFactorState,d,auditedFactorMass,0.);
            };
        }
    }
    auto refreshState=[&](CoupledState state,const Data &d,const std::string &reason,int iteration,bool accepted) {
        auto load=forcing;
        auto adv=Advection(state,d.kRules);
        load.flow-=adv.flow;load.k-=adv.k;load.epsilon-=adv.epsilon;
        SteadyPacket packet;
        if(accepted)packet.meritController=m_currentMerit;
        auto &original=packet.original;
        original.residual=Residual(state,d,load,0,0,timing);
        original.minimumK=d.minK;original.minimumEpsilon=d.minE;
        original.minimumCovarianceOverK=d.minCovariance;original.minimumPrincipalOverMu=d.minPrincipal;
        packet.modifiedMargin=d.minModifiedPrincipal;
        packet.wholeCellLower={m_k.PolynomialBounds(state.k)[0],m_e.PolynomialBounds(state.epsilon)[0]};
        auto solved=original.residual;
        if(modified)
        {
            const Vec extra=regularization*state.flow.head(nv);
            packet.extraWork=state.flow.head(nv).dot(extra);
            load.flow.head(nv)-=extra;
            solved=Residual(state,d,load,0,0,timing);
        }
        packet.candidate={DiagnosticRole::SteadySearch,std::move(state),0,0.,reason,solved,iteration,accepted};
        return std::make_pair(std::move(packet),std::move(adv));
    };
    Vec lastObserved;
    if(diagnostic)
    {
        callbacks.ResidualBlocks=blockNorms;
        if(!options.evaluateInitialOnly)callbacks.ObserveAccepted=[&](const Vec &x,double norm) {
            if(x.size()==lastObserved.size()&&(x.array()==lastObserved.array()).all())return;
            const double start=Clock();
            auto s=unpack(x);
            const auto observed=Evaluate(s,0,PrincipalPolicy::Observe,timing,diagnostic);
            const auto k=InterfaceBounds(s);
            const auto kin=m_reference->Parameters().KIn;
            m_steadyBaseObservations.push_back({norm,observed.minPrincipal,observed.minCovariance,
                                               k[0]/kin,k[1]/kin,Clock()-start});
            lastObserved=x;
            std::cout<<std::setprecision(12)<<"WJ-H27 steady base "<<m_steadyBaseObservations.size()-1
                     <<" residual "<<norm<<" original principal "<<observed.minPrincipal
                     <<" K/Kin ["<<k[0]/kin<<','<<k[1]/kin<<']'<<std::endl;
            const auto sequence=m_steadyBaseObservations.size()-1;
            if(sequence && !diagnostic->acceptedCheckpointPrefix.empty())
            {
                const double refreshStart=Clock();
                auto packet=refreshState(std::move(s),observed,"accepted correction during diagnostic solve",
                                         static_cast<int>(sequence)-1,true).first;
                packet.checkpointSequence=sequence;
                packet.checkpointRefreshSeconds=Clock()-refreshStart;
                std::ostringstream path;
                path<<diagnostic->acceptedCheckpointPrefix<<std::setfill('0')<<std::setw(6)<<sequence<<".json";
                const double writeStart=Clock();
                WriteDiagnostic(path.str(),packet.candidate,&packet);
                const double writeSeconds=Clock()-writeStart;
                timing.io+=writeSeconds;
                std::cout<<"WJ-H27 accepted checkpoint "<<sequence<<" "<<StateHash(packet.candidate.state)
                         <<" refresh "<<packet.checkpointRefreshSeconds<<" s, publication "<<writeSeconds<<" s"
                         <<std::endl;
            }
        };
    }
    // Analyze a bounded sample of the actual line-search trials. The fifth
    // vector measures the effect of varying the dependent trace by difference
    // from a frozen-trace evaluation. It overlaps the four physical terms.
    using EpsilonActions=std::array<Vec,5>;
    EpsilonActions epsilonBase,epsilonSlope;
    Vec baseTrace,traceDirection;
    std::optional<double> halfSpacingDifference;
    std::string epsilonUnavailable;
    if(diagnostic&&options.measureDirection)
    {
        auto terms=[&](const CoupledState &s,const Data &d) {
            const auto &rules=m_rules.Static();
            std::array<Vec,4> value;
            value[0]=m_e.Restrict(m_e.Advection(rules,s.flow,s.epsilon));
            value[1]=m_e.Restrict(m_e.DiffusionAction(rules,d.De,s.epsilon));
            value[2]=-m_e.Restrict(EVolumeLoad(s,d));
            value[3]=m_e.Restrict(m_e.DiffusionAction(rules,d.De,s.epsilon,d.reactionE))-value[1];
            return value;
        };
        auto samples=[&,terms](const Vec &x,bool isolateTrace) {
            ++m_steadyEpsilonEvaluations;
            auto s=unpack(x);
            auto value=terms(s,*currentData(x,s));
            EpsilonActions rows{value[0],value[1],value[2],value[3],Vec::Zero(ne)};
            if(isolateTrace)
            {
                s.epsilon+=baseTrace-m_e.Lift(s.epsilon);
                ++m_steadyEpsilonEvaluations;
                const auto frozen=terms(s,Evaluate(s,0,PrincipalPolicy::Skip,timing,diagnostic));
                for(int j=0;j<4;++j)rows[4]+=value[j]-frozen[j];
            }
            return rows;
        };
        callbacks.InspectDirection=[&,samples](int iteration,const Vec &x,const Vec &dx) {
            if(iteration>=3)return;
            const double started=Clock();
            epsilonUnavailable.clear();halfSpacingDifference.reset();
            baseTrace=m_e.Lift(pcState.epsilon);
            traceDirection=DependentTraceDirection(pcState,
                m_k.Expand(dx.segment(nf,nk),Vec::Zero(m_k.Size())));
            const double h=options.fdRelativeStep*std::max(1.,x.cwiseQuotient(sx).cwiseAbs().maxCoeff())/
                dx.cwiseQuotient(sx).cwiseAbs().maxCoeff();
            try
            {
                epsilonBase=samples(x,false);
                const auto plus=samples(x+h*dx,true),minus=samples(x-h*dx,true);
                const auto halfPlus=samples(x+.5*h*dx,false),halfMinus=samples(x-.5*h*dx,false);
                Vec difference=Vec::Zero(ne);
                for(int j=0;j<5;++j)
                {
                    epsilonSlope[j]=(plus[j]-minus[j])/(2.*h);
                    if(j<4)difference+=epsilonSlope[j]-(halfPlus[j]-halfMinus[j])/h;
                }
                halfSpacingDifference=m_e.FreeResidualNorm(difference);
            }
            catch(const NativeInadmissibleState &error){epsilonUnavailable=error.what();}
            m_steadyEpsilonSeconds+=Clock()-started;
        };
        callbacks.InspectTrial=[&,samples](int iteration,int trial,double step,const Vec &trialX,const Vec &rawRemainder) {
            if(iteration>=3||trial>=6)return;
            const double started=Clock();
            auto &record=m_steadyEpsilonRemainders.emplace_back();
            record.iteration=iteration;record.trial=trial;record.step=step;
            record.halfSpacingDifference=halfSpacingDifference;
            record.unavailable=epsilonUnavailable;
            if(record.unavailable.empty())
            {
                try
                {
                    const auto actual=samples(trialX,true);
                    Vec sum=Vec::Zero(ne);
                    for(int j=0;j<4;++j)
                    {
                        const Vec term=actual[j]-epsilonBase[j]-step*epsilonSlope[j];
                        sum+=term;
                        record.terms[j]=m_e.FreeResidualNorm(term);
                    }
                    const auto tail=rawRemainder.tail(ne);
                    record.closure=m_e.FreeResidualNorm(sum-tail);
                    record.liftResidual=m_e.FreeResidualNorm(actual[4]-step*epsilonSlope[4]);
                    const auto state=unpack(trialX);
                    record.liftCoefficient=(m_e.Lift(state.epsilon)-baseTrace-step*traceDirection).lpNorm<Eigen::Infinity>();
                    Eigen::Index row;
                    tail.cwiseProduct(m_e.ResidualRowScale()).cwiseAbs().maxCoeff(&row);
                    record.epsilonId=m_e.FreeIds()[row];
                    for(const auto &cell:m_flow.Elements())
                        if(m_e.Owns(cell.geometry))
                        {
                            const auto &ids=m_e.CellIds(cell.geometry);
                            if(std::find(ids.begin(),ids.end(),record.epsilonId)!=ids.end())
                                record.owners.push_back(cell.geometry);
                        }
                }
                catch(const NativeInadmissibleState &error){record.unavailable=error.what();}
            }
            m_steadyEpsilonSeconds+=Clock()-started;
        };
    }
    const auto checkedResidual=callbacks.Residual;
    auto auditInverse=[&](const Vec &x,const Vec &raw){
        const double start=Clock();
        std::set<int> matching;
        for(const auto &[id,trace]:m_hKTrace)
            matching.insert(trace.ids.begin(),trace.ids.end());
        std::vector<std::pair<std::string,Vec>> directions;
        for(const char *name:{"velocity","pressure/gauge","bulk k","matching k","outer epsilon","mixed"})
            directions.emplace_back(name,Vec::Zero(size));
        for(int i=0;i<size;++i)
        {
            const double value=sx[i]*(.6*std::sin(.47*(i+1))+.4*std::cos(.21*(i+1)));
            const int block=i<npOffset?0:i<nf?1:i<nf+nk?
                (matching.count(m_k.FreeIds()[i-nf])?3:2):4;
            directions[block].second[i]=value;
            directions.back().second[i]=value;
        }
        const double scale=std::max(1.,x.cwiseQuotient(sx).cwiseAbs().maxCoeff());
        for(const auto &[label,direction]:directions)
        {
            const Vec assembled=coupledInverse->OriginalAction(direction);
            const Vec shift=shifted?applyShift(direction,correctionMass,correctionDiffusion):Vec::Zero(size).eval();
            for(double multiplier:{8.,2.,.5,.125})
            {
                deadline();
                SteadyActionAudit record;
                record.direction=label;
                record.step=multiplier*options.fdRelativeStep*scale;
                const double h=record.step;
                try
                {
                    const Vec plus=x+h*direction,minus=x-h*direction;
                    ++m_steadyAuditResiduals;
                    const Vec fp=checkedResidual(plus);
                    ++m_steadyAuditResiduals;
                    const Vec fm=checkedResidual(minus);
                    const Vec centered=(fp-fm)/(2.*h);
                    const Vec error=assembled-centered;
                    record.differenceBlocks=blockNorms(error);
                    record.centeredRelativeError=error.cwiseProduct(sr).norm()/
                        std::max(1e-20,centered.cwiseProduct(sr).norm());
                    record.taylorRemainder=(fp-raw-h*assembled).cwiseProduct(sr).norm();
                    if(shifted)
                    {
                        const auto delta=Combination(unpack(plus),.5/h,unpack(minus),-.5/h);
                        const Vec reference=differentialShift(delta,correctionMass,correctionDiffusion);
                        record.shiftRelativeError=(shift-reference).cwiseProduct(sr).norm()/
                            std::max(1e-20,reference.cwiseProduct(sr).norm());
                    }
                }
                catch(const SteadyInadmissibleTrial &error)
                {
                    record.rejection=error.what();
                }
                m_steadyActionAudits.push_back(std::move(record));
            }
        }
        m_steadyAuditSeconds+=Clock()-start;
    };
    double observedPreparationSeconds=0.;
    callbacks.PreparePreconditioner=[&](const Vec&x,const Vec&raw){
        const double preparationStart=Clock();
        if(diagnostic&&!diagnostic->acceptedCheckpointPrefix.empty())
        {
            const double remaining=std::min(began+options.maximumSeconds,
                diagnostic->checkpointProcessDeadline)-preparationStart;
            const double estimate=std::max(diagnostic->checkpointInverseSeconds,observedPreparationSeconds);
            if(remaining<=estimate+diagnostic->checkpointReserveSeconds)
            {
                std::ostringstream why;
                why<<"DEADLINE: remaining "<<remaining<<" s before inverse; estimated "
                   <<estimate<<" s plus "<<diagnostic->checkpointReserveSeconds<<" s for output";
                throw SteadyCorrectionFailure(why.str());
            }
            std::cout<<"WJ-H27 inverse preparation start; estimated "<<estimate
                     <<" s, remaining "<<remaining<<" s"<<std::endl;
        }
        const double inverseMass=preconditionerMass+correctionMass;
        if(sparse)
        {
            SteadySparseInverse::Action addShift;
            if(inverseMass>0.||correctionDiffusion>0.)
                addShift=[&](const Vec &dx){return applyShift(dx,inverseMass,correctionDiffusion);};
            auto &work=m_sparseSteadyWork.emplace_back();
            work.baseNorm=norms(raw).Maximum();
            work.baseSha256=StateHash(pcState);
            work.inverseMass=inverseMass;
            work.stiffness=correctionDiffusion;
            auto factorOptions=options;
            if(diagnostic->factorFDStep>0.)factorOptions.fdRelativeStep=diagnostic->factorFDStep;
            coupledInverse->Prepare(x,checkedResidual,addShift,factorOptions,deadline,work);
            if(!diagnostic->actionAuditPrefix.empty())
            {
                auditedFactorState=pcState;
                auditedFactorMass=inverseMass;
                const double auditStart=Clock();
                const Vec rhs=raw/raw.cwiseProduct(sr).stableNorm();
                const Vec solved=coupledInverse->Solve(rhs);
                const Vec unshifted=coupledInverse->OriginalAction(solved)-rhs;
                const Vec shiftedError=unshifted+applyShift(solved,inverseMass,correctionDiffusion);
                work.backsolveAudit=std::array<double,2>{sr.cwiseProduct(shiftedError).stableNorm(),
                                                       sr.cwiseProduct(unshifted).stableNorm()};
                work.backsolveAuditSeconds=Clock()-auditStart;
                std::cout<<"WJ-H27 factor backsolve scaled check shifted "<<(*work.backsolveAudit)[0]
                         <<", unshifted "<<(*work.backsolveAudit)[1]<<std::endl;
            }
            std::cout<<"WJ-H27 coupled inverse "<<work.colors<<" colors, "<<work.nonzeros
                     <<" entries, assembly "<<work.assemblySeconds<<" s, factor "
                     <<work.factorSeconds<<" s, correction mass "<<correctionMass
                     <<", inverse-only mass "<<preconditionerMass<<", total "<<inverseMass
                     <<", factor FD "<<work.fdRelativeStep<<", current FD "<<options.fdRelativeStep<<std::endl;
            if(diagnostic->auditCorrection&&m_sparseSteadyWork.size()==1)auditInverse(x,raw);
            observedPreparationSeconds=std::max(observedPreparationSeconds,Clock()-preparationStart);
            return;
        }
        auto owned=currentData(x,pcState);
        const auto &d=*owned;
        Sparse flowPC=(m_options.steadyVelocityPC=="FullGradientV1"?VelocityTangent(d,timing):Strain(d))+MomentumTransportJacobian(pcState);
        if(modified)flowPC+=regularization;
        if(correctionDiffusion>0.)flowPC+=correctionDiffusion*positiveDiffusion;
        Sparse kPC=m_k.Diffusion(d.kRules,d.Dk)+m_k.TransportMatrix(d.kRules,pcState.flow);
        Sparse ePC=m_e.Diffusion(m_rules.Static(),d.De,d.reactionE)+m_e.TransportMatrix(m_rules.Static(),pcState.flow);
        // preconditionerMass belongs solely to the inverse approximation.
        epsilonPCFull=ePC+inverseMass*m_e.Mass();
        m_flow.PrepareLinearPreconditioner(m_flow.System(inverseMass,flowPC,1));
        m_k.PrepareLinearPreconditioner(m_k.System(inverseMass,kPC));
        m_e.PrepareLinearPreconditioner(m_e.System(inverseMass,ePC));
        observedPreparationSeconds=std::max(observedPreparationSeconds,Clock()-preparationStart);
    };
    callbacks.Precondition=[&](const Vec&rhs){
        if(sparse)return coupledInverse->Solve(rhs);
        Vec result(size);result.segment(nf,nk)=m_k.ApplyLinearPreconditioner(rhs.segment(nf,nk));
        Vec dk=m_k.Expand(result.segment(nf,nk),Vec::Zero(m_k.Size()));
        Vec cross=Vec::Zero(2*m_flow.VelocitySize()+m_flow.PressureSize());cross.head(2*m_flow.VelocitySize())=m_flow.KnownScalarKLoad(dk);
        result.head(nf)=m_flow.ApplyLinearPreconditioner(rhs.head(nf)+m_flow.Restrict(cross));
        const auto trace=DependentTraceDirection(pcState,dk);
        result.tail(ne)=m_e.ApplyLinearPreconditioner(rhs.tail(ne)-m_e.Restrict(epsilonPCFull*trace));return result;
    };
    if(m_options.checkCurrentActions){
        // Uncommitted audit only: explicitly exercise a weakly divergent
        // free-velocity state and the otherwise often-zero gauge multiplier.
        Vec probe=pack(seed);probe[0]+=1e-7;
        if(m_flow.Gauge())probe[nf-1]=.003;
        (void)equations(probe);
    }
    SteadyKrylov solver(std::move(callbacks),sx,sr,options);
    const Vec initial=pack(seed);
    if(options.frozenSizeColumns)
        Require(diagnostic&&StateHash(unpack(initial))==m_initialSourceStateSha,
                "frozen correction basis requires an exact signed base without projection");
    if(diagnostic&&!diagnostic->actionProbePath.empty())
    {
        Require(StateHash(unpack(initial))==m_initialSourceStateSha,
                "fixed direction requires an exact signed base without projection");
        std::ifstream file(diagnostic->actionProbePath,std::ios::binary);
        Mat directions(size,diagnostic->actionProbeCount);
        const auto bytes=directions.size()*sizeof(double);
        file.read(reinterpret_cast<char*>(directions.data()),bytes);
        Require(file.gcount()==static_cast<std::streamsize>(bytes)&&file.peek()==std::char_traits<char>::eof()&&
                D::Hash(directions.data(),bytes)==diagnostic->actionProbeSha,"fixed direction size or digest differs");
        if(options.frozenFiniteTrialLimit)
        {
            D::Sha256 ordering;
            for(const auto *ids:{&m_flow.FreeVelocityIds()[0],&m_flow.FreeVelocityIds()[1],&m_k.FreeIds(),&m_e.FreeIds()})
                ordering.Update(ids->data(),ids->size()*sizeof(int));
            const auto binding=D::Hash("WJH27FrozenFiniteBindingV1|"+m_initialSourceStateSha+"|"+
                ordering.Finish()+"|"+D::Hash(sx.data(),sx.size()*sizeof(double))+"|"+
                D::Hash(sr.data(),sr.size()*sizeof(double))+"|"+diagnostic->actionProbeSha);
            Require(binding==diagnostic->actionProbeBinding,"finite-trial input base, free ordering or scales differ");
        }
        const Vec direction=directions.col(0);
        if(diagnostic->actionProbeCount==1&&!options.frozenFiniteTrialLimit)solver.AuditDirection(initial,direction);
        else solver.AuditDirections(initial,directions);
        if(diagnostic->auditColorColumns)
        {
            auto pattern=MakeSteadySparseInverse(sx,sr);
            const double foot=m_geometry.normalMap->lengthX;
            std::array<std::set<int>,3> selected;
            auto addFree=[](std::set<int> &out,const auto &ids,const std::vector<int> &free,int offset) {
                for(int id:ids)
                {
                    const auto pos=std::lower_bound(free.begin(),free.end(),id);
                    if(pos!=free.end()&&*pos==id)out.insert(offset+int(pos-free.begin()));
                }
            };
            for(const auto &cell:m_flow.Elements())
                if(cell.xl<foot+.05&&cell.xr>foot-.05&&cell.rl<m_geometry.H+.04)
                {
                    addFree(selected[0],m_hKTrace.at(HOwner(cell).geometry).ids,m_k.FreeIds(),nf);
                    if(m_e.Owns(cell.geometry))
                        addFree(selected[1],m_e.CellIds(cell.geometry),m_e.FreeIds(),nf+nk);
                    if(cell.rl<m_geometry.H&&cell.rr>m_geometry.h)
                    {
                        const auto &free=m_flow.FreeVelocityIds();
                        addFree(selected[2],cell.u,free[0],0);
                        addFree(selected[2],cell.u,free[1],int(free[0].size()));
                    }
                }
            const std::array<std::string,3> labels{"color-matching-k","color-outer-epsilon","color-foot-velocity"};
            for(std::size_t i=0;i<selected.size();++i)
            {
                Require(!selected[i].empty(),"selected foot column inventory is empty");
                const int col=*std::max_element(selected[i].begin(),selected[i].end(),[&](int a,int b) {
                    return std::abs(direction[a]/sx[a])<std::abs(direction[b]/sx[b]);
                });
                solver.AuditColumn(initial,direction,pattern->ColumnProbe(col),labels[i]);
            }
        }
    }
    auto result=solver.Solve(initial);
    auto refresh=[&](const Vec &x,int iteration,bool accepted) {
        auto state=unpack(x);
        const auto d=Evaluate(state,0,diagnostic?PrincipalPolicy::Observe:PrincipalPolicy::Enforce,timing,diagnostic);
        return refreshState(std::move(state),d,result.reason,iteration,accepted);
    };
    if(result.converged||(diagnostic&&result.hasEvaluatedState)){
        auto [packet,adv]=refresh(result.lastAcceptedX,result.lastAcceptedNewtonIteration,true);
        final.residual=*packet.candidate.residual;
        final.minimumK=packet.original.minimumK;final.minimumEpsilon=packet.original.minimumEpsilon;
        final.minimumCovarianceOverK=packet.original.minimumCovarianceOverK;
        final.minimumPrincipalOverMu=packet.original.minimumPrincipalOverMu;
        if(result.converged&&final.residual.Maximum()>options.acceptanceNorm){
            result.converged=false;result.status="NOT_CONVERGED";
            result.reason="fresh steady residual failed after eliminating the gauge multiplier";
        }
        if(result.converged&&physicalHistory){
            CoupledHistory fresh;fresh.dt=m_options.dt;fresh.current=packet.candidate.state;fresh.previous=fresh.current;fresh.advection=adv;fresh.previousAdvection=adv;CheckHistory(fresh);*physicalHistory=std::move(fresh);
        }
        if(diagnostic)
        {
            packet.candidate.reason=result.reason;
            m_steadyDiagnosticCandidate=std::move(packet);
            if(!result.bestFromAcceptedIterate)
                m_steadyUncommitted=refresh(result.x,result.bestNewtonIteration,false).first;
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
    auto payload=stream.str();payload+=D::Hash(payload);auto tmp=path.string()+".tmp";int fd=::open(tmp.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);D::Check(fd>=0,"coupled checkpoint create");D::WriteAll(fd,payload.data(),payload.size());D::Check(::fsync(fd)==0,"coupled checkpoint sync");::close(fd);std::filesystem::rename(tmp,path);D::SyncDirectory(std::filesystem::absolute(path));
}
CoupledHistory Coupling::ReadHistory(const std::filesystem::path &path,
                                    CoupledRecoveryTiming *recoveryTiming) const
{
    const double readBegan = Clock();
    std::ifstream file(path,std::ios::binary);std::string bytes((std::istreambuf_iterator<char>(file)),{});Require(bytes.size()>64&&bytes.size()<128*1024*1024,"invalid coupled checkpoint size");auto payload=bytes.substr(0,bytes.size()-64);Require(D::Hash(payload)==bytes.substr(bytes.size()-64),"coupled checkpoint checksum mismatch");std::istringstream stream(payload,std::ios::binary);std::string magic,identity;std::getline(stream,magic);std::getline(stream,identity);Require(magic=="WJ_H27_NATIVE_COUPLED_V2"&&identity==m_signature,"coupled checkpoint scientific identity differs");
    CoupledHistory h;stream.read(reinterpret_cast<char*>(&h.step),sizeof(h.step));stream.read(reinterpret_cast<char*>(&h.time),sizeof(double));stream.read(reinterpret_cast<char*>(&h.dt),sizeof(double));
    for(auto*s:{&h.current,&h.previous,&h.advection,&h.previousAdvection}){int which=0;for(auto*v:{&s->flow,&s->k,&s->epsilon}){std::uint64_t n=0;stream.read(reinterpret_cast<char*>(&n),sizeof(n));const int expected=which==0?2*m_flow.VelocitySize()+m_flow.PressureSize():which==1?m_k.Size():m_e.Size();++which;Require(n==std::uint64_t(expected),"coupled checkpoint vector inventory differs");v->resize(n);stream.read(reinterpret_cast<char*>(v->data()),n*sizeof(double));}}
    Require(stream.good()&&stream.peek()==std::char_traits<char>::eof(),"coupled checkpoint has truncated/trailing bytes");CheckHistory(h);
    const double rebuildBegan = Clock();
    // Reconstruct scientific derived data and compare the explicit histories.
    // Normal factory cold restarts use a new Coupling instance. External
    // same-instance reloads call ResetPreconditioners for cold inverse setup.
    CoupledTiming timing;auto data=Evaluate(h.current,h.time,true,timing);Require(data.traceError<1e-10,"restored dependent trace differs");auto adv=Advection(h.current,data.kRules);
    for(const auto pair:{std::make_pair(&adv.flow,&h.advection.flow),std::make_pair(&adv.k,&h.advection.k),std::make_pair(&adv.epsilon,&h.advection.epsilon)})Require((*pair.first-*pair.second).norm()<=1e-12*std::max(pair.first->norm(),1e-30),"restored explicit history differs from current state");
    auto previous=Evaluate(h.previous,h.step>0?(h.step-1)*h.dt:0,true,timing);auto oldAdv=Advection(h.previous,previous.kRules);
    for(const auto pair:{std::make_pair(&oldAdv.flow,&h.previousAdvection.flow),std::make_pair(&oldAdv.k,&h.previousAdvection.k),std::make_pair(&oldAdv.epsilon,&h.previousAdvection.epsilon)})Require((*pair.first-*pair.second).norm()<=1e-12*std::max(pair.first->norm(),1e-30),"restored previous explicit history differs");
    if (recoveryTiming)
    {
        recoveryTiming->readSeconds = rebuildBegan - readBegan;
        recoveryTiming->cacheRebuildSeconds = Clock() - rebuildBegan;
    }
    return h;
}

std::array<double,5>Coupling::MMSErrors(const CoupledState&s,double time)const
{
    Require(bool(m_mms),"MMS errors requested for application");std::array<double,5>error{};double outerVolume=0;const int n=m_flow.VelocitySize();
    double pressureMean=0.;
    if (m_geometry.normalMap)
    {
        double volume=0.;
        for (const auto &rule : m_rules.Static())
            for (int q=0;q<rule->w.size();++q)
            {
                volume+=rule->w[q];
                pressureMean+=rule->w[q]*m_mms->Fields(rule->x[q],rule->z[q],time).Pi.value;
            }
        pressureMean/=volume;
    }
    for(const auto&rp:m_rules.Static()){
        const auto&r=*rp;Vec u=r.B*Local(s.flow,r.u),v=r.B*Local(s.flow,r.u,n),p=r.P*Local(s.flow,r.p,2*n),k=m_k.Evaluate(r,s.k),e;
        bool outer=m_e.Owns(r.geometry);if(outer)e=m_e.Evaluate(r,s.epsilon);
        for(int q=0;q<r.w.size();++q){auto exact=m_mms->Fields(r.x[q],r.z[q],time);error[0]+=r.w[q]*(std::pow(u[q]-exact.u.value,2)+std::pow(v[q]-exact.w.value,2));error[1]+=r.w[q]*std::pow(p[q]-(exact.Pi.value-pressureMean),2);error[2]+=r.w[q]*std::pow((k[q]-exact.k.value)/m_reference->Parameters().KIn,2);
            if(outer){const double ex=m_mms->OuterEpsilon(r.x[q],r.z[q],time).value;error[3]+=r.w[q]*std::pow((e[q]-ex)/ex,2);outerVolume+=r.w[q];}
            Vec ul=Local(s.flow,r.u),vl=Local(s.flow,r.u,n);Tensor G{};G[0][0]=r.Dx.row(q).dot(ul);G[0][2]=r.Dz.row(q).dot(ul);G[2][0]=r.Dx.row(q).dot(vl);G[2][2]=r.Dz.row(q).dot(vl);
            for(int i:{0,2})for(int j:{0,2})error[4]+=r.w[q]*std::pow(G[i][j]-exact.gradient[i][j],2);
        }
    }
    for(int i=0;i<5;++i)error[i]=std::sqrt(error[i]/(i==3?outerVolume:m_volume));return error;
}
void Coupling::WriteOutletActions(std::ostream&out,const CoupledHistory&h)const
{
    Require(h.step==0&&h.time==0,"outlet action export requires a zero-time state");
    CoupledTiming timing;
    const auto data=Evaluate(h.current,h.time,false,timing);
    const auto forcing=Forcing(h.time,timing);
    const auto advection=Advection(h.current,data.kRules);
    const auto current=CurrentActions(h.current,data);
    const int n=m_flow.VelocitySize();
    out<<",\"outlet_actions\":{\"normal_applied\":";
    VectorJSON(out,m_outletLoad.head(n));
    out<<",\"tangent_applied\":";
    VectorJSON(out,m_outletLoad.tail(n));
    out<<",\"tangent_reaction\":[";
    bool first=true;
    // The unconstrained steady LHS minus prescribed RHS is the additional
    // outward traction action. A shared corner row belongs to both traces
    // and is exported once, without assigning it to either face.
    for(auto [id,shared]:m_outletReactionRows){
        if(!first)out<<',';
        first=false;
        const double reaction=current.flow[n+id]+advection.flow[n+id]-forcing.flow[n+id];
        out<<"{\"dof\":"<<id<<",\"value\":"<<reaction
           <<",\"scope\":\""<<(shared?"shared_corner":"outlet")<<"\"}";
    }
    out<<"]}";
}

void WriteFlowSolverFields(std::ostream &out, const CoupledOptions &options,
                           const CoupledStepInfo &info)
{
    using RANSWallSolveJson::Quote;
    const auto &work = info.flowWork;
    out << ",\"flow_preconditioner\":" << Quote(options.flowPreconditioner)
        << ",\"pressure_correctors\":"
        << (options.flowPreconditioner == "PressureCorrection" ? options.pressureCorrectors : 0)
        << ",\"picard_stop\":" << Quote(options.picardStop)
        << ",\"picard_advection\":" << Quote(options.picardAdvection)
        << ",\"picard_stress_linearization\":" << Quote(options.picardStressLinearization)
        << ",\"picard_outer_correctors\":" << options.picardOuterCorrectors
        << ",\"accepted_by_fixed_budget\":" << (info.acceptedByFixedBudget ? "true" : "false")
        << ",\"mixed_flow_work\":{\"setup_seconds\":" << work.setup
        << ",\"velocity_factor_seconds\":" << work.velocityFactor
        << ",\"pressure_assembly_seconds\":" << work.pressureAssembly
        << ",\"pressure_factor_seconds\":" << work.pressureFactor
        << ",\"predictor_seconds\":" << work.predictor
        << ",\"pressure_solve_seconds\":" << work.pressureSolve
        << ",\"velocity_correction_seconds\":" << work.correction
        << ",\"pressure_corrections\":" << work.pressureCorrections
        << ",\"pressure_krylov_iterations\":" << work.pressureKrylovIterations
        << ",\"fixed_frozen_residuals\":[" << work.fixedMomentum << ','
        << work.fixedContinuity << ',' << work.fixedGauge << ']'
        << ",\"krylov_iterations\":" << work.krylovIterations << '}'
        << ",\"scalar_factor_seconds\":{\"k\":" << info.timingDetail.kFactor
        << ",\"epsilon\":" << info.timingDetail.epsilonFactor << '}';
    if (options.advancement == "CoupledJFNKV1")
    {
        const auto &newton = options.stepNewton;
        out << ",\"step_advection\":"
            << Quote(options.UsesImplicitAdvection() ? "BDF2ImplicitAdvectionV1" : "BDF2EX2")
            << ",\"newton_forcing\":{\"initial\":" << newton.forcingInitial
            << ",\"minimum\":" << newton.forcingMinimum
            << ",\"maximum\":" << newton.forcingMaximum
            << ",\"gamma\":" << newton.forcingGamma
            << ",\"exponent\":" << newton.forcingExponent << '}';
    }
}

void Coupling::Report(std::ostream&out,const CoupledHistory&h,const CoupledStepInfo&info,bool coefficients)const
{
    out<<std::setprecision(17)<<"{\"schema\":\"WJH27NativeCoupledStateV2\",\"model\":\""<<WallModelV2<<"\",\"signature\":\""<<m_signature<<"\",\"source_identity\":\""<<m_options.sourceIdentity<<"\",\"reference_identity\":\""<<m_reference->Identity()<<"\",\"manufactured\":"<<(m_mms?"true":"false")<<",\"step\":"<<h.step<<",\"time\":"<<h.time<<",\"dt\":"<<h.dt<<",\"state_sha256\":\""<<StateHash(h.current)<<"\",\"previous_sha256\":\""<<StateHash(h.previous)<<"\",\"advection_sha256\":\""<<StateHash(h.advection)<<"\",\"previous_advection_sha256\":\""<<StateHash(h.previousAdvection)<<"\",\"sizes\":["<<m_flow.VelocitySize()<<','<<m_flow.PressureSize()<<','<<m_k.Size()<<','<<m_e.Size()<<"],\"gauge\":"<<(m_flow.Gauge()?"true":"false")<<",\"nonlinear_iterations\":"<<info.nonlinearIterations<<",\"linear_iterations\":["<<info.flowIterations<<','<<info.kIterations<<','<<info.epsilonIterations<<"],\"rebuilds\":"<<info.rebuilds<<",\"retries\":"<<info.retries<<",\"seconds\":"<<info.seconds<<",\"coefficient_drift\":"<<info.coefficientDrift<<",\"residuals\":["<<info.residual.momentum<<','<<info.residual.continuity<<','<<info.residual.gauge<<','<<info.residual.k<<','<<info.residual.epsilon<<','<<info.residual.trace<<"],\"minimum_k\":"<<info.minimumK<<",\"minimum_epsilon\":"<<info.minimumEpsilon<<",\"minimum_covariance_over_k\":"<<info.minimumCovarianceOverK<<",\"minimum_principal_over_mu\":"<<info.minimumPrincipalOverMu<<",\"disjoint_seconds\":{\"constitutive\":"<<info.timing.constitutive<<",\"wall\":"<<info.timing.wall<<",\"assembly\":"<<info.timing.assembly<<",\"flow\":"<<info.timing.flow<<",\"scalar\":"<<info.timing.scalar<<",\"rebuild\":"<<info.timing.rebuild<<",\"guards\":"<<info.timing.guards<<",\"forcing\":"<<info.timing.forcing<<",\"residual\":"<<info.timing.residual<<",\"nonlinear\":"<<info.timing.nonlinear<<",\"io\":"<<info.timing.io<<",\"other\":"<<info.timing.other<<'}';
    out<<",\"quadratureX\":"<<m_options.quadratureX<<",\"quadratureR\":"<<m_options.quadratureR<<",\"linear_tolerance\":"<<m_options.linearTolerance<<",\"nonlinear_tolerance\":"<<m_options.nonlinearTolerance<<",\"maximum_iterations\":"<<m_options.maximumIterations<<",\"mms_body_amplitude\":"<<m_options.mmsBodyAmplitude<<",\"retry_reason\":"<<RANSWallSolveJson::Quote(info.retryReason);
    out << ",\"timing_detail_seconds\":{\"momentum_advection\":" << info.timingDetail.momentumAdvection
        << ",\"scalar_advection\":" << info.timingDetail.scalarAdvection
        << ",\"scalar_momentum_coupling\":" << info.timingDetail.scalarMomentumCoupling
        << ",\"history_movement\":" << info.timingDetail.historyMovement << '}';
    WriteFlowSolverFields(out, m_options, info);
    out<<",\"inner_k_transport\":"<<RANSWallSolveJson::Quote(m_reference->TransportName())
       <<",\"boundary_identity\":"<<RANSWallSolveJson::Quote(m_reference->BoundaryIdentity());
    out<<",\"outlet_policy\":"<<RANSWallSolveJson::Quote(m_options.outletPolicy)
       <<",\"outlet_identity\":"<<RANSWallSolveJson::Quote(m_outletIdentity);
    const auto &p=m_reference->Parameters();
    const auto K=InterfaceBounds(h.current);
    out<<",\"interface_K_over_Kin\":["<<K[0]/p.KIn<<','<<K[1]/p.KIn<<']'
       <<",\"interface_K_table_margin_over_Kin\":["<<(K[0]-p.KMinimum)/p.KIn
       <<','<<(p.KMaximum-K[1])/p.KIn<<']'
       <<",\"interface_cutoff_range\":["<<m_reference->Cutoff(K[1])
       <<','<<m_reference->Cutoff(K[0])<<']';
    if(m_options.checkOutletActions)WriteOutletActions(out,h);
    const auto slopes=m_reference->LogCutoffDerivativeBounds();
    out<<",\"scalar_quadrature_policy\":\"CubicCutoffFacePanelsV1\",\"cutoff_radius_range\":["<<m_cutoffRange[0]<<','<<m_cutoffRange[1]
       <<"],\"cutoff_log_slope_bounds\":["<<slopes[0]<<','<<slopes[1]<<']';
    out<<",\"mms_variant\":"<<RANSWallSolveJson::Quote(m_options.mmsVariant);
    out<<",\"implicit_solver\":"<<RANSWallSolveJson::Quote(m_options.advancement)
        <<",\"coupled_newton\":{\"krylov_iterations\":"<<info.coupledKrylovIterations
        <<",\"jacobian_actions\":"<<info.coupledJacobianActions
        <<",\"residual_evaluations\":"<<info.coupledResidualEvaluations
        <<",\"preconditioner_applications\":"<<info.coupledPreconditionerApplications
        <<",\"rejected_trials\":"<<info.coupledRejectedTrials
        <<",\"counts_complete\":"<<(info.coupledStatisticsComplete?"true":"false")
        <<",\"step_seconds_cap\":"<<m_options.stepNewton.maximumSeconds
        <<",\"uncommitted_guard\":\"positivity/root/covariance\",\"physical_acceptance_guard\":"
        <<RANSWallSolveJson::Quote(m_options.picardStop == "FixedCorrectors"
             ? "fresh full planar principal, continuity, gauge and trace at fixed budget"
             : "fresh full planar principal and complete current residual");
    WriteStepPreconditionerFields(out,m_options,info);out<<'}';
    out<<",\"nonlinear_trace_scope\":"<<RANSWallSolveJson::Quote(m_options.advancement=="CoupledJFNKV1"?
        "Uncommitted Newton bases and Armijo candidates; FD excluded; records are not accepted updates; times are inclusive observations":
        "Completed Picard iterations and guarded acceleration proposals; times are inclusive observations");
    out<<",\"steady_velocity_preconditioner\":"<<RANSWallSolveJson::Quote(m_options.steadyVelocityPC);
    out<<",\"stress_split\":"<<RANSWallSolveJson::Quote(m_options.stressSplit)<<",\"all_row_action_crosschecks\":"<<(m_options.checkCurrentActions?"true":"false");
    out<<",\"current_action_checks\":{\"all_row_calls\":"<<m_actionChecks<<",\"maximum_scaled_differences_flow_k_epsilon\":["<<m_actionDifference[0]<<','<<m_actionDifference[1]<<','<<m_actionDifference[2]<<"],\"augmented_calls\":"<<m_augmentedChecks<<",\"maximum_augmented_scaled_difference\":"<<m_augmentedDifference<<",\"maximum_weak_divergence_load_norm\":"<<m_checkedWeakDivergence<<",\"maximum_gauge_multiplier\":"<<m_checkedGaugeMultiplier<<'}';
    if (m_geometry.normalMap && m_options.checkCurrentActions)
        out<<",\"dependent_trace_direction_relative_error\":"<<m_dependentTraceDerivativeError;
    if (m_geometry.normalMap && h.step==0)
        out<<",\"initial_continuity_before_projection\":"<<m_initialContinuityBefore
           <<",\"initial_continuity_after_projection\":"<<m_initialContinuityAfter;
    out<<",\"basis\":\"Modified_A\",\"anderson_depth\":"<<m_options.andersonDepth<<",\"acceleration\":{\"proposals\":"<<info.accelerationProposals<<",\"accepted\":"<<info.accelerationAccepted<<",\"inadmissible_trials\":"<<info.accelerationInadmissible<<'}';
    out<<",\"initial_field_transfer\":{\"packet_sha256\":"<<RANSWallSolveJson::Quote(m_initialPacketSha)<<",\"source_state_sha256\":"<<RANSWallSolveJson::Quote(m_initialSourceStateSha)<<",\"projected\":"<<(m_initialProjected?"true":"false")<<",\"source_was_candidate\":"<<(m_initialCandidate?"true":"false")<<'}';
    if (m_geometrySeedHeight)
    {
        out<<",\"geometry_seed\":{\"from_height\":"<<*m_geometrySeedHeight
           <<",\"to_height\":"<<m_geometry.normalMap->amplitude
           <<",\"nested_refinement\":"<<(m_geometrySeedRefined?"true":"false")
           <<",\"source_step\":"<<m_geometrySeedStep
           <<",\"source_time\":"<<m_geometrySeedTime
           <<",\"piola_flux_roundtrip_error\":"<<m_initialFluxRoundtripError
           <<",\"outer_taper_end_radius\":"<<(m_geometry.H+m_options.geometrySeedTaper)
           <<",\"matching_tangent_correction\":"<<m_initialTangentialTraceCorrection
           <<",\"continuity_before\":"<<m_initialContinuityBefore
           <<",\"continuity_after\":"<<m_initialContinuityAfter
           <<",\"foot_gradient_length\":"<<m_options.geometrySeedFootLength
           <<",\"collar_gradient_length\":"<<m_options.geometrySeedCollarLength
           <<",\"collar_gradient_reach\":"<<m_options.geometrySeedCollarReach
           <<",\"source_identity\":"<<RANSWallSolveJson::Quote(m_options.initialSourceIdentity.empty()
                ? m_options.sourceIdentity : m_options.initialSourceIdentity);
        if (m_initialGradientFit)
            out<<",\"weighted_gradient_rms_before_after\":["<<m_initialGradientFit->front()
               <<','<<m_initialGradientFit->back()<<']';
        out<<'}';
    }
    if(m_options.recordNonlinearTrace){out<<",\"nonlinear_trace\":[";for(int i=0;i<int(info.nonlinearTrace.size());++i){const auto&t=info.nonlinearTrace[i];if(i)out<<',';out<<"{\"residuals\":["<<t.residual.momentum<<','<<t.residual.continuity<<','<<t.residual.gauge<<','<<t.residual.k<<','<<t.residual.epsilon<<','<<t.residual.trace<<"],\"accelerated\":"<<(t.accelerated?"true":"false")<<",\"seconds\":"<<t.seconds<<'}';}out<<']';}

    out<<",\"constitutive_backend\":\""<<(m_options.backend==ConstitutiveBackend::Historical?"Historical":"FastFull3DV1")<<"\",\"constitutive_calls\":{\"values\":"<<info.constitutiveStats.values<<",\"directions\":"<<info.constitutiveStats.directions<<",\"fallbacks\":"<<info.constitutiveStats.fallbacks<<",\"root_iterations\":"<<info.constitutiveStats.rootIterations<<",\"maximum_normalized_root_residual\":"<<info.constitutiveStats.maximumRootResidual<<'}';
    if(m_mms){auto errors=MMSErrors(h.current,h.time);out<<",\"mms_rms_errors_U_Pi_k_epsilon_gradU\":[";for(int i=0;i<5;++i){if(i)out<<',';out<<errors[i];}out<<']';}
    if(coefficients){
        auto state=[&](const CoupledState&s){out<<"{\"flow\":";VectorJSON(out,s.flow);out<<",\"k\":";VectorJSON(out,s.k);out<<",\"epsilon\":";VectorJSON(out,s.epsilon);out<<'}';};
        WriteCoefficientFields(out,h.current,h.time);
        out<<",\"current_advection\":";state(h.advection);out<<",\"previous_advection\":";state(h.previousAdvection);
        WriteInputHistory(out);
    }
    out<<'}';
}

void Coupling::WriteCoefficientFields(std::ostream&out,const CoupledState&s,double time,bool refreshLift)const
{
    auto state=[&](const CoupledState&s){out<<"{\"flow\":";VectorJSON(out,s.flow);out<<",\"k\":";VectorJSON(out,s.k);out<<",\"epsilon\":";VectorJSON(out,s.epsilon);out<<'}';};
    auto integers=[&](const std::vector<int>&v){out<<'[';for(int i=0;i<int(v.size());++i){if(i)out<<',';out<<v[i];}out<<']';};
    auto numbers=[&](const std::vector<double>&v){out<<'[';for(int i=0;i<int(v.size());++i){if(i)out<<',';out<<v[i];}out<<']';};
    out<<",\"coefficients\":";state(s);
    if(refreshLift)
    {
        auto essential=DirichletState(s,time);
        out<<",\"essential_flow_lift\":";VectorJSON(out,m_flow.Lift(essential.flow));
        out<<",\"essential_k_lift\":";VectorJSON(out,m_k.Lift(essential.k));
        out<<",\"essential_epsilon_lift\":";VectorJSON(out,m_e.Lift(essential.epsilon));
    }
    else
        out<<",\"essential_flow_lift\":null,\"essential_k_lift\":null,\"essential_epsilon_lift\":null";
    out<<",\"mms_bulk_body_rhs\":";VectorJSON(out,m_bulkBodyAction);out<<",\"mms_sampled_body_rhs\":";VectorJSON(out,m_sampledBodyAction);
    out<<",\"geometry\":{\"d0\":"<<m_geometry.d0<<",\"H\":"<<m_geometry.H<<",\"h\":"<<m_geometry.h<<",\"left\":"<<m_left<<",\"right\":"<<m_right<<",\"top\":"<<m_top;
    if(m_geometry.normalMap)
    {
        const auto &map=*m_geometry.normalMap;
        out<<",\"mapping\":\"NormalCosineHillV1\",\"normal_map_signature\":"
           <<RANSWallSolveJson::Quote(map.Signature())
           <<",\"hill_height\":"<<map.amplitude<<",\"hill_half_length\":"<<map.lengthX
           <<",\"tunnel_top\":"<<map.tunnelTop
           <<",\"initial_slowdown\":"<<m_options.hillInitialSlowdown
           <<",\"initial_epsilon_bump\":"<<m_options.hillInitialEpsilonBump
           <<",\"coordinates\":\"affine chart (a,d0+rho); velocity components Cartesian\"";
    }
    out<<"},\"free_velocity_ids\":[";integers(m_flow.FreeVelocityIds()[0]);out<<',';integers(m_flow.FreeVelocityIds()[1]);out<<"],\"k_free_ids\":";integers(m_k.FreeIds());out<<",\"epsilon_free_ids\":";integers(m_e.FreeIds());
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

void Coupling::WriteInputHistory(std::ostream&out)const
{
    auto state=[&](const CoupledState&s){out<<"{\"flow\":";VectorJSON(out,s.flow);out<<",\"k\":";VectorJSON(out,s.k);out<<",\"epsilon\":";VectorJSON(out,s.epsilon);out<<'}';};
    if(m_hasInputHistory){const auto&a=m_lastInputHistory;out<<",\"input_history\":{\"step\":"<<a.step<<",\"time\":"<<a.time<<",\"dt\":"<<a.dt<<",\"current\":";state(a.current);out<<",\"previous\":";state(a.previous);out<<",\"advection\":";state(a.advection);out<<",\"previous_advection\":";state(a.previousAdvection);out<<'}';}
}

void Coupling::WriteRejection(const std::filesystem::path&path)const
{
    if(m_rejectedCandidate)WriteDiagnostic(path,*m_rejectedCandidate);
}

void Coupling::WriteSteadyDiagnostic(const std::filesystem::path&path,
                                    const std::filesystem::path&uncommittedPath)const
{
    if(m_steadyDiagnosticCandidate)
        WriteDiagnostic(path,m_steadyDiagnosticCandidate->candidate,&*m_steadyDiagnosticCandidate);
    if(m_steadyUncommitted)
        WriteDiagnostic(uncommittedPath,m_steadyUncommitted->candidate,&*m_steadyUncommitted);
}

void Coupling::WriteDiagnostic(const std::filesystem::path&path,const DiagnosticCandidate&candidate,
    const SteadyPacket *steady)const
{
    namespace J=RANSWallSolveJson;
    const auto *options=steady?&m_diagnosticOptions:nullptr;
    const bool timeStep=candidate.role==DiagnosticRole::RejectedTimeStep;
    const char *role=timeStep?"REJECTED_TIME_STEP":
        candidate.role==DiagnosticRole::RejectedInitialization?"REJECTED_INITIALIZATION":"STEADY_SEARCH";
    auto norms=[](std::ostream &out,const CoupledResidual&r){
        out<<'['<<r.momentum<<','<<r.continuity<<','<<r.gauge<<','<<r.k<<','<<r.epsilon<<','<<r.trace<<']';
    };
    std::ostringstream out;
    out<<std::setprecision(17)<<"{\"schema\":\"WJH27DiagnosticCandidateV1\",\"physical_history\":false"
       <<",\"role\":"<<J::Quote(role)<<",\"reason\":"<<J::Quote(candidate.reason)
       <<",\"attempted_step\":"<<candidate.attemptedStep<<",\"attempted_time\":"<<candidate.attemptedTime
       <<",\"attempted_iteration\":"<<candidate.attemptedIteration
       <<",\"residual_kind\":"<<J::Quote(timeStep?(candidate.residual?
            "refreshed original time step":"original time residual unavailable"):
            options?(options->eta>0.?"refreshed modified steady":"refreshed original steady"):"unevaluated preflight")
       <<",\"residuals\":";
    if(candidate.residual)norms(out,*candidate.residual);else out<<"null";
    out<<",\"initial_packet_sha256\":"<<J::Quote(m_initialPacketSha)
       <<",\"initial_coefficient_sha256\":"<<J::Quote(m_initialSourceStateSha);
    if(options)
    {
        const auto *original=&steady->original;
        out<<",\"iterate_custody\":"<<J::Quote(candidate.acceptedSteady?"LAST_ACCEPTED":"UNCOMMITTED_TRIAL");
        if(steady->checkpointSequence)
        {
            out<<",\"accepted_checkpoint\":{\"sequence\":"<<steady->checkpointSequence
               <<",\"refresh_seconds\":"<<steady->checkpointRefreshSeconds<<",\"session_sha256\":[";
            for(std::size_t i=0;i<options->checkpointSessionDigests.size();++i)
            {
                if(i)out<<',';
                out<<J::Quote(options->checkpointSessionDigests[i]);
            }
            out<<"]}";
        }
        std::ostringstream identity;
        identity<<std::setprecision(17)<<m_signature<<"|steady";
        if(options->eta>0.)
            identity<<"|CartesianBoundaryPolynomialV1|"<<options->eta<<'|'<<options->regularizationViscosity
                    <<'|'<<options->regularizationStart<<'|'<<options->regularizationFull;
        const std::string equationIdentity=D::Hash(identity.str());
        if(steady->meritController)
            out<<",\"merit_controller\":{\"policy_sha256\":"<<J::Quote(m_meritPolicy)
               <<",\"coefficient_sha256\":"<<J::Quote(StateHash(candidate.state))
               <<",\"equation_identity\":"<<J::Quote(equationIdentity)
               <<",\"entry\":"<<J::Quote(options->meritEntry)
               <<",\"state\":"<<steady->meritController->ToJson()<<'}';
        out<<",\"equation_identity\":"<<J::Quote(equationIdentity)
           <<",\"regularization\":{\"eta\":"<<options->eta
           <<",\"viscosity\":"<<options->regularizationViscosity
           <<",\"shape\":\"CartesianBoundaryPolynomialV1\",\"radius_ramp\":["
           <<options->regularizationStart<<','<<options->regularizationFull
           <<"],\"mean_energy_loss\":"<<steady->extraWork<<'}'
           <<",\"correction\":{\"initial_pseudo_time\":"<<options->pseudoTime
           <<",\"initial_diffusion\":"<<options->correctionDiffusion<<",\"rows\":[";
        for(std::size_t i=0;i<m_diagnosticCorrections.size();++i)
        {
            const auto &row=m_diagnosticCorrections[i];
            if(i)out<<',';
            out<<'['<<row[0]<<','<<row[1]<<','<<row[2]<<']';
        }
        out<<"]},\"original_residuals\":";
        norms(out,original->residual);
        out<<",\"base_coefficient_sha256\":"<<J::Quote(m_steadyBaseSha256)
           <<",\"steady_inverse\":"<<J::Quote(options->sparseInverse?"ColoredSparseV1":"BlockV1")
           <<",\"scale_extrema\":[";
        for(std::size_t i=0;i<m_steadyScaleExtrema.size();++i)
        {
            if(i)out<<',';
            out<<J::Numbers(m_steadyScaleExtrema[i]);
        }
        out<<"],\"accepted_bases\":[";
        for(std::size_t i=0;i<m_steadyBaseObservations.size();++i)
        {
            if(i)out<<',';
            out<<J::Numbers(m_steadyBaseObservations[i]);
        }
        out<<"],\"sparse_pattern_seconds\":"<<m_steadyPatternSeconds<<",\"sparse_inverse_work\":[";
        for(std::size_t i=0;i<m_sparseSteadyWork.size();++i)
        {
            const auto &work=m_sparseSteadyWork[i];
            if(i)out<<',';
            out<<"{\"base_norm\":"<<work.baseNorm<<",\"base_coefficient_sha256\":"<<J::Quote(work.baseSha256)
               <<",\"inverse_mass\":"<<work.inverseMass<<",\"stiffness\":"<<work.stiffness<<",\"colors\":"<<work.colors
               <<",\"factor_fd_relative_step\":"<<work.fdRelativeStep
               <<",\"backsolve_audit_seconds\":"<<work.backsolveAuditSeconds
               <<",\"backsolve_scaled_shifted_then_unshifted\":"
               <<(work.backsolveAudit?J::Numbers(*work.backsolveAudit):"null")
               <<",\"nonzeros\":"<<work.nonzeros<<",\"factor_nonzeros\":"<<work.factorNonzeros
               <<",\"residual_evaluations\":"<<work.residualEvaluations
               <<",\"rejected_samples\":"<<work.rejectedSamples
               <<",\"assembly_seconds\":"<<work.assemblySeconds<<",\"factor_seconds\":"<<work.factorSeconds
               <<",\"maximum_unassigned_action\":"<<work.maximumUnassignedAction
               <<",\"last_rejection\":"<<J::Quote(work.lastRejection)<<'}';
        }
        out<<"],\"direction_audit_seconds\":"<<m_steadyAuditSeconds
           <<",\"direction_audit_residual_evaluations\":"<<m_steadyAuditResiduals<<",\"direction_audit\":[";
        for(std::size_t i=0;i<m_steadyActionAudits.size();++i)
        {
            const auto &record=m_steadyActionAudits[i];
            if(i)out<<',';
            out<<"{\"direction\":"<<J::Quote(record.direction)<<",\"step\":"<<record.step;
            if(record.differenceBlocks)
                out<<",\"taylor_remainder\":"<<record.taylorRemainder
                   <<",\"centered_relative_error\":"<<record.centeredRelativeError
                   <<",\"difference_blocks\":"<<J::Numbers(*record.differenceBlocks)
                   <<",\"shift_relative_error\":"<<(options->pseudoTime>0.||options->correctionDiffusion>0.?
                       J::Number(record.shiftRelativeError):"null");
            else out<<",\"rejection\":"<<J::Quote(record.rejection);
            out<<'}';
        }
        out<<"],\"epsilon_direction_seconds\":"<<m_steadyEpsilonSeconds
           <<",\"epsilon_direction_evaluations\":"<<m_steadyEpsilonEvaluations
           <<",\"epsilon_remainders\":[";
        for(std::size_t i=0;i<m_steadyEpsilonRemainders.size();++i)
        {
            const auto &r=m_steadyEpsilonRemainders[i];
            if(i)out<<',';
            out<<"{\"iteration\":"<<r.iteration<<",\"trial\":"<<r.trial<<",\"step\":"<<r.step;
            if(r.unavailable.empty())
                out<<",\"dominant_scaled_epsilon_row\":"<<r.epsilonId<<",\"owners\":"<<J::Numbers(r.owners)
                   <<",\"advection_diffusion_production_destruction\":"<<J::Numbers(r.terms)
                   <<",\"sum_minus_full_remainder\":"<<r.closure
                   <<",\"dependent_lift_residual_overlap\":"<<r.liftResidual
                   <<",\"dependent_lift_coefficient_remainder\":"<<r.liftCoefficient
                   <<",\"half_spacing_derivative_difference\":"<<(r.halfSpacingDifference?J::Number(*r.halfSpacingDifference):"null");
            else out<<",\"unavailable\":"<<J::Quote(r.unavailable);
            out<<'}';
        }
        out<<']';
        out<<",\"minimum_k\":"<<original->minimumK<<",\"minimum_epsilon\":"<<original->minimumEpsilon
           <<",\"minimum_covariance_over_k\":"<<original->minimumCovarianceOverK
           <<",\"original_principal_over_mu\":"<<original->minimumPrincipalOverMu
           <<",\"modified_principal_over_original_mu\":"<<steady->modifiedMargin
           <<",\"scalar_bernstein_depth\":"<<options->scalarSubdivisionDepth
           <<",\"physical_whole_cell_scalar_lower_bounds\":["
           <<steady->wholeCellLower[0]<<','<<steady->wholeCellLower[1]<<']'
           <<",\"dependent_trace_direction_relative_error\":"<<m_dependentTraceDerivativeError;
    }
    if(timeStep)
    {
        const auto &h=m_lastInputHistory;
        out<<",\"input_history_identity\":{\"signature\":"<<J::Quote(m_signature)
           <<",\"step\":"<<h.step<<",\"time\":"<<h.time<<",\"dt\":"<<h.dt
           <<",\"current_sha256\":"<<J::Quote(StateHash(h.current))
           <<",\"previous_sha256\":"<<J::Quote(StateHash(h.previous))
           <<",\"advection_sha256\":"<<J::Quote(StateHash(h.advection))
           <<",\"previous_advection_sha256\":"<<J::Quote(StateHash(h.previousAdvection))<<'}';
    }
    const auto &s=candidate.state;
    const auto bounds=InterfaceBounds(s,false);
    const auto &p=m_reference->Parameters();
    const bool inTable=bounds[0]>=p.KMinimum&&bounds[1]<=p.KMaximum;
    out<<",\"state\":{\"schema\":\"WJH27NativeDiagnosticStateV1\",\"model\":"<<J::Quote(WallModelV2)
       <<",\"signature\":"<<J::Quote(m_signature)<<",\"source_identity\":"<<J::Quote(m_options.sourceIdentity)
       <<",\"reference_identity\":"<<J::Quote(m_reference->Identity())
       <<",\"boundary_identity\":"<<J::Quote(m_reference->BoundaryIdentity())
       <<",\"state_sha256\":"<<J::Quote(StateHash(s))<<",\"basis\":\"Modified_A\""
       <<",\"manufactured\":"<<(m_mms?"true":"false")<<",\"mms_variant\":"<<J::Quote(m_options.mmsVariant)
       <<",\"gauge\":"<<(m_flow.Gauge()?"true":"false")
       <<",\"sizes\":["<<m_flow.VelocitySize()<<','<<m_flow.PressureSize()<<','<<m_k.Size()<<','<<m_e.Size()<<']'
       <<",\"outlet_policy\":"<<J::Quote(m_options.outletPolicy)<<",\"outlet_identity\":"<<J::Quote(m_outletIdentity)
       <<",\"quadratureX\":"<<m_options.quadratureX<<",\"quadratureR\":"<<m_options.quadratureR
       <<",\"nonlinear_tolerance\":"<<m_options.nonlinearTolerance
       <<",\"dt\":"<<m_options.dt<<",\"step_advection\":"
       <<J::Quote(m_options.UsesImplicitAdvection()?"BDF2ImplicitAdvectionV1":"BDF2EX2")
       <<",\"interface_K_over_Kin\":["<<bounds[0]/p.KIn<<','<<bounds[1]/p.KIn<<']'
       <<",\"interface_K_table_margin_over_Kin\":["<<(bounds[0]-p.KMinimum)/p.KIn
       <<','<<(p.KMaximum-bounds[1])/p.KIn<<']'
       <<",\"interface_derived_status\":"<<J::Quote(inTable?"AVAILABLE":"UNAVAILABLE_OUTSIDE_TABLE")
       <<",\"interface_cutoff_range\":";
    if(inTable)out<<'['<<m_reference->Cutoff(bounds[1])<<','<<m_reference->Cutoff(bounds[0])<<']';
    else out<<"null";
    WriteCoefficientFields(out,s,candidate.attemptedTime,inTable);
    if(timeStep)WriteInputHistory(out);
    out<<"}}\n";

    const auto bytes=out.str();
    PublishDiagnostic(path,[&](int fd){D::WriteAll(fd,bytes.data(),bytes.size());});
}

}
