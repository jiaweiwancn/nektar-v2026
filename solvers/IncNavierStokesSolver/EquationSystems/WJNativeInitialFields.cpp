// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeInitialFields.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <boost/property_tree/json_parser.hpp>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace Nektar::WJNative
{
namespace
{
using Tree=boost::property_tree::ptree;
void Require(bool value,const std::string&reason)
{if(!value)throw std::invalid_argument("H27 initial fields: "+reason);}
double Tolerance(double a,double b,double scale)
{return 64*std::numeric_limits<double>::epsilon()*std::max({std::abs(a),std::abs(b),scale});}
bool Same(double a,double b,double scale)
{return std::isfinite(a)&&std::isfinite(b)&&std::abs(a-b)<=Tolerance(a,b,scale);}
bool SHA(const std::string&s)
{return s.size()==64&&std::all_of(s.begin(),s.end(),[](char c){return (c>='0'&&c<='9')||(c>='a'&&c<='f');});}
double Number(const Tree&t,const std::string&path)
{const double v=t.get<double>(path);Require(std::isfinite(v),"nonfinite "+path);return v;}
void UniqueKeys(const Tree&t)
{
    std::set<std::string>keys;
    for(const auto&item:t){if(!item.first.empty())Require(keys.insert(item.first).second,"duplicate JSON key "+item.first);UniqueKeys(item.second);}
}
template<class T>std::vector<T>Array(const Tree&t,const std::string&path)
{
    const auto&child=t.get_child(path);std::vector<T>out;out.reserve(child.size());
    for(const auto&item:child){Require(item.first.empty()&&item.second.empty(),"non-array/non-scalar "+path);out.push_back(item.second.get_value<T>());}
    return out;
}
Eigen::VectorXd Coefficients(const Tree&t,const std::string&path,Eigen::Index size)
{
    const auto values=Array<double>(t,path);
    Require(values.size()==static_cast<std::size_t>(size),"coefficient size differs: "+path);
    Eigen::VectorXd result(size);for(Eigen::Index i=0;i<size;++i)result[i]=values[i];
    Require(result.allFinite(),"nonfinite coefficients: "+path);return result;
}
InitialFieldGeometry Geometry(const Tree&t)
{
    InitialFieldGeometry g;g.d0=Number(t,"d0");g.H=Number(t,"H");g.h=Number(t,"h");
    g.left=Number(t,"left");g.right=Number(t,"right");g.top=Number(t,"top");
    g.curve=t.get<double>("curve",0.);
    g.hillHalfLength=t.get<double>("hill_half_length",0.);
    g.normalMap=t.get<std::string>("normal_map_signature","");
    Require(std::isfinite(g.curve)&&g.curve==0&&g.H>0&&g.h>=0&&g.h<g.H&&g.top>g.H&&g.right>g.left,
            "unsupported nonflat/invalid geometry");return g;
}
bool SameGeometry(const InitialFieldGeometry&a,const InitialFieldGeometry&b)
{
    const double scale=std::max(a.H,b.H);
    return scale>0&&Same(a.d0,b.d0,scale)&&Same(a.H,b.H,scale)&&Same(a.h,b.h,scale)&&
        Same(a.left,b.left,scale)&&Same(a.right,b.right,scale)&&Same(a.top,b.top,scale)&&
        std::isfinite(a.curve)&&std::isfinite(b.curve)&&a.curve==0&&b.curve==0&&
        a.normalMap==b.normalMap;
}
std::array<double,4>Modified(double t)
{
    const double b=t*(1-t);
    return {1-t,t,b,b*(4*t-2)};
}
const std::vector<int>&Ids(const InitialFieldElement&e,unsigned field)
{return field<2?e.u:field==2?e.k:field==3?e.epsilon:e.p;}
const std::vector<double>&Signs(const InitialFieldElement&e,unsigned field)
{return field<2?e.us:field==2?e.ks:field==3?e.epsilonSigns:e.ps;}
}

InitialFields::InitialFields(const std::string &path, bool allowAcceptedTimeState,
    bool allowDiagnostic,bool allowUncommittedGuess)
{
    Require(std::filesystem::is_regular_file(path),"input is not a regular coefficient packet");
    Require(std::filesystem::file_size(path)<=256ULL*1024*1024,"coefficient packet exceeds declared 256 MiB input allowance");
    std::ifstream file(path,std::ios::binary);Require(bool(file),"cannot open coefficient packet");
    const std::string bytes((std::istreambuf_iterator<char>(file)),{});
    Require(file.eof()||bool(file),"coefficient packet read failed");
    Tree root;std::istringstream input(bytes);boost::property_tree::read_json(input,root);UniqueKeys(root);
    auto&m=m_metadata;m.packetSha256=RANSWallR3::Hash(bytes);m.schema=root.get<std::string>("schema");
    if(m.schema=="WJH27RegisteredCoupledRunV2"){
        Require(root.get<bool>("normal_main")&&root.get<std::string>("factory")=="WJNativeCoupledV2",
                "registered packet is not from the declared normal-main factory");
        if (root.get<bool>("steady_mode"))
        {
            Require(root.get<bool>("steady_solver.converged"),
                    "registered initial reference has not converged in steady mode; use explicit candidate schema for a candidate");
            Require(root.get<std::string>("steady_solver.status")=="CONVERGED",
                    "contradictory registered steady status");
            m.role=InitialFieldRole::QualifiedReference;
            m.qualifiedReference=true;
        }
        else
        {
            Require(allowAcceptedTimeState,"positive-time coefficients require a geometry-continuation import");
            m.role=InitialFieldRole::AcceptedTimeState;
        }
    }else if(m.schema=="WJH27DiagnosticCandidateV1"){
        Require(allowDiagnostic,"diagnostic coefficient packets require explicit diagnostic steady preparation");
        const auto role=root.get<std::string>("role");
        Require(role=="REJECTED_INITIALIZATION"||role=="REJECTED_TIME_STEP"||role=="STEADY_SEARCH",
                "unsupported diagnostic coefficient role");
        Require(!root.get<bool>("physical_history"),"diagnostic packet falsely claims a physical history");
        m.role=InitialFieldRole::Diagnostic;
        const auto custody=root.get<std::string>("iterate_custody","");
        Require(custody.empty()||custody=="LAST_ACCEPTED"||
                (custody=="UNCOMMITTED_TRIAL"&&allowUncommittedGuess&&!allowAcceptedTimeState),
                "uncommitted steady trial requires an explicit fixed-height diagnostic guess");
        const auto checkpoint=root.get_child_optional("accepted_checkpoint");
        if(checkpoint)
        {
            const auto digests=Array<std::string>(*checkpoint,"session_sha256");
            Require(role=="STEADY_SEARCH"&&custody=="LAST_ACCEPTED"&&
                    checkpoint->get<std::size_t>("sequence")>0&&!digests.empty()&&
                    std::all_of(digests.begin(),digests.end(),SHA),"invalid accepted diagnostic checkpoint");
        }
        m.originalDiagnosticSteady=!checkpoint&&role=="STEADY_SEARCH"&&custody!="UNCOMMITTED_TRIAL"&&
            root.get<std::string>("residual_kind","")=="refreshed original steady"&&
            root.get<double>("regularization.eta",-1.)==0.;
        m.equationIdentity=root.get<std::string>("equation_identity","");
    }else{
        Require(m.schema=="WJH27SteadyCandidateV1","unsupported initial coefficient packet schema");
        Require(!root.get<bool>("qualified_reference"),"candidate packet falsely claims reference qualification");
        m.role=InitialFieldRole::Candidate;m.qualifiedReference=false;
    }
    const auto&state=root.get_child("state");m.stateSchema=state.get<std::string>("schema");
    const bool diagnostic=m.role==InitialFieldRole::Diagnostic;
    Require(m.stateSchema==(diagnostic?"WJH27NativeDiagnosticStateV1":"WJH27NativeCoupledStateV2"),
            "unsupported native state schema");
    m.step=diagnostic?root.get<std::uint64_t>("attempted_step"):state.get<std::uint64_t>("step");
    m.time=diagnostic?Number(root,"attempted_time"):Number(state,"time");
    if (m.role==InitialFieldRole::AcceptedTimeState)
    {
        const double dt=Number(state,"dt");
        Require(m.step>0 && dt>0. && m.time>0. &&
                std::abs(m.time-m.step*dt)<=2e-13*std::max(1.,m.time),
                "geometry seed has inconsistent positive-time counters");
        Require(!state.get<bool>("manufactured") && !state.get<bool>("accepted_by_fixed_budget") &&
                state.get<std::string>("step_advection")=="BDF2ImplicitAdvectionV1" &&
                bool(state.get_child_optional("input_history")),
                "geometry seed requires an accepted physical implicit BE/BDF2 state");
        const double tolerance=Number(state,"nonlinear_tolerance");
        const auto residuals=Array<double>(state,"residuals");
        Require(tolerance>0. && residuals.size()==6 &&
                std::all_of(residuals.begin(),residuals.end(),[tolerance](double value) {
                    return std::isfinite(value) && value>=0. && value<=std::min(tolerance,1e-10);
                }),"geometry seed has an unqualified equation residual");
        Require(Number(state,"minimum_k")>0. && Number(state,"minimum_epsilon")>0. &&
                Number(state,"minimum_covariance_over_k")>1e-10 &&
                Number(state,"minimum_principal_over_mu")>1e-8,
                "geometry seed has unqualified admissibility margins");
    }
    else if(!diagnostic)
    {
        Require(m.step==0 && m.time==0.,"steady initial fields require step=time=0");
    }
    m.basis=state.get<std::string>("basis");Require(m.basis=="Modified_A","only native Modified_A Q3/Q2 initial fields are supported");
    m.outletPolicy=state.get<std::string>("outlet_policy","FullTractionV1");
    m.model=state.get<std::string>("model");m.sourceIdentity=state.get<std::string>("source_identity");
    m.referenceIdentity=state.get<std::string>("reference_identity");m.mmsVariant=state.get<std::string>("mms_variant");
    m.manufactured=state.get<bool>("manufactured");m.gauge=state.get<bool>("gauge");
    m.stateSha256=state.get<std::string>("state_sha256");
    m.signature=state.get<std::string>("signature","");
    Require(!m.model.empty()&&!m.mmsVariant.empty()&&SHA(m.sourceIdentity)&&SHA(m.referenceIdentity)&&SHA(m.stateSha256),
            "missing/invalid scientific metadata or identities");
    if(const auto controller=root.get_child_optional("merit_controller"))
    {
        auto &c=*controller;
        Require(diagnostic&&root.get<std::string>("role")=="STEADY_SEARCH"&&
                root.get<std::string>("iterate_custody","")=="LAST_ACCEPTED"&&
                root.get<std::string>("residual_kind","")=="refreshed original steady"&&
                root.get<double>("regularization.eta",-1.)==0.&&m.step==0&&m.time==0.&&
                c.get<std::string>("coefficient_sha256")==m.stateSha256&&
                c.get<std::string>("equation_identity")==m.equationIdentity,
                "invalid original-merit checkpoint coefficients or equation binding");
        m.meritPolicy=c.get<std::string>("policy_sha256");
        Require(SHA(m.meritPolicy)&&c.get<std::string>("state.schema")==SteadyMeritState::Identity,
                "invalid original-merit policy binding or schema");
        SteadyMeritState r;
        r.usedMass=Number(c,"state.used_mass");
        r.nextMass=Number(c,"state.next_mass");
        r.referenceNorm=Number(c,"state.reference_norm");
        r.previousNorm=Number(c,"state.previous_norm");
        r.acceptedNorm=Number(c,"state.accepted_norm");
        r.usedForcing=Number(c,"state.used_forcing");
        r.nextForcing=Number(c,"state.next_forcing");
        r.step=Number(c,"state.step");
        if(c.get<std::string>("state.actual_to_predicted_decrease")!="null")
            r.ratio=Number(c,"state.actual_to_predicted_decrease");
        r.accepted=c.get<int>("state.accepted");
        r.goodSteps=c.get<int>("state.good_steps");
        r.decision=c.get<std::string>("state.decision");
        if(c.get_child_optional("state.radius"))
        {
            const auto schema=c.get<std::string>("state.radius.schema");
            Require(schema==SteadyBoundedOptions::Identity||schema==SteadyBoundedOptions::RatioIdentity,
                    "invalid bounded-direction checkpoint schema");
            auto &a=r.radius.emplace();
            a.modelRatio=schema==SteadyBoundedOptions::RatioIdentity;
            a.used=Number(c,"state.radius.used");
            a.next=Number(c,"state.radius.next");
            a.displacement=Number(c,"state.radius.scaled_displacement");
            a.windowBase=Number(c,"state.radius.window_base");
            a.windowUpdates=c.get<int>("state.radius.window_updates");
            a.firstRate=Number(c,"state.radius.first_rate");
            a.lastRate=Number(c,"state.radius.last_rate");
            a.stagnated=c.get<bool>("state.radius.stagnated");
            if(r.accepted>0)
            {
                const auto delta=Array<double>(c,"accepted_displacement");
                Require(!delta.empty()&&RANSWallR3::Hash(delta.data(),delta.size()*sizeof(double))==
                        c.get<std::string>("accepted_displacement_sha256"),"bounded accepted displacement identity differs");
                a.acceptedDisplacement=Eigen::Map<const Eigen::VectorXd>(delta.data(),delta.size());
                Require(a.acceptedDisplacement.allFinite(),"nonfinite bounded accepted displacement");
            }
        }
        m.meritController=r;
    }
    m_geometry=Geometry(state.get_child("geometry"));
    const auto sizes=Array<int>(state,"sizes");Require(sizes.size()==4&&std::all_of(sizes.begin(),sizes.end(),[](int n){return n>0;}),"invalid coefficient sizes");
    std::copy(sizes.begin(),sizes.end(),m_sizes.begin());
    const auto flowSize=2LL*m_sizes[0]+m_sizes[1];
    Require(flowSize<=std::numeric_limits<Eigen::Index>::max(),"coefficient size overflows native vector index");
    m_flow=Coefficients(state,"coefficients.flow",static_cast<Eigen::Index>(flowSize));
    m_k=Coefficients(state,"coefficients.k",m_sizes[2]);m_epsilon=Coefficients(state,"coefficients.epsilon",m_sizes[3]);
    RANSWallR3::Sha256 hash;for(const auto*v:{&m_flow,&m_k,&m_epsilon})hash.Update(v->data(),v->size()*sizeof(double));
    Require(hash.Finish()==m.stateSha256,"native coefficient state SHA256 differs");
    for(const auto&item:state.get_child("elements")){
        Require(item.first.empty(),"elements is not an array");const auto&t=item.second;InitialFieldElement e;
        e.geometry=t.get<int>("geometry");Require(e.geometry>=0,"negative geometry ID");
        const auto bounds=Array<double>(t,"bounds");Require(bounds.size()==4,"invalid element bounds inventory");
        e.xl=bounds[0];e.xr=bounds[1];e.rl=bounds[2];e.rr=bounds[3];
        Require(std::all_of(bounds.begin(),bounds.end(),[](double v){return std::isfinite(v);}),"nonfinite element bounds");
        if(auto basis=t.get_optional<std::string>("basis"))Require(*basis==m.basis,"element basis differs from declared native basis");
        e.u=Array<int>(t,"u");e.us=Array<double>(t,"us");e.p=Array<int>(t,"p");e.ps=Array<double>(t,"ps");
        e.k=Array<int>(t,"k");e.ks=Array<double>(t,"ks");e.epsilon=Array<int>(t,"epsilon");e.epsilonSigns=Array<double>(t,"epsilon_signs");
        m_elements.push_back(std::move(e));
    }
    ValidateInventory();ValidateContinuity();
}

void InitialFields::ValidateOriginalSteadyRelease(const std::string &receiptPath,
    const std::string &receiptSha,const std::string &source,const std::string &reference,
    const std::string &signature,int quadratureX,int quadratureR,double hillHalfLength)const
{
    const auto &m=m_metadata;
    Require(m.originalDiagnosticSteady&&!m.manufactured&&m.step==0&&m.time==0.&&
            m.sourceIdentity==source&&m.referenceIdentity==reference&&m.signature==signature&&
            m.equationIdentity==RANSWallR3::Hash(signature+"|steady"),
            "physical release requires the exact current zero-regularization steady equations");
    auto pinned=[](const std::string &path,const std::string &sha) {
        Require(SHA(sha)&&std::filesystem::is_regular_file(path)&&
                std::filesystem::file_size(path)<=4*1024*1024,"invalid original steady audit receipt");
        std::ifstream input(path,std::ios::binary);
        const std::string bytes((std::istreambuf_iterator<char>(input)),{});
        Require(!input.bad()&&RANSWallR3::Hash(bytes)==sha,"original steady audit receipt digest differs");
        Tree tree;
        std::istringstream json(bytes);
        boost::property_tree::read_json(json,tree);
        UniqueKeys(tree);
        return tree;
    };
    const Tree receipt=pinned(receiptPath,receiptSha);
    Require(receipt.get<std::string>("schema")=="WJH27OriginalSteadyReleaseV1"&&
            receipt.get<std::string>("packet.sha256")==m.packetSha256&&
            receipt.get<std::string>("coefficient_sha256")==m.stateSha256,
            "release receipt binds different candidate coefficients");
    auto audit=[&](const std::string &key) {
        const Tree tree=pinned(receipt.get<std::string>(key+".path"),receipt.get<std::string>(key+".sha256"));
        Require(tree.get<std::string>("native.sha256")==m.packetSha256&&
                tree.get<std::string>("table.sha256")==reference,
                "independent audit binds a different packet or table");
        return tree;
    };
    auto limit=[](const Tree &tree,const std::string &key,double tolerance) {
        Require(std::abs(Number(tree,key))<=tolerance,"independent release gate failed: "+key);
    };
    auto original=[&](const Tree &tree) {
        Require(tree.get<std::string>("candidate_role")=="STEADY_SEARCH"&&
                tree.get<std::string>("regularization")=="null"&&tree.get<int>("step")==0&&
                Number(tree,"time")==0.&&tree.get<int>("order")>quadratureX&&
                tree.get<int>("radial_order")>quadratureR,
                "release requires richer independent original steady audits");
    };
    const Tree equations=audit("equations"),response=audit("response"),points=audit("points");
    original(equations);original(response);
    Require(equations.get<std::string>("source_identity")==source&&equations.get<bool>("residual_pass"),
            "independent original equation receipt did not pass");
    limit(equations,"residuals.momentum",1e-9);
    limit(equations,"residuals.continuity",1e-9);
    limit(equations,"residuals.k",1e-8);
    limit(equations,"residuals.epsilon",1e-8);
    limit(equations,"epsilon_trace.interior_projection_defect",1e-10);
    limit(equations,"interface.physical_face_measure_max_difference",1e-10);
    limit(equations,"interface.shared_trace_max",1e-10);
    limit(equations,"power.algebraic_error",1e-10);
    limit(equations,"mass_integration_identity_error",1e-10);
    Require(Number(equations,"integration_minimum_jacobian")>0.&&
            Number(equations,"raw_boundary.outlet_minimum_normal_velocity")>=-1e-12,
            "independent geometry or outlet receipt failed");
    Require(response.get<bool>("one_sided_faces"),"release response audit must include one-sided faces");
    Require(response.get<std::string>("source_identity")==source&&
            response.get<std::string>("coefficient_sha256")==m.stateSha256,
            "response audit source or coefficient digest differs");
    for(const std::string group:{"inner","outer","matching_inner","matching_outer",
                                 "faces_inner","faces_outer","cutoff","forcing_support"})
    {
        Require(bool(response.get_child_optional("regions."+group)),
                "missing independent original response for "+group);
        const auto &region=response.get_child("regions."+group);
        const double uncertainty=Number(region,"max_directional_difference");
        Require(region.get<int>("samples")>0&&region.get<bool>("sample_pass")&&uncertainty>=0.&&
                Number(region,"principal_over_mu")-uncertainty>1e-8&&
                Number(region,"covariance_over_k")>1e-10&&Number(region,"minimum_k")>0.&&
                Number(region,"minimum_epsilon")>0.,
                "independent original response failed for "+group);
    }
    Require(points.get<std::string>("source_identity")==source&&
            points.get<std::string>("coefficient_sha256")==m.stateSha256&&
            points.get<int>("order")>std::max(quadratureX,quadratureR),
            "point audit source, coefficient digest or rule differs");
    const double H=m_geometry.H,L=m_geometry.hillHalfLength,h=m_geometry.h,left=m_geometry.left;
    const double witnessScale=std::max(1.,H);
    Require(Same(L,hillHalfLength,H)&&L>0.&&left<-L&&m_geometry.right>L,
            "different physical hill geometry for release controls");
    // Mirror the named independent controls, using the packet's exact support radius.
    const std::map<std::string,std::array<double,2>> controls={
        {"upstream_collar",{left+.228*(-L-left),H*(1111./1500.)}},
        {"downstream_outer",{L,H+std::min(.957899,.5*(m_geometry.top-H))}},
        {"support_original",{.75*L,h}},
        {"support_refined",{.5625*L,h}},
        {"support_downstream_foot",{L,h}}};
    // The XML chart faces and the exact reference h can differ by a few ulps
    // after decimal I/O and subtracting d0. Keep both one-sided witnesses at h.
    auto containsWitness=[&](const InitialFieldElement &cell,double a,double r) {
        return a>=cell.xl-Tolerance(a,cell.xl,witnessScale)&&a<=cell.xr+Tolerance(a,cell.xr,witnessScale)&&
               r>=cell.rl-Tolerance(r,cell.rl,witnessScale)&&r<=cell.rr+Tolerance(r,cell.rr,witnessScale);
    };
    auto owner=[&](const Tree &point) {
        const int i=point.get<int>("geometry");
        const auto found=m_index.find(i);
        Require(found!=m_index.end()&&containsWitness(m_elements[found->second],
                Number(point,"a"),Number(point,"rho")),"response point has a different owner");
        return i;
    };
    auto positive=[](const Tree &point,const std::string &label) {
        const double uncertainty=Number(point,"differentiation_spread");
        const auto &minima=point.get_child("minima");
        Require(uncertainty>=0.&&minima.size()>=3,"missing independent differentiation spacings for "+label);
        Require(Number(point,"k")>0.&&Number(point,"epsilon")>0.&&
                Number(point,"covariance_over_k")>1e-10&&Number(point,"mu")>0.,
                "independent point domain failed for "+label);
        for(const auto &minimum:minima)
            Require(Number(minimum.second,"over_mu")-uncertainty>1e-8,
                    "independent original point response failed for "+label);
    };
    std::set<std::pair<std::string,int>> fixed;
    for(const auto &entry:points.get_child("fixed_controls"))
    {
        const auto &point=entry.second;
        const int i=owner(point);
        const auto label=point.get_optional<std::string>("control");
        if(label)
        {
            const auto required=controls.find(*label);
            Require(required!=controls.end()&&Same(Number(point,"a"),required->second[0],witnessScale)&&
                    Same(Number(point,"rho"),required->second[1],witnessScale)&&fixed.emplace(*label,i).second,
                    "fixed original response location or owner differs for "+*label);
        }
        positive(point,label?*label:"additional fixed control");
    }
    for(const auto &[label,coordinate]:controls)
    {
        bool hasOwner=false;
        for(const auto &e:m_elements)
            if(containsWitness(e,coordinate[0],coordinate[1]))
            {
                hasOwner=true;
                Require(fixed.count({label,e.geometry})==1,"missing fixed original response owner for "+label);
            }
        Require(hasOwner,"fixed original response location has no owner for "+label);
    }
    std::set<std::string> checked;
    for(const auto &report:points.get_child("local_minima"))
    {
        Require(report.second.get<std::string>("receipt.sha256")==receipt.get<std::string>("response.sha256"),
                "local response audit has a different volume scan");
        const auto region=report.second.get<std::string>("region");
        Require(checked.insert(region).second,"duplicate local response audit for "+region);
        const auto &values=report.second.get_child("values");
        Require(!values.empty(),"empty local response audit");
        bool hasReportedWitness=false;
        for(const auto &entry:values)
        {
            const auto &point=entry.second;
            const int i=owner(point);
            positive(point,region);
            if(region=="forcing_support"&&point.get<std::string>("side")=="reported-minimum")
            {
                const auto &witness=response.get_child("regions.forcing_support.witness");
                hasReportedWitness|=i==witness.get<int>("geometry")&&
                    Same(Number(point,"a"),Number(witness,"a"),witnessScale)&&
                    Same(Number(point,"rho"),Number(witness,"rho"),witnessScale);
            }
        }
        Require(region!="forcing_support"||hasReportedWitness,
                "missing reported one-sided forcing support witness");
    }
    for(const std::string region:{"inner","outer","faces_inner","faces_outer","forcing_support"})
        Require(checked.count(region)>0,"missing local original response audit for "+region);
}

void InitialFields::ValidateInventory()
{
    Require(!m_elements.empty(),"empty element inventory");
    std::sort(m_elements.begin(),m_elements.end(),[](const auto&a,const auto&b){return a.geometry<b.geometry;});
    std::array<std::set<int>,4>visited;
    const auto&g=m_geometry;
    for(std::size_t j=0;j<m_elements.size();++j){const auto&e=m_elements[j];
        Require(m_index.emplace(e.geometry,j).second,"duplicate geometry ID");
        Require(e.xr-e.xl>Tolerance(e.xl,e.xr,g.H)&&e.rr-e.rl>Tolerance(e.rl,e.rr,g.H),"degenerate element bounds");
        Require(e.xl>=g.left-Tolerance(e.xl,g.left,g.H)&&e.xr<=g.right+Tolerance(e.xr,g.right,g.H)&&
                e.rl>=-Tolerance(e.rl,0,g.H)&&e.rr<=g.top+Tolerance(e.rr,g.top,g.H),"element outside declared rectangle");
        const bool outer=e.rl>=g.H-Tolerance(e.rl,g.H,g.H);
        Require(outer||e.rr<=g.H+Tolerance(e.rr,g.H,g.H),"element straddles dependent H interface");
        Require(outer==!e.epsilon.empty(),"epsilon inventory is not exactly the outer domain");
        for(unsigned f:{0U,4U,2U,3U}){
            const auto&ids=Ids(e,f);const auto&signs=Signs(e,f);
            const int slot=f==0?0:f==4?1:f==2?2:3;
            const std::size_t expected=f==3&&!outer?0:f==4?9:16;
            Require(ids.size()==expected&&signs.size()==expected,"local Q3/Q2 signed map dimension differs");
            std::set<int>local;
            for(std::size_t i=0;i<ids.size();++i){Require(ids[i]>=0&&ids[i]<m_sizes[slot]&&local.insert(ids[i]).second,
                    "local signed map has repeated/out-of-range coefficient ID");
                Require(signs[i]==1.||signs[i]==-1.,"signed map sign is not exactly +/-1");visited[slot].insert(ids[i]);}
        }
        m_xBreaks.push_back(e.xl);m_xBreaks.push_back(e.xr);
    }
    for(int f=0;f<4;++f)Require(visited[f].size()==static_cast<std::size_t>(m_sizes[f]),"global coefficient inventory has unused/missing IDs");
    std::sort(m_xBreaks.begin(),m_xBreaks.end());std::vector<double>unique;
    for(double x:m_xBreaks)if(unique.empty()||!Same(x,unique.back(),g.H))unique.push_back(x);
    m_xBreaks=std::move(unique);
    Require(m_xBreaks.size()>=2&&Same(m_xBreaks.front(),g.left,g.H)&&Same(m_xBreaks.back(),g.right,g.H),"x cover differs from declared rectangle");
    for(std::size_t i=0;i+1<m_xBreaks.size();++i){
        const double x=m_xBreaks[i]+(m_xBreaks[i+1]-m_xBreaks[i])/2;
        std::vector<std::size_t>slab;
        for(std::size_t j=0;j<m_elements.size();++j)if(x>m_elements[j].xl&&x<m_elements[j].xr)slab.push_back(j);
        Require(!slab.empty(),"uncovered x slab");
        std::sort(slab.begin(),slab.end(),[&](auto a,auto b){return m_elements[a].rl<m_elements[b].rl;});
        Require(Same(m_elements[slab.front()].rl,0,g.H)&&Same(m_elements[slab.back()].rr,g.top,g.H),"radial slab does not cover full rectangle");
        for(std::size_t j=1;j<slab.size();++j)Require(Same(m_elements[slab[j-1]].rr,m_elements[slab[j]].rl,g.H),"rectangular cells overlap or leave a gap");
        m_slabs.push_back(std::move(slab));
    }
}

bool InitialFields::Contains(const InitialFieldElement&e,double x,double r)const
{
    const double scale=m_geometry.H;
    return x>=e.xl-Tolerance(x,e.xl,scale)&&x<=e.xr+Tolerance(x,e.xr,scale)&&
           r>=e.rl-Tolerance(r,e.rl,scale)&&r<=e.rr+Tolerance(r,e.rr,scale);
}

double InitialFields::Value(const InitialFieldElement&e,unsigned field,double x,double r)const
{
    const auto&ids=Ids(e,field);const auto&signs=Signs(e,field);
    Require(!ids.empty(),"requested owner has no transported epsilon field");
    const auto bx=Modified((x-e.xl)/(e.xr-e.xl)),br=Modified((r-e.rl)/(e.rr-e.rl));
    const int width=field==4?3:4;
    const Eigen::VectorXd&coefficient=field==2?m_k:field==3?m_epsilon:m_flow;
    const Eigen::Index offset=field==1?m_sizes[0]:field==4?2LL*m_sizes[0]:0;
    double value=0;
    for(int j=0;j<width;++j)for(int i=0;i<width;++i){const int local=i+width*j;value+=bx[i]*br[j]*signs[local]*coefficient[offset+ids[local]];}
    Require(std::isfinite(value),"nonfinite evaluated coefficient field");return value;
}

double InitialFields::Evaluate(unsigned field,double x,double physicalZ,int ownerGeometry)const
{
    Require(field<5&&std::isfinite(x)&&std::isfinite(physicalZ),"invalid requested field or coordinate");
    const auto&g=m_geometry;double r=physicalZ-g.d0;
    Require(std::isfinite(r),"nonfinite physical radius");
    for(double edge:{g.left,g.right})if(Same(x,edge,g.H))x=edge;
    for(double edge:{0.,g.H,g.top})if(Same(r,edge,g.H))r=edge;
    Require(x>=g.left&&x<=g.right&&r>=0&&r<=g.top,"requested point lies outside declared domain");
    Require(field!=3||r>=g.H,"transported epsilon requested in the inner domain");
    if(ownerGeometry>=0){auto owner=m_index.find(ownerGeometry);Require(owner!=m_index.end(),"requested geometry owner is absent");
        const auto&e=m_elements[owner->second];Require(Contains(e,x,r),"requested point lies outside its explicit owner");return Value(e,field,x,r);}
    Require(ownerGeometry==-1,"invalid optional owner ID");
    const InitialFieldElement*selected=nullptr;
    auto include=[&](std::size_t slab){for(auto index:m_slabs[slab]){const auto&e=m_elements[index];
        if((field!=3||!e.epsilon.empty())&&Contains(e,x,r)&&(!selected||e.geometry<selected->geometry))selected=&e;}};
    auto upper=std::upper_bound(m_xBreaks.begin(),m_xBreaks.end(),x);
    std::size_t slab=upper==m_xBreaks.begin()?0:std::min<std::size_t>(upper-m_xBreaks.begin()-1,m_slabs.size()-1);
    include(slab);
    if(slab>0&&Same(x,m_xBreaks[slab],g.H))include(slab-1);
    if(slab+1<m_slabs.size()&&Same(x,m_xBreaks[slab+1],g.H))include(slab+1);
    Require(selected!=nullptr,"no native element owns the requested point/field");
    return Value(*selected,field,x,r);
}

void InitialFields::ValidateContinuity()const
{
    const double eps=4096*std::numeric_limits<double>::epsilon();
    std::array<double,5>scale{m_flow.head(m_sizes[0]).cwiseAbs().maxCoeff(),
        m_flow.segment(m_sizes[0],m_sizes[0]).cwiseAbs().maxCoeff(),m_k.cwiseAbs().maxCoeff(),
        m_epsilon.cwiseAbs().maxCoeff(),m_flow.tail(m_sizes[1]).cwiseAbs().maxCoeff()};
    for(std::size_t ia=0;ia<m_elements.size();++ia)for(std::size_t ib=ia+1;ib<m_elements.size();++ib){
        const auto&a=m_elements[ia];const auto&b=m_elements[ib];
        bool vertical=false,horizontal=false;double af=0,bf=0,lo=0,hi=0;
        if(Same(a.xr,b.xl,m_geometry.H)){vertical=true;af=a.xr;bf=b.xl;}
        else if(Same(a.xl,b.xr,m_geometry.H)){vertical=true;af=a.xl;bf=b.xr;}
        if(vertical){lo=std::max(a.rl,b.rl);hi=std::min(a.rr,b.rr);}
        else{
            if(Same(a.rr,b.rl,m_geometry.H)){horizontal=true;af=a.rr;bf=b.rl;}
            else if(Same(a.rl,b.rr,m_geometry.H)){horizontal=true;af=a.rl;bf=b.rr;}
            if(horizontal){lo=std::max(a.xl,b.xl);hi=std::min(a.xr,b.xr);}
        }
        if((!vertical&&!horizontal)||hi-lo<=Tolerance(lo,hi,m_geometry.H))continue;
        for(unsigned field=0;field<5;++field){if(field==3&&(a.epsilon.empty()||b.epsilon.empty()))continue;
            // Four distinct points determine equality of a Q3 face trace.
            for(double t:{0.,1./3,2./3,1.}){const double q=lo+(hi-lo)*t;
                const double av=Value(a,field,vertical?af:q,vertical?q:af),bv=Value(b,field,vertical?bf:q,vertical?q:bf);
                Require(std::abs(av-bv)<=eps*std::max({scale[field],std::abs(av),std::abs(bv)}),
                    "discontinuous signed coefficient trace between geometry "+std::to_string(a.geometry)+" and "+std::to_string(b.geometry)+" field "+std::to_string(field));}
        }
    }
}

void InitialFields::ValidateBindings(const std::string&source,const std::string&reference,const std::string&model,
    bool manufactured,const std::string&variant,const InitialFieldGeometry&geometry)const
{
    Require(m_metadata.sourceIdentity==source,"source identity binding differs");
    Require(m_metadata.referenceIdentity==reference,"reference identity binding differs");
    Require(m_metadata.model==model,"model binding differs");
    Require(m_metadata.manufactured==manufactured,"manufactured/physical binding differs");
    Require(m_metadata.mmsVariant==variant,"MMS variant binding differs");
    Require(SameGeometry(m_geometry,geometry),"physical geometry binding differs");
}

bool InitialFields::CoversRefinement(const std::vector<InitialFieldElement>&elements)const
{
    if(elements.size()<=m_elements.size())return false;
    for(const auto &cell:elements)
    {
        const auto owner=std::find_if(m_elements.begin(),m_elements.end(),[&](const auto &source){
            return cell.epsilon.empty()==source.epsilon.empty() &&
                Contains(source,cell.xl,cell.rl) && Contains(source,cell.xr,cell.rr);
        });
        if(owner==m_elements.end())return false;
    }
    return true;
}

bool InitialFields::SameInventory(const std::vector<InitialFieldElement>&elements,const std::array<int,4>&sizes)const
{
    if(sizes!=m_sizes||elements.size()!=m_elements.size())return false;
    std::set<int>seen;
    for(const auto&a:elements){auto position=m_index.find(a.geometry);if(position==m_index.end()||!seen.insert(a.geometry).second)return false;
        const auto&b=m_elements[position->second];
        if(!Same(a.xl,b.xl,m_geometry.H)||!Same(a.xr,b.xr,m_geometry.H)||!Same(a.rl,b.rl,m_geometry.H)||!Same(a.rr,b.rr,m_geometry.H)||
            a.u!=b.u||a.us!=b.us||a.p!=b.p||a.ps!=b.ps||a.k!=b.k||a.ks!=b.ks||a.epsilon!=b.epsilon||a.epsilonSigns!=b.epsilonSigns)return false;
    }
    return true;
}
}
