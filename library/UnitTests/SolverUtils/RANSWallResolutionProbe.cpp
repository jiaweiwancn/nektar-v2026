///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallPhysicsSignature.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallKBoundProjection.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <boost/property_tree/json_parser.hpp>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <sstream>

using namespace Nektar;
using namespace Nektar::SolverUtils;
using namespace Nektar::RANSWallAlgebra;
namespace J=Nektar::RANSWallSolveJson;
using Tree=boost::property_tree::ptree;
namespace
{
Vector Read(const Tree &packet,const std::string &name)
{Vector r;for(const auto &entry:packet.get_child(name)){r.push_back(entry.second.get_value<double>());}return r;}
std::string Maps(const RANSWallEvolutionPatch::Diagnostic &values)
{std::vector<std::pair<std::string,std::string>> r;for(const auto &[name,v]:values){r.push_back({name,J::Numbers(v)});}return J::Object(r);}
std::string Wide(const std::map<std::string,std::vector<long double>> &values)
{
    std::vector<std::pair<std::string,std::string>> r;
    for(const auto &[name,v]:values)
    {
        std::vector<std::string> row;
        for(long double x:v){std::ostringstream s;s<<std::setprecision(std::numeric_limits<long double>::max_digits10)<<x;row.push_back(J::Quote(s.str()));}
        r.push_back({name,J::List(row)});
    }
    return J::Object(r);
}
std::string Field(const RANSWallEvolutionPatch &patch,int block)
{
    std::vector<std::string> elements;const auto field=patch.fields[block];const auto map=field->Field()->GetLocalToGlobalMap();
    for(int e=0;e<field->Elements().size();++e)
    {
        const auto &element=field->Elements()[e];Vector ids,signs;
        const int start=field->Field()->GetCoeff_Offset(e),count=field->Field()->GetExp(e)->GetNcoeffs();
        for(int j=0;j<count;++j){ids.push_back(map->GetLocalToGlobalMap(start+j));signs.push_back(map->GetLocalToGlobalSign(start+j));}
        elements.push_back(J::Object({{"geometry",J::Number(element.geometry)},{"degree",J::Number(element.degree)},
            {"lower",J::Numbers(element.lower)},{"upper",J::Numbers(element.upper)},{"ids",J::Numbers(ids)},{"sign",J::Numbers(signs)}}));
    }
    return J::Object({{"size",J::Number(patch.size[block])},{"mesh",J::List(elements)}});
}
Vector MapParameters(const RANSWallNormalMap &map)
{return {double(map.kind),double(map.dimension),map.height,map.heightSlopeX,map.heightSlopeY,map.heightAmplitude,
    map.heightWaveX,map.heightWaveY,map.slopeX,map.slopeY,map.radius,map.amplitude,map.lengthX,map.lengthY};}

#include "RANSWallOpenWeakChecks.h"
#include "RANSWallFlatInletChecks.h"

Vector ConstructContact(const RANSWallEvolutionPatch &patch,const Vector &history)
{
    Require(patch.Dimension()==2 && !patch.Mapped(),"constructed contact requires flat2D");
    const auto field=patch.fields[patch.K()];const auto map=field->Field()->GetLocalToGlobalMap();
    const double K=patch.SampleK(history,{0,0,1});std::map<int,double> values;
    Require(patch.lower[0]==-4 && patch.upper[0]==2 && patch.Parameters().geometry.sampleHeight==1,
        "constructed contact is preregistered for the original x bounds and sample height");
    for(int e=0;e<field->Elements().size();++e)
    {
        const auto &element=field->Elements()[e];Require(element.degree==3,"constructed polynomial requires p3");
        const int start=field->Field()->GetCoeff_Offset(e);std::array<int,4> f{},one{1,1,0,0},h{};
        if(element.lower[0]==-4 && element.upper[0]==-3)f={0,1,1,0};
        else if(element.lower[0]==1 && element.upper[0]==2)f={1,0,1,0};
        else{Require(element.lower[0]>=-3 && element.upper[0]<=1,"unexpected horizontal contact cell");f=one;}
        if(element.lower[2]==0){Require(std::abs(element.upper[2]-.01)<1e-16,"contact shape must use the existing first inner interface");h={1,0,1,0};}
        else{Require(element.lower[2]>=.01-1e-16,"nonconforming contact shape interval");}
        for(int j=0;j<16;++j)
        {
            const int ix=j%4,ir=j/4,id=map->GetLocalToGlobalMap(start+j);
            const double sign=map->GetLocalToGlobalSign(start+j);
            const double value=sign*K*(one[ix]*one[ir]-f[ix]*h[ir]);
            const auto old=values.find(id);Require(old==values.end() || old->second==value,"inconsistent exact contact coefficient");values[id]=value;
        }
    }
    Vector result=history;for(auto [id,value]:values){result[patch.offset[patch.K()]+id]=value;}
    return result;
}

std::string Rate(const RANSWallEvolutionPatch &patch,const RANSWallScalarMassRate &rate,int block,
    const std::vector<RANSWallCoordinate> &points,const Vector &effectiveState)
{
    std::vector<std::string> samples;
    for(const auto &point:points)
    {
        const auto row=patch.Point(block,point,block==patch.Epsilon() && std::abs(point[2]-patch.Parameters().geometry.height)<1e-11?1:0);
        std::vector<std::pair<std::string,std::string>> terms;
        for(const auto &[name,v]:rate.componentRates){terms.push_back({name,J::Number(row.Value(v))});}
        const double value=row.Value(effectiveState,patch.offset[block]);
        samples.push_back(J::Object({{"coordinate",J::Numbers(point)},{"value",J::Number(value)},
            {"exact_native_zero",J::Boolean(value==0)},{"rate",J::Number(row.Value(rate.rate))},{"component_rates",J::Object(terms)}}));
    }
    std::vector<std::pair<std::string,std::string>> components;
    for(const auto &[name,v]:rate.componentRates){components.push_back({name,J::Numbers(v)});}
    return J::Object({{"block",J::Number(block)},{"rate",J::Numbers(rate.rate)},{"boundary_rate",J::Numbers(rate.boundaryRate)},
        {"spatial_load",J::Numbers(rate.spatialLoad)},{"free_rows",J::Numbers(rate.freeRows)},
        {"fixed_rows",J::Numbers(rate.fixedRows)},{"free_mass_residual",J::Numbers(rate.freeMassResidual)},
        {"component_coefficient_rates",J::Object(components)},{"physical_samples",J::List(samples)}});
}

struct TraceRates
{
    bool verified=true;Vector epsilonBoundary;std::vector<std::string> nodes;
    std::string reason;
};
TraceRates QuietTraceRate(const RANSWallEvolutionPatch &patch,const Vector &state,const Vector &kRate,
    double time,const RANSWallEvolutionData &data)
{
    Require(!patch.AcceptedAWF(),"legacy quiet-chi trace derivative is unsupported for accepted AWF");
    TraceRates result;result.epsilonBoundary.assign(patch.size[patch.Epsilon()],0.);
    const double H=patch.Parameters().geometry.height,delta=patch.RobustParameters()[3];
    const auto field=patch.fields[patch.Epsilon()];
    const auto face=[&](int axis,int side,const RANSWallGlobalField::Element &e)
        {return axis==2 && side==0 && std::abs(e.lower[2]-H)<1e-11;};
    auto derivative=[&](const RANSWallCoordinate &p,bool baseline)
    {
        const auto current=patch.StateProfile(state,p,time,data);const auto transition=patch.Transition(state,p,time,data);
        // Reconstruct the actual chi=0 profile from the same constructor
        // inputs as Profile(K,p,0), including its physical panel interfaces.
        // A reweighted active profile would retain an extra cutoff and would
        // therefore not be the baseline discrete normalization operator.
        auto options=current->Options();options.requestedWeight=0;
        const auto zero=std::make_shared<EquilibriumAWFProfiles>(current->Geometry(),current->K(),
            patch.Parameters().model,current->Metric(),options);
        const auto profile=baseline?zero:current;
        const auto norm=patch.WallNormalization(state,p,*profile);
        const bool quiet=transition.state==1 && transition.admissibility==1 && norm.b>delta && norm.lambda==1 &&
            H>profile->Original().DissipationCutoff() && H>profile->Reference().cutoff;
        if(!quiet){result.verified=false;result.reason="Quiet/unit/uncapped trace assumptions not established at all projection nodes";}
        auto sample=p;sample[2]=patch.Parameters().geometry.sampleHeight;
        const double dotK=patch.SampleKRow(sample).Value(kRate);
        const double E=profile->Dissipation(H);
        result.nodes.push_back(J::Object({{"coordinate",J::Numbers(p)},{"baseline",J::Boolean(baseline)},
            {"chi",J::Number(transition.chi)},{"state_gate",J::Number(transition.state)},
            {"admissibility_gate",J::Number(transition.admissibility)},{"d2",J::Number(transition.d2)},
            {"b",J::Number(norm.b)},{"lambda",J::Number(norm.lambda)},{"K",J::Number(profile->K())},
            {"pre_E_H",J::Number(E)},{"sampled_K_rate",J::Number(dotK)}}));
        return 1.5*E/profile->K()*dotK;
    };
    const auto baseline=field->Boundary([&](const RANSWallCoordinate &p){return derivative(p,true);},face);
    auto correction=field->Boundary([&](const RANSWallCoordinate &p){return derivative(p,false)-derivative(p,true);},face);
    const auto map=field->Field()->GetLocalToGlobalMap();
    if(!patch.UniformEPFlat())for(int e=0;e<field->Elements().size();++e)
    {
        const auto &element=field->Elements()[e];
        if(std::abs(element.lower[2]-H)>1e-11 || element.lower[0]<patch.TransitionDistances()[2]-1e-11)continue;
        const int start=field->Field()->GetCoeff_Offset(e);
        for(int j=0;j<=element.degree;++j){correction.at(map->GetLocalToGlobalMap(start+j))=0.;}
    }
    for(auto [id,value]:baseline){result.epsilonBoundary[id]=value+correction.at(id);}
    return result;
}
}
int main(int argc,char **argv)
{
    try
    {
        Require(argc==4,"usage: RANSWallResolutionProbe SESSION.xml INPUT.json OUTPUT.json");
        Tree packet;boost::property_tree::read_json(argv[2],packet);
        auto world=LibUtilities::SessionReader::CreateInstance(2,argv);world->InitSession();
        char program[]="wall-resolution-probe";char *args[]{program,nullptr};
        std::vector<std::string> filenames{argv[1]};
        auto serial=LibUtilities::SessionReader::CreateInstance(1,args,filenames,LibUtilities::CommSerial::create(1,args));serial->InitSession();
        serial->SetSolverInfo("GlobalSysSoln","DirectFull");const auto graph=SpatialDomains::MeshGraphIO::Read(serial);
        RANSWallPatchParameters p;
        serial->LoadParameter("RANSWallChi",p.chi,1.);serial->LoadParameter("RANSWallHeight",p.geometry.height,2.);
        serial->LoadParameter("RANSWallRoughness",p.geometry.roughness,.1);serial->LoadParameter("Kinvis",p.geometry.molecularNu,1.5e-5);
        serial->LoadParameter("RANSWallSampleHeight",p.geometry.sampleHeight,1.);serial->LoadParameter("RANSWallTailStart",p.geometry.tailStart,1.);
        serial->LoadParameter("RANSWallTangentialOrder",p.tangentialOrder,8);serial->LoadParameter("RANSWallNormalOrder",p.normalOrder,16);
        std::string model;serial->LoadSolverInfo("TurbulenceModel",model,"CLSHighRe");
        if(model=="WJHighRe")LoadWJNativeConfiguration(serial,p);
        RANSWallEvolutionPatch patch(serial,graph,world->GetComm(),p,LoadRANSWallNormalMap(serial,2,p.geometry.height));
        const auto operation=packet.get<std::string>("operation","");
        Require(patch.Dimension()==2 && (!patch.Mapped() || (patch.OpenBoundary() && (operation=="open_checks" || operation=="initializer_checks"))) && world->GetComm()->GetSize()==1 && patch.MatchingEnabled(),"resolution scope is a supported serial flat2D wall formulation");
        Require(patch.InnerStressName()==packet.get<std::string>("inner_stress_model"),"packet law differs from session");
        for(const std::string name:{"BodyForceInner","BodyForceOuter","ScalarSourceInner","ScalarSourceOuter"})
            Require(!serial->DefinesFunction(name),"R3 resolution input must retain unforced physical data");
        double us=0;serial->LoadParameter("RANSWallFrictionVelocity",us,.4);
        RANSWallEvolutionData data(2);Vector reference,formedHistory;std::string replayContext="null";
        if(operation=="initializer_checks")
        {
            std::ofstream out(argv[3]);out<<FlatInletChecks(patch,serial,graph,world->GetComm(),argv[3])<<'\n';
            Require(out.good(),"initializer fixture write failed");world->Finalise();return 0;
        }
        if(operation=="open_checks")
        {
            Require(patch.OpenBoundary(),"open checks require OpenTunnelV1");
            patch.LoadOpenBoundaryData(serial,data);
            const auto reportPath=packet.get<std::string>("native_report");
            const auto reportHash=RANSWallR3::HashFile(reportPath);
            Require(reportHash==packet.get<std::string>("native_report_sha256"),"open producing report SHA256 differs");
            Tree native;boost::property_tree::read_json(reportPath,native);
            Require(native.get<std::string>("failure").empty() || packet.get<bool>("allow_failed_producer",false),
                "open readback of an accepted predecessor in a failed report requires an explicit diagnostic role");
            const auto signature=RANSWallPhysicsSignature(patch,serial,0.,false,us,{1.,0.},false);
            Require(signature==native.get<std::string>("signature"),"open native physical signature differs");
            const auto cpPath=packet.get<std::string>("native_checkpoint");
            const auto checkpointHash=RANSWallR3::HashFile(cpPath+"/State.bin");
            Require(checkpointHash==packet.get<std::string>("checkpoint_sha256"),"open checkpoint bytes differ");
            double dt;serial->LoadParameter("TimeStep",dt);
            const auto cp=RANSWallTimeIntegration::ReadCheckpoint(cpPath,signature,1,dt,patch.Size(),patch.DifferentialSize());
            const auto &state=cp.coupledState;const double time=cp.integration.time;
            const auto producerSteps=native.get<std::size_t>("steps");
            const auto &accepted=native.get_child("accepted_states"),&packets=native.get_child("wj_audit_packets");
            Require(!accepted.empty() && accepted.size()==packets.size() && producerSteps+1>=accepted.size(),"open accepted lineage shape differs");
            const auto first=producerSteps+1-accepted.size();
            Require(first<=cp.integration.steps && cp.integration.steps<=producerSteps,"open checkpoint outside producer accepted lineage");
            auto ai=accepted.begin(),pi=packets.begin();std::advance(ai,cp.integration.steps-first);std::advance(pi,cp.integration.steps-first);
            Vector acceptedState;for(const auto &v:ai->second)acceptedState.push_back(v.second.get_value<double>());
            Require(RANSWallR3::Hash(state)==RANSWallR3::Hash(acceptedState) && time==pi->second.get<double>("time"),"open selected state/time producer differs");
            Require(RANSWallR3::Hash(cp.fixedReference)==RANSWallR3::Hash(Read(native,"fixed_reference")),"open fixed reference differs");
            Require(RANSWallR3::Hash(Vector(state.begin(),state.begin()+patch.DifferentialSize()))==RANSWallR3::Hash(cp.integration.state),
                "open native integrator primary state differs from coupled fields");
            std::vector<std::string> checkpointHistory;
            for(const auto &slot:cp.integration.history)checkpointHistory.push_back(J::Quote(RANSWallR3::Hash(slot)));
            if(cp.integration.steps==producerSteps)
            {
                Require(time==native.get<double>("time") && RANSWallR3::Hash(state)==RANSWallR3::Hash(Read(native,"state")),"open terminal state differs");
                std::size_t slot=0;for(const auto &v:native.get_child("native_history"))
                {Vector h;for(const auto &x:v.second)h.push_back(x.second.get_value<double>());
                    Require(slot<cp.integration.history.size() && RANSWallR3::Hash(h)==RANSWallR3::Hash(cp.integration.history[slot++]),"open native GLM history differs");}
                Require(slot==cp.integration.history.size() && Read(native,"native_history_times")==cp.integration.times,"open native history time labels differ");
            }
            patch.RestoreReference(cp.fixedReference,cp.fixedDefect,cp.eta,data,time);
            Vector history=state;double alpha=0;
            if(cp.integration.steps>0)
            {
                const auto &audit=pi->second;history=Read(audit,"callback_history");alpha=audit.get<double>("alpha_mass");
                Require(history.size()==state.size() && alpha==1/dt && native.get<int>("implicit_callback_count")>0 &&
                    audit.get<std::string>("state_sha256")==RANSWallR3::Hash(state) &&
                    audit.get<std::string>("stage_history_sha256")==RANSWallR3::Hash(history),"open accepted audit/callback binding differs");
                if(cp.integration.steps==producerSteps)Require(native.get<double>("last_stage_time")==time &&
                    native.get<double>("last_stage_alpha")==alpha && RANSWallR3::Hash(Read(native,"last_stage_history"))==RANSWallR3::Hash(history),
                    "open completed callback must come from its producing report, not a zero-step restart");
            }
            std::string diagnostics,assembly,reuse="null",domain,weak="null";
            Vector spatial,stage,spatialNorms,stageNorms;
            double amplitudeK=0,amplitudeEpsilon=0;
            {
                RANSWallEvolutionPatch::EvaluationScope evaluation(patch,data,time);
                patch.ValidateOutlet(state,"immutable_readback");
                domain=patch.ScalarDomainJSON(patch.ScalarDomain(state,time,data,true),true);
                std::vector<RANSWallCoordinate> locations;
                if(auto requested=packet.get_child_optional("physical_observations"))for(const auto &entry:*requested)
                {RANSWallCoordinate q{};int i=0;for(const auto &v:entry.second){Require(i<3,"physical observation shape");q[i++]=v.second.get_value<double>();}Require(i==3,"physical observation shape");locations.push_back(q);}
                diagnostics=patch.OpenBoundaryDiagnostic(state,history,alpha,time,data,locations);
                spatial=patch.Residual(state,state,0.,time,data,false);
                stage=alpha>0?patch.Residual(state,history,alpha,time,data,false):Vector{};
                const double H=p.geometry.height,Q=p.model.Q(),kappa=p.model.Kappa();
                amplitudeK=patch.SampleK(cp.fixedReference,{patch.lower[0],0,H});
                amplitudeEpsilon=std::pow(amplitudeK/Q,1.5)/(kappa*H);
                const auto constraints=patch.Constraints(state,time,data);
                auto physicalNorms=[&](const Vector &r){auto v=patch.Norms(r,constraints,1.);
                    v[0]/=H;v[1]/=H;v[2]*=Q/(H*amplitudeK);v[3]/=(kappa*H*H*amplitudeEpsilon);v[4]/=H;v[5]/=H;
                    v.back()=0;for(auto [row,target]:constraints){(void)target;int b=0;while(b+1<patch.offset.size() && row>=patch.offset[b+1])++b;
                        v.back()=std::max(v.back(),std::abs(r[row])/(b<2?1.:b==patch.K()?amplitudeK:amplitudeEpsilon));}return v;};
                spatialNorms=physicalNorms(spatial);if(alpha>0)stageNorms=physicalNorms(stage);
                if(packet.get<bool>("weak",false))weak=Weak(patch,state,history,alpha,time,data,packet.get<int>("weak_outer_order",24));
                assembly=packet.get<bool>("assembly",true)?patch.WJTerrainAssemblyChecks(state,history,alpha>0?alpha:1/dt,time,data):"null";
            }
            if(packet.get<bool>("reuse",false))reuse=patch.EvaluationReuseChecks(state,history,time,data);
            Require(patch.InitializationWork()==Vector({0.,0.}) && patch.InitialVelocityCalls()==0 && reportHash==RANSWallR3::HashFile(reportPath) &&
                checkpointHash==RANSWallR3::HashFile(cpPath+"/State.bin"),"open passive reader changed saved bytes or lifted state");
            std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote(operation)},
                {"report_sha256",J::Quote(reportHash)},{"checkpoint_sha256",J::Quote(checkpointHash)},
                {"selected_absolute_step",J::Number(cp.integration.steps)},{"producer_absolute_step",J::Number(producerSteps)},
                {"producer_failure",J::Quote(native.get<std::string>("failure"))},
                {"checkpoint_native_history_sha256",J::List(checkpointHistory)},
                {"checkpoint_native_history_times",J::Numbers(cp.integration.times)},
                {"state_sha256",J::Quote(RANSWallR3::Hash(state))},{"callback_history_sha256",J::Quote(RANSWallR3::Hash(history))},
                {"open_boundary",diagnostics},{"weak_checks",weak},{"scalar_domain",domain},{"assembly_checks",assembly},{"reuse_checks",reuse},
                {"spatial_residual",J::Numbers(spatial)},{"actual_completed_stage_residual",J::Numbers(stage)},
                {"physical_mass_dual_norms",J::Object({{"uncancelled_spatial",J::Numbers(spatialNorms)},
                    {"completed_stage",J::Numbers(stageNorms)},{"fixed_amplitudes_U_k_epsilon",J::Numbers(Vector{1.,amplitudeK,amplitudeEpsilon})},
                    {"time_scale_H_h_over_U_infinity","1"},{"order",J::Quote("u,v,k,entire_outer_epsilon,inner_continuity,outer_continuity,absent_gauge,constraints")},
                    {"scope",J::Quote("consistent-mass dual passive algebra; no governing solve, synthetic callback or state projection")}})},
                {"initialization_work",J::Numbers(patch.InitializationWork())},
                {"initial_velocity_mode",J::Quote(patch.InitialVelocityMode())},{"initial_velocity_calls",J::Number(patch.InitialVelocityCalls())},{"governing_solves","0"}})<<'\n';
            Require(out.good(),"open immutable checks write failed");world->Finalise();return 0;
        }
        if(operation=="archived_report_stages")
        {
            // Older short runs retained only periodic checkpoints. This is a
            // distinct read-only report/packet lineage contract, never a
            // restart or a waiver of the exact-checkpoint spatial reader.
            Require(patch.AcceptedAWF(),"archived native-report scope is accepted AWF");
            const auto nativeFile=packet.get<std::string>("native_report");
            Require(RANSWallR3::HashFile(nativeFile)==packet.get<std::string>("native_report_sha256"),"archived native report SHA256 differs");
            Tree native,attempt;boost::property_tree::read_json(nativeFile,native);
            const auto packetFile=packet.get<std::string>("archived_packet");
            Require(RANSWallR3::HashFile(packetFile)==packet.get<std::string>("archived_packet_sha256"),"archived packet SHA256 differs");
            boost::property_tree::read_json(packetFile,attempt);
            Require(!native.get<std::string>("failure").empty(),"report bridge requires a terminal native report");
            const auto &last=native.get_child("pimple_attempts").back().second;
            Require(last.get<std::string>("k_snapshot_sha256")==packet.get<std::string>("archived_packet_sha256") &&
                last.get<bool>("scheduled_complete"),"native report does not identify the completed attempted packet");
            Require(RANSWallR3::HashFile(attempt.get<std::string>("matrix_file"))==attempt.get<std::string>("matrix_sha256"),"archived matrix SHA256 differs");
            double dt,eta,dx,dy;serial->LoadParameter("TimeStep",dt);serial->LoadParameter("RANSWallEta",eta,0.);
            serial->LoadParameter("RANSWallDirectionX",dx,1.);serial->LoadParameter("RANSWallDirectionY",dy,0.);
            std::string mode,coordinates;serial->LoadSolverInfo("RANSWallReference",mode,"Neutral");serial->LoadSolverInfo("RANSWallFunctionCoordinates",coordinates,"Physical");
            Require(mode=="Neutral" && eta==0,"archived report requires unforced native neutral data");
            const auto signature=RANSWallPhysicsSignature(patch,serial,eta,true,us,{dx,dy},coordinates=="Chart");
            Require(signature==native.get<std::string>("signature") && signature==packet.get<std::string>("native_signature"),"archived native-report physical signature differs");
            const auto cpPath=packet.get<std::string>("native_checkpoint");
            Require(RANSWallR3::HashFile(cpPath+"/State.bin")==packet.get<std::string>("checkpoint_sha256"),"parent native checkpoint SHA256 differs");
            const auto cp=RANSWallTimeIntegration::ReadCheckpoint(cpPath,signature,1,dt,patch.Size(),patch.DifferentialSize());
            const auto steps=native.get<std::size_t>("steps");const auto &states=native.get_child("accepted_states");
            Require(!states.empty() && steps+1>=states.size(),"archived accepted-state lineage is incomplete");
            const auto first=steps+1-states.size();Require(first<=cp.integration.steps && cp.integration.steps<=steps,"parent checkpoint is outside native report lineage");
            auto parent=states.begin();std::advance(parent,cp.integration.steps-first);Vector parentState;
            for(const auto &v:parent->second)parentState.push_back(v.second.get_value<double>());
            Require(RANSWallR3::Hash(parentState)==RANSWallR3::Hash(cp.coupledState),"parent checkpoint differs from immutable native accepted-state lineage");
            const auto committed=Read(attempt,"committed");Require(RANSWallR3::Hash(committed)==RANSWallR3::Hash(Read(native,"state")),"terminal native accepted state differs from packet committed state");
            const double t=attempt.get<double>("time"),alpha=attempt.get<double>("alpha");
            Require(std::abs(t-(native.get<double>("time")+dt))<1e-14 && std::abs(alpha-1/dt)<1e-10,"archived report attempted time/alpha differs");
            Require(RANSWallR3::Hash(cp.fixedReference)==RANSWallR3::Hash(Read(attempt,"fixed_reference")),"packet fixed reference differs from parent checkpoint");
            data.fixedBoundary=patch.FixedReferenceBoundary(cp.fixedReference);std::map<int,double> external;
            for(const auto &item:attempt.get_child("fixed_external_coefficients"))
            {auto it=item.second.begin();const int row=int(it->second.get_value<double>());++it;Require(external.emplace(row,it->second.get_value<double>()).second,"duplicate archived external row");}
            Require(external==data.fixedBoundary,"archived report external data differ");data.fixedBoundary=external;
            patch.RestoreReference(cp.fixedReference,cp.fixedDefect,cp.eta,data,0.);
            auto snapshot=cp.integration;snapshot.steps=steps;snapshot.time=native.get<double>("time");
            snapshot.state.assign(committed.begin(),committed.begin()+patch.DifferentialSize());snapshot.history.clear();
            std::vector<std::string> historyHashes;
            for(const auto &slot:native.get_child("native_history"))
            {Vector values;for(const auto &v:slot.second)values.push_back(v.second.get_value<double>());snapshot.history.push_back(values);historyHashes.push_back(J::Quote(RANSWallR3::Hash(values)));}
            snapshot.times=Read(native,"native_history_times");
            Require(J::List(historyHashes)==packet.get<std::string>("native_history_hashes_json") &&
                RANSWallR3::Hash(snapshot.times)==packet.get<std::string>("native_history_times_sha256"),"archived report native histories differ");
            struct Captured {};
            RANSWallTimeIntegration integration(1,dt,[&](const Vector &h,double callTime,double lambda)->Vector
                {formedHistory=committed;std::copy(h.begin(),h.end(),formedHistory.begin());Require(std::abs(callTime-t)<1e-14 && std::abs(1/lambda-alpha)<1e-10,"report callback time differs");throw Captured{};},
                [&](const Vector &v,double callTime){auto full=committed;std::copy(v.begin(),v.end(),full.begin());return patch.ExplicitBody(full,callTime,data);});
            integration.Restore(snapshot);try{integration.Advance();}catch(const Captured&){}
            Require(RANSWallR3::Hash(formedHistory)==RANSWallR3::Hash(Read(attempt,"history")),"archived report GLM history does not reproduce packet callback history");
            const auto after=integration.Capture();Require(after.steps==snapshot.steps && after.time==snapshot.time &&
                RANSWallR3::Hash(after.state)==RANSWallR3::Hash(snapshot.state) && RANSWallR3::Hash(after.times)==RANSWallR3::Hash(snapshot.times),"report history capture changed accepted state/time");
            for(int slot=0;slot<snapshot.history.size();++slot)Require(RANSWallR3::Hash(after.history[slot])==RANSWallR3::Hash(snapshot.history[slot]),"report history capture changed native slots");
            std::vector<RANSWallCoordinate> witnesses;
            for(const auto &entry:packet.get_child("witnesses")){RANSWallCoordinate q{};int a=0;for(const auto &v:entry.second){Require(a<3,"report witness shape");q[a++]=v.second.get_value<double>();}Require(a==3,"report witness shape");witnesses.push_back(q);}
            std::vector<std::pair<std::string,std::string>> diagnostics,audits;
            for(const std::string name:{"committed","post_flow","pre_k","raw","corrected"})
            {const auto values=Read(attempt,name);Require(values.size()==patch.Size() && RANSWallR3::Hash(values)==packet.get<std::string>("stage_hashes."+name),"archived report stage hash differs: "+name);
                diagnostics.push_back({name,patch.ConstitutiveStateDiagnostic(values,t,data,witnesses)});audits.push_back({name,Maps(patch.CompleteRestrictedAudit(values,formedHistory,alpha,t,data))});}
            Require(patch.InitializationWork()==Vector({0.,0.}),"report readback performed initialization lifting");
            std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote(operation)},{"audit_revision",J::Quote(patch.AuditRevision())},
                {"parent_checkpoint_step",J::Number(cp.integration.steps)},{"last_native_accepted_step",J::Number(steps)},
                {"lineage",J::Quote("immutable terminal native report accepted states, native GLM slots and callback packet; periodic parent checkpoint is not represented as the final checkpoint")},
                {"diagnostics",J::Object(diagnostics)},{"complete_audits",J::Object(audits)},{"initialization_work",J::Numbers(patch.InitializationWork())},
                {"native_history_sha256",J::List(historyHashes)},{"governing_solves","0"},{"physical_commits","0"},{"native_history_advanced","false"}})<<'\n';
            Require(out.good(),"archived native-report diagnostic write failed");world->Finalise();return 0;
        }
        if(operation=="interface_provenance")
        {RANSWallEvolutionPatch::EvaluationScope evaluation(patch,data,0.);std::ofstream out(argv[3]);out<<patch.InterfaceProvenanceChecks()<<'\n';Require(out.good(),"interface fixture write failed");world->Finalise();return 0;}
        const bool archivedStages=operation=="archived_stages";
        Tree archived;
        if(archivedStages)
        {
            Require(packet.get_optional<std::string>("native_checkpoint").is_initialized(),"archived stages require a compatible native checkpoint");
            const auto file=packet.get<std::string>("archived_packet");
            Require(RANSWallR3::HashFile(file)==packet.get<std::string>("archived_packet_sha256"),"archived stage packet SHA256 differs");
            boost::property_tree::read_json(file,archived);
            for(const std::string name:{"committed","post_flow","pre_k","raw","corrected"})
            {
                const auto values=Read(archived,name);
                Require(values.size()==patch.Size() && RANSWallR3::Hash(values)==packet.get<std::string>("stage_hashes."+name),"archived stage vector SHA256 differs: "+name);
                Maximum(values);
            }
            Require(RANSWallR3::Hash(Read(packet,"state"))==RANSWallR3::Hash(Read(archived,"committed")),"archived diagnostic state placeholder differs from committed checkpoint");
            for(const std::string name:{"committed","history","fixed_reference"})
                Require(RANSWallR3::Hash(Read(packet,name))==RANSWallR3::Hash(Read(archived,name)),"archived stage context differs: "+name);
            Require(RANSWallR3::Hash(Read(packet,"history"))==packet.get<std::string>("expected_history_sha256"),"archived callback history SHA256 differs");
            Require(packet.get<double>("time")==archived.get<double>("time") && packet.get<double>("alpha")==archived.get<double>("alpha"),"archived attempt time/alpha differs");
            Require(RANSWallR3::HashFile(archived.get<std::string>("matrix_file"))==archived.get<std::string>("matrix_sha256"),"archived frozen transport matrix SHA256 differs");
        }
        const bool strictSpatialReadback=packet.get_optional<std::string>("native_checkpoint") &&
            ((patch.AcceptedAWF() && (operation=="accepted_awf" || operation=="height_comparison" || operation=="selected_original_operator")) || operation=="evaluation_reuse_checks" || (patch.WJ() && (operation=="wj_checks" || operation=="wj_measure")));
        if(packet.get_optional<std::string>("native_checkpoint"))
        {
            Require((operation=="k_bound" || operation=="branch_coverage" || archivedStages || strictSpatialReadback) && !packet.get<bool>("projection_only",false),
                "native replay requires the ordinary production entry guard");
            Require(!(strictSpatialReadback || archivedStages) || (!packet.get_optional<std::string>("transfer_source_xml") &&
                !packet.get_optional<bool>("require_history_h0") && !packet.get<bool>("form_native_callback_history",false)),
                "strict native spatial readback forbids transfer, history waivers and callback construction");
            const auto checkpointPath=packet.get<std::string>("native_checkpoint");
            Require(RANSWallR3::HashFile(checkpointPath+"/State.bin")==packet.get<std::string>("checkpoint_sha256"),"archived checkpoint SHA256 differs");
            double dt,eta,dx,dy;serial->LoadParameter("TimeStep",dt);serial->LoadParameter("RANSWallEta",eta,0.);
            serial->LoadParameter("RANSWallDirectionX",dx,1.);serial->LoadParameter("RANSWallDirectionY",dy,0.);
            std::string mode,coordinates;serial->LoadSolverInfo("RANSWallReference",mode,"Neutral");
            serial->LoadSolverInfo("RANSWallFunctionCoordinates",coordinates,"Physical");
            Require(mode=="Neutral" && eta==0,"native replay requires the archived unforced neutral eta=0 configuration");
            const auto signature=RANSWallPhysicsSignature(patch,serial,eta,true,us,{dx,dy},coordinates=="Chart");
            Require(signature==packet.get<std::string>("native_signature"),"native replay physical signature differs from archived report");
            const auto checkpoint=RANSWallTimeIntegration::ReadCheckpoint(checkpointPath,signature,1,dt,patch.Size(),patch.DifferentialSize());
            bool tailIdentityRejected=false;
            if(packet.get<bool>("test_tail_identity_only",false))
            {
                // Pure reader-contract test. No altered profile, simulation,
                // reference or history is constructed; all real cases keep z_m=H.
                double tail;serial->LoadParameter("RANSWallTailStart",tail);
                double changed=tail*1.01;serial->SetParameter("RANSWallTailStart",changed);
                const auto incompatible=RANSWallPhysicsSignature(patch,serial,eta,true,us,{dx,dy},coordinates=="Chart");
                serial->SetParameter("RANSWallTailStart",tail);
                Require(RANSWallPhysicsSignature(patch,serial,eta,true,us,{dx,dy},coordinates=="Chart")==signature,
                    "reader identity test did not restore session parameter");
                try{RANSWallTimeIntegration::ReadCheckpoint(checkpointPath,incompatible,1,dt,patch.Size(),patch.DifferentialSize());}
                catch(const std::exception &e){tailIdentityRejected=std::string(e.what()).find("signature mismatch")!=std::string::npos;}
                Require(tailIdentityRejected,"tail parameter did not protect checkpoint compatibility");
            }
            reference=checkpoint.fixedReference;
            if(archivedStages)
            {
                Require(RANSWallR3::Hash(reference)==RANSWallR3::Hash(Read(archived,"fixed_reference")),"archived fixed reference differs from checkpoint");
                std::map<int,double> saved;
                for(const auto &item:archived.get_child("fixed_external_coefficients"))
                {auto it=item.second.begin();const int row=int(it->second.get_value<double>());++it;Require(!saved.count(row),"duplicate archived external row");saved[row]=it->second.get_value<double>();}
                Require(saved==patch.FixedReferenceBoundary(reference),"archived external coefficients differ from saved native reference");
            }
            const auto keys=patch.FixedReferenceBoundary(reference);
            if(strictSpatialReadback || archivedStages)data.fixedBoundary=keys;
            if(packet.get_child_optional("fixed_external_coefficients"))
            {
                std::map<int,double> archived;
                for(const auto &item:packet.get_child("fixed_external_coefficients"))
                {auto it=item.second.begin();const int row=int(it->second.get_value<double>());++it;Require(!archived.count(row),"duplicate archived essential row");archived[row]=it->second.get_value<double>();}
                if(strictSpatialReadback || archivedStages)Require(archived==keys,"saved-reference boundary data differ from archived external data");
                else data.fixedBoundary=archived;
            }
            Require(keys.size()==data.fixedBoundary.size(),"archived essential row inventory differs");
            double boundaryDefect=0.;const auto matching=patch.DissipationTrace(reference,0,data);
            for(auto [row,value]:data.fixedBoundary)
            {
                Require(keys.count(row),"archived essential row absent from native geometry");
                const bool trace=row>=patch.offset[patch.Epsilon()] && row<patch.offset[patch.InnerPressure()] && matching[0].count(row-patch.offset[patch.Epsilon()]);
                if(!trace)boundaryDefect=std::max(boundaryDefect,std::abs(checkpoint.coupledState[row]-value));
            }
            Require(boundaryDefect==0.,"archived accepted native state disagrees with external data");
            patch.RestoreReference(reference,checkpoint.fixedDefect,checkpoint.eta,data,0.);
            std::vector<std::string> nativeHistoryHashes;
            if(strictSpatialReadback)
            {
                // This is an alpha=0 observation of exactly the retained
                // native checkpoint, not a new callback or a restart.  Its
                // complete GLM history is independently checked byte for byte.
                Require(packet.get<double>("alpha")==0. && packet.get<double>("time")==checkpoint.integration.time,
                    "strict native spatial readback requires alpha=0 and exact checkpoint time");
                Require(RANSWallR3::Hash(Read(packet,"state"))==RANSWallR3::Hash(checkpoint.coupledState),
                    "strict native spatial state differs from checkpoint");
                Require(RANSWallR3::Hash(Read(packet,"history"))==RANSWallR3::Hash(checkpoint.coupledState),
                    "alpha=0 API storage placeholder differs from retained coupled state");
            }
            if(strictSpatialReadback || archivedStages)
            {
                Require(RANSWallR3::Hash(Read(packet,"native_differential_state"))==RANSWallR3::Hash(checkpoint.integration.state),
                    "native differential state differs from checkpoint");
                const auto &savedHistory=packet.get_child("native_history");
                Require(savedHistory.size()==checkpoint.integration.history.size(),"native history slot count differs from checkpoint");
                int slot=0;for(const auto &item:savedHistory)
                {
                    Vector values;for(const auto &value:item.second)values.push_back(value.second.get_value<double>());
                    const auto hash=RANSWallR3::Hash(values);
                    Require(hash==RANSWallR3::Hash(checkpoint.integration.history[slot++]),"native history slot differs from checkpoint");
                    nativeHistoryHashes.push_back(J::Quote(hash));
                }
                Require(RANSWallR3::Hash(Read(packet,"native_history_times"))==RANSWallR3::Hash(checkpoint.integration.times),
                    "native history time labels differ from checkpoint");
            }
            if(!strictSpatialReadback)
            {
                Require(RANSWallR3::Hash(Read(packet,"committed"))==RANSWallR3::Hash(checkpoint.coupledState),"raw packet is not based on archived accepted state");
                Require(std::abs(packet.get<double>("time")-(checkpoint.integration.time+dt))<1e-14,"native replay physical time differs");
                Require(std::abs(packet.get<double>("alpha")-1/dt)<1e-10,"native replay alpha differs");
            }
            if(archivedStages || packet.get<bool>("form_native_callback_history",false))
            {
                Require(archivedStages || operation=="branch_coverage","history-only replay is a forensic operation");
                struct HistoryCaptured {};
                RANSWallTimeIntegration integration(1,dt,
                    [&](const Vector &h,double t,double lambda)->Vector
                    {
                        Require(std::abs(t-packet.get<double>("time"))<1e-14 && std::abs(1/lambda-packet.get<double>("alpha"))<1e-10,
                            "formed native callback has different time/alpha");
                        formedHistory=checkpoint.coupledState;std::copy(h.begin(),h.end(),formedHistory.begin());
                        throw HistoryCaptured{}; // No governing solve or physical commit.
                    },
                    [&](const Vector &v,double t)
                    {auto full=checkpoint.coupledState;std::copy(v.begin(),v.end(),full.begin());return patch.ExplicitBody(full,t,data);});
                integration.Restore(checkpoint.integration);
                try{integration.Advance();}catch(const HistoryCaptured &){}
                Require(!formedHistory.empty() && RANSWallR3::Hash(formedHistory)==packet.get<std::string>("expected_history_sha256"),
                    "native archived GLM slots do not reproduce recorded failed callback history");
                const auto restored=integration.Capture();
                Require(restored.time==checkpoint.integration.time && restored.steps==checkpoint.integration.steps &&
                    restored.history==checkpoint.integration.history,"forensic history formation did not roll back");
                Require(RANSWallR3::Hash(restored.state)==RANSWallR3::Hash(checkpoint.integration.state) &&
                    RANSWallR3::Hash(restored.times)==RANSWallR3::Hash(checkpoint.integration.times),
                    "forensic callback changed native state/time-label bits");
                for(int slot=0;slot<restored.history.size();++slot)
                    Require(RANSWallR3::Hash(restored.history[slot])==RANSWallR3::Hash(checkpoint.integration.history[slot]),
                        "forensic callback changed native history bits");
            }
            Require(patch.InitializationWork()==Vector({0.,0.}),"native immutable replay called Lift");
            replayContext=J::Object({{"physical_signature_matches","true"},{"checkpoint_sha256",J::Quote(packet.get<std::string>("checkpoint_sha256"))},
                {"accepted_steps",J::Number(checkpoint.integration.steps)},{"accepted_time",J::Number(checkpoint.integration.time)},
                {"fixed_reference_sha256",J::Quote(RANSWallR3::Hash(reference))},{"external_boundary_defect",J::Number(boundaryDefect)},
                {"tail_identity_only_rejected",J::Boolean(tailIdentityRejected)},
                {"initialization_lift_calls_iterations",J::Numbers(patch.InitializationWork())},{"new_initial_state_constructed","false"},
                {"strict_spatial_readback",J::Boolean(strictSpatialReadback)},
                {"native_history_sha256",J::List(nativeHistoryHashes)},
                {"native_history_times_sha256",(strictSpatialReadback || archivedStages)?J::Quote(RANSWallR3::Hash(checkpoint.integration.times)):"null"},
                {"native_history_advanced","false"},{"governing_solves","0"},
                {"formed_callback_history_sha256",formedHistory.empty()?"null":J::Quote(RANSWallR3::Hash(formedHistory))},
                {"history_policy",J::Quote(strictSpatialReadback?"complete native GLM state, slots and time labels verified; alpha=0 API storage placeholder is retained coupled state":formedHistory.empty()?"exact archived full implicit-callback history":
                    "native archived GLM slots; capture callback then abort before governing solves; recorded SHA256 matched; local snapshot restored")}});
        }
        else
        {
            const auto rawReference=patch.NeutralReference(us,{1.,0.});data.fixedBoundary=patch.FixedReferenceBoundary(rawReference);
            reference=patch.Lift(rawReference,0,data);patch.SetReference(reference,data,0,0);
        }
        if(operation=="dense_capacity_preflight")
        {
            std::ofstream out(argv[3]);out<<patch.DenseEvolutionPreflight(reference,0.,data,
                packet.get<bool>("condensed",false),packet.get<bool>("enforce",true))<<'\n';
            Require(out.good(),"capacity preflight write failed");world->Finalise();return 0;
        }
        auto history=formedHistory.empty()?Read(packet,"history"):formedHistory;
        auto state=Read(packet,"state");std::string transfer="null";
        if(packet.get_optional<std::string>("transfer_source_xml"))
        {
            const auto sourceXml=packet.get<std::string>("transfer_source_xml");
            std::vector<std::string> sourceNames{sourceXml};
            auto sourceSession=LibUtilities::SessionReader::CreateInstance(1,args,sourceNames,LibUtilities::CommSerial::create(1,args));sourceSession->InitSession();
            sourceSession->SetSolverInfo("GlobalSysSoln","DirectFull");const auto sourceGraph=SpatialDomains::MeshGraphIO::Read(sourceSession);
            RANSWallEvolutionPatch source(sourceSession,sourceGraph,world->GetComm(),p,LoadRANSWallNormalMap(sourceSession,2,p.geometry.height));
            Require(state.size()==source.Size(),"transfer source shape differs");const auto original=state;state.assign(patch.Size(),0.);Vector error;
            for(int b=0;b<patch.Gauge();++b)
            {
                const auto values=patch.fields[b]->Project([&](const RANSWallCoordinate &x)
                    {return source.Point(b,x,b==patch.Epsilon() && std::abs(x[2]-p.geometry.height)<1e-12?1:0).Value(original,source.offset[b]);});
                std::copy(values.begin(),values.end(),state.begin()+patch.offset[b]);double maximum=0.;
                for(const auto &q:patch.Quadrature(b))maximum=std::max(maximum,std::abs(q.Value(state,patch.offset[b])-source.Point(b,q.coordinate).Value(original,source.offset[b])));
                error.push_back(maximum);
            }
            state[patch.offset[patch.Gauge()]]=original[source.offset[source.Gauge()]];history=state;
            // A representation-only local test keeps the *represented* old
            // essential traces, rather than freshly projecting analytic data.
            for(auto &[id,value]:data.fixedBoundary)value=state[id];
            Require(Maximum(error)<1e-11,"local mesh transfer changed the physical polynomial");
            transfer=J::Object({{"source_xml",J::Quote(sourceXml)},{"source_state_sha256",J::Quote(RANSWallR3::Hash(original))},
                {"target_state_sha256",J::Quote(RANSWallR3::Hash(state))},{"physical_block_errors",J::Numbers(error)},
                {"essential_policy",J::Quote("exact restriction of original represented boundary traces")},{"native_history_advanced","false"}});
        }
        Require(state.size()==patch.Size() && history.size()==patch.Size(),"packet coefficient size differs");
        const bool historyMatch=std::memcmp(history.data(),reference.data(),history.size()*sizeof(double))==0;
        Require(strictSpatialReadback || archivedStages || !packet.get<bool>("require_history_h0",true) || historyMatch,"packet history differs from rebuilt native h0");
        const double time=packet.get<double>("time"),alpha=packet.get<double>("alpha");
        if(operation=="wj_measure")
        {
            Require(strictSpatialReadback && patch.WJ(),"WJ measurements require strict alpha=0 native checkpoint readback");
            const auto reportPath=packet.get<std::string>("native_report");
            const auto reportHash=RANSWallR3::HashFile(reportPath);
            Require(reportHash==packet.get<std::string>("native_report_sha256"),"WJ accepted report SHA256 differs");
            Tree report;boost::property_tree::read_json(reportPath,report);
            Require(report.get<std::string>("signature")==packet.get<std::string>("native_signature"),
                "WJ accepted report signature differs from checkpoint");
            const auto cp=RANSWallTimeIntegration::ReadCheckpoint(packet.get<std::string>("native_checkpoint"),
                packet.get<std::string>("native_signature"),1,serial->GetParameter("TimeStep"),patch.Size(),patch.DifferentialSize());
            const auto reportSteps=report.get<std::size_t>("steps");
            const auto &accepted=report.get_child("accepted_states"),&packets=report.get_child("wj_audit_packets");
            Require(!accepted.empty() && accepted.size()==packets.size() && reportSteps+1>=accepted.size(),"WJ accepted-stage lineage shape differs");
            const auto first=reportSteps+1-accepted.size();
            Require(first<=cp.integration.steps && cp.integration.steps<=reportSteps,"WJ checkpoint step outside accepted-stage lineage");
            auto ai=accepted.begin();auto pi=packets.begin();std::advance(ai,cp.integration.steps-first);std::advance(pi,cp.integration.steps-first);
            Vector acceptedState;for(const auto &v:ai->second)acceptedState.push_back(v.second.get_value<double>());
            Require(RANSWallR3::Hash(acceptedState)==RANSWallR3::Hash(state) && pi->second.get<double>("time")==time,
                "WJ accepted-stage state/time differs from checkpoint");
            Vector stageHistory;double stageAlpha=0;
            if(cp.integration.steps>0)
            {
                const auto &binding=pi->second;
                stageHistory=Read(binding,"callback_history");stageAlpha=binding.get<double>("alpha_mass");
                Require(stageHistory.size()==state.size() && stageAlpha==1/serial->GetParameter("TimeStep"),
                    "WJ completed-stage context differs from accepted checkpoint");
                if(cp.integration.steps==reportSteps)
                    Require(report.get<double>("last_stage_time")==time && report.get<double>("last_stage_alpha")==stageAlpha &&
                        RANSWallR3::Hash(Read(report,"last_stage_history"))==RANSWallR3::Hash(stageHistory),"WJ last callback differs from immutable accepted-stage packet");
                Require(binding.get<std::string>("state_sha256")==RANSWallR3::Hash(state) &&
                    binding.get<std::string>("stage_history_sha256")==RANSWallR3::Hash(stageHistory) &&
                    binding.get<double>("time")==time && binding.get<double>("alpha_mass")==stageAlpha,
                    "WJ completed-stage audit binding differs");
            }
            std::vector<std::array<double,8>> observations;
            for(const auto &entry:packet.get_child("observations"))
            {std::array<double,8> q{};q[5]=-1;int i=0;for(const auto &v:entry.second){Require(i<8,"WJ observation shape");q[i++]=v.second.get_value<double>();}Require(i>=5 && i<=8,"WJ observation shape");observations.push_back(q);}
            const auto checkpointHash=RANSWallR3::HashFile(packet.get<std::string>("native_checkpoint")+"/State.bin");
            std::ofstream out(argv[3]);out<<J::Object({{"replay",replayContext},{"native_report_sha256",J::Quote(reportHash)},
                {"measurements",patch.WJMeasurements(state,stageHistory,stageAlpha,time,data,observations,packet.get<bool>("export_wj5_diagnostics",false))},
                {"assembly_checks",packet.get<bool>("export_assembly_checks",false)?patch.WJNativeChecks(state,data,time,stageHistory,stageAlpha):"null"}})<<'\n';
            Require(out.good() && checkpointHash==RANSWallR3::HashFile(packet.get<std::string>("native_checkpoint")+"/State.bin") &&
                reportHash==RANSWallR3::HashFile(reportPath),"WJ measurement changed saved checkpoint/report bytes");
            world->Finalise();return 0;
        }
        if(operation=="wj_checks")
        {
            Require(strictSpatialReadback && patch.WJ(),"WJ checks require immutable native checkpoint readback");
            std::ofstream out(argv[3]);out<<J::Object({{"replay",replayContext},{"checks",patch.WJNativeChecks(state,data,time)},
                {"source_balance",patch.SourceBalanceDiagnostic(state,time,data,{.5*(patch.lower[0]+patch.upper[0])},{0.,.1,.5,1.})}})<<'\n';
            Require(out.good(),"WJ native check write failed");world->Finalise();return 0;
        }
        if(operation=="evaluation_reuse_checks")
        {std::ofstream out(argv[3]);out<<patch.EvaluationReuseChecks(state,history,time,data)<<'\n';Require(out.good(),"reuse verification write failed");world->Finalise();return 0;}
        RANSWallEvolutionPatch::EvaluationScope evaluation(patch,data,time);
        if(archivedStages || operation=="selected_original_operator")
        {
            Require(patch.InitializationWork()==Vector({0.,0.}),"immutable diagnostic called initialization Lift");
            std::vector<RANSWallCoordinate> witnesses;
            for(const auto &q:packet.get_child("witnesses")){RANSWallCoordinate point{};int axis=0;for(const auto &x:q.second){Require(axis<3,"witness shape differs");point[axis++]=x.second.get_value<double>();}Require(axis==3,"witness shape differs");witnesses.push_back(point);}
            std::vector<std::pair<std::string,std::string>> diagnostics,completeAudits;
            if(archivedStages)for(const std::string name:{"committed","post_flow","pre_k","raw","corrected"})
            {
                const auto values=Read(archived,name);const auto hash=RANSWallR3::Hash(values);
                diagnostics.push_back({name,patch.ConstitutiveStateDiagnostic(values,time,data,witnesses)});
                completeAudits.push_back({name,Maps(patch.CompleteRestrictedAudit(values,Read(archived,"history"),alpha,time,data))});
                Require(hash==RANSWallR3::Hash(values),"archived stage coefficients changed");
            }
            else
            {
                diagnostics.push_back({"selected_state",patch.ConstitutiveStateDiagnostic(state,time,data,witnesses)});
                completeAudits.push_back({"selected_state",Maps(patch.CompleteRestrictedAudit(state,history,alpha,time,data))});
                diagnostics.push_back({"operators",patch.CompareOuterOperators(state,time,data)});
                const auto panels=packet.get_child_optional("reporting_breakpoints")?Read(packet,"reporting_breakpoints"):Vector{};
                diagnostics.push_back({"physical",patch.HeightComparison(state,state,time,data,panels)});
                if(packet.get<bool>("export_k_mass",false))
                {
                    const auto constraints=patch.Constraints(state,time,data);std::vector<int> free;
                    for(int i=0;i<patch.size[patch.K()];++i)
                        if(!constraints.count(patch.offset[patch.K()]+i))free.push_back(patch.offset[patch.K()]+i);
                    diagnostics.push_back({"native_free_k_mass",J::Object({{"rows",J::Numbers(free)},
                        {"matrix",J::Rows(patch.PicardKMass(free))},{"governing_solves","0"}})});
                }
            }
            std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote(operation)},
                {"native_replay_context",replayContext},{"audit_revision",J::Quote(patch.AuditRevision())},
                {"diagnostics",J::Object(diagnostics)},{"complete_audits",J::Object(completeAudits)},
                {"physical_time",J::Number(time)},{"governing_solves","0"},{"physical_commits","0"},
                {"native_history_advanced","false"},{"initialization_work",J::Numbers(patch.InitializationWork())}})<<'\n';
            Require(out.good(),"immutable diagnostic write failed");world->Finalise();return 0;
        }
        if(packet.get<std::string>("operation","")=="branch_coverage")
        {
            std::vector<RANSWallCoordinate> anchors;
            const auto xs=packet.get_child_optional("stations")?Read(packet,"stations"):Vector{-3.5,-2.625,-1.,1.6};
            for(double x:xs)anchors.push_back({x,0,0});
            std::vector<std::pair<std::string,std::string>> report{{"state",patch.BranchCoverageDiagnostic(state,time,data,anchors)},
                {"native_replay_context",replayContext}};
            if(packet.get<bool>("compare_prospective_inventory_only",false))
            {
                // A read-only counterfactual inventory check, not a checkpoint
                // reader or migration. The archived patch above has already
                // checked its own signature, boundary values and native slots.
                // This second patch only evaluates the same saved coefficients:
                // no Lift, normalization, transport solve or history advance.
                Require(packet.get_optional<std::string>("native_checkpoint").is_initialized(),
                    "inventory counterfactual requires the verified native archive");
                std::string saved;serial->LoadSolverInfo("RANSWallPimpleKInventory",saved,"LegacyV1");
                serial->SetSolverInfo("RANSWallPimpleKInventory","ConstitutiveV3");
                RANSWallEvolutionPatch prospective(serial,graph,world->GetComm(),p,LoadRANSWallNormalMap(serial,2,p.geometry.height));
                serial->SetSolverInfo("RANSWallPimpleKInventory",saved);
                const auto checked=prospective.ScalarDomain(state,time,data,true,true);
                Require(prospective.InitializationWork()==Vector({0.,0.}),"inventory counterfactual called Lift");
                report.push_back({"prospective_inventory_only",RANSWallEvolutionPatch::ScalarDomainJSON(checked)});
                report.push_back({"prospective_inventory_scope",J::Quote(
                    "read-only V3 constraint enumeration on the original archived coefficients/data; no compatible restart claimed")});
            }
            if(packet.get_child_optional("compare_state"))
            {
                const auto compare=Read(packet,"compare_state");Require(compare.size()==state.size(),"comparison shape differs");
                double difference=0;for(int i=patch.offset[patch.K()];i<patch.offset[patch.Epsilon()];++i)
                    difference=std::max(difference,std::abs(state[i]-compare[i]));
                report.push_back({"compare",patch.BranchCoverageDiagnostic(compare,time,data,anchors)});
                report.push_back({"k_coefficient_difference",J::Number(difference)});
            }
            std::ofstream out(argv[3]);out<<J::Object(report)<<'\n';Require(out.good(),"branch readback write failed");world->Finalise();return 0;
        }
        if(packet.get<std::string>("operation","")=="k_inventory")
        {
            std::ofstream out(argv[3]);out<<patch.ScalarDomainJSON(patch.ScalarDomain(state,time,data,true,true),true)<<'\n';
            Require(out.good(),"inventory write failed");world->Finalise();return 0;
        }
        if(packet.get<std::string>("operation","")=="coverage")
        {
            std::ofstream out(argv[3]);out<<patch.CheckConstitutiveCoverage(state,time,data,{1.23456789,0,0})<<'\n';
            Require(out.good(),"coverage check write failed");world->Finalise();return 0;
        }
        if(packet.get<std::string>("operation","")=="audit_witness")
        {
            RANSWallEvolutionPatch::Diagnostic audit;
            std::string failure;try{audit=patch.CompleteRestrictedAudit(state,history,alpha,time,data);}
            catch(const std::exception &e){failure=e.what();}
            std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote("audit_witness")},
                {"failure",J::Quote(failure)},{"audit",Maps(audit)},{"context",patch.ConstitutiveFailureContext().empty()?"null":patch.ConstitutiveFailureContext()},
                {"governing_solves","0"}})<<'\n';Require(out.good(),"audit witness write failed");world->Finalise();return 0;
        }
        if(packet.get<std::string>("operation","")=="source_balance")
        {
            const auto stations=packet.get_child_optional("stations")?Read(packet,"stations"):Vector{-3.5,-2.625,-1.,1.6};
            const auto heights=packet.get_child_optional("heights")?Read(packet,"heights"):Vector{0.,.00001,.00005,.0001,.0002,.0005,.001,.0025,.005,.01,.05,.1,1.,2.};
            std::ofstream out(argv[3]);out<<patch.SourceBalanceDiagnostic(state,time,data,stations,heights)<<'\n';
            Require(out.good(),"source balance write failed");world->Finalise();return 0;
        }
        if(packet.get<std::string>("operation","")=="sampling_checks")
        {
            if(patch.AcceptedAWF())
            {
                const double H=p.geometry.height;auto polynomial=state;
                const auto coefficients=patch.fields[patch.K()]->Project([&](const RANSWallCoordinate &x)
                    {return 1+.02*x[0]+.1*x[2]+.03*x[2]*x[2]/(H*H);});
                std::copy(coefficients.begin(),coefficients.end(),polynomial.begin()+patch.offset[patch.K()]);
                double endpointError=0,rowError=0,profileError=0,traceError=0;
                for(double x:{-3.5,-1.,1.6})
                {
                    const double expected=1+.02*x+.1*H+.03;const RANSWallCoordinate anchor{x,0,0};
                    const auto row=patch.SampleKRow(anchor);const double got=patch.SampleK(polynomial,anchor);
                    endpointError=std::max(endpointError,std::abs(got-expected));
                    long double direct=0;for(int i=0;i<row.ids.size();++i)direct+=(long double)row.phi[i]*polynomial[patch.offset[patch.K()]+row.ids[i]];
                    rowError=std::max(rowError,std::abs(got-double(direct)));
                    profileError=std::max(profileError,std::abs(patch.StateProfile(polynomial,anchor,0,data)->K()-expected));
                    traceError=std::max(traceError,std::abs(patch.DependentEpsilon(polynomial,{x,0,H},0,data)-
                        std::pow(expected/p.model.Q(),1.5)/(p.model.Kappa()*H)));
                }
                Require(std::max({endpointError,rowError,profileError,traceError})<1e-11,"accepted endpoint functional/dependent consumer mismatch");
                std::ofstream out(argv[3]);out<<J::Object({{"endpoint_sampling_error",J::Number(endpointError)},
                    {"exact_global_row_error",J::Number(rowError)},{"nonuniform_profile_K_error",J::Number(profileError)},
                    {"dependent_epsilon_H_error",J::Number(traceError)},{"chi","null"},{"governing_solves","0"}})<<'\n';
                Require(out.good(),"accepted sampler check write failed");world->Finalise();return 0;
            }
            const double H=p.geometry.height;auto polynomial=reference;
            const auto values=patch.fields[patch.K()]->Project([&](const RANSWallCoordinate &x){return 1+.02*x[0]+.1*x[2]+.03*x[2]*x[2]/(H*H)+.005*std::pow(x[2]/H,3);});
            std::copy(values.begin(),values.end(),polynomial.begin()+patch.offset[patch.K()]);
            double error=0,rowError=0,constantError=0,tangentialError=0,consumerError=0,traceConsumerError=0;
            std::vector<std::string> samples;
            for(double x:{-3.5,-1.,1.6})
            {
                const auto row=patch.SampleKRow({x,0,0});const double got=patch.SampleK(polynomial,{x,0,0});
                constantError=std::max(constantError,std::abs(patch.SampleK(reference,{x,0,0})-p.model.Q()*us*us));
                tangentialError=std::max(tangentialError,std::abs(row.Gradient(polynomial,patch.offset[patch.K()])[0]-.02));
                const double expected=patch.ColumnAverageK()?1+.02*x+.05*H+.01+.00125:
                    1+.02*x+.1*p.geometry.sampleHeight+.03*std::pow(p.geometry.sampleHeight/H,2)+.005*std::pow(p.geometry.sampleHeight/H,3);
                // Nonuniform Upstream fixture: independently known K feeds
                // the transition's intrinsic input AND dependent epsilon.
                // Uniform chi=1 alone cannot exercise these state inputs.
                auto anchor=RANSWallCoordinate{x,0,0};const auto current=patch.StateProfile(polynomial,anchor,0,data);
                auto options=current->Options();options.requestedWeight=0;
                const EquilibriumAWFProfiles baseline(current->Geometry(),expected,p.model,{},options);
                const auto input=patch.TransitionInput(polynomial,anchor,0,data);
                consumerError=std::max({consumerError,std::abs(current->K()-expected),
                    std::abs(input.reynolds-baseline.Reference().reynolds)/baseline.Reference().reynolds,
                    std::abs(input.budgetMargin-baseline.Reference().margin)});
                options.requestedWeight=current->Reference().weight;
                const EquilibriumAWFProfiles selected(current->Geometry(),expected,p.model,{},options);
                anchor[2]=H;
                traceConsumerError=std::max(traceConsumerError,std::abs(patch.DependentEpsilon(polynomial,anchor,0,data)-
                    patch.WallDissipation(polynomial,anchor,selected)));
                long double linear=0;for(int j=0;j<row.ids.size();++j)linear+=(long double)row.phi[j]*polynomial[patch.offset[patch.K()]+row.ids[j]];
                error=std::max(error,std::abs(got-expected));rowError=std::max(rowError,std::abs(got-double(linear)));
                samples.push_back(J::Object({{"x",J::Number(x)},{"value",J::Number(got)},{"expected",J::Number(expected)},
                    {"row_ids",J::Numbers(row.ids)},{"row_phi",J::Numbers(row.phi)}}));
            }
            Require(error<1e-11 && rowError<1e-12 && constantError<1e-12,"representative K constant/polynomial/row verification failed");
            Require(tangentialError<1e-12,"column functional lost tangential variation");
            Require(consumerError<1e-11 && traceConsumerError<1e-10,"transition/dependent epsilon lost selected K functional");
            Vector zero(patch.Size());Tree bands;std::istringstream bandJson(patch.PhysicalKBands(polynomial,0,data,&zero,true));
            boost::property_tree::read_json(bandJson,bands);double bandError=0,bandRowError=0;
            for(const auto &entry:bands.get_child("regions"))
            {
                const auto &region=entry.second;const double z=region.get<double>("upper");
                const auto bounds=Read(region,"x_bounds");const double length=bounds[1]-bounds[0];
                const double expected=length*((1+.01*(bounds[0]+bounds[1]))*z+.05*z*z+.01*std::pow(z,3)/(H*H)+.00125*std::pow(z,4)/(H*H*H));
                const auto row=Read(region,"coefficient_integral_row");long double value=0;
                for(int i=0;i<row.size();++i)value+=(long double)row[i]*polynomial[patch.offset[patch.K()]+i];
                bandError=std::max(bandError,std::abs(region.get<double>("k_integral")-expected));
                bandRowError=std::max(bandRowError,std::abs(double(value)-expected));
            }
            Require(bandError<1e-10 && bandRowError<1e-10,"physical band integral/row missed a cut element");
            const auto domain=patch.ScalarDomain(polynomial,0,data,true,true);int means=0;
            Require(domain.admissible,"polynomial sampler fixture left scalar domain");
            if(patch.ColumnAverageK())for(const auto &o:domain.observations)if(o.inventory=="physical_sampled_K")
            {Require(o.functional && o.functional->columnAverage,"mean observation lost its row");++means;}
            Require(!patch.ColumnAverageK() || means>0,"mean observations were not exercised");
            double meanProjectionError=0,meanStationarity=0;
            if(patch.ColumnAverageK())
            {
                // Exercise the very observation row carried to the production
                // separation wrapper, with essential coefficients eliminated.
                const auto found=std::find_if(domain.observations.begin(),domain.observations.end(),
                    [](const auto &o){return o.inventory=="physical_sampled_K" && std::abs(o.coordinate[0]+1)<.3;});
                Require(found!=domain.observations.end() && found->functional,"missing mean projection test row");
                const auto &row=*found->functional;std::vector<int> free,index(patch.size[patch.K()],-1);
                for(int i=0;i<index.size();++i)if(!data.fixedBoundary.count(patch.offset[patch.K()]+i))
                {index[i]=free.size();free.push_back(patch.offset[patch.K()]+i);}
                Vector evaluation(free.size());
                for(int i=0;i<row.ids.size();++i)if(index[row.ids[i]]>=0)evaluation[index[row.ids[i]]]+=row.phi[i];
                RANSWallKBound::Projection projection(patch.PicardKMass(free));
                Require(projection.Add(evaluation,.125,0),"mean-row consistent-mass projection failed");
                auto corrected=polynomial;for(int i=0;i<free.size();++i)corrected[free[i]]+=projection.delta[i];
                meanProjectionError=std::abs(row.Value(corrected,patch.offset[patch.K()])-row.Value(polynomial,patch.offset[patch.K()])-.125);
                meanStationarity=projection.Stationarity();
                Require(meanProjectionError<1e-12 && meanStationarity<1e-12,"mean constraint became a point constraint");
                for(int i=0;i<corrected.size();++i)if(i<patch.offset[patch.K()] || i>=patch.offset[patch.Epsilon()] || data.fixedBoundary.count(i))
                    Require(corrected[i]==polynomial[i],"mean projection changed protected coefficients");
            }
            auto negative=reference;const auto k=patch.fields[patch.K()]->Project([&](const RANSWallCoordinate &x){return -.1+2*x[2]/H;});
            std::copy(k.begin(),k.end(),negative.begin()+patch.offset[patch.K()]);
            const double mean=patch.SampleK(negative,{-1,0,0});const auto bad=patch.ScalarDomain(negative,0,data,true,true);
            Require(mean>0 && !bad.admissible,"positive representative K waived negative local k");
            std::ofstream out(argv[3]);out<<J::Object({{"polynomial_error",J::Number(error)},{"coefficient_row_error",J::Number(rowError)},
                {"constant_error",J::Number(constantError)},
                {"nonuniform_transition_consumer_error",J::Number(consumerError)},
                {"nonuniform_dependent_epsilon_error",J::Number(traceConsumerError)},
                {"tangential_derivative_error",J::Number(tangentialError)},
                {"physical_band_integral_error",J::Number(bandError)},{"physical_band_row_error",J::Number(bandRowError)},
                {"mean_observations_with_exact_rows",J::Number(means)},{"samples",J::List(samples)},
                {"mean_projection_constraint_error",J::Number(meanProjectionError)},{"mean_projection_stationarity",J::Number(meanStationarity)},
                {"positive_K_negative_local_rejected","true"},{"local_failure",J::Quote(bad.failure)}})<<'\n';
            Require(out.good(),"sampling check write failed");world->Finalise();return 0;
        }
        if(packet.get<std::string>("operation","")=="height_comparison" || packet.get<std::string>("operation","")=="accepted_awf")
        {
            const double time=packet.get<double>("time");
            Require(patch.ScalarDomain(state,time,data).admissible,"height replay outside hard scalar domain");
            const auto reporting=packet.get_child_optional("reporting_breakpoints")?Read(packet,"reporting_breakpoints"):Vector{};
            const auto before=RANSWallR3::Hash(state);auto result=patch.HeightComparison(state,history,time,data,reporting);
            Require(before==RANSWallR3::Hash(state),"spatial readback changed retained coefficients");
            Require(!result.empty() && result.back()=='}',"invalid height comparison JSON");result.pop_back();
            result+=",\"native_replay_context\":"+replayContext+
                ",\"reference_initialization_lift_calls_iterations\":"+J::Numbers(patch.InitializationWork())+
                ",\"readback_selected_state_sha256\":"+J::Quote(before)+
                ",\"readback_selected_state_unchanged\":true,\"readback_governing_solves\":0,\"native_history_advanced\":false}";
            std::ofstream out(argv[3]);out<<result<<'\n';
            Require(out.good(),"height replay write failed");world->Finalise();return 0;
        }
        if(packet.get<std::string>("operation","")=="contact")
        {
            Vector steps=packet.get_child_optional("h")?Read(packet,"h"):Vector{};
            std::vector<RANSWallCoordinate> witnesses;
            if(packet.get_child_optional("witnesses"))for(const auto &item:packet.get_child("witnesses"))
            {RANSWallCoordinate p{};int i=0;for(const auto &v:item.second){Require(i<3,"witness coordinate shape");p[i++]=v.second.get_value<double>();}Require(i==3,"witness coordinate shape");witnesses.push_back(p);}
            const auto &selected=packet.get<bool>("reference_state",false)?reference:state;
            const auto domain=patch.ScalarDomain(selected,time,data,true);
            std::string result;
            if(!domain.admissible)result=J::Object({{"diagnostic_domain_admissible","false"},
                {"scalar_domain",patch.ScalarDomainJSON(domain)},{"state_sha256",J::Quote(RANSWallR3::Hash(selected))},{"production_governing_solves","0"}});
            else result=patch.PimpleContactDiagnostic(selected,time,data,steps,witnesses,
                packet.get<bool>("export_mass",false),packet.get<bool>("isolate",false));
            Require(result.back()=='}',"invalid diagnostic JSON");result.pop_back();result+=",\"transfer\":"+transfer+"}";
            std::ofstream out(argv[3]);out<<result<<'\n';
            Require(out.good(),"contact diagnostic write failed");world->Finalise();return 0;
        }
        if(packet.get<std::string>("operation","")=="samples")
        {
            Require(patch.ScalarDomain(state,time,data).admissible,"sample replay outside scalar domain");
            std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote("samples")},
                {"time",J::Number(time)},{"state",J::Numbers(state)},
                {"physical_samples",Maps(patch.PhysicalSamples(state,time,data))},{"governing_solves","0"}})<<'\n';
            Require(out.good(),"sample replay write failed");world->Finalise();return 0;
        }
        if(packet.get<std::string>("operation","")=="k_budget")
        {
            std::ofstream out(argv[3]);out<<patch.PimpleKRateBudget(state,history,time,data)<<'\n';
            Require(out.good(),"budget replay write failed");world->Finalise();return 0;
        }
        if(packet.get<std::string>("operation","")=="k_bound")
        {
            Vector pre=packet.get_child_optional("pre_k")?Read(packet,"pre_k"):history;
            const auto fixture=packet.get<std::string>("fixture","");
            double nodeMinimum=1e300;std::string fixtureDefinition="null";
            if(fixture=="accepted_negative")
            {
                Require(patch.AcceptedAWF(),"accepted negative fixture needs its own formulation");
                pre=state;const RANSWallCoordinate point{0,0,0};const auto row=patch.Point(patch.K(),point);
                Vector direction(patch.size[patch.K()],0.);long double square=0.;
                for(int j=0;j<row.ids.size();++j)if(!data.fixedBoundary.count(patch.offset[patch.K()]+row.ids[j]))
                    direction[row.ids[j]]+=row.phi[j];
                for(double x:direction)square+=(long double)x*x;
                Require(square>0.,"localized accepted negative fixture has no free evaluation direction");
                const double target=-1e-5*p.model.Q()*us*us;
                const double factor=(row.Value(state,patch.offset[patch.K()])-target)/double(square);
                for(int j=0;j<direction.size();++j)state[patch.offset[patch.K()]+j]-=factor*direction[j];
                fixtureDefinition=J::Object({{"name",J::Quote(fixture)},{"coordinate",J::Numbers(point)},
                    {"definition",J::Quote("single exact physical evaluation row restricted to free global k coefficients")},
                    {"target_raw_k",J::Number(target)},{"actual_raw_k",J::Number(row.Value(state,patch.offset[patch.K()]))},
                    {"free_row_squared_norm",J::Number(double(square))}});
            }
            else if(fixture=="accepted_large_negative")
            {
                Require(patch.AcceptedAWF(),"accepted negative fixture needs its own formulation");
                pre=state;for(int i=0;i<60;++i)state=patch.AcceptedPerturbation(state,-.05);
                fixtureDefinition=J::Object({{"name",J::Quote(fixture)},{"definition",J::Quote("historical 60 successive -0.05 global accepted perturbations")}});
            }
            else if(!fixture.empty())
            {
                const double K0=p.model.Q()*us*us;auto field=patch.fields[patch.K()];
                int chosen=-1;for(int i=0;i<field->Elements().size();++i)if(field->Elements()[i].geometry==4)chosen=i;
                Require(chosen>=0,"missing geometry4 fixture");const auto el=field->Elements()[chosen];
                const auto coefficients=field->Project([&](const RANSWallCoordinate &x)
                {
                    if(fixture=="smooth")return K0*(1+.01*(x[0]-patch.lower[0])*(patch.upper[0]-x[0])*
                        x[2]*(patch.upper[2]-x[2])/std::pow(patch.upper[0]-patch.lower[0],2)/std::pow(patch.upper[2],2));
                    if(fixture=="between_nodes" && x[0]>=el.lower[0] && x[0]<=el.upper[0] && x[2]<=el.upper[2])
                    {const double z=2*(x[0]-el.lower[0])/(el.upper[0]-el.lower[0])-1;return K0*(1-(.997+.16*z)*(1-z*z)*(1-x[2]/el.upper[2]));}
                    return K0;
                });
                pre=reference;state=reference;
                if(fixture!="reference")for(int i=0;i<coefficients.size();++i)state[patch.offset[patch.K()]+i]=coefficients[i];
                for(auto [id,value]:data.fixedBoundary)if(id>=patch.offset[patch.K()] && id<patch.offset[patch.Epsilon()])state[id]=value;
                const auto exp=field->Field()->GetExp(chosen);const auto z=exp->GetBasis(0)->GetZ();
                for(int i=0;i<exp->GetNumPoints(0);++i)
                    nodeMinimum=std::min(nodeMinimum,field->ElementPoint(chosen,{el.lower[0]+(z[i]+1)*(el.upper[0]-el.lower[0])/2,0,0}).Value(state,patch.offset[patch.K()]));
            }
            const auto preDomain=patch.ScalarDomain(pre,time,data,false,true);
            if(packet.get_optional<std::string>("native_checkpoint"))
            {
                Require(preDomain.admissible,"archived pre-k state is incompatible with current production entry: "+preDomain.failure);
                std::vector<int> free;for(int i=0;i<patch.size[patch.K()];++i)
                    if(!data.fixedBoundary.count(patch.offset[patch.K()]+i))free.push_back(patch.offset[patch.K()]+i);
                const auto expected=Read(packet,"transport_rows");Require(free.size()==expected.size(),"archived free-row size differs");
                for(int i=0;i<free.size();++i)Require(free[i]==expected[i],"archived free-row order differs");
                Require(RANSWallR3::Hash(state)==packet.get<std::string>("raw_sha256"),"archived raw coefficient vector differs");
            }
            const auto bounded=patch.ApplyPimpleKBound(state,pre,time,data,packet.get<bool>("projection_only",false));
            double nonK=0,essential=0;
            for(int i=0;i<patch.Size();++i)if(i<patch.offset[patch.K()] || i>=patch.offset[patch.Epsilon()])nonK=std::max(nonK,std::abs(state[i]-bounded.state[i]));
            for(auto [id,value]:data.fixedBoundary)if(id>=patch.offset[patch.K()] && id<patch.offset[patch.Epsilon()])essential=std::max(essential,std::abs(state[id]-bounded.state[id]));
            Vector dk(patch.size[patch.K()]);for(int i=0;i<dk.size();++i)dk[i]=bounded.state[patch.offset[patch.K()]+i]-state[patch.offset[patch.K()]+i];
            double square=0,integral=0,jump=0;
            for(const auto &q:patch.Quadrature(patch.K())){const double x=q.Value(dk);square+=q.weight*x*x;integral+=q.weight*x;}
            const auto f=patch.fields[patch.K()];
            for(int i=0;i<f->Elements().size();++i)for(int j=i+1;j<f->Elements().size();++j)
            {
                const auto &a=f->Elements()[i],&b=f->Elements()[j];
                if(a.upper[0]!=b.lower[0] || a.lower[2]!=b.lower[2] || a.upper[2]!=b.upper[2])continue;
                for(int k=0;k<=6;++k){const RANSWallCoordinate x{a.upper[0],0,a.lower[2]+k*(a.upper[2]-a.lower[2])/6};
                    jump=std::max(jump,std::abs(f->ElementPoint(i,x).Value(bounded.state,patch.offset[patch.K()])-f->ElementPoint(j,x).Value(bounded.state,patch.offset[patch.K()])));}
            }
            std::vector<std::string> witnessValues;
            if(packet.get_child_optional("witnesses"))for(const auto &item:packet.get_child("witnesses"))
            {
                RANSWallCoordinate x{};int i=0;for(const auto &v:item.second)x[i++]=v.second.get_value<double>();
                const auto row=patch.Point(patch.K(),x);
                witnessValues.push_back(J::Object({{"coordinate",J::Numbers(x)},
                    {"raw",J::Number(row.Value(state,patch.offset[patch.K()]))},
                    {"corrected",J::Number(row.Value(bounded.state,patch.offset[patch.K()]))}}));
            }
            const auto domain=patch.ScalarDomain(bounded.state,time,data,true);
            std::string auditFailure;RANSWallEvolutionPatch::Diagnostic completeAudit;
            bool hardAdmissible=false;std::vector<std::string> hardReasons;
            if(bounded.success)try
            {
                completeAudit=patch.CompleteRestrictedAudit(bounded.state,history,alpha,time,data);
                const double K0=p.model.Q()*us*us;
                auto reject=[&](bool valid,const std::string &reason){if(!valid)hardReasons.push_back(J::Quote(reason));};
                reject(completeAudit.at(patch.AcceptedAWF()?"accepted_awf_scalar_domain":"r2_scalar_domain")[0]==1,"complete_scalar_domain_exit");
                reject(completeAudit.at("admissible")[0]==1,"sampled_scalar_or_sampler_inadmissibility");
                for(const std::string name:{"source_volume_inner","source_volume_outer","physical_sample_K","physical_sample_outer",
                    "prescribed_source_inner","mass_volume","interface_inner","interface_outer","reconstruction_dense","reconstruction_native"})
                {const auto &row=completeAudit.at(name);reject(!row[0] || row[3]>=-1e-10*K0,"covariance:"+name);}
                const auto repair=completeAudit.find("robust_source");
                reject(repair==completeAudit.end() || repair->second[6]>=.5,"lambda_applicability_exit");
                const auto &cov=completeAudit.at(patch.AcceptedAWF()?"accepted_awf_covariance":"r2_covariance");
                const bool sourceValid=patch.AcceptedAWF()?completeAudit.at("accepted_awf_hard_admissible")[0]==1 && completeAudit.at("accepted_awf_scalar_domain")[8]<=1e-10:
                    completeAudit.at("r2_scalar_domain")[8]>=.5;
                reject(cov[0]>=-1e-10*K0 && cov[1]>=-1e-10*K0 && sourceValid,
                    "PIMPLE_hard_covariance_or_normalization_exit");
                hardAdmissible=hardReasons.empty();
            }
                catch(const std::exception &error){auditFailure=error.what();}
            std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote("k_bound")},
                {"fixture_definition",fixtureDefinition},
                {"replay_context",replayContext},{"history_sha256",J::Quote(RANSWallR3::Hash(history))},
                {"initialization_lift_calls_iterations",J::Numbers(patch.InitializationWork())},{"pre_k_inventory",patch.ScalarDomainJSON(preDomain)},{"witness_values",J::List(witnessValues)},
                {"scalar_domain",patch.ScalarDomainJSON(domain,true)},{"complete_audit_failure",J::Quote(auditFailure)},
                {"complete_hard_admissible",J::Boolean(hardAdmissible)},{"complete_hard_reasons",J::List(hardReasons)},
                {"complete_audit",Maps(completeAudit)},
                {"complete_audit_context",patch.ConstitutiveFailureContext().empty()?"null":patch.ConstitutiveFailureContext()},{"offsets",J::Numbers(patch.offset)},{"k_bound",bounded.report},
                {"quadrature_correction_mass_norm",J::Number(std::sqrt(square))},{"quadrature_correction_integral",J::Number(integral)},
                {"shared_cg_jump",J::Number(jump)},{"raw_native_ground_node_minimum",J::Number(nodeMinimum)},
                {"raw",J::Numbers(state)},{"corrected",J::Numbers(bounded.state)},
                {"non_k_change",J::Number(nonK)},{"essential_k_change",J::Number(essential)}})<<'\n';
            Require(out.good(),"bound replay write failed");world->Finalise();return 0;
        }
        const bool constructed=packet.get<bool>("construct_contact",false);
        Vector beforeTrace;
        if(constructed)
        {
            Require(patch.InnerStressName()=="SaturatedR2","constructed contact is the selected research law only");
            state=ConstructContact(patch,history);beforeTrace=state;
            for(auto [row,value]:data.fixedBoundary)
                if(row>=patch.offset[patch.K()] && row<patch.offset[patch.Epsilon()])
                    Require(state[row]==value,"constructed polynomial changes an external K condition");
            state=patch.Lift(state,time,data);
            for(int i=0;i<patch.size[patch.K()];++i)
                Require(state[patch.offset[patch.K()]+i]==beforeTrace[patch.offset[patch.K()]+i],"trace construction replaced k");
        }
        const auto domain=patch.ScalarDomain(state,time,data,true);Require(domain.admissible,"resolution state outside native scalar evaluation domain");
        const auto resolution=patch.ResidualResolution(state,history,alpha,time,data,false);
        const RANSWallCoordinate original{1.6046496089512394,0,0},zero{1,0,0};
        std::vector<std::pair<std::string,std::string>> rateEntries;
        if(packet.get<bool>("compute_scalar_rate",false) || constructed)
        {
            const auto kr=RANSWallResolutionMassRate(patch,resolution,patch.K(),original);
            rateEntries.push_back({"k",Rate(patch,kr,patch.K(),{original,zero},resolution.effectiveState)});
            const auto trace=QuietTraceRate(patch,state,kr.rate,time,data);
            rateEntries.push_back({"trace_rate_assumptions_verified",J::Boolean(trace.verified)});
            rateEntries.push_back({"trace_rate_reason",J::Quote(trace.reason)});
            rateEntries.push_back({"trace_projection_node_checks",J::List(trace.nodes)});
            if(trace.verified)
            {
                const auto er=RANSWallResolutionMassRate(patch,resolution,patch.Epsilon(),{1,0,p.geometry.height},trace.epsilonBoundary);
                rateEntries.push_back({"epsilon_with_induced_boundary_mass",Rate(patch,er,patch.Epsilon(),{{1,0,p.geometry.height}},resolution.effectiveState)});
            }
        }
        std::vector<std::string> fields,keys;for(int b=0;b<patch.Gauge();++b){fields.push_back(Field(patch,b));}
        for(const auto &key:RANSWallResolutionSpatialTerms()){keys.push_back(J::Quote(key));}
        const auto &control=patch.TransitionControls();const auto distances=patch.TransitionDistances();
        const auto output=J::Object({{"schema","1"},{"scope",J::Quote("R3 fixed-state arithmetic/component/mass-contact diagnostic only; no solver, trajectory, field repair or root inference")},
            {"dimension","2"},{"ranks","1"},{"mapped","false"},{"wall_policy",J::Quote(patch.WallPolicyName())},{"upstream_transition",J::Boolean(patch.TransitionEnabled() && !patch.UniformEPFlat())},{"chi",patch.AcceptedAWF()?"null":J::Number(p.chi)},{"eta","0"},
            {"inner_stress_model",J::Quote(patch.InnerStressName())},{"model_variant",J::Quote(patch.VariantName())},
            {"robust_parameters",J::Numbers(patch.RobustParameters())},{"robust_reference_flow",J::Numbers(patch.RobustReferenceFlow())},
            {"robust_reference_normal",J::Numbers(patch.RobustReferenceNormal())},{"map_parameters",J::Numbers(MapParameters(patch.Mapping()))},
            {"parameters",J::Numbers(Vector{p.geometry.molecularNu,p.geometry.height,p.geometry.roughness,p.geometry.sampleHeight,p.geometry.tailStart,
                patch.SampleK(history,{0,0,1}),p.model.Q(),p.model.Kappa()})},{"ustar",J::Number(us)},
            {"time",J::Number(time)},{"alpha",J::Number(alpha)},{"state",J::Numbers(state)},{"history",J::Numbers(history)},
            {"initial",J::Numbers(reference)},{"reference_constraints",J::Constraints(patch.ReferenceConstraints())},
            {"packet_history_matches_rebuilt_h0",J::Boolean(historyMatch)},{"construct_contact",J::Boolean(constructed)},
            {"constructed_before_trace",J::Numbers(beforeTrace)},{"constructed_definition",J::Quote(constructed?
                "k=K0*(1-f(x)*(1-g(r)));g=(r/.01)^2 below existing .01 interface and1above;f=1 on[-3,1], quadratic ramps to0at−4/2; exact0/+-K0 tensor coefficients; zero contactx1 differsfromoriginalfixedrow":"none")},
            {"offsets",J::Numbers(patch.offset)},{"fields",J::List(fields)},{"quadrature_orders",J::Numbers(Vector{double(p.tangentialOrder),double(p.normalOrder)})},
            {"transition_controls",J::Numbers(Vector{control.start,control.end,distances[0],distances[1],distances[2],control.state.on,control.state.off,
                control.state.marginOff,control.state.marginOn,control.state.reynoldsOff,control.state.reynoldsOn})},
            {"transition_weights",J::Numbers(control.state.weights)},{"transition_reference",J::Numbers(patch.TransitionReference())},
            {"scalar_domain",RANSWallEvolutionPatch::ScalarDomainJSON(domain,true)},
            {"effective_state",J::Numbers(resolution.effectiveState)},{"constraints",J::Constraints(resolution.constraints)},
            {"spatial_sum_keys",J::List(keys)},{"terms",Maps(resolution.terms)},{"native",Maps(resolution.native)},
            {"compensated",Maps(resolution.compensated)},{"extended_decimal",Wide(resolution.extended)},
            {"absolute_accumulation",Maps(resolution.absoluteAccumulation)},{"addition_counts",Maps(resolution.additionCounts)},
            {"norms",Maps(resolution.norms)},{"extended_precision_bits",J::Number(std::numeric_limits<long double>::digits)},
            {"resolution_scope",J::Quote(resolution.scope)},{"scalar_rates",J::Object(rateEntries)}});
        Require(!std::ifstream(argv[3]).good(),"refusing to overwrite resolution report");
        std::ofstream stream(argv[3]);Require(stream.good(),"cannot open resolution report");stream<<output<<'\n';stream.close();Require(stream.good(),"cannot finish resolution report");
        world->Finalise();return 0;
    }
    catch(const std::exception &error){std::cerr<<error.what()<<std::endl;return 1;}
}
