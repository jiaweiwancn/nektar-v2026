///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <LibUtilities/Communication/CommSerial.h>
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
        const double dotK=patch.Point(patch.K(),sample).Value(kRate);
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
    for(int e=0;e<field->Elements().size();++e)
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
        RANSWallEvolutionPatch patch(serial,graph,world->GetComm(),p,LoadRANSWallNormalMap(serial,2,p.geometry.height));
        Require(patch.Dimension()==2 && !patch.Mapped() && world->GetComm()->GetSize()==1 && patch.TransitionEnabled(),"bounded R3 resolution scope is serial flat2D Upstream");
        Require(patch.InnerStressName()==packet.get<std::string>("inner_stress_model"),"packet law differs from session");
        for(const std::string name:{"BodyForceInner","BodyForceOuter","ScalarSourceInner","ScalarSourceOuter"})
            Require(!serial->DefinesFunction(name),"R3 resolution input must retain unforced physical data");
        double us=0;serial->LoadParameter("RANSWallFrictionVelocity",us,.4);
        const auto rawReference=patch.NeutralReference(us,{1.,0.});RANSWallEvolutionData data(2);
        data.fixedBoundary=patch.FixedReferenceBoundary(rawReference);
        const auto reference=patch.Lift(rawReference,0,data);patch.SetReference(reference,data,0,0);
        auto history=Read(packet,"history");auto state=Read(packet,"state");std::string transfer="null";
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
        Require(!packet.get<bool>("require_history_h0",true) || historyMatch,"packet history differs from rebuilt native h0");
        const double time=packet.get<double>("time"),alpha=packet.get<double>("alpha");
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
            std::string failure;try{patch.CompleteRestrictedAudit(state,history,alpha,time,data);}
            catch(const std::exception &e){failure=e.what();}
            std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote("audit_witness")},
                {"failure",J::Quote(failure)},{"context",patch.ConstitutiveFailureContext().empty()?"null":patch.ConstitutiveFailureContext()},
                {"governing_solves","0"}})<<'\n';Require(out.good(),"audit witness write failed");world->Finalise();return 0;
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
            double nodeMinimum=1e300;
            if(!fixture.empty())
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
            std::string auditFailure;
            if(bounded.success)try{patch.CompleteRestrictedAudit(bounded.state,history,alpha,time,data);}
                catch(const std::exception &error){auditFailure=error.what();}
            std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote("k_bound")},
                {"pre_k_inventory",patch.ScalarDomainJSON(preDomain)},{"witness_values",J::List(witnessValues)},
                {"scalar_domain",patch.ScalarDomainJSON(domain,true)},{"complete_audit_failure",J::Quote(auditFailure)},
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
            {"dimension","2"},{"ranks","1"},{"mapped","false"},{"upstream_transition","true"},{"chi",J::Number(p.chi)},{"eta","0"},
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
