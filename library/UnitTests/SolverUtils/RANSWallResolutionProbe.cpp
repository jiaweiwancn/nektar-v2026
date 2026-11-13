///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallPhysicsSignature.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallModelSpec.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallKBoundProjection.h>
#include "RANSWallScalarFeasibilityContext.h"
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
std::string MappedSeedWitnesses(const RANSWallEvolutionPatch &patch,const Vector &state,
    const RANSWallEvolutionPatch::Diagnostic &audit)
{
    if(!patch.Composite() || patch.Mapping().amplitude==0)return "null";
    std::vector<std::string> result;
    const auto tensor=[](const RANSTensor &g){std::vector<std::string> rows;for(const auto &r:g)rows.push_back(J::Numbers(r));return J::List(rows);};
    for(bool inner:{true,false})
    {
        const auto &point=audit.at(inner?"accepted_awf_inner_minimum_point":"accepted_awf_outer_minimum_point");
        const auto &provenance=audit.at(inner?"native_inner_minimum_geometry_element_side_block":"native_outer_minimum_geometry_element_side_block");
        const RANSWallCoordinate p{point[0],point[1],point[2]};
        for(int element=0;element<patch.fields[0]->Elements().size();++element)
        {
        const auto &el=patch.fields[0]->Elements().at(element);const int geometry=el.geometry;
        if((el.upper[2]<=patch.Parameters().geometry.height)!=inner)continue;
        bool contains=true;for(int a:{0,2})contains=contains && p[a]>=el.lower[a]-1e-12 && p[a]<=el.upper[a]+1e-12;
        if(!contains)continue;
        int foot=0;
        for(double f:{-patch.Mapping().lengthX,patch.Mapping().lengthX})if(std::abs(p[0]-f)<1e-12)foot=(el.lower[0]+el.upper[0])/2<f?-1:1;
        const auto metric=patch.GeometryAt(p,inner?-1:1,foot);
        const auto gradient=[&](double scale)
        {
            RANSTensor native{},physical{};const auto f0=patch.CompositeSeedVelocity(p,foot);
            for(int axis:{0,2})
            {
                const double h=scale*std::min(axis==0?1e-7:1e-8,(el.upper[axis]-el.lower[axis])*1e-4);
                auto a=p,b=p;
                if(p[axis]-h<el.lower[axis])
                {a[axis]+=h;b[axis]+=2*h;const auto fa=patch.CompositeSeedVelocity(a),fb=patch.CompositeSeedVelocity(b);for(int j=0;j<3;++j)native[j][axis]=(-3*f0[j]+4*fa[j]-fb[j])/(2*h);}
                else if(p[axis]+h>el.upper[axis])
                {a[axis]-=h;b[axis]-=2*h;const auto fa=patch.CompositeSeedVelocity(a),fb=patch.CompositeSeedVelocity(b);for(int j=0;j<3;++j)native[j][axis]=(3*f0[j]-4*fa[j]+fb[j])/(2*h);}
                else
                {a[axis]+=h;b[axis]-=h;const auto fa=patch.CompositeSeedVelocity(a),fb=patch.CompositeSeedVelocity(b);for(int j=0;j<3;++j)native[j][axis]=(fa[j]-fb[j])/(2*h);}
            }
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)physical[i][j]+=native[i][k]*metric.inverse[k][j];
            return physical;
        };
        const auto first=gradient(1),target=gradient(.5);RANSTensor represented{};double error=0;
        for(int b=0;b<2;++b)represented[patch.PhysicalAxis(b)]=patch.GeometryPoint(0,geometry,p).Gradient(state,patch.offset[b]);
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)error=std::max(error,std::abs(first[i][j]-target[i][j]));
        const double k=point[3],epsilon=inner?0:point[4];
        const auto stress=inner?patch.InnerStress(k,point[5],target).reynoldsStress:patch.OuterClosure(k,epsilon,target,p).reynoldsStress;
        const auto actualStress=inner?patch.InnerStress(k,point[5],represented).reynoldsStress:patch.OuterClosure(k,epsilon,represented,p).reynoldsStress;
        result.push_back(J::Object({{"inner",J::Boolean(inner)},{"coordinate",J::Numbers(p)},{"physical_coordinate",J::Numbers(patch.Physical(p))},
            {"geometry",J::Number(geometry)},{"element",J::Number(element)},{"represented_gradient",tensor(represented)},
            {"target_gradient",tensor(target)},{"target_gradient_step_halving_max",J::Number(error)},
            {"target_divergence",J::Number(target[0][0]+target[1][1]+target[2][2])},{"target_tensor",tensor(stress)},
            {"target_eigenvalues",J::Numbers(J::Eigenvalues(stress))},{"native_eigenvalues",J::Numbers(J::Eigenvalues(actualStress))},
            {"native_tensor",tensor(actualStress)},{"audit_minimum_eigenvalue",J::Number(point.back())},
            {"is_audit_minimizer",J::Boolean(geometry==int(provenance[0]))},
            {"k",J::Number(k)},{"epsilon",J::Number(epsilon)},
            {"scope",J::Quote("unprojected physical streamfunction velocity at the SAME represented scalar fields and recorded one-sided element; independent step-halved chart differentiation; not an equilibrium residual")}}));
        }
    }
    return J::List(result);
}
std::string ExtraSeedAudit(const RANSWallEvolutionPatch &patch,const Vector &s,const RANSWallEvolutionData &data,const Vector *base=nullptr,int order=12)
{
    double minimum=1e300,innerMinimum=1e300,outerMinimum=1e300,strong=0,volume=0,energy=0,baseEnergy=0,deltaNorm=0,maximumIncrement=0;double tke=0,baseTke=0,outerTke=0,baseOuterTke=0,tmin=1e300,tmax=0,lmin=1e300,lmax=0,trmin=1e300,trmax=0,lrmin=1e300,lrmax=0,maxDivergence=0;std::string witness="null";int count=0;
    double krmin=1e300,krmax=-1e300,ermin=1e300,ermax=-1e300,minK=1e300,minE=1e300;int invalidScalars=0;std::string scalarWitness="null",envelopeWitness="null";
    for(const auto &q:patch.Quadrature(0,order))
    {
        const bool inner=patch.fields[0]->Elements()[q.element].upper[2]<=patch.Parameters().geometry.height;
        RANSTensor g{};RANSVector u{};for(int a=0;a<2;++a){g[patch.PhysicalAxis(a)]=q.Gradient(s,patch.offset[a]);u[patch.PhysicalAxis(a)]=q.Value(s,patch.offset[a]);}
        const double k=patch.GeometryPoint(patch.K(),q.geometry,q.coordinate).Value(s,patch.offset[patch.K()]);
        const auto profile=patch.StateProfile(s,q.coordinate,0,data);
        const double e=inner?profile->Dissipation(q.coordinate[2]):patch.GeometryPoint(patch.Epsilon(),q.geometry,q.coordinate).Value(s,patch.offset[patch.Epsilon()]);
        minK=std::min(minK,k);minE=std::min(minE,e);
        if(!(std::isfinite(k) && std::isfinite(e) && k>0 && e>0))
        {++invalidScalars;scalarWitness=J::Object({{"coordinate",J::Numbers(q.coordinate)},{"physical",J::Numbers(patch.Physical(q.coordinate))},{"geometry",J::Number(q.geometry)},{"element",J::Number(q.element)},{"inner",J::Boolean(inner)},{"k",J::Number(k)},{"epsilon",J::Number(e)}});continue;}
        if(!inner && base)
        {
            const double kb=patch.GeometryPoint(patch.K(),q.geometry,q.coordinate).Value(*base,patch.offset[patch.K()]);
            const double eb=patch.GeometryPoint(patch.Epsilon(),q.geometry,q.coordinate).Value(*base,patch.offset[patch.Epsilon()]);
            const double rk=k/kb,re=e/eb;krmin=std::min(krmin,rk);krmax=std::max(krmax,rk);ermin=std::min(ermin,re);ermax=std::max(ermax,re);
            if(std::min(rk,re)<.5-1e-12 || std::max(rk,re)>2+1e-12)envelopeWitness=J::Object({{"coordinate",J::Numbers(q.coordinate)},{"physical",J::Numbers(patch.Physical(q.coordinate))},{"geometry",J::Number(q.geometry)},{"element",J::Number(q.element)},{"k_over_cold",J::Number(rk)},{"epsilon_over_cold",J::Number(re)}});
        }
        const auto R=inner?patch.InnerStress(k,profile->Original().EddyViscosity(q.coordinate[2]),g).reynoldsStress:patch.OuterClosure(k,e,g,q.coordinate).reynoldsStress;
        tke+=q.weight*k;const double bk=base?patch.GeometryPoint(patch.K(),q.geometry,q.coordinate).Value(*base,patch.offset[patch.K()]):k;baseTke+=q.weight*bk;
        if(!inner){outerTke+=q.weight*k;baseOuterTke+=q.weight*bk;const double ts=k/e,ls=std::pow(k,1.5)/e;tmin=std::min(tmin,ts);tmax=std::max(tmax,ts);lmin=std::min(lmin,ls);lmax=std::max(lmax,ls);
            const double be=base?patch.GeometryPoint(patch.Epsilon(),q.geometry,q.coordinate).Value(*base,patch.offset[patch.Epsilon()]):e;trmin=std::min(trmin,ts/(bk/be));trmax=std::max(trmax,ts/(bk/be));lrmin=std::min(lrmin,ls/(std::pow(bk,1.5)/be));lrmax=std::max(lrmax,ls/(std::pow(bk,1.5)/be));}

        if(base){double b2=0,d2=0;for(int a=0;a<2;++a){const double b=q.Value(*base,patch.offset[a]),du=u[patch.PhysicalAxis(a)]-b;b2+=b*b;d2+=du*du;maximumIncrement=std::max(maximumIncrement,std::abs(du));}baseEnergy+=q.weight*b2/2;deltaNorm+=q.weight*d2;}
        const auto ev=J::Eigenvalues(R);if(inner)innerMinimum=std::min(innerMinimum,ev[0]);else outerMinimum=std::min(outerMinimum,ev[0]);double div=g[0][0]+g[2][2];maxDivergence=std::max(maxDivergence,std::abs(div));strong+=q.weight*div*div;volume+=q.weight;energy+=q.weight*(u[0]*u[0]+u[2]*u[2])/2;++count;
        if(ev[0]<minimum){minimum=ev[0];witness=J::Object({{"coordinate",J::Numbers(q.coordinate)},{"physical",J::Numbers(patch.Physical(q.coordinate))},
            {"geometry",J::Number(q.geometry)},{"element",J::Number(q.element)},{"inner",J::Boolean(inner)},{"k",J::Number(k)},{"epsilon",J::Number(e)},
            {"gradient",J::Rows(g)},{"tensor",J::Rows(R)},{"eigenvalues",J::Numbers(ev)}});}
    }
    return J::Object({{"minimum_covariance",J::Number(minimum)},{"inner_minimum_covariance",J::Number(innerMinimum)},{"outer_minimum_covariance",J::Number(outerMinimum)},{"witness",witness},{"count",J::Number(count)},
        {"quadrature_order",J::Number(order)},{"state_sha256",J::Quote(RANSWallR3::Hash(s))},{"cold_base_sha256",base?J::Quote(RANSWallR3::Hash(*base)):"null"},{"model_specification",patch.ModelMetadata()},
        {"minimum_k",J::Number(minK)},{"minimum_epsilon",J::Number(minE)},{"invalid_scalar_count",J::Number(invalidScalars)},{"scalar_positive_finite",J::Boolean(invalidScalars==0)},{"invalid_scalar_witness",scalarWitness},{"outer_k_over_cold",J::Numbers(Vector{krmin,krmax})},{"outer_epsilon_over_cold",J::Numbers(Vector{ermin,ermax})},{"scalar_envelopes_passed",J::Boolean(base && invalidScalars==0 && krmin>=.5-1e-12 && krmax<=2+1e-12 && ermin>=.5-1e-12 && ermax<=2+1e-12)},{"scalar_envelope_witness",envelopeWitness},{"integrals_complete",J::Boolean(invalidScalars==0)},{"strong_divergence_RMS",J::Number(std::sqrt(strong/volume))},{"maximum_absolute_divergence",J::Number(maxDivergence)},{"TKE_integral",J::Number(tke)},{"base_TKE_integral",J::Number(baseTke)},{"outer_TKE_integral",J::Number(outerTke)},{"base_outer_TKE_integral",J::Number(baseOuterTke)},{"outer_time_scale",J::Numbers(Vector{tmin,tmax})},{"outer_length_scale",J::Numbers(Vector{lmin,lmax})},{"outer_time_scale_ratio",J::Numbers(Vector{trmin,trmax})},{"outer_length_scale_ratio",J::Numbers(Vector{lrmin,lrmax})},{"strong_divergence_L2",J::Number(std::sqrt(strong))},{"volume",J::Number(volume)},{"kinetic_energy",J::Number(energy)},{"base_kinetic_energy",base?J::Number(baseEnergy):"null"},{"velocity_change_L2",base?J::Number(std::sqrt(deltaNorm)):"null"},{"maximum_sampled_velocity_increment",base?J::Number(maximumIncrement):"null"},
        {"scope",J::Quote("independent panelwise physical tensor sampling; complements complete inventory; no domain-wide covariance theorem")}});
}
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
        LoadRANSWallModel(serial,p);
        RANSWallEvolutionPatch patch(serial,graph,world->GetComm(),p,LoadRANSWallNormalMap(serial,2,p.geometry.height));
        {std::ofstream readback(std::string(argv[3])+".model.json");readback<<patch.ModelMetadata()<<std::endl;Require(readback.good(),"native model readback output");}
        const auto operation=packet.get<std::string>("operation","");
        if(auto expected=packet.get_optional<std::string>("expected_model"))Require(*expected==(p.model.Model()==SolverUtils::RANSModel::CLSHighRe?"CLSHighRe":"StandardKEpsilon"),"probe/model-manifest mismatch");
        Require(patch.Dimension()==2 && (!patch.Mapped() || (patch.Composite() && (operation=="h2b_operators" || operation=="h2e_contract_checks" || operation=="h2h_scalar_context" || operation=="h2h_candidate" || operation=="h2g_frozen" || operation=="h2f_context" || operation=="h2f_native" || operation=="h2e_context" || operation=="h2e_candidate" || operation=="h2d_setup" || operation=="h2d_candidate" || operation=="h2d_flat" || operation=="h2c_seed" || operation=="h2b_checkpoint" || operation=="h2c_checkpoint"))) && world->GetComm()->GetSize()==1 && patch.MatchingEnabled(),"resolution scope is a supported serial flat2D wall formulation");
        Require(patch.InnerStressName()==packet.get<std::string>("inner_stress_model"),"packet law differs from session");
        for(const std::string name:{"BodyForceInner","BodyForceOuter","ScalarSourceInner","ScalarSourceOuter"})
            Require(!serial->DefinesFunction(name),"R3 resolution input must retain unforced physical data");
        double us=0;serial->LoadParameter("RANSWallFrictionVelocity",us,.4);
        RANSWallEvolutionData data(2);Vector reference,formedHistory;std::string replayContext="null";
        if(operation=="h2b_operators" || operation=="h2c_seed" || operation=="h2e_contract_checks" || operation=="h2h_scalar_context" || operation=="h2h_candidate" || operation=="h2g_frozen" || operation=="h2f_context" || operation=="h2f_native" || operation=="h2e_context" || operation=="h2e_candidate" || operation=="h2d_setup" || operation=="h2d_candidate" || operation=="h2d_flat")
        {
            Require(patch.PrescribedInflow() && !packet.get_optional<std::string>("native_checkpoint") && RANSWallR3::HashFile(argv[1])==packet.get<std::string>("xml_sha256"),"H2b operator fixture XML/initialization contract differs");
            auto state=patch.FiniteInflowSeed();const auto baseState=state;
            if(operation=="h2e_contract_checks")
            {
                data.fixedBoundary=patch.FixedReferenceBoundary(baseState);patch.SetReference(baseState,data,0,0);
                const auto rows=patch.InitializerRows(baseState);Vector scalar=Read(packet,"state");
                const auto zero=patch.ValidateOuterInitializer(baseState,baseState,data,true);
                const auto scalarContract=patch.ValidateOuterInitializer(scalar,baseState,data,true);
                std::vector<std::string> rejected;
                auto reject=[&](Vector v,const std::string &name,bool coupled=true){bool caught=false;try{patch.ValidateOuterInitializer(v,baseState,data,coupled);}catch(const std::exception&){caught=true;}Require(caught,"H2e contract fixture accepted "+name);rejected.push_back(J::Quote(name));};
                auto v=baseState;v[patch.offset[patch.K()]]+=1e-7;reject(v,"inner_k");
                v=baseState;v[patch.offset[patch.Epsilon()]]+=1e-7;reject(v,"interface_epsilon");
                v=baseState;v[patch.offset[patch.InnerPressure()]]+=1e-7;reject(v,"pressure");
                const auto fixed=patch.FixedReferenceBoundary(baseState);const int f=fixed.begin()->first;v=baseState;v[f]+=1e-7;reject(v,"inlet");
                std::set<int> allowed(rows.velocity.begin(),rows.velocity.end());allowed.insert(rows.k.begin(),rows.k.end());allowed.insert(rows.epsilon.begin(),rows.epsilon.end());bool found=false;
                for(int i=0;i<baseState.size();++i)if(!allowed.count(i) && baseState[i]==0.){v=baseState;v[i]=std::copysign(0.,std::signbit(v[i])?1.:-1.);reject(v,"frozen_signed_zero");found=true;break;}
                Require(found,"no fixed zero for bit fixture");
                reject(scalar,"old_frozen_scalar_contract",false);
                std::ofstream out(argv[3]);out<<J::Object({{"passed","true"},{"zero_amplitude_contract",zero},{"scalar_contract",scalarContract},{"rejected",J::List(rejected)},{"governing_solves","0"},{"physical_commits","0"}})<<std::endl;world->Finalise();return 0;
            }
            if(operation=="h2e_context" || operation=="h2f_context" || operation=="h2f_native")
            {
                Require(RANSWallR3::Hash(state)==packet.get<std::string>("base_sha256"),"H2e saved basis base differs");
                const auto path=packet.get<std::string>("basis_file");Require(RANSWallR3::HashFile(path)==packet.get<std::string>("basis_sha256"),"H2e saved basis checksum");
                const int count=packet.get<int>("velocity_direction_count",32);Require(count>=1 && count<=48,"initializer direction count outside bounded scope");
                RANSWallEvolutionPatch::WallMatrix basis(count,Vector(state.size()));std::ifstream in(path,std::ios::binary);
                for(auto &row:basis){in.read(reinterpret_cast<char*>(row.data()),row.size()*sizeof(double));}
                Require(in.good() && in.peek()==std::char_traits<char>::eof(),"H2e basis dimensions");
                data.fixedBoundary=patch.FixedReferenceBoundary(baseState);patch.SetReference(baseState,data,0,0);
                std::string context;
                if(operation=="h2f_native")
                {
                    std::vector<RANSWallEvolutionPatch::InitializerWitness> witnesses;
                    for(const auto &entry:packet.get_child("witnesses")){const auto &w=entry.second;const auto c=Read(w,"chart"),v=Read(w,"eigenvector");Require(c.size()==3 && v.size()==3,"native witness dimensions");witnesses.push_back({{c[0],c[1],c[2]},w.get<int>("geometry"),{v[0],v[1],v[2]}});}
                    context=patch.InitializerNativeData(state,Read(packet,"start_state"),basis,witnesses,data,packet.get<std::string>("context_directory"));
                }
                else context=patch.CoupledInitializerContext(state,data,basis,packet.get<std::string>("context_directory"),operation=="h2f_context");
                std::ofstream out(argv[3]);out<<J::Object({{"context",context},{"patch_signature",J::Quote(patch.Signature())},{"geometry_identity",J::Quote(patch.GeometrySignature())},{"inventory_identity",J::Quote(patch.PimpleKBoundSignature())},{"reference_sha256",J::Quote(RANSWallR3::Hash(patch.Reference()))},{"defect_sha256",J::Quote(RANSWallR3::Hash(patch.Defect()))}})<<std::endl;
                Require(out.good(),"H2e context manifest output");world->Finalise();return 0;
            }
            if(operation=="h2g_frozen")
            {
                if(packet.get_optional<std::string>("state_sha256")){state=Read(packet,"state");Require(state.size()==baseState.size() && RANSWallR3::Hash(state)==packet.get<std::string>("state_sha256"),"frozen state identity");}
                data.fixedBoundary=patch.FixedReferenceBoundary(baseState);patch.SetReference(baseState,data,0,0);
                std::ofstream out(argv[3]);out<<patch.FrozenConstitutiveData(state,data,packet.get<std::string>("context_directory"))<<std::endl;Require(out.good(),"frozen output");world->Finalise();return 0;
            }
            if(operation=="h2h_scalar_context")
            {
                Require(RANSWallR3::Hash(baseState)==packet.get<std::string>("base_sha256"),"H2h cold base differs");
                const auto b=Read(packet,"state");Require(RANSWallR3::Hash(b)==packet.get<std::string>("state_sha256"),"H2h archived B differs");
                data.fixedBoundary=patch.FixedReferenceBoundary(baseState);patch.SetReference(baseState,data,0,0);
                const auto context=RANSH2h::Export(patch,b,baseState,data,packet);
                std::ofstream out(argv[3]);out<<J::Object({{"context",context},{"patch_signature",J::Quote(patch.Signature())},{"geometry_identity",J::Quote(patch.GeometrySignature())},{"inventory_identity",J::Quote(patch.PimpleKBoundSignature())},{"reference_sha256",J::Quote(RANSWallR3::Hash(patch.Reference()))},{"defect_sha256",J::Quote(RANSWallR3::Hash(patch.Defect()))}})<<std::endl;
                Require(out.good(),"H2h scalar context manifest write");world->Finalise();return 0;
            }
            const bool scalarCandidate=operation=="h2e_candidate" || operation=="h2h_candidate";
            std::string outerContract="null";
            if(scalarCandidate)
            {
                Require(RANSWallR3::Hash(baseState)==packet.get<std::string>("base_sha256"),"H2e candidate base mismatch");state=Read(packet,"state");
                Require(RANSWallR3::Hash(state)==packet.get<std::string>("state_sha256"),"H2e candidate checksum");
                data.fixedBoundary=patch.FixedReferenceBoundary(baseState);patch.SetReference(baseState,data,0,0);
                if(operation=="h2h_candidate")
                {
                    const auto b=Read(packet,"B_state");Require(RANSWallR3::Hash(b)==packet.get<std::string>("B_sha256"),"H2h exact B checksum");
                    RANSH2h::FixedB(patch,state,b,baseState);
                }
                outerContract=patch.ValidateOuterInitializer(state,baseState,data,packet.get<bool>("coupled",false));
            }
            if(operation=="h2d_candidate")
            {
                const auto supplied=Read(packet,"state");Require(supplied.size()==state.size() && RANSWallR3::Hash(state)==packet.get<std::string>("base_sha256"),"candidate base identity mismatch");
                for(int i=patch.offset[patch.K()];i<state.size();++i)Require(std::memcmp(&supplied[i],&state[i],sizeof(double))==0,"velocity-only candidate changed scalar or pressure bits");
                for(auto [row,value]:patch.FixedReferenceBoundary(state))Require(std::memcmp(&supplied[row],&state[row],sizeof(double))==0,"candidate changed essential coefficient bits");
                state=supplied;Require(RANSWallR3::Hash(state)==packet.get<std::string>("state_sha256"),"candidate hash mismatch");
            }
            if(operation=="h2d_flat")
            {data.fixedBoundary=patch.FixedReferenceBoundary(state);patch.SetReference(state,data,0,0);std::ofstream out(argv[3]);out<<patch.FlatSpatialDiagnostic(state,data)<<std::endl;world->Finalise();return 0;}
            const auto hash=RANSWallR3::Hash(state);data.fixedBoundary=patch.FixedReferenceBoundary(scalarCandidate?baseState:state);
            const auto domain=patch.ScalarDomain(state,0,data,true);
            Require(domain.admissible,"H2b setup scalar/budget gate: "+patch.ScalarDomainJSON(domain,true));
            std::string h2hContract="null";
            if(operation=="h2h_candidate")
            {
                const auto b=Read(packet,"B_state");RANSH2h::FixedB(patch,state,b,baseState);
                const auto coldDomain=patch.ScalarDomain(baseState,0,data,true);
                const auto ph=RANSH2h::ProfileHash(domain);Require(ph==RANSH2h::ProfileHash(coldDomain),"H2h candidate changed fixed AWF profile");
                h2hContract=J::Object({{"non_scalar_bits_exact","true"},{"profile_invariant","true"},{"profile_sha256",J::Quote(ph)},{"B_sha256",J::Quote(RANSWallR3::Hash(b))},{"cold_sha256",J::Quote(RANSWallR3::Hash(baseState))}});
            }
            patch.SetReference(scalarCandidate?baseState:state,data,0,0);
            RANSWallEvolutionPatch::EvaluationScope evaluation(patch,data,0);
            const auto checks=operation!="h2b_operators"?patch.FiniteInflowChecks(state):patch.CompositeChecks(state,0,data);const auto audit=patch.CompleteRestrictedAudit(state,state,0,0,data);
            std::string budget;
            try{budget=patch.FiniteTurbulenceBudget(state,state,0,0,data);}
            catch(const std::exception &error)
            {budget=J::Object({{"status",J::Quote("additional_diagnostic_failed")},{"failure",J::Quote(error.what())},
                {"native_operator_and_hard_audit_preserved","true"},{"governing_solves","0"},{"physical_commits","0"}});}
            std::string footFixture="null";
            if(operation=="h2c_seed" && patch.Composite() && patch.Mapping().amplitude!=0)
            {
                try{footFixture=patch.MovingFootChecks(state,data);}
                catch(const std::exception &error)
                {footFixture=J::Object({{"passed","false"},{"failure",J::Quote(error.what())},{"native_seed_audit_preserved","true"}});}
            }
            std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote(operation)},{"checks",checks},{"audit",Maps(audit)},{"feasibility",operation=="h2d_setup"?patch.CompositeFeasibilityData(state,data):"null"},
                {"state",J::Numbers(state)},{"state_sha256",J::Quote(hash)},{"outer_initializer_contract",outerContract},{"H2h_fixed_B_contract",h2hContract},{"fixed_reference",scalarCandidate?J::Numbers(patch.Reference()):"null"},{"fixed_defect",scalarCandidate?J::Numbers(patch.Defect()):"null"},{"reference_sha256",J::Quote(RANSWallR3::Hash(patch.Reference()))},{"defect_sha256",J::Quote(RANSWallR3::Hash(patch.Defect()))},{"extra_sampling",(operation=="h2d_candidate" || scalarCandidate)?ExtraSeedAudit(patch,state,data,&baseState,scalarCandidate?20:12):"null"},{"fresh_sampling",operation=="h2h_candidate"?ExtraSeedAudit(patch,state,data,&baseState,24):"null"},{"seed_setup",patch.CompositeSeedInfo()},{"signed_derivative_checks",patch.Composite() && patch.AnalyticalNeutralInflow() && patch.Mapping().amplitude!=0?patch.CompositeSeedDerivativeChecks():"null"},
                {"scalar_inventory",patch.ScalarDomainJSON(domain)},
                {"seed_target_witnesses",patch.Composite() && patch.Mapping().amplitude!=0?MappedSeedWitnesses(patch,state,audit):"[]"},
                {"moving_foot_fixture",footFixture},
                {"seed_constitutive",patch.ConstitutiveStateDiagnostic(state,0,data,{{.059553093889933906,0,patch.Parameters().geometry.height},{.1,0,.1}})},{"physical_samples",Maps(patch.PhysicalSamples(state,0,data))},
                {"native_budget",budget},
                {"evaluation_reuse",patch.EvaluationReuseJSON()},
                {"patch_signature",J::Quote(patch.Signature())},{"inventory_identity",J::Quote(patch.PimpleKBoundSignature())},{"geometry_identity",J::Quote(patch.GeometrySignature())},{"governing_solves","0"},{"physical_commits","0"}})<<'\n';
            Require(out.good(),"cannot write H2b native operator fixture");world->Finalise();return 0;
        }
        if(operation=="finite_initial_domain")
        {
            Require(patch.FiniteInflowReference() && !packet.get_optional<std::string>("native_checkpoint") && packet.get<double>("time")==0,
                "initial domain operation is an uncommitted finite-profile setup candidate only");
            Require(RANSWallR3::HashFile(argv[1])==packet.get<std::string>("xml_sha256"),"initial candidate XML hash differs");
            const auto state=Read(packet,"state");const auto hash=RANSWallR3::Hash(state);
            Require(state.size()==patch.Size() && hash==packet.get<std::string>("state_sha256"),"initial candidate bytes differ");
            data.fixedBoundary=patch.FixedReferenceBoundary(state);
            const auto domain=patch.ScalarDomain(state,0,data);
            Require(hash==RANSWallR3::Hash(state) && patch.InitializationWork()==Vector({0.,0.}),"initial domain diagnostic modified candidate or called Lift");
            std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote(operation)},
                {"state_sha256",J::Quote(hash)},{"xml_sha256",J::Quote(packet.get<std::string>("xml_sha256"))},
                {"scalar_domain",patch.ScalarDomainJSON(domain,true)},{"governing_solves","0"},{"physical_commits","0"},
                {"native_history_formed","false"},{"selected_state_unchanged","true"},{"initialization_work",J::Numbers(patch.InitializationWork())},
                {"scope",J::Quote("exact retained uncommitted setup candidate; scalar/sampling applicability only; not an accepted checkpoint or migration")}})<<'\n';
            Require(out.good(),"initial domain diagnostic write failed");world->Finalise();return 0;
        }
        const bool finiteH1=operation=="finite_inflow_h1" || operation=="h2b_checkpoint" || operation=="h2c_checkpoint";
        Require(!finiteH1 || ((operation=="h2c_checkpoint"?patch.AnalyticalNeutralInflow():patch.FiniteInflowReference()) && packet.get_optional<std::string>("native_checkpoint")),"finite H1 readback requires its immutable compatible checkpoint");
        Require(operation!="native_epsilon_metric" || packet.get_optional<std::string>("native_checkpoint"),
            "native epsilon metric requires an immutable native checkpoint");
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
            ((patch.AcceptedAWF() && (operation=="accepted_awf" || operation=="height_comparison" || operation=="selected_original_operator")) ||
                operation=="evaluation_reuse_checks" || operation=="native_epsilon_metric" || finiteH1);
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
            Require(((mode=="Neutral" && !patch.FiniteInflowReference()) || (finiteH1 && mode==(operation=="h2c_checkpoint"?"AnalyticalNeutralInflowV1":"FiniteInflowV1"))) && eta==0,
                "native replay requires neutral eta=0, or the explicit finite-H1 immutable operation");
            const bool neutral=mode=="Neutral";
            const auto signature=RANSWallPhysicsSignature(patch,serial,eta,neutral,us,{dx,dy},coordinates=="Chart");
            Require(signature==packet.get<std::string>("native_signature"),"native replay physical signature differs from archived report");
            const auto checkpoint=RANSWallTimeIntegration::ReadCheckpoint(checkpointPath,signature,1,dt,patch.Size(),patch.DifferentialSize());
            bool tailIdentityRejected=false;
            if(packet.get<bool>("test_tail_identity_only",false))
            {
                // Pure reader-contract test. No altered profile, simulation,
                // reference or history is constructed; all real cases keep z_m=H.
                double tail;serial->LoadParameter("RANSWallTailStart",tail);
                double changed=tail*1.01;serial->SetParameter("RANSWallTailStart",changed);
                const auto incompatible=RANSWallPhysicsSignature(patch,serial,eta,neutral,us,{dx,dy},coordinates=="Chart");
                serial->SetParameter("RANSWallTailStart",tail);
                Require(RANSWallPhysicsSignature(patch,serial,eta,neutral,us,{dx,dy},coordinates=="Chart")==signature,
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
        if(operation=="native_epsilon_metric")
        {
            Require(strictSpatialReadback && patch.InitializationWork()==Vector({0.,0.}),
                "native epsilon metric cannot initialize, lift or transfer a selected state");
            const auto before=RANSWallR3::Hash(state);const int block=patch.Epsilon();
            const auto mass=patch.fields[block]->Mass();Vector coefficients(patch.size[block]),constant(patch.size[block]);
            std::vector<bool> assigned(patch.size[block],false);std::vector<std::string> rows;
            for(int i=0;i<coefficients.size();++i)coefficients[i]=state[patch.offset[block]+i];
            double volume=0.;
            for(const auto &q:patch.Quadrature(block,4))
            {
                volume+=q.weight;
                for(int j=0;j<q.ids.size();++j)
                {
                    const int id=q.ids[j];Require(!assigned[id] || constant[id]==q.constant[j],"inconsistent native constant representation");
                    assigned[id]=true;constant[id]=q.constant[j];
                }
                rows.push_back(J::Object({{"geometry",J::Number(q.geometry)},{"element",J::Number(q.element)},
                    {"coordinate",J::Numbers(q.coordinate)},{"weight",J::Number(q.weight)},
                    {"ids",J::Numbers(Vector(q.ids.begin(),q.ids.end()))},{"phi",J::Numbers(q.phi)}}));
            }
            double measure=0.;for(int i=0;i<constant.size();++i)for(int j=0;j<constant.size();++j)measure+=constant[i]*mass[i][j]*constant[j];
            Require(std::all_of(assigned.begin(),assigned.end(),[](bool x){return x;}) &&
                std::abs(measure-volume)<1e-11*std::max(1.,volume),"native epsilon mass/physical measure disagree");
            Require(before==RANSWallR3::Hash(state),"native epsilon metric modified selected coefficients");
            std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote(operation)},{"native_replay_context",replayContext},
                {"state_sha256",J::Quote(before)},{"block",J::Number(block)},{"offset",J::Number(patch.offset[block])},
                {"coefficients",J::Numbers(coefficients)},{"matrix",J::Rows(mass)},{"constant_coefficients",J::Numbers(constant)},
                {"physical_volume",J::Number(volume)},{"constant_mass_measure",J::Number(measure)},{"quadrature_rows",J::List(rows)},
                {"field",Field(patch,block)},{"physical_time",J::Number(time)},{"governing_solves","0"},{"physical_commits","0"},
                {"native_history_advanced","false"},{"selected_state_unchanged","true"},
                {"initialization_work",J::Numbers(patch.InitializationWork())},
                {"scope",J::Quote("actual full physical consistent outer-epsilon mass, including prescribed/matching trace coefficients; no constitutive substitution or state correction")}})<<'\n';
            Require(out.good(),"native epsilon metric write failed");world->Finalise();return 0;
        }
        if(operation=="evaluation_reuse_checks")
        {std::ofstream out(argv[3]);out<<patch.EvaluationReuseChecks(state,history,time,data)<<'\n';Require(out.good(),"reuse verification write failed");world->Finalise();return 0;}
        RANSWallEvolutionPatch::EvaluationScope evaluation(patch,data,time);
        if(finiteH1)
        {
            Require(strictSpatialReadback && patch.InitializationWork()==Vector({0.,0.}),"finite H1 readback performed setup Lift");
            const auto before=RANSWallR3::Hash(state),referenceHash=RANSWallR3::Hash(patch.Reference());
            const auto domain=patch.ScalarDomain(state,time,data);
            Require(domain.admissible,"finite H1 readback outside scalar/budget domain");
            const auto audit=patch.CompleteRestrictedAudit(state,history,0,time,data);
            std::string turbulence="null",boundaryFixture="null",boundaryStatus="null",stepCorrection="null";
            if(packet.get_optional<std::string>("native_report_file"))
            {
                const auto file=packet.get<std::string>("native_report_file");
                Require(RANSWallR3::HashFile(file)==packet.get<std::string>("native_report_sha256"),"budget native report SHA256 differs");
                Tree saved;boost::property_tree::read_json(file,saved);
                Require(saved.get<std::string>("signature")==packet.get<std::string>("native_signature") && saved.get<double>("time")==time &&
                    RANSWallR3::Hash(Read(saved,"state"))==before && RANSWallR3::Hash(Read(saved,"fixed_reference"))==referenceHash,
                    "budget report does not bind selected native state/reference/signature/time");
                Require(RANSWallR3::Hash(Read(saved,"native_history_times"))==RANSWallR3::Hash(Read(packet,"native_history_times")),"budget native report history times differ");
                const auto &a=saved.get_child("native_history"),&b=packet.get_child("native_history");Require(a.size()==b.size(),"budget history shape");
                auto bi=b.begin();for(const auto &slot:a){Vector x,y;for(const auto &v:slot.second)x.push_back(v.second.get_value<double>());for(const auto &v:bi++->second)y.push_back(v.second.get_value<double>());Require(RANSWallR3::Hash(x)==RANSWallR3::Hash(y),"budget report history differs");}
                const auto h=Read(saved,"last_stage_history");const double aa=h.empty()?0:saved.get<double>("last_stage_alpha");
                Require(h.empty() || saved.get<double>("last_stage_time")==time,"budget callback time differs");
                if(const auto attempts=saved.get_child_optional("pimple_attempts"))if(!attempts->empty())
                {const auto &last=attempts->back().second;
                 Require(last.get<double>("time")==time && last.get<double>("alpha")==aa &&
                    RANSWallR3::Hash(h)==last.get<std::string>("history_sha256"),"budget callback history/alpha/time differs from native attempted stage");
                 if(const auto correction=last.get_child_optional("k_bound"))
                    {std::ostringstream text;boost::property_tree::write_json(text,*correction,false);stepCorrection=text.str();}}
                turbulence=patch.FiniteTurbulenceBudget(state,h.empty()?state:h,aa,time,data);
            }
            if(patch.OpenBoundary())
            {
                boundaryStatus=patch.DevelopingChecks(state,time,data);
                if(packet.get<bool>("boundary_fixture",false))boundaryFixture=patch.DevelopingBoundaryFixture(state,data);
            }
            const auto components=patch.ResidualResolutionComponents(state,state,0,time,data);
            const auto native=patch.Residual(state,state,0,time,data,false,true);
            double residualError=0;
            const auto &resolved=components.terms.at("total_constrained");
            for(int i=0;i<patch.Size();++i)residualError=std::max(residualError,std::abs(resolved[i]-native[i]));
            Require(residualError<2e-9,"finite-reference stable decomposition differs from native residual");
            if(operation=="h2b_checkpoint" || operation=="h2c_checkpoint")
            {
                std::string analytical="null";
                if(operation=="h2c_checkpoint" && !patch.Mapped())
                {
                    const auto exact=patch.ResidualResolutionComponents(state,state,0,time,data,false,false,true);
                    auto free=exact.terms.at("total_constrained");
                    for(auto [row,value]:exact.constraints){(void)value;free[row]=0.;}
                    const auto norms=patch.Norms(free,exact.constraints,us);
                    analytical=J::Object({{"first_four_scaled_free_test_load_norms",J::Numbers(Vector(norms.begin(),norms.begin()+4))},
                        {"definition",J::Quote("analytical AWF/log values AND analytical gradients under native tests/quadrature; essential rows zeroed; not a residual of the retained coefficient field")},
                        {"state_unchanged",J::Boolean(before==RANSWallR3::Hash(state))}});
                }
                std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote(operation)},{"native_replay_context",replayContext},
                    {"state_sha256",J::Quote(before)},{"reference_sha256",J::Quote(referenceHash)},{"physical_samples",Maps(patch.PhysicalSamples(state,time,data))},
                    {"audit",Maps(audit)},{"native_turbulence_budget",turbulence},{"developing_boundary",boundaryStatus},
                    {"selected_spatial_block_norms",J::Numbers(patch.Norms(native,components.constraints,us))},
                    {"analytical_free_test_loads",analytical},
                    {"patch_signature",J::Quote(patch.Signature())},{"inventory_identity",J::Quote(patch.PimpleKBoundSignature())},{"geometry_identity",J::Quote(patch.GeometrySignature())},{"provider",patch.FiniteInflowChecks(state)},
                    {"governing_solves","0"},{"physical_commits","0"},{"initialization_work",J::Numbers(patch.InitializationWork())},
                    {"native_history_advanced","false"},{"selected_state_unchanged","true"}})<<'\n';world->Finalise();return 0;
            }
            const double H=p.geometry.height;const auto endpoint=patch.FixedInflow().AtHeight(H);
            AWFProfileOptions targetOptions;targetOptions.formulation=AWFProfileFormulation::CLSCalibratedBudgetAWF;
            const auto wall=std::make_shared<EquilibriumAWFProfiles>(p.geometry,endpoint.k,p.model,AWFNormalMetric{},targetOptions,CalibratedAWFParameters(p.model));
            std::vector<std::string> metric;
            for(int order:{16,32})
            {
                Vector errors(7); // separate units: U,gradU,k,epsilon, seed U/k/epsilon drift
                auto xRule=LibUtilities::PointsManager()[LibUtilities::PointsKey(8,LibUtilities::eGaussGaussLegendre)];
                auto zRule=LibUtilities::PointsManager()[LibUtilities::PointsKey(order,LibUtilities::eGaussGaussLegendre)];
                const auto xs=xRule->GetZ(),xw=xRule->GetW(),zs=zRule->GetZ(),zw=zRule->GetW();
                for(int block:{0,patch.K(),patch.Epsilon()})for(const auto &element:patch.fields[block]->Elements())
                {
                    std::vector<double> panels{element.lower[2],element.upper[2]};
                    for(const auto &q:patch.FixedInflow().Knots())if(q.z>element.lower[2] && q.z<element.upper[2])panels.push_back(q.z);
                    for(double z:{wall->Reference().cutoff,wall->Original().ViscosityInterceptHeight(),p.geometry.roughness})
                        if(z>element.lower[2] && z<element.upper[2])panels.push_back(z);
                    std::sort(panels.begin(),panels.end());panels.erase(std::unique(panels.begin(),panels.end()),panels.end());
                    for(std::size_t panel=1;panel<panels.size();++panel)for(int ix=0;ix<8;++ix)for(int iz=0;iz<order;++iz)
                    {
                        const double x=(element.lower[0]+element.upper[0])/2+xs[ix]*(element.upper[0]-element.lower[0])/2;
                        const double z=(panels[panel-1]+panels[panel])/2+zs[iz]*(panels[panel]-panels[panel-1])/2;
                        const double weight=xw[ix]*zw[iz]*(element.upper[0]-element.lower[0])*(panels[panel]-panels[panel-1])/4;
                        const auto q=patch.GeometryPoint(block,element.geometry,{x,0,z});
                        const auto f=patch.FixedInflow().At(patch.Physical(q.coordinate));
                        const bool inner=z<H;double exact=0;
                        if(block==0)exact=inner?wall->Original().Velocity(z,{endpoint.U,0},{0,0})[0]:f.U;
                        else if(block==patch.K())exact=inner?wall->ReferenceK(z):f.k;
                        else exact=f.epsilon;
                        const int j=block==0?0:(block==patch.K()?2:3);
                        errors[j]+=weight*std::pow(q.Value(state,patch.offset[block])-exact,2);
                        errors[block==0?4:(block==patch.K()?5:6)]+=weight*std::pow(q.Value(state,patch.offset[block])-q.Value(patch.Reference(),patch.offset[block]),2);
                        if(block==0)
                        {
                            const double grad=inner?wall->Original().VelocityGradient(z,{endpoint.U,0},{0,0})[0]:f.dU;
                            const auto actual=q.Gradient(state,patch.offset[block]);
                            errors[1]+=weight*(std::pow(actual[2]-grad,2)+actual[0]*actual[0]);
                        }
                    }
                }
                for(double &v:errors)v=std::sqrt(std::max(0.,v));
                metric.push_back(J::Object({{"Gauss_order",J::Number(order)},{"U_gradU_k_epsilon_seedU_seedK_seedEpsilon_L2",J::Numbers(errors)}}));
            }
            std::vector<std::string> profiles;
            for(double x:{-.15,0.,.15})for(double z:{0.,.00025,.001,.003,.005,.01,.015,.02,.03,.05,.1,.2,.3,.4,.8})
            {
                const int radialSide=std::abs(z-H)<1e-13?-1:0;
                const auto row=patch.Point(0,{x,0,z},radialSide);const auto f=patch.FixedInflow().AtHeight(z);
                profiles.push_back(J::Numbers(Vector{x,z,double(row.geometry),double(radialSide),row.Value(state,patch.offset[0]),row.Gradient(state,patch.offset[0])[2],
                    patch.Point(patch.K(),{x,0,z}).Value(state,patch.offset[patch.K()]),f.U,f.dU,f.k,f.epsilon}));
            }
            std::vector<std::string> dataQueries;
            if(packet.get_optional<std::string>("permitted_queries_file"))
            {
                const auto file=packet.get<std::string>("permitted_queries_file");
                Require(RANSWallR3::HashFile(file)==packet.get<std::string>("permitted_queries_sha256"),"permitted flat query hash differs");
                std::ifstream input(file);double z,target,weight,bound;std::string name;
                while(input>>z>>name>>target>>weight>>bound)
                {
                    const RANSWallCoordinate point{packet.get<double>("observation_x",-.15),0,z};const auto v=patch.Point(0,point);
                    const double k=patch.Point(patch.K(),point).Value(state,patch.offset[patch.K()]);
                    const bool inner=z<H;const double epsilon=inner?0:patch.Point(patch.Epsilon(),point).Value(state,patch.offset[patch.Epsilon()]);
                    const bool valid=std::isfinite(k) && (inner?k>=0:k>0) && (inner || (epsilon>0 && std::isfinite(epsilon)));
                    double predicted=0,production=0;bool tensorValid=valid;
                    if(valid)
                    {
                        RANSTensor gradient{};gradient[0]=v.Gradient(state,patch.offset[0]);gradient[2]=v.Gradient(state,patch.offset[1]);
                        RANSTensor stress{};
                        if(inner)
                        {
                            const double K=patch.SampleKRow(point).Value(state,patch.offset[patch.K()]);
                            const AnalyticalWallFunction auxiliary({p.geometry.molecularNu,H,p.geometry.roughness,H/2},K,CalibratedAWFParameters(p.model));
                            const auto law=patch.InnerStress(k,auxiliary.EddyViscosity(z),gradient);stress=law.reynoldsStress;production=law.production;
                        }
                        else
                        {
                            const auto law=EvaluateRANSRobustClosure(p.model.Model(),k,epsilon,gradient,z,p.model.Parameters(),0,patch.RegularizedReferenceAt(point));
                            stress=law.reynoldsStress;production=law.production;
                        }
                        tensorValid=stress[0][0]>0 && J::Eigenvalues(stress)[0]>=-1e-12;
                        predicted=name=="U"?v.Value(state,patch.offset[0]):(name=="minusRxz"?-stress[0][2]:(tensorValid?std::sqrt(stress[0][0]):0));
                    }
                    Require(name=="U" || name=="minusRxz" || name=="sigmaU","unregistered experimental query quantity");
                    dataQueries.push_back(J::Object({{"z",J::Number(z)},{"x",J::Number(point[0])},{"velocity_geometry",J::Number(v.geometry)},
                        {"quantity",J::Quote(name)},{"target",J::Number(target)},{"weight",J::Number(weight)},{"extraction_envelope",J::Number(bound)},
                        {"valid",J::Boolean(valid && tensorValid)},{"native_prediction",valid && tensorValid?J::Number(predicted):"null"},
                        {"production",valid && tensorValid?J::Number(production):"null"},{"inside_AWF",J::Boolean(inner)},
                        {"failure",valid && tensorValid?"null":J::Quote("additional diagnostic scalar/covariance query invalid; production acceptance retained")}}));
                }
            }
            Require(before==RANSWallR3::Hash(state) && referenceHash==RANSWallR3::Hash(patch.Reference()),"finite readback modified state/reference");
            std::ofstream out(argv[3]);out<<J::Object({{"operation",J::Quote(operation)},{"native_replay_context",replayContext},
                {"state_sha256",J::Quote(before)},{"reference_sha256",J::Quote(referenceHash)},
                {"selected_state_unchanged","true"},{"provider",patch.FiniteInflowChecks(state)},
                {"physical_samples",Maps(patch.PhysicalSamples(state,time,data))},
                {"complete_audit",Maps(audit)},{"native_turbulence_budget",turbulence},{"native_step_correction",stepCorrection},{"developing_boundary",boundaryStatus},{"boundary_fixture",boundaryFixture},{"residual_resolution_native_max_error",J::Number(residualError)},
                {"selected_spatial_block_norms",J::Numbers(patch.Norms(native,components.constraints,us))},
                {"physical_projection_metrics",J::List(metric)},{"physical_profile_columns",J::Quote("x,z,velocity_geometry,radial_side,U_h,dU_h/dz,k_h,U_ref,dU_ref/dz,k_ref,epsilon_ref")},
                {"physical_profiles",J::List(profiles)},{"permitted_flat_native_queries",J::List(dataQueries)},
                {"projection_metric_columns",J::Quote("U,gradU,k,epsilon,U_from_native_seed,k_from_native_seed,epsilon_from_native_seed; dimensional norms kept separate")},{"initialization_work",J::Numbers(patch.InitializationWork())},
                {"physical_time",J::Number(time)},{"governing_solves","0"},{"physical_commits","0"},{"native_history_advanced","false"},
                {"metric_scope",J::Quote("affine integration split at every native interface, finite-profile knot, AWF cutoff, viscosity intercept and roughness; tangential GL8 and normal GL16/32; represented field versus finite completion; inner target is AWF no-slip U and budget k; epsilon uses raw finite completion, so matching changes remain visible; not curved polynomial exactness")}})<<'\n';
            Require(out.good(),"finite H1 readback write failed");world->Finalise();return 0;
        }
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
                if(packet.get<bool>("neutral_reference_balance",false))
                {
                    Require(strictSpatialReadback,"neutral reference balance requires strict immutable checkpoint provenance");
                    diagnostics.push_back({"neutral_reference_balance",patch.NeutralReferenceBalance(state,time,data)});
                }
                const auto panels=packet.get_child_optional("reporting_breakpoints")?Read(packet,"reporting_breakpoints"):Vector{};
                diagnostics.push_back({"physical",patch.HeightComparison(state,state,time,data,panels,packet.get<bool>("complete_outer_epsilon",false))});
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
