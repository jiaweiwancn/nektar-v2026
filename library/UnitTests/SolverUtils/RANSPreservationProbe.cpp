// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallPhysicsSignature.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <LibUtilities/Communication/CommSerial.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <boost/property_tree/json_parser.hpp>
#include <chrono>
#include <fstream>
#include <iostream>
#include <set>
using namespace Nektar;
using namespace Nektar::SolverUtils;
namespace J=Nektar::RANSWallSolveJson;
using Tree=boost::property_tree::ptree;
using Vector=RANSWallAlgebra::Vector;
namespace
{
Vector Values(const Tree &p,const std::string &key)
{Vector v;for(const auto &x:p.get_child(key))v.push_back(x.second.get_value<double>());return v;}
void Require(bool ok,const std::string &message){if(!ok)throw std::runtime_error(message);}
std::string Maps(const RANSWallEvolutionPatch::Diagnostic &values)
{std::vector<std::pair<std::string,std::string>> r;for(const auto &[name,v]:values)r.push_back({name,J::Numbers(v)});return J::Object(r);}
std::string Row(const RANSWallGlobalPoint &p)
{return J::Object({{"ids",J::Numbers(p.ids)},{"phi",J::Numbers(p.phi)},
    {"derivative",J::Rows(p.derivative)},{"geometry",J::Number(p.geometry)}});}
}
int main(int argc,char **argv)
{
    try
    {
        Require(argc==4,"usage: RANSPreservationProbe SESSION.xml REQUEST.json RESPONSE.json");
        const auto start=std::chrono::steady_clock::now();
        Tree request;boost::property_tree::read_json(argv[2],request);
        auto world=LibUtilities::SessionReader::CreateInstance(2,argv);world->InitSession();
        Require(world->GetComm()->GetSize()==1,"preservation verification requires serial execution");
        char name[]="preservation-probe";char *args[]{name,nullptr};
        std::vector<std::string> files{argv[1]};
        auto session=LibUtilities::SessionReader::CreateInstance(1,args,files,LibUtilities::CommSerial::create(1,args));session->InitSession();
        session->SetSolverInfo("GlobalSysSoln","DirectFull");auto graph=SpatialDomains::MeshGraphIO::Read(session);
        RANSWallPatchParameters p;
        session->LoadParameter("RANSWallChi",p.chi,1.);session->LoadParameter("RANSWallHeight",p.geometry.height,2.);
        session->LoadParameter("RANSWallRoughness",p.geometry.roughness,.1);session->LoadParameter("Kinvis",p.geometry.molecularNu,1.5e-5);
        session->LoadParameter("RANSWallSampleHeight",p.geometry.sampleHeight,p.geometry.height);session->LoadParameter("RANSWallTailStart",p.geometry.tailStart,p.geometry.height);
        session->LoadParameter("RANSWallNormalOrder",p.normalOrder,32);session->LoadParameter("RANSWallTangentialOrder",p.tangentialOrder,8);
        std::string model,policy,wall;session->LoadSolverInfo("TurbulenceModel",model,"MISSING");
        session->LoadSolverInfo("RANSKappaPolicy",policy,"Natural");session->LoadSolverInfo("RANSWallFormulation",wall,"MISSING");
        Require(policy=="Natural" || policy=="Prescribed","unknown kappa policy");
        Require(wall=="EPAWF" || wall=="AcceptedBudgetAWF","new-model wall formulation required");
        double kappa;session->LoadParameter("RANSKappa",kappa,.41);
        p.model=RANSResolvedModel::Resolve(ParseRANSModel(model),policy=="Natural"?RANSKappaPolicy::Natural:RANSKappaPolicy::Prescribed,
            kappa,wall=="EPAWF"?RANSWallFormulation::EPAWF:RANSWallFormulation::AcceptedBudgetAWF,p.geometry.molecularNu);
        Require(p.model.NewModel(),"new-model descriptor required; no CLS-default replay");
        for(const auto &name:{"BodyForceInner","BodyForceOuter","ScalarSourceInner","ScalarSourceOuter"})
            Require(!session->DefinesFunction(name),"neutral preservation verifier requires unforced source definitions");
        std::string referenceMode;session->LoadSolverInfo("RANSWallReference",referenceMode,"Neutral");
        Require(referenceMode=="Neutral","neutral verification requires the Neutral reference/boundary contract");
        double eta;session->LoadParameter("RANSWallEta",eta,0.);Require(eta==0,"unanchored verification requires eta zero");
        RANSWallEvolutionPatch patch(session,graph,world->GetComm(),p,LoadRANSWallNormalMap(session,2,p.geometry.height));
        Require(patch.Dimension()==2 && !patch.Mapped() && (patch.AcceptedAWF() || patch.UniformEPFlat()),
            "neutral preservation probe supports accepted or EP plateau references, serial flat2D only");
        RANSWallEvolutionData data(2);Vector state,reference;double time=request.get<double>("time",0.);
        const std::string operation=request.get<std::string>("operation","neutral");
        if(operation=="neutral")
        {
            double us,dx,dy;session->LoadParameter("RANSWallFrictionVelocity",us,.4);
            session->LoadParameter("RANSWallDirectionX",dx,1.);session->LoadParameter("RANSWallDirectionY",dy,0.);
            state=patch.NeutralReference(us,{dx,dy});data.fixedBoundary=patch.FixedReferenceBoundary(state);
            state=patch.Lift(state,time,data);reference=state;
        }
        else
        {
            Require(operation=="replay","unknown read-only verification operation");
            Require(request.get<std::string>("identity")==patch.PreservationIdentity(),"verification physics/mesh/space identity mismatch");
            Require(request.get<std::string>("descriptor")==p.model.Signature(),"verification resolved descriptor mismatch");
            state=Values(request,"state");reference=Values(request,"reference");
            Require(state.size()==std::size_t(patch.Size()) && reference.size()==state.size(),"verification coefficient size mismatch");
            for(const auto &item:request.get_child("fixed_boundary"))
            {auto q=item.second.begin();const int row=q->second.get_value<int>();++q;const double value=q->second.get_value<double>();Require(row>=0 && row<patch.Size(),"invalid fixed boundary row");Require(data.fixedBoundary.emplace(row,value).second,"duplicate fixed boundary row");}
            Require(data.fixedBoundary==patch.FixedReferenceBoundary(reference),"verification external boundary identity mismatch");
        }
        patch.SetReference(reference,data,time,0.);
        std::vector<RANSWallCoordinate> grid;
        if(auto points=request.get_child_optional("grid"))for(const auto &item:*points)
        {RANSWallCoordinate q{};int i=0;for(const auto &v:item.second){Require(i<3,"invalid verification coordinate");q[i++]=v.second.get_value<double>();}Require(i==3,"invalid verification coordinate");grid.push_back(q);}
        std::vector<int> sides;if(auto selected=request.get_child_optional("grid_sides"))for(const auto &v:*selected)sides.push_back(v.second.get_value<int>());
        double receiptUs,receiptDx,receiptDy;session->LoadParameter("RANSWallFrictionVelocity",receiptUs,.4);session->LoadParameter("RANSWallDirectionX",receiptDx,1.);session->LoadParameter("RANSWallDirectionY",receiptDy,0.);
        const auto targetReference=patch.NeutralReference(receiptUs,{receiptDx,receiptDy});
        const auto receipt=patch.InitializationReceipt(state,reference,time,data,data,patch.FixedReferenceBoundary(targetReference));
        const auto stateCopy=state;
        const auto observation=patch.PreservationPacket(state,time,data,request.get<bool>("export_basis",false),grid,sides);
        std::string resolution="null",testBasis="null",constitutive="null";
        if(request.get<bool>("export_constitutive",false))
        {
            std::vector<std::string> rows;
            for(int index=0;index<grid.size();++index)
            {
                const auto &q=grid[index];const int side=sides.empty()?0:sides[index];
                const bool outer=q[2]>p.geometry.height || (q[2]==p.geometry.height && side==1);
                if(!outer)continue;
                const auto v=patch.Point(0,q,side),e=patch.Point(patch.Epsilon(),q,side);
                const double k=v.Value(state,patch.offset[patch.K()]),epsilon=e.Value(state,patch.offset[patch.Epsilon()]);
                RANSTensor gradient{};gradient[0]=v.Gradient(state,patch.offset[0]);gradient[2]=v.Gradient(state,patch.offset[1]);
                const auto law=patch.OuterClosure(k,epsilon,gradient,q);
                const auto pressure=patch.Point(patch.OuterPressure(),q,side).Value(state,patch.offset[patch.OuterPressure()]);
                rows.push_back(J::Object({{"sample_index",J::Number(index)},{"coordinate",J::Numbers(q)},
                    {"side",J::Number(side)},{"geometry",J::Number(v.geometry)},{"raw_gradient",J::Rows(gradient)},
                    {"k",J::Number(k)},{"epsilon",J::Number(epsilon)},{"Pi",J::Number(pressure)},
                    {"R",J::Rows(law.reynoldsStress)},{"Q",J::Rows(law.quadraticStress)},
                    {"P",J::Number(law.production)},{"c_mu",J::Number(law.cMu)},
                    {"nu_t",J::Number(law.eddyViscosity)},{"raw_divergence",J::Number(law.divergence)},
                    {"epsilon_production",J::Number(law.epsilonProduction)},{"epsilon_destruction",J::Number(law.epsilonDestruction)},
                    {"a95_radicand",law.a95RadicandDefined?J::Number(law.a95Radicand):"null"},
                    {"covariance_eigenvalues",J::Numbers(J::Eigenvalues(law.reynoldsStress))}}));
            }
            constitutive=J::Object({{"schema",J::Quote("NativeForwardConstitutiveV1")},
                {"descriptor",J::Quote(p.model.Signature())},{"samples",J::List(rows)},
                {"gradient_convention",J::Quote("A_ij=d_j U_i; DeviatoricStrainV1 for invariants and Q; full R contracted with raw A")}});
        }
        if(request.get<bool>("export_resolution",false))
        {
            const auto r=patch.ResidualResolution(state,state,0.,time,data,false,true);
            std::vector<std::pair<std::string,std::string>> local;
            for(const auto &[key,terms]:r.elementTerms)local.push_back({key,Maps(terms)});
            std::vector<std::string> masses;
            for(int b=0;b<patch.Gauge();++b)masses.push_back(J::Rows(patch.fields[b]->Mass()));
            std::vector<std::string> keys;for(const auto &key:RANSWallResolutionSpatialTerms())keys.push_back(J::Quote(key));
            std::set<double> cuts;for(const auto &e:patch.fields[0]->Elements())
            {cuts.insert(e.lower[2]);cuts.insert(e.upper[2]);}
            std::vector<std::string> fluxes;
            for(double z:cuts)for(int side:{-1,1})for(const auto &f:patch.Surface(0,2,z,side,12))
            {
                const auto &q=f.test;const auto &element=patch.fields[0]->Elements().at(q.element);
                const bool outer=element.lower[2]>=p.geometry.height;
                const double k=q.Value(state,patch.offset[patch.K()]);
                const double pi=patch.GeometryPoint(patch.InnerPressure()+outer,q.geometry,q.coordinate).Value(state,patch.offset[patch.InnerPressure()+outer]);
                RANSTensor g{},stress{};for(int a=0;a<2;++a)g[a==0?0:2]=q.Gradient(state,patch.offset[a]);
                double dk=0,de=0,epsilon=0,epsilonFlux=0;
                if(outer)
                {
                    const auto erow=patch.GeometryPoint(patch.Epsilon(),q.geometry,q.coordinate);
                    epsilon=erow.Value(state,patch.offset[patch.Epsilon()]);
                    const auto c=patch.OuterClosure(k,epsilon,g,q.coordinate);dk=c.kDiffusion;de=c.epsilonDiffusion;
                    epsilonFlux=de*erow.Gradient(state,patch.offset[patch.Epsilon()])[2];
                    for(int a=0;a<3;++a)for(int b=0;b<3;++b)stress[a][b]=-c.reynoldsStress[a][b]+(a==b?2*k/3-pi:0.);
                }
                else
                {
                    const auto profile=patch.StateProfile(state,q.coordinate,time,data);
                    const double nt=profile->Original().EddyViscosity(z);
                    dk=p.geometry.molecularNu+nt/p.model.Parameters().sigmaK;
                    stress=RANSWallInnerMeanStress(patch.InnerStress(k,nt,g),p.geometry.molecularNu,pi);
                }
                const Vector traction{stress[0][2],stress[1][2],stress[2][2]};
                const double vertical=q.Value(state,patch.offset[1]);
                fluxes.push_back(J::Object({{"coordinate",J::Numbers(q.coordinate)},
                    {"weight",J::Number(f.weight)},{"geometry",J::Number(q.geometry)},
                    {"region",outer?"1":"0"},{"side",J::Number(-side)},
                    {"outward_normal_z",J::Number(side)},
                    {"velocity",J::Numbers(Vector{q.Value(state,patch.offset[0]),vertical})},
                    {"k",J::Number(k)},{"epsilon",outer?J::Number(epsilon):"null"},
                    {"Pi",J::Number(pi)},{"velocity_gradient",J::Rows(g)},
                    {"k_gradient",J::Numbers(q.Gradient(state,patch.offset[patch.K()]))},
                    {"full_traction",J::Numbers(traction)},
                    {"k_diffusive_flux",J::Number(dk*q.Gradient(state,patch.offset[patch.K()])[2])},
                    {"k_conservative_flux",J::Number(dk*q.Gradient(state,patch.offset[patch.K()])[2]-vertical*k)},
                    {"epsilon_diffusive_flux",outer?J::Number(epsilonFlux):"null"},
                    {"epsilon_conservative_flux",outer?J::Number(epsilonFlux-vertical*epsilon):"null"},
                    {"test",Row(q)}}));
            }
            resolution=J::Object({{"schema",J::Quote("ModelAwareSpatialResolutionV1")},
                {"terms",Maps(r.terms)},{"native",Maps(r.native)},{"norms",Maps(r.norms)},
                {"constraints",J::Constraints(r.constraints)},{"effective_state",J::Numbers(r.effectiveState)},
                {"physical_term_names",J::List(keys)},{"element_terms",J::Object(local)},
                {"one_sided_faces",J::List(fluxes)},
                {"mass",J::List(masses)},{"alpha","0"},{"scope",J::Quote(r.scope)}});
        }
        if(request.get<bool>("export_test_basis",false))
        {
            std::vector<int> geometries;
            if(auto selected=request.get_child_optional("grid_test_geometries"))
                for(const auto &v:*selected)geometries.push_back(v.second.get_value<int>());
            Require(geometries.empty() || geometries.size()==grid.size(),"test geometry selector count differs");
            std::vector<std::string> rows;
            for(int i=0;i<grid.size();++i)
            {
                const auto &q=grid[i];const int side=sides.empty()?0:sides[i];
                bool outer=q[2]>p.geometry.height || (q[2]==p.geometry.height && side==1);
                if(!geometries.empty())
                {
                    const auto &elements=patch.fields[0]->Elements();
                    const auto parent=std::find_if(elements.begin(),elements.end(),[&](const auto &e){return e.geometry==geometries[i];});
                    Require(parent!=elements.end(),"unknown explicit native test geometry");
                    const bool selectedOuter=parent->lower[2]>=p.geometry.height;
                    if(std::abs(q[2]-p.geometry.height)<1e-12)
                        Require(selectedOuter==(side==1),"test geometry conflicts with declared interface side");
                    outer=selectedOuter;
                }
                auto selected=[&](int block){return geometries.empty()?patch.Point(block,q,side):patch.GeometryPoint(block,geometries[i],q);};
                rows.push_back(J::Object({{"velocity_k",Row(selected(0))},
                    {"pressure",Row(selected(patch.InnerPressure()+outer))},
                    {"epsilon",outer?Row(selected(patch.Epsilon())):"null"}}));
            }
            testBasis=J::List(rows);
        }
        Require(state==stateCopy,"read-only verification modified coefficients");
        const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        Tree observed;std::istringstream input(observation);boost::property_tree::read_json(input,observed);
        const auto &audit=observed.get_child("audit_alpha_zero");
        const bool physical=audit.get_child("accepted_awf_scalar_domain").front().second.get_value<double>()==1 &&
            audit.get_child_optional("accepted_awf_hard_admissible") &&
            audit.get_child("accepted_awf_hard_admissible").front().second.get_value<double>()==1;
        std::string initialHistory;
        if(auto target=request.get_optional<std::string>("write_initial_history"))
        {
            Require(operation=="replay" && time==0. && physical &&
                std::memcmp(state.data(),reference.data(),state.size()*sizeof(double))==0,
                "initial history export requires an admissible identical state/reference at time zero");
            const auto scheme=session->GetTimeIntScheme();
            Require(scheme.method=="IMEX" && scheme.variant.empty() && scheme.order==1,
                "initial history export supports native one-sweep IMEX1 only");
            double dt,us,dx,dy;session->LoadParameter("TimeStep",dt);
            session->LoadParameter("RANSWallFrictionVelocity",us,.4);session->LoadParameter("RANSWallDirectionX",dx,1.);session->LoadParameter("RANSWallDirectionY",dy,0.);
            RANSWallTimeIntegration integrator(1,dt,
                [](const Vector&,double,double)->Vector{throw std::runtime_error("initial history export attempted a governing solve");},
                [&](const Vector &v,double t){auto full=state;std::copy(v.begin(),v.end(),full.begin());return patch.ExplicitBody(full,t,data);});
            integrator.Initialise(Vector(state.begin(),state.begin()+patch.DifferentialSize()),0.);
            RANSWallTimeIntegration::Checkpoint checkpoint;
            std::string coordinates;session->LoadSolverInfo("RANSWallFunctionCoordinates",coordinates,"Physical");
            Require(coordinates=="Chart" || coordinates=="Physical","unknown initialization function coordinates");
            checkpoint.signature=RANSWallPhysicsSignature(patch,session,0.,true,us,{dx,dy},coordinates=="Chart");
            checkpoint.integration=integrator.Capture();checkpoint.coupledState=state;
            checkpoint.fixedReference=reference;checkpoint.fixedDefect=patch.Defect();checkpoint.eta=0.;
            initialHistory=RANSWallTimeIntegration::WriteCheckpoint(*target,checkpoint);
            Require(state==stateCopy,"initial history export modified coefficients");
        }
        const std::string output=J::Object({{"status",J::Quote(physical?"pass":"fail")},{"execution",J::Quote("completed")},
            {"physical_admissibility",J::Quote(physical?"pass":"fail")},
            {"accuracy",J::Quote("not_assessed_by_sampling_adapter")},{"performed_time_steps","0"},{"initial_history",initialHistory.empty()?"null":J::Quote(initialHistory)},
            {"initialization_receipt",receipt},{"reference",J::Numbers(reference)},{"fixed_boundary",J::Constraints(data.fixedBoundary)},
            {"elapsed_seconds",J::Number(elapsed)},{"observation",observation},
            {"resolution",resolution},{"test_basis",testBasis},{"native_constitutive",constitutive}});
        std::ofstream file(argv[3]);Require(bool(file),"cannot write verification output");file<<output<<'\n';
        return physical?0:1;
    }
    catch(const std::exception &error)
    {
        if(argc==4){std::ofstream file(argv[3]);file<<J::Object({{"status",J::Quote("fail")},{"failure",J::Quote(error.what())}})<<'\n';}
        std::cerr<<"preservation verification: "<<error.what()<<std::endl;return 1;
    }
}
