///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <boost/property_tree/json_parser.hpp>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace Nektar;
using namespace Nektar::SolverUtils;
using Vector=RANSWallAlgebra::Vector;
using RANSWallAlgebra::Require;
namespace J=RANSWallSolveJson;
namespace
{
Vector Read(const boost::property_tree::ptree &input,const std::string &key)
{
    Vector values;for(const auto &item:input.get_child(key))
    {const double value=item.second.get_value<double>();Require(std::isfinite(value),"nonfinite input coefficient");values.push_back(value);}
    return values;
}
Vector Binary(const std::string &path,int size)
{
    std::ifstream input(path,std::ios::binary);Require(input.good(),"cannot open input direction");
    Vector values(size);input.read(reinterpret_cast<char*>(values.data()),size*sizeof(double));
    Require(input.gcount()==static_cast<std::streamsize>(size*sizeof(double)) && input.peek()==EOF,
        "input direction byte length differs");
    for(double value:values){Require(std::isfinite(value),"nonfinite direction");}
    return values;
}
std::string Constraints(const std::map<int,double> &values)
{
    std::vector<std::string> rows;for(auto [row,value]:values)rows.push_back(J::Numbers(Vector{double(row),value}));
    return J::List(rows);
}
std::string Field(const RANSWallEvolutionPatch &patch,int block)
{
    std::vector<std::string> elements;const auto field=patch.fields[block];
    const auto map=field->Field()->GetLocalToGlobalMap();
    for(int e=0;e<static_cast<int>(field->Elements().size());++e)
    {
        const auto &element=field->Elements()[e];Vector ids,signs;
        const int start=field->Field()->GetCoeff_Offset(e),count=field->Field()->GetExp(e)->GetNcoeffs();
        for(int j=0;j<count;++j){ids.push_back(map->GetLocalToGlobalMap(start+j));signs.push_back(map->GetLocalToGlobalSign(start+j));}
        elements.push_back(J::Object({{"geometry",J::Number(element.geometry)},{"degree",J::Number(element.degree)},
            {"lower",J::Numbers(element.lower)},{"upper",J::Numbers(element.upper)},
            {"ids",J::Numbers(ids)},{"sign",J::Numbers(signs)}}));
    }
    return J::Object({{"size",J::Number(patch.size[block])},{"mesh",J::List(elements)}});
}
Vector MassAction(const RANSWallEvolutionPatch &patch,const Vector &direction,
    double alpha,const std::map<int,double> &constraints)
{
    Vector result(patch.Size());
    for(int block=0;block<patch.InnerPressure();++block)
    {
        auto part=direction;for(auto [row,value]:constraints)
        {
            (void)value;
            // Retained matching epsilon coefficients are independent. The
            // actual row action is already in J; mass columns at external
            // conditions vanish, while matching columns remain active.
            if(row>=patch.offset[block] && row<patch.offset[block]+patch.size[block] && block!=patch.Epsilon())part[row]=0;
        }
        Vector local(part.begin()+patch.offset[block],part.begin()+patch.offset[block]+patch.size[block]);
        const auto load=RANSWallAlgebra::Multiply(patch.PhysicalMass(block),local);
        for(int i=0;i<patch.size[block];++i)if(!constraints.count(patch.offset[block]+i))
            result[patch.offset[block]+i]=alpha*load[i];
    }
    return result;
}
}

int main(int argc,char **argv)
{
    try
    {
        Require(argc==4 || argc==5,"usage: RANSNonpressureActionProbe SESSION.xml STATE.json REPORT.json [OUTER_K_DIRECTION.bin]");
        Require(!std::ifstream(argv[3]).good(),"refusing to overwrite an action report");
        auto world=LibUtilities::SessionReader::CreateInstance(2,argv);world->InitSession();
        char name[]="nonpressure-action-probe";char *args[]{name,nullptr};
        auto session=LibUtilities::SessionReader::CreateInstance(1,args,std::vector<std::string>{argv[1]},
            LibUtilities::CommSerial::create(1,args));
        session->InitSession();session->SetSolverInfo("GlobalSysSoln","DirectFull");
        session->SetSolverInfo("RANSWallNonlinearPolicy","ScalarFeasibleR3");
        session->SetSolverInfo("RANSWallActionPolicy","VerifiedR3Action");
        const auto graph=SpatialDomains::MeshGraphIO::Read(session);
        Require(world->GetComm()->GetSize()==1 && graph->GetMeshDimension()==2,"R3 action probe is serial flat2D");
        RANSWallPatchParameters parameters;
        session->LoadParameter("RANSWallChi",parameters.chi,1.);
        session->LoadParameter("RANSWallHeight",parameters.geometry.height,2.);
        session->LoadParameter("RANSWallRoughness",parameters.geometry.roughness,.1);
        session->LoadParameter("Kinvis",parameters.geometry.molecularNu,1.5e-5);
        session->LoadParameter("RANSWallSampleHeight",parameters.geometry.sampleHeight,1.);
        session->LoadParameter("RANSWallTailStart",parameters.geometry.tailStart,1.);
        session->LoadParameter("RANSWallTangentialOrder",parameters.tangentialOrder,8);
        session->LoadParameter("RANSWallNormalOrder",parameters.normalOrder,16);
        double us,dx,dy;session->LoadParameter("RANSWallFrictionVelocity",us,.4);
        session->LoadParameter("RANSWallDirectionX",dx,1.);session->LoadParameter("RANSWallDirectionY",dy,0.);
        const double H=parameters.geometry.height,K0=parameters.model.Q()*us*us,E0=us*us*us/(parameters.model.Kappa()*H);
        RANSWallEvolutionPatch patch(session,graph,world->GetComm(),parameters,LoadRANSWallNormalMap(session,2,H));
        Require(!patch.Mapped() && patch.TransitionEnabled() && patch.Size()<=3000,"unsupported R3 action scope");
        auto reference=patch.NeutralReference(us,{dx,dy});RANSWallEvolutionData data(2);
        data.fixedBoundary=patch.FixedReferenceBoundary(reference);reference=patch.Lift(reference,0,data);
        patch.SetReference(reference,data,0,0);
        boost::property_tree::ptree input;boost::property_tree::read_json(argv[2],input);
        const auto state=Read(input,"state"),history=Read(input,"history");
        Require(state.size()==static_cast<std::size_t>(patch.Size()) && history.size()==state.size(),"input field size differs");
        const double time=input.get<double>("time"),alpha=input.get<double>("alpha");
        const auto domain=patch.ScalarDomain(state,time,data,true);
        Require(domain.admissible,"action base is outside common scalar inventory: "+domain.failure);
        const auto constraints=patch.Constraints(state,time,data);
        const auto residual=patch.Residual(state,history,alpha,time,data,false);
        std::vector<std::string> metadata;
        const auto start=std::chrono::steady_clock::now();
        const auto matrix=patch.ScalarFeasibleJacobian(state,history,alpha,time,data,&metadata);
        const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        const std::string matrixFile=std::string(argv[3])+".J.bin";
        Require(!std::ifstream(matrixFile).good(),"refusing to overwrite an action matrix");
        std::ofstream outputMatrix(matrixFile,std::ios::binary);
        for(const auto &row:matrix)outputMatrix.write(reinterpret_cast<const char*>(row.data()),row.size()*sizeof(double));
        outputMatrix.close();Require(outputMatrix.good(),"cannot preserve full action matrix");
        std::vector<std::pair<std::string,Vector>> directions;
        for(const std::string kind:{"velocity","inner_k","sampled_K","epsilon_with_trace","joint"})
        {
            Vector direction(patch.Size());
            if(kind=="inner_k" || kind=="sampled_K")
            {
                const auto values=patch.fields[patch.K()]->Project([&](const RANSWallCoordinate &p)
                {
                    const double x=(p[0]-patch.lower[0])/(patch.upper[0]-patch.lower[0]);
                    if(kind=="inner_k")return p[2]<=.01?4*K0*x*(1-x)*std::pow(1-p[2]/.01,2):0.;
                    const double r=p[2]/H;return p[2]<=H?16*K0*x*(1-x)*r*(1-r):0.;
                });
                std::copy(values.begin(),values.end(),direction.begin()+patch.offset[patch.K()]);
            }
            else for(int block=0;block<patch.InnerPressure();++block)
                if(kind=="joint" || (kind=="velocity" && block<patch.K()) ||
                    (kind=="epsilon_with_trace" && block==patch.Epsilon()))
                    for(int i=0;i<patch.size[block];++i)
                        direction[patch.offset[block]+i]=(block<patch.K()?us:block==patch.K()?K0:E0)*std::sin(1.317*(patch.offset[block]+i+1));
            // Apply the actual external lift to this direction, preserving
            // retained epsilon-H modes. Difference is affine, never a state repair.
            auto moved=state;for(int i=0;i<patch.Size();++i)moved[i]+=direction[i];
            const auto a=patch.ScalarDomain(moved,time,data),b=domain;
            if(!a.effectiveState.empty())for(int i=0;i<patch.Size();++i)direction[i]=a.effectiveState[i]-b.effectiveState[i];
            directions.push_back({kind,std::move(direction)});
        }
        if(argc==5)directions.push_back({"retained_outer_k_witness",Binary(argv[4],patch.Size())});
        std::vector<std::string> actionRecords;
        for(const auto &[name,direction]:directions)
            actionRecords.push_back(J::Object({{"name",J::Quote(name)},{"direction",J::Numbers(direction)},
                {"stage_action",J::Numbers(RANSWallAlgebra::Multiply(matrix,direction))},
                {"exact_mass_action",J::Numbers(MassAction(patch,direction,alpha,constraints))}}));
        const auto polynomial=patch.fields[patch.K()]->Project([](const RANSWallCoordinate &p)
            {return 1+p[0]*p[2]*p[2]+p[2]*p[2]*p[2];});
        double hessianError=0;int hessianPoints=0;
        const auto points=patch.fields[patch.K()]->Quadrature(4);
        for(std::size_t i=0;i<points.size();i+=31)
        {
            const auto &q=points[i];const auto rows=patch.fields[patch.K()]->RadialGradientDerivative(q);
            const double centre=polynomial[q.ids[0]]/q.constant[0];
            for(int axis=0;axis<3;++axis)
            {
                long double value=0;for(std::size_t j=0;j<q.ids.size();++j)value+=rows[axis][j]*(static_cast<long double>(polynomial[q.ids[j]])-centre*q.constant[j]);
                const double expected=axis==0?2*q.coordinate[2]:axis==2?2*q.coordinate[0]+6*q.coordinate[2]:0;
                hessianError=std::max(hessianError,std::abs(static_cast<double>(value)-expected));
            }
            ++hessianPoints;
        }
        std::vector<std::string> fields;for(int block=0;block<patch.Gauge();++block)fields.push_back(Field(patch,block));
        const auto &control=patch.TransitionControls();const auto distances=patch.TransitionDistances();
        const auto report=J::Object({{"schema","1"},{"scope",J::Quote("R3 fixed-state action verification; no trajectory")},
            {"dimension","2"},{"mapped","false"},{"ranks","1"},{"upstream_transition","true"},
            {"chi",J::Number(parameters.chi)},
            {"inner_stress_model",J::Quote(patch.InnerStressName())},{"model_variant",J::Quote(patch.VariantName())},
            {"action_policy",J::Quote(patch.ActionPolicyName())},{"nonlinear_policy",J::Quote(patch.NonlinearPolicyName())},
            {"time",J::Number(time)},{"alpha",J::Number(alpha)},{"eta","0"},{"ustar",J::Number(us)},
            {"parameters",J::Numbers(Vector{parameters.geometry.molecularNu,H,parameters.geometry.roughness,parameters.geometry.sampleHeight,
                parameters.geometry.tailStart,K0,parameters.model.Q(),parameters.model.Kappa()})},
            {"quadrature_orders",J::Numbers(Vector{double(parameters.tangentialOrder),double(parameters.normalOrder)})},
            {"offsets",J::Numbers(Vector(patch.offset.begin(),patch.offset.end()))},{"fields",J::List(fields)},
            {"state",J::Numbers(state)},{"history",J::Numbers(history)},{"initial",J::Numbers(reference)},
            {"constraints",Constraints(constraints)},{"residual",J::Numbers(residual)},
            {"scaled_block_norms",J::Numbers(patch.Norms(residual,constraints,us))},
            {"scalar_domain",RANSWallEvolutionPatch::ScalarDomainJSON(domain,false)},
            {"transition_reference",J::Numbers(patch.TransitionReference())},
            {"transition_controls",J::Numbers(Vector{control.start,control.end,distances[0],distances[1],distances[2],
                control.state.on,control.state.off,control.state.marginOff,control.state.marginOn,control.state.reynoldsOff,control.state.reynoldsOn})},
            {"transition_weights",J::Numbers(control.state.weights)},{"robust_parameters",J::Numbers(patch.RobustParameters())},
            {"robust_reference_flow",J::Numbers(patch.RobustReferenceFlow())},{"robust_reference_normal",J::Numbers(patch.RobustReferenceNormal())},
            {"matrix_file",J::Quote(matrixFile)},{"matrix_size",J::Number(patch.Size())},{"assembly_seconds",J::Number(seconds)},
            {"assembly_metadata",J::List(metadata)},{"directions",J::List(actionRecords)},
            {"radial_hessian_polynomial_error",J::Number(hessianError)},{"radial_hessian_points",J::Number(hessianPoints)}});
        std::ofstream output(argv[3]);output<<report<<'\n';output.close();Require(output.good(),"cannot write action report");
        world->Finalise();return 0;
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 1;}
}
