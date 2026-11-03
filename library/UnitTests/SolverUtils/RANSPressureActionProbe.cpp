///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>

using namespace Nektar;
using namespace Nektar::SolverUtils;
using Vector=RANSWallAlgebra::Vector;
using Vec=RANSVector;
using RANSWallAlgebra::Require;
namespace J=Nektar::RANSWallSolveJson;
namespace
{
std::string Constraints(const std::map<int,double> &values)
{
    std::vector<std::string> rows;
    for(auto [row,value]:values)rows.push_back(J::Numbers(Vector{double(row),value}));
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
        elements.push_back(J::Object({{"geometry",J::Number(element.geometry)},
            {"degree",J::Number(element.degree)},{"lower",J::Numbers(element.lower)},
            {"upper",J::Numbers(element.upper)},{"ids",J::Numbers(ids)},
            {"sign",J::Numbers(signs)}}));
    }
    return J::Object({{"size",J::Number(patch.size[block])},{"mesh",J::List(elements)}});
}
double Step(double x)
{
    if(x<=0)return 0;
    if(x>=1)return 1;
    return x*x*(3-2*x);
}
double Depletion(const Vec &p,double H,double rK,double width)
{
    const double horizontal=Step((p[0]+4)/.1)*Step((-3-p[0])/.1);
    const double distance=std::abs(p[2]-rK)/width;
    const double peak=distance>=1?0:1-Step(distance);
    const double restore=Step((p[2]-(H-width))/width);
    return horizontal*(1-std::max(peak,restore));
}
Vector MapParameters(const RANSWallNormalMap &map)
{
    return {double(map.kind),double(map.dimension),map.height,map.heightSlopeX,map.heightSlopeY,
        map.heightAmplitude,map.heightWaveX,map.heightWaveY,map.slopeX,map.slopeY,map.radius,
        map.amplitude,map.lengthX,map.lengthY};
}
std::string Samples(const RANSWallEvolutionPatch &patch,const Vector &state,double time,
    const RANSWallEvolutionData &data)
{
    std::vector<std::string> records;
    for(const auto &face:patch.Surface(0,2,patch.Parameters().geometry.height,1,8))
    {
        const auto point=face.test.coordinate;const auto input=patch.TransitionInput(state,point,time,data);
        const auto transition=patch.Transition(state,point,time,data);
        const auto profile=patch.StateProfile(state,point,time,data);
        const auto normalization=patch.WallNormalization(state,point,*profile);
        records.push_back(J::Object({{"chart",J::Numbers(point)},
            {"pressure_sensor",J::Numbers(input.sensors.pressure)},
            {"gates",J::Numbers(Vector{transition.chi,transition.spatial,transition.state,transition.admissibility,transition.d2})},
            {"K",J::Number(profile->K())},{"b",J::Number(normalization.b)},
            {"m",J::Number(normalization.m)},{"lambda",J::Number(normalization.lambda)},
            {"pre_trace",J::Number(profile->Dissipation(patch.Distance(point)))},
            {"actual_trace",J::Number(patch.DependentEpsilon(state,point,time,data))}}));
    }
    return J::List(records);
}
} // namespace

int main(int argc,char **argv)
{
    try
    {
        Require(argc==3,"usage: RANSPressureActionProbe SESSION.xml REPORT.json");
        auto world=LibUtilities::SessionReader::CreateInstance(2,argv);world->InitSession();
        char name[]="pressure-action-probe";char *args[]{name,nullptr};
        std::vector<std::string> filenames{argv[1]};
        auto session=LibUtilities::SessionReader::CreateInstance(1,args,filenames,LibUtilities::CommSerial::create(1,args));
        session->InitSession();session->SetSolverInfo("GlobalSysSoln","DirectFull");
        const auto graph=SpatialDomains::MeshGraphIO::Read(session);const int dimension=graph->GetMeshDimension();
        Require(dimension==2 && world->GetComm()->GetSize()==1,"bounded pressure fixture is serial 2D");
        RANSWallPatchParameters parameters;
        session->LoadParameter("RANSWallChi",parameters.chi,1.);
        session->LoadParameter("RANSWallHeight",parameters.geometry.height,2.);
        session->LoadParameter("RANSWallRoughness",parameters.geometry.roughness,.1);
        session->LoadParameter("Kinvis",parameters.geometry.molecularNu,1.5e-5);
        session->LoadParameter("RANSWallSampleHeight",parameters.geometry.sampleHeight,1.);
        session->LoadParameter("RANSWallTailStart",parameters.geometry.tailStart,1.);
        session->LoadParameter("RANSWallTangentialOrder",parameters.tangentialOrder,8);
        session->LoadParameter("RANSWallNormalOrder",parameters.normalOrder,16);
        const auto mapping=LoadRANSWallNormalMap(session,dimension,parameters.geometry.height);
        RANSWallEvolutionPatch patch(session,graph,world->GetComm(),parameters,mapping);
        Require(!patch.Mapped() && patch.TransitionEnabled(),"pressure fixture requires flat Upstream mode");
        int active=0;session->LoadParameter("RANSWallPressureProbeActive",active,0);
        double dose=.06,depth=.99,width=.03;
        session->LoadParameter("RANSWallPressureProbeDose",dose,.06);
        session->LoadParameter("RANSWallPressureProbeDepth",depth,.99);
        session->LoadParameter("RANSWallPressureProbeWidth",width,.03);
        const double us=.4,K0=parameters.model.Q()*us*us,H=parameters.geometry.height,time=.06,alpha=300.;
        const auto reference=patch.NeutralReference(us,{1.,0.});auto state=reference;
        RANSWallEvolutionData data(dimension);data.fixedBoundary=patch.FixedReferenceBoundary(reference);
        if(active)
        {
            const auto k=patch.fields[patch.K()]->Project([&](const Vec &point)
            {return point[2]>H?K0:K0*(1-depth*Depletion(point,H,parameters.geometry.sampleHeight,width));});
            std::copy(k.begin(),k.end(),state.begin()+patch.offset[patch.K()]);
            const double gradient=dose*(K0/parameters.model.Q()+std::pow(parameters.geometry.molecularNu/H,2))/H;
            const auto pressure=patch.fields[patch.OuterPressure()]->Project([&](const Vec &point)
            {return gradient*(point[0]-(patch.lower[0]+patch.upper[0])/2);});
            std::copy(pressure.begin(),pressure.end(),state.begin()+patch.offset[patch.OuterPressure()]);
        }
        state=patch.Lift(state,time,data);auto history=state;
        for(int i=0;i<patch.DifferentialSize();++i)history[i]*=.97;
        const auto constraints=patch.Constraints(state,time,data);
        const auto base=patch.Residual(state,history,alpha,time,data,false);
        std::vector<std::string> fields;
        for(int block=0;block<patch.Gauge();++block)fields.push_back(Field(patch,block));
        std::vector<std::string> directions;
        for(const std::string kind:{"inner_pressure","outer_pressure","gauge","joint_pressure"})
        {
            Vector direction(patch.Size());
            for(int i=patch.DifferentialSize();i<patch.Size();++i)
            {
                const bool inner=i<patch.offset[patch.OuterPressure()];
                const bool gauge=i==patch.offset[patch.Gauge()];
                if(kind=="joint_pressure" || (kind=="inner_pressure" && inner) ||
                    (kind=="outer_pressure" && !inner && !gauge) || (kind=="gauge" && gauge))
                    direction[i]=(gauge?us/H:us*us)*std::sin(1.317*(i+1));
            }
            const auto affine=patch.AffinePressureAction(direction,constraints);
            const auto complete=patch.PressureDirectionalAction(state,direction,time,data);
            std::vector<std::string> levels;
            for(double h:{2e-4,5e-5,1e-5,2e-6,5e-7})
            {
                auto plus=state,minus=state;
                for(int i=0;i<patch.Size();++i){plus[i]+=h*direction[i];minus[i]-=h*direction[i];}
                const auto a=patch.Residual(plus,history,alpha,time,data,false);
                const auto b=patch.Residual(minus,history,alpha,time,data,false);
                levels.push_back(J::Object({{"step",J::Number(h)},
                    {"plus_residual",J::Numbers(a)},{"minus_residual",J::Numbers(b)},
                    {"plus_constraints",Constraints(patch.Constraints(plus,time,data))},
                    {"minus_constraints",Constraints(patch.Constraints(minus,time,data))},
                    {"plus_samples",Samples(patch,plus,time,data)},{"minus_samples",Samples(patch,minus,time,data)}}));
            }
            directions.push_back(J::Object({{"kind",J::Quote(kind)},
                {"direction",J::Numbers(direction)},{"affine_action",J::Numbers(affine)},
                {"complete_action",J::Numbers(complete)},{"levels",J::List(levels)}}));
        }
        std::vector<std::string> affineColumns;
        if(!active)for(int column=patch.DifferentialSize();column<patch.Size();++column)
        {
            Vector direction(patch.Size());direction[column]=1.;
            affineColumns.push_back(J::Object({{"column",J::Number(column)},
                {"action",J::Numbers(patch.AffinePressureAction(direction,constraints))}}));
        }
        const auto &control=patch.TransitionControls();const auto distances=patch.TransitionDistances();
        const auto report=J::Object({{"schema","1"},{"fixture",J::Quote(active?"pressure_sensor_and_capped_trace_active":"affine_r1_geometry")},
            {"operator_verification_only","true"},{"physical_qualification","false"},
            {"dimension",J::Number(dimension)},{"ranks","1"},{"chi",J::Number(parameters.chi)},
            {"eta","0"},{"upstream_transition","true"},{"mapped","false"},
            {"time",J::Number(time)},{"alpha",J::Number(alpha)},{"ustar",J::Number(us)},
            {"parameters",J::Numbers(Vector{parameters.geometry.molecularNu,H,parameters.geometry.roughness,
                parameters.geometry.sampleHeight,parameters.geometry.tailStart,K0,parameters.model.Q(),parameters.model.Kappa()})},
            {"robust_parameters",J::Numbers(patch.RobustParameters())},{"model_variant",J::Quote(patch.VariantName())},
            {"inner_stress_model",J::Quote(patch.InnerStressName())},
            {"robust_reference_flow",J::Numbers(patch.RobustReferenceFlow())},{"robust_reference_normal",J::Numbers(patch.RobustReferenceNormal())},
            {"map_parameters",J::Numbers(MapParameters(mapping))},{"quadrature_orders",J::Numbers(Vector{double(parameters.tangentialOrder),double(parameters.normalOrder)})},
            {"transition_reference",J::Numbers(patch.TransitionReference())},
            {"transition_controls",J::Numbers(Vector{control.start,control.end,distances[0],distances[1],distances[2],
                control.state.on,control.state.off,control.state.marginOff,control.state.marginOn,control.state.reynoldsOff,control.state.reynoldsOn})},
            {"transition_weights",J::Numbers(control.state.weights)},{"offsets",J::Numbers(patch.offset)},
            {"fields",J::List(fields)},{"state",J::Numbers(state)},{"history",J::Numbers(history)},
            {"constraints",Constraints(constraints)},{"base_residual",J::Numbers(base)},
            {"samples",Samples(patch,state,time,data)},{"directions",J::List(directions)},
            {"affine_columns",J::List(affineColumns)},
            {"data",J::Object({{"dose",J::Number(dose)},{"depth",J::Number(depth)},{"width",J::Number(width)},
                {"body_and_scalar_sources",J::Quote("zero")},{"initial_velocity",J::Quote("actual stored native neutral reference")}})}});
        std::ofstream output(argv[2]);Require(output.good(),"cannot open pressure report");output<<report<<'\n';
        Require(output.good(),"cannot finish pressure report");world->Finalise();return 0;
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 1;}
}
