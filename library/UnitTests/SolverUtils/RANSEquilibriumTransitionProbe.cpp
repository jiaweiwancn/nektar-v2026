///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/EquilibriumAWFTransition.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Nektar::SolverUtils;
namespace
{
using V = std::vector<double>;
int checks = 0;
void Check(bool ok, const std::string &message)
{
    ++checks;
    if (!ok) throw std::runtime_error(message);
}
void Near(double a, double b, double tolerance = 2e-12)
{
    Check(std::isfinite(a) && std::isfinite(b) &&
        std::abs(a-b) <= tolerance*std::max({1.0, std::abs(a), std::abs(b)}), "Transition native identity");
}
void Array(std::ostream &out, const V &v)
{
    out << '[';
    for (std::size_t j = 0; j < v.size(); ++j)
    {
        if (j) out << ',';
        out << v[j];
    }
    out << ']';
}
std::string Quote(const std::string &text)
{
    std::string result = "\"";
    for (char c : text)
    {
        if (c == '\\' || c == '"') result += '\\';
        if (c == '\n') result += "\\n";
        else result += c;
    }
    return result+'"';
}
V Flat(const AWFTransitionSensors &s)
{
    return {s.velocity[0],s.velocity[1],s.velocity[2],s.velocityDifference[0],s.velocityDifference[1],
        s.velocityDifference[2],s.normal,s.k,s.epsilon,s.pressure[0],s.pressure[1],s.pressure[2],
        s.convection[0],s.convection[1],s.convection[2]};
}
AWFTransitionSensors Sensors(const V &v)
{
    if (v.size() != 15) throw std::logic_error("Sensor vector shape");
    return {{v[0],v[1],v[2]},{v[3],v[4],v[5]},{v[9],v[10],v[11]},
        {v[12],v[13],v[14]},v[6],v[7],v[8]};
}
V Flat(const AWFTransitionInput &i)
{
    V result{i.coordinate,i.budgetMargin,i.reynolds};
    const auto s = Flat(i.sensors);
    result.insert(result.end(),s.begin(),s.end());
    return result;
}
AWFTransitionInput Input(const V &v)
{
    if (v.size() != 18) throw std::logic_error("Transition input vector shape");
    return {v[0],v[1],v[2],Sensors(V(v.begin()+3,v.end()))};
}
V Flat(const AWFTransitionControls &c)
{
    V result{c.start,c.end,c.state.on,c.state.off,c.state.marginOff,c.state.marginOn,
        c.state.reynoldsOff,c.state.reynoldsOn};
    result.insert(result.end(),c.state.weights.begin(),c.state.weights.end());
    return result;
}
AWFTransitionControls Controls(const V &v)
{
    if (v.size() != 15) throw std::logic_error("Transition controls vector shape");
    AWFTransitionControls c;
    c.start=v[0]; c.end=v[1]; c.state.on=v[2]; c.state.off=v[3];
    c.state.marginOff=v[4]; c.state.marginOn=v[5]; c.state.reynoldsOff=v[6]; c.state.reynoldsOn=v[7];
    std::copy(v.begin()+8,v.end(),c.state.weights.begin());
    return c;
}
AWFTransitionControlDirection ControlDirection(const V &v)
{
    AWFTransitionControlDirection c;
    c.start=v[0]; c.end=v[1]; c.on=v[2]; c.off=v[3]; c.marginOff=v[4]; c.marginOn=v[5];
    c.reynoldsOff=v[6]; c.reynoldsOn=v[7];
    std::copy(v.begin()+8,v.end(),c.weights.begin());
    return c;
}
V Flat(const AWFTransitionValues &v)
{
    return {v.chi,v.spatial,v.state,v.admissibility,v.marginGate,v.reynoldsGate,v.d2};
}
V Flat(const AWFTransitionProfileValues &v) { return {v.dissipation,v.redistribution,v.kSource}; }
const std::array<RANSVector AWFTransitionPhysicalInput::*,11> physicalVectors{
    &AWFTransitionPhysicalInput::normal,&AWFTransitionPhysicalInput::matchingVelocity,
    &AWFTransitionPhysicalInput::matchingReference,&AWFTransitionPhysicalInput::velocityA,
    &AWFTransitionPhysicalInput::velocityAReference,&AWFTransitionPhysicalInput::velocityB,
    &AWFTransitionPhysicalInput::velocityBReference,&AWFTransitionPhysicalInput::pressure,
    &AWFTransitionPhysicalInput::pressureReference,&AWFTransitionPhysicalInput::convection,
    &AWFTransitionPhysicalInput::convectionReference};
const std::array<double AWFTransitionPhysicalInput::*,7> physicalScalars{
    &AWFTransitionPhysicalInput::sampledK,&AWFTransitionPhysicalInput::height,
    &AWFTransitionPhysicalInput::molecularNu,&AWFTransitionPhysicalInput::Q,
    &AWFTransitionPhysicalInput::kA,&AWFTransitionPhysicalInput::epsilonA,
    &AWFTransitionPhysicalInput::epsilonAReference};
V Flat(const AWFTransitionPhysicalInput &p)
{
    V result;
    for (auto member:physicalVectors)
        for (double value:p.*member) result.push_back(value);
    for (auto member:physicalScalars) result.push_back(p.*member);
    return result;
}
AWFTransitionPhysicalInput Physical(const V &v)
{
    if (v.size()!=40) throw std::logic_error("Physical sensor vector shape");
    AWFTransitionPhysicalInput p;
    std::size_t j=0;
    for (auto member:physicalVectors)
        for (double &value:p.*member) value=v[j++];
    for (auto member:physicalScalars) p.*member=v[j++];
    return p;
}
V Map(const RANSWallNormalMap &m)
{
    return {static_cast<double>(m.kind),static_cast<double>(m.dimension),m.height,m.heightSlopeX,
        m.heightSlopeY,m.heightAmplitude,m.heightWaveX,m.heightWaveY,m.slopeX,m.slopeY,
        m.radius,m.amplitude,m.lengthX,m.lengthY};
}
void Transition(std::ostream &out,const std::string &name,const AWFTransitionInput &i,
    const AWFTransitionControls &c,const AWFTransitionInput &d={},const V &dc=V(15,0))
{
    const auto result=EvaluateEquilibriumAWFTransition(i,c,d,ControlDirection(dc));
    const auto &s=i.sensors;
    const auto norm=[](const RANSVector &a){return std::hypot(a[0],a[1],a[2]);};
    const double legacy=EquilibriumAWFProfileWeight(result.value.spatial,
        {norm(s.velocity),norm(s.velocityDifference),s.normal,s.k,s.epsilon,norm(s.pressure),norm(s.convection)},
        i.budgetMargin,i.reynolds,c.state);
    Near(result.value.chi,legacy);
    Near(result.value.chi,result.value.spatial*result.value.state*result.value.admissibility);
    Check(result.value.chi>=0 && result.value.chi<=1,"Transition range");
    out << "{\"name\":" << Quote(name) << ",\"input\":"; Array(out,Flat(i));
    out << ",\"controls\":"; Array(out,Flat(c));
    out << ",\"direction\":"; Array(out,Flat(d));
    out << ",\"control_direction\":"; Array(out,dc);
    out << ",\"value\":"; Array(out,Flat(result.value));
    out << ",\"first\":"; Array(out,Flat(result.derivative));
    out << ",\"second\":"; Array(out,Flat(result.secondDerivative));
    out << ",\"legacy_weight\":" << legacy << '}';
}
void Coordinate(std::ostream &out,const std::string &name,const RANSWallNormalMap &m,
    const RANSVector &q,const RANSVector &d,bool affine=false)
{
    const RANSVector origin{.1,-.2,.3}, covector{1,.3,m.slopeX+.3*m.slopeY};
    const auto v=affine?EquilibriumAWFAffineWallCoordinate(m,q,origin,covector,d):
        EquilibriumAWFWallCoordinateDifferential(m,q,d);
    auto ray=q; ray[2]+=3.1;
    const auto same=affine?EquilibriumAWFAffineWallCoordinate(m,ray,origin,covector,d):
        EquilibriumAWFWallCoordinateDifferential(m,ray,d);
    Near(v.coordinate,same.coordinate); Near(v.derivative,same.derivative); Near(v.chartGradient[2],0);
    out << "{\"name\":" << Quote(name) << ",\"map_parameters\":"; Array(out,Map(m));
    out << ",\"chart\":"; Array(out,{q[0],q[1],q[2]});
    out << ",\"direction\":"; Array(out,{d[0],d[1],d[2]});
    out << ",\"affine\":" << (affine?"true":"false") << ",\"origin\":";
    Array(out,{origin[0],origin[1],origin[2]});
    out << ",\"covector\":"; Array(out,{covector[0],covector[1],covector[2]});
    out << ",\"coordinate\":" << v.coordinate << ",\"first\":" << v.derivative
        << ",\"second\":" << v.secondDerivative << ",\"chart_gradient\":";
    Array(out,{v.chartGradient[0],v.chartGradient[1],v.chartGradient[2]});
    out << ",\"ray_coordinate\":" << same.coordinate << '}';
}
struct ProfilePoint
{
    double K=0,cutoff=0,margin=0,reynolds=0,productionBudget=0,minimumBudget=0;
    AWFTransitionInput input;
    AWFTransitionResult transition;
    std::vector<AWFTransitionProfileInput> profiles;
    std::vector<AWFTransitionProfileValues> values;
    double flux0=0,fluxH=0,redistributionIntegral=0;
};
ProfilePoint Profile(double K,double K0,const AWFProfileGeometry &geometry,
    const RANSNeutralReference &model,const AWFTransitionControls &controls,const V &heights)
{
    AWFProfileOptions options; options.quadratureOrder=32; options.requestedWeight=0;
    EquilibriumAWFProfiles baseline(geometry,K,model,{},options);
    const auto &b=baseline.Reference();
    ProfilePoint point;
    point.K=K; point.margin=b.margin; point.reynolds=b.reynolds;
    point.productionBudget=b.productionBudget; point.minimumBudget=b.minimumBudget;
    point.input={-.5*(std::abs(controls.start)+std::abs(controls.end)),b.margin,b.reynolds,{}};
    point.input.sensors.velocity={.045+.02*(K/K0-1),.012,0};
    point.transition=EvaluateEquilibriumAWFTransition(point.input,controls);
    options.requestedWeight=1;
    EquilibriumAWFProfiles corrected(geometry,K,model,{},options);
    point.cutoff=corrected.Reference().cutoff;
    options.requestedWeight=point.transition.value.chi;
    EquilibriumAWFProfiles blend(geometry,K,model,{},options);
    for (double r:heights)
    {
        const double production=.27+.03*r+.02*K;
        AWFTransitionProfileInput p{baseline.Dissipation(r),corrected.Dissipation(r),
            corrected.ReferenceProduction(r),production};
        const auto v=BlendEquilibriumAWFTransition(point.transition,p);
        const auto actual=blend.AtProduction(r,production);
        Near(v.value.dissipation,actual.dissipation); Near(v.value.redistribution,actual.redistribution);
        Near(v.value.kSource,actual.kSource);
        point.profiles.push_back(p); point.values.push_back(v.value);
    }
    long double net=0;
    for (const auto &q:blend.Quadrature()) net+=q.weight*q.jacobian*blend.AtProduction(q.r,0).redistribution;
    point.redistributionIntegral=static_cast<double>(net);
    Near(point.redistributionIntegral,0,2e-10);
    point.flux0=blend.RedistributionFlux(0); point.fluxH=blend.RedistributionFlux(geometry.height);
    Near(point.flux0,0); Near(point.fluxH,0,2e-10);
    return point;
}
void WriteProfile(std::ostream &out,const ProfilePoint &p)
{
    out << "{\"K\":" << p.K << ",\"cutoff\":" << p.cutoff << ",\"margin\":" << p.margin
        << ",\"reynolds\":" << p.reynolds << ",\"production_budget\":" << p.productionBudget
        << ",\"minimum_budget\":" << p.minimumBudget << ",\"input\":";
    Array(out,Flat(p.input)); out << ",\"transition\":"; Array(out,Flat(p.transition.value));
    out << ",\"samples\":[";
    for (std::size_t j=0;j<p.values.size();++j)
    {
        if(j)out << ',';
        const auto &a=p.profiles[j];
        out << "{\"profile\":"; Array(out,{a.originalDissipation,a.correctedDissipation,a.referenceProduction,a.production});
        out << ",\"value\":"; Array(out,Flat(p.values[j])); out << '}';
    }
    out << "],\"flux_endpoints\":"; Array(out,{p.flux0,p.fluxH});
    out << ",\"redistribution_integral\":" << p.redistributionIntegral << '}';
}
} // namespace

int main(int argc,char **argv)
{
    std::string path;
    for(int j=1;j<argc;++j)
    {
        const std::string option=argv[j];
        if(option=="--output" && j+1<argc)path=argv[++j];
        else {std::cerr << "Usage: RANSEquilibriumTransitionProbe --output REPORT.json\n";return 2;}
    }
    if(path.empty()){std::cerr << "Missing --output\n";return 2;}
    std::ofstream out(path);
    if(!out){std::cerr << "Cannot write transition report\n";return 2;}
    out << std::setprecision(17);
    try
    {
        AWFTransitionControls controls;
        controls.state.weights={1,1.1,1.2,1.3,1.4,1.5,1.6};
        AWFTransitionInput middle{-1.5,.075,150,{}};
        middle.sensors={{.04,.01,0},{.02,.01,0},{.01,0,.012},{.012,.008,0},.015,.01,.02};
        const AWFTransitionInput plateau{-3,.2,300,{}};
        out << "{\"schema\":1,\"transition_value_order\":[\"chi\",\"spatial\",\"state\",\"admissibility\",\"margin_gate\",\"reynolds_gate\",\"d2\"],\"transition_cases\":[";
        Transition(out,"equal_reference",plateau,controls);
        out << ','; Transition(out,"interior",middle,controls);
        for(int j=0;j<18;++j)
        {
            V direction(18,0); direction[j]=1;
            out << ','; Transition(out,"input_direction_"+std::to_string(j),middle,controls,Input(direction));
        }
        for(int j=0;j<15;++j)
        {
            V direction(15,0); direction[j]=1;
            out << ','; Transition(out,"control_direction_"+std::to_string(j),middle,controls,{},direction);
        }
        for(int j=0;j<7;++j)
        {
            AWFTransitionInput active=plateau;
            V sensor(15,0);
            const std::array<int,7> indices{1,4,6,7,8,10,13};
            sensor[indices[j]]=2*controls.state.off/std::sqrt(controls.state.weights[j]);
            active.sensors=Sensors(sensor);
            Near(EvaluateEquilibriumAWFTransition(active,controls).value.chi,0);
            out << ','; Transition(out,"independent_sensor_"+std::to_string(j),active,controls);
        }
        auto cancel=plateau;
        cancel.sensors.pressure={.2,0,0};cancel.sensors.convection={-.2,0,0};
        Near(EvaluateEquilibriumAWFTransition(cancel,controls).value.chi,0);
        out << ',';Transition(out,"pressure_convection_no_cancellation",cancel,controls);
        for(int gate=0;gate<4;++gate)for(int endpoint=0;endpoint<2;++endpoint)
        {
            auto i=plateau; AWFTransitionInput d;
            if(gate==0){i.coordinate=endpoint?controls.end:controls.start;d.coordinate=1;}
            if(gate==1){i.sensors.velocity[0]=endpoint?controls.state.off:controls.state.on;d.sensors.velocity[0]=1;}
            if(gate==2){i.budgetMargin=endpoint?controls.state.marginOn:controls.state.marginOff;d.budgetMargin=1;}
            if(gate==3){i.reynolds=endpoint?controls.state.reynoldsOn:controls.state.reynoldsOff;d.reynolds=1;}
            const auto r=EvaluateEquilibriumAWFTransition(i,controls,d);
            Near(r.derivative.chi,0);Near(r.secondDerivative.chi,0);
            out << ',';Transition(out,"endpoint_"+std::to_string(gate)+"_"+std::to_string(endpoint),i,controls,d);
        }
        for(int gate=0;gate<4;++gate)for(int endpoint=0;endpoint<2;++endpoint)
        {
            int level=0;
            for(double fraction:{.01,.005,.0025})
            {
                auto i=plateau;AWFTransitionInput d;
                const double q=endpoint?1-fraction:fraction;
                if(gate==0){i.coordinate=controls.start+q*(controls.end-controls.start);d.coordinate=1;}
                if(gate==1)
                {
                    i.sensors.velocity[0]=std::sqrt(controls.state.on*controls.state.on+
                        q*(controls.state.off*controls.state.off-controls.state.on*controls.state.on));
                    d.sensors.velocity[0]=1;
                }
                if(gate==2){i.budgetMargin=controls.state.marginOff+q*(controls.state.marginOn-controls.state.marginOff);d.budgetMargin=1;}
                if(gate==3){i.reynolds=controls.state.reynoldsOff+q*(controls.state.reynoldsOn-controls.state.reynoldsOff);d.reynolds=1;}
                out << ',';Transition(out,"approach_"+std::to_string(gate)+"_"+
                    std::to_string(endpoint)+"_"+std::to_string(level++),i,controls,d);
            }
        }
        auto inactive=plateau;inactive.budgetMargin=-1;inactive.reynolds=0;
        out << ',';Transition(out,"inadmissible_inactive",inactive,controls);
        for(double coordinate:{0.,1.,100.})
        {
            auto i=plateau;i.coordinate=coordinate;
            Near(EvaluateEquilibriumAWFTransition(i,controls).value.chi,0);
            out << ',';Transition(out,"downstream_"+std::to_string(coordinate),i,controls);
        }
        out << "],\"physical_cases\":[";
        AWFTransitionPhysicalInput physical;
        physical.normal={0,0,1};physical.matchingVelocity={1.03,.13,.01};physical.matchingReference={1,.1,0};
        physical.velocityA={1.4,.11,.02};physical.velocityAReference={1.38,.1,0};
        physical.velocityB={1.6,.08,-.01};physical.velocityBReference={1.56,.1,0};
        physical.pressure={.015,-.01,0};physical.pressureReference={.01,-.012,0};
        physical.convection={-.012,.006,0};physical.convectionReference={-.01,.005,0};
        physical.sampledK=.4;physical.height=2;physical.molecularNu=1.5e-5;physical.Q=3.4;
        physical.kA=.41;physical.epsilonA=.18;physical.epsilonAReference=.17;
        for(int j=-1;j<40;++j)
        {
            V direction(40,0);
            if(j>=0)direction[j]=1;
            const auto derivative=Physical(direction);
            const auto r=EvaluateEquilibriumAWFTransitionSensors(physical,derivative);
            if(j!=-1)out << ',';
            out << "{\"name\":" << Quote(j<0?"physical_base":"physical_direction_"+std::to_string(j)) << ",\"input\":";
            Array(out,Flat(physical));out << ",\"direction\":";Array(out,direction);
            out << ",\"value\":";Array(out,Flat(r.value));
            out << ",\"first\":";Array(out,Flat(r.derivative));
            out << ",\"second\":";Array(out,Flat(r.secondDerivative));out << '}';
        }
        auto equal=physical;
        equal.matchingVelocity=equal.matchingReference;equal.velocityA=equal.velocityAReference;
        equal.velocityB=equal.velocityBReference;equal.pressure=equal.pressureReference;
        equal.convection=equal.convectionReference;equal.kA=equal.sampledK;equal.epsilonA=equal.epsilonAReference;
        const auto equalSensors=EvaluateEquilibriumAWFTransitionSensors(equal);
        for(double value:Flat(equalSensors.value))Near(value,0);
        out << "],\"equal_physical_input\":";Array(out,Flat(equal));
        out << ",\"equal_physical_sensors\":";Array(out,Flat(equalSensors.value));
        out << ",\"coordinates\":[";
        RANSWallNormalMap wall;
        Coordinate(out,"flat",wall,{-2.1,0,.3},{.7,0,1});
        wall.kind=RANSWallNormalMap::Kind::Plane;wall.slopeX=.3;
        out << ',';Coordinate(out,"plane2d",wall,{1.3,0,.8},{.5,0,.7});
        wall.dimension=3;wall.slopeY=-.2;
        out << ',';Coordinate(out,"plane3d",wall,{1.1,.4,.8},{.5,.2,.7});
        out << ',';Coordinate(out,"affine3d",wall,{1.1,.4,.8},{.5,.2,.7},true);
        wall={};wall.kind=RANSWallNormalMap::Kind::CosineHill;wall.amplitude=.2;wall.lengthX=2;
        for(double x:{-4.,-2.,-.7,0.,.9,2.,4.})
        {
            out << ',';Coordinate(out,"cosine_"+std::to_string(x),wall,{x,0,.4},{.8,0,1});
        }
        AWFTransitionControls upstream;upstream.start=-5;upstream.end=-3;
        ValidateEquilibriumAWFTransitionSupport(wall,upstream);
        out << "],\"blends\":[";
        const auto transition=EvaluateEquilibriumAWFTransition(middle,controls,
            { .2,.01,10,{{.1,-.03,.02},{-.02,.01,0},{.01,.02,0},{-.03,.01,0},.02,.01,-.02}});
        const AWFTransitionProfileInput p{.5,.8,.3,.6},dp{.1,-.2,.07,.09},ddp{-.03,.04,.01,-.02};
        const auto blend=BlendEquilibriumAWFTransition(transition,p,dp,ddp);
        out << "{\"transition_value\":";Array(out,Flat(transition.value));
        out << ",\"transition_first\":";Array(out,Flat(transition.derivative));
        out << ",\"transition_second\":";Array(out,Flat(transition.secondDerivative));
        out << ",\"profile\":";Array(out,{p.originalDissipation,p.correctedDissipation,p.referenceProduction,p.production});
        out << ",\"first_profile\":";Array(out,{dp.originalDissipation,dp.correctedDissipation,dp.referenceProduction,dp.production});
        out << ",\"second_profile\":";Array(out,{ddp.originalDissipation,ddp.correctedDissipation,ddp.referenceProduction,ddp.production});
        out << ",\"value\":";Array(out,Flat(blend.value));out << ",\"first\":";Array(out,Flat(blend.derivative));
        out << ",\"second\":";Array(out,Flat(blend.secondDerivative));out << "}]";
        const double nan=std::numeric_limits<double>::quiet_NaN();
        const auto safe=BlendEquilibriumAWFTransition(EvaluateEquilibriumAWFTransition(inactive,controls),{.5,nan,nan,.6});
        Near(safe.value.dissipation,.5);Near(safe.value.redistribution,0);Near(safe.value.kSource,.1);
        out << ",\"safe_inactive\":{\"corrected_profile_evaluated\":false,\"value\":";
        Array(out,Flat(safe.value));out << ",\"first\":";Array(out,Flat(safe.derivative));
        out << ",\"second\":";Array(out,Flat(safe.secondDerivative));out << '}';
        const auto model=RANSNeutralReference::CLS();
        const double K0=model.Q()*.16;
        const AWFProfileGeometry geometry{1.5e-5,2,.1,.5,1.5};
        AWFProfileOptions budgetOptions;budgetOptions.requestedWeight=0;budgetOptions.quadratureOrder=32;
        const EquilibriumAWFProfiles budget(geometry,K0,model,{},budgetOptions);
        AWFTransitionControls pc;pc.state.marginOff=.5*budget.Reference().margin;
        pc.state.marginOn=1.5*budget.Reference().margin;pc.state.reynoldsOff=.8*budget.Reference().reynolds;
        pc.state.reynoldsOn=1.2*budget.Reference().reynolds;
        auto activeOptions=budgetOptions;activeOptions.requestedWeight=1;
        const EquilibriumAWFProfiles active(geometry,K0,model,{},activeOptions);
        const V heights{.25*active.Reference().cutoff,.7,geometry.height};
        out << ",\"profiles\":{\"geometry\":";Array(out,{geometry.molecularNu,geometry.height,geometry.roughness,geometry.sampleHeight,geometry.tailStart});
        out << ",\"K0\":" << K0 << ",\"Q\":" << model.Q() << ",\"kappa\":" << model.Kappa()
            << ",\"controls\":";Array(out,Flat(pc));out << ",\"heights\":";Array(out,heights);
        out << ",\"base\":";WriteProfile(out,Profile(K0,K0,geometry,model,pc,heights));
        out << ",\"levels\":[";
        int level=0;
        for(double fraction:{.001,.0005,.00025})
        {
            if(level++)out << ',';
            const double h=K0*fraction;
            out << "{\"h\":" << h << ",\"minus\":";WriteProfile(out,Profile(K0-h,K0,geometry,model,pc,heights));
            out << ",\"plus\":";WriteProfile(out,Profile(K0+h,K0,geometry,model,pc,heights));out << '}';
        }
        out << "]},\"rejections\":[";
        int rejected=0;
        const auto reject=[&](const std::string &name,const auto &function)
        {
            std::string error;
            try{function();}catch(const std::exception &e){error=e.what();}
            Check(!error.empty(),"Invalid transition input accepted: "+name);
            if(rejected++)out << ',';
            out << "{\"name\":" << Quote(name) << ",\"message\":" << Quote(error) << '}';
        };
        for(int index:{0,2,4,6})
        {
            auto c=Flat(controls);c[index]=c[index+1];
            reject("unordered_"+std::to_string(index),[&](){EvaluateEquilibriumAWFTransition(middle,Controls(c));});
        }
        for(int j=2;j<15;++j)
        {
            auto c=Flat(controls);c[j]=0;
            reject("nonpositive_control_"+std::to_string(j),[&](){EvaluateEquilibriumAWFTransition(middle,Controls(c));});
        }
        auto invalidInput=middle;invalidInput.reynolds=-1;
        reject("negative_reynolds",[&](){EvaluateEquilibriumAWFTransition(invalidInput,controls);});
        invalidInput=middle;invalidInput.sensors.epsilon=nan;
        reject("nonfinite_sensor",[&](){EvaluateEquilibriumAWFTransition(invalidInput,controls);});
        auto invalidControls=controls;invalidControls.end=nan;
        reject("nonfinite_interval",[&](){EvaluateEquilibriumAWFTransition(middle,invalidControls);});
        for(int j=33;j<40;++j)
        {
            auto values=Flat(physical);values[j]=0;
            reject("nonpositive_physical_"+std::to_string(j),[&](){EvaluateEquilibriumAWFTransitionSensors(Physical(values));});
        }
        auto badNormal=physical;badNormal.normal={0,0,2};
        reject("nonunit_normal",[&](){EvaluateEquilibriumAWFTransitionSensors(badNormal);});
        auto overlap=upstream;overlap.end=-1;
        reject("transition_overlaps_hill",[&](){ValidateEquilibriumAWFTransitionSupport(wall,overlap);});
        auto noFlat=wall;noFlat.kind=RANSWallNormalMap::Kind::Gaussian;
        reject("no_flat_gaussian",[&](){ValidateEquilibriumAWFTransitionSupport(noFlat,upstream);});
        noFlat.kind=RANSWallNormalMap::Kind::Cylinder;
        reject("no_flat_cylinder",[&](){ValidateEquilibriumAWFTransitionSupport(noFlat,upstream);});
        RANSWallNormalMap plane;plane.dimension=3;
        reject("affine_normal_covector",[&](){EquilibriumAWFAffineWallCoordinate(plane,{0,0,0},{0,0,0},{0,0,1});});
        reject("affine_zero_covector",[&](){EquilibriumAWFAffineWallCoordinate(plane,{0,0,0},{0,0,0},{0,0,0});});
        reject("active_undefined_profile",[&](){BlendEquilibriumAWFTransition(EvaluateEquilibriumAWFTransition(plateau,controls),{.5,nan,nan,.6});});
        out << "],\"checks\":" << checks << ",\"rejected\":" << rejected << ",\"success\":true}\n";
        std::cout << "EP-AWF transition kernel checks=" << checks << " rejected=" << rejected << '\n';
        return 0;
    }
    catch(const std::exception &error)
    {
        std::cerr << "EP-AWF transition probe: " << error.what() << '\n';
        return 1;
    }
}
