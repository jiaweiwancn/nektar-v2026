///////////////////////////////////////////////////////////////////////////////
// File: RANSKernelProbe.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/AnalyticalWallFunction.h>
#include <SolverUtils/Turbulence/RANSWallStrip.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Nektar::SolverUtils;
namespace
{
int checks=0,invalidChecks=0;
double maximumTensorSymmetry=0.0,maximumTensorTrace=0.0,maximumFrameError=0.0;
double maximumQuadratureError=0.0,maximumEnergyError=0.0;
void Check(bool condition,const std::string &message)
{
    ++checks; if (!condition) { throw std::runtime_error(message); }
}
void Near(double a,double b,double tolerance=3e-12,double scale=1.0)
{
    Check(std::isfinite(a) && std::isfinite(b) &&
          std::abs(a-b)<=tolerance*std::max({scale,std::abs(a),std::abs(b)}),
          "Kernel equality failed: "+std::to_string(a)+" versus "+std::to_string(b));
}
template<class Function> void Reject(Function function)
{
    bool rejected=false;
    try { function(); } catch (const std::exception &) { rejected=true; }
    ++invalidChecks; Check(rejected,"Invalid kernel input was accepted");
}
template<std::size_t N> void Array(std::ostream &out,const std::array<double,N> &values)
{
    out << '[';
    for (std::size_t i=0;i<N;++i) { if (i) { out << ','; } out << values[i]; }
    out << ']';
}
void Tensor(std::ostream &out,const RANSTensor &matrix)
{
    out << '[';
    for (int i=0;i<3;++i) { if (i) { out << ','; } Array(out,matrix[i]); }
    out << ']';
}
void Parameters(std::ostream &out,const RANSParameters &p)
{
    Array(out,std::array<double,7>{p.cMu,p.cEpsilon1,p.cEpsilon2,p.sigmaK,
                                 p.sigmaEpsilon,p.yapCoefficient,p.yapLength});
}
void Values(std::ostream &out,const AWFValues &v)
{
    out << "{\"traction\":"; Array(out,v.traction);
    out << ",\"production\":" << v.meanProduction << ",\"epsilon_mean\":" << v.meanDissipation
        << ",\"epsilon_P\":" << v.dissipationP << ",\"epsilon_H\":" << v.dissipationH << '}';
}
RANSTensor Rotate(const RANSTensor &matrix)
{
    const double angle=0.73,c=std::cos(angle),s=std::sin(angle);
    const RANSVector axis{1/std::sqrt(14.0),2/std::sqrt(14.0),3/std::sqrt(14.0)};
    RANSTensor rotation{},result{};
    const RANSTensor cross{{{0,-axis[2],axis[1]},{axis[2],0,-axis[0]},{-axis[1],axis[0],0}}};
    for (int i=0;i<3;++i) for (int j=0;j<3;++j)
        { rotation[i][j]=(i==j?c:0)+(1-c)*axis[i]*axis[j]+s*cross[i][j]; }
    for (int i=0;i<3;++i) for (int j=0;j<3;++j)
        for (int k=0;k<3;++k) for (int l=0;l<3;++l)
            { result[i][j]+=rotation[i][k]*matrix[k][l]*rotation[j][l]; }
    return result;
}

struct BulkCase
{
    std::string name;
    RANSModel model;
    RANSTensor gradient;
    double k,epsilon,distance,nu;
    RANSParameters parameters;
    RANSClosureValues result;
};
BulkCase Bulk(std::string name,RANSModel model,RANSTensor g,double k,double epsilon,
              double distance=1.0,double nu=0.0,RANSParameters p={})
{
    auto v=EvaluateRANSClosure(model,k,epsilon,g,distance,p,nu);
    double scale=k,trace=0.0;
    for (const auto &row:v.reynoldsStress) for (double x:row) { scale=std::max(scale,std::abs(x)); }
    for (int i=0;i<3;++i)
    {
        trace+=v.reynoldsStress[i][i];
        for (int j=0;j<3;++j)
        {
            const double error=std::abs(v.reynoldsStress[i][j]-v.reynoldsStress[j][i])/scale;
            maximumTensorSymmetry=std::max(maximumTensorSymmetry,error); Check(error<2e-12,"Nonsymmetric stress");
        }
    }
    if (std::abs(v.divergence)<1e-13)
    {
        const double error=std::abs(trace-2*k)/scale;
        maximumTensorTrace=std::max(maximumTensorTrace,error); Check(error<3e-12,"Incorrect covariance trace");
    }
    const auto rotated=EvaluateRANSClosure(model,k,epsilon,Rotate(g),distance,p,nu);
    const auto expected=Rotate(v.reynoldsStress);
    for (int i=0;i<3;++i) for (int j=0;j<3;++j)
    {
        const double error=std::abs(rotated.reynoldsStress[i][j]-expected[i][j])/scale;
        maximumFrameError=std::max(maximumFrameError,error); Check(error<3e-12,"Stress is not frame covariant");
    }
    Near(rotated.production,v.production,3e-12,std::max(1.0,std::abs(v.production)));
    return {name,model,g,k,epsilon,distance,nu,p,v};
}
void WriteBulk(std::ostream &out,const BulkCase &c)
{
    const auto &v=c.result;
    out << "{\"name\":\"" << c.name << "\",\"model\":\""
        << (c.model==RANSModel::CLSHighRe?"cls":"standard") << "\",\"k\":" << c.k
        << ",\"epsilon\":" << c.epsilon << ",\"distance\":" << c.distance << ",\"nu\":" << c.nu
        << ",\"parameters\":"; Parameters(out,c.parameters);
    out << ",\"gradient\":"; Tensor(out,c.gradient);
    out << ",\"stress\":"; Tensor(out,v.reynoldsStress);
    out << ",\"quadratic\":"; Tensor(out,v.quadraticStress);
    out << ",\"cubic\":"; Tensor(out,v.cubicStress);
    out << ",\"cMu\":" << v.cMu << ",\"nut\":" << v.eddyViscosity
        << ",\"strain\":" << v.strainInvariant << ",\"rotation\":" << v.rotationInvariant
        << ",\"divergence\":" << v.divergence << ",\"production\":" << v.production
        << ",\"k_source\":" << v.kSource << ",\"epsilon_production\":" << v.epsilonProduction
        << ",\"epsilon_destruction\":" << v.epsilonDestruction << ",\"yap\":" << v.yap
        << ",\"epsilon_source\":" << v.epsilonSource << ",\"diffusion\":";
    Array(out,std::array<double,3>{v.momentumDiffusion,v.kDiffusion,v.epsilonDiffusion}); out << '}';
}

// Independent 16-point Gauss integration, split at physical breakpoints and
// geometric subdivisions. No production primitive or moment routine is reused.
double Integrate(const std::function<double(double)> &f,const AnalyticalWallFunction &wall,
                 std::vector<double> extra={})
{
    constexpr double x[8]={.09501250983763744,.2816035507792589,.4580167776572274,.6178762444026438,
                          .755404408355003,.8656312023878318,.9445750230732326,.9894009349916499};
    constexpr double w[8]={.1894506104550685,.1826034150449236,.1691565193950025,.1495959888165767,
                          .1246289712555339,.0951585116824928,.0622535239386479,.027152459411754095};
    const auto &g=wall.Geometry();
    std::vector<double> breaks{0,g.height};
    extra.insert(extra.end(),{g.roughness,wall.ViscosityInterceptHeight(),wall.DissipationCutoff()});
    for (double z:extra) { if (z>0 && z<g.height) { breaks.push_back(z); } }
    for (double z=std::min(wall.DissipationCutoff(),g.molecularNu/
         (wall.Parameters().cMu*wall.Parameters().cLength*std::sqrt(wall.K())));z<g.height;z*=2)
        { breaks.push_back(z); }
    std::sort(breaks.begin(),breaks.end());
    breaks.erase(std::unique(breaks.begin(),breaks.end()),breaks.end());
    long double integral=0;
    for (std::size_t i=1;i<breaks.size();++i)
    {
        const double center=(breaks[i]+breaks[i-1])/2,half=(breaks[i]-breaks[i-1])/2;
        for (int j=0;j<8;++j) { integral+=half*w[j]*(f(center-half*x[j])+f(center+half*x[j])); }
    }
    return static_cast<double>(integral);
}

struct WallCase
{
    std::string name;
    AWFGeometry geometry;
    double k;
    AWFParameters parameters;
    AWFVector velocity,forcing;
};
void WriteWall(std::ostream &out,const WallCase &c,const RANSNeutralReference &model)
{
    const EquilibriumWallFunction ep(c.geometry,c.k,model,c.parameters);
    const auto &wall=ep.Original(); const auto &reference=ep.Reference();
    const auto v=wall.Evaluate(c.velocity,c.forcing);
    const auto j=wall.Integrals();
    const double numericP=Integrate([&](double z){ const auto g=wall.VelocityGradient(z,c.velocity,c.forcing);
        return wall.EddyViscosity(z)*(g[0]*g[0]+g[1]*g[1]); },wall)/c.geometry.height;
    const double error=std::abs(v.meanProduction-numericP)/std::max(1e-100,std::abs(numericP));
    maximumQuadratureError=std::max(maximumQuadratureError,error);
    Check(error<2e-10,"Closed production does not agree with independent quadrature");
    Check(v.meanProduction>=0,"Negative AWF production");
    const auto zero=wall.Velocity(0,c.velocity,c.forcing),match=wall.Velocity(c.geometry.height,c.velocity,c.forcing);
    for (int i=0;i<2;++i) { Check(zero[i]==0,"Auxiliary no slip failed"); Near(match[i],c.velocity[i]); }
    Near(wall.MeanOriginalDissipation(),Integrate([&](double z){return wall.OriginalDissipation(z);},wall)/c.geometry.height,
         3e-12,std::abs(v.meanDissipation));
    const AWFVector otherForce{c.forcing[0]+.07,c.forcing[1]-.03};
    const auto other=wall.Evaluate(c.velocity,otherForce);
    for (int i=0;i<2;++i) { Near(other.traction[i]-v.traction[i],-j.j1/j.j0*(otherForce[i]-c.forcing[i])); }

    if (reference.admissible)
    {
        const AWFVector refU{reference.matchingVelocity,0},noForce{0,0};
        const auto target=ep.EvaluateFixedWeight(refU,noForce,1,AWFFailurePolicy::RequirePreservation);
        Near(target.values.traction[0],reference.traction);
        Near(target.values.meanProduction,target.values.meanDissipation);
        Near(target.values.dissipationP,reference.numerator/c.geometry.pointP);
        const double power=reference.traction*reference.matchingVelocity;
        const double balance=power-c.geometry.height*reference.production-
            reference.traction*reference.traction*c.geometry.molecularNu*j.j2;
        maximumEnergyError=std::max(maximumEnergyError,std::abs(balance)/power);
        Near(balance,0,3e-12,power);
    }

    out << "{\"name\":\"" << c.name << "\",\"q\":" << model.Q()
        << ",\"kappa\":" << model.Kappa() << ",\"model_parameters\":";
    Parameters(out,model.Parameters());
    out << ",\"geometry\":";
    Array(out,std::array<double,4>{c.geometry.molecularNu,c.geometry.height,c.geometry.roughness,c.geometry.pointP});
    out << ",\"k\":" << c.k << ",\"parameters\":";
    Array(out,std::array<double,3>{c.parameters.cMu,c.parameters.cLength,c.parameters.smoothIntercept});
    out << ",\"velocity\":"; Array(out,c.velocity); out << ",\"forcing\":"; Array(out,c.forcing);
    out << ",\"intercept\":" << wall.Intercept() << ",\"zv\":" << wall.ViscosityInterceptHeight()
        << ",\"ze\":" << wall.DissipationCutoff() << ",\"integrals\":";
    Array(out,std::array<double,4>{j.j0,j.j1,j.j2,j.productionWeight});
    out << ",\"original\":"; Values(out,v);
    out << ",\"quadrature_production\":" << numericP
        << ",\"reference\":{\"ustar\":" << reference.frictionVelocity << ",\"traction\":" << reference.traction
        << ",\"U_H\":" << reference.matchingVelocity << ",\"A\":" << reference.numerator
        << ",\"P0\":" << reference.production << ",\"lambda\":" << reference.lambda
        << ",\"minimum_lambda\":" << reference.minimumLambda << ",\"cutoff\":" << reference.cutoff
        << ",\"log_z0\":" << reference.logAerodynamicRoughness << ",\"admissible\":"
        << (reference.admissible?"true":"false") << "},\"blends\":[";
    int count=0;
    for (double weight:{0.0,.25,.7,1.0})
    {
        const auto result=ep.EvaluateFixedWeight(c.velocity,c.forcing,weight);
        Check(result.values.dissipationP==ep.Dissipation(c.geometry.pointP,result.weight) &&
              result.values.dissipationH==ep.Dissipation(c.geometry.height,result.weight),
              "Point and trace dissipation do not use the same auxiliary profile");
        Check(result.values.traction==v.traction && result.values.meanProduction==v.meanProduction,
              "Equilibrium correction changed original momentum or production");
        if (weight==0 || !reference.admissible)
        {
            Check(result.values.meanDissipation==v.meanDissipation && result.values.dissipationP==v.dissipationP &&
                  result.values.dissipationH==v.dissipationH && result.weight==0,"Baseline recovery is not exact");
        }
        const double mean=Integrate([&](double z){ return ep.Dissipation(z,result.weight); },wall,
                                    reference.admissible?std::vector<double>{reference.cutoff}:std::vector<double>{})/c.geometry.height;
        Near(mean,result.values.meanDissipation,5e-12,std::abs(mean));
        double previous=std::numeric_limits<double>::infinity();
        for (int i=0;i<=32;++i)
        {
            const double value=ep.Dissipation(c.geometry.height*i/32.0,result.weight);
            Check(value>0 && value<=previous*(1+1e-14),"Dissipation profile lost positivity or monotonicity");
            previous=value;
        }
        if (count++) { out << ','; }
        out << "{\"requested\":" << weight << ",\"weight\":" << result.weight
            << ",\"reason\":" << static_cast<int>(result.reason) << ",\"values\":";
        Values(out,result.values); out << ",\"quadrature_mean\":" << mean << '}';
    }
    out << "],\"samples\":[";
    count=0;
    std::vector<double> heights{0,c.geometry.pointP,c.geometry.height};
    if (reference.admissible) { heights.push_back(reference.cutoff); }
    for (double z:heights)
    {
        if (count++) { out << ','; }
        out << "{\"z\":" << z << ",\"nu_t\":" << wall.EddyViscosity(z) << ",\"D\":" << wall.TotalViscosity(z)
            << ",\"epsilon_A\":" << wall.OriginalDissipation(z) << ",\"U\":";
        Array(out,wall.Velocity(z,c.velocity,c.forcing)); out << ",\"gradient\":";
        Array(out,wall.VelocityGradient(z,c.velocity,c.forcing)); out << '}';
    }
    out << "]}";
}

void ReferenceTests(const RANSNeutralReference &model)
{
    const double u=.4,H=2,U=5.7;
    for (double z:{.1,1.0,2.0,10.0,100.0})
    {
        auto r=model.AtHeight(z,u,H,U);
        RANSTensor gradient{}; gradient[0][2]=r.velocityGradient;
        const auto value=EvaluateRANSClosure(model.Model(),r.k,r.epsilon,gradient,z,model.Parameters());
        Near(value.eddyViscosity,r.eddyViscosity); Near(value.production,r.epsilon);
        Near(value.epsilonSource+std::pow(u,4)/(model.Parameters().sigmaEpsilon*z*z),0);
        Near(r.epsilonGradient,-r.epsilon/z);
        Near(r.epsilonGradientFlux,r.eddyViscosity/model.Parameters().sigmaEpsilon*r.epsilonGradient);
        Near(value.reynoldsStress[0][2],-u*u);
        if (model.Model()==RANSModel::CLSHighRe)
        {
            Near(value.reynoldsStress[0][0]/r.k,.92);
            Near(value.reynoldsStress[1][1]/r.k,.56);
            Near(value.reynoldsStress[2][2]/r.k,.52);
        }
    }
    Reject([&]{model.AtHeight(0,u,H,U);});
}
void SwitchTests(const EquilibriumWallFunction &ep,std::ostream &out,
                 const RANSNeutralReference &model=RANSNeutralReference::CLS())
{
    const auto &r=ep.Reference(); const auto &wall=ep.Original();
    AWFSwitchParameters parameters{.05,.2,{1,1,1,1,1,1,1},{1,0}};
    const auto atP=model.AtHeight(wall.Geometry().pointP,r.frictionVelocity,wall.Geometry().height,r.matchingVelocity);
    const auto atN=model.AtHeight(4,r.frictionVelocity,wall.Geometry().height,r.matchingVelocity);
    const AWFIndependentSamples base{{r.matchingVelocity,0},{atP.velocity,0},{atN.velocity,0},4,wall.K(),atN.epsilon,0,{0,0},{0,0}};
    const auto equilibrium=ep.Evaluate(base,parameters);
    Check(equilibrium.weight==1,"Switch is not on its exact equilibrium plateau");
    double previous=1;
    for (int i=0;i<=100;++i)
    {
        double value=AWFEquilibriumWeight(i*.003,parameters.on,parameters.off);
        Check(value>=0 && value<=previous && value<=1,"Switch is not monotone in [0,1]"); previous=value;
    }
    for (double point:{parameters.on,parameters.off})
    {
        const double h=1e-6,a=AWFEquilibriumWeight(point-h,parameters.on,parameters.off),
                     b=AWFEquilibriumWeight(point,parameters.on,parameters.off),
                     c=AWFEquilibriumWeight(point+h,parameters.on,parameters.off);
        Check(std::abs((c-a)/(2*h))<2e-6 && std::abs((c-2*b+a)/(h*h))<.05,"Switch endpoint is not C2");
    }
    out << '[';
    for (int index=0;index<9;++index)
    {
        auto sample=base;
        if (index==1) { sample.matchingVelocity[0]+=1; }
        if (index==2) { sample.velocityN[0]+=1; }
        if (index==3) { sample.kN*=2; }
        if (index==4) { sample.epsilonN*=2; }
        if (index==5) { sample.pressureForcing[0]=1; }
        if (index==6) { sample.convectionForcing[0]=1; }
        if (index==7) { sample.normalVelocityP=1; }
        if (index==8) { sample.pressureForcing[0]=1; sample.convectionForcing[0]=-1; }
        const auto value=ep.Evaluate(sample,parameters);
        if (index>0) { Check(value.weight==0,"Independent departure did not recover the baseline"); }
        if (index) { out << ','; }
        out << "{\"case\":" << index << ",\"weight\":" << value.weight << ",\"departure\":" << value.departure
            << ",\"sensors\":"; Array(out,value.sensors); out << '}';
    }
    out << ']';
    auto rotated=base; const double c=std::cos(.9),s=std::sin(.9);
    auto rotate=[&](AWFVector value){ return AWFVector{c*value[0]-s*value[1],s*value[0]+c*value[1]}; };
    rotated.matchingVelocity=rotate(base.matchingVelocity); rotated.velocityP=rotate(base.velocityP);
    rotated.velocityN=rotate(base.velocityN); parameters.incomingDirection={c,s};
    Check(ep.Evaluate(rotated,parameters).weight==1,"Tangential frame rotation broke equilibrium sensing");
    const auto baseline=wall.Evaluate(base.matchingVelocity,{0,0});
    const auto disabled=ep.EvaluateFixedWeight(base.matchingVelocity,{0,0},1,AWFFailurePolicy::UseOriginal,false);
    Check(disabled.weight==0 && disabled.reason==AWFFallbackReason::OutsideHighReRegime &&
          disabled.values.meanDissipation==baseline.meanDissipation,"Wall-regime fallback failed");
    Reject([&]{ep.EvaluateFixedWeight(base.matchingVelocity,{0,0},1,AWFFailurePolicy::RequirePreservation,false);});
}
} // namespace

int main(int argc,char *argv[])
{
    try
    {
        std::string path="rans-kernels.json";
        if (argc==3 && std::string(argv[1])=="--output") { path=argv[2]; }
        else if (argc!=1) { throw std::invalid_argument("Use --output PATH"); }
        const auto standard=RANSNeutralReference::Standard(),cls=RANSNeutralReference::CLS();
        Near(cls.Q(),3.190675081624288); Near(cls.Kappa(),.440553905186659);
        Near(standard.Kappa(),std::sqrt(.48*1.3*.3));
        ReferenceTests(standard); ReferenceTests(cls);
        const auto shifted=standard.StandardShifted(0,.4,.02);
        Check(shifted.velocity==0 && shifted.epsilon>0,"Shifted standard profile failed");
        Reject([&]{cls.StandardShifted(1,.4,.02);});
        const double k=standard.Q()*.16,up=.4/standard.Kappa()*std::log(2/.02);
        const auto traction=standard.StandardLogWallTraction(k,{up,-up},2,.02);
        Near(traction[0],.16); Near(traction[1],-.16);
        Near(RANSTotalDissipation(1,.1,1.5e-5,{2,-3,.5}),.1+1.5e-5*13.25/2);
        Check(RANSTotalDissipation(1,.1,0,{2,-3,.5})==.1,"Zero-nu dissipation conversion changed epsilon");
        const auto coefficients=CLSCubicCoefficients(.1);
        for (int i=0;i<7;++i) { Near(coefficients[i],std::array<double,7>{-.1,.1,.26,-.1,0,-.05,.05}[i]); }
        Check(CLSYapFunction(0)==0 && CLSYapFunction(1)==0,"Yap lower branch failed");
        Near(CLSYapFunction(2),.83*4);

        std::vector<BulkCase> bulk;
        const std::vector<RANSTensor> gradients{
            {}, {{{0,0,1},{0,0,0},{0,0,0}}}, {{{0,-2,0},{2,0,0},{0,0,0}}},
            {{{-1,0,0},{0,.5,0},{0,0,.5}}}, {{{1,.3,-.7},{-.6,-.4,.9},{.4,-.2,-.6}}}};
        for (int i=0;i<gradients.size();++i)
            for (auto model:{RANSModel::StandardKEpsilon,RANSModel::CLSHighRe})
                { bulk.push_back(Bulk("canonical_"+std::to_string(i),model,gradients[i],1,i==3?.01:.7)); }
        Check(bulk[7].result.production<0,"Negative CLS production was silently clipped");
        std::mt19937_64 generator(20261029);
        std::uniform_real_distribution<double> random(-2,2);
        for (int sample=0;sample<32;++sample)
        {
            RANSTensor g{};
            for (auto &row:g) for (auto &v:row) { v=random(generator); }
            g[2][2]=-g[0][0]-g[1][1];
            for (auto model:{RANSModel::StandardKEpsilon,RANSModel::CLSHighRe})
                { bulk.push_back(Bulk("general_"+std::to_string(sample),model,g,.5+sample*.07,.03+sample*.05,.2+sample*.1,1.5e-5)); }
        }
        RANSTensor divergence{{{1,0,0},{0,1,0},{0,0,1}}};
        RANSParameters custom; custom.cMu=.075; custom.cEpsilon1=1.4;custom.cEpsilon2=1.9;
        custom.sigmaK=.85;custom.sigmaEpsilon=1.17;custom.yapCoefficient=.71;custom.yapLength=2.8;
        for (auto model:{RANSModel::StandardKEpsilon,RANSModel::CLSHighRe})
        {
            bulk.push_back(Bulk("raw_divergence",model,divergence,.6,.4));
            Near(bulk.back().result.divergence,3);
            bulk.push_back(Bulk("custom_coefficients",model,gradients[4],.8,.2,1.3,2e-5,custom));
        }
        for (double invariant:{0.0,1.0,8.0,1e3,1e8,1e100})
            { Check(std::isfinite(CLSCmu(invariant)) && CLSCmu(invariant)>0,"Unstable CLS coefficient"); }

        const double unit=1.5e-5/std::sqrt(.5),smooth=10.7*unit;
        std::vector<WallCase> walls;
        struct GeometryCase { std::string name; AWFGeometry geometry; double k; AWFParameters p; };
        const std::vector<GeometryCase> geometries{
            {"note",{1.5e-5,2,.1,1},cls.Q()*.16,{}},
            {"smooth",{1.5e-5,40*unit,0,20*unit},.5,{}},
            {"rough_below_viscous",{1.5e-5,40*unit,unit,20*unit},.5,{}},
            {"rough_above_viscous",{1.5e-5,40*unit,10*unit,20*unit},.5,{}},
            {"fully_rough",{1.5e-5,1000*unit,100*unit,500*unit},.5,{}},
            {"rough_threshold",{1.5e-5,200*unit,70*unit,100*unit},.5,{}},
            {"molecular_inner_epsilon",{1.5e-5,.25*smooth,0,.125*smooth},.5,{}},
            {"molecular_tail_epsilon",{1.5e-5,.75*smooth,0,.375*smooth},.5,{}},
            {"weak_turbulent",{1.5e-5,1.0001*smooth,0,.50005*smooth},.5,{}},
            {"forcing_thin_layer",{1.5e-5,2,2*(1-1e-8),1},.5,{}},
            {"custom_wall",{2e-5,1.7,.05,.8},.7,{.08,2.8,11.0}},
            {"current_K_changed",{1.5e-5,2,.1,1},4*cls.Q()*.16,{}}};
        const CLSEquilibriumWallFunction note(geometries[0].geometry,geometries[0].k,cls);
        for (const auto &g:geometries)
        {
            const double U=g.name=="note" ? note.Reference().matchingVelocity : 2.3;
            for (const auto &f:std::vector<AWFVector>{{-.02,0},{0,0},{.03,0},{10,-3}})
                { walls.push_back({g.name,g.geometry,g.k,g.p,{U,g.name=="note"?0.0:.2*U},f}); }
            walls.push_back({g.name+"_reversed",g.geometry,g.k,g.p,{-U,.3*U},{0,.07}});
        }
        auto edge=geometries[0].geometry; edge.pointP=note.Reference().cutoff;
        const CLSEquilibriumWallFunction equal(edge,geometries[0].k,cls);
        Check(equal.Reference().admissible,"Cutoff equality was rejected");
        edge.pointP=std::nextafter(edge.pointP,0.0);
        const CLSEquilibriumWallFunction fail(edge,geometries[0].k,cls);
        Check(!fail.Reference().admissible,"Failed cutoff was silently clipped");
        Reject([&]{fail.EvaluateFixedWeight({1,0},{0,0},1,AWFFailurePolicy::RequirePreservation);});
        Check(fail.EvaluateFixedWeight({1,0},{0,0},1).weight==0,"Inadmissible reference did not fall back");

        const double nan=std::numeric_limits<double>::quiet_NaN();
        Reject([&]{CLSCmu(-1);}); Reject([&]{CLSCmu(nan);});
        Reject([&]{EvaluateRANSClosure(RANSModel::CLSHighRe,0,1,{},1);});
        Reject([&]{EvaluateRANSClosure(RANSModel::CLSHighRe,1,-1,{},1);});
        Reject([&]{EvaluateRANSClosure(RANSModel::CLSHighRe,1,1,{},0);});
        Reject([&]{EvaluateRANSClosure(RANSModel::StandardKEpsilon,1,1,{},1,{},-1);});
        Reject([&]{auto g=gradients[1];g[2][0]=nan;EvaluateRANSClosure(RANSModel::CLSHighRe,1,1,g,1);});
        Reject([&]{auto p=RANSParameters{};p.sigmaK=0;EvaluateRANSClosure(RANSModel::CLSHighRe,1,1,{},1,p);});
        Reject([&]{CLSCubicCoefficients(0);}); Reject([&]{CLSYapFunction(-1);});
        Reject([&]{auto p=RANSParameters{};p.cEpsilon2=p.cEpsilon1;RANSNeutralReference::CLS(p);});
        Reject([&]{AnalyticalWallFunction({0,1,0,.5},1);});
        Reject([&]{AnalyticalWallFunction({1e-5,1,1,.5},1);});
        Reject([&]{AnalyticalWallFunction({1e-5,1,0,1},1);});
        Reject([&]{AnalyticalWallFunction({1e-5,1,0,.5},-1);});
        Reject([&]{note.Original().Velocity(-1,{1,0},{0,0});});
        Reject([&]{note.EvaluateFixedWeight({1,0},{0,0},1.01);});
        Reject([&]{EquilibriumWallFunction({1.5e-5,2,.1,2},1,standard);});
        Reject([&]{AWFEquilibriumWeight(.1,.2,.1);});

        std::ofstream out(path); out << std::setprecision(17);
        out << "{\"schema\":1,\"bulk\":[\n";
        for (int i=0;i<bulk.size();++i) { if (i) { out << ",\n"; } WriteBulk(out,bulk[i]); }
        out << "],\"walls\":[\n";
        for (int i=0;i<walls.size();++i) { if (i) { out << ",\n"; } WriteWall(out,walls[i],cls); }
        out << "],\"switch_cases\":"; SwitchTests(note,out);
        // The standard closure has a different Q and Kappa, including when
        // its configured bulk constants differ from the wall-kernel constants.
        // Independent Gauss integration in WriteWall verifies both families.
        out << ",\"standard_walls\":[";
        int standardCases=0;
        for (const auto &model : {standard,RANSNeutralReference::Standard(custom)})
        {
            ReferenceTests(model);
            for (const auto &g : geometries)
            {
                const double energy=g.name=="note" ? model.Q()*.16 : g.k;
                const EquilibriumWallFunction ep(g.geometry,energy,model,g.p);
                for (const auto &forcing : {AWFVector{0,0},AWFVector{.03,-.02}})
                {
                    if (standardCases++) out << ',';
                    WriteWall(out,{g.name,g.geometry,energy,g.p,
                                   {ep.Reference().matchingVelocity,0},forcing},model);
                }
            }
            const EquilibriumWallFunction ep(geometries[0].geometry,model.Q()*.16,model);
            auto boundary=geometries[0].geometry;
            boundary.pointP=ep.Reference().cutoff;
            Check(EquilibriumWallFunction(boundary,model.Q()*.16,model).Reference().admissible,
                  "Standard cutoff equality rejected");
            boundary.pointP=std::nextafter(boundary.pointP,0.0);
            const EquilibriumWallFunction invalid(boundary,model.Q()*.16,model);
            Check(!invalid.Reference().admissible,"Standard cutoff silently clipped");
            Reject([&]{invalid.EvaluateFixedWeight({1,0},{0,0},1,AWFFailurePolicy::RequirePreservation);});
        }
        out << "],\"standard_switch_cases\":";
        SwitchTests(EquilibriumWallFunction(geometries[0].geometry,standard.Q()*.16,standard),out,standard);
        out << ",\"configured_standard_switch_cases\":";
        const auto configuredStandard=RANSNeutralReference::Standard(custom);
        SwitchTests(EquilibriumWallFunction(geometries[0].geometry,configuredStandard.Q()*.16,configuredStandard),
                    out,configuredStandard);
        out << ",\"references\":[";
        RANSParameters noYap;noYap.yapCoefficient=0;
        const std::vector<RANSNeutralReference> references{standard,cls,RANSNeutralReference::Standard(custom),
            RANSNeutralReference::CLS(custom),RANSNeutralReference::CLS(noYap)};
        for (int i=0;i<references.size();++i)
        {
            if (i) { out << ','; } const auto &r=references[i];
            out << "{\"model\":\"" << (r.Model()==RANSModel::CLSHighRe?"cls":"standard")
                << "\",\"q\":" << r.Q() << ",\"kappa\":" << r.Kappa() << ",\"cMu\":" << r.Cmu()
                << ",\"parameters\":"; Parameters(out,r.Parameters()); out << '}';
        }
        out << "],\"summary\":{\"native_checks\":" << checks << ",\"rejected_inputs\":" << invalidChecks
            << ",\"bulk_cases\":" << bulk.size() << ",\"wall_cases\":" << walls.size()
            << ",\"standard_wall_cases\":" << standardCases
            << ",\"maximum_tensor_symmetry_error\":" << maximumTensorSymmetry
            << ",\"maximum_tensor_trace_error\":" << maximumTensorTrace
            << ",\"maximum_frame_error\":" << maximumFrameError
            << ",\"maximum_production_quadrature_relative_error\":" << maximumQuadratureError
            << ",\"maximum_reference_energy_relative_error\":" << maximumEnergyError << "}}\n";
        out.close();
        if (!out) { throw std::runtime_error("Cannot write kernel report"); }
        std::cout << "RANS kernel checks passed: " << checks << " checks, " << bulk.size()
                  << " tensor cases, " << walls.size() << " wall cases\n";
        return 0;
    }
    catch (const std::exception &error)
    {
        std::cerr << "RANS kernels: " << error.what() << '\n'; return 1;
    }
}
