// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallPhysicsSignature.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallManufactured.h>
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
#include "RANSTensorValidation.h"
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
    {"constant",J::Numbers(p.constant)},{"derivative",J::Rows(p.derivative)},
    {"geometry",J::Number(p.geometry)},{"endpoint_roundoff_axes",J::Number(p.endpointRoundoffAxes)},
    {"raw_local_coordinate",J::Numbers(p.rawLocalCoordinate)}});}

// Static variable-height geometry/trace qualification. These polynomial
// fields are probes, not a claimed equilibrium or an evolving solution.
std::string StaticCT6Transmission(RANSWallEvolutionPatch &patch)
{
    Require(patch.StaticTransmissionOnly(),"static transmission requires CT6StaticTransmissionV1");
    const auto &p=patch.Parameters();const double H0=p.geometry.height,K0=p.model.Q()*.16;
    std::vector<RANSWallEvolutionFunction> functions{
        [](const auto &x,double){return .15*x[2]+.01*x[0];},
        [](const auto &x,double){return .002*x[0]*x[0];},
        [=](const auto &x,double){return K0*(1+.02*x[0]+.03*x[2]/H0);},
        [=](const auto &,double){return std::pow(K0/p.model.Q(),1.5)/(p.model.Kappa()*H0);},
        [](const auto &x,double){return .02*x[0]+.01*x[2];},
        [](const auto &x,double){return .03*x[0]-.01*x[2];}};
    const auto state=patch.Project(functions,0.);
    const auto zeroNormal=patch.fields[patch.K()]->Project([=](const auto &x){return K0*(1+.02*x[0]);});
    RANSWallEvolutionData data(2);std::vector<std::string> rows;
    double maximum=0,tiltedSignal=0,omittedTiltSignal=0,fixedChartSignal=0,minimumCovariance=1e100;
    const auto dot=[](const RANSVector &a,const RANSVector &b)
    {double x=0;for(int j=0;j<3;++j)x+=a[j]*b[j];return x;};
    for(const auto &face:patch.Surface(0,2,H0,1,12))
    {
        const auto q=face.test.coordinate;const auto g=patch.Mapping().At(q);
        const auto &m=patch.Mapping();
        const double Hp=m.heightSlopeX+m.heightAmplitude*m.heightWaveX*std::cos(m.heightWaveX*q[0]);
        const double height=H0+m.heightAmplitude*std::sin(m.heightWaveX*q[0]);
        const double graph=m.amplitude*std::exp(-q[0]*q[0]/(m.lengthX*m.lengthX));
        const double graphDerivative=-2*q[0]*graph/(m.lengthX*m.lengthX),a=std::hypot(1.,graphDerivative);
        const double gamma=-(4*q[0]*q[0]/std::pow(m.lengthX,4)-2/(m.lengthX*m.lengthX))*graph/(a*a*a);
        const double J=1+height*gamma;
        maximum=std::max({maximum,std::abs(g.height-height),std::abs(g.baseArea-a),std::abs(g.J-J)});
        const double qH=Hp/(g.baseArea*g.J);
        RANSVector N{},area=patch.Mapping().MatchingArea(q);
        for(int j=0;j<3;++j)N[j]=(g.normal[j]-qH*g.tangent1[j])/std::sqrt(1+qH*qH);
        const double areaNorm=std::hypot(area[0],area[1],area[2]);
        for(int j=0;j<3;++j)maximum=std::max(maximum,std::abs(N[j]-area[j]/areaNorm));
        const auto profile=patch.StateProfile(state,q,0.,data);const auto &pg=profile->Geometry();
        Require(pg.height==g.height && pg.sampleHeight==g.height && pg.tailStart==g.height,
            "Accepted local endpoint metadata mismatch");
        const double sampled=face.test.Value(state,patch.offset[patch.K()]);
        const double dependent=patch.DependentEpsilon(state,q,0.,data);
        const double expected=std::pow(sampled/p.model.Q(),1.5)/(p.model.Kappa()*g.height);
        maximum=std::max(maximum,std::abs(dependent-expected));
        const auto forcing=patch.ForcingParts(state,q,0.,data);
        // p_o=.03q-.01rho and k=K0(1+.02q+.03rho/H0).
        // Differentiate at fixed PHYSICAL r, then take the outer trace.
        const double gp=.03+.01*H0*Hp/height,gk=.02*K0-.03*K0*Hp/height;
        const double expectedMechanical=(gp-2*gk/3)/(a*J);
        const double wrongChart=(.03-2*(.02*K0)/3)/(a*J);
        const double mechanicalDefect=std::abs(dot(forcing[0],g.tangent1)-expectedMechanical);
        maximum=std::max(maximum,mechanicalDefect);
        fixedChartSignal=std::max(fixedChartSignal,std::abs(wrongChart-expectedMechanical));
        RANSVector u{face.test.Value(state,patch.offset[0]),0,face.test.Value(state,patch.offset[1])};
        const double mH=g.J*(dot(u,g.normal)-dot(u,g.tangent1)*qH);
        maximum=std::max(maximum,std::abs(mH-dot(u,area)/g.baseArea));
        std::vector<std::string> regions;
        for(int side:{-1,1})
        {
            const bool outer=side>0;const auto v=patch.Point(0,q,side);
            const auto kg=patch.Point(patch.K(),q,side).Gradient(state,patch.offset[patch.K()]);
            const double k=v.Value(state,patch.offset[patch.K()]);
            RANSTensor G{},T{},R{};G[0]=v.Gradient(state,patch.offset[0]);G[2]=v.Gradient(state,patch.offset[1]);
            const double pi=patch.Point(patch.InnerPressure()+outer,q,side).Value(state,patch.offset[patch.InnerPressure()+outer]);
            double D=0,production=0,radicand=0;bool radicandDefined=false;
            if(outer)
            {
                const double eps=patch.Point(patch.Epsilon(),q,side).Value(state,patch.offset[patch.Epsilon()]);
                const auto law=patch.OuterClosure(k,eps,G,q);D=law.kDiffusion;R=law.reynoldsStress;production=law.production;
                radicand=law.a95Radicand;radicandDefined=law.a95RadicandDefined;
                for(int i=0;i<3;++i)for(int j=0;j<3;++j)T[i][j]=-law.reynoldsStress[i][j]+(i==j?2*k/3-pi:0);
            }
            else
            {
                const double nt=profile->Original().EddyViscosity(g.height);D=p.geometry.molecularNu+nt/p.model.Parameters().sigmaK;
                const auto law=patch.InnerStress(k,nt,G);R=law.reynoldsStress;production=law.production;
                T=RANSWallInnerMeanStress(law,p.geometry.molecularNu,pi);
            }
            const auto covariance=J::Eigenvalues(R);minimumCovariance=std::min(minimumCovariance,covariance[0]);
            Require(J::CovarianceAdmissible(covariance[0],K0),"CT6 static sampled covariance failure");
            const double dr=dot(kg,g.normal),dt=dot(kg,g.tangent1);
            const double normalFlux=mH*k-g.J*D*(dr-qH*dt);
            const double physicalFlux=dot(u,area)*k/g.baseArea-D*dot(kg,area)/g.baseArea;
            maximum=std::max(maximum,std::abs(normalFlux-physicalFlux));
            const double chain=g.baseArea*g.J*dt+Hp*dr;
            maximum=std::max(maximum,std::abs(chain-.02*K0));
            RANSVector a{},b{};
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            {a[i]+=T[i][j]*area[j]/g.baseArea;b[i]+=g.J*T[i][j]*(g.normal[j]-qH*g.tangent1[j]);}
            for(int i=0;i<3;++i)maximum=std::max(maximum,std::abs(a[i]-b[i]));
            const auto zn=patch.Point(patch.K(),q,side).Gradient(zeroNormal,0);
            const double tilt=-D*dot(zn,area)/g.baseArea;
            // A zero base-normal derivative does not annul tilted diffusion.
            maximum=std::max(maximum,std::abs(dot(zn,g.normal)));
            tiltedSignal=std::max(tiltedSignal,std::abs(tilt));
            omittedTiltSignal=std::max(omittedTiltSignal,std::abs(g.J*D*qH*dt));
            regions.push_back(J::Object({{"side",J::Number(side)},{"geometry",J::Number(v.geometry)},
                {"k_gradient",J::Numbers(kg)},{"normal_flux",J::Number(normalFlux)},
                {"Cartesian_flux",J::Number(physicalFlux)},{"traction_per_ground_area",J::Numbers(a)},
                {"Reynolds_stress",J::Rows(R)},{"covariance_eigenvalues",J::Numbers(covariance)},
                {"raw_divergence",J::Number(G[0][0]+G[2][2])},{"signed_production",J::Number(production)},
                {"A95_radicand",radicandDefined?J::Number(radicand):"null"},
                {"zero_base_normal_gradient_tilted_diffusion",J::Number(tilt)}}));
        }
        rows.push_back(J::Object({{"chart",J::Numbers(q)},{"physical",J::Numbers(g.point)},
            {"H",J::Number(g.height)},{"qH",J::Number(qH)},{"J",J::Number(g.J)},
            {"normal",J::Numbers(N)},{"cofactor",J::Numbers(area)},{"mH",J::Number(mH)},
            {"mechanical_sampler_defect",J::Number(mechanicalDefect)},
            {"fixed_chart_wrong_sampler_difference",J::Number(wrongChart-expectedMechanical)},
            {"dependent_epsilon",J::Number(dependent)},{"regions",J::List(regions)}}));
    }
    const auto ground=patch.ContinuousGroundDiagnostic(state);
    Tree trace;std::istringstream traceInput(ground);boost::property_tree::read_json(traceInput,trace);
    Require(trace.get<bool>("admissible"),"CT6 ground chart polynomial check rejected the positive static polynomial");
    // The independently prescribed ground polynomial is K0*(1+.02*q).
    // Its minimum is the left endpoint regardless of varying curvature/H.
    const double exactMinimum=K0*(1+.02*patch.lower[0]);
    maximum=std::max(maximum,std::abs(trace.get<double>("minimum")-exactMinimum));
    auto negativeFunctions=functions;
    negativeFunctions[patch.K()]=[=](const auto &x,double)
        {return K0*((x[0]+.23)*(x[0]+.23)-1e-6);};
    const auto negativeState=patch.Project(negativeFunctions,0.);
    const auto negativeGround=patch.ContinuousGroundDiagnostic(negativeState);
    Tree negativeTrace;std::istringstream negativeInput(negativeGround);boost::property_tree::read_json(negativeInput,negativeTrace);
    maximum=std::max(maximum,std::abs(negativeTrace.get<double>("minimum")+K0*1e-6));
    Require(!negativeTrace.get<bool>("admissible"),"CT6 continuous ground check missed a narrow interior negative interval");
    Require(maximum<=2e-10 && tiltedSignal>1e-8 && omittedTiltSignal>1e-8 && fixedChartSignal>1e-8,"CT6 static transmission identity or discriminating signal failed");
    return J::Object({{"schema",J::Quote("CT6StaticTransmissionV1")},{"status",J::Quote("pass")},
        {"stored_slots",J::Number(patch.Size())},{"maximum_primitive_defect",J::Number(maximum)},
        {"tilted_diffusion_signal",J::Number(tiltedSignal)},{"omitted_tilt_signal",J::Number(omittedTiltSignal)},
        {"fixed_chart_derivative_signal",J::Number(fixedChartSignal)},
        {"minimum_sampled_covariance",J::Number(minimumCovariance)},
        {"continuous_ground_trace",ground},{"independent_ground_minimum",J::Number(exactMinimum)},
        {"negative_ground_control",negativeGround},
        {"rows",J::List(rows)},{"scope",J::Quote("native geometry/profile/trace identity; full traction and conservative k-flux evaluated on both sides; arbitrary polynomial probe fields are not a transmission solution; no variable-height evolution")}});
}

std::string ContinuousCT6Samples(const RANSWallEvolutionPatch &patch,const Tree &request)
{
    Require(patch.CurvedManufactured(),"continuous CT6 samples require the manufactured verification capability");
    const double time=request.get<double>("time",0.);
    Require(time>=0 && time<=.016,"CT6 continuous data check exceeds the contracted horizon");
    const RANSWallManufactured mms(patch);const auto &p=patch.Parameters();
    std::vector<int> sides;for(const auto &s:request.get_child("grid_sides"))sides.push_back(s.second.get_value<int>());
    Require(sides.size()==request.get_child("grid").size(),"continuous sample side count differs");
    std::vector<std::string> rows;int index=0;
    for(const auto &entry:request.get_child("grid"))
    {
        Vector values;for(const auto &x:entry.second)values.push_back(x.second.get_value<double>());
        Require(values.size()==3,"continuous sample coordinate size differs");
        const RANSWallCoordinate q{values[0],values[1],values[2]};const int side=sides[index++];
        Require(q[2]==p.geometry.height?(side==-1 || side==1):side==0,"ambiguous continuous interface sample");
        const bool outer=q[2]>p.geometry.height || (q[2]==p.geometry.height && side>0);
        const auto v=mms.At(q,time,outer);RANSTensor nativeR{},nativeQ{};double nativeP=0,nativeD=0,ep=0,ed=0;
        double nativeCutoff=0;
        if(outer)
        {
            const auto law=EvaluateRANSClosure(p.model.Model(),v.value[2],v.value[3],v.velocityGradient,q[2],p.model.Parameters(),0);
            nativeR=law.reynoldsStress;nativeQ=law.quadraticStress;nativeP=law.production;nativeD=law.kDiffusion;
            ep=law.epsilonProduction;ed=law.epsilonDestruction;
        }
        else
        {
            auto options=AWFProfileOptions{};options.formulation=AWFProfileFormulation::AcceptedBudgetAWF;
            options.curvedCapability=AWFCurvedCapability::AcceptedVerificationV1;options.quadratureOrder=96;
            const auto metric=patch.Mapping().At(q);EquilibriumAWFProfiles f(p.geometry,v.K,p.model,{metric.kappa1,metric.kappa2},options,CalibratedAWFParameters(p.model));
            const double nt=f.Original().EddyViscosity(q[2]);const auto law=EvaluateRANSWallInnerStress(v.value[2],nt,v.velocityGradient,RANSWallInnerStressModel::Linear);
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)nativeR[i][j]=(i==j?2*v.value[2]/3:0)-law.turbulentStress[i][j];
            nativeP=law.production;nativeD=p.geometry.molecularNu+nt/p.model.Parameters().sigmaK;nativeCutoff=f.Reference().cutoff;
        }
        rows.push_back(J::Object({{"coordinate",J::Numbers(q)},{"outer",J::Boolean(outer)},
            {"values",J::Numbers(v.value)},{"gradients",J::Rows(v.gradient)},{"velocity_gradient",J::Rows(v.velocityGradient)},
            {"required_body",J::Numbers(v.body)},{"k_source",J::Number(v.kSource)},{"epsilon_source",J::Number(v.epsilonSource)},
            {"R",J::Rows(v.reynolds)},{"Q",J::Rows(v.quadratic)},{"T",J::Rows(v.stress)},
            {"production",J::Number(v.production)},{"cutoff",J::Number(v.cutoff)},
            {"native_R",J::Rows(nativeR)},{"native_Q",J::Rows(nativeQ)},
            {"native_P",J::Number(nativeP)},{"native_Dk",J::Number(nativeD)},
            {"native_epsilon_production",J::Number(ep)},{"native_epsilon_destruction",J::Number(ed)},
            {"native_cutoff",outer?"null":J::Number(nativeCutoff)}}));
    }
    return J::Object({{"schema",J::Quote("CT6ContinuousNativeCrossCheckV1")},{"status",J::Quote("execution_completed")},
        {"time",J::Number(time)},{"stored_slots",J::Number(patch.Size())},{"samples",J::List(rows)},
        {"scope",J::Quote("independent continuous-data provider and production constitutive/profile evaluations; zero governing solves; primitive agreement classified by external oracle")}});
}

std::string ZeroCurvatureCT6Comparison(RANSWallEvolutionPatch &curved,
    const SpatialDomains::MeshGraphSharedPtr &graph,const std::vector<std::string> &files)
{
    Require(curved.CurvedManufactured() && !curved.Mapped() && curved.Size()<500,
        "zero-curvature comparison requires the small CT6 amplitude-zero fixture");
    char name[]="ct6-flat-comparison";char *args[]{name,nullptr};
    auto flatFiles=files;
    auto flatSession=LibUtilities::SessionReader::CreateInstance(1,args,flatFiles,LibUtilities::CommSerial::create(1,args));
    flatSession->InitSession();flatSession->SetSolverInfo("GlobalSysSoln","DirectFull");
    flatSession->SetSolverInfo("RANSWallCurvedCapability","Off");
    flatSession->SetSolverInfo("RANSWallManufactured","Off");
    flatSession->SetSolverInfo("RANSWallPilotAccuracy","Off");
    flatSession->SetSolverInfo("RANSWallBoundaryPolicy","AcceptedOpenFlatV1");
    flatSession->SetSolverInfo("RANSWallFunctionCoordinates","Physical");
    auto map=curved.Mapping();map.kind=RANSWallNormalMap::Kind::Flat;map.amplitude=0;
    RANSWallEvolutionPatch flat(flatSession,graph,curved.Comm(),curved.Parameters(),map);
    const RANSWallManufactured provider(curved);auto cd=provider.Data();
    const auto state=curved.Lift(curved.Project(provider.InitialFunctions(),0.),0.,cd);
    auto fd=cd;fd.openBoundaryIdentity=flat.OpenBoundaryIdentity();fd.Touch();
    Require(flat.Size()==curved.Size(),"zero-curvature native dimensions differ");
    const auto a=curved.Residual(state,state,0.,0.,cd,false,false),b=flat.Residual(state,state,0.,0.,fd,false,false);
    double raw=0,scale=1.;for(int i=0;i<a.size();++i){raw=std::max(raw,std::abs(a[i]-b[i]));scale=std::max({scale,std::abs(a[i]),std::abs(b[i])});}
    Tree ca,fa;std::istringstream ci(curved.OpenFrozenSystemsPacket(state,0.,cd)),fi(flat.OpenFrozenSystemsPacket(state,0.,fd));
    boost::property_tree::read_json(ci,ca);boost::property_tree::read_json(fi,fa);
    Require(ca.get_child("blocks").size()==fa.get_child("blocks").size(),"zero-curvature block count differs");
    double matrix=0;auto ib=ca.get_child("blocks").begin(),jb=fa.get_child("blocks").begin();
    for(;ib!=ca.get_child("blocks").end();++ib,++jb)
    {
        Require(jb!=fa.get_child("blocks").end() && ib->second.get<int>("block")==jb->second.get<int>("block"),"zero-curvature block identities differ");
        Require(Values(ib->second,"rows")==Values(jb->second,"rows"),"zero-curvature constrained rows differ");
        Require(ib->second.get_child("matrix").size()==jb->second.get_child("matrix").size(),"zero-curvature matrix row count differs");
        auto ir=ib->second.get_child("matrix").begin(),jr=jb->second.get_child("matrix").begin();
        for(;ir!=ib->second.get_child("matrix").end();++ir,++jr)
        {
            Require(jr!=jb->second.get_child("matrix").end(),"zero-curvature matrix shape differs");
            Require(ir->second.size()==jr->second.size(),"zero-curvature matrix column count differs");
            auto x=ir->second.begin(),y=jr->second.begin();
            for(;x!=ir->second.end();++x,++y)
            {Require(y!=jr->second.end(),"zero-curvature matrix width differs");matrix=std::max(matrix,std::abs(x->second.get_value<double>()-y->second.get_value<double>()));}
        }
    }
    double profile=0,constitutive=0;
    const RANSWallCoordinate H{.1,0,curved.Parameters().geometry.height};
    const auto cp=curved.StateProfile(state,H,0.,cd),fp=flat.StateProfile(state,H,0.,fd);
    for(double r:{0.,.1,.25,.5})
        for(const auto &pair:std::vector<std::pair<double,double>>{{cp->Dissipation(r),fp->Dissipation(r)},
            {cp->ReferenceK(r),fp->ReferenceK(r)},{cp->ReferenceKGradient(r),fp->ReferenceKGradient(r)},
            {cp->ReferenceKFlux(r),fp->ReferenceKFlux(r)},{cp->ReferenceKDiffusivity(r),fp->ReferenceKDiffusivity(r)}})
            profile=std::max(profile,std::abs(pair.first-pair.second));
    for(int side:{-1,1})
    {
        RANSTensor cg{},fg{};for(int component=0;component<2;++component)
        {cg[component?2:0]=curved.Point(0,H,side).Gradient(state,curved.offset[component]);fg[component?2:0]=flat.Point(0,H,side).Gradient(state,flat.offset[component]);}
        const auto kv=curved.Point(0,H,side).Value(state,curved.offset[curved.K()]);
        if(side>0)
        {
            const double e=curved.Point(curved.Epsilon(),H,1).Value(state,curved.offset[curved.Epsilon()]);
            const auto x=curved.OuterClosure(kv,e,cg,H),y=flat.OuterClosure(kv,e,fg,H);
            constitutive=std::max(constitutive,std::abs(x.production-y.production));
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)constitutive=std::max(constitutive,std::abs(x.reynoldsStress[i][j]-y.reynoldsStress[i][j]));
        }
        else
        {
            const auto x=curved.InnerStress(kv,cp->Original().EddyViscosity(H[2]),cg),y=flat.InnerStress(kv,fp->Original().EddyViscosity(H[2]),fg);
            constitutive=std::max(constitutive,std::abs(x.production-y.production));
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)constitutive=std::max(constitutive,std::abs(x.reynoldsStress[i][j]-y.reynoldsStress[i][j]));
        }
    }
    Require(std::max({raw/scale,matrix,profile,constitutive})<=2e-10,"zero-curvature flat/CT6 operator disagreement");
    return J::Object({{"schema",J::Quote("CT6ZeroCurvatureFlatComparisonV1")},{"status",J::Quote("pass")},
        {"stored_slots",J::Number(curved.Size())},{"raw_scaled_defect",J::Number(raw/scale)},
        {"frozen_matrix_defect",J::Number(matrix)},{"profile_defect",J::Number(profile)},
        {"constitutive_defect",J::Number(constitutive)},
        {"curved_identity",J::Quote(curved.PreservationIdentity())},{"flat_identity",J::Quote(flat.PreservationIdentity())},
        {"scope",J::Quote("same native state and same physical C++ functions passed explicitly to ordinary flat and CT6-zero operators; policy identities remain distinct; no checkpoint relabelling or governing solve")}});
}

// Physical volume/face integration, deliberately independent of assembled
// constant-test loads and prescribed-row reactions. IMEX1, unforced one sweep.
std::string PhysicalScalarBudget(RANSWallEvolutionPatch &patch,
    const RANSWallPatchParameters &parameters,const Vector &state,const Vector &previous,
    double dt,double time,const RANSWallEvolutionData &data,int faceOrder,
    bool reconcile=false,const Vector &stageHistory={},double diffusionIBPFraction=0.)
{
    Require(dt>0 && previous.size()==state.size(),"invalid physical budget history");
    Require(!patch.Mapped() || diffusionIBPFraction==0,"Cartesian finite-difference IBP diagnostic is not a mapped derivative");
    Require(faceOrder==16 || faceOrder==24 || faceOrder==32,"unsupported physical budget face order");
    auto flow=previous;
    for(int b:{0,1,patch.InnerPressure(),patch.OuterPressure()})
        std::copy_n(state.begin()+patch.offset[b],patch.size[b],flow.begin()+patch.offset[b]);
    flow.back()=state.back();auto preK=flow;
    std::copy_n(state.begin()+patch.offset[patch.Epsilon()],patch.size[patch.Epsilon()],preK.begin()+patch.offset[patch.Epsilon()]);
    Tree packet;std::istringstream input(patch.PreservationPacket(preK,time,data,true,{},{}));
    boost::property_tree::read_json(input,packet);
    std::array<Vector,3> sum{Vector(7),Vector(7),Vector(7)};
    std::array<std::array<Vector,8>,3> loads;
    std::array<Vector,3> advectiveFaces,diffusiveFaces;
    if(reconcile)for(int region=0;region<3;++region)
    {
        const int block=region==2?patch.Epsilon():patch.K();
        for(auto &v:loads[region])v.resize(patch.size[block]);
        advectiveFaces[region].resize(patch.size[block]);
        diffusiveFaces[region].resize(patch.size[block]);
    }
    const auto model=parameters.model.Parameters();
    auto gradient=[&](const RANSWallGlobalPoint &v)
    {RANSTensor g{};g[0]=v.Gradient(state,patch.offset[0]);g[2]=v.Gradient(state,patch.offset[1]);return g;};
    auto coefficients=[&](const RANSWallGlobalPoint &v,bool outer)
    {
        const auto g=gradient(v);const double oldK=v.Value(previous,patch.offset[patch.K()]),newK=v.Value(state,patch.offset[patch.K()]);
        Vector c(12); // original eight totals, then independent frozen RHS/rates
        if(!outer)
        {
            const auto frozen=patch.StateProfile(preK,v.coordinate,time,data),current=patch.StateProfile(state,v.coordinate,time,data);
            const double nt=frozen->Original().EddyViscosity(patch.Distance(v.coordinate)),now=current->Original().EddyViscosity(patch.Distance(v.coordinate));
            const auto law=patch.InnerStress(oldK,nt,g),actual=patch.InnerStress(newK,now,g);
            double rate=0;const auto source=patch.WallSource(preK,v.coordinate,law.production,*frozen,nullptr,&rate);
            Require(rate==0,"prescribed wall sink unexpectedly has local reaction");
            c[0]=parameters.geometry.molecularNu+nt/model.sigmaK;c[1]=parameters.geometry.molecularNu+now/model.sigmaK;
            c[4]=source.kSource;c[5]=patch.WallSource(state,v.coordinate,actual.production,*current).kSource;
            c[8]=source.kSource;c[9]=0;
        }
        else
        {
            const auto e=patch.GeometryPoint(patch.Epsilon(),v.geometry,v.coordinate);
            const double oldE=e.Value(previous,patch.offset[patch.Epsilon()]),newE=e.Value(state,patch.offset[patch.Epsilon()]);
            const auto kLaw=patch.OuterClosure(oldK,newE,g,v.coordinate),eLaw=patch.OuterClosure(oldK,oldE,g,v.coordinate),current=patch.OuterClosure(newK,newE,g,v.coordinate);
            c[0]=kLaw.kDiffusion;c[1]=current.kDiffusion;c[2]=eLaw.epsilonDiffusion;c[3]=current.epsilonDiffusion;
            c[4]=std::max(0.,kLaw.production)-(newE-std::min(0.,kLaw.production))*newK/oldK;
            c[5]=current.kSource;c[6]=eLaw.epsilonProduction-model.cEpsilon2*oldE/oldK*newE;c[7]=current.epsilonSource;
            c[8]=std::max(0.,kLaw.production);c[9]=(newE-std::min(0.,kLaw.production))/oldK;
            c[10]=eLaw.epsilonProduction;c[11]=model.cEpsilon2*oldE/oldK;
        }
        if(patch.OpenBoundary())
        {
            const double kSource=(outer?data.outerKSource:data.innerKSource)(v.coordinate,time);
            c[4]+=kSource;c[5]+=kSource;c[8]+=kSource;
            if(outer){const double eSource=data.outerEpsilonSource(v.coordinate,time);c[6]+=eSource;c[7]+=eSource;c[10]+=eSource;}
        }
        return c;
    };
    auto weak=[&](int region,const RANSWallGlobalPoint &q,const RANSWallGlobalPoint &v,
                  const Vector &c,double weight)
    {
        if(!reconcile)return;
        const bool epsilon=region==2;const int block=epsilon?patch.Epsilon():patch.K();
        const double f=q.Value(state,patch.offset[block]),old=q.Value(previous,patch.offset[block]);
        const auto grad=q.Gradient(state,patch.offset[block]);
        const auto ux=v.Gradient(state,patch.offset[0]),uz=v.Gradient(state,patch.offset[1]);
        const double u=v.Value(state,patch.offset[0]),w=v.Value(state,patch.offset[1]);
        const double D=c[epsilon?2:0],rhs=c[epsilon?10:8],rate=c[epsilon?11:9];
        double divDiffusion=0;
        if(diffusionIBPFraction)
        {
            Require(diffusionIBPFraction==std::ldexp(1.,-11) || diffusionIBPFraction==std::ldexp(1.,-12),
                    "undeclared read-only diffusion IBP stencil");
            const auto &el=patch.fields[0]->Elements()[v.element];
            for(int a:{0,2})
            {
                const double h=(el.upper[a]-el.lower[a])*diffusionIBPFraction;
                auto flux=[&](int shift)
                {
                    auto p=v.coordinate;p[a]+=shift*h;
                    const auto vv=patch.GeometryPoint(0,v.geometry,p);
                    const auto cc=coefficients(vv,region!=0);
                    const auto scalar=patch.GeometryPoint(block,v.geometry,p);
                    return cc[epsilon?2:0]*scalar.Gradient(state,patch.offset[block])[a];
                };
                if(v.coordinate[a]-2*h>=el.lower[a] && v.coordinate[a]+2*h<=el.upper[a])
                    divDiffusion+=(-flux(2)+8*flux(1)-8*flux(-1)+flux(-2))/(12*h);
                else
                {
                    const int side=v.coordinate[a]+4*h<=el.upper[a]?1:-1;
                    Require(v.coordinate[a]+side*4*h>=el.lower[a] && v.coordinate[a]+side*4*h<=el.upper[a],
                            "diffusion IBP stencil cannot fit its generating owner");
                    divDiffusion+=side*(-25*flux(0)+48*flux(side)-36*flux(2*side)+16*flux(3*side)-3*flux(4*side))/(12*h);
                }
            }
        }
        for(int i=0;i<q.ids.size();++i)
        {
            const int id=q.ids[i];const double test=q.phi[i];
            loads[region][0][id]+=weight*test*(f-old)/dt;
            loads[region][1][id]-=weight*test*rhs;
            loads[region][2][id]+=weight*test*rate*f;
            loads[region][3][id]-=weight*(q.derivative[0][i]*u+q.derivative[2][i]*w)*f;
            loads[region][4][id]+=weight*D*(q.derivative[0][i]*grad[0]+q.derivative[2][i]*grad[2]);
            // Conservative strong advection, used only for an IBP check.
            // The f*div(U) term is not added to the physical budget or operator.
            loads[region][5][id]+=weight*test*(u*grad[0]+w*grad[2]+f*(ux[0]+uz[2]));
            if(diffusionIBPFraction)loads[region][6][id]+=weight*test*divDiffusion;
        }
    };
    for(const auto &item:packet.get_child("volume"))
    {
        const auto &q=item.second;const auto coordinate=Values(q,"coordinate");RANSWallCoordinate p{};std::copy_n(coordinate.begin(),3,p.begin());
        const int region=q.get<int>("region"),geometry=q.get<int>("velocity_k.geometry");const auto v=patch.GeometryPoint(0,geometry,p);const auto c=coefficients(v,region==1);const double w=q.get<double>("weight");
        sum[region][0]+=w*(v.Value(state,patch.offset[patch.K()])-v.Value(previous,patch.offset[patch.K()]))/dt;
        sum[region][1]+=w*c[4];sum[region][2]+=w*c[5];sum[region][5]+=w;
        weak(region,v,v,c,w);
        if(region==1)
        {
            const auto e=patch.GeometryPoint(patch.Epsilon(),geometry,p);
            sum[2][0]+=w*(e.Value(state,patch.offset[patch.Epsilon()])-e.Value(previous,patch.offset[patch.Epsilon()]))/dt;
            sum[2][1]+=w*c[6];sum[2][2]+=w*c[7];sum[2][5]+=w;
            weak(2,e,v,c,w);
        }
    }
    std::string reconciliation="null";
    if(reconcile)
    {
        Require(stageHistory.size()==state.size(),"frozen audit requires actual native stage history");
        std::vector<std::pair<std::string,std::string>> regional;
        const auto fixed=patch.Constraints(flow,time,data);
        if(patch.OpenBoundary())for(const auto &face:patch.Surface(0,0,patch.upper[0],1,faceOrder))
        {
            const auto &v=face.test;const int side=patch.fields[0]->Elements()[v.element].lower[2]>=parameters.geometry.height?1:0;
            const double un=face.area[0]*v.Value(state,patch.offset[0])+face.area[2]*v.Value(state,patch.offset[1]);
            const double area=std::hypot(face.area[0],face.area[2]);
            auto add=[&](int region,int block,const RANSWallGlobalPoint &q,double datum)
            {
                const double load=un*q.Value(state,patch.offset[block])-area*datum;
                for(int i=0;i<q.ids.size();++i)loads[region][7][q.ids[i]]+=face.weight*q.phi[i]*load;
            };
            add(side,patch.K(),v,data.outletKDiffusion[side](v.coordinate,time));
            if(side)add(2,patch.Epsilon(),patch.GeometryPoint(patch.Epsilon(),v.geometry,v.coordinate),data.outletEpsilonDiffusion(v.coordinate,time));
        }
        for(int region=0;region<3;++region)
        {
            const int block=region==2?patch.Epsilon():patch.K();
            std::array<std::set<double>,2> cuts;
            for(const auto &e:patch.fields[block]->Elements())
            {
                if(region<2 && int(e.lower[2]>=parameters.geometry.height)!=region)continue;
                cuts[0].insert(e.lower[0]);cuts[0].insert(e.upper[0]);
                cuts[1].insert(e.lower[2]);cuts[1].insert(e.upper[2]);
            }
            for(int d=0;d<2;++d)for(double plane:cuts[d])for(int side:{-1,1})
                for(const auto &face:patch.Surface(block,d==0?0:2,plane,side,faceOrder))
                {
                    const auto &q=face.test;const auto &e=patch.fields[block]->Elements()[q.element];
                    if(region<2 && int(e.lower[2]>=parameters.geometry.height)!=region)continue;
                    const auto v=patch.GeometryPoint(0,q.geometry,q.coordinate);
                    const double adv=(face.area[0]*v.Value(state,patch.offset[0])+face.area[2]*v.Value(state,patch.offset[1]))*
                        q.Value(state,patch.offset[block]);
                    const auto c=coefficients(v,region!=0);const auto g=q.Gradient(state,patch.offset[block]);
                    const double diff=c[region==2?2:0]*(face.area[0]*g[0]+face.area[2]*g[2]);
                    for(int i=0;i<q.ids.size();++i)
                    {
                        advectiveFaces[region][q.ids[i]]+=face.weight*q.phi[i]*adv;
                        diffusiveFaces[region][q.ids[i]]+=face.weight*q.phi[i]*diff;
                    }
                }
            Vector ibp(patch.size[block]),diffusionIBP(patch.size[block]);for(int j=0;j<ibp.size();++j)
            {
                ibp[j]=loads[region][5][j]-loads[region][3][j]-advectiveFaces[region][j];
                diffusionIBP[j]=loads[region][6][j]+loads[region][4][j]-diffusiveFaces[region][j];
            }
            regional.push_back({region==0?"inner_k":region==1?"outer_k":"outer_epsilon",J::Object({
                {"backward_storage",J::Numbers(loads[region][0])},{"negative_frozen_RHS",J::Numbers(loads[region][1])},
                {"reaction",J::Numbers(loads[region][2])},{"conservative_advection",J::Numbers(loads[region][3])},
                {"diffusion",J::Numbers(loads[region][4])},{"strong_conservative_advection",J::Numbers(loads[region][5])},
                {"outlet_weak_load",J::Numbers(loads[region][7])},
                {"all_element_advective_faces",J::Numbers(advectiveFaces[region])},{"advective_IBP_defect",J::Numbers(ibp)},
                {"all_element_diffusive_faces",J::Numbers(diffusiveFaces[region])},
                {"independent_strong_diffusion",diffusionIBPFraction?J::Numbers(loads[region][6]):"null"},
                {"diffusion_IBP_defect",diffusionIBPFraction?J::Numbers(diffusionIBP):"null"}})});
        }
        std::vector<std::pair<std::string,std::string>> assembled;
        for(int block:{patch.K(),patch.Epsilon()})
        {
            const auto &frozen=block==patch.K()?preK:flow;
            const auto raw=patch.FrozenScalarWeakLoad(block,frozen,state,stageHistory,1/dt,time,data);
            Vector reconstructed(raw.size()),storageDifference(raw.size());std::vector<int> free,prescribed;
            const int first=block==patch.K()?0:2,last=block==patch.K()?2:3;
            for(int region=first;region<last;++region)for(int term=0;term<5;++term)
                for(int j=0;j<raw.size();++j)reconstructed[j]+=loads[region][term][j];
            if(patch.OpenBoundary())for(int region=first;region<last;++region)
                for(int j=0;j<raw.size();++j)reconstructed[j]+=loads[region][7][j];
            Vector defect=raw;for(int j=0;j<raw.size();++j)
            {
                defect[j]-=reconstructed[j];
                (fixed.count(patch.offset[block]+j)?prescribed:free).push_back(j);
                storageDifference[j]=previous[patch.offset[block]+j]-stageHistory[patch.offset[block]+j];
            }
            assembled.push_back({block==patch.K()?"k":"outer_epsilon",J::Object({
                {"native_raw_frozen_load",J::Numbers(raw)},{"independent_term_sum",J::Numbers(reconstructed)},
                {"reconstruction_defect",J::Numbers(defect)},{"previous_minus_stage_history",J::Numbers(storageDifference)},
                {"free_rows",J::Numbers(free)},{"prescribed_or_dependent_rows",J::Numbers(prescribed)}})});
        }
        reconciliation=J::Object({{"schema",J::Quote("ConservativeFrozenScalarReconciliationV1")},
            {"regions",J::Object(regional)},{"assembled",J::Object(assembled)},
            {"transport_velocity",J::Quote("actual corrected flow, unchanged by subsequent scalar solves")},
            {"epsilon_coefficients",J::Quote("old k and old epsilon; transported epsilon is newly solved, including dependent trace storage")},
            {"k_coefficients",J::Quote("old k and newly solved epsilon; transported k is newly solved")},
            {"diffusion_IBP_element_step_fraction",J::Number(diffusionIBPFraction)},
            {"diffusion_IBP_method",J::Quote("optional fourth-order physical flux differences within each known generating owner, central or declared one-sided near its endpoint; h/h2 control; forward constitutive values only, no tangent or ellipticity claim")},
            {"scope",J::Quote("read-only actual unmasked frozen scalar matrix and independently rebuilt weak terms; prescribed reactions are algebraic diagnostics, not physical flux substitutes; no scalar*div(U) correction")}});
    }
    std::vector<std::string> faces;
    for(int region=0;region<2;++region)for(int axis:{0,2})for(int side:{-1,1})
    {
        const double z=axis==0?(side<0?patch.lower[0]:patch.upper[0]):
            (side<0?(region==0?0.:parameters.geometry.height):(region==0?parameters.geometry.height:patch.upper[2]));
        for(const auto &face:patch.Surface(0,axis,z,side,faceOrder))
        {
            const auto &v=face.test;const auto &element=patch.fields[0]->Elements().at(v.element);
            if((element.lower[2]>=parameters.geometry.height?1:0)!=region)continue;
            const auto c=coefficients(v,region==1);const auto gk=v.Gradient(state,patch.offset[patch.K()]);
            const double k=v.Value(state,patch.offset[patch.K()]);RANSVector u{};u[0]=v.Value(state,patch.offset[0]);u[2]=v.Value(state,patch.offset[1]);
            double frozen=0,current=0;
            for(int a=0;a<3;++a){frozen+=face.area[a]*(c[0]*gk[a]-u[a]*k);current+=face.area[a]*(c[1]*gk[a]-u[a]*k);}
            // Diagnostic components only. Keep the original net-flux and
            // budget accumulation above unchanged; no row reaction is used.
            double ka=0,kdf=0,kdc=0;
            for(int a=0;a<3;++a){ka+=face.area[a]*u[a]*k;kdf+=face.area[a]*c[0]*gk[a];kdc+=face.area[a]*c[1]*gk[a];}
            sum[region][3]+=face.weight*frozen;sum[region][4]+=face.weight*current;
            double ef=0,ec=0,ea=0,edf=0,edc=0;
            if(region==1)
            {
                const auto e=patch.GeometryPoint(patch.Epsilon(),v.geometry,v.coordinate);const auto ge=e.Gradient(state,patch.offset[patch.Epsilon()]);const double value=e.Value(state,patch.offset[patch.Epsilon()]);
                for(int a=0;a<3;++a){ef+=face.area[a]*(c[2]*ge[a]-u[a]*value);ec+=face.area[a]*(c[3]*ge[a]-u[a]*value);}
                for(int a=0;a<3;++a){ea+=face.area[a]*u[a]*value;edf+=face.area[a]*c[2]*ge[a];edc+=face.area[a]*c[3]*ge[a];}
                sum[2][3]+=face.weight*ef;sum[2][4]+=face.weight*ec;
            }
            faces.push_back(J::Object({{"coordinate",J::Numbers(v.coordinate)},{"geometry",J::Number(v.geometry)},
                {"region",J::Number(region)},{"axis",J::Number(axis)},{"outward_side",J::Number(side)},
                {"weight",J::Number(face.weight)},{"area",J::Numbers(face.area)},
                {"component_schema",J::Quote("SeparatePhysicalFaceFluxV1")},
                {"k_advective_outward",J::Number(ka)},
                {"k_frozen_diffusive_outward",J::Number(kdf)},
                {"k_current_diffusive_outward",J::Number(kdc)},
                {"epsilon_advective_outward",region?J::Number(ea):"null"},
                {"epsilon_frozen_diffusive_outward",region?J::Number(edf):"null"},
                {"epsilon_current_diffusive_outward",region?J::Number(edc):"null"},
                {"k_frozen_flux",J::Number(frozen)},{"k_current_flux",J::Number(current)},
                {"epsilon_frozen_flux",region?J::Number(ef):"null"},{"epsilon_current_flux",region?J::Number(ec):"null"}}));
        }
    }
    std::vector<std::pair<std::string,std::string>> rows;
    for(int i=0;i<3;++i)
    {
        const auto &s=sum[i];rows.push_back({i==0?"inner_k":i==1?"outer_k":"outer_epsilon",J::Object({
            {"storage_backward_difference",J::Number(s[0])},{"stage_source",J::Number(s[1])},
            {"instantaneous_source",J::Number(s[2])},{"stage_physical_boundary_flux",J::Number(s[3])},
            {"instantaneous_physical_boundary_flux",J::Number(s[4])},{"volume",J::Number(s[5])},
            {"frozen_stage_physical_defect",J::Number(s[0]-s[1]-s[3])},
            {"instantaneous_physical_defect",J::Number(s[0]-s[2]-s[4])},
            {"frozen_operator_splitting",J::Number(s[1]+s[3]-s[2]-s[4])}})});
    }
    return J::Object({{"schema",J::Quote("PhysicalFrozenStageScalarBudgetV1")},{"time",J::Number(time)},
        {"dt",J::Number(dt)},{"equations",J::Object(rows)},{"physical_faces",J::List(faces)},
        {"face_order",J::Number(faceOrder)},{"wall_local_reaction_rate","0"},
        {"reconciliation",reconciliation},
        {"scope",J::Quote("independently integrated storage, frozen one-sweep scalar source and physical gradient/advective face flux; no prescribed-row reaction used. Dependent epsilon boundary storage is retained, with no inner epsilon equation. Finite discretization defects are observations, not forced balances.")}});
}
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
        std::string boundaryPolicy;session->LoadSolverInfo("RANSWallBoundaryPolicy",boundaryPolicy,"Closed");
        const bool open=boundaryPolicy=="AcceptedOpenFlatV1" || boundaryPolicy=="AcceptedManufacturedOpenV1";
        for(const auto &name:{"BodyForceInner","BodyForceOuter","ScalarSourceInner","ScalarSourceOuter"})
            Require(open || !session->DefinesFunction(name),"neutral preservation verifier requires unforced source definitions");
        std::string referenceMode;session->LoadSolverInfo("RANSWallReference",referenceMode,"Neutral");
        std::string curved;session->LoadSolverInfo("RANSWallCurvedCapability",curved,"Off");
        const bool staticCT6=curved=="CT6StaticTransmissionV1";
        Require(referenceMode=="Neutral" || ((open || staticCT6) && referenceMode=="Functions"),"neutral verification requires the Neutral reference/boundary contract");
        double eta;session->LoadParameter("RANSWallEta",eta,0.);Require(eta==0,"unanchored verification requires eta zero");
        RANSWallEvolutionPatch patch(session,graph,world->GetComm(),p,LoadRANSWallNormalMap(session,2,p.geometry.height));
        if(staticCT6)
        {
            Require(request.get<std::string>("operation")=="ct6_static_transmission","static CT6 route cannot replay or evolve a history");
            std::ofstream out(argv[3]);out<<StaticCT6Transmission(patch)<<'\n';return 0;
        }
        if(request.get<std::string>("operation","")=="ct6_continuous_samples")
        {std::ofstream out(argv[3]);out<<ContinuousCT6Samples(patch,request)<<'\n';return 0;}
        if(request.get<std::string>("operation","")=="ct6_zero_flat")
        {std::ofstream out(argv[3]);out<<ZeroCurvatureCT6Comparison(patch,graph,files)<<'\n';return 0;}
        Require(patch.Dimension()==2 && (!patch.Mapped() || patch.CurvedVerification()) && (patch.AcceptedAWF() || patch.UniformEPFlat()),
            "neutral preservation probe supports accepted or EP plateau references, serial flat2D only");
        RANSWallEvolutionData data(2);Vector state,reference;double time=request.get<double>("time",0.);
        if(open)data=LoadRANSAcceptedOpenData(session,patch);
        const std::string operation=request.get<std::string>("operation","neutral");
        const bool legacy=operation=="legacy_operator_replay";
        const std::string sourceOperator=legacy?"LegacyCoordinateToleranceV1":RANSWallNumericalOperatorIdentity();
        if(operation=="neutral")
        {
            Require(!open,"open-boundary sampling requires actual Functions or a bound replay, not a closed neutral initializer");
            double us,dx,dy;session->LoadParameter("RANSWallFrictionVelocity",us,.4);
            session->LoadParameter("RANSWallDirectionX",dx,1.);session->LoadParameter("RANSWallDirectionY",dy,0.);
            state=patch.NeutralReference(us,{dx,dy});data.fixedBoundary=patch.FixedReferenceBoundary(state);
            state=patch.Lift(state,time,data);reference=state;
        }
        else if(operation=="functions")
        {
            Require(open && referenceMode=="Functions" && time==0.,"Functions verification is a fresh open-policy initialization only");
            state=patch.Lift(data.projectedReference.empty()?patch.Project(RANSAcceptedOpenInitialFunctions(session,patch),time):data.projectedReference,time,data);reference=state;
            if(!patch.CurvedManufactured())data.fixedBoundary=patch.FixedReferenceBoundary(state);
        }
        else
        {
            Require(operation=="replay" || legacy,"unknown read-only verification operation");
            if(legacy)
            {
                Require(request.get<std::string>("source_numerical_operator")==sourceOperator &&
                    request.get<std::string>("target_numerical_operator")==RANSWallNumericalOperatorIdentity(),
                    "legacy read-only comparison requires both numerical operator identities");
                Require(!request.get_optional<std::string>("write_initial_history"),
                    "legacy operator coefficients cannot create or continue a history");
            }
            Require(request.get<std::string>("identity")==patch.PreservationIdentity(legacy),"verification physics/mesh/space/operator identity mismatch");
            Require(request.get<std::string>("descriptor")==p.model.Signature(),"verification resolved descriptor mismatch");
            state=Values(request,"state");reference=Values(request,"reference");
            Require(state.size()==std::size_t(patch.Size()) && reference.size()==state.size(),"verification coefficient size mismatch");
            for(const auto &item:request.get_child("fixed_boundary"))
            {auto q=item.second.begin();const int row=q->second.get_value<int>();++q;const double value=q->second.get_value<double>();Require(row>=0 && row<patch.Size(),"invalid fixed boundary row");Require(data.fixedBoundary.emplace(row,value).second,"duplicate fixed boundary row");}
            Require(patch.CurvedManufactured()?data.fixedBoundary.empty():data.fixedBoundary==patch.FixedReferenceBoundary(reference),"verification external boundary identity mismatch");
        }
        // CT6 has time-dependent live boundary functions and a fixed initial
        // reference. Validate/store that reference at its actual time zero,
        // just as native restart does before restoring the evolved time.
        // Current-state rows, boundary functions and residuals below still
        // use the requested observation time. This is diagnostic-only; no
        // field or history is reset and no residual enters evolution.
        const double referenceTime=patch.CurvedManufactured()?0.:time;
        patch.SetReference(reference,data,referenceTime,0.);
        std::vector<RANSWallCoordinate> grid;
        if(auto points=request.get_child_optional("grid"))for(const auto &item:*points)
        {RANSWallCoordinate q{};int i=0;for(const auto &v:item.second){Require(i<3,"invalid verification coordinate");q[i++]=v.second.get_value<double>();}Require(i==3,"invalid verification coordinate");grid.push_back(q);}
        std::vector<int> sides;if(auto selected=request.get_child_optional("grid_sides"))for(const auto &v:*selected)sides.push_back(v.second.get_value<int>());
        double receiptUs,receiptDx,receiptDy;session->LoadParameter("RANSWallFrictionVelocity",receiptUs,.4);session->LoadParameter("RANSWallDirectionX",receiptDx,1.);session->LoadParameter("RANSWallDirectionY",receiptDy,0.);
        const auto targetReference=open?reference:patch.NeutralReference(receiptUs,{receiptDx,receiptDy});
        const auto receipt=patch.InitializationReceipt(state,reference,time,data,data,patch.FixedReferenceBoundary(targetReference));
        const auto stateCopy=state;
        const auto observation=patch.PreservationPacket(state,time,data,request.get<bool>("export_basis",false),grid,sides,
            request.get<bool>("export_ct6_primitives",false));
        std::string resolution="null",testBasis="null",constitutive="null";
        if(request.get<bool>("export_constitutive",false))
        {
            std::vector<std::string> rows;
            for(int index=0;index<grid.size();++index)
            {
                const auto &q=grid[index];const int side=sides.empty()?0:sides[index];
                const bool atInterface=q[2]==p.geometry.height;
                const bool outer=atInterface?side==1:q[2]>p.geometry.height;
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
                RANSVector normal{0,0,1};double measure=1;
                if(patch.CurvedVerification())
                {
                    measure=std::hypot(f.area[0],f.area[1],f.area[2]);
                    for(int a=0;a<3;++a)normal[a]=side*f.area[a]/measure;
                }
                const auto normalDerivative=[&](const auto &gradient)
                {double v=0;for(int a=0;a<3;++a)v+=gradient[a]*normal[a];return v;};
                double dk=0,de=0,epsilon=0,epsilonFlux=0;
                if(outer)
                {
                    const auto erow=patch.GeometryPoint(patch.Epsilon(),q.geometry,q.coordinate);
                    epsilon=erow.Value(state,patch.offset[patch.Epsilon()]);
                    const auto c=patch.OuterClosure(k,epsilon,g,q.coordinate);dk=c.kDiffusion;de=c.epsilonDiffusion;
                    epsilonFlux=de*normalDerivative(erow.Gradient(state,patch.offset[patch.Epsilon()]));
                    for(int a=0;a<3;++a)for(int b=0;b<3;++b)stress[a][b]=-c.reynoldsStress[a][b]+(a==b?2*k/3-pi:0.);
                }
                else
                {
                    const auto profile=patch.StateProfile(state,q.coordinate,time,data);
                    const double nt=profile->Original().EddyViscosity(patch.Distance(q.coordinate));
                    dk=p.geometry.molecularNu+nt/p.model.Parameters().sigmaK;
                    stress=RANSWallInnerMeanStress(patch.InnerStress(k,nt,g),p.geometry.molecularNu,pi);
                }
                Vector traction(3);for(int a=0;a<3;++a)for(int b=0;b<3;++b)traction[a]+=stress[a][b]*normal[b];
                const double vertical=q.Value(state,patch.offset[1]);
                const double un=q.Value(state,patch.offset[0])*normal[0]+vertical*normal[2];
                const double kFlux=dk*normalDerivative(q.Gradient(state,patch.offset[patch.K()]));
                fluxes.push_back(J::Object({{"coordinate",J::Numbers(q.coordinate)},
                    {"weight",J::Number(f.weight*measure)},{"chart_weight",J::Number(f.weight)},
                    {"positive_chart_normal",J::Numbers(normal)},{"area",J::Numbers(f.area)},
                    {"geometry",J::Number(q.geometry)},
                    {"region",outer?"1":"0"},{"side",J::Number(-side)},
                    {"outward_normal_z",J::Number(side)},
                    {"velocity",J::Numbers(Vector{q.Value(state,patch.offset[0]),vertical})},
                    {"k",J::Number(k)},{"epsilon",outer?J::Number(epsilon):"null"},
                    {"Pi",J::Number(pi)},{"velocity_gradient",J::Rows(g)},
                    {"k_gradient",J::Numbers(q.Gradient(state,patch.offset[patch.K()]))},
                    {"full_traction",J::Numbers(traction)},
                    {"k_diffusive_flux",J::Number(kFlux)},
                    {"k_conservative_flux",J::Number(kFlux-un*k)},
                    {"epsilon_diffusive_flux",outer?J::Number(epsilonFlux):"null"},
                    {"epsilon_conservative_flux",outer?J::Number(epsilonFlux-un*epsilon):"null"},
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
                const bool atInterface=q[2]==p.geometry.height;
                bool outer=atInterface?side==1:q[2]>p.geometry.height;
                if(!geometries.empty())
                {
                    const auto &elements=patch.fields[0]->Elements();
                    const auto parent=std::find_if(elements.begin(),elements.end(),[&](const auto &e){return e.geometry==geometries[i];});
                    Require(parent!=elements.end(),"unknown explicit native test geometry");
                    const bool selectedOuter=parent->lower[2]>=p.geometry.height;
                    if(q[2]==p.geometry.height)
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
        std::string tensorValidation="null",mixedSpace="null",openSystems="null";
        if(request.get<bool>("export_tensor_validation",false))tensorValidation=RANSTensorValidation(patch);
        if(request.get<bool>("export_mixed_space",false))mixedSpace=patch.MixedSpaceRankPacket(state,time,data);
        if(request.get<bool>("export_open_verification",false))openSystems=patch.OpenFrozenSystemsPacket(state,time,data);
        std::string physicalBudget="null";
        if(request.get<bool>("export_physical_scalar_budget",false))
        {
            const auto scheme=session->GetTimeIntScheme();
            Require(scheme.method=="IMEX" && scheme.variant.empty() && scheme.order==1,"physical stage budget requires IMEX1");
            double dt;session->LoadParameter("TimeStep",dt);Require(time>0,"physical stage budget needs an advanced state");
            const bool reconciliation=request.get<bool>("export_frozen_scalar_reconciliation",false);
            physicalBudget=PhysicalScalarBudget(patch,p,state,Values(request,"previous_state"),dt,time,data,request.get<int>("physical_face_order",16),
                reconciliation,reconciliation?Values(request,"stage_history"):Vector{},request.get<double>("diffusion_IBP_element_step_fraction",0.));
        }
        Require(state==stateCopy,"read-only verification modified coefficients");
        const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
        Tree observed;std::istringstream input(observation);boost::property_tree::read_json(input,observed);
        const auto &audit=observed.get_child("audit_alpha_zero");
        bool physical=audit.get_child("accepted_awf_scalar_domain").front().second.get_value<double>()==1 &&
            audit.get_child_optional("accepted_awf_hard_admissible") &&
            audit.get_child("accepted_awf_hard_admissible").front().second.get_value<double>()==1;
        std::string openPacket;
        if(open)
        {
            bool outflow=false;patch.OpenBoundaryStatus(state,time,data,outflow);physical=physical && outflow;
            openPacket=patch.OpenBoundaryPacket(state,time,data);
        }
        std::string initialHistory;
        if(auto target=request.get_optional<std::string>("write_initial_history"))
        {
            Require(!open && operation=="replay" && time==0. && physical &&
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
        std::string output=J::Object({{"status",J::Quote(physical?"pass":"fail")},{"execution",J::Quote("completed")},
            {"source_numerical_operator",J::Quote(sourceOperator)},
            {"numerical_operator",J::Quote(RANSWallNumericalOperatorIdentity())},
            {"physical_admissibility",J::Quote(physical?"pass":"fail")},
            {"accuracy",J::Quote("not_assessed_by_sampling_adapter")},{"performed_time_steps","0"},{"initial_history",initialHistory.empty()?"null":J::Quote(initialHistory)},
            {"initialization_receipt",receipt},{"reference",J::Numbers(reference)},{"fixed_boundary",J::Constraints(data.fixedBoundary)},
            {"elapsed_seconds",J::Number(elapsed)},{"observation",observation},
            {"resolution",resolution},{"test_basis",testBasis},{"native_constitutive",constitutive},{"physical_scalar_budget",physicalBudget},
            {"tensor_validation",tensorValidation},{"mixed_space",mixedSpace}});
        if(open){output.pop_back();output+=",\"open_boundary\":"+openPacket+",\"open_frozen_systems\":"+openSystems+"}";}
        if(request.get<bool>("export_manufactured",false))
        {
            Require(patch.CurvedManufactured(),"manufactured export requires its verification capability");
            const RANSWallManufactured mms(patch);std::vector<std::string> samples;
            for(int i=0;i<grid.size();++i)
            {
                const auto &q=grid[i];const bool outer=q[2]==p.geometry.height?sides[i]>0:q[2]>p.geometry.height;const auto v=mms.At(q,time,outer);
                samples.push_back(J::Object({{"coordinate",J::Numbers(q)},{"outer",J::Boolean(outer)},
                    {"values",J::Numbers(v.value)},{"gradients",J::Rows(v.gradient)},{"velocity_gradient",J::Rows(v.velocityGradient)},
                    {"velocity_time",J::Numbers(v.velocityTime)},{"k_time",J::Number(v.kTime)},{"epsilon_time",J::Number(v.epsilonTime)},
                    {"required_body",J::Numbers(v.body)},{"k_source",J::Number(v.kSource)},{"epsilon_source",J::Number(v.epsilonSource)},
                    {"R",J::Rows(v.reynolds)},{"Q",J::Rows(v.quadratic)},{"T",J::Rows(v.stress)},
                    {"k_diffusion",J::Number(v.kDiffusion)},{"epsilon_diffusion",J::Number(v.epsilonDiffusion)},
                    {"production",J::Number(v.production)},{"epsilon_production",J::Number(v.epsilonProduction)},{"epsilon_destruction",J::Number(v.epsilonDestruction)},
                    {"K",J::Number(v.K)},{"cutoff",J::Number(v.cutoff)},{"J",J::Number(v.J)},{"a",J::Number(v.a)},{"gamma",J::Number(v.gamma)}}));
            }
            output.pop_back();output+=",\"manufactured\":"+J::Object({{"version",J::Quote("CT6ContinuousV1")},{"samples",J::List(samples)}})+"}";
        }
        if(request.get<bool>("export_neutral_composite_reference",false))
        {
            output.pop_back();output+=",\"neutral_composite_reference\":"+RANSAcceptedCompositeReferencePacket(session,patch,grid,sides)+"}";
        }
        std::ofstream file(argv[3]);Require(bool(file),"cannot write verification output");file<<output<<'\n';
        return physical?0:1;
    }
    catch(const std::exception &error)
    {
        if(argc==4){std::ofstream file(argv[3]);file<<J::Object({{"status",J::Quote("fail")},{"failure",J::Quote(error.what())}})<<'\n';}
        std::cerr<<"preservation verification: "<<error.what()<<std::endl;return 1;
    }
}
