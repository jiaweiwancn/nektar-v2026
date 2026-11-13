// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallManufactured.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>

namespace Nektar
{
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;
bool RANSUsesNeutralComposite(const LibUtilities::SessionReaderSharedPtr &s)
{
    std::string mode;s->LoadSolverInfo("RANSOpenReference",mode,"Functions");
    Require(mode=="Functions" || mode=="NeutralCompositeV1","unknown open reference provider version");
    return mode=="NeutralCompositeV1";
}
namespace
{
// Immutable continuous physical reference, shared by initialization, live
// inlet/top, absolute outlet traction and read-only exports. This never sets
// the evolving wall traction or replaces the current sampled endpoint K.
class NeutralComposite
{
public:
    const RANSWallEvolutionPatch &patch;
    double us,H,Hh,K,A,UH,PiOuter,PiInner=0.;
    std::shared_ptr<SolverUtils::EquilibriumAWFProfiles> profile;
    SolverUtils::RANSTensor Q{};
    mutable std::map<std::pair<double,bool>,std::array<double,9>> cache;
    NeutralComposite(const LibUtilities::SessionReaderSharedPtr &s,const RANSWallEvolutionPatch &p):patch(p)
    {
        const auto &m=p.Parameters().model;const auto &g=p.Parameters().geometry;
        double displacement;s->LoadParameter("RANSOpenDisplacement",displacement,0.);
        Require(displacement==0.,
            "original Accepted/NeutralCompositeV1 cannot silently use aerodynamic displacement; a separately qualified displaced wall policy is required");
        s->LoadParameter("RANSOpenPressureLevel",PiInner,0.);
        Require(std::isfinite(PiInner) && (p.RichardsNorrisTop() || PiInner==0.),
            "a changed absolute pressure/traction level is qualified only for RichardsNorris2011V1");
        if(p.RichardsNorrisTop())for(const std::string name:{"TopTotalTraction","TopDiffusiveFlux"})
            Require(!s->DefinesFunction(name),"RichardsNorris2011V1 forbids competing top functions "+name);
        Require(p.OpenBoundary() && p.NewModels() && p.AcceptedAWF() && !p.Mapped() && p.Dimension()==2 &&
            m.KappaPolicy()==SolverUtils::RANSKappaPolicy::Natural,
            "NeutralCompositeV1 requires original Natural Accepted serial flat2D open route");
        s->LoadParameter("RANSWallFrictionVelocity",us);s->LoadParameter("RANSOpenReferenceLength",Hh);
        H=g.height;K=m.Q()*us*us;A=us*us*us/m.Kappa();
        Require(std::isfinite(us) && us>0 && std::isfinite(Hh) && Hh>0 && g.sampleHeight==H && g.tailStart==H,
            "invalid neutral physical scales or endpoint sampler");
        double q,kappa;s->LoadParameter("RANSOpenEquilibriumQ",q,m.Q());s->LoadParameter("RANSOpenNaturalKappa",kappa,m.Kappa());
        Require(std::abs(q-m.Q())<=2e-13*m.Q() && std::abs(kappa-m.Kappa())<=2e-13*m.Kappa(),
            "NeutralCompositeV1 wrong model equilibrium data");
        for(const std::string name:{"InitialConditionsInner","InitialConditionsOuter","BoundaryValuesInner","BoundaryValuesOuter",
            "OutletTotalTractionInner","OutletTotalTractionOuter","OutletDiffusiveFluxInner","OutletDiffusiveFluxOuter",
            "BodyForceInner","BodyForceOuter","ScalarSourceInner","ScalarSourceOuter"})
            Require(!s->DefinesFunction(name),"NeutralCompositeV1 forbids a competing physical function "+name);
        SolverUtils::AWFProfileOptions o;o.formulation=SolverUtils::AWFProfileFormulation::AcceptedBudgetAWF;o.quadratureOrder=96;
        profile=std::make_shared<SolverUtils::EquilibriumAWFProfiles>(g,K,m,SolverUtils::AWFNormalMetric{},o,SolverUtils::CalibratedAWFParameters(m));
        UH=profile->Reference().matchingVelocity;
        double B;s->LoadParameter("RANSOpenMeanIntercept",B);
        Require(std::abs(UH-(B+us/m.Kappa()*std::log(H/Hh)))<=2e-10*std::max(1.,std::abs(UH)),
            "NeutralCompositeV1 mean-fit and finite-nu wall mismatch");
        SolverUtils::RANSTensor G{};G[0][2]=us/(m.Kappa()*H);
        const auto law=SolverUtils::EvaluateRANSClosure(m.Model(),K,A/H,G,H,m.Parameters(),0);
        Q=law.quadraticStress;
        PiOuter=PiInner==0.?-Q[2][2]:PiInner-Q[2][2];
    }
    std::array<double,9> At(double z,bool outer) const
    {
        // Native FwdTrans samples GetCoords() at an outer-field endpoint.
        // It can differ from the declared H vertex by an arithmetic ulp.
        // Evaluate the analytic outer continuation at the ORIGINAL positive
        // coordinate, never snap coordinates or relax native element ownership.
        const bool endpoint=outer && z<H && H-z<=8*std::numeric_limits<double>::epsilon()*H;
        std::ostringstream where;where<<std::setprecision(17)<<"neutral reference coordinate/side mismatch: z="<<z<<",H="<<H<<",outer="<<outer;
        Require(std::isfinite(z) && z>=0 && (outer?(z>=H || endpoint):z<=H),where.str());
        const auto key=std::make_pair(z,outer);const auto found=cache.find(key);if(found!=cache.end())return found->second;
        if(endpoint)std::cerr<<std::setprecision(17)<<"NeutralCompositeV1: outer-field endpoint analytic continuation at original z="<<z<<", H="<<H<<", delta="<<H-z<<"; no coordinate change\n";
        std::array<double,9> v{}; // U,k,epsilon,Pi,U',k',epsilon',Q0,Dk
        const double kap=patch.Parameters().model.Kappa();
        if(outer)v={UH+us/kap*std::log(z/H),K,A/z,PiOuter,us/(kap*z),0,-A/(z*z),0,kap*us*z};
        else
        {
            const double cutoff=profile->Reference().cutoff;
            v={profile->Original().Velocity(z,{UH,0},{0,0})[0],profile->ReferenceK(z),profile->Dissipation(z),PiInner,
                us*us/profile->Original().TotalViscosity(z),profile->ReferenceKGradient(z),
                z<cutoff?0.:-A/(z*z),profile->ReferenceKFlux(z),profile->ReferenceKDiffusivity(z)};
        }
        cache.emplace(key,v);return v;
    }
    RANSWallEvolutionPatch::Vector Project() const
    {
        using V=RANSWallEvolutionPatch::Vector;const double Z=patch.upper[2];V state(patch.Size());
        const auto at=[&](double z){return At(z,z>H);};
        const auto u=patch.fields[0]->RadialProject([&](double z){return at(z)[0];},
            {{0,0,0},{0,1,At(0,false)[4]},{H,0,UH},{H,-1,At(H,false)[4]},
             {H,1,At(H,true)[4]},{Z,0,At(Z,true)[0]},{Z,-1,At(Z,true)[4]}},
            {profile->Reference().cutoff,profile->Original().ViscosityInterceptHeight(),profile->Geometry().roughness});
        const auto k=patch.fields[patch.K()]->RadialProject([&](double z){return at(z)[1];},
            {{0,1,0},{H,0,K},{H,-1,0},{H,1,0},{Z,0,K},{Z,-1,0}},
            {profile->Reference().cutoff,profile->Original().ViscosityInterceptHeight(),profile->Geometry().roughness});
        // Only an outer epsilon representation. H value is dependent; no
        // independent inner equation or two-sided epsilon-flux condition.
        const auto e=patch.fields[patch.Epsilon()]->RadialProject([&](double z){return A/z;},
            {{H,0,A/H},{Z,0,A/Z},{Z,-1,-A/(Z*Z)}});
        std::copy(u.begin(),u.end(),state.begin()+patch.offset[0]);
        std::copy(k.begin(),k.end(),state.begin()+patch.offset[patch.K()]);
        std::copy(e.begin(),e.end(),state.begin()+patch.offset[patch.Epsilon()]);
        const auto p=patch.fields[patch.OuterPressure()]->Project([&](const RANSWallCoordinate &){return PiOuter;});
        std::copy(p.begin(),p.end(),state.begin()+patch.offset[patch.OuterPressure()]);
        if(PiInner!=0.)
        {
            const auto inner=patch.fields[patch.InnerPressure()]->Project([&](const RANSWallCoordinate &){return PiInner;});
            std::copy(inner.begin(),inner.end(),state.begin()+patch.offset[patch.InnerPressure()]);
        }
        return state;
    }
};
RANSWallEvolutionFunction Load(const LibUtilities::SessionReaderSharedPtr &s,
    const std::string &name,const std::string &variable,bool required)
{
    Require(!required || s->DefinesFunction(name,variable),"missing open-boundary function "+name+"/"+variable);
    if(!s->DefinesFunction(name,variable))return [](const RANSWallCoordinate&,double){return 0.;};
    Require(s->GetFunctionType(name,variable)==LibUtilities::eFunctionTypeExpression,
            "AcceptedOpenFlatV1 requires analytic session functions");
    const auto equation=s->GetFunction(name,variable);
    return [equation](const RANSWallCoordinate &p,double t){return equation->Evaluate(p[0],p[2],0.,t);};
}
}

std::string RANSAcceptedOpenBoundaryIdentity(const LibUtilities::SessionReaderSharedPtr &s)
{
    if(RANSUsesAcceptedManufactured(s))return RANSAcceptedManufacturedIdentity(s);
    std::string policy;s->LoadSolverInfo("RANSWallBoundaryPolicy",policy,"Closed");
    const bool driven=policy=="AcceptedRichardsNorrisFlatV1";
    std::ostringstream out;out<<std::setprecision(17)<<(driven?
        "AcceptedRichardsNorrisFlatV1:top=RichardsNorris2011V1:w=0:t_x=ustar2:g_k=0:g_epsilon=prescribed-negative:backflow=RejectV1:pressure=absolute-total-traction:unused-gauge-zero":
        "AcceptedOpenFlatV1:top=DirichletCompatibleV1:backflow=RejectV1:pressure=absolute-total-traction:unused-gauge-zero");
    if(RANSUsesNeutralComposite(s))
    {
        std::string projection;s->LoadSolverInfo("RANSOpenReferenceProjection",projection,"FunctionsL2V1");
        Require(projection=="FunctionsL2V1" || projection=="ConstrainedRadialV1","unknown neutral reference projection");
        out<<":reference=NeutralCompositeV1:projection="<<projection<<":coordinate=physical-z-unshifted";
        if(driven)
        {
            double level;s->LoadParameter("RANSOpenPressureLevel",level,0.);
            out<<":Pi_inner="<<level<<":Pi_outer=level-Qzz:no-finite-depth-cap";
        }
        else out<<":Pi_inner=0:Pi_outer=-Qzz:no-finite-depth-cap";
    }
    const std::set<std::string> work{"NUMSTEPS","TIMESTEP","TIME","IO_CHECKSTEPS","IO_INFOSTEPS","IO_TIMER_LEVEL",
        "RANSWALLNORMALORDER","RANSWALLTANGENTIALORDER","RANSWALLMASSORDER","RANSWALLVOLUMEORDER","RANSWALLINTERFACEORDER",
        "RANSWALLPIMPLEREJECTSTEP","RANSWALLPIMPLEHARDAUDITREJECTSTEP"};
    for(const auto &[name,value]:s->GetParameters())if(!work.count(name))out<<'|'<<name<<'='<<value;
    for(const std::string name:{"BoundaryValuesInner","BoundaryValuesOuter","OutletTotalTractionInner","OutletTotalTractionOuter",
        "OutletDiffusiveFluxInner","OutletDiffusiveFluxOuter","BodyForceInner","BodyForceOuter",
        "ScalarSourceInner","ScalarSourceOuter"})
        for(const auto &variable:s->GetVariables())if(s->DefinesFunction(name,variable))
            out<<'|'<<name<<'/'<<variable<<'='<<s->GetFunction(name,variable)->GetExpression();
    return out.str();
}

RANSWallEvolutionData LoadRANSAcceptedOpenData(const LibUtilities::SessionReaderSharedPtr &s,
    const RANSWallEvolutionPatch &patch)
{
    Require(patch.OpenBoundary() && patch.Dimension()==2 && (!patch.Mapped() || patch.CurvedManufactured()),"open data loader requires qualified Accepted policy");
    if(patch.CurvedManufactured())return RANSWallManufactured(patch).Data();
    RANSWallEvolutionData data(2);data.openBoundaryIdentity=RANSAcceptedOpenBoundaryIdentity(s);
    if(RANSUsesNeutralComposite(s))
    {
        const auto ref=std::make_shared<NeutralComposite>(s,patch);
        data.boundary[0]=[ref](const RANSWallCoordinate &p,double){return ref->At(p[2],p[2]>ref->H)[0];};
        data.boundary[1]=[](const RANSWallCoordinate&,double){return 0.;};
        data.boundary[patch.K()]=[ref](const RANSWallCoordinate &p,double){return ref->At(p[2],p[2]>ref->H)[1];};
        data.boundary[patch.Epsilon()]=[ref](const RANSWallCoordinate &p,double){return ref->At(p[2],true)[2];};
        std::string projection;s->LoadSolverInfo("RANSOpenReferenceProjection",projection,"FunctionsL2V1");
        if(projection=="ConstrainedRadialV1")
        {
            data.projectedReference=ref->Project();
            data.projectedReferenceBoundary=patch.FixedReferenceBoundary(data.projectedReference);
            // Apply the existing dependent-trace initialization once before
            // defining the immutable live corner rows. This records the same
            // numerical corner in initialization and replay, including ulps.
            // No such reference operation occurs after advancing k.
            data.projectedReference=patch.Lift(data.projectedReference,0.,data);
            data.projectedReferenceBoundary=patch.FixedReferenceBoundary(data.projectedReference);
        }
        for(int region=0;region<2;++region)
        {
            data.outletTraction[region][0]=[ref,region](const RANSWallCoordinate&,double){
                const double base=region?ref->Q[2][2]-ref->Q[0][0]:0.;return ref->PiInner==0.?base:base-ref->PiInner;};
            data.outletTraction[region][1]=[ref](const RANSWallCoordinate&,double){return ref->us*ref->us;};
        }
        if(patch.RichardsNorrisTop())
        {
            const double shear=ref->us*ref->us;
            const double sigma=patch.Parameters().model.Parameters().sigmaEpsilon;
            const double flux=-shear*shear/(sigma*patch.upper[2]);
            data.topTraction[0]=[shear](const RANSWallCoordinate&,double){return shear;};
            data.topEpsilonDiffusion=[flux](const RANSWallCoordinate&,double){return flux;};
            // Physical flux data: deliberately independent of evolving nu_t.
            // Normal velocity is essential; no normal traction/pressure datum.
        }
        return data; // default zero body/source/diffusive data; advection stays assembled
    }
    const double H=patch.Parameters().geometry.height;
    for(int b=0;b<patch.InnerPressure();++b)
    {
        const std::string variable=b<2?(b==0?"u":"v"):(b==patch.K()?"k":"eps_tilde");
        const auto outer=Load(s,"BoundaryValuesOuter",variable,true);
        if(b==patch.Epsilon())data.boundary[b]=outer;
        else
        {
            const auto inner=Load(s,"BoundaryValuesInner",variable,true);
            data.boundary[b]=[inner,outer,H](const RANSWallCoordinate &p,double t){return p[2]<=H?inner(p,t):outer(p,t);};
        }
    }
    for(int region=0;region<2;++region)
    {
        const std::string suffix=region?"Outer":"Inner";
        for(int a=0;a<2;++a)data.outletTraction[region][a]=Load(s,"OutletTotalTraction"+suffix,a==0?"u":"v",true);
        data.outletKDiffusion[region]=Load(s,"OutletDiffusiveFlux"+suffix,"k",true);
    }
    data.outletEpsilonDiffusion=Load(s,"OutletDiffusiveFluxOuter","eps_tilde",true);
    for(int a=0;a<2;++a)
    {
        const std::string variable=a==0?"u":"v";
        data.innerBody[a]=Load(s,"BodyForceInner",variable,false);
        data.outerBody[a]=Load(s,"BodyForceOuter",variable,false);
    }
    data.innerKSource=Load(s,"ScalarSourceInner","k",false);
    data.outerKSource=Load(s,"ScalarSourceOuter","k",false);
    data.outerEpsilonSource=Load(s,"ScalarSourceOuter","eps_tilde",false);
    return data;
}

std::vector<RANSWallEvolutionFunction> RANSAcceptedOpenInitialFunctions(
    const LibUtilities::SessionReaderSharedPtr &s,const RANSWallEvolutionPatch &patch)
{
    Require(patch.OpenBoundary(),"open initial-function loader requires its policy");
    if(patch.CurvedManufactured())return RANSWallManufactured(patch).InitialFunctions();
    std::vector<RANSWallEvolutionFunction> functions;
    if(RANSUsesNeutralComposite(s))
    {
        const auto ref=std::make_shared<NeutralComposite>(s,patch);
        for(int b=0;b<patch.Gauge();++b)functions.push_back([ref,b](const RANSWallCoordinate &p,double){
            if(b==1)return 0.;
            if(b==ref->patch.InnerPressure())return ref->PiInner;
            if(b==ref->patch.OuterPressure())return ref->PiOuter;
            const bool outer=b==ref->patch.Epsilon() || b==ref->patch.OuterPressure() || p[2]>ref->H;
            const auto v=ref->At(p[2],outer);
            return b==0?v[0]:b==1?0.:b==ref->patch.K()?v[1]:b==ref->patch.Epsilon()?v[2]:v[3];});
        return functions;
    }
    const double H=patch.Parameters().geometry.height;
    for(int b=0;b<patch.Gauge();++b)
    {
        const std::string variable=b<2?(b==0?"u":"v"):(b==patch.K()?"k":b==patch.Epsilon()?"eps_tilde":"p");
        if(b==patch.Epsilon() || b==patch.OuterPressure())functions.push_back(Load(s,"InitialConditionsOuter",variable,true));
        else if(b==patch.InnerPressure())functions.push_back(Load(s,"InitialConditionsInner",variable,true));
        else
        {
            const auto inner=Load(s,"InitialConditionsInner",variable,true),outer=Load(s,"InitialConditionsOuter",variable,true);
            functions.push_back([inner,outer,H](const RANSWallCoordinate &p,double t){return p[2]<=H?inner(p,t):outer(p,t);});
        }
    }
    return functions;
}

std::string RANSAcceptedCompositeReferencePacket(const LibUtilities::SessionReaderSharedPtr &s,
    const RANSWallEvolutionPatch &patch,const std::vector<RANSWallCoordinate> &grid,const std::vector<int> &sides)
{
    Require(RANSUsesNeutralComposite(s) && grid.size()==sides.size(),"missing neutral provider or side inventory");
    const NeutralComposite ref(s,patch);std::vector<std::string> rows;
    for(std::size_t i=0;i<grid.size();++i)
    {
        const double z=grid[i][2];Require(z!=ref.H || sides[i]!=0,"ambiguous exact-interface neutral reference request");
        const bool outer=z>ref.H || (z==ref.H && sides[i]>0);const auto v=ref.At(z,outer);
        rows.push_back(J::Object({{"coordinate",J::Numbers(grid[i])},{"side",J::Number(sides[i])},
            {"region",J::Number(outer?1:0)},{"values_U_k_epsilon_Pi_gradients_Q0_D",J::Numbers(v)}}));
    }
    return J::Object({{"provider",J::Quote("NeutralCompositeV1")},{"identity",J::Quote(patch.OpenBoundaryIdentity())},
        {"K",J::Number(ref.K)},{"A",J::Number(ref.A)},{"UH",J::Number(ref.UH)},
        {"cutoff",J::Number(ref.profile->Reference().cutoff)},{"Q",J::Rows(ref.Q)},{"samples",J::List(rows)}});
}

RANSWallEvolutionPatch::OpenBoundaryTerms RANSWallEvolutionPatch::OpenBoundaryContributions(
    const Vector &state,double time,const RANSWallEvolutionData &data,int selectedRegion) const
{
    RANSWallTiming::Scope timing("assembly.open_faces",true,true);
    Require(OpenBoundary() && state.size()==Size() && selectedRegion>=-1 && selectedRegion<=1,
            "invalid open-boundary contribution request");
    Require(data.openBoundaryIdentity==m_openBoundaryIdentity,"open live boundary definition differs from its immutable policy identity");
    OpenBoundaryTerms result{Vector(Size()),Vector(Size()),Vector(Size())};
    auto add=[&](Vector &load,int block,const RANSWallGlobalPoint &q,double value,double weight)
    {
        Require(std::isfinite(value),"nonfinite open-boundary physical datum or flux");
        for(int j=0;j<q.ids.size();++j)load[offset[block]+q.ids[j]]+=weight*q.phi[j]*value;
    };
    for(const auto &face:Surface(0,0,upper[0],1))
    {
        const auto &q=face.test;const auto &element=fields[0]->Elements()[q.element];
        if(Owner(element.geometry)!=m_comm->GetRank())continue;
        const int region=element.lower[2]>=m_parameters.geometry.height?1:0;
        if(selectedRegion>=0 && selectedRegion!=region)continue;
        SolverUtils::RANSVector u{};for(int a=0;a<2;++a)u[PhysicalAxis(a)]=q.Value(state,offset[a]);
        const double un=u[0]*face.area[0]+u[2]*face.area[2];
        const double area=std::hypot(face.area[0],face.area[2]);
        for(int a=0;a<2;++a)
        {
            add(result.advection,a,q,u[PhysicalAxis(a)]*un,face.weight);
            add(result.traction,a,q,-area*data.outletTraction[region][a](q.coordinate,time),face.weight);
        }
        add(result.advection,K(),q,un*q.Value(state,offset[K()]),face.weight);
        add(result.diffusion,K(),q,-area*data.outletKDiffusion[region](q.coordinate,time),face.weight);
        if(region)
        {
            const auto e=GeometryPoint(Epsilon(),q.geometry,q.coordinate);
            add(result.advection,Epsilon(),e,un*e.Value(state,offset[Epsilon()]),face.weight);
            add(result.diffusion,Epsilon(),e,-area*data.outletEpsilonDiffusion(q.coordinate,time),face.weight);
        }
    }
    if(RichardsNorrisTop() && selectedRegion!=0)
        for(const auto &face:Surface(0,2,upper[2],1))
        {
            const auto &v=face.test;if(Owner(v.geometry)!=m_comm->GetRank())continue;
            const auto e=GeometryPoint(Epsilon(),v.geometry,v.coordinate);
            const double un=v.Value(state,offset[0])*face.area[0]+v.Value(state,offset[1])*face.area[2];
            const double area=std::hypot(face.area[0],face.area[2]);
            for(int a=0;a<2;++a)add(result.advection,a,v,un*v.Value(state,offset[a]),face.weight);
            add(result.traction,0,v,-area*data.topTraction[0](v.coordinate,time),face.weight);
            add(result.advection,K(),v,un*v.Value(state,offset[K()]),face.weight);
            add(result.advection,Epsilon(),e,un*e.Value(state,offset[Epsilon()]),face.weight);
            add(result.diffusion,K(),v,-area*data.topKDiffusion(v.coordinate,time),face.weight);
            add(result.diffusion,Epsilon(),e,-area*data.topEpsilonDiffusion(v.coordinate,time),face.weight);
        }
    return result;
}

std::string RANSWallEvolutionPatch::OpenBoundaryStatus(const Vector &state,double time,
    const RANSWallEvolutionData &data,bool &admissible) const
{
    Require(OpenBoundary() && data.openBoundaryIdentity==m_openBoundaryIdentity,"invalid open-boundary audit identity");
    double minimum=1e300;RANSWallCoordinate point{};int geometry=-1,count=0;
    for(const auto &face:Surface(0,0,upper[0],1))
    {
        const double un=face.test.Value(state,offset[0])*face.area[0]+
            face.test.Value(state,offset[1])*face.area[2];++count;
        Require(std::isfinite(un),"nonfinite outlet normal velocity");
        if(un<minimum){minimum=un;point=face.test.coordinate;geometry=face.test.geometry;}
    }
    const double tolerance=1e-12*m_openVelocityScale;
    admissible=minimum>=-tolerance;
    return J::Object({{"policy",J::Quote(m_boundaryPolicy)},{"admissible",J::Boolean(admissible)},
        {"minimum_outward_velocity",J::Number(minimum)},{"roundoff_tolerance",J::Number(tolerance)},
        {"point",J::Numbers(point)},{"geometry",J::Number(geometry)},{"side","1"},{"time",J::Number(time)},
        {"samples",J::Number(count)},{"reason",J::Quote(admissible?"outflow":"outlet_backflow_unsupported")},
        {"scope",J::Quote("declared Gaussian outlet face inventory; signed flux retained; no rejection of internal reversed flow or negative wall shear")}});
}

std::string RANSWallEvolutionPatch::OpenBoundaryPacket(const Vector &state,double time,
    const RANSWallEvolutionData &data) const
{
    Require(OpenBoundary(),"open-boundary packet requires its explicit policy");
    const auto terms=OpenBoundaryContributions(state,time,data);std::vector<std::string> faces,top;
    for(const auto &face:Surface(0,0,upper[0],1))
    {
        const auto &q=face.test;const int region=fields[0]->Elements()[q.element].lower[2]>=m_parameters.geometry.height?1:0;
        std::vector<double> velocity,traction;
        for(int a=0;a<2;++a){velocity.push_back(q.Value(state,offset[a]));traction.push_back(data.outletTraction[region][a](q.coordinate,time));}
        auto row=[&](const RANSWallGlobalPoint &p){return J::Object({{"ids",J::Numbers(p.ids)},{"phi",J::Numbers(p.phi)},{"derivative",J::Rows(p.derivative)},{"geometry",J::Number(p.geometry)}});};
        const auto e=region?GeometryPoint(Epsilon(),q.geometry,q.coordinate):RANSWallGlobalPoint{};
        faces.push_back(J::Object({{"coordinate",J::Numbers(q.coordinate)},{"area",J::Numbers(face.area)},
            {"weight",J::Number(face.weight)},{"region",J::Number(region)},{"velocity",J::Numbers(velocity)},
            {"traction",J::Numbers(traction)},{"k",J::Number(q.Value(state,offset[K()]))},
            {"epsilon",region?J::Number(e.Value(state,offset[Epsilon()])):"null"},
            {"k_diffusive_datum",J::Number(data.outletKDiffusion[region](q.coordinate,time))},
            {"epsilon_diffusive_datum",region?J::Number(data.outletEpsilonDiffusion(q.coordinate,time)):"null"},
            {"velocity_k_row",row(q)},{"epsilon_row",region?row(e):"null"}}));
    }
    if(RichardsNorrisTop())for(const auto &face:Surface(0,2,upper[2],1))
    {
        const auto &v=face.test;const auto e=GeometryPoint(Epsilon(),v.geometry,v.coordinate);
        auto row=[](const RANSWallGlobalPoint &p){return J::Object({{"ids",J::Numbers(p.ids)},{"phi",J::Numbers(p.phi)},{"derivative",J::Rows(p.derivative)},{"geometry",J::Number(p.geometry)}});};
        top.push_back(J::Object({{"coordinate",J::Numbers(v.coordinate)},{"area",J::Numbers(face.area)},
            {"weight",J::Number(face.weight)},{"velocity_k_row",row(v)},{"epsilon_row",row(e)},
            {"tangential_total_traction",J::Number(data.topTraction[0](v.coordinate,time))},
            {"k_diffusive_datum",J::Number(data.topKDiffusion(v.coordinate,time))},
            {"epsilon_diffusive_datum",J::Number(data.topEpsilonDiffusion(v.coordinate,time))},
            {"normal_velocity_essential","true"},{"independent_normal_traction","false"}}));
    }
    bool admissible=false;const auto status=OpenBoundaryStatus(state,time,data,admissible);
    return J::Object({{"schema",J::Quote("AcceptedOpenBoundaryPacketV1")},{"definition",J::Quote(m_openBoundaryIdentity)},
        {"pressure_level",J::Quote("absolute total traction; physical pressure rows free; unused gauge fixed zero")},
        {"advection_load",J::Numbers(terms.advection)},{"traction_load",J::Numbers(terms.traction)},
        {"diffusive_load",J::Numbers(terms.diffusion)},{"faces",J::List(faces)},{"top_faces",J::List(top)},{"outlet_status",status}});
}
}
