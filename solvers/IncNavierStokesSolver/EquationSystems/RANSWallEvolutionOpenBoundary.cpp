// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <cmath>
#include <iomanip>
#include <sstream>

namespace Nektar
{
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;
namespace
{
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
    std::ostringstream out;out<<std::setprecision(17)<<"AcceptedOpenFlatV1:top=DirichletCompatibleV1:backflow=RejectV1:pressure=absolute-total-traction:unused-gauge-zero";
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
    Require(patch.OpenBoundary() && patch.Dimension()==2 && !patch.Mapped(),"open data loader requires AcceptedOpenFlatV1");
    RANSWallEvolutionData data(2);data.openBoundaryIdentity=RANSAcceptedOpenBoundaryIdentity(s);
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
    std::vector<RANSWallEvolutionFunction> functions;
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
    const auto terms=OpenBoundaryContributions(state,time,data);std::vector<std::string> faces;
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
    bool admissible=false;const auto status=OpenBoundaryStatus(state,time,data,admissible);
    return J::Object({{"schema",J::Quote("AcceptedOpenBoundaryPacketV1")},{"definition",J::Quote(m_openBoundaryIdentity)},
        {"pressure_level",J::Quote("absolute total traction; physical pressure rows free; unused gauge fixed zero")},
        {"advection_load",J::Numbers(terms.advection)},{"traction_load",J::Numbers(terms.traction)},
        {"diffusive_load",J::Numbers(terms.diffusion)},{"faces",J::List(faces)},{"outlet_status",status}});
}
}
