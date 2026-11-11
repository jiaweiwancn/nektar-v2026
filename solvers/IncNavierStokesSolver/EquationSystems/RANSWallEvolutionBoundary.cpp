// SPDX-License-Identifier: MIT
// Physical, conservative exterior terms for the bounded serial WJ flat tunnel.
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace Nektar
{
using namespace SolverUtils;
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;
namespace
{
double Dot(const RANSVector &a,const RANSVector &b)
{double r=0;for(int i=0;i<3;++i)r+=a[i]*b[i];return r;}
bool Same(double a,double b){return std::abs(a-b)<1e-11;}
struct Values
{
    RANSVector u{},normal{},gk{},ge{},traction{},applied{};
    RANSTensor G{},R{},T{};
    double k=0,e=0,pi=0,un=0,dk=0,de=0,weight=0,nt=0,K=0,rc=0;
    RANSWallGlobalPoint epsilon,pressure;
};
Values At(const RANSWallEvolutionPatch &p,const RANSWallEvolutionPatch::ExteriorPoint &f,
    const Vector &s,double time,const RANSWallEvolutionData &data)
{
    Values v;const auto &q=f.face.test;const auto a=q.coordinate;
    const double length=std::hypot(f.face.area[0],f.face.area[1],f.face.area[2]);
    Require(length>0,"open_boundary_invalid_surface");
    for(int d=0;d<3;++d)v.normal[d]=f.face.area[d]/length;
    v.weight=f.face.weight*length;
    for(int b=0;b<p.Dimension();++b)
    {const int d=p.PhysicalAxis(b);v.u[d]=q.Value(s,p.offset[b]);v.G[d]=q.Gradient(s,p.offset[b]);}
    v.k=q.Value(s,p.offset[p.K()]);v.gk=q.Gradient(s,p.offset[p.K()]);
    v.pressure=p.Point(f.inner?p.InnerPressure():p.OuterPressure(),a);
    v.pi=v.pressure.Value(s,p.offset[f.inner?p.InnerPressure():p.OuterPressure()]);
    if(f.inner)
    {
        const auto profile=p.StateProfile(s,a,time,data);v.K=profile->K();v.rc=profile->Reference().cutoff;
        v.nt=profile->Original().EddyViscosity(a[2]);v.e=profile->Dissipation(a[2]);
        const auto law=p.InnerStress(v.k,v.nt,v.G);v.R=law.reynoldsStress;
        v.dk=p.Parameters().geometry.molecularNu+v.nt/p.SigmaK();
        v.T=RANSWallInnerMeanStress(law,p.Parameters().geometry.molecularNu,v.pi);
    }
    else
    {
        v.epsilon=p.Point(p.Epsilon(),a);v.e=v.epsilon.Value(s,p.offset[p.Epsilon()]);
        v.ge=v.epsilon.Gradient(s,p.offset[p.Epsilon()]);const auto law=p.OuterClosure(v.k,v.e,v.G,a);
        v.R=law.reynoldsStress;v.dk=law.kDiffusion;v.de=law.epsilonDiffusion;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)v.T[i][j]=-v.R[i][j]+(i==j?2*v.k/3-v.pi:0);
    }
    for(int i=0;i<3;++i)v.traction[i]=Dot(v.T[i],v.normal);
    v.applied=v.traction;v.un=Dot(v.u,v.normal);
    if(f.role=="outlet")for(int b=0;b<p.Dimension();++b)
        v.applied[p.PhysicalAxis(b)]=(f.inner?data.outletInner:data.outletOuter)[b](a,time);
    if(f.role=="top")
    {
        RANSVector given{};for(int b=0;b<p.Dimension();++b)given[p.PhysicalAxis(b)]=data.topTraction[b](a,time);
        const double normal=Dot(v.traction,v.normal),remove=Dot(given,v.normal);
        for(int i=0;i<3;++i)v.applied[i]=given[i]+(normal-remove)*v.normal[i];
    }
    for(double x:v.applied)Require(std::isfinite(x),"open_boundary_nonfinite_traction");
    return v;
}
Vector Flux(const Values &v,const RANSWallEvolutionPatch::ExteriorPoint &f,
    double time,const RANSWallEvolutionData &data)
{
    double jk=0,je=0;
    if(f.role=="inlet"){jk=-v.dk*Dot(v.gk,v.normal);if(!f.inner)je=-v.de*Dot(v.ge,v.normal);}
    if(f.role=="top"){jk=data.topKFlux(f.face.test.coordinate,time);je=data.topEpsilonFlux(f.face.test.coordinate,time);}
    return {v.u[0]*v.un-v.applied[0],v.u[2]*v.un-v.applied[2],v.k*v.un+jk,v.e*v.un+je};
}
}

void RANSWallEvolutionPatch::LoadOpenBoundaryData(const LibUtilities::SessionReaderSharedPtr &session,RANSWallEvolutionData &data) const
{
    if(!OpenBoundary())return;
    auto function=[&](const std::string &name,const std::string &variable,bool required)
    {
        Require(!required || session->DefinesFunction(name,variable),"missing open boundary function "+name+"/"+variable);
        if(!session->DefinesFunction(name,variable))return RANSWallEvolutionFunction([](const RANSWallCoordinate&,double){return 0.;});
        Require(session->GetFunctionType(name,variable)==LibUtilities::eFunctionTypeExpression,"open boundaries require expression Functions");
        auto equation=session->GetFunction(name,variable);
        return RANSWallEvolutionFunction([this,equation](const RANSWallCoordinate &q,double t){const auto p=Physical(q);return equation->Evaluate(p[0],p[2],0.,t);});
    };
    for(int b=0;b<InnerPressure();++b)
    {
        const std::string variable=b==0?"u":b==1?"v":b==K()?"k":"epsilon";
        const auto inner=function("BoundaryValuesInner",variable,b!=Epsilon());
        const auto outer=function("BoundaryValuesOuter",variable,true);
        data.boundary[b]=[inner,outer,b,ep=Epsilon(),H=m_parameters.geometry.height](const RANSWallCoordinate &q,double t){
            if(b!=ep && Same(q[2],H))Require(std::abs(inner(q,t)-outer(q,t))<1e-10*std::max(1.,std::abs(outer(q,t))),
                "open_inlet_shared_value_corner_conflict");
            return q[2]<H?inner(q,t):outer(q,t);};
    }
    for(int b=0;b<m_dimension;++b)
    {
        const std::string variable=b==0?"u":"v";
        data.outletInner[b]=function("OutletTractionInner",variable,true);
        data.outletOuter[b]=function("OutletTractionOuter",variable,true);
    }
    Require(!session->DefinesFunction("BoundaryValuesInner","p") && !session->DefinesFunction("BoundaryValuesOuter","p"),
        "open outlet fixes mean traction; an extra pressure Dirichlet datum is unsupported");
    if(m_openTop=="ReferenceTractionV1")
    {
        Require(!session->DefinesFunction("TopTractionOuter","v"),"normal top traction is a reaction, not extra prescribed data");
        data.topTraction[0]=function("TopTractionOuter","u",true);
        data.topKFlux=function("TopDiffusiveFluxOuter","k",true);
        data.topEpsilonFlux=function("TopDiffusiveFluxOuter","epsilon",true);
    }
    else Require(!session->DefinesFunction("TopTractionOuter") && !session->DefinesFunction("TopDiffusiveFluxOuter"),
        "SlipV1 specifies zero tangential traction and zero scalar flux; conflicting top Functions");
    data.Touch();
}

std::map<int,double> RANSWallEvolutionPatch::OpenConstraints(double time,const RANSWallEvolutionData &data) const
{
    Require(OpenBoundary() && data.fixedBoundary.empty(),"open boundaries cannot use fixed closed-reference traces");
    auto corner=lower;corner[2]=m_parameters.geometry.height;
    const double K=data.boundary[this->K()](corner,time),epsilon=data.boundary[Epsilon()](corner,time);
    Require(K>0 && std::isfinite(K),"open_inlet_nonpositive_K");
    const double target=std::pow(K/m_parameters.model.Q(),1.5)/(m_parameters.model.Kappa()*m_parameters.geometry.height);
    Require(std::isfinite(epsilon) && std::abs(epsilon-target)<=1e-10*target,
        "open_inlet_epsilon_corner_conflict: matching epsilon belongs to A(K)/H");
    std::map<int,double> rows;
    for(int b=0;b<InnerPressure();++b)
    {
        const auto values=fields[b]->Boundary([&](const RANSWallCoordinate &q)
        {return b<m_dimension && (Same(q[2],0) || (b==1 && Same(q[2],upper[2])))?0.:data.boundary[b](q,time);},
        [&](int axis,int side,const RANSWallGlobalField::Element &el)
        {
            if(axis==0)return !side && Same(el.lower[0],lower[0]);
            return axis==2 && b<m_dimension && ((!side && Same(el.lower[2],0)) ||
                (side && b==1 && Same(el.upper[2],upper[2])));
        });
        for(auto [id,value]:values)if(b!=Epsilon() || !m_traceIds.count(id))rows[offset[b]+id]=value;
    }
    m_boundaryCorrection=0.;return rows;
}

std::vector<RANSWallEvolutionPatch::ExteriorPoint> RANSWallEvolutionPatch::Exterior(int order) const
{
    std::vector<ExteriorPoint> rows;if(!OpenBoundary())return rows;
    if(!order)order=std::max(m_parameters.normalOrder,m_parameters.tangentialOrder);
    for(const auto &[axis,position,side,role]:std::vector<std::tuple<int,double,int,std::string>>{
        {0,lower[0],-1,"inlet"},{0,upper[0],1,"outlet"},{2,0.,-1,"ground"},{2,upper[2],1,"top"}})
        for(const auto &f:Surface(0,axis,position,side,order))
            rows.push_back({f,role,f.test.coordinate[2]<m_parameters.geometry.height});
    return rows;
}

void RANSWallEvolutionPatch::ValidateOutlet(const Vector &s,const std::string &phase) const
{
    if(!OpenBoundary())return;
    for(const auto &f:Exterior())if(f.role=="outlet")
    {
        RANSVector u{};for(int b=0;b<m_dimension;++b)u[PhysicalAxis(b)]=f.face.test.Value(s,offset[b]);
        const double un=Dot(u,f.face.area)/std::hypot(f.face.area[0],f.face.area[1],f.face.area[2]);
        if(!std::isfinite(un) || un < -1e-12)
        {
            m_solveInvalidContext=J::Object({{"operation",J::Quote("open_boundary_backflow")},{"phase",J::Quote(phase)},
                {"physical",J::Numbers(Physical(f.face.test.coordinate))},{"geometry",J::Number(f.face.test.geometry)},
                {"normal_velocity",J::Number(un)}});
            throw std::runtime_error("open_boundary_backflow: "+m_solveInvalidContext);
        }
    }
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::ExteriorResidual(const Vector &s,double time,
    const RANSWallEvolutionData &data,int region) const
{
    Vector r(Size());if(!OpenBoundary())return r;
    for(const auto &f:Exterior())if(region<0 || int(!f.inner)==region)
    {
        const auto v=At(*this,f,s,time,data);const auto flux=Flux(v,f,time,data);
        for(int b=0;b<InnerPressure();++b)
        {
            if(f.inner && b==Epsilon())continue;
            const auto &q=b==Epsilon()?v.epsilon:f.face.test;
            for(int j=0;j<q.ids.size();++j)r[offset[b]+q.ids[j]]+=v.weight*q.phi[j]*flux[b];
        }
    }
    return r;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::ExteriorPressureAction(const Vector &direction,int region) const
{
    Vector r(Size());if(!OpenBoundary())return r;
    for(const auto &f:Exterior())if(f.role!="outlet" && (region<0 || int(!f.inner)==region))
    {
        const auto &q=f.face.test;const auto p=Point(f.inner?InnerPressure():OuterPressure(),q.coordinate);
        const double pressure=p.Value(direction,offset[f.inner?InnerPressure():OuterPressure()]);
        for(int b=0;b<m_dimension;++b)for(int j=0;j<q.ids.size();++j)
            r[offset[b]+q.ids[j]]+=f.face.weight*q.phi[j]*pressure*f.face.area[PhysicalAxis(b)];
    }
    return r;
}

void RANSWallEvolutionPatch::ExteriorPicard(int block,const Vector &s,double time,const RANSWallEvolutionData &data,
    PicardBlock &system,const std::vector<int> &index) const
{
    if(!OpenBoundary())return;
    const auto residual=ExteriorResidual(s,time,data);
    for(int i=0;i<system.rows.size();++i)system.residual[i]+=residual[system.rows[i]];
    for(const auto &f:Exterior())
    {
        if(block==Epsilon() && f.inner)continue;
        const auto v=At(*this,f,s,time,data);const auto &q=block==Epsilon()?v.epsilon:f.face.test;
        for(int b=block==0?0:block;b<(block==0?m_dimension:block+1);++b)
            for(int i=0;i<q.ids.size();++i)
            {
                const int row=index[offset[b]+q.ids[i]];if(row<0)continue;
                for(int j=0;j<q.ids.size();++j)
                {
                    const int col=index[offset[b]+q.ids[j]];if(col<0)continue;
                    double value=v.un*q.phi[j];
                    if(block!=0 && f.role=="inlet")
                        for(int d=0;d<3;++d)value-=(block==K()?v.dk:v.de)*v.normal[d]*q.derivative[d][j];
                    system.matrix[row][col]+=v.weight*q.phi[i]*value;
                }
                if(block==0 && f.role!="outlet")
                {
                    const int pb=f.inner?InnerPressure():OuterPressure();
                    for(int j=0;j<v.pressure.ids.size();++j)
                        system.matrix[row][index[offset[pb]+v.pressure.ids[j]]]+=
                            v.weight*q.phi[i]*v.normal[PhysicalAxis(b)]*v.pressure.phi[j];
                }
            }
    }
}

std::string RANSWallEvolutionPatch::OpenBoundaryDiagnostic(const Vector &s,const Vector &history,
    double alpha,double time,const RANSWallEvolutionData &data,
    const std::vector<RANSWallCoordinate> &physicalLocations) const
{
    if(!OpenBoundary())return "null";
    const auto hash=RANSWallR3::Hash(s);const auto work=InitializationWork();
    std::vector<std::string> rows,stations;double totalMass=0;
    for(const auto &f:Exterior())
    {
        const auto v=At(*this,f,s,time,data);const auto flux=Flux(v,f,time,data);totalMass+=v.weight*v.un;
        rows.push_back(J::Object({{"role",J::Quote(f.role)},{"inner",J::Boolean(f.inner)},
            {"chart",J::Numbers(f.face.test.coordinate)},{"physical",J::Numbers(Physical(f.face.test.coordinate))},
            {"geometry",J::Number(f.face.test.geometry)},{"element",J::Number(f.face.test.element)},
            {"weight",J::Number(v.weight)},{"normal",J::Numbers(v.normal)},{"U",J::Numbers(v.u)},
            {"G",J::Rows(v.G)},{"R",J::Rows(v.R)},{"k",J::Number(v.k)},{"epsilon",J::Number(v.e)},
            {"pi",J::Number(v.pi)},{"grad_k",J::Numbers(v.gk)},{"grad_epsilon",J::Numbers(v.ge)},
            {"K",J::Number(v.K)},{"cutoff",J::Number(v.rc)},{"actual_traction",J::Numbers(v.traction)},
            {"applied_traction",J::Numbers(v.applied)},{"flux",J::Numbers(flux)},
            {"ids",J::Numbers(f.face.test.ids)},{"phi",J::Numbers(f.face.test.phi)},
            {"epsilon_ids",J::Numbers(v.epsilon.ids)},{"epsilon_phi",J::Numbers(v.epsilon.phi)}}));
    }
    const double H=m_parameters.geometry.height;
    std::set<double> heights{Mapping().displacement,Mapping().displacement+H,Mapping().tunnelTop};
    for(int i=0;i<=160;++i)heights.insert(Mapping().displacement+(Mapping().tunnelTop-Mapping().displacement)*i/160.);
    std::vector<RANSWallCoordinate> locations=physicalLocations;
    if(locations.empty())for(double fraction:{0.,.25,.5,.75,1.})for(double z:heights)
        locations.push_back({lower[0]+fraction*(upper[0]-lower[0]),0,z});
    for(const auto &physicalPoint:locations)
    {
        Require(physicalPoint[0]>=lower[0] && physicalPoint[0]<=upper[0] && physicalPoint[1]==0 &&
            physicalPoint[2]>=Mapping().displacement && physicalPoint[2]<=Mapping().tunnelTop,"open observation outside physical tunnel");
        const auto q=Mapping().Inverse(physicalPoint);
        for(int side:std::abs(q[2]-H)<1e-11?std::vector<int>{-1,1}:std::vector<int>{0})
        {
            const auto v=Point(0,q,side);const bool inner=q[2]<H-1e-11 || (std::abs(q[2]-H)<1e-11 && side<0);
            RANSVector u{};RANSTensor G{},R{};
            for(int b=0;b<m_dimension;++b){u[PhysicalAxis(b)]=v.Value(s,offset[b]);G[PhysicalAxis(b)]=v.Gradient(s,offset[b]);}
            const double k=v.Value(s,offset[K()]);double epsilon=0,cutoff=0,K=0,P=0,P0=0,E=0,nt=0,nud=0;
            auto anchor=q;anchor[2]=H;const auto profile=StateProfile(s,anchor,time,data);K=profile->K();cutoff=profile->Reference().cutoff;
            if(inner){epsilon=profile->Dissipation(q[2]);nt=profile->Original().EddyViscosity(q[2]);
                const auto law=InnerStress(k,nt,G);R=law.reynoldsStress;P=law.production;P0=profile->ReferenceProduction(q[2]);E=epsilon;}
            else{epsilon=GeometryPoint(Epsilon(),v.geometry,q).Value(s,offset[Epsilon()]);const auto law=OuterClosure(k,epsilon,G,q);
                R=law.reynoldsStress;P=law.production;E=epsilon;nt=law.eddyViscosity;nud=law.scalarViscosity;}
            const auto pressure=GeometryPoint(inner?InnerPressure():OuterPressure(),v.geometry,q);
            const double pi=pressure.Value(s,offset[inner?InnerPressure():OuterPressure()]);
            stations.push_back(J::Object({{"physical",J::Numbers(Physical(q))},{"chart",J::Numbers(q)},
                {"inner",J::Boolean(inner)},{"normal_side",J::Number(side)},{"geometry",J::Number(v.geometry)},
                {"element",J::Number(v.element)},{"U",J::Numbers(u)},
                {"G",J::Rows(G)},{"R",J::Rows(R)},{"covariance_eigenvalues",J::Numbers(J::Eigenvalues(R))},
                {"k",J::Number(k)},{"epsilon",J::Number(epsilon)},{"pi",J::Number(pi)},
                {"mechanical_pressure",J::Number(pi-2*k/3)},{"actual_production",J::Number(P)},
                {"counterfactual_production",inner?J::Number(P0):"null"},{"destruction",J::Number(E)},
                {"nu_eff",J::Number(nt)},{"nu_d",inner?"null":J::Number(nud)},
                {"K",J::Number(K)},{"cutoff",J::Number(cutoff)}}));
        }
    }
    std::vector<std::string> interfaces;
    for(const auto &f:Surface(0,2,H,-1,8))
    {
        const auto q=f.test.coordinate;const auto vi=Point(0,q,-1),vo=Point(0,q,1);
        const auto profile=StateProfile(s,q,time,data);const auto n=Mapping().At(q).normal;
        RANSVector ui{},uo{};RANSTensor gi{},go{};
        for(int b=0;b<m_dimension;++b){const int d=PhysicalAxis(b);ui[d]=vi.Value(s,offset[b]);uo[d]=vo.Value(s,offset[b]);gi[d]=vi.Gradient(s,offset[b]);go[d]=vo.Gradient(s,offset[b]);}
        const double k=vi.Value(s,offset[K()]),ko=vo.Value(s,offset[K()]),e=Point(Epsilon(),q,1).Value(s,offset[Epsilon()]);
        const double pii=Point(InnerPressure(),q,-1).Value(s,offset[InnerPressure()]),pio=Point(OuterPressure(),q,1).Value(s,offset[OuterPressure()]);
        const double nt=profile->Original().EddyViscosity(H),D=m_parameters.geometry.molecularNu+nt;
        const auto li=InnerStress(k,nt,gi);const auto ti=RANSWallInnerMeanStress(li,m_parameters.geometry.molecularNu,pii);
        const auto lo=OuterClosure(ko,e,go,q);RANSVector traction{},innerTraction{},outerTraction{};
        double targetJump=SolverUtils::WJInterfacePressureJump(k,D,li.strain,lo.reynoldsStress,n);
        for(int a=0;a<3;++a)for(int b=0;b<3;++b){innerTraction[a]+=ti[a][b]*n[b];outerTraction[a]+=(lo.reynoldsStress[a][b]+(a==b?pio-2*ko/3:0))*n[b];}
        for(int a=0;a<3;++a)traction[a]=innerTraction[a]+outerTraction[a];
        const auto gki=vi.Gradient(s,offset[K()]),gko=vo.Gradient(s,offset[K()]);
        const double ji=-(m_parameters.geometry.molecularNu+nt/SigmaK())*Dot(gki,n),jo=lo.kDiffusion*Dot(gko,n);
        interfaces.push_back(J::Object({{"physical",J::Numbers(Physical(q))},{"chart",J::Numbers(q)},
            {"geometry_inner_outer",J::Numbers(std::array<int,2>{vi.geometry,vo.geometry})},
            {"outward_inner_normal",J::Numbers(n)},{"inner_outward_traction",J::Numbers(innerTraction)},
            {"outer_outward_traction",J::Numbers(outerTraction)},{"traction_defect",J::Numbers(traction)},
            {"mass_flux_defect",J::Number(Dot(ui,n)-Dot(uo,n))},{"k_flux_inner_outer",J::Numbers(std::array<double,2>{ji,jo})},
            {"shared_k_defect",J::Number(k-ko)},{"K",J::Number(profile->K())},{"cutoff",J::Number(profile->Reference().cutoff)},
            {"epsilon_trace_target",J::Numbers(std::array<double,2>{e,profile->Dissipation(H)})},
            {"mechanical_jump_actual_target",J::Numbers(std::array<double,2>{pii-2*k/3-pio+2*ko/3,targetJump})}}));
    }
    Require(hash==RANSWallR3::Hash(s) && work==InitializationWork(),"open diagnostic mutated state or Lift count");
    return J::Object({{"schema",J::Quote("WJOpenTunnelV1")},{"state_sha256",J::Quote(hash)},
        {"gauge_unknowns",J::Number(size[Gauge()])},{"mass_flux",J::Number(totalMass)},
        {"boundary_flux_correction",J::Number(m_boundaryCorrection)},{"boundary_rows",J::List(rows)},
        {"stations",J::List(stations)},{"interface",J::List(interfaces)},{"boundary_load",J::Numbers(ExteriorResidual(s,time,data))},
        {"alpha_mass",J::Number(alpha)},{"time",J::Number(time)},{"history_sha256",J::Quote(RANSWallR3::Hash(history))},
        {"scope",J::Quote("physical fluxes and region-owned fields; prescribed natural and actual pointwise tractions distinct; no governing solve/Lift")}});
}
} // namespace Nektar
