// SPDX-License-Identifier: MIT
// Explicit developing boundary policy: affine or the qualified composite
// chart only. Unrelated legacy mapped-production paths remain guarded.
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <cmath>
#include <cstring>
#include <numeric>

namespace Nektar
{
using namespace SolverUtils;
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;
bool RANSRetainAttemptPacket(const std::string&,bool,bool,long long,int);

std::string RANSWallEvolutionPatch::BoundaryPolicyJSON() const
{
    return J::Object({{"mode",J::Quote(m_boundaryPolicy)},
        {"inlet",J::Quote("native finite/AWF seed Dirichlet; dependent epsilon corner takes exact E_C(Hw)")},
        {"top",J::Quote("impermeable full-selected-stress free slip; natural scalar flux; floor-only Yap distance")},
        {"outlet",J::Quote("upwind conservative flux: un+*q+un-*q_seed; zero outgoing scalar diffusion; sigma*n=-Pout*n+un-*(u-u_seed)")},
        {"pressure",J::Quote("absolute total normal traction datum, not a Dirichlet mechanical or modified pressure; free regional pressure; unused gauge fixed zero")},
        {"Pout",J::Number(m_outletTotalPressure)},
        {"backflow",J::Quote("immutable native fixed-reference velocity/k/outer-epsilon; positive scalar data required; no clipping")}});
}

RANSWallEvolutionPatch::ExternalFluxTerms RANSWallEvolutionPatch::DevelopingFlux(
    const Vector &s,int selectedRegion) const
{
    Require(OpenBoundary() && (!Mapped() || Composite()) && m_dimension==2 && s.size()==Size() && m_backflowSeed.size()==Size(),
        "developing flux needs affine state and immutable native backflow seed");
    RANSWallTiming::Scope timing("boundary.assembly");
    ExternalFluxTerms f{Vector(Size()),Vector(Size()),Vector(Size())};
    for(const auto &face:SurfaceView(0,0,upper[0],1,8))
    {
        const auto &v=face.test;const auto &el=fields[0]->Elements()[v.element];
        const int region=el.lower[2]>=m_parameters.geometry.height?1:0;
        if(selectedRegion>=0 && selectedRegion!=region)continue;
        const double un=v.Value(s,offset[0])*face.area[0]+v.Value(s,offset[1])*face.area[2];
        const double positive=std::max(0.,un),negative=std::min(0.,un);
        for(int b=0;b<InnerPressure();++b)
        {
            if(b==Epsilon() && !region)continue;
            const auto q=b==Epsilon()?GeometryPoint(b,v.geometry,v.coordinate):v;
            const double prescribed=q.Value(m_backflowSeed,offset[b]);
            Require(std::isfinite(prescribed) && (b<K() || prescribed>0),"invalid prescribed backflow turbulence data");
            const double flux=positive*q.Value(s,offset[b])+negative*prescribed;
            for(int i=0;i<q.ids.size();++i)
            {
                f.advection[offset[b]+q.ids[i]]+=face.weight*q.phi[i]*flux;
                if(b<K())f.traction[offset[b]+q.ids[i]]+=face.weight*q.phi[i]*m_outletTotalPressure*face.area[PhysicalAxis(b)];
            }
        }
    }
    return f;
}

void RANSWallEvolutionPatch::AddDevelopingMatrix(int block,const Vector &s,
    const std::vector<int> &rows,WallMatrix &a) const
{
    std::vector<int> index(Size(),-1);for(int i=0;i<rows.size();++i)index[rows[i]]=i;
    for(const auto &face:SurfaceView(0,0,upper[0],1,8))
    {
        const auto &v=face.test;const bool outer=fields[0]->Elements()[v.element].lower[2]>=m_parameters.geometry.height;
        if(block==Epsilon() && !outer)continue;
        const double un=std::max(0.,v.Value(s,offset[0])*face.area[0]+v.Value(s,offset[1])*face.area[2]);
        for(int b=block==0?0:block;b<(block==0?K():block+1);++b)
        {
            const auto q=b==Epsilon()?GeometryPoint(b,v.geometry,v.coordinate):v;
            for(int i=0;i<q.ids.size();++i)if(index[offset[b]+q.ids[i]]>=0)
                for(int j=0;j<q.ids.size();++j)if(index[offset[b]+q.ids[j]]>=0)
                    a[index[offset[b]+q.ids[i]]][index[offset[b]+q.ids[j]]]+=face.weight*q.phi[i]*un*q.phi[j];
        }
    }
}

std::string RANSWallEvolutionPatch::DevelopingChecks(const Vector &s,double time,
    const RANSWallEvolutionData &data) const
{
    std::vector<std::string> rows;double net=0,minimum=1e300,strong[2]={};
    for(int axis:{0,2})for(int side:{-1,1})
    {
        double flux=0,area=0;
        for(const auto &f:SurfaceView(0,axis,side<0?lower[axis]:upper[axis],side,8))
        {
            const double un=f.test.Value(s,offset[0])*f.area[0]+f.test.Value(s,offset[1])*f.area[2];
            flux+=f.weight*un;area+=f.weight*std::hypot(f.area[0],f.area[2]);
            if(axis==0 && side==1)minimum=std::min(minimum,un);
        }
        net+=flux;rows.push_back(J::Numbers(Vector{double(axis),double(side),area,flux}));
    }
    for(const auto &q:Quadrature(0,8))
    {const double d=q.Gradient(s,offset[0])[0]+q.Gradient(s,offset[1])[2];strong[q.coordinate[2]<m_parameters.geometry.height?0:1]+=q.weight*d*d;}
    const auto r=Residual(s,s,0,time,data,false,true);Vector weak;
    for(int b=InnerPressure();b<=Gauge();++b)weak.push_back(Maximum(Vector(r.begin()+offset[b],r.begin()+offset[b]+size[b])));
    return J::Object({{"definition",BoundaryPolicyJSON()},{"face_columns",J::Quote("axis,side,physical_area,outward_volume_flux")},
        {"faces",J::List(rows)},{"net_volume_flux",J::Number(net)},{"minimum_outlet_un",J::Number(minimum)},
        {"strong_divergence_inner_outer_L2",J::Numbers(Vector{std::sqrt(strong[0]),std::sqrt(strong[1])})},
        {"weak_regional_continuity_unused_gauge",J::Numbers(weak)},{"top_net_flux_adjustment",J::Number(m_boundaryCorrection)}});
}

std::string RANSWallEvolutionPatch::DevelopingBoundaryFixture(const Vector &s,
    const RANSWallEvolutionData &data) const
{
    Require(OpenBoundary(),"boundary fixture requires the explicit developing mode");
    std::vector<RANSWallEvolutionFunction> functions;
    for(int b=0;b<Gauge();++b)functions.push_back([b,this](const RANSWallCoordinate&,double){return b>=InnerPressure()?.3:(b>=K()?.1:0.);});
    const auto projected=Project(functions,0.);double absoluteProjection=0.;
    for(int b=InnerPressure();b<=OuterPressure();++b)for(const auto &q:Quadrature(b,4))absoluteProjection=std::max(absoluteProjection,std::abs(q.Value(projected,offset[b])-.3));
    Require(absoluteProjection<1e-12,"open projection silently removed absolute pressure level");
    // Independent divergence theorem with a frozen FULL selected tensor,
    // including its nonlinear/reference/pressure pieces. This is a flux-sign
    // fixture, not a claim that the manufactured constant tensor satisfies
    // the slip/outlet boundary data. Scalar diffusion is checked separately.
    RANSTensor manufacturedG{};manufacturedG[0][2]=.1;
    const RANSWallCoordinate anchor{(lower[0]+upper[0])/2,0,.04};
    const auto law=OuterClosure(.1,.1,manufacturedG,anchor);
    RANSTensor stress{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)stress[i][j]=-law.reynoldsStress[i][j]+(i==j?2*.1/3-.3:0.);
    Vector stressErrors;
    for(int b:{0,1,K(),Epsilon()})
    {
        Vector test(size[b]);RANSVector flux{};
        if(b<K())flux=stress[PhysicalAxis(b)];
        else{const double d=b==K()?law.kDiffusion:law.epsilonDiffusion;flux={.1*d,0,.2*d};}
        for(const auto &q:Quadrature(b,8))for(int j=0;j<q.ids.size();++j)
            test[q.ids[j]]+=q.weight*(q.derivative[0][j]*flux[0]+q.derivative[2][j]*flux[2]);
        for(int axis:{0,2})for(int side:{-1,1})
        {
            const double lo=b==Epsilon()?m_parameters.geometry.height:0.;
            const double pos=axis==0?(side<0?lower[0]:upper[0]):(side<0?lo:upper[2]);
            for(const auto &f:SurfaceView(b,axis,pos,side,8))for(int j=0;j<f.test.ids.size();++j)
                test[f.test.ids[j]]-=f.weight*f.test.phi[j]*(flux[0]*f.area[0]+flux[2]*f.area[2]);
        }
        stressErrors.push_back(Maximum(test));
    }
    Require(Maximum(stressErrors)<2e-12,"selected tensor/scalar diffusive boundary identity failed");
    auto negative=s;const auto ground=Point(K(),{upper[0],0,0});const auto prescribed=ExternalConstraints(0,data);
    int chosen=-1;for(int j=0;j<ground.ids.size();++j)if(!prescribed.count(offset[K()]+ground.ids[j]) && std::abs(ground.phi[j])>.5){chosen=j;break;}
    Require(chosen>=0,"outlet ground k must remain free");
    negative[offset[K()]+ground.ids[chosen]]-=(ground.Value(s,offset[K()])+1e-4)/ground.phi[chosen];
    const auto repair=ApplyPimpleKBound(negative,s,0.,data);
    Require(repair.success && repair.activated && ScalarDomain(repair.state,0,data).admissible,"new free-outlet k inventory/projection fixture failed");
    for(int i=0;i<Size();++i)if(i<offset[K()] || i>=offset[Epsilon()] || prescribed.count(i))
        Require(std::memcmp(&s[i],&repair.state[i],sizeof(double))==0,"synthetic correction changed non-k or essential data");
    Require(RANSRetainAttemptPacket("SelectedV1",repair.activated,false,9,32) &&
        !RANSRetainAttemptPacket("SelectedV1",false,false,9,32),"active/unselected attempt retention policy failed");
    Vector identityError;std::vector<std::string> fluxChecks;
    for(int b=0;b<InnerPressure();++b)
    {
        Vector identity(size[b]);
        // Independent conservative integration by parts for every native test
        // row, including constrained modes. Polynomial physical fields make
        // GL8 exact on this affine fixture.
        auto scalar=[&](const RANSWallCoordinate&p){return b==0?1+.2*p[0]+.1*p[2]:(b==1?.1*p[0]:.5+.1*p[0]+.2*p[2]);};
        const RANSVector dq=b==0?RANSVector{.2,0,.1}:(b==1?RANSVector{.1,0,0}:RANSVector{.1,0,.2});
        for(const auto &q:Quadrature(b,8))
        {
            const auto p=q.coordinate;const double u=1+.2*p[0]+.1*p[2],w=.1*p[0],v=scalar(p),div=.2*v+u*dq[0]+w*dq[2];
            for(int j=0;j<q.ids.size();++j)identity[q.ids[j]]+=q.weight*(q.derivative[0][j]*u*v+q.derivative[2][j]*w*v+q.phi[j]*div);
        }
        for(int axis:{0,2})for(int side:{-1,1})
        {
            const double z0=b==Epsilon()?m_parameters.geometry.height:0.;
            const double pos=axis==0?(side<0?lower[0]:upper[0]):(side<0?z0:upper[2]);
            for(const auto &f:SurfaceView(b,axis,pos,side,8))
            {const auto p=f.test.coordinate;const double un=(1+.2*p[0]+.1*p[2])*f.area[0]+.1*p[0]*f.area[2];
             for(int j=0;j<f.test.ids.size();++j)identity[f.test.ids[j]]-=f.weight*f.test.phi[j]*un*scalar(p);}
        }
        identityError.push_back(Maximum(identity));
    }
    for(double sign:{1.,-1.})
    {
        auto q=s;for(int b=0;b<InnerPressure();++b)for(int j=0;j<size[b];++j)
            q[offset[b]+j]=(b==0?sign:(b==1?0.:(b==K()?2.:.3)))*m_constants[b][j];
        const auto f=DevelopingFlux(q);Vector error;
        for(int b=0;b<InnerPressure();++b)
        {
            double actual=0,independent=0;for(int j=0;j<size[b];++j)actual+=m_constants[b][j]*f.advection[offset[b]+j];
            for(const auto &face:SurfaceView(b,0,upper[0],1,16))
                independent+=face.weight*(sign>0?sign*(b==0?sign:(b==1?0.:(b==K()?2.:.3))):sign*face.test.Value(m_backflowSeed,offset[b]));
            error.push_back(std::abs(actual-independent));
        }
        Require(Maximum(error)<2e-12,"outflow/backflow upwind flux identity failed");
        fluxChecks.push_back(J::Object({{"un",J::Number(sign)},{"independent_integral_errors",J::Numbers(error)}}));
    }
    const auto external=ExternalConstraints(0,data);Vector direction(Size());
    for(int b=InnerPressure();b<=OuterPressure();++b)for(int j=0;j<size[b];++j)direction[offset[b]+j]=m_constants[b][j];
    const auto pressure=AffinePressureAction(direction,external);
    const auto flow=AssemblePicardBlock(0,s,s,1.,0,data,external);Vector error(flow.rows.size());
    for(int i=0;i<flow.rows.size();++i)
    {double v=0;for(int j=0;j<flow.rows.size();++j)v+=flow.matrix[i][j]*direction[flow.rows[j]];error[i]=v-pressure[flow.rows[i]];}
    Require(Maximum(pressure)>1e-6 && Maximum(error)<1e-10 && external.at(offset[Gauge()])==0.,"open absolute pressure/full native G/gauge identity failed");
    double corner=0;const auto trace=Constraints(s,0,data);
    for(int id:m_traceIds)if(external.count(offset[Epsilon()]+id))corner=std::max(corner,std::abs(s[offset[Epsilon()]+id]-trace.at(offset[Epsilon()]+id)));
    Require(corner<1e-12 && Maximum(identityError)<2e-12,"boundary corner/conservative identity failed");
    return J::Object({{"status",J::Quote("passed")},{"projection_only_fixture_bound_calls","1"},{"projection_only_fixture_governing_solves","0"},{"projection_only_fixture",repair.report},{"active_packet_policy_passed","true"},{"all_native_mode_flux_identity_errors",J::Numbers(identityError)},
        {"outflow_backflow",J::List(fluxChecks)},{"selected_stress_and_scalar_diffusion_flux_errors",J::Numbers(stressErrors)},
        {"common_constant_pressure_action",J::Number(Maximum(pressure))},{"absolute_pressure_projection_error",J::Number(absoluteProjection)},{"fixture_setup_projections","6"},
        {"complete_G_action_error",J::Number(Maximum(error))},{"epsilon_inlet_interface_corner_error",J::Number(corner)},
        {"gauge",J::Quote("unused slot fixed zero; pressure level fixed by total traction")},{"governing_solves","0"}});
}

std::string RANSWallEvolutionPatch::FiniteTurbulenceBudget(const Vector &s,const Vector &history,
    double alpha,double time,const RANSWallEvolutionData &data) const
{
    Require(PrescribedInflow() && AcceptedAWF() && (!Mapped() || Composite()),"finite native budget is the affine accepted-AWF diagnostic");
    EvaluationScope evaluation(*this,data,time);Refresh(s,time,&data);
    std::array<Vector,2> rates{Vector(8),Vector(8)};
    // columns: k production, k sink, epsilon production, destruction, Yap,
    // counterfactual p_A0, prescribed k source, prescribed epsilon source.
    for(const auto &ray:m_rays)for(const auto &q:ray.points)
    {
        RANSTensor g{};for(int a=0;a<2;++a)g[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);
        const auto law=InnerStress(q.v.Value(s,offset[K()]),ray.profile->Original().EddyViscosity(Distance(q.v.coordinate)),g);
        rates[0][0]+=q.weight*law.production;rates[0][1]+=q.weight*ray.profile->Dissipation(Distance(q.v.coordinate));
        rates[0][5]+=q.weight*ray.profile->ReferenceProduction(Distance(q.v.coordinate));
        rates[0][6]+=q.weight*data.innerKSource(q.v.coordinate,time);
    }
    for(const auto &q:m_outer)
    {
        RANSTensor g{};for(int a=0;a<2;++a)g[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);
        const double e=q.e.Value(s,offset[Epsilon()]);const auto c=OuterClosure(q.v.Value(s,offset[K()]),e,g,q.v.coordinate);
        rates[1][0]+=q.weight*c.production;rates[1][1]+=q.weight*e;
        rates[1][2]+=q.weight*c.epsilonProduction;rates[1][3]+=q.weight*c.epsilonDestruction;rates[1][4]+=q.weight*c.yap;
        rates[1][6]+=q.weight*data.outerKSource(q.v.coordinate,time);rates[1][7]+=q.weight*data.outerEpsilonSource(q.v.coordinate,time);
    }
    std::vector<std::string> independentRates;
    for(int order:Composite()?std::vector<int>{8,16,32}:std::vector<int>{16,32})
    {
        std::array<Vector,2> values{Vector(8),Vector(8)};
        const int horizontal=Composite()?order:8;
        const auto gx=LibUtilities::PointsManager()[LibUtilities::PointsKey(horizontal,LibUtilities::eGaussGaussLegendre)];
        const auto gz=LibUtilities::PointsManager()[LibUtilities::PointsKey(order,LibUtilities::eGaussGaussLegendre)];
        for(const auto &el:fields[0]->Elements())for(int ix=0;ix<horizontal;++ix)
        {
            const double x=el.lower[0]+(gx->GetZ()[ix]+1)*(el.upper[0]-el.lower[0])/2;
            const bool inner=el.upper[2]<=m_parameters.geometry.height;const int region=inner?0:1;
            const auto profile=StateProfile(s,{x,0,0},time,data);
            std::vector<double> panels{el.lower[2],el.upper[2]};
            // Finite reference knots are incoming PHYSICAL heights. Each
            // regional chart is linear in its vertical coordinate at fixed s;
            // never treat an outer zHat knot as physical reference height.
            const double physicalLo=Physical({x,0,el.lower[2]})[2],physicalHi=Physical({x,0,el.upper[2]})[2];
            for(double knot:InflowBreakpoints())
            {
                const double chart=Composite()?el.lower[2]+(knot-physicalLo)*(el.upper[2]-el.lower[2])/(physicalHi-physicalLo):knot;
                if(chart>el.lower[2] && chart<el.upper[2])panels.push_back(chart);
            }
            for(double z:{profile->Reference().cutoff,profile->Original().ViscosityInterceptHeight(),m_parameters.geometry.roughness})
                if(z>el.lower[2] && z<el.upper[2])panels.push_back(z);
            std::sort(panels.begin(),panels.end());panels.erase(std::unique(panels.begin(),panels.end()),panels.end());
            for(int panel=1;panel<panels.size();++panel)for(int iz=0;iz<order;++iz)
            {
                const double z=panels[panel-1]+(gz->GetZ()[iz]+1)*(panels[panel]-panels[panel-1])/2;
                double w=gx->GetW()[ix]*gz->GetW()[iz]*(el.upper[0]-el.lower[0])*(panels[panel]-panels[panel-1])/4;
                const auto q=GeometryPoint(0,el.geometry,{x,0,z});w*=GeometryAt(q.coordinate).determinant;const double k=q.Value(s,offset[K()]);
                Require(k>=0 && (inner || k>0),"additional source-budget quadrature finds inadmissible k at chart="+J::Numbers(q.coordinate));
                RANSTensor g{};for(int a=0;a<2;++a)g[PhysicalAxis(a)]=q.Gradient(s,offset[a]);
                if(inner)
                {const auto c=InnerStress(k,profile->Original().EddyViscosity(z),g);values[0][0]+=w*c.production;
                 values[0][1]+=w*profile->Dissipation(z);values[0][5]+=w*profile->ReferenceProduction(z);values[0][6]+=w*data.innerKSource(q.coordinate,time);}
                else
                {const double e=GeometryPoint(Epsilon(),el.geometry,q.coordinate).Value(s,offset[Epsilon()]);Require(e>0,"additional source-budget epsilon invalid at chart="+J::Numbers(q.coordinate));
                 const auto c=OuterClosure(k,e,g,q.coordinate);values[1][0]+=w*c.production;values[1][1]+=w*e;
                 values[1][2]+=w*c.epsilonProduction;values[1][3]+=w*c.epsilonDestruction;values[1][4]+=w*c.yap;
                 values[1][6]+=w*data.outerKSource(q.coordinate,time);values[1][7]+=w*data.outerEpsilonSource(q.coordinate,time);}
                (void)region;
            }
        }
        independentRates.push_back(J::Object({{"Gauss_order",J::Number(order)},{"tangential_Gauss_order",J::Number(horizontal)},
            {"inner",J::Numbers(values[0])},{"outer",J::Numbers(values[1])}}));
    }
    const auto constraints=Constraints(s,time,data);
    std::vector<std::string> regions,faces,components;
    for(int region=0;region<2;++region)
    {
        const auto native=Residual(s,history,alpha,time,data,false,false,region);
        for(int b:{K(),Epsilon()})
        {
            if(b==Epsilon() && !region)continue;
            const auto mass=b==K()?PhysicalMass(b,region):m_mass[b];Vector delta(size[b]);
            for(int j=0;j<size[b];++j)delta[j]=s[offset[b]+j]-history[offset[b]+j];
            const auto storage=Multiply(mass,delta);double stored=0,free=0,reaction=0,interface=0,total=0;
            const auto atH=fields[b]->Boundary([](const RANSWallCoordinate&){return 0.;},[&](int axis,int side,const RANSWallGlobalField::Element &e){return axis==2 && std::abs((side?e.upper[2]:e.lower[2])-m_parameters.geometry.height)<1e-13;});
            for(int j=0;j<size[b];++j)
            {stored+=alpha*m_constants[b][j]*storage[j];const double v=m_constants[b][j]*native[offset[b]+j];total+=v;
             if(constraints.count(offset[b]+j))reaction+=v;else if(atH.count(j))interface+=v;else free+=v;}
            const double source=b==K()?rates[region][0]-rates[region][1]+rates[region][6]:rates[region][2]-rates[region][3]+rates[region][4]+rates[region][7];
            double open=0;if(OpenBoundary()){const auto f=DevelopingFlux(s,region);for(int j=0;j<size[b];++j)open+=m_constants[b][j]*(f.advection[offset[b]+j]+f.diffusion[offset[b]+j]);}
            regions.push_back(J::Object({{"region",J::Number(region)},{"block",J::Number(b)},
                {"physical_storage_rate",J::Number(stored)},{"source_rate",J::Number(source)},
                {"explicit_open_outward_flux",J::Number(open)},{"unmasked_constant_test_residual",J::Number(total)},
                {"free_interior_row_projection",J::Number(free)},{"constrained_row_reaction_projection",J::Number(reaction)},
                {"free_interface_row_projection",J::Number(interface)},
                {"exact_weak_balance_error",J::Number(stored-source+open-total)},
                {"row_partition_error",J::Number(free+reaction+interface-total)}}));
            const double lo=region?m_parameters.geometry.height:0.,hi=region?upper[2]:m_parameters.geometry.height;
            for(int axis:{0,2})for(int side:{-1,1})
            {
                const double pos=axis==0?(side<0?lower[0]:upper[0]):(side<0?lo:hi);double adv=0,diff=0;
                for(const auto &f:SurfaceView(0,axis,pos,side,16))
                {
                    const auto &el=fields[0]->Elements()[f.test.element];if((el.lower[2]>=m_parameters.geometry.height)!=bool(region))continue;
                    const auto q=b==Epsilon()?GeometryPoint(b,f.test.geometry,f.test.coordinate):f.test;
                    const double k=f.test.Value(s,offset[K()]);double d=0;
                    if(!region){const auto p=StateProfile(s,q.coordinate,time,data);d=m_parameters.geometry.molecularNu+p->Original().EddyViscosity(Distance(q.coordinate))/m_parameters.model.Parameters().sigmaK;}
                    else{RANSTensor g{};for(int a=0;a<2;++a)g[PhysicalAxis(a)]=f.test.Gradient(s,offset[a]);const auto e=GeometryPoint(Epsilon(),q.geometry,q.coordinate);const auto c=OuterClosure(k,e.Value(s,offset[Epsilon()]),g,q.coordinate);d=b==K()?c.kDiffusion:c.epsilonDiffusion;}
                    const auto grad=q.Gradient(s,offset[b]);
                    adv+=f.weight*q.Value(s,offset[b])*(f.test.Value(s,offset[0])*f.area[0]+f.test.Value(s,offset[1])*f.area[2]);
                    diff+=f.weight*d*(grad[0]*f.area[0]+grad[2]*f.area[2]);
                }
                faces.push_back(J::Numbers(Vector{double(region),double(b),double(axis),double(side),adv,diff}));
            }
        }
    }
    const auto r=ResidualResolutionComponents(s,history,alpha,time,data);
    const auto native=Residual(s,history,alpha,time,data,false,false);double sumError=0;
    for(int j=0;j<Size();++j)sumError=std::max(sumError,std::abs(r.terms.at("total_unmasked")[j]-native[j]));
    for(int b:{K(),Epsilon()})
    {
        std::vector<int> free;for(int i=0;i<size[b];++i)if(!constraints.count(offset[b]+i))free.push_back(i);
        WallMatrix mass(free.size(),Vector(free.size()));Vector total;
        for(int i=0;i<free.size();++i){total.push_back(native[offset[b]+free[i]]);for(int j=0;j<free.size();++j)mass[i][j]=m_mass[b][free[i]][free[j]];}
        const Factor factor(mass);const auto dual=factor.Apply(total);double norm=std::sqrt(std::max(0.,std::inner_product(total.begin(),total.end(),dual.begin(),0.)));
        auto names=RANSWallResolutionSpatialTerms();names.push_back("time_mass");
        for(const auto &name:names)
        {Vector load;double integral=0;for(int i=0;i<size[b];++i)integral+=m_constants[b][i]*r.terms.at(name)[offset[b]+i];for(int i:free)load.push_back(r.terms.at(name)[offset[b]+i]);
         const auto inverse=factor.Apply(load);const double own=std::sqrt(std::max(0.,std::inner_product(load.begin(),load.end(),inverse.begin(),0.)));
         components.push_back(J::Object({{"block",J::Number(b)},{"term",J::Quote(name)},{"constant_test_projection",J::Number(integral)},
             {"free_mass_dual_norm",J::Number(own)},{"signed_projection_on_total_residual",norm>1e-30?J::Number(std::inner_product(load.begin(),load.end(),dual.begin(),0.)/norm):"null"}}));}
    }
    // Independent analytical parallel seed, using physical GL64 panels.
    const double H=m_parameters.geometry.height,K=InflowAtHeight(H).k;
    const auto p=Profile(K,{lower[0],0,0},1.);const double UH=m_robustUstar*m_robustUstar*p->Original().Integrals().j0;
    std::vector<double> panels{0.,H,p->Reference().cutoff};for(double v:{p->Original().ViscosityInterceptHeight(),m_parameters.geometry.roughness})if(v>0 && v<H)panels.push_back(v);
    std::sort(panels.begin(),panels.end());panels.erase(std::unique(panels.begin(),panels.end()),panels.end());
    const auto gl=LibUtilities::PointsManager()[LibUtilities::PointsKey(64,LibUtilities::eGaussGaussLegendre)];double P=0,E=0;
    for(int j=1;j<panels.size();++j)for(int i=0;i<64;++i)
    {const double z=panels[j-1]+(gl->GetZ()[i]+1)*(panels[j]-panels[j-1])/2,w=gl->GetW()[i]*(panels[j]-panels[j-1])/2;
     RANSTensor g{};g[0][2]=p->Original().VelocityGradient(z,{UH,0},{0,0})[0];P+=w*InnerStress(p->ReferenceK(z),p->Original().EddyViscosity(z),g).production;E+=w*p->Dissipation(z);}
    const double ideal=std::pow(m_parameters.model.Q()*m_robustUstar*m_robustUstar/K,2);
    Require(std::abs(P/E-ideal)<1e-11 && sumError<1e-9,"native turbulence budget decomposition or independent ideal ratio failed");
    return J::Object({{"rate_columns",J::Quote("Pk,Ek,P_epsilon,D_epsilon,Yap,p_A0,prescribed_k,prescribed_epsilon")},
        {"inner_rates",J::Numbers(rates[0])},{"outer_rates",J::Numbers(rates[1])},{"regional_weak_balances",J::List(regions)},
        {"independent_source_integration",J::List(independentRates)},
        {"physical_face_columns",J::Quote("region,block,axis,side,outward_advection,outward_diffusion")},{"physical_faces",J::List(faces)},
        {"weak_components",J::List(components)},{"native_component_sum_max_error",J::Number(sumError)},
        {"ideal_seed_ratio",J::Number(ideal)},{"independent_analytical_P_over_E",J::Number(P/E)},
        {"native_inner_P_over_E",J::Number(rates[0][0]/rates[0][1])},{"history_alpha",J::Number(alpha)},
        {"ground_reaction",FiniteGroundReaction(s,time,data)},
        {"physical_time",J::Number(time)},{"governing_solves","0"},{"physical_commits","0"},
        {"scope",J::Quote("Native physical quadrature and saved callback history. Exact weak row partition is separate from gradient-based face flux estimates; interface CG row projections are not unique conservative control-volume fluxes. Source and diffusion are separate; net destruction counted once. No Lift or epsilon substitution.")}});
}

std::string RANSWallEvolutionPatch::FiniteGroundReaction(const Vector &s,double time,
    const RANSWallEvolutionData &data) const
{
    // The same physical comparison lifting is used on both flat meshes. It
    // vanishes at side corners and above the common .00025 native interface.
    // This is a weak-load diagnostic, not a new wall-traction condition.
    // The neutral family has its own shared fitted physical panel. Preserve
    // the old finite-family lifting exactly; neither is a first-cell budget.
    const double height=AnalyticalNeutralInflow()?.0002:.00025;
    const auto weight=[&](const RANSWallCoordinate &p)
    {const double x=(p[0]-lower[0])/(upper[0]-lower[0]);return 4*x*(1-x);};
    const auto function=[&](const RANSWallCoordinate &p)
    {return weight(p)*std::max(0.,1-p[2]/height);};
    const auto lift=fields[0]->Project(function);
    double error=0;for(const auto &q:Quadrature(0,8))error=std::max(error,std::abs(q.Value(lift)-function(q.coordinate)));
    Require(error<1e-11,"finite ground reaction requires the declared common lifting interface");
    const auto residual=ResidualResolutionComponents(s,s,0,time,data);
    long double total=0,free=0;std::vector<std::pair<std::string,std::string>> components;
    for(int i=0;i<size[0];++i)
    {const double v=lift[i]*residual.terms.at("spatial_sum")[i];total+=v;if(!residual.constraints.count(i))free+=v;}
    for(const auto &name:RANSWallResolutionSpatialTerms())
    {long double v=0;for(int i=0;i<size[0];++i)v+=(long double)lift[i]*residual.terms.at(name)[i];components.push_back({name,J::Number(v)});}
    double sampled=0,measure=0;
    for(const auto &face:SurfaceView(0,2,0.,-1,16))
    {
        const auto &q=face.test;const auto p=StateProfile(s,q.coordinate,time,data);RANSTensor G{};
        for(int b=0;b<m_dimension;++b)G[PhysicalAxis(b)]=q.Gradient(s,offset[b]);
        const auto stress=RANSWallInnerMeanStress(InnerStress(q.Value(s,offset[K()]),p->Original().EddyViscosity(0),G),m_parameters.geometry.molecularNu,0.);
        const double scale=std::hypot(face.area[0],face.area[2]),w=face.weight*scale*weight(q.coordinate);
        measure+=w;sampled+=w*(stress[0][0]*(-face.area[0]/scale)+stress[0][2]*(-face.area[2]/scale));
    }
    return J::Object({{"lifting",J::Quote("global-x Q2(s) times max(0,1-r/lifting_height), zero side corners")},
        {"lifting_height",J::Number(height)},
        {"weighted_surface_measure",J::Number(measure)},{"weighted_sampled_global_x_traction",J::Number(sampled/measure)},
        {"full_spatial_lifting_work",J::Number(total)},{"free_row_lifting_work",J::Number(free)},
        {"prescribed_row_reaction",J::Number(total-free)},{"inferred_ground_traction",J::Number(-(total-free)/measure)},
        {"signed_volume_components",J::Object(components)},{"lifting_representation_error",J::Number(error)},
        {"scope",J::Quote("spatial weak contraction; free/unsteady residual and other loads remain explicit; not an exact physical stress estimator")}});
}
}
