///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>

using namespace Nektar;
using namespace Nektar::SolverUtils;
using Vector=RANSWallAlgebra::Vector;
using WallMatrix=RANSWallAlgebra::DenseMatrix;
using Vec=RANSVector;
using Tensor=RANSTensor;
using RANSWallAlgebra::Require;
namespace
{
double Dot(const Vec &a,const Vec &b){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
double Norm(const Vec &a){return std::hypot(a[0],a[1],a[2]);}
double Max(const Vector &v){double r=0;for(double x:v){Require(std::isfinite(x),"nonfinite probe load");r=std::max(r,std::abs(x));}return r;}
bool Same(double a,double b){return std::abs(a-b)<1e-11;}
void Write(std::ostream &o,const Vector &v)
{o<<'[';bool first=true;for(double x:v){if(!first)o<<',';first=false;o<<x;}o<<']';}
void Write(std::ostream &o,const Vec &v){Write(o,Vector(v.begin(),v.end()));}
void Rows(std::ostream &o,const WallMatrix &v)
{o<<'[';bool first=true;for(const auto &row:v){if(!first)o<<',';first=false;Write(o,row);}o<<']';}
void TensorRows(std::ostream &o,const Tensor &v)
{o<<'[';for(int i=0;i<3;++i){if(i)o<<',';Write(o,v[i]);}o<<']';}
void Reduce(const RANSWallEvolutionPatch &p,Vector &v){p.Comm()->AllReduce(v,LibUtilities::ReduceSum);}
Vector MapParameters(const RANSWallNormalMap &m)
{return {static_cast<double>(m.kind),static_cast<double>(m.dimension),m.height,m.heightSlopeX,m.heightSlopeY,m.heightAmplitude,
         m.heightWaveX,m.heightWaveY,m.slopeX,m.slopeY,m.radius,m.amplitude,m.lengthX,m.lengthY};}
WallMatrix RowMatrix(const Vector &v,int n)
{WallMatrix a(n,Vector(3));for(int i=0;i<n;++i)for(int d=0;d<3;++d)a[i][d]=v[3*i+d];return a;}
void Describe(std::ostream &o,const RANSWallEvolutionPatch &p,int block)
{
    const auto field=p.fields[block];const auto native=field->Field()->GetLocalToGlobalMap();
    o<<"{\"block\":"<<block<<",\"size\":"<<field->Size()<<",\"constant\":";
    Write(o,field->Project([](const Vec&){return 1.;}));o<<",\"mesh\":[";
    for(int e=0;e<static_cast<int>(field->Elements().size());++e)
    {
        if(e){o<<',';}
        const auto &g=field->Elements()[e];const int count=field->Field()->GetExp(e)->GetNcoeffs(),offset=field->Field()->GetCoeff_Offset(e);
        o<<"{\"geometry\":"<<g.geometry<<",\"owner\":"<<p.Owner(g.geometry)<<",\"degree\":"<<g.degree<<",\"lower\":";Write(o,g.lower);
        o<<",\"upper\":";Write(o,g.upper);o<<",\"ids\":[";
        for(int j=0;j<count;++j){if(j)o<<',';o<<native->GetLocalToGlobalMap(offset+j);}o<<"],\"sign\":[";
        for(int j=0;j<count;++j){if(j)o<<',';o<<native->GetLocalToGlobalSign(offset+j);}o<<"]}";
    }
    o<<"]}";
}
void Row(std::ostream &o,const RANSWallEvolutionPatch &p,const RANSWallGlobalPoint &q)
{
    o<<"{\"chart\":";Write(o,q.coordinate);o<<",\"point\":";Write(o,p.Physical(q.coordinate));
    o<<",\"geometry\":"<<q.geometry<<",\"owner\":"<<p.Owner(q.geometry)<<",\"weight\":"<<q.weight<<",\"ids\":[";
    for(int j=0;j<static_cast<int>(q.ids.size());++j){if(j)o<<',';o<<q.ids[j];}o<<"],\"phi\":";Write(o,q.phi);
    o<<",\"constant\":";Write(o,q.constant);o<<",\"gradient\":[";
    for(int d=0;d<3;++d){if(d)o<<',';Write(o,q.derivative[d]);}
    const auto g=p.Mapping().At(q.coordinate);
    o<<"],\"volume_metric\":"<<g.determinant<<",\"normal_J\":"<<g.J<<",\"base_area\":"<<g.baseArea
     <<",\"physical_r\":"<<g.distance<<",\"H\":"<<g.height<<'}';
}
void BoundedQuadrature(std::ostream &o,const RANSWallEvolutionPatch &p)
{
    o<<'[';for(int block=0;block<p.Gauge();++block)
    {
        if(block){o<<',';}
        const auto points=p.Quadrature(block,8);std::map<int,int> index,count;
        for(const auto &q:points)++count[q.geometry];
        o<<"{\"block\":"<<block<<",\"order\":8,\"total_points\":"<<points.size()<<",\"rows\":[";int emitted=0;
        for(const auto &q:points)
        {
            const int i=index[q.geometry]++;if(i!=0 && i!=count[q.geometry]/2)continue;
            if(emitted++){o<<',';}
            Row(o,p,q);
        }
        o<<"]}";
    }o<<']';
}
void Masses(std::ostream &o,const RANSWallEvolutionPatch &p)
{
    o<<'[';bool first=true;
    for(int block:{0,p.Epsilon()})
    {
        const int n=p.size[block];Vector mass(n*n);
        for(const auto &q:p.Quadrature(block))if(p.Owner(q.geometry)==p.Comm()->GetRank())
            for(int i=0;i<static_cast<int>(q.ids.size());++i)for(int j=0;j<static_cast<int>(q.ids.size());++j)
                mass[q.ids[i]*n+q.ids[j]]+=q.weight*q.phi[i]*q.phi[j];
        Reduce(p,mass);if(!first){o<<',';}first=false;o<<"{\"block\":"<<block<<",\"shared_blocks\":[";
        if(block==0){for(int a=0;a<p.Dimension();++a){if(a)o<<',';o<<a;}o<<','<<p.K();}else o<<block;
        o<<"],\"values\":[";for(int i=0;i<n;++i){if(i)o<<',';o<<'[';for(int j=0;j<n;++j){if(j)o<<',';o<<mass[i*n+j];}o<<']';}o<<"]}";
    }o<<']';
}
// An actual all-mode weak Piola check. Volume derivatives and every physical
// cell face are integrated independently; prescribed source is identically 0.
void WeakMetric(std::ostream &o,const RANSWallEvolutionPatch &p)
{
    o<<'[';int level=0;
    for(int order:{8,12,16})
    {
        if(level++){o<<',';}
        o<<"{\"order\":"<<order<<",\"prescribed_volume_source\":[0,0,0],\"fields\":[";
        for(int block=0;block<p.Gauge();++block)
        {
            if(block){o<<',';}
            const int n=p.size[block];Vector volume(3*n),external(3*n),seam(3*n),flux(3),owned(p.Comm()->GetSize());double measure=0;
            Vec lower{1e100,1e100,1e100},upper{-1e100,-1e100,-1e100};std::array<std::set<double>,3> positions;
            for(const auto &e:p.fields[block]->Elements())for(int a=0;a<3;++a)
            {lower[a]=std::min(lower[a],e.lower[a]);upper[a]=std::max(upper[a],e.upper[a]);positions[a].insert(e.lower[a]);positions[a].insert(e.upper[a]);}
            for(const auto &q:p.Quadrature(block,order))if(p.Owner(q.geometry)==p.Comm()->GetRank())
            {
                measure+=q.weight;owned[p.Comm()->GetRank()]+=1;
                for(int j=0;j<static_cast<int>(q.ids.size());++j)for(int a=0;a<3;++a)volume[3*q.ids[j]+a]+=q.weight*q.derivative[a][j];
            }
            for(int axis=0;axis<3;++axis)if(p.Dimension()==3 || axis!=1)
                for(double position:positions[axis])for(int side:{-1,1})
                    for(const auto &face:p.Surface(block,axis,position,side,order))if(p.Owner(face.test.geometry)==p.Comm()->GetRank())
                    {
                        const bool boundary=Same(position,side<0?lower[axis]:upper[axis]);auto &load=boundary?external:seam;
                        for(int a=0;a<3;++a)
                        {
                            if(boundary)flux[a]+=face.weight*face.area[a];
                            for(int j=0;j<static_cast<int>(face.test.ids.size());++j)
                                load[3*face.test.ids[j]+a]+=face.weight*face.test.phi[j]*face.area[a];
                        }
                    }
            Reduce(p,volume);Reduce(p,external);Reduce(p,seam);Reduce(p,flux);Reduce(p,owned);p.Comm()->AllReduce(measure,LibUtilities::ReduceSum);
            Vector weak(3*n),allCells(3*n);for(int j=0;j<3*n;++j){weak[j]=volume[j]-external[j];allCells[j]=weak[j]-seam[j];}
            o<<"{\"block\":"<<block<<",\"size\":"<<n<<",\"weak_residual\":";Rows(o,RowMatrix(weak,n));o<<",\"seam_load\":";Rows(o,RowMatrix(seam,n));
            o<<",\"all_cell_residual\":";Rows(o,RowMatrix(allCells,n));o<<",\"maximum\":"<<Max(weak)<<",\"seam_maximum\":"<<Max(seam)
             <<",\"all_cell_maximum\":"<<Max(allCells)<<",\"constant_flux_boundary\":";Write(o,flux);
            o<<",\"owned_quadrature_counts\":";Write(o,owned);o<<",\"volume\":"<<measure<<'}';
        }o<<"]}";
    }o<<']';
}
void Velocities(const RANSWallEvolutionPatch &p,const Vector &s,const RANSWallGlobalPoint &q,Vec &u,Tensor &G)
{for(int a=0;a<p.Dimension();++a){const int d=p.PhysicalAxis(a);u[d]=q.Value(s,p.offset[a]);G[d]=q.Gradient(s,p.offset[a]);}}
void Interface(std::ostream &o,const RANSWallEvolutionPatch &p,const Vector &s)
{
    const int order=12;const auto &param=p.Parameters();Vector integral(20),sharedDefect(5*p.size[0]);
    double heightControl=0,velocityJump=0,tractionJump=0,epsilonTraceError=0;
    std::ostringstream samples;samples<<std::setprecision(17);int emitted=0;
    for(const auto &face:p.Surface(0,2,param.geometry.height,1,order))
    {
        const auto &q=face.test.coordinate;const auto metric=p.Mapping().At(q);const auto base=p.Mapping().At({q[0],q[1],0});
        const auto vi=p.Point(0,q,-1),vo=p.Point(0,q,1),pi=p.Point(p.InnerPressure(),q,-1),po=p.Point(p.OuterPressure(),q,1),ep=p.Point(p.Epsilon(),q,1);
        Vec u{},uo{};Tensor Gi{},Go{},Ti{},To{};Velocities(p,s,vi,u,Gi);Velocities(p,s,vo,uo,Go);
        const double k=vi.Value(s,p.offset[p.K()]),e=ep.Value(s,p.offset[p.Epsilon()]),K=p.SampleK(s,q),Ew=p.DependentEpsilon(s,q);
        const auto ki=vi.Gradient(s,p.offset[p.K()]),ko=vo.Gradient(s,p.offset[p.K()]),ge=ep.Gradient(s,p.offset[p.Epsilon()]);
        AnalyticalWallFunction awf({param.geometry.molecularNu,metric.height,param.geometry.roughness,param.geometry.sampleHeight},K);
        const double nt=awf.EddyViscosity(metric.height);
        const auto closure=p.OuterClosure(k,e,Go,q);
        Ti=RANSWallInnerMeanStress(p.InnerStress(k,nt,Gi),param.geometry.molecularNu,
            pi.Value(s,p.offset[p.InnerPressure()]));
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            To[a][b]=-closure.reynoldsStress[a][b]+(a==b?2*k/3-po.Value(s,p.offset[p.OuterPressure()]):0.);
        Vec tractionI{},tractionO{},N{},ut{},v{},gradH{},advMomentum{},fkI{},fkO{},feO{},feI{};
        const double areaNorm=Norm(face.area),mass=Dot(u,face.area),w=Dot(u,metric.normal);
        const double H1=p.Mapping().heightSlopeX+p.Mapping().heightAmplitude*p.Mapping().heightWaveX*std::cos(p.Mapping().heightWaveX*q[0])*std::cos(p.Mapping().heightWaveY*q[1]);
        const double H2=p.Mapping().heightSlopeY-p.Mapping().heightAmplitude*p.Mapping().heightWaveY*std::sin(p.Mapping().heightWaveX*q[0])*std::sin(p.Mapping().heightWaveY*q[1]);
        const auto t1=Vec{base.jacobian[0][0],base.jacobian[1][0],base.jacobian[2][0]},t2=Vec{base.jacobian[0][1],base.jacobian[1][1],base.jacobian[2][1]};
        const double g11=Dot(t1,t1),g12=Dot(t1,t2),g22=Dot(t2,t2),det=g11*g22-g12*g12;
        for(int a=0;a<3;++a)
        {
            N[a]=face.area[a]/areaNorm;ut[a]=u[a]-w*metric.normal[a];v[a]=ut[a];
            gradH[a]=((g22*H1-g12*H2)*t1[a]+(g11*H2-g12*H1)*t2[a])/det;
        }
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)v[a]-=metric.height*metric.Br[a][b]*ut[b];
        for(int a=0;a<3;++a)
        {
            tractionI[a]=Dot(Ti[a],face.area);tractionO[a]=Dot(To[a],face.area);advMomentum[a]=u[a]*mass;
            fkI[a]=u[a]*k-(param.geometry.molecularNu+nt/param.model.Parameters().sigmaK)*ki[a];
            fkO[a]=u[a]*k-closure.kDiffusion*ko[a];feO[a]=u[a]*e-closure.epsilonDiffusion*ge[a];feI[a]=u[a]*Ew;
        }
        const double m=mass/metric.baseArea,mFormula=metric.J*(w-Dot(v,gradH)),naive=metric.J*w;
        Require(std::abs(m-mFormula)<1e-10*std::max(1.,std::abs(m)),"matching mass-flux metric identity failed");
        Vec jump{};for(int a=0;a<3;++a){jump[a]=(tractionI[a]-tractionO[a])/areaNorm;velocityJump=std::max(velocityJump,std::abs(u[a]-uo[a]));}
        tractionJump=std::max(tractionJump,Norm(jump));heightControl=std::max(heightControl,std::abs(m-naive));epsilonTraceError=std::max(epsilonTraceError,std::abs(e-Ew));
        if(p.Owner(face.test.geometry)==p.Comm()->GetRank())
        {
            const double weight=face.weight;integral[0]+=weight*areaNorm;integral[1]+=weight*mass;integral[2]+=weight*k*mass;
            const std::array<double,5> shared{mass,advMomentum[0],advMomentum[1],advMomentum[2],k*mass};
            for(int a=0;a<5;++a)
            {
                for(int j=0;j<static_cast<int>(vi.ids.size());++j)sharedDefect[5*vi.ids[j]+a]+=weight*vi.phi[j]*shared[a];
                for(int j=0;j<static_cast<int>(vo.ids.size());++j)sharedDefect[5*vo.ids[j]+a]-=weight*vo.phi[j]*shared[a];
            }
            for(int a=0;a<3;++a){integral[3+a]+=weight*advMomentum[a];integral[6+a]+=weight*tractionI[a];integral[9+a]+=weight*tractionO[a];}
            integral[12]+=weight*Dot(tractionI,N);integral[13]+=weight*Dot(tractionO,N);
            integral[14]+=weight*Dot(fkI,face.area);integral[15]+=weight*Dot(fkO,face.area);
            integral[16]+=weight*Dot(feO,face.area);integral[17]+=weight*Dot(feI,face.area);
            integral[18]+=weight*(Dot(feO,face.area)-Dot(feI,face.area));integral[19]+=weight*metric.baseArea*naive;
        }
        if(emitted<8)
        {
            if(emitted++){samples<<',';}
            samples<<"{\"chart\":";Write(samples,q);samples<<",\"point\":";Write(samples,p.Physical(q));
            samples<<",\"area\":";Write(samples,face.area);samples<<",\"weight\":"<<face.weight<<",\"N_H\":";Write(samples,N);
            samples<<",\"normal\":";Write(samples,metric.normal);samples<<",\"U\":";Write(samples,u);samples<<",\"G_inner\":";TensorRows(samples,Gi);
            samples<<",\"G_outer\":";TensorRows(samples,Go);samples<<",\"traction_inner_area\":";Write(samples,tractionI);
            samples<<",\"traction_outer_area\":";Write(samples,tractionO);samples<<",\"shared_advective_momentum\":";Write(samples,advMomentum);
            samples<<",\"k\":"<<k<<",\"sampled_k\":"<<K<<",\"epsilon\":"<<e<<",\"dependent_epsilon\":"<<Ew
                   <<",\"physical_r\":"<<metric.distance<<",\"J\":"<<metric.J<<",\"base_area\":"<<metric.baseArea
                   <<",\"m_H\":"<<m<<",\"m_H_formula\":"<<mFormula<<",\"m_H_naive\":"<<naive
                   <<",\"shared_advective_mass\":"<<mass<<",\"shared_advective_k\":"<<k*mass
                   <<",\"k_flux_inner\":"<<Dot(fkI,face.area)<<",\"k_flux_outer\":"<<Dot(fkO,face.area)
                   <<",\"epsilon_flux_outer\":"<<Dot(feO,face.area)<<",\"epsilon_inner_advective_flux\":"<<Dot(feI,face.area)
                   <<",\"epsilon_flux_reaction\":"<<Dot(feO,face.area)-Dot(feI,face.area)<<'}';
        }
    }
    Reduce(p,integral);Reduce(p,sharedDefect);
    for(double *v:{&heightControl,&velocityJump,&tractionJump,&epsilonTraceError})p.Comm()->AllReduce(*v,LibUtilities::ReduceMax);
    o<<"{\"quadrature_order\":"<<order<<",\"integral_labels\":[\"area\",\"shared_advective_mass\",\"shared_advective_k\",\"adv_momentum_x\",\"adv_momentum_y\",\"adv_momentum_z\",\"traction_inner_x\",\"traction_inner_y\",\"traction_inner_z\",\"traction_outer_x\",\"traction_outer_y\",\"traction_outer_z\",\"normal_traction_inner\",\"normal_traction_outer\",\"k_flux_inner\",\"k_flux_outer\",\"epsilon_flux_outer\",\"epsilon_inner_advective_flux\",\"epsilon_flux_reaction\",\"naive_mass\"],\"integrals\":";
    Write(o,integral);o<<",\"maximum_velocity_trace_jump\":"<<velocityJump<<",\"maximum_traction_trace_jump\":"<<tractionJump
     <<",\"maximum_epsilon_projection_difference\":"<<epsilonTraceError<<",\"height_slope_omission_control\":"<<heightControl
     <<",\"shared_advective_seam_maximum\":"<<Max(sharedDefect)<<",\"shared_advective_seam_load\":[";
    for(int i=0;i<p.size[0];++i){if(i)o<<',';o<<'[';for(int a=0;a<5;++a){if(a)o<<',';o<<sharedDefect[5*i+a];}o<<']';}
    o<<"],\"samples\":["<<samples.str()<<"],\"interpretation\":\"Unsolved generic state: shared advective traces; regional traction/diffusive jumps are diagnostics, not zero constraints\"}";
}
void InnerEpsilonReaction(std::ostream &o,const RANSWallEvolutionPatch &p,const Vector &s,const Vector &history,double alpha)
{
    o<<'[';int count=0;
    for(const auto &e:p.innerVelocity->Elements())
    {
        Vec q{e.lower[0]+.37*(e.upper[0]-e.lower[0]),p.Dimension()==3?e.lower[1]+.41*(e.upper[1]-e.lower[1]):0.,e.lower[2]+.43*(e.upper[2]-e.lower[2])};
        const auto row=p.Point(0,q);Vec u{};Tensor G{};Velocities(p,s,row,u,G);Vec gradient{};const auto metric=p.Mapping().At(q);
        const double E=p.DependentEpsilon(s,q),old=p.DependentEpsilon(history,q);const double h=1e-5*std::min(1.,e.upper[2]-e.lower[2]);
        for(int a=0;a<3;++a)if(p.Dimension()==3 || a!=1)
        {
            auto plus=q,minus=q;plus[a]+=h;minus[a]-=h;const double derivative=(p.DependentEpsilon(s,plus)-p.DependentEpsilon(s,minus))/(2*h);
            for(int b=0;b<3;++b)gradient[b]+=derivative*metric.inverse[a][b];
        }
        const double divergence=G[0][0]+G[1][1]+G[2][2],reaction=alpha*(E-old)+E*divergence+Dot(u,gradient);
        if(count++){o<<',';}
        o<<"{\"chart\":";Write(o,q);o<<",\"epsilon\":"<<E<<",\"history_epsilon\":"<<old<<",\"gradient\":";Write(o,gradient);
        o<<",\"U\":";Write(o,u);o<<",\"divergence_U\":"<<divergence<<",\"diagnostic_source\":0,\"reaction\":"<<reaction<<",\"chart_fd_h\":"<<h<<'}';
    }o<<']';
}
}
namespace
{
void TransitionAudit(std::ostream &o,const RANSWallEvolutionPatch &p,const Vector &s,const Vector &history,
    double time,double alpha,const RANSWallEvolutionData &data,const LibUtilities::SessionReaderSharedPtr &serial,
    const SpatialDomains::MeshGraphSharedPtr &graph,const LibUtilities::CommSharedPtr &comm,RANSWallPatchParameters parameters)
{
    const auto &c=p.TransitionControls();const auto &t=c.state;const auto distances=p.TransitionDistances();
    o<<"{\"controls\":";Write(o,Vector{c.start,c.end,distances[0],distances[1],distances[2],
        t.on,t.off,t.marginOff,t.marginOn,t.reynoldsOff,t.reynoldsOn});
    o<<",\"weights\":";Write(o,Vector(t.weights.begin(),t.weights.end()));
    o<<",\"reference\":";Write(o,p.TransitionReference());
    const auto trace=p.DissipationTrace(s,time,data);o<<",\"trace\":[";
    for(int part=0;part<2;++part){if(part)o<<',';o<<'[';int count=0;
        for(auto [id,value]:trace[part]){if(count++)o<<',';Write(o,Vector{static_cast<double>(id),value});}o<<']';}o<<']';
    o<<",\"samples\":[";int count=0;
    for(const auto &e:p.outerVelocity->Elements())if(Same(e.lower[2],parameters.geometry.height))
        for(double fraction:{.137,.481,.863})
    {
        const Vec q{e.lower[0]+fraction*(e.upper[0]-e.lower[0]),p.Dimension()==3?e.lower[1]+.37*(e.upper[1]-e.lower[1]):0.,parameters.geometry.height};
        if(count++){o<<',';}
        RANSWallEvolutionData zero(p.Dimension());const auto input=p.TransitionInput(s,q,time,data);
        const auto v=p.Transition(s,q,time,data),rv=p.Transition(p.TransitionReference(),q,time,zero);
        const auto profile=p.StateProfile(s,q,time,data);const auto &pr=profile->Reference();
        Require(std::abs(p.DependentEpsilon(s,q,time,data)-p.WallDissipation(s,q,*profile))<1e-14,
            "state-aware dependent epsilon evaluation differs");
        o<<"{\"chart\":";Write(o,q);o<<",\"input\":";Write(o,Vector{input.coordinate,input.budgetMargin,input.reynolds});
        o<<",\"gates\":";Write(o,Vector{v.chi,v.spatial,v.state,v.admissibility,v.marginGate,v.reynoldsGate,v.d2});
        o<<",\"reference_d2\":"<<rv.d2<<",\"reference_state\":"<<rv.state<<",\"sensors\":[";
        const auto &a=input.sensors;Write(o,a.velocity);o<<',';Write(o,a.velocityDifference);o<<',';Write(o,Vector{a.normal,a.k,a.epsilon});o<<',';
        Write(o,a.pressure);o<<',';Write(o,a.convection);o<<"],\"profile\":";
        Write(o,Vector{profile->K(),pr.cutoff,static_cast<double>(pr.cutoffIterations),pr.weight,profile->Dissipation(0),profile->Dissipation(parameters.geometry.height)});
        double correction=0,budget=0,expected=0;
        for(const auto &z:profile->Quadrature())
        {const auto value=profile->AtProduction(z.r,0);const auto old=profile->Original().OriginalDissipation(z.r);
            correction+=z.weight*z.jacobian*value.redistribution;budget+=z.weight*z.jacobian*(value.kSource+old);
            expected+=z.weight*z.jacobian*pr.weight*(old-value.referenceProduction);}
        o<<",\"budgets\":";Write(o,Vector{correction,budget,expected,profile->RedistributionFlux(0),profile->RedistributionFlux(parameters.geometry.height)});
        const auto parts=p.ForcingParts(s,q,time,data);o<<",\"forcing_parts\":[";Write(o,parts[0]);o<<',';Write(o,parts[1]);o<<"]}";
    }o<<']';
    auto low=s;const auto constant=p.fields[p.K()]->Project([](const Vec&){return 1.;});
    for(int i=0;i<p.size[p.K()];++i)low[p.offset[p.K()]+i]=1e-10*constant[i];
    const Vec inactive{(p.lower[0]-3.)/2,p.Dimension()==3?.37:0.,parameters.geometry.sampleHeight};
    const auto fallback=p.StateProfile(low,inactive,time,data);
    o<<",\"inactive_low_K\":";Write(o,Vector{fallback->K(),fallback->Reference().weight,fallback->Reference().cutoff,
        static_cast<double>(fallback->Reference().cutoffIterations),fallback->Reference().margin,fallback->Reference().reynolds});
    // Equal-input recovery uses two instances of this same curved CG patch.
    // Both chi and eta are zero; this is not a different upstream trajectory.
    parameters.chi=0;serial->SetSolverInfo("RANSWallTransition","Fixed");
    RANSWallEvolutionPatch baseline(serial,graph,comm,parameters,p.Mapping());
    serial->SetSolverInfo("RANSWallTransition","Upstream");
    RANSWallEvolutionPatch disabled(serial,graph,comm,parameters,p.Mapping());
    // The generic independent outer data are original-AWF values at corners.
    auto common=s;const double K0=parameters.model.Q()*.16;
    std::vector<RANSWallEvolutionFunction> original;
    for(int b=0;b<p.Gauge();++b)original.push_back([&,b](const Vec &q,double){
        if(b==p.Epsilon()){
            const double K=K0+.02*q[0]+.03*q[1];return std::pow(K,1.5)/(2.55*parameters.geometry.height)+.01*(q[2]-parameters.geometry.height);}
        return p.fields[b]->Point(q).Value(s,p.offset[b]);});
    auto baselineData=data;baselineData.boundary[p.Epsilon()]=original[p.Epsilon()];
    common=baseline.Project(original,time);common=baseline.Lift(common,time,baselineData);
    const auto a=baseline.Residual(common,history,alpha,time,baselineData,false,false),b=disabled.Residual(common,history,alpha,time,baselineData,false,false);
    Vector delta(a.size());for(int i=0;i<static_cast<int>(a.size());++i)delta[i]=a[i]-b[i];
    o<<",\"equal_input_baseline\":{\"maximum\":"<<Max(delta)<<",\"state\":";Write(o,common);
    o<<",\"fixed_operator\":";Write(o,a);o<<",\"disabled_transition_operator\":";Write(o,b);
    const auto ba=baseline.Constraints(common,time,baselineData),bb=disabled.Constraints(common,time,baselineData);
    double traceError=0;Require(ba.size()==bb.size(),"baseline constraint size differs");
    for(auto [row,value]:ba)traceError=std::max(traceError,std::abs(value-bb.at(row)));
    o<<",\"constraint_maximum\":"<<traceError<<'}';
    // Full coupled actions include K/cutoff, every sensor and implicit trace
    // columns. The independent Python oracle also reconstructs weak loads.
    o<<",\"directions\":[";
    const auto essential=p.Constraints(s,time,data);int emitted=0;
    for(int block=0;block<=p.OuterPressure();++block)
    {
        if(block==p.InnerPressure())continue;
        Vector direction(s.size());
        for(int i=0;i<p.size[block];++i)
        {
            const int row=p.offset[block]+i;
            if(essential.count(row) && block!=p.Epsilon())continue;
            direction[row]=(block==p.Epsilon()?.001:.01)*std::sin(.71*(i+1));
        }
        if(emitted++){o<<',';}
        o<<"{\"block\":"<<block<<",\"vector\":";Write(o,direction);o<<",\"actions\":[";int level=0;
        for(double h:{8e-4,4e-4,2e-4,1e-4,5e-5,2e-6})
        {if(level++)o<<',';o<<"{\"step\":"<<h<<",\"value\":";Write(o,p.DirectionalAction(s,direction,history,alpha,time,data,h));o<<'}';}
        o<<"]}";
    }o<<"]}";
}
}
int main(int argc,char **argv)
{
    try
    {
        Require(argc==3,"usage: RANSCurvedWallProbe SESSION.xml REPORT.json");
        auto world=LibUtilities::SessionReader::CreateInstance(2,argv);world->InitSession();
        char program[]="curved-wall-probe";char *args[]{program,nullptr};
        std::vector<std::string> filenames{argv[1]};
        // SessionReader prepends optimiser sidecars to GetFilenames(). Keep
        // the original XML as the clone's first input on repeated runs.
        auto serial=LibUtilities::SessionReader::CreateInstance(1,args,filenames,LibUtilities::CommSerial::create(1,args));serial->InitSession();
        serial->SetSolverInfo("GlobalSysSoln","DirectFull");const auto graph=SpatialDomains::MeshGraphIO::Read(serial);const int d=graph->GetMeshDimension();
        RANSWallPatchParameters parameters;parameters.chi=0;
        serial->LoadParameter("RANSWallChi",parameters.chi,0.);
        serial->LoadParameter("RANSWallHeight",parameters.geometry.height,2.);serial->LoadParameter("RANSWallRoughness",parameters.geometry.roughness,.1);
        serial->LoadParameter("Kinvis",parameters.geometry.molecularNu,1.5e-5);serial->LoadParameter("RANSWallSampleHeight",parameters.geometry.sampleHeight,1.);
        serial->LoadParameter("RANSWallTailStart",parameters.geometry.tailStart,1.);serial->LoadParameter("RANSWallTangentialOrder",parameters.tangentialOrder,12);
        serial->LoadParameter("RANSWallNormalOrder",parameters.normalOrder,16);
        const auto mapping=LoadRANSWallNormalMap(serial,d,parameters.geometry.height);RANSWallEvolutionPatch patch(serial,graph,world->GetComm(),parameters,mapping);
        int constructorOnly=0;serial->LoadParameter("RANSWallProbeConstructorOnly",constructorOnly,0);
        if(constructorOnly)
        {
            Require(patch.TransitionEnabled(),"constructor audit requires Upstream mode");
            double dx,dy;serial->LoadParameter("RANSWallDirectionX",dx,1.);serial->LoadParameter("RANSWallDirectionY",dy,0.);
            const auto frame=mapping.At({patch.lower[0],patch.lower[1],0});Vec flow{};
            for(int a=0;a<3;++a)flow[a]=dx*frame.tangent1[a]+dy*frame.tangent2[a];
            const auto coordinate=EquilibriumAWFWallCoordinateDifferential(mapping,{patch.lower[0],patch.lower[1],0});
            double rate=0;for(int a=0;a<3;++a)for(int b=0;b<3;++b)rate+=coordinate.chartGradient[b]*frame.inverse[b][a]*flow[a];
            double departure=0;RANSWallEvolutionData zero(d);
            for(const auto &q:patch.Surface(0,2,parameters.geometry.height,1,4))if(patch.WallCoordinate(q.test.coordinate)<patch.TransitionControls().end)
                departure=std::max(departure,patch.Transition(patch.TransitionReference(),q.test.coordinate,0,zero).d2);
            if(world->GetComm()->GetRank()==0)
            {
                std::ofstream file(argv[2]);file<<std::setprecision(17)<<"{\"schema\":1,\"constructor_only\":true,\"dimension\":"<<d
                    <<",\"ranks\":"<<world->GetComm()->GetSize()<<",\"incoming_flow_cartesian\":";Write(file,flow);
                file<<",\"actual_wall_coordinate_rate\":"<<rate<<",\"reference_departure_maximum\":"<<departure<<"}\n";
                Require(static_cast<bool>(file),"cannot write constructor audit");
            }
            world->Finalise();return 0;
        }
        int innerFixture=0;serial->LoadParameter("RANSWallProbeInnerStressFixture",innerFixture,0);
        Require(innerFixture>=0 && innerFixture<=3,"unknown R2 inner-stress verification fixture");
        if(innerFixture)Require(d==2 && mapping.Identity() && !patch.TransitionEnabled(),
            "R2 inner-stress verification fixture is explicitly flat 2D with Fixed chi");
        const double K=parameters.model.Q()*.16,H=parameters.geometry.height,Z=patch.upper[2],time=.2,alpha=innerFixture?0.:5.3;
        using F=RANSWallEvolutionFunction;std::vector<F> functions;
        functions.push_back([](const Vec &p,double){return (.1+.03*p[0]+.02*p[1])*p[2];});
        if(d==3)functions.push_back([](const Vec &p,double){return .04*p[2]+.01*p[0]*p[1]*p[2];});
        functions.push_back([Z](const Vec &p,double){return .02*p[2]+.01*(p[0]+p[1])*p[2]*(Z-p[2]);});
        functions.push_back([K,rK=parameters.geometry.sampleHeight](const Vec &p,double){return K+.02*p[0]+.03*p[1]-.02*(p[2]-rK);});
        functions.push_back([K,H,mapping,rK=parameters.geometry.sampleHeight](const Vec &p,double){
            const auto g=mapping.At(p);const double sampled=K+.02*p[0]+.03*p[1]-.02*(rK/g.scale-rK);
            return std::pow(sampled,1.5)/(2.55*g.height)+.01*(p[2]-H);});
        functions.push_back([](const Vec &p,double){return -.01*p[0]+.015*p[1]+.02*p[2];});
        functions.push_back([](const Vec &p,double){return .04*p[0]+.03*p[1]+.05*p[2];});
        int invariant=0;serial->LoadParameter("RANSWallProbeInvariant",invariant,0);
        if(invariant)
        {
            Require(mapping.kind==RANSWallNormalMap::Kind::Cylinder && mapping.heightSlopeX==0 && mapping.heightSlopeY==0 && mapping.heightAmplitude==0,
                    "partition invariant probe requires a constant-H cylinder");
            for(int a=0;a<d;++a)functions[a]=[Z,d,a](const Vec &q,double){
                const double envelope=.02*q[0]*(1-q[0])*(d==3?q[1]*(1-q[1]):1.)*q[2]*(Z-q[2])*(Z-q[2]);
                return (a==0?1.:(a==d-1?-.3:.4))*envelope;};
            functions[d]=[K](const Vec&,double){return K;};
            functions[d+1]=[K,H](const Vec &q,double){return std::pow(K,1.5)/(2.55*H)+.01*(q[2]-H);};
        }
        if(innerFixture)
        {
            const double shear=innerFixture==2?.1:6.;
            functions[0]=[shear](const Vec &p,double){return shear*p[2];};
            if(innerFixture==3)
                functions[d]=[K,rK=parameters.geometry.sampleHeight](const Vec &p,double)
                    {return K*(.001+.999*std::pow(p[2]/rK,3));};
        }
        RANSWallEvolutionData data(d);for(int b=0;b<patch.InnerPressure();++b)data.boundary[b]=functions[b];
        if(!innerFixture)
        {
        data.innerBody[0]=[](const Vec &p,double){return .01+.02*p[0];};data.outerBody[0]=[](const Vec &p,double){return .04*p[0]+.01*p[2];};
        if(d==3){data.innerBody[1]=[](const Vec &p,double){return .015-.01*p[1]+.02*p[2];};data.outerBody[1]=[](const Vec &p,double){return .02*p[1]-.01*p[2];};}
        data.innerBody[d-1]=[](const Vec &p,double){return -.02+.03*p[2];};data.outerBody[d-1]=[](const Vec &p,double){return -.03+.02*p[2];};
        data.innerKSource=[](const Vec &p,double){return .03+.01*p[0]+.04*p[1]+.02*p[2];};
        data.outerKSource=[](const Vec &p,double){return -.02+.01*p[0]+.02*p[1]+.03*p[2];};
        data.outerEpsilonSource=[](const Vec &p,double){return .01*(1+p[0]+p[1]+p[2]);};
        }
        if(innerFixture)
        {
            // Fixed-chi lateral Dirichlet EPS and the dependent matching
            // trace share corner modes. Choose compatible fixed corner data
            // once; interior matching modes still rebuild lambda*Epre(H)
            // from every perturbed current k field in Constraints/Lift.
            const auto seed=patch.Project(functions,time);
            auto endpoint=[&](double x)
            {
                const Vec q{x,0,H};const auto p=patch.StateProfile(seed,q,time,data);
                return patch.WallDissipation(seed,q,*p);
            };
            const double left=endpoint(patch.lower[0]),right=endpoint(patch.upper[0]);
            functions[patch.Epsilon()]=[left,right,H,x0=patch.lower[0],length=patch.upper[0]-patch.lower[0]](const Vec &p,double)
                {return left+(right-left)*(p[0]-x0)/length+.01*(p[2]-H);};
            data.boundary[patch.Epsilon()]=functions[patch.Epsilon()];
        }
        const auto state=patch.Project(functions,time);auto history=state;for(int i=0;i<patch.DifferentialSize();++i)history[i]*=.97;
        const auto constraints=patch.Constraints(state,time,data);const auto lifted=patch.Lift(state,time,data);
        const auto operatorState=patch.TransitionEnabled()?lifted:state;
        const auto raw=patch.Residual(operatorState,history,alpha,time,data,false,false),inner=patch.Residual(operatorState,history,alpha,time,data,false,false,0),outer=patch.Residual(operatorState,history,alpha,time,data,false,false,1);
        Vector delta(raw.size());for(int i=0;i<static_cast<int>(raw.size());++i)delta[i]=inner[i]+outer[i]-raw[i];
        Require(Max(delta)<2e-10*std::max(1.,Max(raw)),"regional raw residuals do not sum");
        const auto roughness=patch.RoughnessBudget(operatorState,time,data),energy=patch.EnergyBudget(operatorState,history,alpha,time,data);
        std::ostringstream out;out<<std::setprecision(17)<<"{\"schema\":1,\"dimension\":"<<d<<",\"ranks\":"<<world->GetComm()->GetSize()
            <<",\"chi\":"<<parameters.chi<<",\"eta\":0,\"invariant_partition_fixture\":"<<(invariant?"true":"false")<<",\"time\":"<<time<<",\"alpha\":"<<alpha<<",\"function_coordinates\":\"Chart\",\"parameters\":";
        Write(out,Vector{parameters.geometry.molecularNu,H,parameters.geometry.roughness,parameters.geometry.sampleHeight,parameters.geometry.tailStart,K,parameters.model.Q(),parameters.model.Kappa()});
        double ustar=0;serial->LoadParameter("RANSWallFrictionVelocity",ustar,.4);
        out<<",\"model_variant\":\""<<patch.VariantName()<<"\",\"robust_parameters\":";Write(out,patch.RobustParameters());
        out<<",\"inner_stress_model\":\""<<patch.InnerStressName()<<"\",\"inner_stress_fixture\":"<<innerFixture
           <<",\"physical_signature\":"<<std::quoted(patch.Signature());
        if(innerFixture)
            out<<",\"upstream_transition\":false,\"unforced_alpha0_fixture\":true"
               <<",\"inner_stress_fixture_definitions\":{\"u\":\""<<(innerFixture==2?".1*r":"6*r")
               <<"\",\"k\":\""<<(innerFixture==3?"K0*(.001+.999*(r/rK)^3)":"base definition")
               <<"\",\"epsilon_seed\":\"linear interpolation of fixed compatible endpoint lambda*Epre(H), plus .01*(r-H); interior matching modes remain dependent\""
               <<",\"all_prescribed_body_and_scalar_sources\":0,\"alpha\":0,\"qualification\":\"operator verification only; outer state box is not required or enlarged\"}";
        out<<",\"ustar\":"<<ustar<<",\"robust_reference_flow\":";Write(out,patch.RobustReferenceFlow());
        out<<",\"robust_reference_normal\":";Write(out,patch.RobustReferenceNormal());
        out<<",\"map_parameters\":";Write(out,MapParameters(mapping));out<<",\"quadrature_orders\":";Write(out,Vector{static_cast<double>(parameters.tangentialOrder),static_cast<double>(parameters.normalOrder)});
        out<<",\"lower\":";Write(out,patch.lower);out<<",\"upper\":";Write(out,patch.upper);out<<",\"offsets\":[";
        for(int b=0;b<=patch.Gauge();++b){if(b)out<<',';out<<patch.offset[b];}out<<"],\"sizes\":[";
        for(int b=0;b<=patch.Gauge();++b){if(b)out<<',';out<<patch.size[b];}out<<"],\"cell_ownership\":[";int count=0;
        for(const auto &e:patch.fields[0]->Elements()){if(count++)out<<',';Write(out,Vector{static_cast<double>(e.geometry),static_cast<double>(patch.Owner(e.geometry))});}
        out<<"],\"fields\":[";for(int b=0;b<patch.Gauge();++b){if(b)out<<',';Describe(out,patch,b);}out<<"],\"constraints\":[";count=0;
        for(auto [row,value]:constraints){if(count++)out<<',';Write(out,Vector{static_cast<double>(row),value});}out<<"],\"free_differential\":[";count=0;
        for(int i=0;i<patch.DifferentialSize();++i)if(!constraints.count(i)){if(count++)out<<',';out<<i;}
        out<<"],\"state\":";Write(out,state);out<<",\"lifted_state\":";Write(out,lifted);out<<",\"history\":";Write(out,history);
        out<<",\"raw_residual\":";Write(out,raw);out<<",\"region_residuals\":[";Write(out,inner);out<<',';Write(out,outer);out<<"],\"regional_sum_error\":"<<Max(delta);
        out<<",\"boundary_flux_correction\":"<<patch.BoundaryFluxCorrection()<<",\"roughness_budget\":";Write(out,roughness);out<<",\"energy_budget\":";Write(out,energy);
        out<<",\"forcing_samples\":[";count=0;for(double x:{.13,.47,.91})
        {
            const Vec q{patch.lower[0]+x*(patch.upper[0]-patch.lower[0]),d==3?patch.lower[1]+.37*(patch.upper[1]-patch.lower[1]):0.,H};
            if(count++){out<<',';}
            out<<"{\"chart\":";Write(out,q);out<<",\"force\":";Write(out,patch.CartesianForcing(lifted,q,time,data));out<<",\"sampled_k\":"<<patch.SampleK(lifted,q)<<'}';
        }
        out<<"],\"bounded_quadrature\":";BoundedQuadrature(out,patch);out<<",\"physical_masses\":";Masses(out,patch);
        out<<",\"weak_metric\":";WeakMetric(out,patch);
        if(innerFixture)
        {
            const auto audit=patch.Audit(operatorState,history,alpha,time,data,ustar,true);
            out<<",\"constitutive_audit\":{";bool first=true;
            for(const auto &[name,value]:audit){if(!first)out<<',';first=false;out<<std::quoted(name)<<':';Write(out,value);}out<<'}';
            std::vector<F> changes(patch.Gauge(),[](const Vec&,double){return 0.;});
            changes[0]=[H](const Vec &p,double){return .002*p[2]*(H-p[2]);};
            changes[1]=[](const Vec &p,double){return .001*p[2];};
            changes[patch.K()]=[](const Vec &p,double){return .001*(1+.1*p[0]+.2*p[2]);};
            changes[patch.Epsilon()]=[H](const Vec &p,double){return .001*(p[2]-H);};
            changes[patch.InnerPressure()]=[](const Vec &p,double){return .002*p[0];};
            changes[patch.OuterPressure()]=[](const Vec &p,double){return -.001*p[0];};
            auto direction=patch.Project(changes,time);direction.back()=.001;
            out<<",\"r2_coupled_action\":{\"direction\":";Write(out,direction);out<<",\"samples\":[";int sample=0;
            for(double step:{1e-4,1e-5})
            {
                if(sample++)out<<',';auto plus=operatorState,minus=operatorState;
                for(int i=0;i<patch.Size();++i){plus[i]+=step*direction[i];minus[i]-=step*direction[i];}
                const auto lp=patch.Lift(plus,time,data),lm=patch.Lift(minus,time,data);
                out<<"{\"step\":"<<step<<",\"input_plus\":";Write(out,plus);out<<",\"input_minus\":";Write(out,minus);
                out<<",\"lifted_plus\":";Write(out,lp);out<<",\"lifted_minus\":";Write(out,lm);
                out<<",\"raw_plus\":";Write(out,patch.Residual(plus,history,alpha,time,data,false,false));
                out<<",\"raw_minus\":";Write(out,patch.Residual(minus,history,alpha,time,data,false,false));
                for(const auto &[name,value]:std::vector<std::pair<std::string,Vector>>{{"plus",plus},{"minus",minus}})
                {
                    out<<",\"constraints_"<<name<<"\":[";int row=0;
                    for(const auto &[id,boundary]:patch.Constraints(value,time,data))
                    {if(row++)out<<',';Write(out,Vector{static_cast<double>(id),boundary});}out<<']';
                }
                out<<",\"native_action\":";Write(out,patch.DirectionalAction(operatorState,direction,history,alpha,time,data,step));out<<'}';
            }out<<"]}";
        }
        if(patch.TransitionEnabled()){out<<",\"transition\":";TransitionAudit(out,patch,lifted,history,time,alpha,data,serial,graph,world->GetComm(),parameters);}
        else{out<<",\"interface\":";Interface(out,patch,lifted);out<<",\"inner_epsilon_diagnostic\":";InnerEpsilonReaction(out,patch,lifted,history,alpha);}
        out<<",\"data_definitions\":{\"coordinates\":\"q1,q2,rho; fixed physical rK uses rhoK=rK/scale\",\"u\":\"(.1+.03*q1+.02*q2)*rho\",\"v_3D\":\".04*rho+.01*q1*q2*rho\",\"w\":\".02*rho+.01*(q1+q2)*rho*(Z-rho)\",\"k\":\"K0+.02*q1+.03*q2-.02*(rho-rK)\",\"epsilon_seed\":\"(K0+.02*q1+.03*q2-.02*(rK/scale-rK))^1.5/(2.55*H(q))+.01*(rho-H0)\",\"pi_inner\":\"-.01*q1+.015*q2+.02*rho\",\"pi_outer\":\".04*q1+.03*q2+.05*rho\",\"body_inner_x\":\".01+.02*q1\",\"body_inner_y_3D\":\".015-.01*q2+.02*rho\",\"body_inner_z\":\"-.02+.03*rho\",\"body_outer_x\":\".04*q1+.01*rho\",\"body_outer_y_3D\":\".02*q2-.01*rho\",\"body_outer_z\":\"-.03+.02*rho\",\"source_inner_k\":\".03+.01*q1+.04*q2+.02*rho\",\"source_outer_k\":\"-.02+.01*q1+.02*q2+.03*rho\",\"source_outer_epsilon\":\".01*(1+q1+q2+rho)\",\"epsilon_diagnostic\":\"inner flux=U*Ew, source=0, reaction=alpha*(Ew-Ew_history)+div(U*Ew); no inner epsilon PDE\"}}\n";
        if(world->GetComm()->GetRank()==0){std::ofstream file(argv[2]);file<<out.str();Require(static_cast<bool>(file),"cannot write curved wall probe JSON");}
        world->Finalise();return 0;
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 1;}
}
