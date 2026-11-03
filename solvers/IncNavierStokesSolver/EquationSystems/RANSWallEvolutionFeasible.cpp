///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>

namespace Nektar
{
using namespace SolverUtils;
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;
namespace
{
bool Same(double a,double b){return std::abs(a-b)<1e-11;}
using Matrix=RANSWallAlgebra::DenseMatrix;
class FeasibleFactor
{
public:
    explicit FeasibleFactor(const Matrix &a):m_n(a.size()),m_lu(m_n*m_n),m_pivots(m_n)
    {
        for(int i=0;i<m_n;++i)for(int j=0;j<m_n;++j)m_lu[i+m_n*j]=a[i][j];
        int info=0;Lapack::Dgetrf(m_n,m_n,m_lu.data(),m_n,m_pivots.data(),info);
        Require(info==0,"ScalarFeasibleR2 linear factorization failure");
    }
    RANSWallAlgebra::Vector Apply(RANSWallAlgebra::Vector rhs) const
    {
        auto pivots=m_pivots;int info=0;
        Lapack::Dgetrs('N',m_n,1,m_lu.data(),m_n,pivots.data(),rhs.data(),m_n,info);
        Require(info==0,"ScalarFeasibleR2 linear solve failure");Maximum(rhs);return rhs;
    }
private:
    int m_n;RANSWallAlgebra::Vector m_lu;std::vector<int> m_pivots;
};
}

void RANSWallEvolutionPatch::InitialiseScalarLocations() const
{
    if(!m_scalarLocations.empty())return;
    std::set<RANSWallCoordinate> samplers,trace;
    auto add=[&](const std::string &name,const RANSWallGlobalPoint &v,int block,bool strict,int side=0)
    {
        m_scalarLocations.push_back({name,v,block,side,strict});
        samplers.insert({v.coordinate[0],v.coordinate[1],0});
    };
    auto pair=[&](const std::string &name,const RANSWallGlobalPoint &v,bool inner,int side=0)
    {
        add(name,v,K(),!inner,side);
        if(!inner)add(name,Point(Epsilon(),v.coordinate,Same(v.coordinate[2],m_parameters.geometry.height)?1:0),Epsilon(),true,side);
    };
    for(const auto &q:Quadrature(0))pair("mass_volume",q,q.coordinate[2]<m_parameters.geometry.height);
    for(const auto &q:m_bodyInner)pair("prescribed_source_inner",q.v,true);
    for(const auto &q:m_outer)pair("source_volume_outer",q.v,false);
    for(const auto &face:Surface(0,2,m_parameters.geometry.height,1,8))
    {
        pair("interface_inner",Point(K(),face.test.coordinate,-1),true,-1);
        pair("interface_outer",Point(K(),face.test.coordinate,1),false,1);
    }
    for(int e=0;e<fields[K()]->Elements().size();++e)
    {
        const auto &el=fields[K()]->Elements()[e];const auto exp=fields[K()]->Field()->GetExp(e);
        const bool inner=el.upper[2]<=m_parameters.geometry.height;
        for(int kind=0;kind<2;++kind)
        {
            const int dense=2*el.degree+(inner?7:3);
            const int nx=kind?dense:exp->GetNumPoints(0);
            const int ny=m_dimension==3?(kind?dense:exp->GetNumPoints(1)):1;
            const int nz=kind?dense:exp->GetNumPoints(m_dimension-1);
            for(int iz=0;iz<nz;++iz)for(int iy=0;iy<ny;++iy)for(int ix=0;ix<nx;++ix)
            {
                RANSWallCoordinate p{};const std::array<int,3> indices{ix,iy,iz};
                for(int a=0;a<3;++a)if(m_dimension==3 || a!=1)
                {
                    const int axis=m_dimension==2 && a==2?1:a;
                    const double z=kind?-1.+2.*indices[a]/(dense-1):exp->GetBasis(axis)->GetZ()[indices[a]];
                    p[a]=el.lower[a]+(z+1)*(el.upper[a]-el.lower[a])/2;
                }
                pair(kind?"reconstruction_dense":"reconstruction_native",
                    Transform(fields[K()]->ElementPoint(e,p)),inner,Same(p[2],m_parameters.geometry.height)?(inner?-1:1):0);
            }
        }
    }
    // EXACT evaluation nodes used by RANSWallGlobalField::Boundary, not the
    // source Gauss anchors or an assumed degree+1 interpolation grid.
    for(const auto &el:fields[Epsilon()]->Elements())if(Same(el.lower[2],m_parameters.geometry.height))
    {
        const int n=std::max(8,el.degree+3);
        const auto rule=LibUtilities::PointsManager()[LibUtilities::PointsKey(n,LibUtilities::eGaussLobattoLegendre)];
        const auto z=rule->GetZ();
        for(int iy=0;iy<(m_dimension==3?n:1);++iy)for(int ix=0;ix<n;++ix)
        {
            RANSWallCoordinate p{el.lower[0]+(z[ix]+1)*(el.upper[0]-el.lower[0])/2,
                m_dimension==3?el.lower[1]+(z[iy]+1)*(el.upper[1]-el.lower[1])/2:0.,m_parameters.geometry.height};
            trace.insert(p);samplers.insert({p[0],p[1],0});
            pair("matching_trace_nodes",Point(K(),p,1),false,1);
        }
    }
    for(const auto &ray:m_rays)samplers.insert({ray.coordinate[0],ray.coordinate[1],0});
    m_traceScalarAnchors.assign(trace.begin(),trace.end());
    m_scalarSamplerAnchors.assign(samplers.begin(),samplers.end());
    for(auto p:m_scalarSamplerAnchors)
    {
        const double scale=m_mapping.At(p).scale;
        p[2]=m_parameters.geometry.height;
        pair("physical_matching_transfer",Point(K(),p,1),false,1);
        p[2]=m_parameters.geometry.sampleHeight/scale;
        add("physical_sampled_K",Point(K(),p,Same(p[2],m_parameters.geometry.height)?-1:0),K(),true);
        if(m_transitionEnabled)for(double r:{m_sensorA,m_sensorB})
        {
            p[2]=r/scale;pair("physical_sensor_transfer",Point(K(),p),false);
        }
    }
}

RANSWallEvolutionPatch::ScalarDomainResult RANSWallEvolutionPatch::ScalarDomain(
    const Vector &state,double time,const RANSWallEvolutionData &data,bool retain) const
{
    ScalarDomainResult result;result.phase="coefficients";
    const double us=m_robustUstar,K0=m_parameters.model.Q()*us*us;
    const double E0=us*us*us/(m_parameters.model.Kappa()*m_parameters.geometry.height);
    auto check=[&](const std::string &name,const RANSWallGlobalPoint &q,int block,bool strict,
                   const Vector &values,int side=0,bool dependent=false,double weight=0.)
    {
        const double value=q.Value(values,offset[block]);
        const double scale=block==K()?K0:E0*m_parameters.geometry.height/Distance(q.coordinate);
        const double margin=value/scale;++result.counts[name];
        if(block==K())
        {
            if(name=="physical_sampled_K")result.minimumSampledK=std::min(result.minimumSampledK,value);
            if(Distance(q.coordinate)<=m_mapping.At(q.coordinate).height && side<=0)
                result.minimumInnerK=std::min(result.minimumInnerK,value);
            else result.minimumOuterK=std::min(result.minimumOuterK,value);
        }
        else if(dependent)result.minimumTraceEpsilon=std::min(result.minimumTraceEpsilon,value);
        else result.minimumOuterEpsilon=std::min(result.minimumOuterEpsilon,value);
        if(margin<result.minimumMargin || !std::isfinite(value))
        {result.minimumMargin=margin;result.limitingPoint=q.coordinate;result.limitingInventory=name;}
        if(retain)result.observations.push_back({name,q.coordinate,block,q.geometry,side,value,scale,weight,strict,dependent});
        if(!std::isfinite(value) || (strict?value<=0:value<0))
        {
            result.failure=block==K()?(strict?"nonpositive_K_or_outer_k":"negative_inner_k"):"nonpositive_epsilon";
            result.limitingPoint=q.coordinate;result.limitingInventory=name;return false;
        }
        return true;
    };
    try
    {
        Require(state.size()==static_cast<std::size_t>(Size()),"scalar-domain state size differs");
        for(double value:state)if(!std::isfinite(value)){result.failure="nonfinite_coefficient";return result;}
        InitialiseScalarLocations();
        result.effectiveState=ExternalLift(state,time,data);result.phase="fixed_fields_and_samplers";
        for(const auto &point:m_scalarLocations)
            if(!check(point.inventory,point.row,point.block,point.strict,result.effectiveState,point.side))return result;
        struct Pending
        {
            RANSWallCoordinate anchor;
            std::shared_ptr<EquilibriumAWFProfiles> profile;
            bool source=false;
        };
        std::vector<Pending> pending;
        result.phase="moving_ray_k_before_normalization";
        auto ray=[&](const RANSWallCoordinate &p,const std::shared_ptr<EquilibriumAWFProfiles> &profile,
                     const std::string &name,bool source,double horizontalWeight)
        {
            const double scale=m_mapping.At(p).scale;
            result.minimumCutoff=std::min(result.minimumCutoff,profile->Reference().cutoff);
            result.maximumCutoff=std::max(result.maximumCutoff,profile->Reference().cutoff);
            result.profileBranches.push_back(profile->Reference().weight==0?'B':profile->Reference().weight==1?'E':'M');
            for(const auto &q:profile->Quadrature())
            {
                auto point=p;point[2]=q.r/scale;
                if(!check(name,Point(K(),point),K(),false,result.effectiveState,0,false,
                    source?horizontalWeight*q.weight*q.jacobian:0.))return false;
            }
            pending.push_back({p,profile,source});return true;
        };
        // Check every physical k BEFORE calling WallNormalization or a
        // dependent trace. Current profiles rebuild trial-specific cutoffs.
        for(const auto &r:m_rays)
            if(!ray(r.coordinate,StateProfile(result.effectiveState,r.coordinate,time,data),
                    "volume_normalization_ray",true,r.weight))return result;
        for(const auto &p:m_traceScalarAnchors)
        {
            if(!ray(p,Profile(SampleK(result.effectiveState,p),p,0),"trace_baseline_normalization_ray",false,0))return result;
            if(!ray(p,StateProfile(result.effectiveState,p,time,data),
                    "trace_current_normalization_ray",false,0))return result;
        }
        result.phase="normalization_and_dependent_trace";
        for(const auto &p:pending)if(m_cappedDestruction)
        {
            const auto n=WallNormalization(result.effectiveState,p.anchor,*p.profile);
            result.minimumB=std::min(result.minimumB,n.b);result.minimumLambda=std::min(result.minimumLambda,n.lambda);
            result.profileBranches.push_back(n.capped?'C':'P');
        }
        // Retain independent epsilon-H unknowns exactly as Residual does for
        // Upstream. A second vector checks the projected target at the trace;
        // it is NEVER substituted into that coupled residual/Newton iterate.
        result.traceState=result.effectiveState;
        for(auto [row,value]:Constraints(state,time,data))result.traceState[row]=value;
        if(!m_transitionEnabled)result.effectiveState=result.traceState;
        for(const auto &point:m_scalarLocations)if(point.block==Epsilon())
        {
            if(!check("actual_lifted_outer_epsilon",point.row,Epsilon(),true,result.effectiveState,point.side))return result;
            if(Same(point.row.coordinate[2],m_parameters.geometry.height) &&
                !check("projected_dependent_trace",point.row,Epsilon(),true,result.traceState,1,true))return result;
        }
        result.phase="complete";result.admissible=true;return result;
    }
    catch(const std::exception &error)
    {result.failure=std::string(J::ExceptionType(error))+": "+error.what();return result;}
}

std::string RANSWallEvolutionPatch::ScalarDomainJSON(const ScalarDomainResult &d,bool locations)
{
    std::vector<std::pair<std::string,std::string>> counts;
    for(auto [name,count]:d.counts)counts.push_back({name,J::Number(count)});
    std::vector<std::pair<std::string,std::string>> fields{
        {"admissible",J::Boolean(d.admissible)},{"scope",J::Quote("declared finite scalar evaluation sets only")},
        {"phase",J::Quote(d.phase)},{"failure",J::Quote(d.failure)},{"minimum_dimensionless_margin",J::Number(d.minimumMargin)},
        {"limiting_inventory",J::Quote(d.limitingInventory)},{"limiting_point",J::Numbers(d.limitingPoint)},
        {"minimum_inner_k",J::Number(d.minimumInnerK)},{"minimum_outer_k",J::Number(d.minimumOuterK)},
        {"minimum_outer_epsilon",J::Number(d.minimumOuterEpsilon)},{"minimum_sampled_K",J::Number(d.minimumSampledK)},
        {"minimum_projected_trace_epsilon",J::Number(d.minimumTraceEpsilon)},
        {"minimum_b",J::Number(d.minimumB)},{"minimum_lambda",J::Number(d.minimumLambda)},
        {"cutoff_range",J::Numbers(Vector{d.minimumCutoff,d.maximumCutoff})},
        {"inventory_counts",J::Object(counts)},{"profile_branch_pattern",J::Quote(d.profileBranches)}};
    if(locations)
    {
        std::vector<std::string> points;
        for(const auto &p:d.observations)points.push_back(J::Object({
            {"inventory",J::Quote(p.inventory)},{"coordinate",J::Numbers(p.coordinate)},
            {"geometry",J::Number(p.geometry)},{"block",J::Number(p.block)},{"side",J::Number(p.side)},
            {"value",J::Number(p.value)},{"scale",J::Number(p.scale)},{"weight",J::Number(p.weight)},
            {"strict",J::Boolean(p.strict)},{"dependent_trace",J::Boolean(p.dependentTrace)}}));
        fields.push_back({"observations",J::List(points)});
        fields.push_back({"effective_state",J::Numbers(d.effectiveState)});
        fields.push_back({"trace_target_state",J::Numbers(d.traceState)});
    }
    return J::Object(fields);
}

RANSWallEvolutionPatch::WallMatrix RANSWallEvolutionPatch::ScalarFeasibleJacobian(
    const Vector &s,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,
    std::vector<std::string> *metadata,bool scalarOnly) const
{
    Require(Size()<=3000,"ScalarFeasibleR2 retains the 3000-unknown dense guard");
    const auto domain=ScalarDomain(s,time,data);Require(domain.admissible,"Jacobian base is outside scalar domain");
    const auto constraints=Constraints(s,time,data);const auto base=Residual(s,history,0,time,data);
    const double us=m_robustUstar,K0=m_parameters.model.Q()*us*us;
    const double E0=us*us*us/(m_parameters.model.Kappa()*m_parameters.geometry.height);
    WallMatrix jac(Size(),Vector(Size()));if(metadata)metadata->clear();
    for(int column=0;column<Size();++column)
    {
        if(column%128==0)std::cout<<"EP-AWF ScalarFeasibleR2 Jacobian column "<<column<<" / "<<Size()<<std::endl;
        int block=0;while(block<Gauge() && column>=offset[block+1])++block;
        if(scalarOnly && (block<K() || block>Epsilon()))
        {
            jac[column][column]=1.;
            if(metadata)metadata->push_back(J::Object({{"column",J::Number(column)},
                {"stencil",J::Quote("frozen_scalar_diagnostic_not_a_coupled_column")}}));
            continue;
        }
        if(block>=InnerPressure())
        {
            Vector direction(Size());direction[column]=1.;
            const auto action=PressureDirectionalAction(s,direction,time,data);
            for(int i=0;i<Size();++i)jac[i][column]=action[i];
            if(metadata)metadata->push_back(J::Object({{"column",J::Number(column)},
                {"stencil",J::Quote("exact_affine_and_full_nonlinear_pressure_action")}}));
            continue;
        }
        if(constraints.count(column) && !(block==Epsilon() && m_transitionEnabled && m_traceIds.count(column-offset[Epsilon()])))
        {
            jac[column][column]=1.;
            if(metadata)metadata->push_back(J::Object({{"column",J::Number(column)},{"stencil",J::Quote("lifted_constraint_identity")}}));
            continue;
        }
        const double scale=block<m_dimension?us:block==K()?K0:E0;
        const double initial=2e-6*std::max(scale,std::abs(s[column]));
        Vector plusValue,minusValue,one;double selected=0;int side=0;bool central=false;
        std::vector<std::string> probes;
        for(double shrink:{1.,.25,.0625,.015625,.00390625})
        {
            const double h=initial*shrink;auto p=s,m=s;p[column]+=h;m[column]-=h;
            auto pd=ScalarDomain(p,time,data),md=ScalarDomain(m,time,data);
            bool pv=pd.admissible,mv=md.admissible;Vector a,b;
            try{if(pv)a=Residual(p,history,0,time,data);}catch(const std::exception &error){pv=false;pd.failure=error.what();}
            try{if(mv)b=Residual(m,history,0,time,data);}catch(const std::exception &error){mv=false;md.failure=error.what();}
            probes.push_back(J::Object({{"step",J::Number(h)},{"plus_domain",ScalarDomainJSON(pd)},
                {"minus_domain",ScalarDomainJSON(md)},{"plus_evaluated",J::Boolean(pv)},{"minus_evaluated",J::Boolean(mv)}}));
            if(pv && mv){plusValue=std::move(a);minusValue=std::move(b);selected=h;central=true;break;}
            if(pv || mv){one=pv?std::move(a):std::move(b);side=pv?1:-1;selected=h;}
        }
        Require(central || side,"ScalarFeasibleR2 no admissible coordinate action at column "+std::to_string(column));
        std::string stencil=central?"scaled_central":"admissible_one_sided_first_order";
        if(central)for(int i=0;i<Size();++i)jac[i][column]=(plusValue[i]-minusValue[i])/(2*selected);
        else
        {
            auto p=s;p[column]+=2*side*selected;Vector second;bool secondValid=false;
            const auto d=ScalarDomain(p,time,data);
            if(d.admissible)try{second=Residual(p,history,0,time,data);secondValid=true;}catch(const std::exception &){}
            for(int i=0;i<Size();++i)jac[i][column]=secondValid?
                side*(-3*base[i]+4*one[i]-second[i])/(2*selected):side*(one[i]-base[i])/selected;
            if(secondValid)stencil="admissible_one_sided_second_order";
        }
        // Time mass is affine in retained Upstream unknowns. Assemble it
        // exactly instead of subtracting large, unchanged time loads.
        const int local=column-offset[block];
        for(int i=0;i<size[block];++i)if(!constraints.count(offset[block]+i))
            jac[offset[block]+i][column]+=alpha*m_mass[block][i][local];
        if(metadata)metadata->push_back(J::Object({{"column",J::Number(column)},
            {"stencil",J::Quote(stencil)},{"step",J::Number(selected)},
            {"side",J::Number(central?0:side)},{"probes",J::List(probes)}}));
    }
    if(scalarOnly)for(int row=0;row<Size();++row)
        if(row<offset[K()] || row>=offset[InnerPressure()])
        {std::fill(jac[row].begin(),jac[row].end(),0.);jac[row][row]=1.;}
    return jac;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::SolveScalarFeasible(
    Vector s,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,
    bool condensed,int &iterations,double tolerance,int maximumIterations,bool scalarOnly) const
{
    Require(!condensed,"ScalarFeasibleR2 currently uses the full coupled dense solve, without condensation");
    Require(m_dimension==2 && !Mapped() && m_transitionEnabled && m_cappedDestruction && m_comm->GetSize()==1,
        "ScalarFeasibleR2 requires serial flat 2D Upstream+Capped");
    Require(Size()<=3000 && alpha>0,"invalid bounded scalar-feasible implicit stage");
    std::vector<std::string> records;std::string terminal="started";
    const int solveNumber=++m_feasibleSolveCount;
    auto domain=ScalarDomain(s,time,data,true);Vector finalResidual,finalNorms;
    const double us=m_robustUstar,K0=m_parameters.model.Q()*us*us;
    const double E0=us*us*us/(m_parameters.model.Kappa()*m_parameters.geometry.height);
    auto publish=[&](bool converged)
    {
        const auto finalConstraints=domain.admissible?Constraints(s,time,data):std::map<int,double>{};
        Require(m_feasibleTraceStarted || !std::ifstream(m_feasibleTraceFile).good(),"refusing to overwrite an existing scalar-feasible trace");
        std::ofstream out(m_feasibleTraceFile,std::ios::app);Require(out.good(),"cannot open scalar-feasible trace");
        out<<J::Object({{"schema","1"},{"policy",J::Quote("ScalarFeasibleR2:v1")},
            {"inner_stress_model",J::Quote(InnerStressName())},{"time",J::Number(time)},{"alpha",J::Number(alpha)},
            {"tolerance",J::Number(tolerance)},{"maximum_iterations",J::Number(maximumIterations)},
            {"converged",J::Boolean(converged)},{"terminal",J::Quote(terminal)},
            {"frozen_velocity_pressure_diagnostic",J::Boolean(scalarOnly)},
            {"history_advanced","false"},{"residual_unchanged","true"},
            {"jacobian_rebuild_policy",J::Quote("fresh every Newton iteration; exact mass and complete analytic pressure actions")},
            {"domain",ScalarDomainJSON(domain,true)},{"state",J::Numbers(s)},{"history",J::Numbers(history)},
            {"constraints",J::Constraints(finalConstraints)},
            {"residual",J::Numbers(finalResidual)},{"scaled_block_norms",J::Numbers(finalNorms)},
            {"iterations",J::List(records)}})<<'\n';out.close();Require(out.good(),"cannot finish scalar-feasible trace");
        m_feasibleTraceStarted=true;
    };
    try
    {
        Require(domain.admissible,"ScalarFeasibleR2 initial iterate is outside declared scalar domain: "+domain.failure);
        int stagnant=0;
        for(iterations=0;iterations<=maximumIterations;++iterations)
        {
            const auto constraints=Constraints(s,time,data);
            const auto r=Residual(s,history,alpha,time,data);const auto norms=Norms(r,constraints,us);
            const double error=scalarOnly?std::max({norms[K()],norms[Epsilon()],norms.back()}):Maximum(norms);
            finalResidual=r;finalNorms=norms;
            if(error<tolerance)
            {
                domain=ScalarDomain(s,time,data,true);Require(domain.admissible,"converged iterate lost its scalar domain");
                terminal=scalarOnly?"frozen_scalar_residual_converged_not_coupled":"true_coupled_residual_converged";
                m_lastFeasibleState=s;publish(true);return s;
            }
            Require(iterations<maximumIterations,"ScalarFeasibleR2 Newton iteration limit");
            std::cout<<"EP-AWF ScalarFeasibleR2 fresh Jacobian iteration "<<iterations<<" alpha "<<alpha<<" residual "<<error<<std::endl;
            std::vector<std::string> columns;const auto jac=ScalarFeasibleJacobian(s,history,alpha,time,data,&columns,scalarOnly);
            const std::string matrixFile=m_feasibleTraceFile+".solve"+std::to_string(solveNumber)+".J"+std::to_string(iterations)+".bin";
            Require(!std::ifstream(matrixFile).good(),"refusing to overwrite a scalar-feasible Jacobian");
            std::ofstream matrixOut(matrixFile,std::ios::binary);
            for(const auto &row:jac)matrixOut.write(reinterpret_cast<const char*>(row.data()),row.size()*sizeof(double));
            matrixOut.close();Require(matrixOut.good(),"cannot preserve scalar-feasible Jacobian");
            auto rhs=r;for(int i=0;i<Size();++i)rhs[i]=scalarOnly && (i<offset[K()] || i>=offset[InnerPressure()])?0.:-rhs[i];
            auto step=FeasibleFactor(jac).Apply(std::move(rhs));auto linear=Multiply(jac,step);
            const auto jacStep=linear;
            for(int i=0;i<Size();++i)linear[i]+=scalarOnly && (i<offset[K()] || i>=offset[InnerPressure()])?0.:r[i];
            const auto linearNorms=Norms(linear,constraints,us);
            Require(Maximum(linearNorms)<1e-8*std::max(1.,error),"ScalarFeasibleR2 inaccurate linear solve");
            std::vector<std::string> actions;
            const auto spatialBase=Residual(s,history,0,time,data);
            for(const std::string kind:{"velocity","free_k","epsilon_with_trace","pressure_and_gauge","joint"})
            {
                if(scalarOnly && kind!="free_k" && kind!="epsilon_with_trace")continue;
                Vector direction(Size());
                for(int b=0;b<=Gauge();++b)for(int i=0;i<size[b];++i)
                {
                    const int column=offset[b]+i;
                    const bool selected=kind=="joint" || (kind=="velocity" && b<K()) ||
                        (kind=="free_k" && b==K()) || (kind=="epsilon_with_trace" && b==Epsilon()) ||
                        (kind=="pressure_and_gauge" && b>=InnerPressure());
                    const double scale=b<m_dimension?us:b==K()?K0:b==Epsilon()?E0:b==Gauge()?us/m_parameters.geometry.height:us*us;
                    if(selected && (!constraints.count(column) || (b==Epsilon() && m_traceIds.count(i))))
                        direction[column]=scale*std::sin(1.317*(column+1));
                }
                if(kind=="free_k" || kind=="joint")
                {
                    const auto bubble=fields[K()]->Project([&](const RANSWallCoordinate &p)
                    {
                        if(p[2]>m_parameters.geometry.height)return 0.;
                        const double x=(p[0]-lower[0])/(upper[0]-lower[0]),r=p[2]/m_parameters.geometry.height;
                        return 16*K0*x*(1-x)*r*(1-r);
                    });
                    for(int i=0;i<size[K()];++i)direction[offset[K()]+i]=constraints.count(offset[K()]+i)?0.:bubble[i];
                }
                Vector massAction(Size());
                for(int b=0;b<InnerPressure();++b)
                {
                    Vector part(direction.begin()+offset[b],direction.begin()+offset[b]+size[b]);
                    const auto load=Multiply(m_mass[b],part);
                    for(int i=0;i<size[b];++i)if(!constraints.count(offset[b]+i))massAction[offset[b]+i]=alpha*load[i];
                }
                std::vector<std::string> levels;
                for(double h:{2e-5,2e-6,5e-7})
                {
                    std::vector<std::string> sides;
                    for(int sign:{1,-1})
                    {
                        auto point=s;for(int i=0;i<Size();++i)point[i]+=sign*h*direction[i];
                        const auto checked=ScalarDomain(point,time,data);
                        std::vector<std::pair<std::string,std::string>> observation{
                            {"domain",ScalarDomainJSON(checked)},{"evaluated",J::Boolean(checked.admissible)}};
                        if(checked.admissible)
                        {
                            try
                            {
                                observation.push_back({"spatial_residual",J::Numbers(Residual(point,history,0,time,data))});
                                observation.push_back({"constraints",J::Constraints(Constraints(point,time,data))});
                            }
                            catch(const std::exception &errorValue)
                            {
                                observation[1].second="false";
                                observation.push_back({"exception_type",J::Quote(J::ExceptionType(errorValue))});
                                observation.push_back({"exception_message",J::Quote(errorValue.what())});
                            }
                        }
                        sides.push_back(J::Object(observation));
                    }
                    levels.push_back(J::Object({{"step",J::Number(h)},{"plus",sides[0]},{"minus",sides[1]}}));
                }
                actions.push_back(J::Object({{"kind",J::Quote(kind)},{"direction",J::Numbers(direction)},
                    {"jacobian_action",J::Numbers(Multiply(jac,direction))},{"exact_mass_action",J::Numbers(massAction)},
                    {"spatial_base",J::Numbers(spatialBase)},{"levels",J::List(levels)}}));
            }
            // External lifting is affine here. The retained matching EPS
            // coefficients stay independent; the target prediction below is
            // only an additional bound and is rebuilt nonlinearly per trial.
            auto moved=s;for(int i=0;i<Size();++i)moved[i]+=step[i];
            moved=ExternalLift(moved,time,data);Vector actualStep(Size());
            for(int i=0;i<Size();++i)actualStep[i]=moved[i]-domain.effectiveState[i];
            auto traceStep=actualStep;
            for(int id:m_traceIds)traceStep[offset[Epsilon()]+id]=step[offset[Epsilon()]+id]-jacStep[offset[Epsilon()]+id];
            double maximum=1;ScalarObservation limiting;bool limited=false;
            for(const auto &p:domain.observations)
            {
                const auto row=Point(p.block,p.coordinate,p.side);
                const double delta=row.Value(p.dependentTrace?traceStep:actualStep,offset[p.block]);
                if(delta<0)
                {
                    const double bound=p.value/(-delta);
                    if(bound<maximum){maximum=std::max(0.,bound);limiting=p;limited=true;}
                }
            }
            double fraction=maximum<1?.995*maximum:1.;
            std::vector<std::string> trials;bool accepted=false;double nextError=error;Vector acceptedState;
            ScalarDomainResult acceptedDomain;
            for(int backtrack=0;backtrack<32 && fraction>=std::ldexp(1.,-40);++backtrack,fraction*=.5)
            {
                auto trial=s;for(int i=0;i<Size();++i)trial[i]+=fraction*step[i];
                auto candidateDomain=ScalarDomain(trial,time,data,true);
                if(!candidateDomain.admissible)
                {
                    trials.push_back(J::Object({{"fraction",J::Number(fraction)},
                        {"classification",J::Quote("scalar_domain_rejection_before_residual")},
                        {"domain",ScalarDomainJSON(candidateDomain)}}));continue;
                }
                try
                {
                    const auto next=Residual(trial,history,alpha,time,data);
                    const auto nextConstraints=Constraints(trial,time,data);
                    const auto nextNorms=Norms(next,nextConstraints,us);
                    nextError=scalarOnly?std::max({nextNorms[K()],nextNorms[Epsilon()],nextNorms.back()}):Maximum(nextNorms);
                    const bool descent=nextError<(1.-1e-4*fraction)*error;
                    auto rechecked=ScalarDomain(trial,time,data,true);
                    const bool accept=descent && rechecked.admissible;
                    trials.push_back(J::Object({{"fraction",J::Number(fraction)},
                        {"classification",J::Quote(accept?"feasible_accepted_descent":"valid_non_descent_or_recheck_failure")},
                        {"state",J::Numbers(trial)},{"residual",J::Numbers(next)},
                        {"constraints",J::Constraints(nextConstraints)},
                        {"scaled_block_norms",J::Numbers(nextNorms)},{"domain",ScalarDomainJSON(rechecked)},
                        {"profile_branches_changed",J::Boolean(domain.profileBranches!=rechecked.profileBranches)}}));
                    if(accept){accepted=true;acceptedState=std::move(trial);acceptedDomain=std::move(rechecked);break;}
                }
                catch(const std::exception &errorValue)
                {
                    trials.push_back(J::Object({{"fraction",J::Number(fraction)},
                        {"classification",J::Quote("evaluable_scalar_state_residual_exception")},
                        {"type",J::Quote(J::ExceptionType(errorValue))},{"message",J::Quote(errorValue.what())}}));
                }
            }
            double stepNorm=0;
            for(int b=0;b<=Gauge();++b)for(int i=0;i<size[b];++i)
            {
                const double scale=b<m_dimension?us:b==K()?K0:b==Epsilon()?E0:b==Gauge()?us/m_parameters.geometry.height:us*us;
                stepNorm=std::max(stepNorm,std::abs(step[offset[b]+i])/scale);
            }
            records.push_back(J::Object({{"iteration",J::Number(iterations)},{"state",J::Numbers(s)},
                {"constraints",J::Constraints(constraints)},
                {"residual",J::Numbers(r)},{"scaled_block_norms",J::Numbers(norms)},
                {"domain",ScalarDomainJSON(domain,true)},{"step",J::Numbers(step)},
                {"dimensionless_step",J::Number(stepNorm)},{"linear_scaled_block_norms",J::Numbers(linearNorms)},
                {"jacobian_binary_file",J::Quote(matrixFile)},
                {"jacobian_binary_shape",J::Numbers(Vector{double(Size()),double(Size())})},
                {"jacobian_binary_encoding",J::Quote("row-major native IEEE754 binary64; byte order in host receipt")},
                {"current_directional_checks",J::List(actions)},
                {"fixed_location_maximum_fraction",J::Number(maximum)},
                {"boundary_limited",J::Boolean(limited)},{"limiting_inventory",J::Quote(limiting.inventory)},
                {"limiting_point",J::Numbers(limiting.coordinate)},
                {"limiting_value",J::Number(limiting.value)},{"accepted_fraction",J::Number(accepted?fraction:0)},
                {"columns",J::List(columns)},{"trials",J::List(trials)},{"accepted",J::Boolean(accepted)}}));
            if(!accepted){terminal=maximum<std::ldexp(1.,-40)?"feasible_step_collapsed":"feasible_line_search_exhausted";throw std::runtime_error(terminal);}
            s=std::move(acceptedState);domain=std::move(acceptedDomain);
            stagnant=nextError>error*(1.-1e-4)?stagnant+1:0;
            if(stagnant>=3){terminal="feasible_residual_stagnation";throw std::runtime_error(terminal);}
        }
        throw std::runtime_error("unreachable scalar-feasible solve");
    }
    catch(const std::exception &error)
    {
        m_lastFeasibleState=s;
        if(terminal=="started")terminal=std::string(J::ExceptionType(error))+": "+error.what();
        domain=ScalarDomain(s,time,data,true);
        if(domain.admissible)
        {
            try{finalResidual=Residual(s,history,alpha,time,data);finalNorms=Norms(finalResidual,Constraints(s,time,data),us);}
            catch(const std::exception &){finalResidual.clear();finalNorms.clear();}
        }
        try{publish(false);}catch(const std::exception &exportError){std::cerr<<"ScalarFeasibleR2 trace export: "<<exportError.what()<<std::endl;}
        throw;
    }
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::SolveFrozenScalars(Vector state,
    const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,
    int &iterations,double tolerance,int maximumIterations) const
{
    Require(m_nonlinearPolicy=="ScalarFeasibleR2","frozen scalar diagnostic requires ScalarFeasibleR2");
    return SolveScalarFeasible(std::move(state),history,alpha,time,data,false,iterations,tolerance,maximumIterations,true);
}

RANSWallEvolutionPatch::Diagnostic RANSWallEvolutionPatch::CompleteRestrictedAudit(
    const Vector &state,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data) const
{
    const auto domain=ScalarDomain(state,time,data,true);Diagnostic audit;
    audit["r2_scalar_domain"]={domain.admissible?1.:0.,domain.minimumMargin,domain.minimumInnerK,
        domain.minimumOuterK,domain.minimumOuterEpsilon,domain.minimumSampledK,
        domain.minimumTraceEpsilon,domain.minimumB,domain.minimumLambda};
    if(!domain.admissible){audit["r2_qualified"]={0};return audit;}
    audit=Audit(domain.effectiveState,history,alpha,time,data,m_robustUstar,true);
    audit["r2_scalar_domain"]={1,domain.minimumMargin,domain.minimumInnerK,domain.minimumOuterK,
        domain.minimumOuterEpsilon,domain.minimumSampledK,domain.minimumTraceEpsilon,domain.minimumB,domain.minimumLambda};
    const double K0=m_parameters.model.Q()*m_robustUstar*m_robustUstar;
    Vector outer{0,1e300,-1e300,1e300,-1e300,domain.minimumSampledK/K0,-1e300};
    double innerCov=1e300,outerCov=1e300;Vector innerPoint,outerPoint;
    for(const auto &p:domain.observations)
    {
        if(p.block!=K() && !p.dependentTrace)continue;
        const double r=Distance(p.coordinate);const bool inner=r<m_mapping.At(p.coordinate).height-1e-11 ||
            (Same(r,m_mapping.At(p.coordinate).height) && p.side<0);
        auto row=Point(0,p.coordinate,p.side);
        // Preserve both one-sided native element gradients at internal
        // interfaces, rather than reselecting the lowest geometry ID.
        for(int e=0;e<fields[0]->Elements().size();++e)
            if(fields[0]->Elements()[e].geometry==p.geometry)
            {row=Transform(fields[0]->ElementPoint(e,p.coordinate));break;}
        RANSTensor G{};
        for(int a=0;a<m_dimension;++a)G[PhysicalAxis(a)]=row.Gradient(domain.effectiveState,offset[a]);
        const double k=Point(K(),p.coordinate,p.side).Value(domain.effectiveState,offset[K()]);
        if(inner)
        {
            const auto profile=StateProfile(domain.effectiveState,p.coordinate,time,data);
            const double nt=profile->Original().EddyViscosity(r);const auto law=InnerStress(k,nt,G);
            const double cov=J::Eigenvalues(law.reynoldsStress)[0];
            if(cov<innerCov){innerCov=cov;innerPoint={p.coordinate[0],p.coordinate[1],p.coordinate[2],k,profile->K(),nt,law.a,law.c,cov};}
            continue;
        }
        const double epsilon=p.dependentTrace?p.value:Point(Epsilon(),p.coordinate,Same(r,m_mapping.At(p.coordinate).height)?1:0).Value(domain.effectiveState,offset[Epsilon()]);
        const auto law=OuterClosure(k,epsilon,G,p.coordinate);const double cov=J::Eigenvalues(law.reynoldsStress)[0];
        if(cov<outerCov){outerCov=cov;outerPoint={p.coordinate[0],p.coordinate[1],p.coordinate[2],k,epsilon,cov};}
        const double g0=m_robustUstar/(m_parameters.model.Kappa()*r);
        const double e0=m_robustUstar*m_robustUstar*m_robustUstar/(m_parameters.model.Kappa()*r);
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            outer[0]=std::max(outer[0],std::abs(G[a][b]/g0-m_referenceFlow[a]*m_referenceNormal[b]));
        outer[1]=std::min(outer[1],k/K0);outer[2]=std::max(outer[2],k/K0);
        outer[3]=std::min(outer[3],epsilon/e0);outer[4]=std::max(outer[4],epsilon/e0);
    }
    for(const auto &p:domain.observations)if(p.inventory=="physical_sampled_K")outer[6]=std::max(outer[6],p.value/K0);
    const bool box=outer[0]<=.001 && outer[1]>=.99 && outer[2]<=1.01 && outer[3]>=.98 && outer[4]<=1.02 &&
        outer[5]>=.99 && outer[6]<=1.01;
    const bool covariance=innerCov>=-1e-10*K0 && outerCov>=-1e-10*K0;
    audit["r2_outer_range"]=outer;audit["r2_covariance"]={innerCov,outerCov};
    audit["r2_inner_minimum_point"]=innerPoint;audit["r2_outer_minimum_point"]=outerPoint;
    audit["r2_qualified"]={box && covariance && domain.minimumLambda>=.5?1.:0.};
    return audit;
}
} // namespace Nektar
