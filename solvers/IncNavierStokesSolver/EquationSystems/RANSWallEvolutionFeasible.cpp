///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>
#include <limits>
#include <numeric>
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
    RANSWallAlgebra::Vector Apply(RANSWallAlgebra::Vector rhs,bool transpose=false) const
    {
        auto pivots=m_pivots;int info=0;
        Lapack::Dgetrs(transpose?'T':'N',m_n,1,m_lu.data(),m_n,pivots.data(),rhs.data(),m_n,info);
        Require(info==0,"ScalarFeasibleR2 linear solve failure");Maximum(rhs);return rhs;
    }
    RANSWallAlgebra::Vector Refine(const Matrix &matrix,const RANSWallAlgebra::Vector &rhs,bool transpose=false) const
    {
        // Keep the original pivoted LU as a correction solver. Recompute
        // the defect and accumulate both corrections in extended precision.
        const auto initial=Apply(rhs,transpose);
        std::vector<long double> solution(initial.begin(),initial.end()),defect(m_n);
        for(int pass=0;pass<2;++pass)
        {
            std::copy(rhs.begin(),rhs.end(),defect.begin());
            if(transpose)
                for(int i=0;i<m_n;++i)for(int j=0;j<m_n;++j)
                    defect[j]-=(long double)matrix[i][j]*solution[i];
            else
                for(int i=0;i<m_n;++i)for(int j=0;j<m_n;++j)
                    defect[i]-=(long double)matrix[i][j]*solution[j];
            const auto correction=Apply(RANSWallAlgebra::Vector(defect.begin(),defect.end()),transpose);
            for(int i=0;i<m_n;++i)solution[i]+=correction[i];
        }
        return RANSWallAlgebra::Vector(solution.begin(),solution.end());
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
    if(m_nonlinearPolicy=="ScalarFeasibleR3")
    {
        Require(!scalarOnly,"VerifiedR3Action is a coupled action, not the frozen R2 scalar diagnostic");
        const auto domainStarted=std::chrono::steady_clock::now();
        const auto checked=ScalarDomain(s,time,data);
        const double domainSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-domainStarted).count();
        Require(checked.admissible,"R3 Jacobian base is outside its scalar evaluation domain");
        const auto constraintsStarted=std::chrono::steady_clock::now();
        const auto constraints=Constraints(s,time,data);
        const double constraintsSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-constraintsStarted).count();
        auto jac=NonpressureJacobianR3(s,time,data,metadata);
        Require(jac.size()==Size(),"R3 spatial action matrix has wrong row count");
        for(const auto &row:jac)Require(row.size()==Size(),"R3 spatial action matrix has wrong column count");
        const auto pressureStarted=std::chrono::steady_clock::now();
        for(int column=offset[InnerPressure()];column<Size();++column)
        {
            Vector direction(Size());direction[column]=1.;
            const auto action=PressureDirectionalAction(s,direction,time,data);
            for(int row=0;row<Size();++row)jac[row][column]=action[row];
        }
        const double pressureSeconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-pressureStarted).count();
        const auto massStarted=std::chrono::steady_clock::now();
        for(int block=0;block<InnerPressure();++block)for(int local=0;local<size[block];++local)
        {
            const int column=offset[block]+local;
            const bool external=constraints.count(column) &&
                !(block==Epsilon() && m_transitionEnabled && m_traceIds.count(local));
            if(external)continue;
            for(int row=0;row<size[block];++row)if(!constraints.count(offset[block]+row))
                jac[offset[block]+row][column]+=alpha*m_mass[block][row][local];
        }
        if(metadata)metadata->push_back(J::Object({{"component",J::Quote("complete_verified_R3_action")},
            {"pressure_and_gauge_columns",J::Number(Size()-offset[InnerPressure()])},
            {"wrapper_scalar_domain_calls","1"},{"wrapper_scalar_domain_seconds",J::Number(domainSeconds)},
            {"wrapper_constraints_calls","1"},{"wrapper_constraints_seconds",J::Number(constraintsSeconds)},
            {"pressure_action_seconds",J::Number(pressureSeconds)},
            {"exact_time_mass_seconds",J::Number(std::chrono::duration<double>(std::chrono::steady_clock::now()-massStarted).count())},
            {"exact_time_mass","true"},{"whole_residual_difference_calls","0"}}));
        return jac;
    }
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

std::vector<RANSWallEvolutionPatch::WallExtremum>
RANSWallEvolutionPatch::WallTraceR4(const Vector &state) const
{
    Require(m_dimension==2 && !Mapped(),"R4 wall trace requires canonical flat 2D");
    const auto &field=*fields[K()];
    const auto native=field.Field();
    const auto map=native->GetLocalToGlobalMap();
    std::vector<WallExtremum> extrema;
    for(int e=0;e<field.Elements().size();++e)
    {
        const auto &el=field.Elements()[e];
        if(el.lower[2]!=0.)continue;
        Require(el.degree==3 && native->GetExp(e)->GetBasis(0)->GetBasisType()==LibUtilities::eModified_A &&
            native->GetExp(e)->GetBasis(1)->GetBasisType()==LibUtilities::eModified_A,
            "R4 wall trace requires the frozen tensor cubic basis");
        const int start=native->GetCoeff_Offset(e);
        std::array<long double,4> c{};
        for(int j=0;j<4;++j)
            c[j]=map->GetLocalToGlobalSign(start+j)*state[offset[K()]+map->GetLocalToGlobalMap(start+j)];
        const long double a0=(c[0]+c[1])/2+c[2]/4, a1=(c[1]-c[0]+c[3])/2;
        const long double a2=-c[2]/4, a3=-c[3]/2;
        const long double A=3*a3,B=2*a2,C=a1;
        std::vector<long double> points{-1,1};
        auto add=[&](long double x){if(std::abs(x)<1)points.push_back(x);};
        if(A==0)
        {
            if(B!=0)add(-C/B);
        }
        else
        {
            const long double discriminant=B*B-4*A*C;
            if(discriminant>=0)
            {
                const long double q=-(B+std::copysign(std::sqrt(discriminant),B))/2;
                if(q==0)add(-B/(2*A));
                else {add(q/A);add(C/q);}
            }
        }
        const long double error=128*std::numeric_limits<long double>::epsilon()*
            (std::abs(c[0])+std::abs(c[1])+std::abs(c[2])+std::abs(c[3]));
        for(long double x:points)
        {
            const double physical=el.lower[0]+(x+1)*(el.upper[0]-el.lower[0])/2;
            const long double value=x==-1?c[0]:x==1?c[1]:a0+x*(a1+x*(a2+x*a3));
            extrema.push_back({field.ElementPoint(e,{physical,0,0}),value,
                x==-1 || x==1?0.L:error});
        }
    }
    Require(!extrema.empty(),"R4 wall trace has no ground elements");
    return extrema;
}

std::string RANSWallEvolutionPatch::WallTraceJSON(const Vector &state,bool &admissible) const
{
    const auto points=WallTraceR4(state);
    const auto worst=std::min_element(points.begin(),points.end(),[](const auto &a,const auto &b)
        {return a.value-a.error<b.value-b.error;});
    admissible=worst->value>=worst->error;
    return J::Object({{"admissible",J::Boolean(admissible)},
        {"minimum",J::Number(worst->value)},{"arithmetic_scale",J::Number(worst->error)},
        {"point",J::Numbers(Vector(worst->row.coordinate.begin(),worst->row.coordinate.end()))},
        {"geometry",J::Number(worst->row.geometry)},{"extrema_checked",J::Number(points.size())}});
}

bool RANSWallEvolutionPatch::RestrictedInner(const ScalarObservation &p) const
{
    const double r=Distance(p.coordinate),H=m_mapping.At(p.coordinate).height;
    return r<H-1e-11 || (Same(r,H) && p.side<0);
}

RANSWallGlobalPoint RANSWallEvolutionPatch::RestrictedGradientRow(const ScalarObservation &p) const
{
    // Keep the one-sided native element gradient even at a shared endpoint.
    for(int e=0;e<fields[0]->Elements().size();++e)
        if(fields[0]->Elements()[e].geometry==p.geometry)
            return Transform(fields[0]->ElementPoint(e,p.coordinate));
    return Point(0,p.coordinate,p.side);
}

double RANSWallEvolutionPatch::OuterBound::Action(const Vector &step,int offset) const
{
    return (axis<0?row.Value(step,offset):row.Gradient(step,offset)[axis])/(scale*width);
}

double RANSWallEvolutionPatch::OuterViolation::Max() const
{
    return *std::max_element(family.begin(),family.end());
}

std::vector<RANSWallEvolutionPatch::OuterBound> RANSWallEvolutionPatch::OuterBoundsR5(
    const ScalarDomainResult &domain) const
{
    Require(domain.admissible,"R5 box requires the complete scalar domain");
    const double us=m_robustUstar,K0=m_parameters.model.Q()*us*us;
    const auto &state=domain.effectiveState;
    std::vector<OuterBound> bounds;
    // Identical audit rows can be reached from several inventories. Keep
    // their inventory identities in the reported aggregate.
    for(const auto &p:domain.observations)
    {
        if(p.inventory=="physical_sampled_K")
            bounds.push_back({Point(K(),p.coordinate,p.side),3,K(),-1,p.side,p.geometry,false,
                K0,p.value/K0,1.,.01,.99,1.01,p.inventory});
        if((p.block!=K() && !p.dependentTrace) || RestrictedInner(p))continue;
        // Normalization quadrature is open and lies strictly inside H.
        // Fail explicitly if a moving point ever enters this frozen box
        // inventory. Moving rays still enter projected trace through J.
        Require(p.inventory.find("normalization_ray")==std::string::npos,
            "R5 outer inventory has a moving normalization ray");
        const double r=Distance(p.coordinate),g0=us/(m_parameters.model.Kappa()*r);
        const auto velocity=RestrictedGradientRow(p);
        for(int a=0;a<m_dimension;++a)
        {
            const auto g=velocity.Gradient(state,offset[a]);
            for(int b=0;b<m_dimension;++b)
            {
                const int axis=PhysicalAxis(b);
                const double reference=m_referenceFlow[PhysicalAxis(a)]*m_referenceNormal[axis];
                bounds.push_back({velocity,0,a,axis,p.side,p.geometry,false,
                    g0,g[axis]/g0-reference,reference,.001,-.001,.001,p.inventory});
            }
        }
        const auto k=Point(K(),p.coordinate,p.side);
        bounds.push_back({k,1,K(),-1,p.side,p.geometry,false,K0,k.Value(state,offset[K()])/K0,
            1.,.01,.99,1.01,p.inventory});
        const double e0=us*us*us/(m_parameters.model.Kappa()*r);
        const int side=Same(r,m_mapping.At(p.coordinate).height)?1:0;
        const auto epsilon=Point(Epsilon(),p.coordinate,side);
        const double value=p.dependentTrace?p.value:epsilon.Value(state,offset[Epsilon()]);
        bounds.push_back({epsilon,2,Epsilon(),-1,p.side,p.geometry,p.dependentTrace,e0,value/e0,
            1.,.02,.98,1.02,p.inventory});
    }
    return bounds;
}

RANSWallEvolutionPatch::OuterViolation RANSWallEvolutionPatch::OuterViolationR5(
    const std::vector<OuterBound> &bounds,const Vector &effective,const Vector &trace,double fraction) const
{
    OuterViolation result;
    std::array<const OuterBound *,5> witness{};
    std::array<double,5> excess,values,signs,actions;
    excess.fill(-std::numeric_limits<double>::infinity());
    for(const auto &b:bounds)
    {
        const double action=effective.empty()?0.:b.Action(b.projected?trace:effective,offset[b.block]);
        const double change=fraction*action;
        const double value=b.value+b.width*change;
        for(double sign:{-1.,1.})
        {
            const double violation=sign*((b.value-(sign<0?b.lower:b.upper))/b.width+change);
            result.family[b.family]=std::max(result.family[b.family],violation);
            result.sum+=std::max(0.,violation);
            const int index=b.projected?4:b.family;
            if(violation>excess[index])
            {
                excess[index]=violation;values[index]=value;signs[index]=sign;actions[index]=action;witness[index]=&b;
            }
        }
    }
    std::vector<std::string> records;
    for(int i=0;i<witness.size();++i)if(const auto *b=witness[i])
        records.push_back(J::Object({{"family",J::Number(b->family)},{"inventory",J::Quote(b->inventory)},
            {"point",J::Numbers(b->row.coordinate)},{"element_geometry",J::Number(b->row.geometry)},
            {"observation_geometry",J::Number(b->observationGeometry)},{"side",J::Number(b->side)},
            {"component",J::Numbers(Vector{double(b->axis<0?-1:PhysicalAxis(b->block)),double(b->axis)})},
            {"projected_trace",J::Boolean(b->projected)},{"sign",J::Number(signs[i])},
            {"reference_scale",J::Number(b->scale)},{"reference_component",J::Number(b->reference)},
            {"signed_excess",J::Number(excess[i])},{"value",J::Number(values[i])},
            {"directional_quantity_change",J::Number(actions[i]*b->width)}}));
    result.report=J::Object({{"maximum",J::Number(result.Max())},{"family_maxima",J::Numbers(result.family)},
        {"sum_positive",J::Number(result.sum)},{"quantity_count",J::Number(bounds.size())},
        {"witnesses",J::List(records)}});
    return result;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::CoupledDirectionR4(
    const ScalarDomainResult &domain,const Vector &newton,const WallMatrix &jac,
    const std::map<int,double> &constraints,const std::function<Vector(Vector,bool)> &inverse,
    std::string &details,const std::vector<OuterBound> &outer,double restoration,bool inward,
    RANSWallR3::DurableProgress *diagnostic,const std::string &tag,const Vector &residual) const
{
    // Parameterize d=d_N+J^{-1} B z. B injects the scaled Cholesky
    // factors on free differential rows, leaving every algebraic row zero.
    // The physical least-squares objective becomes ||z||^2/2.
    struct MetricBlock {std::vector<int> free;WallMatrix lower;double scale;};
    const double us=m_robustUstar,H=m_parameters.geometry.height;
    const double K0=m_parameters.model.Q()*us*us,E0=us*us*us/(m_parameters.model.Kappa()*H);
    const std::array<double,7> scales{us,us,K0,E0,us*us,us*us,us/H};
    std::vector<MetricBlock> metric;
    for(int b=0;b<InnerPressure();++b)
    {
        MetricBlock block;
        block.scale=scales[b]*us/H*std::sqrt(m_volumes[b]);
        for(int i=0;i<size[b];++i)if(!constraints.count(offset[b]+i))block.free.push_back(offset[b]+i);
        const int n=block.free.size();
        block.lower.assign(n,Vector(n));
        for(int i=0;i<n;++i)for(int j=0;j<=i;++j)
        {
            long double value=m_mass[b][block.free[i]-offset[b]][block.free[j]-offset[b]];
            for(int k=0;k<j;++k)value-=block.lower[i][k]*block.lower[j][k];
            Require(i!=j || value>0,"R4 free consistent mass is not SPD");
            block.lower[i][j]=i==j?std::sqrt(value):value/block.lower[j][j];
        }
        metric.push_back(std::move(block));
    }
    auto dot=[](const Vector &a,const Vector &b)
    {
        long double result=0;
        for(std::size_t i=0;i<a.size();++i)result+=(long double)a[i]*b[i];
        return static_cast<double>(result);
    };
    auto axpy=[](Vector &a,double factor,const Vector &b)
    {for(std::size_t i=0;i<a.size();++i)a[i]+=factor*b[i];};
    const bool refined=RefinedDirectionR6();
    auto metricAction=[&](const Vector &x,bool transpose)
    {
        if(refined)
        {
            std::vector<long double> result(Size());
            for(const auto &block:metric)
                for(std::size_t i=0;i<block.free.size();++i)
                    for(std::size_t j=0;j<=i;++j)
                    {
                        const long double entry=(long double)block.scale*block.lower[i][j];
                        if(transpose)result[block.free[j]]+=entry*x[block.free[i]];
                        else result[block.free[i]]+=entry*x[block.free[j]];
                    }
            return Vector(result.begin(),result.end());
        }
        Vector result(Size());
        for(const auto &block:metric)
            for(std::size_t i=0;i<block.free.size();++i)
                for(std::size_t j=0;j<=i;++j)
                    if(transpose)result[block.free[j]]+=block.scale*block.lower[i][j]*x[block.free[i]];
                    else result[block.free[i]]+=block.scale*block.lower[i][j]*x[block.free[j]];
        return result;
    };
    auto external=[&](int col)
    {
        return constraints.count(col) && !(col>=offset[Epsilon()] && col<offset[InnerPressure()] &&
            m_traceIds.count(col-offset[Epsilon()]));
    };
    auto scalarRow=[&](const RANSWallGlobalPoint &point,int block,double scale,bool trace,int axis=-1)
    {
        const auto &basis=axis<0?point.phi:point.derivative[axis];
        const double constantAction=axis<0?1.:0.;
        Vector a(Size());long double constant=0;
        for(std::size_t j=0;j<point.ids.size();++j)constant+=(long double)basis[j]*point.constant[j];
        for(std::size_t j=0;j<point.ids.size();++j)
        {
            const int col=offset[block]+point.ids[j];
            if(external(col))continue;
            const double phi=(basis[j]+(j==0?(constantAction-constant)/point.constant[0]:0))/scale;
            a[col]+=phi;
            if(trace && m_traceIds.count(point.ids[j]))axpy(a,-phi,jac[col]);
        }
        return a;
    };

    constexpr double radius=2., reserve=.05, near=.005, feasibility=5e-12;
    struct Witness
    {
        int block=-1,axis=-1,side=0,coefficient=-1;
        double scale=1;
        bool trace=false;
        RANSWallGlobalPoint point;
    };
    struct Cut
    {
        Vector normal;
        double rhs,length;
        std::string kind;
        RANSWallCoordinate point;
        int id;
    };
    struct OriginalCut
    {
        Witness witness;
        Vector row,transposeSolution;
        double rhs;
        std::string kind;
        RANSWallCoordinate point;
    };
    std::vector<Cut> active;
    std::vector<OriginalCut> original;
    std::vector<std::string> exchanges;
    auto ground=WallTraceR4(domain.effectiveState);
    // Inward motion is allowed at the observed native witness and all
    // currently near-active extrema, including tied and endpoint minima.
    const RANSWallCoordinate witness{1.6046496089512394,0,0};
    const auto witnessRow=Point(K(),witness);
    const double witnessValue=witnessRow.Value(domain.effectiveState,offset[K()]);
    if(witnessValue/K0<=near)ground.push_back({witnessRow,witnessValue,0});

    Vector direction=newton,z(Size()),multipliers;
    int releases=0,dependentExchanges=0;
    auto qr=[&]()
    {
        // Twice reorthogonalized QR of the active residual-space normals.
        WallMatrix Q,R(active.size(),Vector(active.size()));
        for(std::size_t j=0;j<active.size();++j)
        {
            auto v=active[j].normal;
            for(int pass=0;pass<2;++pass)for(std::size_t i=0;i<j;++i)
            {
                const double coefficient=dot(Q[i],v);
                R[i][j]+=coefficient;
                axpy(v,-coefficient,Q[i]);
            }
            R[j][j]=std::sqrt(dot(v,v));
            Require(R[j][j]>1e-10,"R4 active cuts are numerically rank deficient");
            for(double &entry:v)entry/=R[j][j];
            Q.push_back(std::move(v));
        }
        return std::pair{std::move(Q),std::move(R)};
    };
    auto triangular=[](const WallMatrix &R,Vector value)
    {
        for(int i=static_cast<int>(value.size())-1;i>=0;--i)
        {
            for(std::size_t j=i+1;j<value.size();++j)value[i]-=R[i][j]*value[j];
            value[i]/=R[i][i];
        }
        return value;
    };
    std::string lastCut;
    std::array<double,4> ceiling{};
    if(!outer.empty())
    {
        ceiling=OuterViolationR5(outer).family;
        for(double &value:ceiling)value*=1.-restoration;
    }
    const int maxExchange=outer.empty()?40:96;
    int exchange=0;
    bool published=false;
    double worstNative=0,worstExtended=0,maximumScanDifference=0;
    auto record=[&](const Cut *entering,const WallMatrix &R,const Vector &effective,const Vector &trace,
                    const std::string &status)
    {
        if(!diagnostic)return;
        std::vector<std::string> rows;
        Vector diagonal;
        for(int i=0;i<R.size();++i)diagonal.push_back(R[i][i]);
        auto rowRecord=[&](const Cut &cut,bool incoming,double multiplier)
        {
            const auto &o=original[cut.id];
            const auto &w=o.witness;
            const double native=w.coefficient>=0?direction[w.coefficient]/w.scale:
                (w.axis<0?w.point.Value(w.trace?trace:effective,offset[w.block]):
                    w.point.Gradient(effective,offset[w.block])[w.axis])/w.scale;
            const std::string name=tag+".cut"+std::to_string(rows.size());
            auto unscaled=cut.normal;
            for(double &x:unscaled)x*=cut.length;
            rows.push_back(J::Object({{"id",J::Number(cut.id)},{"kind",J::Quote(cut.kind)},
                {"point",J::Numbers(cut.point)},{"geometry",J::Number(w.point.geometry)},
                {"block",J::Number(w.block)},{"axis",J::Number(w.axis)},{"side",J::Number(w.side)},
                {"scale",J::Number(w.scale)},{"projected",J::Boolean(w.trace)},
                {"coefficient",J::Number(w.coefficient)},{"incoming",J::Boolean(incoming)},
                {"rhs",J::Number(o.rhs)},{"transformed_rhs",J::Number(cut.rhs*cut.length)},
                {"normal_length",J::Number(cut.length)},{"multiplier_scaled",J::Number(multiplier)},
                {"native_gap",J::Number(native-o.rhs)},{"coefficient_gap",J::Number(dot(o.row,direction)-o.rhs)},
                {"transformed_gap",J::Number(dot(unscaled,z)-cut.rhs*cut.length)},
                {"a",diagnostic->VectorPayload(name+".row",o.row)},
                {"n",diagnostic->VectorPayload(name+".normal",unscaled)},
                {"transpose_solution",diagnostic->VectorPayload(name+".transpose",o.transposeSolution)}}));
        };
        for(int i=0;i<active.size();++i)rowRecord(active[i],false,multipliers[i]);
        if(entering)rowRecord(*entering,true,0.);
        diagnostic->Event("r6_cut_audit",{{"name",J::Quote(tag)},{"status",J::Quote(status)},
            {"exchange",J::Number(exchange+1)},{"rows",J::List(rows)},{"qr_diagonal",J::Numbers(diagonal)},
            {"original_tolerance",J::Number(feasibility)},{"worst_native_violation",J::Number(worstNative)},
            {"worst_extended_violation",J::Number(worstExtended)},
            {"maximum_scan_difference",J::Number(maximumScanDifference)},
            {"exchanges",J::List(exchanges)},{"z",diagnostic->VectorPayload(tag+".z",z)},
            {"direction",diagnostic->VectorPayload(tag+".direction",direction)},
            {"newton",diagnostic->VectorPayload(tag+".newton",newton)},
            {"Bz",diagnostic->VectorPayload(tag+".Bz",metricAction(z,false))}});
        published=true;
    };
    for(;exchange<maxExchange;++exchange)
    {
        auto [Q,R]=qr();
        Vector qCoordinates(active.size());
        for(std::size_t j=0;j<active.size();++j)
        {
            double rhs=active[j].rhs;
            for(std::size_t i=0;i<j;++i)rhs-=R[i][j]*qCoordinates[i];
            qCoordinates[j]=rhs/R[j][j];
        }
        multipliers=triangular(R,qCoordinates);
        Require(multipliers.empty() || *std::min_element(multipliers.begin(),multipliers.end())>=-1e-8,
            "R4 active solution lost dual feasibility");
        z.assign(Size(),0.);
        for(std::size_t i=0;i<Q.size();++i)axpy(z,qCoordinates[i],Q[i]);
        direction=newton;
        if(!active.empty())
        {
            auto load=metricAction(z,false);
            if(refined)
            {
                Require(residual.size()==load.size(),"R6 direction requires its actual coupled residual");
                axpy(load,-1,residual);
                direction=inverse(std::move(load),false);
            }
            else axpy(direction,1,inverse(std::move(load),false));
        }
        auto effective=direction;
        const auto action=Multiply(jac,direction);
        for(const auto &[col,value]:constraints)if(external(col))effective[col]=0.;
        auto trace=effective;
        for(int id:m_traceIds)trace[offset[Epsilon()]+id]-=action[offset[Epsilon()]+id];
        std::vector<long double> extendedTrace;
        if(refined)
        {
            extendedTrace.assign(effective.begin(),effective.end());
            for(int id:m_traceIds)
            {
                const int i=offset[Epsilon()]+id;
                for(int j=0;j<Size();++j)extendedTrace[i]-=(long double)jac[i][j]*direction[j];
            }
        }

        Vector violated;
        worstNative=worstExtended=maximumScanDifference=0;
        double violation=feasibility,rhs=0;
        std::string kind;
        RANSWallCoordinate location{};
        Witness selected;
        auto inspect=[&](double lower,double observed,const std::string &label,
                         const RANSWallCoordinate &point,auto makeRow,int block,int axis,double scale,
                         bool projected,int side,const RANSWallGlobalPoint *row,int coefficient=-1)
        {
            worstNative=std::max(worstNative,lower-observed);
            if(refined)
            {
                long double value=0;
                if(coefficient>=0)value=(long double)direction[coefficient]/scale;
                else
                {
                    const auto &basis=axis<0?row->phi:row->derivative[axis];
                    auto entry=[&](int i)->long double
                    {
                        const int j=offset[block]+row->ids[i];
                        return projected?extendedTrace[j]:effective[j];
                    };
                    const long double centre=entry(0)/row->constant[0];
                    if(axis<0)value=centre;
                    for(int i=0;i<basis.size();++i)value+=(long double)basis[i]*(entry(i)-centre*row->constant[i]);
                    value/=scale;
                }
                worstExtended=std::max(worstExtended,static_cast<double>(lower-value));
                maximumScanDifference=std::max(maximumScanDifference,static_cast<double>(std::abs(value-observed)));
                // Both the native scan and the separate extended row
                // evaluation must satisfy the unchanged tolerance.
                observed=std::min(observed,static_cast<double>(value));
            }
            if(lower-observed<=violation)return;
            violation=lower-observed;
            rhs=lower;
            kind=label;
            location=point;
            violated=makeRow();
            if(diagnostic)
                selected={block,axis,side,coefficient,scale,projected,row?*row:RANSWallGlobalPoint{}};
        };
        for(int b=0;b<=Gauge();++b)for(int i=offset[b];i<offset[b]+size[b];++i)
            for(double sign:{-1.,1.})
                inspect(-radius,sign*direction[i]/scales[b],"dimensionless_box",{},[&]
                {Vector row(Size());row[i]=sign/scales[b];return row;},b,-1,scales[b]/sign,false,0,nullptr,i);
        for(const auto &observation:domain.observations)
        {
            const auto point=Point(observation.block,observation.coordinate,observation.side);
            const double delta=point.Value(observation.dependentTrace?trace:effective,offset[observation.block]);
            inspect(-(1-reserve)*observation.value/observation.scale,delta/observation.scale,
                observation.inventory,point.coordinate,[&]
                {return scalarRow(point,observation.block,observation.scale,observation.dependentTrace);},
                observation.block,-1,observation.scale,observation.dependentTrace,observation.side,&point);
        }
        for(const auto &p:ground)if(p.value/K0<=near && (inward || p.value<=0))
            inspect(0,p.row.Value(effective,offset[K()])/K0,"inward_ground_row",p.row.coordinate,[&]
                {return scalarRow(p.row,K(),K0,false);},K(),-1,K0,false,0,&p.row);
        auto wall=effective;
        axpy(wall,1-reserve,domain.effectiveState);
        for(const auto &p:WallTraceR4(wall))
            inspect(-(1-reserve)*p.row.Value(domain.effectiveState,offset[K()])/K0,
                p.row.Value(effective,offset[K()])/K0,"continuous_wall_reserve",p.row.coordinate,[&]
                {return scalarRow(p.row,K(),K0,false);},K(),-1,K0,false,0,&p.row);
        for(const auto &bound:outer)
        {
            const double change=bound.Action(bound.projected?trace:effective,offset[bound.block]);
            for(double sign:{-1.,1.})
            {
                const double limit=sign<0?bound.lower:bound.upper;
                const double excess=sign*(bound.value-limit)/bound.width;
                inspect(excess-ceiling[bound.family],-sign*change,
                    "outer_"+bound.inventory+"_"+std::to_string(bound.family),bound.row.coordinate,[&]
                    {return scalarRow(bound.row,bound.block,-sign*bound.scale*bound.width,bound.projected,bound.axis);},
                    bound.block,bound.axis,-sign*bound.scale*bound.width,bound.projected,bound.side,&bound.row);
            }
        }
        if(violated.empty())
        {
            Vector stationarity=z;
            std::vector<std::string> activeRows;
            double activeError=0,complementarity=0;
            for(std::size_t i=0;i<active.size();++i)
            {
                const auto &row=active[i];
                axpy(stationarity,-multipliers[i],row.normal);
                const double gap=dot(row.normal,z)-row.rhs;
                activeError=std::max(activeError,std::abs(gap));
                complementarity=std::max(complementarity,std::abs(gap*multipliers[i]));
                activeRows.push_back(J::Object({{"inventory",J::Quote(row.kind)},
                    {"point",J::Numbers(Vector(row.point.begin(),row.point.end()))},
                    {"multiplier",J::Number(multipliers[i]/row.length)}}));
            }
            const double stationarityNorm=std::sqrt(dot(stationarity,stationarity));
            Require(stationarityNorm<=1e-9 && activeError<=1e-9 && complementarity<=1e-9,
                "R4 constrained solve failed its KKT check");
            details=J::Object({{"method",J::Quote(active.empty()?"ordinary_newton":"residual_space_active_QR")},
                {"active",J::List(activeRows)},{"releases",J::Number(releases)},
                {"dependent_exchanges",J::Number(dependentExchanges)},
                {"exchange_passes",J::Number(exchange+1)},{"trust_radius","2"},
                {"model_feasibility_tolerance",J::Number(feasibility)},{"stationarity_l2",J::Number(stationarityNorm)},
                {"active_equality_error",J::Number(activeError)},{"complementarity",J::Number(complementarity)},
                {"scaled_residual_objective",J::Number(.5*dot(z,z))},
                {"witness_k",J::Number(witnessValue)},
                {"witness_newton_delta",J::Number(witnessRow.Value(newton,offset[K()]))},
                {"witness_direction_delta",J::Number(witnessRow.Value(effective,offset[K()]))},
                {"outer_restoration_fraction",outer.empty()?"null":J::Number(restoration)},
                {"outer_linear_ceiling",J::Numbers(ceiling)},
                {"inward_positive_ground_rule",J::Boolean(inward)}});
            record(nullptr,R,effective,trace,"passed_original_cut_scan");
            return direction;
        }
        auto transposeSolution=inverse(violated,true);
        auto normal=metricAction(transposeSolution,true);
        const double length=std::sqrt(dot(normal,normal));
        Require(length>1e-14,"R4 feasibility cut has no free differential response");
        for(double &value:normal)value/=length;
        const int id=original.size();
        const Cut incoming{std::move(normal),(rhs-dot(violated,newton))/length,length,kind,location,id};
        if(diagnostic)
        {
            Vector aliases,activeIds;
            for(int i=0;i<original.size();++i)if(original[i].row==violated && original[i].rhs==rhs)aliases.push_back(i);
            for(const auto &cut:active)activeIds.push_back(cut.id);
            original.push_back({std::move(selected),violated,std::move(transposeSolution),rhs,kind,location});
            exchanges.push_back(J::Object({{"op",J::Quote("enter")},{"id",J::Number(id)},
                {"kind",J::Quote(kind)},{"point",J::Numbers(location)},{"active",J::Numbers(activeIds)},
                {"exact_previous_alias_ids",J::Numbers(aliases)},{"violation",J::Number(violation)}}));
            if(exchange+1==maxExchange)record(&incoming,R,effective,trace,"exchange_limit_candidate");
            if(refined && std::any_of(active.begin(),active.end(),[&](const Cut &cut)
                {return original[cut.id].row==violated && original[cut.id].rhs==rhs;}))
            {
                if(!published)record(&incoming,R,effective,trace,"active_original_row_accuracy_unresolved");
                throw std::runtime_error("R6 active original row still violates the verified cut scan");
            }
        }
        lastCut=kind+" violation "+J::Number(violation)+" at "+J::Numbers(Vector(location.begin(),location.end()));
        // Dual active-set projection. Keep the incoming constraint through
        // every blocking release. Re-selecting a different cut mid-update
        // can cycle when several evaluation rows share the same polynomial.
        for(;;)
        {
            auto [basis,upper]=qr();
            auto primal=incoming.normal;
            Vector projected(basis.size());
            for(int pass=0;pass<2;++pass)for(std::size_t i=0;i<basis.size();++i)
            {
                const double coefficient=dot(basis[i],primal);
                projected[i]+=coefficient;
                axpy(primal,-coefficient,basis[i]);
            }
            const auto dual=triangular(upper,std::move(projected));
            const double square=dot(primal,primal);
            const double distance=square>1e-20?
                std::max(0.,incoming.rhs-dot(incoming.normal,z))/square:
                std::numeric_limits<double>::infinity();
            double blocker=std::numeric_limits<double>::infinity();
            std::size_t leaving=basis.size();
            for(std::size_t i=0;i<basis.size();++i)if(dual[i]>1e-12)
            {
                const double candidate=std::max(0.,multipliers[i])/dual[i];
                if(candidate<blocker){blocker=candidate;leaving=i;}
            }
            if(std::isfinite(distance) && distance<=blocker)
            {
                if(diagnostic)exchanges.push_back(J::Object({{"op",J::Quote("add")},{"id",J::Number(id)},
                    {"projected_length",J::Number(std::sqrt(square))}}));
                active.push_back(incoming);
                break;
            }
            if(leaving>=basis.size() && !published)record(&incoming,upper,effective,trace,"dual_exchange_unresolved");
            Require(leaving<basis.size(),"R4 cuts have no feasible dual exchange: "+lastCut);
            if(diagnostic)exchanges.push_back(J::Object({{"op",J::Quote("release")},
                {"id",J::Number(active[leaving].id)},{"incoming",J::Number(id)},
                {"projected_length",J::Number(std::sqrt(square))},{"dual_blocker",J::Number(blocker)}}));
            axpy(z,blocker,primal);
            axpy(multipliers,-blocker,dual);
            active.erase(active.begin()+leaving);
            multipliers.erase(multipliers.begin()+leaving);
            ++releases;
            if(square<=1e-20)++dependentExchanges;
        }
    }
    throw std::runtime_error("R4 bounded active cut exchange exhausted: "+lastCut);
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::SolveScalarFeasibleR3(
    Vector s,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,
    int &iterations,double tolerance,int maximumIterations,bool benchmarkOnly,bool localR4,bool localR5,
    const WallMatrix &frozenMatrix,const std::string &frozenAssembly) const
{
    Require(m_nonlinearPolicy=="ScalarFeasibleR3" && m_dimension==2 && !Mapped() &&
        m_transitionEnabled && m_cappedDestruction && m_comm->GetSize()==1,
        "ScalarFeasibleR3 requires its opt-in serial flat 2D Upstream+Capped policy");
    Require(Size()<=3000 && alpha==400. && time==.0025 && tolerance==2e-12 &&
        maximumIterations>0 && maximumIterations<=(JointRestorationR5()?32:16),"R3 requires the unchanged bounded exact target");
    Require(!localR4 || (CoupledR4() && maximumIterations==1),"R4 local experiment uses one bounded update");
    Require(!localR5 || (JointRestorationR5() && maximumIterations==1),
        "R5 local experiment requires the explicit bounded restoration policy");
    Require(frozenMatrix.empty() || ((DirectionReplayR6() || RefinedDirectionR6()) && localR5),"Frozen J is restricted to local R6");
    iterations=0;
    using Clock=std::chrono::steady_clock;
    std::map<std::string,std::pair<int,double>> profile;
    auto timed=[&](const std::string &name,auto function)
    {
        const auto start=Clock::now();++profile[name].first;
        auto value=function();profile[name].second+=std::chrono::duration<double>(Clock::now()-start).count();
        return value;
    };
    auto profiling=[&]()
    {
        std::vector<std::pair<std::string,std::string>> rows;
        for(const auto &[name,value]:profile)rows.push_back({name,J::Object({
            {"calls",J::Number(value.first)},{"seconds",J::Number(value.second)}})});
        rows.push_back({"scope",J::Quote("outer solver calls; complete action assembly includes separately reported wrapper/backend work; line-search time includes its domain/residual calls")});
        return J::Object(rows);
    };
    const int solve=++m_feasibleSolveCount;
    RANSWallR3::DurableProgress durable(m_r3ProgressFile,solve,alpha,time,history,ActionPolicyName()+":v1");
    durable.Event("protocol",{{"inner_stress_model",J::Quote(InnerStressName())},
        {"tolerance",J::Number(tolerance)},{"maximum_Newton_updates",J::Number(maximumIterations)},
        {"benchmark_only",J::Boolean(benchmarkOnly)},{"globalization",J::Quote(m_globalization)},
        {"direction_accuracy",J::Quote(m_directionAccuracy)},
        {"local_direction_experiment",J::Boolean(localR4 || localR5)},{"no_auxiliary_alpha_sequence","true"}});
    const double us=m_robustUstar,K0=m_parameters.model.Q()*us*us;
    const double E0=us*us*us/(m_parameters.model.Kappa()*m_parameters.geometry.height);
    auto domain=timed("scalar_domain",[&]{return ScalarDomain(s,time,data,true);});
    auto constraints=std::map<int,double>{};Vector residual,norms;std::string terminal="started";
    bool alreadySaved=false;int stagnant=0;
    auto domainReport=[&](const ScalarDomainResult &value)
    {
        const auto json=ScalarDomainJSON(value);
        if(!GroundSafeguard())return json;
        bool wallOK=false;
        const auto trace=WallTraceJSON(value.effectiveState,wallOK);
        return json.substr(0,json.size()-1)+",\"r4_wall_trace\":"+trace+"}";
    };
    auto checkWall=[&](const ScalarDomainResult &value)
    {
        if(!GroundSafeguard())return true;
        bool admissible=false;
        WallTraceJSON(value.effectiveState,admissible);
        return admissible;
    };
    auto recordRanges=[&](int iteration,const Vector &state)
    {
        if(!GroundSafeguard())return;
        const auto audit=CompleteRestrictedAudit(state,history,alpha,time,data,true);
        std::vector<std::pair<std::string,std::string>> values;
        for(const auto &[key,row]:audit)values.push_back({key,J::Numbers(row)});
        durable.Event("r4_applicability",{{"iteration",J::Number(iteration)},{"audit",J::Object(values)}});
    };
    auto evaluate=[&]()
    {
        Require(domain.admissible,"R3 iterate is outside its scalar evaluation domain: "+domain.failure);
        constraints=timed("constraints",[&]{return Constraints(s,time,data);});
        residual=timed("residual",[&]{return Residual(s,history,alpha,time,data);});
        norms=timed("mass_dual_norms",[&]{return Norms(residual,constraints,us);});
    };
    try
    {
        evaluate();
        Require(checkWall(domain),"R4 base wall polynomial fails the local trace safeguard");
        for(iterations=0;iterations<=maximumIterations;++iterations)
        {
            // Every accepted iterate is persisted immediately, including one
            // at which stagnation stops the search. No next J starts first.
            if(!alreadySaved)durable.Snapshot(iterations,s,residual,constraints,norms,
                domainReport(domain),profiling());
            if(iterations==0)recordRanges(0,s);
            alreadySaved=false;
            const double error=Maximum(norms);
            if(error<tolerance && !benchmarkOnly)
            {
                terminal="true_coupled_residual_converged";m_lastFeasibleState=s;
                durable.Event("solve_terminal",{{"converged","true"},{"terminal",J::Quote(terminal)},
                    {"iteration",J::Number(iterations)},{"scaled_block_norms",J::Numbers(norms)},{"profiling",profiling()}});
                return s;
            }
            Require(iterations<maximumIterations,JointRestorationR5()?"R5 nonlinear update limit":"R3 Newton iteration limit");
            durable.JacobianStatus(iterations,"incomplete");
            std::vector<std::string> assembly;
            const auto jac=frozenMatrix.empty()?
                timed("complete_action_assembly",[&]{return ScalarFeasibleJacobian(s,history,alpha,time,data,&assembly,false);}):
                frozenMatrix;
            const auto matrix=durable.MatrixPayload(iterations,jac);
            durable.JacobianStatus(iterations,"complete",matrix,J::Object({
                {"shape",J::Numbers(Vector{double(Size()),double(Size())})},
                {"assembly",frozenMatrix.empty()?J::List(assembly):frozenAssembly},
                {"loaded_bound_matrix",J::Boolean(!frozenMatrix.empty())},{"profiling",profiling()}}));
            const auto factor=timed("factorization",[&]{return FeasibleFactor(jac);});
            durable.Event("factorization_completed",{{"iteration",J::Number(iterations)},
                {"matrix",matrix},{"profiling",profiling()}});
            if(benchmarkOnly)
            {
                terminal="one_complete_current_J_and_factorization_benchmark_no_update";m_lastFeasibleState=s;
                durable.Event("solve_terminal",{{"converged","false"},{"benchmark_completed","true"},
                    {"terminal",J::Quote(terminal)},{"iteration",J::Number(iterations)},{"profiling",profiling()}});
                return s;
            }
            const auto inverse=[&](Vector rhs,bool transpose)
            {
                return RefinedDirectionR6()?factor.Refine(jac,rhs,transpose):factor.Apply(std::move(rhs),transpose);
            };
            auto rhs=residual;for(auto &value:rhs)value=-value;
            const auto newton=timed("linear_solve",[&]{return inverse(std::move(rhs),false);});
            auto linear=Multiply(jac,newton);
            for(int row=0;row<Size();++row)linear[row]+=residual[row];
            const auto linearNorms=Norms(linear,constraints,us);
            Require(Maximum(linearNorms)<1e-8*std::max(1.,error),"R3 inaccurate linear solve");
            if(localR4)
            {
                // Compare the ordinary Newton step under the identical wall
                // safeguard. Compute the first ground contact by bisection.
                double hi=1,lo=0;
                auto trial=s;
                for(int pass=0;pass<40;++pass)
                {
                    const double f=(lo+hi)/2;
                    for(int row=0;row<Size();++row)trial[row]=s[row]+f*newton[row];
                    bool admissible=false;
                    WallTraceJSON(ExternalLift(trial,time,data),admissible);
                    if(admissible)lo=f;else hi=f;
                }
                const double fraction=.995*lo;
                for(int row=0;row<Size();++row)trial[row]=s[row]+fraction*newton[row];
                const auto scalar=ScalarDomain(trial,time,data,true);
                Vector values;
                if(scalar.admissible)
                    values=Norms(Residual(trial,history,alpha,time,data),Constraints(trial,time,data),us);
                durable.Event("r4_newton_comparison",{{"fraction",J::Number(fraction)},
                    {"wall_contact_fraction",J::Number(lo)},{"scalar_domain",domainReport(scalar)},
                    {"predicted_norm",J::Number((1-fraction)*error)},{"actual_norms",J::Numbers(values)},
                    {"direction",durable.VectorPayload("comparison.newton",newton)}});
                if(scalar.admissible)recordRanges(-1,trial);
            }
            if(JointRestorationR5())
            {
                const auto bounds=timed("outer_box",[&]{return OuterBoundsR5(domain);});
                const auto baseBox=OuterViolationR5(bounds);
                durable.Event("r5_box",{{"iteration",J::Number(iterations)},{"box",baseBox.report}});
                struct Trial
                {
                    Vector state,residual,norms;
                    ScalarDomainResult domain;
                    std::map<int,double> constraints;
                    OuterViolation box;
                    double fraction=0;
                    std::string stepPin;
                };
                int trialCount=0;
                auto compare=[&](const std::string &name,bool joint,double restoration,bool inward)
                {
                    Trial result;
                    std::string details;
                    Vector step;
                    const auto started=Clock::now();++profile["r5_direction"].first;
                    try
                    {
                        step=CoupledDirectionR4(domain,newton,jac,constraints,inverse,details,
                            joint?bounds:std::vector<OuterBound>{},restoration,inward,
                            (DirectionReplayR6() || RefinedDirectionR6())?&durable:nullptr,
                            "iterate"+std::to_string(iterations)+"."+name,residual);
                    }
                    catch(const std::exception &e)
                    {
                        profile["r5_direction"].second+=std::chrono::duration<double>(Clock::now()-started).count();
                        details=e.what();
                        durable.Event("r5_direction_failure",{{"iteration",J::Number(iterations)},
                            {"name",J::Quote(name)},{"message",J::Quote(details)},{"profiling",profiling()}});
                        return std::pair{std::move(result),details};
                    }
                    profile["r5_direction"].second+=std::chrono::duration<double>(Clock::now()-started).count();
                    const auto action=Multiply(jac,step);
                    auto effective=step;
                    for(const auto &[row,value]:constraints)
                        if(!(row>=offset[Epsilon()] && row<offset[InnerPressure()] &&
                            m_traceIds.count(row-offset[Epsilon()])))effective[row]=0.;
                    auto trace=effective;
                    for(int id:m_traceIds)trace[offset[Epsilon()]+id]-=action[offset[Epsilon()]+id];
                    const auto prefix="iterate"+std::to_string(iterations)+"."+name;
                    result.stepPin=durable.VectorPayload(prefix+".step",step);
                    auto model=[&](double fraction)
                    {
                        auto prediction=residual;
                        for(int i=0;i<Size();++i)prediction[i]+=fraction*action[i];
                        return Norms(prediction,constraints,us);
                    };
                    const auto full=model(1.);
                    Require(*std::max_element(full.begin()+InnerPressure(),full.end())<1e-9,
                        "R5 direction failed full linearized algebraic restoration");
                    double stepNorm=0;
                    const std::array<double,7> scales{us,us,K0,E0,us*us,us*us,us/m_parameters.geometry.height};
                    for(int b=0;b<=Gauge();++b)for(int i=offset[b];i<offset[b]+size[b];++i)
                        stepNorm=std::max(stepNorm,std::abs(step[i])/scales[b]);
                    durable.Event("r5_direction",{{"iteration",J::Number(iterations)},{"name",J::Quote(name)},
                        {"step",result.stepPin},{"action",durable.VectorPayload(prefix+".action",action)},
                        {"dimensionless_step",J::Number(stepNorm)},{"details",details},
                        {"model_block_norms",J::Numbers(full)},
                        {"model_box",OuterViolationR5(bounds,effective,trace,1.).report},
                        {"base_witness_actions",OuterViolationR5(bounds,effective,trace).report},
                        {"profiling",profiling()}});
                    const auto search=Clock::now();++profile["line_search"].first;
                    double fraction=1.;
                    // Reserve half the trial budget for partial restoration
                    // even if the full model constructs a valid direction.
                    const int allowance=RefinedDirectionR6()?
                        std::min(localR5?4:16,(localR5?8:32)-trialCount):
                        (localR5?8:32-trialCount);
                    for(int k=0;k<allowance;++k,fraction*=.5)
                    {
                        ++trialCount;
                        auto trial=s;
                        for(int i=0;i<Size();++i)trial[i]+=fraction*step[i];
                        const auto label=prefix+".trial"+std::to_string(k);
                        auto checked=timed("scalar_domain",[&]{return ScalarDomain(trial,time,data,true);});
                        std::vector<std::pair<std::string,std::string>> event{
                            {"iteration",J::Number(iterations)},{"name",J::Quote(name)},{"backtrack",J::Number(k)},
                            {"fraction",J::Number(fraction)},{"state",durable.VectorPayload(label+".state",trial)},
                            {"scalar_domain",domainReport(checked)},{"model_block_norms",J::Numbers(model(fraction))},
                            {"model_box",OuterViolationR5(bounds,effective,trace,fraction).report},
                            {"profile_branches_changed",J::Boolean(domain.profileBranches!=checked.profileBranches)}};
                        if(!checked.admissible || !checkWall(checked))
                        {
                            event.push_back({"classification",J::Quote(checked.admissible?"ground_trace":"scalar_domain")});
                            durable.Event("r5_trial",event);continue;
                        }
                        try
                        {
                            const auto candidateBox=timed("outer_box",[&]{return OuterViolationR5(OuterBoundsR5(checked));});
                            const auto candidate=timed("residual",[&]{return Residual(trial,history,alpha,time,data);});
                            const auto current=timed("constraints",[&]{return Constraints(trial,time,data);});
                            const auto candidateNorms=Norms(candidate,current,us);
                            const double required=localR5?.01:1e-4*fraction;
                            bool boxDescent=candidateBox.Max()<=(1.-required)*baseBox.Max();
                            for(int family=0;family<4;++family)
                                boxDescent=boxDescent && candidateBox.family[family]<=baseBox.family[family];
                            const bool descent=Maximum(candidateNorms)<(1.-required)*error;
                            const bool accept=descent && (!joint || boxDescent);
                            event.insert(event.end(),{{"classification",J::Quote(accept?"accepted_comparison":
                                !descent?"residual_non_descent":"outer_non_descent_or_family_worsening")},
                                {"box",candidateBox.report},{"scaled_block_norms",J::Numbers(candidateNorms)},
                                {"trace_constraint_norm",J::Number(candidateNorms.back())},
                                {"residual",durable.VectorPayload(label+".residual",candidate)}});
                            durable.Event("r5_trial",event);
                            if(accept)
                            {
                                const auto audit=CompleteRestrictedAudit(trial,history,alpha,time,data,true);
                                durable.Event("r5_trial_audit",{{"iteration",J::Number(iterations)},
                                    {"name",J::Quote(name)},{"fraction",J::Number(fraction)},
                                    {"covariance",J::Numbers(audit.at("r2_covariance"))},
                                    {"qualified",J::Boolean(audit.at("r2_qualified")[0]==1.)}});
                                result.state=std::move(trial);result.residual=candidate;result.norms=candidateNorms;
                                result.domain=std::move(checked);result.constraints=current;result.box=candidateBox;
                                result.fraction=fraction;break;
                            }
                        }
                        catch(const std::exception &e)
                        {
                            event.push_back({"classification",J::Quote("trial_exception")});
                            event.push_back({"message",J::Quote(e.what())});durable.Event("r5_trial",event);
                        }
                    }
                    profile["line_search"].second+=std::chrono::duration<double>(Clock::now()-search).count();
                    return std::pair{std::move(result),details};
                };
                if(localR5 && !(DirectionReplayR6() || RefinedDirectionR6()))compare("r4_control",false,0.,true);
                auto [trial,details]=compare("r5_full",true,1.,true);
                if(!RefinedDirectionR6() && trial.fraction==0. && baseBox.Max()>0. &&
                    details.find("inward_ground_row")!=std::string::npos)
                    std::tie(trial,details)=compare("r5_inward_control",true,1.,false);
                if(trial.fraction==0. && baseBox.Max()>0. && (RefinedDirectionR6() || trial.stepPin.empty()))
                {
                    if(RefinedDirectionR6())
                        durable.Event("r6_partial_fallback",{{"iteration",J::Number(iterations)},
                            {"full_direction_constructed",J::Boolean(!trial.stepPin.empty())},
                            {"full_trials",J::Number(trialCount)},{"reserved_trials",J::Number(localR5?4:16)}});
                    std::tie(trial,details)=compare("r5_partial",true,.25,true);
                }
                if(localR5)
                {
                    m_lastFeasibleState=s;
                    durable.Event("solve_terminal",{{"converged","false"},{"terminal",J::Quote("r5_local_comparison_completed")},
                        {"iteration",J::Number(iterations)},{"useful_joint_descent",J::Boolean(trial.fraction>0.)},
                        {"profiling",profiling()}});
                    return s;
                }
                if(trial.fraction==0.)
                {
                    terminal=trial.stepPin.empty()?
                        (RefinedDirectionR6()?"r6_inner_direction_accuracy_unresolved":"r5_direction_unresolved"):
                        "r5_joint_line_search_exhausted";
                    throw std::runtime_error(terminal);
                }
                const auto update=J::Object({{"iteration",J::Number(iterations)},{"accepted","true"},
                    {"accepted_fraction",J::Number(trial.fraction)},{"matrix",matrix},{"step",trial.stepPin},
                    {"trial_count",J::Number(trialCount)},{"box",trial.box.report}});
                durable.Snapshot(iterations+1,trial.state,trial.residual,trial.constraints,trial.norms,
                    domainReport(trial.domain),profiling(),update);
                durable.Event("newton_update",{{"details",update},{"profiling",profiling()}});
                s=std::move(trial.state);residual=std::move(trial.residual);norms=std::move(trial.norms);
                constraints=std::move(trial.constraints);domain=std::move(trial.domain);alreadySaved=true;
                recordRanges(iterations+1,s);
                continue;
            }

            std::string directionDetails="null";
            const auto step=CoupledR4()?timed("r4_direction",[&]
                {return CoupledDirectionR4(domain,newton,jac,constraints,
                    [&](Vector rhs,bool transpose){return factor.Apply(std::move(rhs),transpose);},directionDetails);}):newton;
            const auto stepPin=durable.VectorPayload("iterate"+std::to_string(iterations)+".step",step);
            const auto jacStep=Multiply(jac,step);
            if(CoupledR4())
            {
                auto prediction=residual;
                for(int row=0;row<Size();++row)prediction[row]+=jacStep[row];
                const auto model=Norms(prediction,constraints,us);
                Require(*std::max_element(model.begin()+InnerPressure(),model.end())<1e-9,
                    "R4 direction failed full linearized algebraic restoration");
                durable.Event("r4_direction",{{"iteration",J::Number(iterations)},{"step",stepPin},
                    {"model_block_norms",J::Numbers(model)},{"details",directionDetails},{"profiling",profiling()}});
            }
            auto moved=s;for(int row=0;row<Size();++row)moved[row]+=step[row];
            moved=ExternalLift(moved,time,data);Vector actualStep(Size());
            for(int row=0;row<Size();++row)actualStep[row]=moved[row]-domain.effectiveState[row];
            auto traceStep=actualStep;
            for(int id:m_traceIds)traceStep[offset[Epsilon()]+id]=step[offset[Epsilon()]+id]-jacStep[offset[Epsilon()]+id];
            double maximum=1;ScalarObservation limiting;bool limited=false;
            for(const auto &point:domain.observations)
            {
                const auto row=Point(point.block,point.coordinate,point.side);
                const double delta=row.Value(point.dependentTrace?traceStep:actualStep,offset[point.block]);
                if(delta<0)
                {
                    const double bound=point.value/(-delta);
                    if(bound<maximum){maximum=std::max(0.,bound);limiting=point;limited=true;}
                }
            }
            double stepNorm=0;
            for(int block=0;block<=Gauge();++block)for(int local=0;local<size[block];++local)
            {
                const double scale=block<m_dimension?us:block==K()?K0:block==Epsilon()?E0:
                    block==Gauge()?us/m_parameters.geometry.height:us*us;
                stepNorm=std::max(stepNorm,std::abs(step[offset[block]+local])/scale);
            }
            double fraction=maximum<1?.995*maximum:1.;bool accepted=false;
            Vector nextState,nextResidual,nextNorms;std::map<int,double> nextConstraints;
            ScalarDomainResult nextDomain;const auto searchStarted=Clock::now();++profile["line_search"].first;
            for(int backtrack=0;backtrack<32 && fraction>=std::ldexp(1.,-40);++backtrack,fraction*=.5)
            {
                auto trial=s;for(int row=0;row<Size();++row)trial[row]+=fraction*step[row];
                const auto trialName="iterate"+std::to_string(iterations)+".trial"+std::to_string(backtrack);
                const auto trialPin=durable.VectorPayload(trialName+".state",trial);
                auto checked=timed("scalar_domain",[&]{return ScalarDomain(trial,time,data,true);});
                if(!checked.admissible || !checkWall(checked))
                {
                    durable.Event("line_search_trial",{{"iteration",J::Number(iterations)},{"backtrack",J::Number(backtrack)},
                        {"fraction",J::Number(fraction)},{"classification",J::Quote(checked.admissible?"r4_wall_trace_rejection":"scalar_domain_rejection_before_residual")},
                        {"state",trialPin},{"scalar_domain",domainReport(checked)}});continue;
                }
                try
                {
                    auto candidate=timed("residual",[&]{return Residual(trial,history,alpha,time,data);});
                    auto current=timed("constraints",[&]{return Constraints(trial,time,data);});
                    auto candidateNorms=timed("mass_dual_norms",[&]{return Norms(candidate,current,us);});
                    const auto candidatePin=durable.VectorPayload(trialName+".residual",candidate);
                    const bool descent=Maximum(candidateNorms)<(1.-1e-4*fraction)*error;
                    auto rechecked=timed("scalar_domain",[&]{return ScalarDomain(trial,time,data,true);});
                    const bool accept=descent && rechecked.admissible && checkWall(rechecked);
                    durable.Event("line_search_trial",{{"iteration",J::Number(iterations)},{"backtrack",J::Number(backtrack)},
                        {"fraction",J::Number(fraction)},{"classification",J::Quote(accept?"feasible_accepted_descent":"valid_non_descent_or_recheck_failure")},
                        {"state",trialPin},{"residual",candidatePin},{"scaled_block_norms",J::Numbers(candidateNorms)},
                        {"scalar_domain",domainReport(rechecked)},
                        {"profile_branches_changed",J::Boolean(domain.profileBranches!=rechecked.profileBranches)}});
                    if(accept)
                    {
                        accepted=true;nextState=std::move(trial);nextResidual=std::move(candidate);
                        nextNorms=std::move(candidateNorms);nextConstraints=std::move(current);nextDomain=std::move(rechecked);break;
                    }
                }
                catch(const std::exception &errorValue)
                {
                    durable.Event("line_search_trial",{{"iteration",J::Number(iterations)},{"backtrack",J::Number(backtrack)},
                        {"fraction",J::Number(fraction)},{"classification",J::Quote("evaluable_scalar_state_residual_exception")},
                        {"state",trialPin},{"exception_type",J::Quote(J::ExceptionType(errorValue))},
                        {"exception_message",J::Quote(errorValue.what())}});
                }
            }
            profile["line_search"].second+=std::chrono::duration<double>(Clock::now()-searchStarted).count();
            const auto update=J::Object({{"iteration",J::Number(iterations)},{"accepted",J::Boolean(accepted)},
                {"matrix",matrix},{"step",stepPin},{"dimensionless_step",J::Number(stepNorm)},
                {"linear_scaled_block_norms",J::Numbers(linearNorms)},
                {"fixed_location_maximum_fraction",J::Number(maximum)},{"boundary_limited",J::Boolean(limited)},
                {"limiting_inventory",J::Quote(limiting.inventory)},{"limiting_point",J::Numbers(limiting.coordinate)},
                {"limiting_value",J::Number(limiting.value)},{"accepted_fraction",J::Number(accepted?fraction:0)}});
            if(!accepted)
            {
                durable.Event("newton_update",{{"details",update},{"profiling",profiling()}});
                terminal=maximum<std::ldexp(1.,-40)?"feasible_step_collapsed":"feasible_line_search_exhausted";
                throw std::runtime_error(terminal);
            }
            const double nextError=Maximum(nextNorms);
            // Publish the completed accepted state before any further action,
            // including the stagnation decision and the next loop's J start.
            durable.Snapshot(iterations+1,nextState,nextResidual,nextConstraints,nextNorms,
                domainReport(nextDomain),profiling(),update);
            durable.Event("newton_update",{{"details",update},{"profiling",profiling()}});
            s=std::move(nextState);residual=std::move(nextResidual);norms=std::move(nextNorms);
            constraints=std::move(nextConstraints);domain=std::move(nextDomain);alreadySaved=true;
            recordRanges(iterations+1,s);
            if(localR4)
            {
                m_lastFeasibleState=s;
                durable.Event("solve_terminal",{{"converged","false"},{"terminal",J::Quote("r4_local_experiment_completed")},
                    {"iteration",J::Number(iterations+1)},{"scaled_block_norms",J::Numbers(norms)},{"profiling",profiling()}});
                return s;
            }
            stagnant=nextError>error*(1.-1e-4)?stagnant+1:0;
            if(stagnant>=3){terminal="feasible_residual_stagnation";throw std::runtime_error(terminal);}
        }
        throw std::runtime_error("unreachable R3 solve");
    }
    catch(const std::exception &error)
    {
        m_lastFeasibleState=s;
        if(terminal=="started")terminal=std::string(J::ExceptionType(error))+": "+error.what();
        durable.Event("solve_terminal",{{"converged","false"},{"terminal",J::Quote(terminal)},
            {"iteration",J::Number(iterations)},{"last_completed_iteration",J::Number(iterations+(alreadySaved?1:0))},
            {"scaled_block_norms",J::Numbers(norms)},{"profiling",profiling()}});
        throw;
    }
}

RANSWallEvolutionPatch::Diagnostic RANSWallEvolutionPatch::CompleteRestrictedAudit(
    const Vector &state,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,bool rangesOnly) const
{
    const auto domain=ScalarDomain(state,time,data,true);Diagnostic audit;
    audit["r2_scalar_domain"]={domain.admissible?1.:0.,domain.minimumMargin,domain.minimumInnerK,
        domain.minimumOuterK,domain.minimumOuterEpsilon,domain.minimumSampledK,
        domain.minimumTraceEpsilon,domain.minimumB,domain.minimumLambda};
    if(!domain.admissible){audit["r2_qualified"]={0};return audit;}
    if(!rangesOnly)audit=Audit(domain.effectiveState,history,alpha,time,data,m_robustUstar,true);
    audit["r2_scalar_domain"]={1,domain.minimumMargin,domain.minimumInnerK,domain.minimumOuterK,
        domain.minimumOuterEpsilon,domain.minimumSampledK,domain.minimumTraceEpsilon,domain.minimumB,domain.minimumLambda};
    const double K0=m_parameters.model.Q()*m_robustUstar*m_robustUstar;
    Vector outer{0,1e300,-1e300,1e300,-1e300,domain.minimumSampledK/K0,-1e300};
    double innerCov=1e300,outerCov=1e300;Vector innerPoint,outerPoint,gradientPoint;
    for(const auto &p:domain.observations)
    {
        if(p.block!=K() && !p.dependentTrace)continue;
        const double r=Distance(p.coordinate);const bool inner=RestrictedInner(p);
        const auto row=RestrictedGradientRow(p);
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
        {
            const double departure=G[a][b]/g0-m_referenceFlow[a]*m_referenceNormal[b];
            if(std::abs(departure)>outer[0])
            {
                outer[0]=std::abs(departure);
                gradientPoint={p.coordinate[0],p.coordinate[1],p.coordinate[2],double(a),double(b),
                    departure,double(p.geometry),double(p.side)};
            }
        }
        outer[1]=std::min(outer[1],k/K0);outer[2]=std::max(outer[2],k/K0);
        outer[3]=std::min(outer[3],epsilon/e0);outer[4]=std::max(outer[4],epsilon/e0);
    }
    for(const auto &p:domain.observations)if(p.inventory=="physical_sampled_K")outer[6]=std::max(outer[6],p.value/K0);
    const bool box=outer[0]<=.001 && outer[1]>=.99 && outer[2]<=1.01 && outer[3]>=.98 && outer[4]<=1.02 &&
        outer[5]>=.99 && outer[6]<=1.01;
    const bool covariance=innerCov>=-1e-10*K0 && outerCov>=-1e-10*K0;
    if(rangesOnly || GroundSafeguard())audit["r4_outer_gradient_witness"]=gradientPoint;
    audit["r2_outer_range"]=outer;audit["r2_covariance"]={innerCov,outerCov};
    audit["r2_inner_minimum_point"]=innerPoint;audit["r2_outer_minimum_point"]=outerPoint;
    audit["r2_qualified"]={box && covariance && domain.minimumLambda>=.5?1.:0.};
    return audit;
}
} // namespace Nektar
