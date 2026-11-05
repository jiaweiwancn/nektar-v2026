///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallCutBudget.h>
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

enum class FullRestorationOutcome { DirectionUnresolved, LineSearchExhausted, Accepted };

constexpr int R6RestorationTrialLimit(bool local)
{
    return local?4:16;
}

constexpr bool R6TryPartialRestoration(FullRestorationOutcome outcome,double violation)
{
    return outcome!=FullRestorationOutcome::Accepted && violation>0.;
}

static_assert(2*R6RestorationTrialLimit(false)==32 && 2*R6RestorationTrialLimit(true)==8);
static_assert(R6TryPartialRestoration(FullRestorationOutcome::DirectionUnresolved,1.));
static_assert(R6TryPartialRestoration(FullRestorationOutcome::LineSearchExhausted,1.));
static_assert(!R6TryPartialRestoration(FullRestorationOutcome::Accepted,1.) &&
              !R6TryPartialRestoration(FullRestorationOutcome::LineSearchExhausted,0.));

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
    RANSWallR3::DurableProgress *diagnostic,const std::string &tag,const Vector &residual,
    double radius,double *guardWork,const std::array<double,4> &interior,
    const std::vector<DirectionRestriction> &extra,bool scanOnly) const
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

    constexpr double reserve=.05, near=.005, feasibility=5e-12;
    struct Witness
    {
        int block=-1,axis=-1,side=0,coefficient=-1;
        double scale=1;
        bool trace=false,general=false;
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
        for(int i=0;i<ceiling.size();++i)ceiling[i]=(1.-restoration)*ceiling[i]-interior[i];
    }
    const bool r12=ExperimentR12();
    const int maxExchange=r12?m_r12CutLimit:(outer.empty()?40:96);
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
            const double native=w.general?std::inner_product(o.row.begin(),o.row.end(),direction.begin(),0.):w.coefficient>=0?direction[w.coefficient]/w.scale:
                (w.axis<0?w.point.Value(w.trace?trace:effective,offset[w.block]):
                    w.point.Gradient(effective,offset[w.block])[w.axis])/w.scale;
            const std::string name=tag+(r12?".pass"+std::to_string(exchange):"")+".cut"+std::to_string(rows.size());
            auto unscaled=cut.normal;
            for(double &x:unscaled)x*=cut.length;
            rows.push_back(J::Object({{"id",J::Number(cut.id)},{"kind",J::Quote(cut.kind)},{"general_affine_row",J::Boolean(w.general)},
                {"point",w.general?"null":J::Numbers(cut.point)},{"geometry",w.general?"null":J::Number(w.point.geometry)},
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
            {"exchanges",J::List(exchanges)},{"z",diagnostic->VectorPayload(tag+(r12?".pass"+std::to_string(exchange):"")+".z",z)},
            {"direction",diagnostic->VectorPayload(tag+(r12?".pass"+std::to_string(exchange):"")+".direction",direction)},
            {"newton",diagnostic->VectorPayload(tag+(r12?".pass"+std::to_string(exchange):"")+".newton",newton)},
            {"Bz",diagnostic->VectorPayload(tag+(r12?".pass"+std::to_string(exchange):"")+".Bz",metricAction(z,false))}});
        published=true;
    };
    for(;exchange<(r12?CutBudgetR12{maxExchange}.Scans():maxExchange);++exchange)
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
                         bool projected,int side,const RANSWallGlobalPoint *row,int coefficient=-1,const Vector *general=nullptr)
        {
            worstNative=std::max(worstNative,lower-observed);
            if(refined)
            {
                long double value=0;
                if(general)for(int j=0;j<Size();++j)value+=(long double)(*general)[j]*direction[j];
                else if(coefficient>=0)value=(long double)direction[coefficient]/scale;
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
                selected={block,axis,side,coefficient,scale,projected,general!=nullptr,row?*row:RANSWallGlobalPoint{}};
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
        std::vector<std::string> extraAudit;
        for(const auto &cut:extra)
        {
            Require(cut.row.size()==Size(),"Invalid R10 general affine row");
            inspect(cut.lower,std::inner_product(cut.row.begin(),cut.row.end(),direction.begin(),0.),
                cut.name,{},[&]{return cut.row;},-1,-1,1.,false,0,nullptr,-1,&cut.row);
            extraAudit.push_back(J::Object({{"name",J::Quote(cut.name)},{"point","null"},
                {"lower",J::Number(cut.lower)},{"native_gap",J::Number(
                    std::inner_product(cut.row.begin(),cut.row.end(),direction.begin(),0.)-cut.lower)},
                {"independent_accumulation_gap",J::Number(dot(cut.row,direction)-cut.lower)}}));
        }
        if(diagnostic && !extra.empty() && (scanOnly || violated.empty()))
            diagnostic->Event("r10_general_row_scan",{{"name",J::Quote(tag)},{"rows",J::List(extraAudit)}});
        if(r12 && diagnostic && !scanOnly)
        {
            Vector nonnegative=multipliers,dualVector(Size());
            double dualLinear=0.;
            for(int i=0;i<active.size();++i)
            {
                nonnegative[i]=std::max(0.,nonnegative[i]);
                dualLinear+=active[i].rhs*nonnegative[i];
                axpy(dualVector,nonnegative[i],active[i].normal);
            }
            Vector activeIds,diagonal;
            for(int i=0;i<active.size();++i){activeIds.push_back(active[i].id);diagonal.push_back(R[i][i]);}
            diagnostic->Event("r12_qp_pass",{{"name",J::Quote(tag)},{"oracle_scan",J::Number(exchange+1)},
                {"admitted",J::Number(original.size())},{"releases",J::Number(releases)},
                {"active_ids",J::Numbers(activeIds)},{"qr_diagonal",J::Numbers(diagonal)},
                {"primal_objective",J::Number(.5*dot(z,z))},
                {"nonnegative_dual_bound",J::Number(dualLinear-.5*dot(dualVector,dualVector))},
                {"worst_native_violation",J::Number(worstNative)},
                {"worst_extended_violation",J::Number(worstExtended)},
                {"direction",diagnostic->VectorPayload(tag+".scan"+std::to_string(exchange)+".d",direction)},
                {"z",diagnostic->VectorPayload(tag+".scan"+std::to_string(exchange)+".z",z)}});
        }
        if(scanOnly)
        {
            details=J::Object({{"passed",J::Boolean(violated.empty())},
                {"worst_native_violation",J::Number(worstNative)},
                {"worst_extended_violation",J::Number(worstExtended)}});
            if(diagnostic)diagnostic->Event("r10_anchor_original_rows",{{"name",J::Quote(tag)},{"scan",details}});
            return violated.empty()?direction:Vector{};
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
                    {"point",row.kind=="native_merit_slope"?"null":J::Numbers(Vector(row.point.begin(),row.point.end()))},
                    {"multiplier",J::Number(multipliers[i]/row.length)}}));
            }
            const double stationarityNorm=std::sqrt(dot(stationarity,stationarity));
            const double objective=.5*dot(z,z);
            if(guardWork)
            {
                *guardWork=0.;
                for(int i=0;i<active.size();++i)
                    if(active[i].kind=="inward_ground_row")
                    {
                        const auto &point=original[active[i].id].witness.point;
                        const double k=point.Value(domain.effectiveState,offset[K()])/K0;
                        if(k>0.)*guardWork+=k*multipliers[i]/active[i].length;
                    }
                *guardWork/=std::max(objective,1e-30);
            }
            Require(stationarityNorm<=1e-9 && activeError<=1e-9 && complementarity<=1e-9,
                "R4 constrained solve failed its KKT check");
            details=J::Object({{"method",J::Quote(active.empty()?"ordinary_newton":"residual_space_active_QR")},
                {"active",J::List(activeRows)},{"releases",J::Number(releases)},
                {"dependent_exchanges",J::Number(dependentExchanges)},
                {"exchange_passes",J::Number(exchange+1)},{"trust_radius",J::Number(radius)},
                {"model_feasibility_tolerance",J::Number(feasibility)},{"stationarity_l2",J::Number(stationarityNorm)},
                {"active_equality_error",J::Number(activeError)},{"complementarity",J::Number(complementarity)},
                {"scaled_residual_objective",J::Number(objective)},
                {"witness_k",J::Number(witnessValue)},
                {"witness_newton_delta",J::Number(witnessRow.Value(newton,offset[K()]))},
                {"witness_direction_delta",J::Number(witnessRow.Value(effective,offset[K()]))},
                {"outer_restoration_fraction",outer.empty()?"null":J::Number(restoration)},
                {"outer_linear_ceiling",J::Numbers(ceiling)},
                {"inward_positive_ground_rule",J::Boolean(inward)}});
            record(nullptr,R,effective,trace,"passed_original_cut_scan");
            return direction;
        }
        // R12 performs one last reconstruction and complete scan after the last
        // admitted cut. This check cannot discover/admit another cut.
        if(r12 && !CutBudgetR12{maxExchange}.CanAdmit(exchange))
        {
            record(nullptr,R,effective,trace,"post_limit_updated_candidate_infeasible");
            throw std::runtime_error("R12 post-limit updated candidate still violates original rows: "+
                kind+" violation "+J::Number(violation));
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
                {"kind",J::Quote(kind)},{"point",kind=="native_merit_slope"?"null":J::Numbers(location)},{"active",J::Numbers(activeIds)},
                {"exact_previous_alias_ids",J::Numbers(aliases)},{"violation",J::Number(violation)}}));
            if(r12)
            {
                const std::string prefix=tag+".original"+std::to_string(id);
                diagnostic->Event("r12_original_cut",{{"name",J::Quote(tag)},{"id",J::Number(id)},
                    {"kind",J::Quote(kind)},{"point",J::Numbers(location)},{"rhs",J::Number(rhs)},
                    {"identity",J::Quote(RANSWallR3::Hash([&]{auto v=violated;v.push_back(rhs);return v;}()))},
                    {"aliases",J::Numbers(aliases)},{"length",J::Number(length)},
                    {"transformed_rhs",J::Number(incoming.rhs*length)},
                    {"row",diagnostic->VectorPayload(prefix+".row",violated)},
                    {"normal",diagnostic->VectorPayload(prefix+".normal",[&]{auto v=incoming.normal;for(auto &x:v)x*=length;return v;}())}});
            }
            if(exchange+1==maxExchange)record(&incoming,R,effective,trace,"exchange_limit_candidate");
            if(refined && std::any_of(active.begin(),active.end(),[&](const Cut &cut)
                {return original[cut.id].row==violated && original[cut.id].rhs==rhs;}))
            {
                if(!published)record(&incoming,R,effective,trace,"active_original_row_accuracy_unresolved");
                throw std::runtime_error("R6 active original row still violates the verified cut scan");
            }
        }
        lastCut=kind+" violation "+J::Number(violation)+
            (kind=="native_merit_slope"?"":" at "+J::Numbers(Vector(location.begin(),location.end())));
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

RANSWallEvolutionPatch::R8Trial RANSWallEvolutionPatch::RunR8Comparison(const Vector &state,const Vector &history,double alpha,double time,
    const RANSWallEvolutionData &data,const ScalarDomainResult &baseDomain,const Vector &residual,
    const std::map<int,double> &constraints,const WallMatrix &jac,const Vector &newton,
    const std::function<Vector(Vector,bool)> &inverse,const R8Ray &ray,RANSWallR3::DurableProgress &durable,
    int iteration,const std::string &comparison,std::map<double,R8Trial> *samples,
    std::map<Vector,R8Trial> *cache,R11TrialWork *work) const
{
    using Clock=std::chrono::steady_clock;
    const bool pilot=PilotR8() || ExperimentR9() || ExperimentR10() || ExperimentR11() || ExperimentR12();
    const std::string prefix="iterate"+std::to_string(iteration)+"."+comparison;
    constexpr double fit=.125,holdout=.0625;
    const double us=m_robustUstar,K0=m_parameters.model.Q()*us*us,H=m_parameters.geometry.height;
    const double E0=us*us*us/(m_parameters.model.Kappa()*H);
    const auto baseNorms=Norms(residual,constraints,us);
    const double baseE=Maximum(baseNorms);
    const auto outer=OuterBoundsR5(baseDomain);
    const auto baseBox=OuterViolationR5(outer);
    Require(baseBox.Max()==0. && ray.step.size()==Size() && Multiply(jac,ray.step)==ray.action &&
        (pilot || (ray.residuals.size()==4 && ray.residuals.count(1.) && ray.residuals.count(.5) &&
        ray.residuals.count(.25) && ray.residuals.count(fit))),"R8 binding or base box differs");
    bool baseWallOK=false;
    const auto baseWall=WallTraceJSON(baseDomain.effectiveState,baseWallOK);
    Require(baseWallOK,"R8 base wall safeguard failed");
    const auto baseExtrema=WallTraceR4(baseDomain.effectiveState);
    const double baseK=std::min_element(baseExtrema.begin(),baseExtrema.end(),
        [](const auto &a,const auto &b){return a.value<b.value;})->value;
    const std::array<double,7> scales{us,us,K0,E0,us*us,us*us,us/H};
    auto normStep=[&](const Vector &step)
    {
        double norm=0.;
        for(int b=0;b<=Gauge();++b)for(int i=offset[b];i<offset[b]+size[b];++i)
            norm=std::max(norm,std::abs(step[i])/scales[b]);
        return norm;
    };
    const double retainedNorm=normStep(ray.step);
    Require(pilot?retainedNorm<=.125+5e-12:std::abs(retainedNorm-.125)<5e-12,"R8 radius differs");
    durable.Event("r8_base",{{"iteration",J::Number(iteration)},{"scaled_block_norms",J::Numbers(baseNorms)},{"box",baseBox.report},
        {"wall",baseWall},{"scalar",ScalarDomainJSON(baseDomain)},{"retained_step_norm",J::Number(retainedNorm)},
        {"maximum_new_native_trials",pilot?"9":"20"},{"fit_fraction","0.125"},{"holdout_fraction","0.0625"}});
    auto along=[&](const Vector &step,double f)
    {
        auto q=state;
        for(int i=0;i<Size();++i)q[i]+=f*step[i];
        return q;
    };
    auto linear=[&](const Vector &action,double f)
    {
        auto r=residual;
        for(int i=0;i<Size();++i)r[i]+=f*action[i];
        return r;
    };
    const auto fullModel=Norms(linear(ray.action,1.),constraints,us);
    Require(*std::max_element(fullModel.begin()+InnerPressure(),fullModel.end())<1e-9,
        "R8 direction failed full linearized algebraic restoration");
    auto difference=[&](const Vector &a,const Vector &b)
    {
        Vector delta(Size());
        for(int i=0;i<Size();++i)delta[i]=a[i]-b[i];
        return delta;
    };
    if(!pilot)for(const auto &[f,r]:ray.residuals)
        durable.Event("r8_retained_remainder",{{"fraction",J::Number(f)},{"actual",J::Numbers(Norms(r,constraints,us))},
            {"linear",J::Numbers(Norms(linear(ray.action,f),constraints,us))},
            {"remainder",J::Numbers(Norms(difference(r,linear(ray.action,f)),constraints,us))}});

    int evaluations=0;
    auto evaluate=[&](const std::string &name,const Vector &step,const Vector &action,double f,const Vector *model=nullptr)
    {
        Require(++evaluations<=(pilot?9:20),"R8 trial budget exhausted");
        const auto started=Clock::now();
        R8Trial t;
        t.name=name;t.f=f;t.state=along(step,f);
        if(cache)
        {
            const auto previous=cache->find(t.state);
            if(previous!=cache->end())
            {
                t=previous->second;t.name=name;t.f=f;
                t.accepted=t.qualified && t.E<(1.-1e-4*f)*baseE;
                bool wallOK=false;
                const auto predicted=Norms(linear(action,f),constraints,us);
                double trace=0.;
                for(int id:m_traceIds)trace=std::max(trace,std::abs(t.residual[offset[Epsilon()]+id])/E0);
                const double seconds=std::chrono::duration<double>(Clock::now()-started).count();
                if(work){++work->reused;work->actual+=seconds;work->charged+=t.evaluationSeconds;}
                durable.Event("r8_trial",{{"iteration",J::Number(iteration)},{"comparison",J::Quote(comparison)},
                    {"name",J::Quote(name)},{"fraction",J::Number(f)},{"cache_reused","true"},
                    {"state",durable.VectorPayload(prefix+name+".state",t.state)},
                    {"residual",durable.VectorPayload(prefix+name+".residual",t.residual)},
                    {"scalar",ScalarDomainJSON(t.domain)},{"wall",WallTraceJSON(t.domain.effectiveState,wallOK)},
                    {"actual",J::Numbers(t.norms)},{"linear",J::Numbers(predicted)},{"box",t.box.report},
                    {"covariance",J::Numbers(t.covariance)},{"matching_trace",J::Number(trace)},
                    {"classification",J::Quote(t.accepted?"accepted":"residual_non_descent")},
                    {"seconds",J::Number(seconds)},{"charged_seconds",J::Number(t.evaluationSeconds)}});
                return t;
            }
        }
        const auto checked=ScalarDomain(t.state,time,data,true);
        std::string verdict=checked.admissible?"ground":"scalar_domain";
        std::vector<std::pair<std::string,std::string>> fields{{"iteration",J::Number(iteration)},
            {"name",J::Quote(name)},{"comparison",J::Quote(comparison)},{"fraction",J::Number(f)},
            {"state",durable.VectorPayload(prefix+name+".state",t.state)},{"scalar",ScalarDomainJSON(checked)},
            {"profile_branches_changed",J::Boolean(checked.profileBranches!=baseDomain.profileBranches)}};
        if(checked.admissible)
        {
            bool wallOK=false;
            fields.emplace_back("wall",WallTraceJSON(checked.effectiveState,wallOK));
            if(wallOK)
            {
                const auto extrema=WallTraceR4(checked.effectiveState);
                t.ground=std::min_element(extrema.begin(),extrema.end(),
                    [](const auto &a,const auto &b){return a.value<b.value;})->value;
                const auto box=OuterViolationR5(OuterBoundsR5(checked));
                const auto actualConstraints=Constraints(t.state,time,data);
                t.residual=Residual(t.state,history,alpha,time,data);
                const auto norms=Norms(t.residual,actualConstraints,us);
                t.E=Maximum(norms);
                const auto predicted=Norms(linear(action,f),constraints,us);
                const auto audit=CompleteRestrictedAudit(t.state,history,alpha,time,data,true);
                const bool descent=t.E<(1.-1e-4*f)*baseE;
                t.qualified=box.Max()==0. && audit.at("r2_qualified")[0]==1.;
                t.accepted=descent && t.qualified;
                verdict=t.accepted?"accepted":!descent?"residual_non_descent":box.Max()!=0.?"outer_box":"physical_inventory";
                double trace=0.;
                for(int id:m_traceIds)trace=std::max(trace,std::abs(t.residual[offset[Epsilon()]+id])/E0);
                fields.insert(fields.end(),{{"box",box.report},{"actual",J::Numbers(norms)},{"linear",J::Numbers(predicted)},
                    {"remainder",J::Numbers(Norms(difference(t.residual,linear(action,f)),constraints,us))},
                    {"linear_agreement",baseE>Maximum(predicted)?J::Number((baseE-t.E)/(baseE-Maximum(predicted))):"null"},
                    {"matching_trace",J::Number(trace)},{"covariance",J::Numbers(audit.at("r2_covariance"))},
                    {"residual",durable.VectorPayload(prefix+name+".residual",t.residual)}});
                if(model)
                {
                    const auto modelNorm=Norms(*model,constraints,us);
                    fields.emplace_back("quadratic",J::Numbers(modelNorm));
                    fields.emplace_back("quadratic_error",J::Numbers(Norms(difference(t.residual,*model),constraints,us)));
                    fields.emplace_back("quadratic_agreement",baseE>Maximum(modelNorm)?
                        J::Number((baseE-t.E)/(baseE-Maximum(modelNorm))):"null");
                }
                if(t.accepted || cache)
                {
                    t.covariance=audit.at("r2_covariance");
                    t.box=box;t.constraints=actualConstraints;t.norms=norms;t.domain=checked;
                }
            }
        }
        t.evaluationSeconds=std::chrono::duration<double>(Clock::now()-started).count();
        if(work){++work->fresh;work->actual+=t.evaluationSeconds;work->charged+=t.evaluationSeconds;}
        if(cache && !t.residual.empty())cache->emplace(t.state,t);
        fields.emplace_back("classification",J::Quote(verdict));
        fields.emplace_back("seconds",J::Number(t.evaluationSeconds));
        if(cache)fields.emplace_back("cache_reused","false");
        durable.Event("r8_trial",fields);
        return t;
    };
    std::map<double,R8Trial> sameRay;
    auto tryRay=[&](double f,const std::string &name,const Vector *model=nullptr)->const R8Trial &
    {
        auto p=sameRay.find(f);
        if(p!=sameRay.end())return p->second;
        return sameRay.emplace(f,evaluate(name,ray.step,ray.action,f,model)).first->second;
    };
    const R8Trial *ordinary=nullptr,*bestRay=nullptr;
    if(!pilot)for(double f:{.0625,.03125,.015625,.0078125})
    {
        const auto &t=tryRay(f,"r8_ordinary"+std::to_string(sameRay.size()));
        if(t.accepted){ordinary=bestRay=&t;break;}
    }
    if(!pilot)durable.Event("r8_ordinary_outcome",{{"accepted",J::Boolean(ordinary)},
        {"fraction",J::Number(ordinary?ordinary->f:0.)},{"new_evaluations",J::Number(evaluations)}});
    if(!pilot && !bestRay)
        for(double f:{.00390625,.001953125})
        {
            const auto &t=tryRay(f,"r8_small"+std::to_string(sameRay.size()));
            if(t.accepted){bestRay=&t;break;}
        }

    const Vector *fitResidual=nullptr,*heldResidual=nullptr;
    if(pilot && ray.residuals.empty())
    {
        const auto &fitTrial=tryRay(fit,"r8_live_fit");
        const auto &heldTrial=tryRay(holdout,"r8_live_holdout");
        if(fitTrial.qualified && heldTrial.qualified)
        {
            fitResidual=&fitTrial.residual;heldResidual=&heldTrial.residual;
        }
        if(fitTrial.accepted)bestRay=&fitTrial;
        if(heldTrial.accepted && (!bestRay || heldTrial.E<bestRay->E))bestRay=&heldTrial;
    }
    else
    {
        fitResidual=&ray.residuals.at(fit);
        if(pilot)heldResidual=&ray.residuals.at(holdout);
        else if(!sameRay.at(holdout).residual.empty())heldResidual=&sameRay.at(holdout).residual;
    }
    if((pilot || bestRay) && fitResidual && heldResidual)
    {
        const auto started=Clock::now();
        auto curvature=difference(*fitResidual,linear(ray.action,fit));
        for(double &v:curvature)v/=fit*fit;
        auto quadratic=[&](double f)
        {
            auto r=linear(ray.action,f);
            for(int i=0;i<Size();++i)r[i]+=f*f*curvature[i];
            return r;
        };
        const auto validation=quadratic(holdout);
        const auto validationNorm=Norms(validation,constraints,us);
        const auto heldNorm=Norms(*heldResidual,constraints,us);
        const double heldE=Maximum(heldNorm);
        const bool agrees=std::abs(Maximum(validationNorm)-heldE)<=std::max(2e-6,.5*std::abs(baseE-heldE));
        // The interval is bracketed by the held-out point and the
        // predeclared small-step range. Poor validation halves the cap.
        double candidate=0.,best=baseE;
        for(int i=16;i<=(agrees?128:64);++i)
        {
            const double f=i/2048.;
            const double E=Maximum(Norms(quadratic(f),constraints,us));
            if(E<best){best=E;candidate=f;}
        }
        const int before=evaluations;
        const double predictorSeconds=std::chrono::duration<double>(Clock::now()-started).count();
        const R8Trial *curved=nullptr;
        if(candidate>0.)
        {
            const auto prediction=quadratic(candidate);
            const auto &t=tryRay(candidate,"r8_curvature",&prediction);
            if(t.accepted)curved=&t;
        }
        // The ordinary control is a cached bounded fallback on this
        // exact ray. It requires no repeated native evaluation here.
        if(!curved)curved=bestRay;
        if(curved && (!bestRay || curved->E<bestRay->E))bestRay=curved;
        durable.Event("r8_curvature_outcome",{{"iteration",J::Number(iteration)},
            {"fit_fraction",J::Number(fit)},{"holdout_fraction",J::Number(holdout)},
            {"heldout_predicted",J::Numbers(validationNorm)},{"heldout_actual",J::Numbers(heldNorm)},
            {"heldout_vector_error",J::Numbers(Norms(difference(*heldResidual,validation),constraints,us))},
            {"comparison",J::Quote(comparison)},
            {"samples",J::Quote(pilot?(ray.residuals.empty()?"current_live":"current_bound"):"fixed_comparison")},
            {"validation_agrees",J::Boolean(agrees)},{"proposed_fraction",J::Number(candidate)},
            {"used_fraction",J::Number(curved?curved->f:0.)},{"new_evaluations",J::Number(evaluations-before)},
            {"standalone_holdout_evaluations","1"},{"model_construction_seconds",J::Number(predictorSeconds)}});
    }

    if(pilot)
    {
        if(!bestRay)
            for(int i=0;i<8;++i)
            {
                const auto &t=tryRay(std::ldexp(fit,-i),"r8_fallback"+std::to_string(i));
                if(t.accepted){bestRay=&t;break;}
            }
        durable.Event("r8_work",{{"iteration",J::Number(iteration)},{"comparison",J::Quote(comparison)},
            {"native_trial_states",J::Number(evaluations)},{"selected",J::Quote(bestRay?bestRay->name:"none")}});
        if(samples)*samples=sameRay;
        return bestRay?*bestRay:R8Trial{};
    }

    const auto &held=sameRay.at(holdout);
    if(heldResidual)
    {
        const auto started=Clock::now();
        const auto zero=ResidualResolutionComponents(state,history,alpha,time,data);
        const auto middle=ResidualResolutionComponents(held.state,history,alpha,time,data);
        const auto fitState=along(ray.step,fit);
        const auto upper=ResidualResolutionComponents(fitState,history,alpha,time,data);
        Vector saturation;
        for(const auto *q:{&state,&held.state,&fitState})
            saturation.push_back(Audit(*q,history,alpha,time,data,us,true).at("inner_stress")[0]);
        Vector reconstruct;
        for(const auto &[components,value]:{std::pair{&zero,&residual},std::pair{&middle,&held.residual},
                                            std::pair{&upper,fitResidual}})
            reconstruct.push_back(Maximum(Norms(difference(components->terms.at("total_constrained"),*value),constraints,us)));
        Require(Maximum(reconstruct)<1e-8,"R8 component residual reconstruction differs");
        std::vector<std::pair<std::string,std::string>> components;
        for(const std::string key:{"time_mass","conservative_advection","k_diffusion","inner_k_production",
            "inner_k_destruction","outer_k_production","outer_k_destruction","total_constrained"})
        {
            const auto &a=zero.terms.at(key),&b=middle.terms.at(key),&c=upper.terms.at(key);
            Vector bend(Size());
            for(int i=0;i<Size();++i)bend[i]=holdout*holdout*((c[i]-a[i])/fit-(b[i]-a[i])/holdout)/(fit-holdout);
            components.emplace_back(key,J::Number(Norms(bend,constraints,us)[K()]));
        }
        durable.Event("r8_components",{{"scaled_k_second_divided_differences_times_holdout_squared",J::Object(components)},
            {"residual_reconstruction_error",J::Numbers(reconstruct)},{"states","[0,0.0625,0.125]"},
            {"active_inner_stress_volume",J::Numbers(saturation)},
            {"seconds",J::Number(std::chrono::duration<double>(Clock::now()-started).count())}});
    }

    R8Trial shorter;
    // A sub-per-mille improvement is a weak same-ray outcome even when
    // the selector reduces trial cost. Re-solve rather than rescaling.
    if(!bestRay || baseE-bestRay->E<.001*baseE)
    {
        double lower=0.;
        for(int i=0;i<Size();++i)if(i>=offset[InnerPressure()] || constraints.count(i))
        {
            long double capacity=0.;
            for(int b=0;b<=Gauge();++b)
                for(int j=offset[b];j<offset[b]+size[b];++j)capacity+=std::abs(jac[i][j])*scales[b];
            if(capacity>0.)lower=std::max(lower,static_cast<double>(std::abs(residual[i])/capacity));
        }
        const double displacement=retainedNorm*(bestRay?bestRay->f:.0078125);
        double radius=std::min(.0625,std::max(4.*displacement,2.*lower));
        for(int attempt=0;attempt<2;++attempt)
        {
            const auto started=Clock::now();
            const std::string tag="r8_radius"+std::to_string(attempt);
            durable.Event("r8_radius_selection",{{"name",J::Quote(tag)},{"radius",J::Number(radius)},
                {"retained_displacement",J::Number(displacement)},{"necessary_algebraic_radius",J::Number(lower)}});
            std::string details;
            Vector step;
            try
            {
                Require(radius>=lower,"R8 radius violates necessary algebraic bound");
                step=CoupledDirectionR4(baseDomain,newton,jac,constraints,inverse,details,
                    outer,1.,false,&durable,tag,residual,radius);
            }
            catch(const std::exception &e)
            {
                durable.Event("r8_direction_failure",{{"name",J::Quote(tag)},{"failure",J::Quote(e.what())},
                    {"seconds",J::Number(std::chrono::duration<double>(Clock::now()-started).count())}});
                radius=std::min(.09375,2.*radius);
                continue;
            }
            const auto action=Multiply(jac,step);
            const auto model=Norms(linear(action,1.),constraints,us);
            Require(*std::max_element(model.begin()+InnerPressure(),model.end())<1e-9,
                "R8 smaller radius failed full algebraic restoration");
            const double norm=normStep(step);
            durable.Event("r8_direction",{{"name",J::Quote(tag)},{"radius",J::Number(radius)},
                {"step_norm",J::Number(norm)},{"details",details},{"linear",J::Numbers(model)},
                {"step",durable.VectorPayload(tag+".step",step)},{"action",durable.VectorPayload(tag+".action",action)},
                {"seconds",J::Number(std::chrono::duration<double>(Clock::now()-started).count())}});
            if(norm<radius*(1.-1e-7) && attempt==0 && .5*radius>=lower)
            {
                radius*=.5;
                continue;
            }
            for(int i=0;i<5;++i)
            {
                auto t=evaluate(tag+".trial"+std::to_string(i),step,action,std::ldexp(1.,-i));
                if(t.accepted){shorter=std::move(t);break;}
            }
            break;
        }
    }
    const R8Trial *selected=bestRay;
    if(shorter.accepted && (!selected || shorter.E<selected->E))selected=&shorter;
    if(selected)
    {
        const auto audit=CompleteRestrictedAudit(selected->state,history,alpha,time,data);
        std::vector<std::pair<std::string,std::string>> fields;
        for(const auto &[key,value]:audit)fields.emplace_back(key,J::Numbers(value));
        durable.Event("r8_selected",{{"name",J::Quote(selected->name)},{"fraction",J::Number(selected->f)},
            {"residual_E",J::Number(selected->E)},{"ground",J::Number(selected->ground)},
            {"ground_ratio",J::Number(selected->ground/baseK)},{"audit",J::Object(fields)}});
    }
    durable.Event("r8_work",{{"native_trial_states",J::Number(evaluations)},
        {"retained_trial_states_reused","4"},{"selected",J::Quote(selected?selected->name:"none")}});
    return selected?*selected:R8Trial{};
}

RANSWallEvolutionPatch::R8Trial RANSWallEvolutionPatch::RunR9Comparison(const Vector &state,const Vector &history,
    double alpha,double time,const RANSWallEvolutionData &data,const ScalarDomainResult &domain,const Vector &residual,
    const std::map<int,double> &constraints,const WallMatrix &jac,const Vector &newton,
    const std::function<Vector(Vector,bool)> &inverse,RANSWallR3::DurableProgress &durable) const
{
    using Clock=std::chrono::steady_clock;
    const double us=m_robustUstar;
    const auto bounds=OuterBoundsR5(domain);
    const auto baseNorms=Norms(residual,constraints,us);
    Require(OuterViolationR5(bounds).Max()==0.,"R9 requires an in-box endpoint");
    auto at=[&](const Vector &d,double f=1.)
    {
        auto s=state;
        for(int i=0;i<Size();++i)s[i]+=f*d[i];
        return s;
    };
    auto difference=[&](const Vector &a,const Vector &b)
    {
        auto v=a;
        for(int i=0;i<Size();++i)v[i]-=b[i];
        return v;
    };
    auto scale=[](Vector v,double f)
    {
        for(double &x:v)x*=f;
        return v;
    };
    auto merit=[&](const Vector &norms)
    {
        double sum=0.;
        for(int i=0;i<InnerPressure();++i)sum+=norms[i]*norms[i];
        return .5*sum;
    };
    auto support=[&](const Vector &v)
    {
        double outside=0.;
        for(int i=0;i<Size();++i)
            if(i<offset[K()] || i>=offset[Epsilon()] || constraints.count(i))
                outside=std::max(outside,std::abs(v[i]));
        return outside;
    };
    int sourceValues=0,sourceActions=0,factorizations=0,refinements=0;
    double sourceSeconds=0.,actionSeconds=0.,factorSeconds=0.;
    auto source=[&](const Vector &q,InnerDestruction &value)
    {
        const auto started=Clock::now();++sourceValues;
        const auto checked=ScalarDomain(q,time,data,true);
        bool wall=false;
        if(checked.admissible)WallTraceJSON(checked.effectiveState,wall);
        if(wall)value=DestructionValueR9(checked,time,data,constraints);
        sourceSeconds+=std::chrono::duration<double>(Clock::now()-started).count();
        return wall;
    };
    auto sourceAction=[&](const Vector &q)
    {
        const auto started=Clock::now();++sourceActions;
        auto a=AssembleNonpressureR3(q,time,data,false,0,nullptr,true);
        actionSeconds+=std::chrono::duration<double>(Clock::now()-started).count();
        return a;
    };
    auto changes=[](const std::string &a,const std::string &b)
    {
        if(a.size()!=b.size())return -1;
        return static_cast<int>(std::inner_product(a.begin(),a.end(),b.begin(),0,std::plus<int>{},
            [](char x,char y){return x!=y;}));
    };
    InnerDestruction D0;
    Require(source(state,D0) && support(D0.value)==0.,"R9 base destruction support differs");
    const auto JD=sourceAction(state);
    double excluded=0.;
    for(int i=0;i<Size();++i)if(i<offset[K()] || i>=offset[Epsilon()] || constraints.count(i))
        excluded=std::max(excluded,Maximum(JD[i]));
    Require(excluded==0.,"R9 source action touches algebraic or overwritten rows");
    durable.Event("r9_component",{{"D0",durable.VectorPayload("r9.D0",D0.value)},
        {"JD0",durable.MatrixPayload(10000,JD)},{"source_only_matrix_label",J::Quote("JD0, separate from physical J")},
        {"native_b_m_lambda_chi_ranges_and_branch_counts",J::Numbers(D0.ranges)},
        {"smoothing",J::Quote(m_destruction.smooth?"C2":"Linear")},{"component_excluded_rows",J::Number(support(D0.value))},
        {"action_excluded_rows",J::Number(excluded)}});
    auto model=[&](const Vector &d,const InnerDestruction &value)
    {
        auto m=Multiply(jac,d);
        const auto action=Multiply(JD,d);
        for(int i=0;i<Size();++i)m[i]+=residual[i]+(value.value[i]-D0.value[i])-action[i];
        return m;
    };
    const auto zero=Norms(difference(model(Vector(Size()),D0),residual),constraints,us);
    Require(Maximum(zero)==0.,"R9 source model does not recover its base residual");
    const auto reference=ResidualResolutionComponents(state,history,alpha,time,data);
    auto inner=reference.terms.at("inner_k_destruction");
    for(const auto &[row,value]:constraints)inner[row]=0.;
    const auto baseDifference=Norms(difference(D0.value,inner),constraints,us);
    Require(Maximum(baseDifference)<1e-8,"R9 source component and residual decomposition differ");
    durable.Event("r9_base_consistency",{{"decomposition_error",J::Numbers(baseDifference)},
        {"model_at_zero_error",J::Numbers(zero)},
        {"first_order_model_identity",J::Quote("J + JD0 - JD0 = J on the selected base branch")}});

    auto directionalCheck=[&](const std::string &name,const Vector &d)
    {
        const auto exact=Multiply(JD,d);
        double error=1e300;
        for(int level:{9,11,13})
        {
            const double h=std::ldexp(1.,-level);
            std::array<InnerDestruction,4> values;
            bool stable=true;
            for(int i=0;i<4;++i)
            {
                const double f=(i%2?-1.:1.)*h/(i<2?1.:2.);
                if(!source(at(d,f),values[i])){stable=false;break;}
                stable=stable && values[i].limiter==D0.limiter && values[i].cap==D0.cap &&
                    values[i].profile==D0.profile && values[i].panels==D0.panels;
            }
            if(!stable)
            {
                durable.Event("r9_component_action_check",{{"name",J::Quote(name)},{"h",J::Number(h)},
                    {"same_ordered_branch_inventory","false"}});
                continue;
            }
            const auto coarse=scale(difference(values[0].value,values[1].value),.5/h);
            const auto fine=scale(difference(values[2].value,values[3].value),1/h);
            auto richardson=fine;
            for(int i=0;i<Size();++i)richardson[i]+=(fine[i]-coarse[i])/3.;
            const auto errors=Norms(difference(richardson,exact),constraints,us);
            const auto derivative=Norms(exact,constraints,us);
            const auto gap=Norms(difference(fine,coarse),constraints,us);
            durable.Event("r9_component_action_check",{{"name",J::Quote(name)},{"h",J::Number(h)},
                {"same_ordered_branch_inventory","true"},{"errors",J::Numbers(errors)},
                {"central_difference_gap",J::Numbers(gap)},{"analytic_norms",J::Numbers(derivative)}});
            error=std::min(error,Maximum(errors));
            if(error<=1e-7 && Maximum(gap)<=1e-6)break;
        }
        return error;
    };

    const auto controlStarted=Clock::now();
    R8Ray linear;
    std::string linearDetails;
    linear.step=CoupledDirectionR4(domain,newton,jac,constraints,inverse,linearDetails,bounds,1.,false,
        &durable,"r9_linear",residual,.125);
    linear.action=Multiply(jac,linear.step);
    const double linearConstructSeconds=std::chrono::duration<double>(Clock::now()-controlStarted).count();
    durable.Event("r9_linear_direction",{{"step",durable.VectorPayload("r9.linear.step",linear.step)},
        {"action",durable.VectorPayload("r9.linear.action",linear.action)},{"details",linearDetails},
        {"seconds",J::Number(linearConstructSeconds)}});
    auto p=linear.step;
    std::fill(p.begin(),p.begin()+offset[InnerPressure()],0.);
    const double linearError=directionalCheck("coupled_control",linear.step);
    const double pressureError=directionalCheck("control_pressure_part",p);
    const auto controlTrialStarted=Clock::now();
    std::map<double,R8Trial> controlSamples;
    auto control=RunR8Comparison(state,history,alpha,time,data,domain,residual,constraints,jac,newton,inverse,
        linear,durable,0,"r9_linear",&controlSamples);
    const double linearTrialSeconds=std::chrono::duration<double>(Clock::now()-controlTrialStarted).count();
    bool supportsModel=false;
    for(double f:{control.accepted?control.f:.0625,.125})
    {
        const auto &t=control.accepted && f==control.f?control:controlSamples.at(f);
        if(!t.qualified || t.residual.empty())continue;
        const auto d=scale(linear.step,f);
        InnerDestruction value;
        if(!source(t.state,value))continue;
        auto simple=residual;
        for(int i=0;i<Size();++i)simple[i]+=f*linear.action[i];
        const auto corrected=model(d,value);
        const auto L=Norms(difference(t.residual,simple),constraints,us);
        const auto M=Norms(difference(t.residual,corrected),constraints,us);
        const double algebraic=support(difference(difference(value.value,D0.value),Multiply(JD,d)));
        const auto decomposition=ResidualResolutionComponents(t.state,history,alpha,time,data);
        auto native=decomposition.terms.at("inner_k_destruction");
        for(const auto &[row,target]:constraints)native[row]=0.;
        const auto valueError=Norms(difference(native,value.value),constraints,us);
        Require(Maximum(valueError)<1e-8 && algebraic==0.,"R9 endpoint destruction component mismatch");
        supportsModel=supportsModel || (M[K()]<.5*L[K()] && Maximum(M)<.5*Maximum(L));
        durable.Event("r9_endpoint_model",{{"fraction",J::Number(f)},{"actual",J::Numbers(t.norms.empty()?Norms(t.residual,constraints,us):t.norms)},
            {"linear",J::Numbers(Norms(simple,constraints,us))},{"source_model",J::Numbers(Norms(corrected,constraints,us))},
            {"linear_vector_error",J::Numbers(L)},{"source_model_vector_error",J::Numbers(M)},
            {"decomposition_error",J::Numbers(valueError)},{"excluded_row_correction",J::Number(algebraic)},
            {"source_ranges",J::Numbers(value.ranges)},{"ordered_limiter_changes",J::Number(changes(D0.limiter,value.limiter))},
            {"ordered_cap_changes",J::Number(changes(D0.cap,value.cap))},{"ordered_profile_changes",J::Number(changes(D0.profile,value.profile))},
            {"ordered_panel_identity_changes",J::Number(changes(D0.panels,value.panels))}});
    }
    const bool actionPass=std::max(linearError,pressureError)<=1e-7;
    durable.Event("r9_model_gate",{{"action_pass",J::Boolean(actionPass)},{"endpoint_supports_model",J::Boolean(supportsModel)},
        {"control_trials_seconds",J::Number(linearTrialSeconds)}});
    if(!control.accepted || !supportsModel || !actionPass)
    {
        durable.Event("r9_decision",{{"selected",J::Quote(control.accepted?"linear":"none")},
            {"promote","false"},{"reason",J::Quote(!actionPass?"component_action_gate":!supportsModel?"endpoint_model_gate":"linear_trial_unresolved")},
            {"source_value_calls",J::Number(sourceValues)},{"source_action_calls",J::Number(sourceActions)},
            {"source_seconds",J::Number(sourceSeconds)},{"source_action_seconds",J::Number(actionSeconds)}});
        return control;
    }

    const auto modelStarted=Clock::now();
    auto anchor=linear.step;
    InnerDestruction anchorValue;
    if(!source(at(anchor),anchorValue))
    {
        anchor=scale(linear.step,control.f);
        Require(source(at(anchor),anchorValue),"R9 accepted control tangent is outside the scalar domain");
    }
    double incumbent=merit(Norms(model(anchor,anchorValue),constraints,us));
    Vector best;
    double bestMerit=1e300;
    for(int i=0;i<3;++i)
    {
        ++refinements;
        const auto ASource=sourceAction(at(anchor));
        auto A=jac;
        for(int row=offset[K()];row<offset[Epsilon()];++row)
            for(int col=0;col<Size();++col)A[row][col]+=ASource[row][col]-JD[row][col];
        const auto dSource=Multiply(ASource,anchor);
        auto r=residual;
        for(int row=offset[K()];row<offset[Epsilon()];++row)
            r[row]+=(anchorValue.value[row]-D0.value[row])-dSource[row];
        const auto factStarted=Clock::now();
        const auto factor=FeasibleFactor(A);++factorizations;
        factorSeconds+=std::chrono::duration<double>(Clock::now()-factStarted).count();
        auto apply=[&](Vector rhs,bool transpose){return factor.Refine(A,rhs,transpose);};
        const auto n=apply(scale(r,-1.),false);
        const std::string tag="r9_model"+std::to_string(i);
        const auto directionStarted=Clock::now();
        std::string details;
        Vector d;
        try
        {
            d=CoupledDirectionR4(domain,n,A,constraints,apply,details,bounds,1.,false,&durable,tag,r,.125);
        }
        catch(const std::exception &e)
        {
            durable.Event("r9_refinement_failure",{{"index",J::Number(i)},{"reason",J::Quote(e.what())},
                {"seconds",J::Number(std::chrono::duration<double>(Clock::now()-directionStarted).count())}});
            break;
        }
        const double directionSeconds=std::chrono::duration<double>(Clock::now()-directionStarted).count();
        const auto affine=Norms([&]{
            auto v=Multiply(A,d);
            for(int j=0;j<Size();++j)v[j]+=r[j];
            return v;
        }(),constraints,us);
        Require(*std::max_element(affine.begin()+InnerPressure(),affine.end())<1e-9,
            "R9 model direction loses frozen-base algebraic restoration");
        InnerDestruction value;
        const bool valid=source(at(d),value);
        const auto norms=valid?Norms(model(d,value),constraints,us):Vector{};
        const double objective=valid?merit(norms):1e300;
        durable.Event("r9_refinement",{{"index",J::Number(i)},{"source_tangent_state",durable.VectorPayload(tag+".anchor",at(anchor))},
            {"direction",durable.VectorPayload(tag+".step",d)},{"affine",J::Numbers(affine)},{"full_source_model",J::Numbers(norms)},
            {"merit",valid?J::Number(objective):"null"},{"previous_merit",J::Number(incumbent)},
            {"valid_source_state",J::Boolean(valid)},{"details",details},{"direction_seconds",J::Number(directionSeconds)}});
        if(valid && objective<bestMerit){best=d;bestMerit=objective;}
        if(valid && objective<incumbent)
        {
            const double previous=incumbent;
            incumbent=objective;anchor=std::move(d);anchorValue=std::move(value);
            if(previous-incumbent<=1e-5*std::max(1.,previous))break;
        }
        else
        {
            auto middle=anchor;
            for(int row=0;row<Size();++row)middle[row]=.5*(anchor[row]+d[row]);
            InnerDestruction candidate;
            if(!source(at(middle),candidate))break;
            anchor=std::move(middle);anchorValue=std::move(candidate);
            incumbent=merit(Norms(model(anchor,anchorValue),constraints,us));
        }
    }
    const double modelSeconds=std::chrono::duration<double>(Clock::now()-modelStarted).count();
    R8Trial nonlinear;
    double nonlinearTrialSeconds=0.;
    if(!best.empty())
    {
        const double error=directionalCheck("nonlinear_direction",best);
        if(error<=1e-7)
        {
            R8Ray ray;ray.step=best;ray.action=Multiply(jac,best);
            const auto trialStarted=Clock::now();
            nonlinear=RunR8Comparison(state,history,alpha,time,data,domain,residual,constraints,jac,newton,
                inverse,ray,durable,0,"r9_source");
            nonlinearTrialSeconds=std::chrono::duration<double>(Clock::now()-trialStarted).count();
            if(nonlinear.accepted)
            {
                InnerDestruction value;
                Require(source(nonlinear.state,value),"R9 qualified native trial lost source domain");
                const auto d=scale(best,nonlinear.f);
                const auto M=Norms(model(d,value),constraints,us);
                durable.Event("r9_source_trial_model",{{"fraction",J::Number(nonlinear.f)},{"model",J::Numbers(M)},
                    {"actual",J::Numbers(nonlinear.norms)},{"source_ranges",J::Numbers(value.ranges)}});
            }
        }
    }
    const double oldGain=Maximum(baseNorms)-control.E;
    const double newGain=nonlinear.accepted?Maximum(baseNorms)-nonlinear.E:0.;
    // Charge the source policy for its extra construction and trial work.
    const double controlCost=linearConstructSeconds+linearTrialSeconds;
    const double sourceCost=linearConstructSeconds+modelSeconds+nonlinearTrialSeconds;
    const bool promote=nonlinear.accepted && newGain>=1.25*oldGain &&
        newGain>=1e-4*Maximum(baseNorms) && nonlinear.ground>=control.ground &&
        nonlinear.norms.back()<=1.1*control.norms.back() && newGain/sourceCost>oldGain/controlCost;
    durable.Event("r9_decision",{{"selected",J::Quote(promote?"source":"linear")},{"promote",J::Boolean(promote)},
        {"reason",J::Quote(best.empty()?"bounded_model_construction_unresolved":!nonlinear.accepted?"source_trial_unresolved":
            promote?"material_matched_gain":"material_gain_gate_not_met")},
        {"linear_gain",J::Number(oldGain)},{"source_gain",J::Number(newGain)},
        {"linear_seconds",J::Number(controlCost)},{"source_seconds",J::Number(sourceCost)},
        {"refinements",J::Number(refinements)},{"source_value_calls",J::Number(sourceValues)},
        {"source_action_calls",J::Number(sourceActions)},{"source_value_seconds",J::Number(sourceSeconds)},
        {"source_action_seconds",J::Number(actionSeconds)},{"model_factorizations",J::Number(factorizations)},
        {"model_factorization_seconds",J::Number(factorSeconds)},{"model_seconds",J::Number(modelSeconds)}});
    return promote?nonlinear:control;
}

namespace
{
constexpr bool R10ModelAccept(double incumbent,double candidate,double slope,double eta,
    double nativeSlope,double bound,bool feasible)
{
    return feasible && slope<0. && nativeSlope<=bound &&
        candidate<=incumbent+1e-4*eta*slope && candidate<incumbent;
}
constexpr double R10Retain(double incumbent,double candidate,bool accepted)
{
    return accepted && candidate<incumbent?candidate:incumbent;
}
constexpr double R10AffineDefect(double C,double Jd){return C+Jd;}
static_assert(R10Retain(4.,10.,false)==4. && R10Retain(4.,3.,false)==4.);
static_assert(!R10ModelAccept(4.,5.,-1.,.5,-1.,-.01,true));
static_assert(R10Retain(4.,5.,R10ModelAccept(4.,5.,-1.,.5,-1.,-.01,true))==4.);
static_assert(R10AffineDefect(.2,.5*(-.2)+.5*(-.2))==0.);
static_assert(R10AffineDefect(.2,.125*(-.2))>.1);
static_assert(!R10ModelAccept(4.,3.,-1.,1.,.2361516653,-.000582328302,true));
}

RANSWallEvolutionPatch::R8Trial RANSWallEvolutionPatch::RunR10Comparison(
    const Vector &state,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,
    const ScalarDomainResult &domain,const Vector &residual,const std::map<int,double> &constraints,
    const WallMatrix &jac,const Vector &newton,const std::function<Vector(Vector,bool)> &inverse,
    const R8Ray &binding,RANSWallR3::DurableProgress &durable) const
{
    using Clock=std::chrono::steady_clock;
    const auto started=Clock::now();
    auto elapsed=[&]{return std::chrono::duration<double>(Clock::now()-started).count();};
    const double us=m_robustUstar;
    const auto baseNorms=Norms(residual,constraints,us);
    const double baseE=Maximum(baseNorms),slopeBound=-2e-4*baseE;
    const auto outer=OuterBoundsR5(domain);
    Require(OuterViolationR5(outer).Max()==0. && !m_destruction.smooth,
        "R10 requires the in-box frozen Linear-smoothing base");
    durable.Event("r10_budget",{{"tangent_limit","4"},{"armijo_limit_per_tangent","10"},
        {"new_source_values_limit","50"},{"native_limit_per_ray","9"},{"comparison_seconds","600"},
        {"source_stop_seconds","540"},{"native_selector_tests_full_endpoint","false"}});
    auto dot=[](const Vector &a,const Vector &b)
    {
        long double sum=0.;for(int i=0;i<a.size();++i)sum+=(long double)a[i]*b[i];
        return static_cast<double>(sum);
    };
    auto minus=[](Vector a,const Vector &b)
    {for(int i=0;i<a.size();++i)a[i]-=b[i];return a;};
    auto add=[](Vector a,const Vector &b,double f)
    {for(int i=0;i<a.size();++i)a[i]+=f*b[i];return a;};
    auto at=[&](const Vector &d){return add(state,d,1.);};
    auto objective=[&](const Vector &value)
    {
        const auto products=DualProductsR10(value,value,constraints);
        return .5*std::accumulate(products.begin(),products.begin()+InnerPressure(),0.);
    };
    const auto meritRows=MeritRowsR10(residual,jac,constraints);
    auto nativeSlope=[&](const Vector &d)
    {
        double value=-std::numeric_limits<double>::infinity();
        for(const auto &row:meritRows)value=std::max(value,dot(row,d));
        return value;
    };
    std::vector<DirectionRestriction> slopeCuts;
    for(auto row:meritRows)
    {
        for(double &x:row)x=-x;
        slopeCuts.push_back({std::move(row),-slopeBound,"native_merit_slope"});
    }
    auto algebraic=[&](const Vector &d)
    {
        const auto norms=Norms(add(residual,Multiply(jac,d),1.),constraints,us);
        return *std::max_element(norms.begin()+InnerPressure(),norms.end());
    };
    int valueCalls=0,cacheHits=0,actionCalls=0,solves=0;
    double valueSeconds=0.,actionSeconds=0.,solveSeconds=0.;
    struct CachedSource {bool valid=false;InnerDestruction source;};
    std::map<Vector,CachedSource> cache;
    auto source=[&](const Vector &d)->const CachedSource *
    {
        auto found=cache.find(d);
        if(found!=cache.end()){++cacheHits;return &found->second;}
        if(valueCalls>=50 || elapsed()>=540.)return nullptr;
        ++valueCalls;const auto begin=Clock::now();
        CachedSource item;
        const auto checked=ScalarDomain(at(d),time,data,true);
        if(checked.admissible)
        {
            bool wall=false;WallTraceJSON(checked.effectiveState,wall);
            if(wall){item.source=DestructionValueR9(checked,time,data,constraints);item.valid=true;}
        }
        valueSeconds+=std::chrono::duration<double>(Clock::now()-begin).count();
        return &cache.emplace(d,std::move(item)).first->second;
    };
    const Vector zero(Size());
    const auto *base=source(zero);
    Require(base && base->valid && base->source.value==binding.destruction &&
        binding.destructionAction.size()==Size(),"R10 bound source value differs from current extraction");
    const auto &D0=binding.destruction;
    const auto &JD=binding.destructionAction;
    auto model=[&](const Vector &d,const InnerDestruction &D)
    {return add(add(add(residual,Multiply(jac,d),1.),minus(D.value,D0),1.),Multiply(JD,d),-1.);};
    auto action=[&](const Vector &d)
    {
        ++actionCalls;const auto begin=Clock::now();
        auto value=AssembleNonpressureR3(at(d),time,data,false,0,nullptr,true);
        actionSeconds+=std::chrono::duration<double>(Clock::now()-begin).count();
        return value;
    };
    auto verified=[&](const Vector &d,const std::string &name)
    {
        std::string scan;
        const auto rows=CoupledDirectionR4(domain,d,jac,constraints,inverse,scan,outer,1.,false,
            &durable,name,residual,.125,nullptr,{},slopeCuts,true);
        const double affine=algebraic(d),slope=nativeSlope(d);
        const auto *value=source(d);
        const bool passed=!rows.empty() && affine<=5e-12 && slope<=slopeBound && value && value->valid;
        durable.Event("r10_anchor_feasibility",{{"name",J::Quote(name)},{"original_rows",scan},
            {"algebraic",J::Number(affine)},{"native_slope",J::Number(slope)},{"slope_bound",J::Number(slopeBound)},
            {"source_valid",J::Boolean(value && value->valid)},{"passed",J::Boolean(passed)}});
        return passed;
    };
    const double oldSlope=nativeSlope(binding.ascentRegression);
    durable.Event("r10_ascent_regression",{{"slope",J::Number(oldSlope)},{"bound",J::Number(slopeBound)},
        {"rejected_before_native_trials",J::Boolean(oldSlope>slopeBound)}});
    Require(std::abs(oldSlope-.2361516653)<1e-8 && oldSlope>slopeBound,"R10 retained ascent regression differs");

    const auto controlStarted=Clock::now();
    R8Ray controlRay;std::string details;
    controlRay.step=CoupledDirectionR4(domain,newton,jac,constraints,inverse,details,outer,1.,false,
        &durable,"r10_linear",residual,.125);
    controlRay.action=Multiply(jac,controlRay.step);
    Require(nativeSlope(controlRay.step)<=slopeBound,"R10 linear control fails native descent");
    durable.Event("r10_control",{{"direction",durable.VectorPayload("r10.linear.step",controlRay.step)},
        {"physical_action",durable.VectorPayload("r10.linear.action",controlRay.action)},
        {"native_slope",J::Number(nativeSlope(controlRay.step))},{"details",details}});
    const double directionSeconds=std::chrono::duration<double>(Clock::now()-controlStarted).count();
    const auto control=RunR8Comparison(state,history,alpha,time,data,domain,residual,constraints,jac,newton,
        inverse,controlRay,durable,0,"r10_linear");
    const double controlSeconds=std::chrono::duration<double>(Clock::now()-controlStarted).count();
    Require(control.accepted,"R10 matched native control unresolved");
    const auto sourcePolicyStarted=Clock::now();
    const double baseValueSeconds=valueSeconds;
    auto scaled=controlRay.step;for(double &x:scaled)x*=control.f;
    const double scaledDefect=algebraic(scaled);
    durable.Event("r10_scaled_anchor_regression",{{"fraction",J::Number(control.f)},
        {"actual_affine_defect",J::Number(scaledDefect)},
        {"expected_trace_defect",J::Number((1.-control.f)*baseNorms.back())},
        {"rejected_as_restored_anchor",J::Boolean(scaledDefect>5e-12)}});

    Vector anchor=controlRay.step,best=anchor;
    if(!verified(anchor,"initial"))
    {
        durable.Event("r10_decision",{{"promote","false"},{"reason",J::Quote("initial_model_anchor_unverified")},
            {"selected",J::Quote("linear")},{"source_values",J::Number(valueCalls)}});
        return control;
    }
    auto incumbent=model(anchor,source(anchor)->source);
    double phi=objective(incumbent),bestPhi=phi;
    const double initialPhi=phi;
    durable.Event("r10_anchor",{{"index","0"},{"direction",durable.VectorPayload("r10.anchor0",anchor)},
        {"Phi",J::Number(phi)},{"model",J::Numbers(Norms(incumbent,constraints,us))},
        {"native_slope",J::Number(nativeSlope(anchor))}});
    bool retryUsed=false,needSlopeCut=false,checkDerivative=true;
    std::string stop="tangent_budget";
    const auto modelStarted=Clock::now();
    for(int index=0;index<4 && elapsed()<500. && valueCalls<45;++index)
    {
        const auto *anchorSource=source(anchor);
        Require(anchorSource && anchorSource->valid,"R10 incumbent lost its cached source");
        WallMatrix JDanchor;
        try{JDanchor=action(anchor);}
        catch(const std::exception &error)
        {
            stop="moved_source_action_unresolved";
            durable.Event("r10_action_failure",{{"message",J::Quote(error.what())}});
            break;
        }
        auto A=jac;
        for(int row=0;row<Size();++row)for(int col=0;col<Size();++col)
            A[row][col]+=JDanchor[row][col]-JD[row][col];
        const auto r=add(model(anchor,anchorSource->source),Multiply(A,anchor),-1.);
        const auto consistency=Norms(minus(add(r,Multiply(A,anchor),1.),incumbent),constraints,us);
        Require(Maximum(consistency)<1e-10,"R10 affine tangent identity failed");
        const auto begin=Clock::now();
        const auto factor=FeasibleFactor(A);
        auto apply=[&](Vector rhs,bool transpose){return factor.Refine(A,rhs,transpose);};
        auto rhs=r;for(double &x:rhs)x=-x;
        const auto n=apply(rhs,false);
        Vector target;
        std::string subproblem;
        ++solves;
        try
        {
            target=CoupledDirectionR4(domain,n,A,constraints,apply,subproblem,outer,1.,false,
                &durable,"r10_tangent"+std::to_string(index),r,.125,nullptr,{},
                needSlopeCut?slopeCuts:std::vector<DirectionRestriction>{});
        }
        catch(const std::exception &error)
        {
            stop="tangent_construction_unresolved";
            durable.Event("r10_tangent_failure",{{"index",J::Number(index)},{"message",J::Quote(error.what())}});
            break;
        }
        const double seconds=std::chrono::duration<double>(Clock::now()-begin).count();solveSeconds+=seconds;
        const auto p=minus(target,anchor),Ap=Multiply(A,p);
        const auto cross=DualProductsR10(incumbent,Ap,constraints),square=DualProductsR10(Ap,Ap,constraints);
        const double slope=std::accumulate(cross.begin(),cross.begin()+InnerPressure(),0.);
        const double curvature=std::accumulate(square.begin(),square.begin()+InnerPressure(),0.);
        const double diagnostic=slope+.5*curvature,tolerance=1e-9*std::max({1.,phi,curvature});
        durable.Event("r10_tangent",{{"index",J::Number(index)},{"slope_cut",J::Boolean(needSlopeCut)},
            {"anchor",durable.VectorPayload("r10.t"+std::to_string(index)+".anchor",anchor)},
            {"target",durable.VectorPayload("r10.t"+std::to_string(index)+".target",target)},
            {"affine_identity_error",J::Numbers(consistency)},{"s_Phi",J::Number(slope)},
            {"Ap_squared",J::Number(curvature)},{"diagnostic",J::Number(diagnostic)},
            {"diagnostic_tolerance",J::Number(tolerance)},{"details",subproblem},{"seconds",J::Number(seconds)}});
        if(diagnostic>tolerance || slope>=0.)
        {stop="tangent_slope_consistency";break;}
        if(checkDerivative)
        {
            bool passed=false;
            for(int level:{11,13})
            {
                const double h=std::ldexp(1.,-level);
                std::array<const CachedSource *,4> values{};
                bool same=true;
                for(int i=0;i<4;++i)
                {
                    values[i]=source(add(anchor,p,(i%2?-1.:1.)*h/(i<2?1.:2.)));
                    if(!values[i] || !values[i]->valid){same=false;continue;}
                    const auto &v=values[i]->source,&a=anchorSource->source;
                    same=same && v.limiter==a.limiter && v.cap==a.cap && v.profile==a.profile && v.panels==a.panels;
                }
                const auto &reference=anchorSource->source;
                bool oneSided=!same && values[0] && values[2] && values[0]->valid && values[2]->valid;
                if(oneSided)for(int id:{0,2})
                {
                    const auto &v=values[id]->source;
                    oneSided=oneSided && v.limiter==reference.limiter && v.cap==reference.cap &&
                        v.profile==reference.profile && v.panels==reference.panels;
                }
                if(!same && !oneSided)
                {
                    durable.Event("r10_moved_action",{{"index",J::Number(index)},{"h",J::Number(h)},
                        {"ordered_branches_match","false"},{"passed","false"}});
                    continue;
                }
                auto coarse=minus(values[0]->source.value,oneSided?reference.value:values[1]->source.value);
                auto fine=minus(values[2]->source.value,oneSided?reference.value:values[3]->source.value);
                for(int row=0;row<Size();++row)
                {
                    coarse[row]/=oneSided?h:2*h;fine[row]/=oneSided?.5*h:h;
                    fine[row]+=(fine[row]-coarse[row])/(oneSided?1.:3.);
                }
                const auto error=Norms(minus(fine,Multiply(JDanchor,p)),constraints,us);
                passed=Maximum(error)<=1e-7;
                durable.Event("r10_moved_action",{{"index",J::Number(index)},{"h",J::Number(h)},
                    {"ordered_branches_match","true"},{"one_sided",J::Boolean(oneSided)},
                    {"error",J::Numbers(error)},{"passed",J::Boolean(passed)}});
                if(passed)break;
            }
            if(!passed){stop="moved_source_action_unresolved";break;}
        }
        bool accepted=false,slopeBlocked=false;
        for(int trial=0;trial<10 && valueCalls<50 && elapsed()<540.;++trial)
        {
            const double eta=std::ldexp(1.,-trial);
            const auto d=add(anchor,p,eta);
            const auto *value=source(d);
            if(!value)break;
            const auto valueM=value->valid?model(d,value->source):Vector{};
            const double nextPhi=value->valid?objective(valueM):phi;
            const double g=nativeSlope(d);
            const bool armijo=value->valid && nextPhi<=phi+1e-4*eta*slope && nextPhi<phi;
            slopeBlocked=slopeBlocked || (armijo && g>slopeBound);
            const bool possible=R10ModelAccept(phi,nextPhi,slope,eta,g,slopeBound,value->valid);
            const bool feasible=possible && verified(d,"t"+std::to_string(index)+".trial"+std::to_string(trial));
            const bool accept=possible && feasible;
            const std::string reason=!value->valid?"source_domain":!armijo?"model_armijo":
                g>slopeBound?"native_slope":!feasible?"original_rows_or_affine":"accepted";
            durable.Event("r10_model_trial",{{"index",J::Number(index)},{"eta",J::Number(eta)},
                {"Phi_before",J::Number(phi)},{"Phi_trial",value->valid?J::Number(nextPhi):"null"},
                {"Armijo_ceiling",J::Number(phi+1e-4*eta*slope)},{"native_slope",J::Number(g)},
                {"model_armijo",J::Boolean(armijo)},{"reason",J::Quote(reason)},{"accepted",J::Boolean(accept)}});
            if(!accept)continue;
            const auto &before=anchorSource->source,&after=value->source;
            checkDerivative=before.limiter!=after.limiter || before.cap!=after.cap ||
                before.profile!=after.profile || before.panels!=after.panels;
            anchor=d;incumbent=valueM;phi=nextPhi;
            if(phi<bestPhi){best=anchor;bestPhi=R10Retain(bestPhi,phi,true);}
            durable.Event("r10_anchor",{{"index",J::Number(index+1)},
                {"direction",durable.VectorPayload("r10.anchor"+std::to_string(index+1),anchor)},
                {"Phi",J::Number(phi)},{"best_Phi",J::Number(bestPhi)},{"eta",J::Number(eta)},
                {"model",J::Numbers(Norms(incumbent,constraints,us))},{"native_slope",J::Number(g)},
                {"branches_changed",J::Boolean(checkDerivative)}});
            accepted=true;break;
        }
        if(slopeBlocked && !retryUsed)
        {retryUsed=true;needSlopeCut=true;continue;}
        if(!accepted){stop=valueCalls>=50?"source_value_budget":"model_backtracking";break;}
    }
    const double modelSeconds=std::chrono::duration<double>(Clock::now()-modelStarted).count();
    if(elapsed()>=500.)stop="cooperative_time";
    else if(valueCalls>=45 && solves<4)stop="source_value_budget";
    durable.Event("r10_model_outcome",{{"stop",J::Quote(stop)},{"Phi_initial",J::Number(initialPhi)},
        {"Phi_best",J::Number(bestPhi)},{"solves",J::Number(solves)},{"retry_used",J::Boolean(retryUsed)},
        {"source_values",J::Number(valueCalls)},{"cache_hits",J::Number(cacheHits)},{"seconds",J::Number(modelSeconds)}});
    R8Trial nonlinear;
    double sourceTrialSeconds=0.;
    const bool distinct=best!=controlRay.step;
    if(distinct && elapsed()<565. && verified(best,"before_native"))
    {
        R8Ray ray;ray.step=best;ray.action=Multiply(jac,best);
        const auto trialStarted=Clock::now();
        nonlinear=RunR8Comparison(state,history,alpha,time,data,domain,residual,constraints,jac,newton,
            inverse,ray,durable,0,"r10_source");
        sourceTrialSeconds=std::chrono::duration<double>(Clock::now()-trialStarted).count();
    }
    else if(!distinct)nonlinear=control;
    const double controlGain=baseE-control.E,sourceGain=nonlinear.accepted?baseE-nonlinear.E:0.;
    const double sourceSeconds=directionSeconds+baseValueSeconds+
        std::chrono::duration<double>(Clock::now()-sourcePolicyStarted).count();
    const bool promote=distinct && nonlinear.accepted && sourceGain>=1.25*controlGain &&
        sourceGain>=1e-4*baseE && nonlinear.ground>=control.ground &&
        nonlinear.norms.back()<=1.1*control.norms.back() && sourceGain/sourceSeconds>controlGain/controlSeconds;
    durable.Event("r10_decision",{{"promote",J::Boolean(promote)},{"distinct",J::Boolean(distinct)},
        {"selected",J::Quote(promote?"source":"linear")},{"control_gain",J::Number(controlGain)},
        {"source_gain",J::Number(sourceGain)},{"control_seconds",J::Number(controlSeconds)},
        {"source_seconds",J::Number(sourceSeconds)},{"source_values",J::Number(valueCalls)},
        {"source_actions",J::Number(actionCalls)},{"source_value_seconds",J::Number(valueSeconds)},
        {"source_action_seconds",J::Number(actionSeconds)},{"tangent_seconds",J::Number(solveSeconds)},{"native_source_trial_seconds",J::Number(sourceTrialSeconds)},
        {"elapsed",J::Number(elapsed())},{"reason",J::Quote(!distinct?"no_distinct_improved_source_candidate":
            promote?"promotion_passed":!nonlinear.accepted?"native_trial_unresolved":"promotion_gain_or_cost")}});
    return promote?nonlinear:control;
}

namespace
{
constexpr double R11MixedRhs(double sampled,double anchorAction,double original,double fraction,bool freeDifferential)
{
    return freeDifferential?sampled/fraction-anchorAction:original;
}
static_assert(R11MixedRhs(2.,3.,.25,.125,false)==.25);
static_assert(R11MixedRhs(2.,3.,.25,.125,true)==13.);
static_assert(.25+.125*(-.25)==.21875);
constexpr double R11Square(double x){return x*x;}
static_assert((R11Square(2.+.125*(2.+.5))-R11Square(2.+.125*(2.-.5)))/(2*.125*.5)==4.5);
static_assert(R10Retain(4.,3.,false)==4.);
}

RANSWallEvolutionPatch::R8Trial RANSWallEvolutionPatch::RunR12Frozen(
    const Vector &state,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,
    const ScalarDomainResult &domain,const Vector &residual,const std::map<int,double> &constraints,
    const WallMatrix &jac,const Vector &newton,const std::function<Vector(Vector,bool)> &inverse,
    const R8Ray &binding,RANSWallR3::DurableProgress &durable) const
{
    const auto started=std::chrono::steady_clock::now();
    const auto outer=OuterBoundsR5(domain);
    const auto merit=MeritRowsR10(residual,jac,constraints);
    const double E=Maximum(Norms(residual,constraints,m_robustUstar)),slopeBound=-2e-4*E;
    std::vector<DirectionRestriction> cuts;
    for(auto row:merit){for(auto &x:row)x=-x;cuts.push_back({row,-slopeBound,"native_merit_slope"});}
    auto r=residual;
    for(int i=0;i<Size();++i)if(i<offset[InnerPressure()] && !constraints.count(i))r[i]/=.125;
    auto inspect=[&](const Vector &d,const std::string &name)
    {
        Require(d.size()==Size(),"R12 candidate binding missing");
        std::string scan;
        const auto passed=CoupledDirectionR4(domain,d,jac,constraints,inverse,scan,outer,1.,false,
            &durable,name,residual,.125,nullptr,{},cuts,true);
        auto jd=Multiply(jac,d),linear=jd,mixed=jd;
        for(int i=0;i<Size();++i){linear[i]+=residual[i];mixed[i]+=r[i];}
        const auto products=DualProductsR10(mixed,mixed,constraints);
        const double psi=.5*std::accumulate(products.begin(),products.begin()+InnerPressure(),0.);
        const auto norms=Norms(linear,constraints,m_robustUstar);
        const double C=*std::max_element(norms.begin()+InnerPressure(),norms.end());
        double slope=-std::numeric_limits<double>::infinity();
        for(const auto &row:merit)
        {long double s=0.;for(int i=0;i<Size();++i)s+=(long double)row[i]*d[i];slope=std::max(slope,double(s));}
        auto stateAt=state;for(int i=0;i<Size();++i)stateAt[i]+=d[i];
        const auto sourceDomain=ScalarDomain(stateAt,time,data,true);
        bool wall=false;if(sourceDomain.admissible)WallTraceJSON(sourceDomain.effectiveState,wall);
        const bool qualified=!passed.empty() && C<=5e-12 && slope<=slopeBound;
        durable.Event("r12_candidate",{{"name",J::Quote(name)},{"original_scan",scan},
            {"qualified_affine",J::Boolean(qualified)},{"full_anchor_source_domain",J::Boolean(sourceDomain.admissible && wall)},
            {"algebraic",J::Number(C)},{"native_slope",J::Number(slope)},{"slope_bound",J::Number(slopeBound)},
            {"Psi_Lf",J::Number(psi)},{"Phi_Lf",J::Number(psi*.125*.125)},
            {"direction",durable.VectorPayload("r12."+name+".direction",d)},
            {"mixed_residual",durable.VectorPayload("r12."+name+".mixed",mixed)}});
        Require(qualified,"R12 frozen candidate failed original feasibility/algebraic/slope checks");
        return psi;
    };
    const double upper=inspect(binding.r12Seed,"L1_witness");
    Require(std::abs(upper-272.9590401766774)<2e-9,"R12 feasible objective bound changed");
    if(m_r12Experiment=="Verify")
    {
        inspect(binding.r12Direction,"native_completed");
        inspect(binding.r12Reference,"independent_reference");
        durable.Event("r12_frozen_outcome",{{"status",J::Quote("independent_direction_original_oracle_verified")}});
        return {};
    }
    auto rhs=r;for(auto &x:rhs)x=-x;
    const auto n=inverse(rhs,false);
    std::string details;
    try
    {
        const auto d=CoupledDirectionR4(domain,n,jac,constraints,inverse,details,outer,1.,false,
            &durable,"r12_Lf",r,.125,nullptr,{},cuts);
        const auto objective=inspect(d,"Lf_completed");
        Require(objective<=upper+2e-9,"R12 result exceeds its feasible L1 upper bound");
        durable.Event("r12_frozen_outcome",{{"status",J::Quote("qualified")},{"details",details},
            {"cut_admission_limit",J::Number(m_r12CutLimit)},{"upper_Psi",J::Number(upper)},
            {"seconds",J::Number(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count())}});
    }
    catch(const std::exception &e)
    {
        durable.Event("r12_frozen_outcome",{{"status",J::Quote("unresolved")},{"reason",J::Quote(e.what())},
            {"cut_admission_limit",J::Number(m_r12CutLimit)},{"upper_Psi",J::Number(upper)},
            {"seconds",J::Number(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count())}});
    }
    return {};
}

RANSWallEvolutionPatch::R8Trial RANSWallEvolutionPatch::RunR11Comparison(
    const Vector &state,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,
    const ScalarDomainResult &domain,const Vector &residual,const std::map<int,double> &constraints,
    const WallMatrix &jac,const Vector &newton,const std::function<Vector(Vector,bool)> &inverse,
    const R8Ray &binding,RANSWallR3::DurableProgress &durable) const
{
    using Clock=std::chrono::steady_clock;
    const bool r12=ExperimentR12();
    const double us=m_robustUstar,referenceFraction=.125;
    const auto baseNorms=Norms(residual,constraints,us);
    const double baseE=Maximum(baseNorms),slopeBound=-2e-4*baseE;
    const auto outer=OuterBoundsR5(domain);
    const auto &D0=binding.destruction;
    const auto &JD=binding.destructionAction;
    Require(OuterViolationR5(outer).Max()==0. && !m_destruction.smooth && JD.size()==Size(),
        "R11 requires its frozen in-box Linear-smoothing source binding");
    auto dot=[](const Vector &a,const Vector &b)
    {long double sum=0.;for(int i=0;i<a.size();++i)sum+=(long double)a[i]*b[i];return double(sum);};
    auto add=[](Vector a,const Vector &b,double f)
    {for(int i=0;i<a.size();++i)a[i]+=f*b[i];return a;};
    auto scaled=[](Vector a,double f){for(auto &x:a)x*=f;return a;};
    auto at=[&](const Vector &d){return add(state,d,1.);};
    auto freeDifferential=[&](int row){return row<offset[InnerPressure()] && !constraints.count(row);};
    auto objective=[&](const Vector &r)
    {
        const auto v=DualProductsR10(r,r,constraints);
        return .5*std::accumulate(v.begin(),v.begin()+InnerPressure(),0.);
    };
    auto affine=[&](const Vector &d,double s){return add(residual,Multiply(jac,d),s);};
    auto model=[&](const Vector &d,const InnerDestruction &value)
    {return add(add(affine(d,1.),value.value,1.),add(D0,Multiply(JD,d),1.),-1.);};
    const auto rows=MeritRowsR10(residual,jac,constraints);
    auto nativeSlope=[&](const Vector &d)
    {
        double result=-std::numeric_limits<double>::infinity();
        for(const auto &r:rows)result=std::max(result,dot(r,d));
        return result;
    };
    std::vector<DirectionRestriction> slopeCuts;
    for(auto row:rows)
    {for(double &x:row)x=-x;slopeCuts.push_back({std::move(row),-slopeBound,"native_merit_slope"});}
    auto algebraic=[&](const Vector &d)
    {
        const auto v=Norms(affine(d,1.),constraints,us);
        return *std::max_element(v.begin()+InnerPressure(),v.end());
    };
    auto sourceValue=[&](const Vector &d,InnerDestruction &out)
    {
        const auto check=ScalarDomain(at(d),time,data,true);
        if(!check.admissible)return false;
        bool wall=false;WallTraceJSON(check.effectiveState,wall);
        if(wall)out=DestructionValueR9(check,time,data,constraints);
        return wall;
    };
    auto branchChanges=[](const InnerDestruction &a,const InnerDestruction &b)
    {
        auto count=[](const std::string &x,const std::string &y)
        {
            if(x.size()!=y.size())return -1;
            return int(std::inner_product(x.begin(),x.end(),y.begin(),0,std::plus<int>{},
                [](char a,char b){return a!=b;}));
        };
        return Vector{double(count(a.limiter,b.limiter)),double(count(a.cap,b.cap)),
            double(count(a.profile,b.profile)),double(count(a.panels,b.panels))};
    };
    const auto replayStarted=Clock::now();
    InnerDestruction baseSource;
    Require(sourceValue(Vector(Size()),baseSource) && baseSource.value==D0,"R11 base source binding differs");
    std::map<Vector,InnerDestruction> replayCache;
    bool replayPass=true;
    int replayValues=1;
    if(!r12)
    {
    for(const auto &item:binding.replay)
    {
        const auto d=scaled(item.direction,item.fraction);
        auto found=replayCache.find(d);
        if(found==replayCache.end())
        {
            InnerDestruction value;++replayValues;
            Require(sourceValue(d,value),"R11 retrospective sample is outside the source domain");
            found=replayCache.emplace(d,std::move(value)).first;
        }
        const auto M=model(d,found->second);
        const auto norms=Norms(M,constraints,us);
        double excluded=0.;
        const auto correction=add(add(found->second.value,D0,-1.),Multiply(JD,d),-1.);
        for(int row=0;row<Size();++row)if(!freeDifferential(row))excluded=std::max(excluded,std::abs(correction[row]));
        Require(excluded==0.,"R11 source correction touches excluded rows");
        std::vector<std::pair<std::string,std::string>> record{{"name",J::Quote(item.name)},
            {"fraction",J::Number(item.fraction)},{"anchor",J::Boolean(item.anchor)},
            {"model",J::Numbers(norms)},{"Phi_fraction",J::Number(objective(M))},
            {"source_ranges",J::Numbers(found->second.ranges)},{"ordered_branch_changes",J::Numbers(branchChanges(baseSource,found->second))},
            {"excluded_correction",J::Number(excluded)}};
        if(!item.actual.empty())
        {
            const auto actual=Norms(item.actual,constraints,us);
            const double error=std::abs(Maximum(actual)-Maximum(norms));
            const double gate=std::max(2e-6,.01*std::abs(baseE-Maximum(actual)));
            replayPass=replayPass && error<=gate;
            record.insert(record.end(),{{"actual",J::Numbers(actual)},{"Phi_actual",J::Number(objective(item.actual))},
                {"vector_error",J::Numbers(Norms(add(item.actual,M,-1.),constraints,us))},
                {"governing_merit_error",J::Number(error)},{"merit_error_gate",J::Number(gate)}});
        }
        durable.Event("r11_replay",record);
    }
    durable.Event("r11_replay_outcome",{{"passed",J::Boolean(replayPass)},{"source_values",J::Number(replayValues)},
        {"seconds",J::Number(std::chrono::duration<double>(Clock::now()-replayStarted).count())}});
    if(!replayPass || binding.replay.size()!=9)
    {
        durable.Event("r11_decision",{{"promote","false"},{"reason",J::Quote("retrospective_model_mismatch")},
            {"selected_policy","null"},{"retained_policy","null"}});
        return {};
    }
    }
    const auto started=Clock::now();
    auto elapsed=[&]{return std::chrono::duration<double>(Clock::now()-started).count();};
    durable.Event("r11_budget",{{"affine_solves_per_arm","1"},{"tangents_per_nonlinear_arm","3"},
        {"armijo_fractions_per_tangent","10"},{"source_values_per_nonlinear_arm","50"},{"native_trials_per_policy","9"},
        {"comparison_seconds","600"},{"reference_fraction","0.125"},{"promotion_comparators",J::Quote("Lf:L1; N1:L1; Nf:L1 and qualified Lf")},
        {"selection",J::Quote("highest gain/second promoted policy; otherwise lowest-E qualified trial")},
        {"shared_cost",J::Quote("physical J/factor shared; each nonlinear arm charged its seed construction and own work; cached native trials charged original evaluation cost")}});
    std::map<Vector,R8Trial> nativeCache;
    struct Policy
    {
        std::string name,status;
        Vector d;
        R8Trial trial;
        double directionSeconds=0.,cost=0.;
    };
    std::vector<Policy> policies;
    auto validate=[&](const Vector &d,const std::string &label)
    {
        std::string scan;
        const auto scanned=CoupledDirectionR4(domain,d,jac,constraints,inverse,scan,outer,1.,false,
            &durable,label,residual,.125,nullptr,{},slopeCuts,true);
        const auto check=ScalarDomain(at(d),time,data,true);
        bool wall=false;if(check.admissible)WallTraceJSON(check.effectiveState,wall);
        const double defect=algebraic(d),g=nativeSlope(d);
        const bool valid=!scanned.empty() && wall && defect<=5e-12 && g<=slopeBound;
        durable.Event("r11_anchor_check",{{"name",J::Quote(label)},{"scan",scan},{"algebraic",J::Number(defect)},
            {"native_slope",J::Number(g)},{"slope_bound",J::Number(slopeBound)},
            {"full_source_domain",J::Boolean(check.admissible && wall)},{"passed",J::Boolean(valid)}});
        return valid;
    };
    auto native=[&](Policy &p)
    {
        if(p.d.empty() || elapsed()>=565.)return;
        R8Ray ray;ray.step=p.d;ray.action=Multiply(jac,p.d);
        R11TrialWork work;
        const auto begin=Clock::now();
        p.trial=RunR8Comparison(state,history,alpha,time,data,domain,residual,constraints,jac,newton,
            inverse,ray,durable,0,p.name,nullptr,&nativeCache,&work);
        const double actual=std::chrono::duration<double>(Clock::now()-begin).count();
        p.cost+=actual+work.charged-work.actual;
        durable.Event("r11_native_cost",{{"name",J::Quote(p.name)},{"fresh",J::Number(work.fresh)},
            {"reused",J::Number(work.reused)},{"actual_seconds",J::Number(actual)},
            {"charged_seconds",J::Number(actual+work.charged-work.actual)}});
    };
    auto publish=[&](const Policy &p)
    {
        std::vector<std::pair<std::string,std::string>> r{{"name",J::Quote(p.name)},
            {"status",J::Quote(p.status)},{"cost",J::Number(p.cost)},
            {"direction_seconds",J::Number(p.directionSeconds)},{"accepted",J::Boolean(p.trial.accepted)}};
        if(!p.d.empty())r.insert(r.end(),{{"direction",durable.VectorPayload("r11."+p.name+".direction",p.d)},
            {"native_slope",J::Number(nativeSlope(p.d))}});
        if(p.trial.accepted)r.insert(r.end(),{{"selected_name",J::Quote(p.trial.name)},
            {"fraction",J::Number(p.trial.f)},{"E",J::Number(p.trial.E)},
            {"gain",J::Number(baseE-p.trial.E)},{"ground",J::Number(p.trial.ground)},
            {"trace_block",J::Number(p.trial.norms.back())},
            {"selected_state",durable.VectorPayload("r11."+p.name+".selected.state",p.trial.state)},
            {"selected_residual",durable.VectorPayload("r11."+p.name+".selected.residual",p.trial.residual)}});
        durable.Event("r11_policy",r);
        if(p.trial.accepted && (p.name=="L1" || p.name=="Lf"))
        {
            const auto prediction=affine(p.d,p.trial.f);
            durable.Event("r11_selected_prediction",{{"name",J::Quote(p.name)},{"fraction",J::Number(p.trial.f)},
                {"model",J::Numbers(Norms(prediction,constraints,us))},{"actual",J::Numbers(p.trial.norms)},
                {"vector_error",J::Numbers(Norms(add(p.trial.residual,prediction,-1.),constraints,us))}});
        }
    };
    for(const auto &name:{"L1","Lf"})
    {
        Policy p;p.name=name;
        const double s=p.name=="L1"?1.:referenceFraction;
        const auto begin=Clock::now();
        auto r=residual;
        for(int row=0;row<Size();++row)r[row]=R11MixedRhs(residual[row],0.,residual[row],s,freeDifferential(row));
        auto rhs=r;for(auto &x:rhs)x=-x;
        const auto n=s==1.?newton:inverse(rhs,false);
        std::string detail;
        if(r12)
        {
            p.d=p.name=="L1"?binding.r12Seed:binding.r12Direction;
            p.status=!p.d.empty() && validate(p.d,p.name+".initial")?"bound_qualified_affine":"bound_affine_invalid";
            p.directionSeconds=p.name=="L1"?binding.r12SeedCost:binding.r12SeedCost+binding.r12DirectionCost;
            p.directionSeconds+=std::chrono::duration<double>(Clock::now()-begin).count();
            p.cost=p.directionSeconds;
            if(p.status=="bound_qualified_affine")native(p);
            publish(p);policies.push_back(std::move(p));continue;
        }
        try
        {
            p.d=CoupledDirectionR4(domain,n,jac,constraints,inverse,detail,outer,1.,false,&durable,
                "r11_"+p.name,r,.125,nullptr,{},p.name=="Lf"?slopeCuts:std::vector<DirectionRestriction>{});
            p.status=validate(p.d,p.name+".initial")?"affine_constructed":"affine_anchor_invalid";
        }
        catch(const std::exception &e){p.status=std::string("affine_unresolved: ")+e.what();}
        p.directionSeconds=std::chrono::duration<double>(Clock::now()-begin).count();
        p.cost=p.directionSeconds;
        if(p.status=="affine_constructed")native(p);
        publish(p);policies.push_back(std::move(p));
    }
    for(int arm=0;arm<2;++arm)
    {
        Policy p;p.name=arm?"Nf":"N1";
        if(r12 && !arm)
        {p.status="historical_N1_not_repeated";publish(p);policies.push_back(std::move(p));continue;}
        const double s=arm?referenceFraction:1.;
        const auto modelStarted=Clock::now();
        int calls=0,hits=0,actions=0,solves=0,fullChecks=0;
        double actionSeconds=0.,sourceSeconds=0.;
        struct Value {bool valid=false;InnerDestruction D;};
        std::map<Vector,Value> cache;
        auto source=[&](const Vector &physicalDisplacement)->const Value *
        {
            const auto found=cache.find(physicalDisplacement);
            if(found!=cache.end()){++hits;return &found->second;}
            if(calls>=50 || elapsed()>=540.)return nullptr;
            ++calls;const auto begin=Clock::now();Value v;
            v.valid=sourceValue(physicalDisplacement,v.D);
            sourceSeconds+=std::chrono::duration<double>(Clock::now()-begin).count();
            return &cache.emplace(physicalDisplacement,std::move(v)).first->second;
        };
        const int seed=arm?1:0;
        double seedCost=policies[seed].directionSeconds;
        Vector anchor=policies[seed].d;
        ++fullChecks;
        if(anchor.empty() || !validate(anchor,p.name+".seed"))
        {
            if(arm)
            {
                anchor=policies[0].d;seedCost+=policies[0].directionSeconds;++fullChecks;
                durable.Event("r11_seed_fallback",{{"name",J::Quote(p.name)},{"fallback",J::Quote("L1 total displacement")}});
            }
            if(anchor.empty() || !validate(anchor,p.name+".fallback"))
            {p.status="no_verified_seed";p.cost=seedCost;publish(p);policies.push_back(std::move(p));continue;}
        }
        const auto *sample=source(scaled(anchor,s));
        if((!sample || !sample->valid) && arm && anchor!=policies[0].d && !policies[0].d.empty())
        {
            anchor=policies[0].d;seedCost+=policies[0].directionSeconds;++fullChecks;
            durable.Event("r11_seed_fallback",{{"name",J::Quote(p.name)},{"fallback",J::Quote("L1 after fractional sample rejection")}});
            if(validate(anchor,p.name+".sample_fallback"))sample=source(scaled(anchor,s));
            else sample=nullptr;
        }
        if(!sample || !sample->valid)
        {p.status="seed_sample_domain";p.cost=seedCost;publish(p);policies.push_back(std::move(p));continue;}
        auto M=model(scaled(anchor,s),sample->D);
        double psi=objective(M)/(s*s);
        Vector best=anchor;
        double bestPsi=psi;
        auto fractionalC=[&](const Vector &m)
        {
            const auto error=Norms(add(m,residual,-(1.-s)),constraints,us);
            const double result=*std::max_element(error.begin()+InnerPressure(),error.end());
            Require(result<=5e-12,"R11 physical prediction lost its fractional algebraic residual");
            return result;
        };
        bool checkAction=true;
        p.status="model_tangent_budget";
        auto normalized=[&](const Vector &m)
        {
            auto r=residual;
            for(int row=0;row<Size();++row)if(freeDifferential(row))r[row]=m[row]/s;
            return r;
        };
        durable.Event("r11_model_anchor",{{"name",J::Quote(p.name)},{"index","0"},
            {"direction",durable.VectorPayload("r11."+p.name+".anchor0",anchor)},
            {"fraction",J::Number(s)},{"Phi",J::Number(psi*s*s)},{"Psi",J::Number(psi)},
            {"native_prediction",J::Numbers(Norms(M,constraints,us))},{"native_slope",J::Number(nativeSlope(anchor))},
            {"fractional_C_identity_error",J::Number(fractionalC(M))}});
        for(int iteration=0;iteration<3 && elapsed()<500. && calls<44;++iteration)
        {
            const auto *base=source(scaled(anchor,s));
            if(!base || !base->valid){p.status="model_sample_domain";break;}
            ++actions;
            const auto actionStarted=Clock::now();
            WallMatrix sampleJD;
            try{sampleJD=AssembleNonpressureR3(at(scaled(anchor,s)),time,data,false,0,nullptr,true);}
            catch(const std::exception &e){p.status=std::string("source_action: ")+e.what();break;}
            actionSeconds+=std::chrono::duration<double>(Clock::now()-actionStarted).count();
            auto A=jac;
            double excluded=0.;
            for(int row=0;row<Size();++row)for(int col=0;col<Size();++col)
            {
                const double change=sampleJD[row][col]-JD[row][col];
                if(freeDifferential(row))A[row][col]+=change;
                else excluded=std::max(excluded,std::abs(change));
            }
            Require(excluded==0.,"R11 source action touches excluded rows");
            const auto Aa=Multiply(A,anchor);
            auto r=residual;
            for(int row=0;row<Size();++row)r[row]=R11MixedRhs(M[row],Aa[row],residual[row],s,freeDifferential(row));
            auto expected=normalized(M);
            const auto physicalAffine=affine(anchor,1.);
            for(int row=0;row<Size();++row)if(!freeDifferential(row))expected[row]=physicalAffine[row];
            const auto identity=Norms(add(add(r,Aa,1.),expected,-1.),constraints,us);
            Require(Maximum(identity)<1e-9 && algebraic(anchor)<=5e-12,"R11 mixed tangent/restoration identity failed");
            const auto solveStarted=Clock::now();
            const auto factor=FeasibleFactor(A);
            auto apply=[&](Vector rhs,bool transpose){return factor.Refine(A,rhs,transpose);};
            auto rhs=r;for(double &x:rhs)x=-x;
            const auto n=apply(rhs,false);
            ++solves;Vector target;std::string details;
            try{target=CoupledDirectionR4(domain,n,A,constraints,apply,details,outer,1.,false,&durable,
                "r11_"+p.name+".t"+std::to_string(iteration),r,.125,nullptr,{},slopeCuts);}
            catch(const std::exception &e){p.status=std::string("model_direction: ")+e.what();break;}
            const auto delta=add(target,anchor,-1.),Ap=Multiply(A,delta);
            const auto cross=DualProductsR10(normalized(M),Ap,constraints),square=DualProductsR10(Ap,Ap,constraints);
            const double slope=std::accumulate(cross.begin(),cross.begin()+InnerPressure(),0.);
            const double norm2=std::accumulate(square.begin(),square.begin()+InnerPressure(),0.);
            const double diagnostic=slope+.5*norm2;
            durable.Event("r11_tangent",{{"name",J::Quote(p.name)},{"index",J::Number(iteration)},
                {"source_sample",durable.VectorPayload("r11."+p.name+".sample"+std::to_string(iteration),at(scaled(anchor,s)))},
                {"s",J::Number(s)},{"identity_error",J::Numbers(identity)},{"excluded_source_action",J::Number(excluded)},
                {"slope_Psi",J::Number(slope)},{"norm2",J::Number(norm2)},{"diagnostic",J::Number(diagnostic)},
                {"target",durable.VectorPayload("r11."+p.name+".target"+std::to_string(iteration),target)},
                {"details",details},{"seconds",J::Number(std::chrono::duration<double>(Clock::now()-solveStarted).count())}});
            if(diagnostic>1e-9*std::max({1.,psi,norm2}) || slope>=0.)
            {p.status="model_slope_consistency";break;}
            if(checkAction)
            {
                bool passed=false;
                for(int level:{11,13})
                {
                    const double h=std::ldexp(1.,-level)/s;
                    std::array<const Value *,4> v{};bool same=true;
                    for(int i=0;i<4;++i)
                    {
                        v[i]=source(scaled(add(anchor,delta,(i%2?-1.:1.)*h/(i<2?1.:2.)),s));
                        if(!v[i] || !v[i]->valid){same=false;continue;}
                        same=same && v[i]->D.limiter==base->D.limiter && v[i]->D.cap==base->D.cap &&
                            v[i]->D.profile==base->D.profile && v[i]->D.panels==base->D.panels;
                    }
                    bool oneSided=!same && v[0] && v[2] && v[0]->valid && v[2]->valid;
                    if(oneSided)for(int i:{0,2})oneSided=oneSided && v[i]->D.limiter==base->D.limiter &&
                        v[i]->D.cap==base->D.cap && v[i]->D.profile==base->D.profile && v[i]->D.panels==base->D.panels;
                    if(!same && !oneSided)
                    {
                        durable.Event("r11_action",{{"name",J::Quote(p.name)},{"index",J::Number(iteration)},
                            {"optimizer_h",J::Number(h)},{"physical_h",J::Number(s*h)},{"ordered_branches_match","false"},{"passed","false"}});
                        continue;
                    }
                    auto coarse=add(v[0]->D.value,oneSided?base->D.value:v[1]->D.value,-1.);
                    auto fine=add(v[2]->D.value,oneSided?base->D.value:v[3]->D.value,-1.);
                    for(int row=0;row<Size();++row)
                    {
                        coarse[row]/=oneSided?s*h:2*s*h;fine[row]/=oneSided?.5*s*h:s*h;
                        fine[row]+=(fine[row]-coarse[row])/(oneSided?1.:3.);
                    }
                    const auto error=Norms(add(fine,Multiply(sampleJD,delta),-1.),constraints,us);
                    passed=Maximum(error)<=1e-7;
                    durable.Event("r11_action",{{"name",J::Quote(p.name)},{"index",J::Number(iteration)},
                        {"optimizer_h",J::Number(h)},{"physical_h",J::Number(s*h)},{"one_sided",J::Boolean(oneSided)},
                        {"error",J::Numbers(error)},{"passed",J::Boolean(passed)}});
                    if(passed)break;
                }
                if(!passed){p.status=calls>=50?"source_value_budget":"moved_action_unresolved";break;}
            }
            bool accepted=false;
            for(int trial=0;trial<10 && calls<48 && elapsed()<540.;++trial)
            {
                const double eta=std::ldexp(1.,-trial);
                const auto d=add(anchor,delta,eta);
                const auto *value=source(scaled(d,s));
                if(!value)break;
                const auto nextM=value->valid?model(scaled(d,s),value->D):Vector{};
                const double nextPsi=value->valid?objective(nextM)/(s*s):psi;
                const double g=nativeSlope(d);
                const bool armijo=value->valid && nextPsi<psi && nextPsi<=psi+1e-4*eta*slope;
                bool valid=false;
                if(armijo && g<=slopeBound){++fullChecks;valid=validate(d,p.name+".t"+std::to_string(iteration)+".a"+std::to_string(trial));}
                durable.Event("r11_model_trial",{{"name",J::Quote(p.name)},{"index",J::Number(iteration)},
                    {"eta",J::Number(eta)},{"Psi_before",J::Number(psi)},{"Psi_trial",value->valid?J::Number(nextPsi):"null"},
                    {"native_slope",J::Number(g)},{"accepted",J::Boolean(valid)},
                    {"reason",J::Quote(!value->valid?"sample_domain":!armijo?"model_armijo":g>slopeBound?"native_slope":!valid?"full_anchor":"accepted")}});
                if(!valid)continue;
                const auto changes=branchChanges(base->D,value->D);
                checkAction=std::any_of(changes.begin(),changes.end(),[](double x){return x!=0.;});
                anchor=d;M=nextM;psi=nextPsi;
                if(psi<bestPsi){best=anchor;bestPsi=R10Retain(bestPsi,psi,true);}
                durable.Event("r11_model_anchor",{{"name",J::Quote(p.name)},{"index",J::Number(iteration+1)},
                    {"eta",J::Number(eta)},{"direction",durable.VectorPayload("r11."+p.name+".anchor"+std::to_string(iteration+1),anchor)},
                    {"Phi",J::Number(psi*s*s)},{"Psi",J::Number(psi)},{"native_slope",J::Number(g)},
                    {"native_prediction",J::Numbers(Norms(M,constraints,us))},{"branch_changes",J::Numbers(changes)},
                    {"fractional_C_identity_error",J::Number(fractionalC(M))}});
                accepted=true;break;
            }
            if(!accepted){p.status=calls>=48?"source_value_budget":"model_backtracking";break;}
        }
        if(elapsed()>=500.)p.status="comparison_cooperative_limit";
        else if(solves<3 && calls>=44)p.status="source_value_budget";
        p.d=best;p.directionSeconds=seedCost+std::chrono::duration<double>(Clock::now()-modelStarted).count();
        p.cost=p.directionSeconds;
        ++fullChecks;
        if(elapsed()<565. && validate(best,p.name+".before_native"))native(p);
        const auto predictionStarted=Clock::now();
        if(p.trial.accepted)
        {
            const auto *value=source(scaled(best,p.trial.f));
            if(value && value->valid)
            {
                const auto prediction=model(scaled(best,p.trial.f),value->D);
                durable.Event("r11_selected_prediction",{{"name",J::Quote(p.name)},{"fraction",J::Number(p.trial.f)},
                    {"model",J::Numbers(Norms(prediction,constraints,us))},{"actual",J::Numbers(p.trial.norms)},
                    {"vector_error",J::Numbers(Norms(add(p.trial.residual,prediction,-1.),constraints,us))}});
            }
        }
        p.cost+=std::chrono::duration<double>(Clock::now()-predictionStarted).count();
        durable.Event("r11_model_work",{{"name",J::Quote(p.name)},{"source_values",J::Number(calls)},{"cache_hits",J::Number(hits)},
            {"source_actions",J::Number(actions)},{"tangent_solves",J::Number(solves)},{"full_anchor_checks",J::Number(fullChecks)},
            {"source_seconds",J::Number(sourceSeconds)},{"action_seconds",J::Number(actionSeconds)},
            {"elapsed_model_and_trials",J::Number(std::chrono::duration<double>(Clock::now()-modelStarted).count())}});
        publish(p);policies.push_back(std::move(p));
    }
    auto promote=[&](int candidate,int comparator)
    {
        const auto &p=policies[candidate],&c=policies[comparator];
        if(!p.trial.accepted || !c.trial.accepted)return false;
        const double gain=baseE-p.trial.E,baseGain=baseE-c.trial.E;
        const bool passed=gain>=1.25*baseGain && gain>=1e-4*baseE && p.trial.ground>=c.trial.ground &&
            p.trial.norms.back()<=1.1*c.trial.norms.back() && gain/p.cost>baseGain/c.cost;
        durable.Event("r11_promotion",{{"name",J::Quote(p.name)},{"comparator",J::Quote(c.name)},
            {"gain",J::Number(gain)},{"comparator_gain",J::Number(baseGain)},
            {"gain_gate",J::Boolean(gain>=1.25*baseGain)},{"absolute_gate",J::Boolean(gain>=1e-4*baseE)},
            {"ground_gate",J::Boolean(p.trial.ground>=c.trial.ground)},
            {"trace_gate",J::Boolean(p.trial.norms.back()<=1.1*c.trial.norms.back())},
            {"efficiency_gate",J::Boolean(gain/p.cost>baseGain/c.cost)},{"passed",J::Boolean(passed)}});
        return passed;
    };
    int chosen=-1;
    for(int i=1;i<4;++i)
    {
        bool passed=promote(i,0);
        if(i==3 && policies[1].trial.accepted)passed=promote(i,1) && passed;
        if(passed && (chosen<0 || (baseE-policies[i].trial.E)/policies[i].cost>
            (baseE-policies[chosen].trial.E)/policies[chosen].cost))chosen=i;
    }
    int retained=chosen;
    if(retained<0)for(int i=0;i<4;++i)
        if(policies[i].trial.accepted && (retained<0 || policies[i].trial.E<policies[retained].trial.E))retained=i;
    durable.Event("r11_decision",{{"promote",J::Boolean(chosen>=0)},
        {"selected_policy",chosen<0?"null":J::Quote(policies[chosen].name)},
        {"retained_policy",retained<0?"null":J::Quote(policies[retained].name)},
        {"elapsed_seconds",J::Number(elapsed())},{"distinct_native_states",J::Number(nativeCache.size())}});
    return retained<0?R8Trial{}:policies[retained].trial;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::SolveScalarFeasibleR3(
    Vector s,const Vector &history,double alpha,double time,const RANSWallEvolutionData &data,
    int &iterations,double tolerance,int maximumIterations,bool benchmarkOnly,bool localR4,bool localR5,
    const WallMatrix &frozenMatrix,const std::string &frozenAssembly,const R8Ray &ray) const
{
    Require(m_nonlinearPolicy=="ScalarFeasibleR3" && m_dimension==2 && !Mapped() &&
        m_transitionEnabled && m_cappedDestruction && m_comm->GetSize()==1,
        "ScalarFeasibleR3 requires its opt-in serial flat 2D Upstream+Capped policy");
    Require(Size()<=3000 && alpha==400. && time==.0025 && tolerance==2e-12 &&
        maximumIterations>0 && maximumIterations<=(JointRestorationR5()?32:16),"R3 requires the unchanged bounded exact target");
    Require(!localR4 || (CoupledR4() && maximumIterations==1),"R4 local experiment uses one bounded update");
    Require(!localR5 || (JointRestorationR5() && maximumIterations==1),
        "R5 local experiment requires the explicit bounded restoration policy");
    Require(frozenMatrix.empty() || ((DirectionReplayR6() || RefinedDirectionR6()) && localR5) || PilotR8(),
        "Frozen J is restricted to the bound comparisons and the first R8 pilot state");
    iterations=0;
    using Clock=std::chrono::steady_clock;
    const auto solveStarted=Clock::now();
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
    double adaptiveRadius=m_r7Radius;
    std::array<double,4> interior{};
    try
    {
        evaluate();
        Require(checkWall(domain),"R4 base wall polynomial fails the local trace safeguard");
        if(m_r7Interior)
        {
            const auto base=OuterViolationR5(OuterBoundsR5(domain));
            // Half-step entry needs m >= e for an affine inequality.
            // Add 25% headroom and keep the target fixed until box entry.
            for(int i=0;i<interior.size();++i)
                if(base.family[i]>0.)interior[i]=std::max(1e-8,1.25*base.family[i]);
            durable.Event("r7_fixed_interior_target",{{"margin",J::Numbers(interior)}});
        }
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
            if(PilotR8() && std::chrono::duration<double>(Clock::now()-solveStarted).count()>=300.)
            {
                terminal="r8_cooperative_completed_state_limit";
                throw std::runtime_error(terminal);
            }
            Require(iterations<maximumIterations,JointRestorationR5()?"R5 nonlinear update limit":"R3 Newton iteration limit");
            durable.JacobianStatus(iterations,"incomplete");
            std::vector<std::string> assembly;
            const bool useFrozen=!frozenMatrix.empty() && (!PilotR8() || iterations==0);
            const auto jac=useFrozen?frozenMatrix:
                timed("complete_action_assembly",[&]{return ScalarFeasibleJacobian(s,history,alpha,time,data,&assembly,false);});
            const auto matrix=durable.MatrixPayload(iterations,jac);
            durable.JacobianStatus(iterations,"complete",matrix,J::Object({
                {"shape",J::Numbers(Vector{double(Size()),double(Size())})},
                {"assembly",useFrozen?frozenAssembly:J::List(assembly)},
                {"loaded_bound_matrix",J::Boolean(useFrozen)},{"profiling",profiling()}}));
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
            if(ExperimentR12() && m_r12Experiment!="Compare")
            {
                Require(localR5 && useFrozen,"R12 requires its bound frozen problem");
                RunR12Frozen(s,history,alpha,time,data,domain,residual,constraints,jac,newton,inverse,ray,durable);
                m_lastFeasibleState=s;
                durable.Event("solve_terminal",{{"converged","false"},{"terminal",J::Quote("r12_frozen_completed")}});
                return s;
            }
            if(ExperimentR11() || ExperimentR12())
            {
                Require(localR5 && useFrozen,"R11 must use its bound current physical matrix");
                const auto trial=RunR11Comparison(s,history,alpha,time,data,domain,residual,constraints,
                    jac,newton,inverse,ray,durable);
                if(trial.accepted)
                {
                    durable.Snapshot(1,trial.state,trial.residual,trial.constraints,trial.norms,
                        domainReport(trial.domain),profiling(),J::Object({{"trial",J::Quote(trial.name)},
                            {"fraction",J::Number(trial.f)},{"physical_matrix",matrix},{"box",trial.box.report}}));
                    iterations=1;
                }
                m_lastFeasibleState=trial.accepted?trial.state:s;
                durable.Event("solve_terminal",{{"converged","false"},
                    {"terminal",J::Quote(ExperimentR12()?"r12_fixed_comparison_completed":"r11_fixed_comparison_completed")},{"iteration",J::Number(iterations)},
                    {"profiling",profiling()}});
                return m_lastFeasibleState;
            }
            if(ExperimentR10())
            {
                Require(localR5 && useFrozen,"R10 must use its bound current physical matrix");
                const auto trial=RunR10Comparison(s,history,alpha,time,data,domain,residual,constraints,
                    jac,newton,inverse,ray,durable);
                if(trial.accepted)
                {
                    durable.Snapshot(1,trial.state,trial.residual,trial.constraints,trial.norms,
                        domainReport(trial.domain),profiling(),J::Object({{"trial",J::Quote(trial.name)},
                            {"fraction",J::Number(trial.f)},{"physical_matrix",matrix},{"box",trial.box.report}}));
                    iterations=1;
                }
                m_lastFeasibleState=trial.accepted?trial.state:s;
                durable.Event("solve_terminal",{{"converged","false"},
                    {"terminal",J::Quote("r10_fixed_comparison_completed")},{"iteration",J::Number(iterations)},
                    {"profiling",profiling()}});
                return m_lastFeasibleState;
            }
            if(ExperimentR9())
            {
                Require(localR5 && !useFrozen,"R9 requires a fresh current endpoint matrix");
                auto trial=RunR9Comparison(s,history,alpha,time,data,domain,residual,constraints,jac,newton,inverse,durable);
                if(trial.accepted)
                {
                    durable.Snapshot(1,trial.state,trial.residual,trial.constraints,trial.norms,
                        domainReport(trial.domain),profiling(),J::Object({{"trial",J::Quote(trial.name)},
                            {"accepted_fraction",J::Number(trial.f)},{"matrix",matrix},{"box",trial.box.report}}));
                    iterations=1;
                }
                m_lastFeasibleState=trial.accepted?trial.state:s;
                durable.Event("solve_terminal",{{"converged","false"},
                    {"terminal",J::Quote("r9_fixed_comparison_completed")},{"iteration",J::Number(iterations)},
                    {"profiling",profiling()}});
                return m_lastFeasibleState;
            }
            if(ExperimentR8())
            {
                Require(localR5 && !frozenMatrix.empty(),"R8 comparison requires its bound current matrix");
                RunR8Comparison(s,history,alpha,time,data,domain,residual,constraints,jac,newton,inverse,ray,durable);
                m_lastFeasibleState=s;
                durable.Event("solve_terminal",{{"converged","false"},
                    {"terminal",J::Quote("r8_fixed_comparison_completed")},{"iteration",J::Number(iterations)},
                    {"profiling",profiling()}});
                return s;
            }
            if(PilotR8())
            {
                R8Ray currentRay;
                if(useFrozen)
                {
                    currentRay=ray;
                    Require(currentRay.residuals.size()==2 && currentRay.residuals.count(.125) &&
                        currentRay.residuals.count(.0625),"R8 pilot needs its two current first-state samples");
                }
                else
                {
                    std::string details;
                    currentRay.step=timed("r8_direction",[&]
                    {
                        return CoupledDirectionR4(domain,newton,jac,constraints,inverse,details,
                            OuterBoundsR5(domain),1.,false,&durable,"iterate"+std::to_string(iterations)+".r8_full",
                            residual,.125);
                    });
                    currentRay.action=Multiply(jac,currentRay.step);
                }
                auto trial=timed("r8_model_and_trials",[&]{
                    return RunR8Comparison(s,history,alpha,time,data,domain,residual,constraints,jac,newton,
                        inverse,currentRay,durable,iterations);
                });
                if(!trial.accepted)
                {
                    terminal="r8_nonlinear_trial_allocation_exhausted";
                    throw std::runtime_error(terminal);
                }
                const auto update=J::Object({{"iteration",J::Number(iterations)},{"accepted","true"},
                    {"accepted_fraction",J::Number(trial.f)},{"matrix",matrix},{"trial",J::Quote(trial.name)},
                    {"box",trial.box.report}});
                durable.Snapshot(iterations+1,trial.state,trial.residual,trial.constraints,trial.norms,
                    domainReport(trial.domain),profiling(),update);
                durable.Event("newton_update",{{"details",update},{"profiling",profiling()}});
                s=std::move(trial.state);residual=std::move(trial.residual);norms=std::move(trial.norms);
                constraints=std::move(trial.constraints);domain=std::move(trial.domain);alreadySaved=true;
                recordRanges(iterations+1,s);
                continue;
            }
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
                    double fraction=0,stepNorm=0,guardWork=0,modelAgreement=0;
                    std::string stepPin;
                };
                int trialCount=0;
                const std::array<double,7> scales{us,us,K0,E0,us*us,us*us,us/m_parameters.geometry.height};
                double algebraicRadiusLower=0;
                if(CompareR7() || InteriorProbeR7() || AdaptiveR7())
                    for(int i=0;i<Size();++i)
                        if(i>=offset[InnerPressure()] || constraints.count(i))
                        {
                            long double capacity=0;
                            for(int b=0;b<=Gauge();++b)
                                for(int j=offset[b];j<offset[b]+size[b];++j)
                                    capacity+=std::abs(jac[i][j])*scales[b];
                            if(capacity>0.)algebraicRadiusLower=std::max(algebraicRadiusLower,
                                static_cast<double>(std::abs(residual[i])/capacity));
                        }
                auto compare=[&](const std::string &name,bool joint,double restoration,bool inward,
                                 double radius=2.,int trialLimit=0,const std::array<double,4> &margin={})
                {
                    Trial result;
                    std::string details;
                    Vector step;
                    const auto started=Clock::now();++profile["r5_direction"].first;
                    try
                    {
                        Require(radius>=algebraicRadiusLower,"Radius below necessary algebraic restoration bound");
                        step=CoupledDirectionR4(domain,newton,jac,constraints,inverse,details,
                            joint?bounds:std::vector<OuterBound>{},restoration,inward,
                            (DirectionReplayR6() || RefinedDirectionR6())?&durable:nullptr,
                            "iterate"+std::to_string(iterations)+"."+name,residual,radius,
                            CompareR7()?&result.guardWork:nullptr,margin);
                    }
                    catch(const std::exception &e)
                    {
                        profile["r5_direction"].second+=std::chrono::duration<double>(Clock::now()-started).count();
                        details=e.what();
                        durable.Event("r5_direction_failure",{{"iteration",J::Number(iterations)},
                            {"name",J::Quote(name)},{"message",J::Quote(details)},{"radius",J::Number(radius)},
                            {"algebraic_radius_lower_bound",J::Number(algebraicRadiusLower)},{"profiling",profiling()}});
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
                    for(int b=0;b<=Gauge();++b)for(int i=offset[b];i<offset[b]+size[b];++i)
                        result.stepNorm=std::max(result.stepNorm,std::abs(step[i])/scales[b]);
                    durable.Event("r5_direction",{{"iteration",J::Number(iterations)},{"name",J::Quote(name)},
                        {"step",result.stepPin},{"action",durable.VectorPayload(prefix+".action",action)},
                        {"dimensionless_step",J::Number(result.stepNorm)},{"details",details},
                        {"guard_work_over_objective",J::Number(result.guardWork)},
                        {"algebraic_radius_lower_bound",J::Number(algebraicRadiusLower)},
                        {"model_block_norms",J::Numbers(full)},
                        {"model_box",OuterViolationR5(bounds,effective,trace,1.).report},
                        {"base_witness_actions",OuterViolationR5(bounds,effective,trace).report},
                        {"profiling",profiling()}});
                    const auto search=Clock::now();++profile["line_search"].first;
                    double fraction=1.;
                    // Reserve half the trial budget for partial restoration
                    // even if the full model constructs a valid direction.
                    const int allowance=trialLimit>0?trialLimit:RefinedDirectionR6()?
                        std::min(R6RestorationTrialLimit(localR5),2*R6RestorationTrialLimit(localR5)-trialCount):
                        (localR5?8:32-trialCount);
                    for(int k=0;k<allowance;++k,fraction*=.5)
                    {
                        ++trialCount;
                        auto trial=s;
                        for(int i=0;i<Size();++i)trial[i]+=fraction*step[i];
                        const auto label=prefix+".trial"+std::to_string(k);
                        auto checked=timed("scalar_domain",[&]{return ScalarDomain(trial,time,data,true);});
                        const auto modelNorms=model(fraction);
                        std::vector<std::pair<std::string,std::string>> event{
                            {"iteration",J::Number(iterations)},{"name",J::Quote(name)},{"backtrack",J::Number(k)},
                            {"fraction",J::Number(fraction)},{"state",durable.VectorPayload(label+".state",trial)},
                            {"scalar_domain",domainReport(checked)},{"model_block_norms",J::Numbers(modelNorms)},
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
                            const double required=localR5 && !(CompareR7() || InteriorProbeR7())?.01:1e-4*fraction;
                            bool boxDescent=candidateBox.Max()<=(1.-required)*baseBox.Max();
                            for(int family=0;family<4;++family)
                                boxDescent=boxDescent && candidateBox.family[family]<=baseBox.family[family];
                            const bool descent=Maximum(candidateNorms)<(1.-required)*error;
                            const bool accept=descent && (!joint || boxDescent);
                            const double predicted=error-Maximum(modelNorms);
                            const double ratio=predicted>0.?(error-Maximum(candidateNorms))/predicted:0.;
                            double traceNorm=0.;
                            for(int id:m_traceIds)traceNorm=std::max(traceNorm,std::abs(candidate[offset[Epsilon()]+id])/E0);
                            event.insert(event.end(),{{"classification",J::Quote(accept?"accepted_comparison":
                                !descent?"residual_non_descent":"outer_non_descent_or_family_worsening")},
                                {"box",candidateBox.report},{"scaled_block_norms",J::Numbers(candidateNorms)},
                                {"trace_constraint_norm",J::Number(candidateNorms.back())},
                                {"actual_matching_trace_norm",J::Number(traceNorm)},
                                {"required_relative_decrease",J::Number(required)},
                                {"model_agreement",predicted>0.?J::Number(ratio):"null"},
                                {"residual",durable.VectorPayload(label+".residual",candidate)}});
                            durable.Event("r5_trial",event);
                            if(accept || CompareR7() || InteriorProbeR7())
                            {
                                const auto audit=CompleteRestrictedAudit(trial,history,alpha,time,data,true);
                                durable.Event("r5_trial_audit",{{"iteration",J::Number(iterations)},
                                    {"name",J::Quote(name)},{"fraction",J::Number(fraction)},
                                    {"covariance",J::Numbers(audit.at("r2_covariance"))},
                                    {"qualified",J::Boolean(audit.at("r2_qualified")[0]==1.)}});
                            }
                            if(accept)
                            {
                                result.modelAgreement=ratio;
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
                    if(CompareR7() || InteriorProbeR7())
                        durable.Event("r7_local_trial_outcome",{{"name",J::Quote(name)},{"fraction",J::Number(result.fraction)},
                            {"trial_allocation",J::Number(allowance)},
                            {"allocation_exhausted",J::Boolean(result.fraction==0.)},{"profiling",profiling()}});
                    return std::pair{std::move(result),details};
                };
                if(CompareR7())
                {
                    Require(localR5,"R7 comparison cannot advance nonlinear iterates");
                    const auto baseline=compare("r7_r6_radius2",true,1.,true,2.,8).first;
                    Require(baseline.stepNorm>0.,"R7 could not measure the radius-2 direction");
                    // Twice the last accepted displacement is still close
                    // to the observed finite-step regime at the R6 endpoint.
                    double radius=std::min(.5*baseline.stepNorm,2*m_r7PreviousDisplacement);
                    durable.Event("r7_initial_radius",{{"baseline_norm",J::Number(baseline.stepNorm)},
                        {"last_accepted_displacement",J::Number(m_r7PreviousDisplacement)},
                        {"radius",J::Number(radius)}});
                    double restoration=1.;
                    auto contracted=compare("r7_contracted",true,restoration,true,radius,4).first;
                    if(contracted.stepPin.empty() && baseBox.Max()>0.)
                    {
                        restoration=.25;
                        contracted=compare("r7_contracted_partial",true,restoration,true,radius,4).first;
                    }
                    if(contracted.stepPin.empty() || contracted.stepNorm<radius*(1.-1e-7))
                    {
                        radius=contracted.stepPin.empty()?std::min(2*radius,.75*baseline.stepNorm):.5*radius;
                        restoration=1.;
                        contracted=compare("r7_adjusted_radius",true,restoration,true,radius,4).first;
                        if(contracted.stepPin.empty() && baseBox.Max()>0.)
                        {
                            restoration=.25;
                            contracted=compare("r7_adjusted_partial",true,restoration,true,radius,4).first;
                        }
                    }
                    // A material multiplier only triggers the matched
                    // control. Trial outcomes decide whether to use it.
                    if(!contracted.stepPin.empty() && contracted.guardWork>.01)
                        compare("r7_positive_ground_control",true,restoration,false,radius,8);
                    m_lastFeasibleState=s;
                    durable.Event("solve_terminal",{{"converged","false"},
                        {"terminal",J::Quote("r7_local_comparison_completed")},{"iteration",J::Number(iterations)},
                        {"profiling",profiling()}});
                    return s;
                }
                if(InteriorProbeR7())
                {
                    const auto trial=compare("r7_interior_probe",true,1.,false,m_r7Radius,8,interior).first;
                    m_lastFeasibleState=s;
                    durable.Event("solve_terminal",{{"converged","false"},
                        {"terminal",J::Quote("r7_interior_probe_completed")},{"iteration",J::Number(iterations)},
                        {"useful_joint_descent",J::Boolean(trial.fraction>0.)},{"profiling",profiling()}});
                    return s;
                }
                if(AdaptiveR7())
                {
                    const auto margin=baseBox.Max()>0.?interior:std::array<double,4>{};
                    const double startRadius=std::min(2.,std::max(adaptiveRadius,1.05*algebraicRadiusLower));
                    double radius=startRadius;
                    auto trial=compare("r7_full",true,1.,false,radius,8,margin).first;
                    bool constructed=!trial.stepPin.empty();
                    if(trial.fraction==0.)
                    {
                        const double retry=constructed?std::max({.125,.5*radius,1.05*algebraicRadiusLower}):
                            std::min(2.,2.*radius);
                        if(retry!=radius && retry<=2.)
                        {
                            radius=retry;
                            trial=compare("r7_retry_full",true,1.,false,radius,8,margin).first;
                            constructed=constructed || !trial.stepPin.empty();
                        }
                    }
                    if(trial.fraction==0. && baseBox.Max()>0.)
                    {
                        radius=std::max(radius,startRadius);
                        trial=compare("r7_partial",true,.25,false,radius,16).first;
                        constructed=constructed || !trial.stepPin.empty();
                    }
                    if(trial.fraction==0.)
                    {
                        terminal=constructed?"r7_nonlinear_trial_allocation_exhausted":"r7_direction_construction_unresolved";
                        throw std::runtime_error(terminal);
                    }
                    adaptiveRadius=radius;
                    if(trial.fraction==1. && trial.modelAgreement>=.75)adaptiveRadius=std::min(2.,1.5*radius);
                    else if(trial.fraction<=.25 || trial.modelAgreement<.25)adaptiveRadius=std::max(.125,.5*radius);
                    const auto update=J::Object({{"iteration",J::Number(iterations)},{"accepted","true"},
                        {"accepted_fraction",J::Number(trial.fraction)},{"matrix",matrix},{"step",trial.stepPin},
                        {"trial_count",J::Number(trialCount)},{"box",trial.box.report},{"radius",J::Number(radius)},
                        {"next_radius",J::Number(adaptiveRadius)},{"model_agreement",J::Number(trial.modelAgreement)}});
                    durable.Snapshot(iterations+1,trial.state,trial.residual,trial.constraints,trial.norms,
                        domainReport(trial.domain),profiling(),update);
                    durable.Event("newton_update",{{"details",update},{"profiling",profiling()}});
                    s=std::move(trial.state);residual=std::move(trial.residual);norms=std::move(trial.norms);
                    constraints=std::move(trial.constraints);domain=std::move(trial.domain);alreadySaved=true;
                    recordRanges(iterations+1,s);
                    continue;
                }
                if(localR5 && !(DirectionReplayR6() || RefinedDirectionR6()))compare("r4_control",false,0.,true);
                auto [trial,details]=compare("r5_full",true,1.,true);
                if(!RefinedDirectionR6() && trial.fraction==0. && baseBox.Max()>0. &&
                    details.find("inward_ground_row")!=std::string::npos)
                    std::tie(trial,details)=compare("r5_inward_control",true,1.,false);
                const auto fullOutcome=trial.fraction>0.?FullRestorationOutcome::Accepted:
                    trial.stepPin.empty()?FullRestorationOutcome::DirectionUnresolved:
                        FullRestorationOutcome::LineSearchExhausted;
                const bool tryPartial=RefinedDirectionR6()?R6TryPartialRestoration(fullOutcome,baseBox.Max()):
                    trial.fraction==0. && baseBox.Max()>0. && trial.stepPin.empty();
                if(tryPartial)
                {
                    if(RefinedDirectionR6())
                        durable.Event("r6_partial_fallback",{{"iteration",J::Number(iterations)},
                            {"full_direction_constructed",J::Boolean(!trial.stepPin.empty())},
                            {"full_trials",J::Number(trialCount)},
                            {"reserved_trials",J::Number(R6RestorationTrialLimit(localR5))}});
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
