// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <cmath>
#include <algorithm>

namespace Nektar
{
using namespace SolverUtils;
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;

std::string RANSWallEvolutionPatch::SourceBalanceDiagnostic(const Vector &s,double time,
    const RANSWallEvolutionData &data,const Vector &stations,const Vector &heights) const
{
    Require(!NewModels(),"CLS-specific diagnostic/action unsupported for A93/A95/E94E95; use native model-aware audits");
    Require(!AcceptedAWF(),"legacy chi/cap source diagnostic is unsupported for the accepted wall formulation; use accepted_awf readback");
    Require(Dimension()==2 && !Mapped(),"source diagnostic is restricted to the flat 2D pilot");
    const auto hash=RANSWallR3::Hash(s);const auto domain=ScalarDomain(s,time,data);
    Require(domain.admissible,"selected source readback is outside its declared scalar inventory: "+domain.failure);
    auto resolution=ResidualResolutionComponents(s,s,0.,time,data,false);
    resolution.native["spatial_unmasked"]=Residual(s,s,0.,time,data,false,false);
    const auto rate=RANSWallResolutionMassRate(*this,resolution,K(),{1.6,0,0});
    double componentLoadDefect=0.,componentRateDefect=0.;
    for(int i=0;i<size[K()];++i)
    {
        double load=0.,sum=0.;
        for(const auto &name:RANSWallResolutionSpatialTerms())load+=resolution.terms.at(name)[offset[K()]+i];
        for(const auto &[name,value]:rate.componentRates)sum+=value[i];
        componentLoadDefect=std::max(componentLoadDefect,std::abs(load-resolution.native.at("spatial_unmasked")[offset[K()]+i]));
        componentRateDefect=std::max(componentRateDefect,std::abs(sum-rate.rate[i]));
    }
    Require(componentLoadDefect<1e-9 && componentRateDefect<1e-7,"source component/native sum mismatch");
    auto tensor=[](const RANSTensor &g){Vector v;for(const auto &row:g)for(double x:row)v.push_back(x);return v;};
    auto checkedNormalization=[&](const RANSWallCoordinate &p,const EquilibriumAWFProfiles &profile)
    {
        // Safe for read-only extra queries: no closure/normalization at an
        // invalid transported k. Do not add constraints or repair this state.
        for(const auto &row:NormalizationRows(p,profile))
        {
            const double k=row.Value(s,offset[K()]);
            if(!std::isfinite(k) || k<0)
            {
                m_solveInvalidContext=J::Object({{"query",J::Quote("source_balance_ray")},{"coordinate",J::Numbers(row.coordinate)},
                    {"geometry",J::Number(row.geometry)},{"k",J::Number(k)},{"anchor",J::Numbers(p)},
                    {"chi",J::Number(profile.Reference().weight)}});
                throw std::runtime_error("extra source-query ray has inadmissible k");
            }
        }
        return WallNormalization(s,p,profile);
    };
    std::vector<std::string> points,budgets,weak;
    int normalizations=0,massApplications=rate.componentRates.size()+1;
    for(double x:stations)
    {
        const RANSWallCoordinate anchor{x,0,0};
        // Only profile/sampler construction is needed, never Lift. This
        // branch has its own invalid-query result, separate from production.
        const auto query=ScalarDomain(s,time,data,false,true,{anchor});
        if(!query.admissible)
        {
            points.push_back(J::Object({{"x",J::Number(x)},{"valid","false"},{"reason",J::Quote(query.failure)},
                {"witness",J::Numbers(query.limitingPoint)},{"inventory",J::Quote(query.limitingInventory)}}));continue;
        }
        const auto selected=StateProfile(s,anchor,time,data);const auto gates=Transition(s,anchor,time,data);
        for(int variant:{-1,0,1})
        {
            m_solveInvalidContext.clear();
            try
            {
                const auto profile=variant<0?selected:Profile(SampleK(s,anchor),anchor,variant);
                if(variant==1)Require(profile->Reference().weight==1 && profile->Reference().reason==AWFFallbackReason::Active,
                    "fixed-state chi=1 profile is intrinsically inadmissible");
                const auto normal=checkedNormalization(anchor,*profile);++normalizations;
                double sumD=0,sumPre=0;
                for(const auto &q:profile->Quadrature())
                {
                    auto p=anchor;p[2]=q.r;const double k=Point(K(),p).Value(s,offset[K()]);
                    const auto pre=profile->AtProduction(q.r,0.);const double d=pre.dissipation-pre.redistribution;
                    const auto value=EvaluateRANSCappedDestructionPoint(k,profile->K(),pre.dissipation,d,normal.b,normal.m,m_destruction);
                    sumD+=q.weight*q.jacobian*value.destruction;sumPre+=q.weight*q.jacobian*d;
                }
                budgets.push_back(J::Object({{"x",J::Number(x)},{"variant",J::Number(variant)},
                    {"chi",J::Number(profile->Reference().weight)},{"sum_w_D",J::Number(sumD)},
                    {"lambda_sum_w_Dpre",J::Number(normal.lambda*sumPre)},
                    {"budget_identity_defect",J::Number(sumD-normal.lambda*sumPre)},
                    {"B",J::Number(normal.b)},{"m",J::Number(normal.m)},{"lambda",J::Number(normal.lambda)},
                    {"quadrature_count",J::Number(profile->Quadrature().size())}}));
                const double ze=profile->Original().DissipationCutoff();Vector zs=heights;
                for(double f:{.25,.5,1.,2.,5.})zs.push_back(f*ze);
                std::sort(zs.begin(),zs.end());zs.erase(std::unique(zs.begin(),zs.end()),zs.end());
                for(double z:zs)if(z>=0 && z<=m_parameters.geometry.height)
                {
                    RANSWallCoordinate p{x,0,z};const auto row=Point(K(),p);const double k=row.Value(s,offset[K()]);
                    if(!std::isfinite(k) || k<0)
                    {points.push_back(J::Object({{"x",J::Number(x)},{"z",J::Number(z)},{"variant",J::Number(variant)},
                        {"valid","false"},{"k",J::Number(k)},{"reason",J::Quote("extra point has inadmissible k")}}));continue;}
                    RANSTensor g{},reference{};for(int a=0;a<m_dimension;++a)g[PhysicalAxis(a)]=Point(a,p).Gradient(s,offset[a]);
                    const auto ag=profile->Original().VelocityGradient(z,{profile->Reference().matchingVelocity,0},{0,0});
                    for(int a=0;a<3;++a)for(int b=0;b<3;++b)reference[a][b]=m_referenceFlow[a]*m_referenceNormal[b]*ag[0];
                    Vector velocity;for(int a=0;a<m_dimension;++a)velocity.push_back(Point(a,p).Value(s,offset[a]));
                    const double nt=profile->Original().EddyViscosity(z);const auto actual=InnerStress(k,nt,g);
                    const auto linear=EvaluateRANSWallInnerStress(k,nt,g,RANSWallInnerStressModel::Linear);
                    const auto analyticalLinear=EvaluateRANSWallInnerStress(k,nt,reference,RANSWallInnerStressModel::Linear);
                    const auto analyticalSelected=InnerStress(k,nt,reference);
                    const auto pre=profile->AtProduction(z,actual.production);const double dpre=pre.dissipation-pre.redistribution;
                    const double ea=profile->Original().OriginalDissipation(z),p0=profile->ReferenceProduction(z),chi=profile->Reference().weight;
                    const auto value=EvaluateRANSCappedDestructionPoint(k,profile->K(),pre.dissipation,dpre,normal.b,normal.m,m_destruction);
                    std::vector<std::pair<std::string,std::string>> components;
                    for(const auto &[name,v]:rate.componentRates)components.push_back({name,J::Number(row.Value(v))});
                    points.push_back(J::Object({{"x",J::Number(x)},{"z",J::Number(z)},{"valid","true"},{"variant",J::Number(variant)},
                        {"chi",J::Number(chi)},{"spatial_gate",J::Number(gates.spatial)},{"state_gate",J::Number(gates.state)},
                        {"admissibility_gate",J::Number(gates.admissibility)},{"margin_gate",J::Number(gates.marginGate)},
                        {"reynolds_gate",J::Number(gates.reynoldsGate)},{"sensor_d2",J::Number(gates.d2)},
                        {"local_k",J::Number(k)},{"sampled_K",J::Number(profile->K())},
                        {"K_sampling",J::Quote(m_kSampling)},{"sample_height",ColumnAverageK()?"null":J::Number(m_parameters.geometry.sampleHeight)},
                        {"z_epsilon",J::Number(ze)},{"c_length",J::Number(profile->Original().Parameters().cLength)},
                        {"molecular_nu",J::Number(m_parameters.geometry.molecularNu)},
                        {"represented_velocity",J::Numbers(velocity)},{"represented_gradient",J::Numbers(tensor(g))},{"analytical_gradient",J::Numbers(tensor(reference))},
                        {"production",J::Number(actual.production)},{"linear_production_same_gradient",J::Number(linear.production)},
                        {"analytical_gradient_linear_production",J::Number(analyticalLinear.production)},
                        {"analytical_gradient_saturated_production",J::Number(analyticalSelected.production)},
                        {"reference_production",J::Number(p0)},{"epsilon_A",J::Number(ea)},{"D_pre",J::Number(dpre)},
                        {"D_pre_identity_defect",J::Number(dpre-((1-chi)*ea+chi*p0))},
                        {"limiter",J::Number(value.limiter)},{"B",J::Number(normal.b)},{"m",J::Number(normal.m)},{"lambda",J::Number(normal.lambda)},
                        {"dependent_epsilon",J::Number(value.dissipation)},{"redistribution",J::Number(value.redistribution)},
                        {"D",J::Number(value.destruction)},{"D_over_k",J::Number(value.rate)},{"source",J::Number(actual.production-value.destruction)},
                        {"auxiliary_nu",J::Number(nt)},{"k_diffusivity",J::Number(m_parameters.geometry.molecularNu+nt/m_parameters.model.Parameters().sigmaK)},
                        {"saturation_active",J::Boolean(actual.active)},{"saturation_a_c",J::Numbers(Vector{actual.a,actual.c})},
                        {"selected_native_spatial_rate",J::Number(row.Value(rate.rate))},{"selected_native_component_rates",J::Object(components)}}));
                }
            }
            catch(const std::exception &error)
            {points.push_back(J::Object({{"x",J::Number(x)},{"variant",J::Number(variant)},{"valid","false"},
                {"reason",J::Quote(error.what())},{"witness",m_solveInvalidContext.empty()?"null":m_solveInvalidContext}}));}
        }
    }
    // Source-only counterfactuals on the SAME selected native volume rows.
    // Normalization uses each profile's required quadrature. No boundary,
    // matching coefficient, pressure, velocity or history is replaced.
    for(int variant:{-1,0,1})
    {
        m_solveInvalidContext.clear();
        try
        {
            Vector r(Size());
            for(const auto &ray:m_rays)
            {
                const auto profile=variant<0?ray.profile:Profile(SampleK(s,ray.coordinate),ray.coordinate,variant);
                if(variant==1)Require(profile->Reference().weight==1,"chi=1 frozen-source profile inactive");
                const auto normal=checkedNormalization(ray.coordinate,*profile);++normalizations;
                for(const auto &q:ray.points)
                {
                    RANSTensor g{};for(int a=0;a<m_dimension;++a)g[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);
                    const double k=q.v.Value(s,offset[K()]);const auto stress=InnerStress(k,profile->Original().EddyViscosity(Distance(q.v.coordinate)),g);
                    const auto value=WallSource(s,q.v.coordinate,stress.production,*profile,&normal);
                    for(int j=0;j<q.v.ids.size();++j)
                        r[offset[K()]+q.v.ids[j]]-=q.weight*q.v.phi[j]*value.kSource;
                }
            }
            auto only=resolution;for(auto &[name,value]:only.terms)value.assign(Size(),0.);
            only.terms.at("inner_k_production")=r;only.native["spatial_unmasked"]=r;
            const auto projected=RANSWallResolutionMassRate(*this,only,K(),{1.6,0,0});massApplications+=projected.componentRates.size()+1;
            std::vector<std::string> samples;
            for(double x:stations)for(double z:{0.,.0001,.0025,.01})
            {const auto row=Point(K(),{x,0,z});samples.push_back(J::Object({{"coordinate",J::Numbers(row.coordinate)},
                {"source_only_mass_rate",J::Number(row.Value(projected.rate))}}));}
            Vector load(size[K()]);for(int i=0;i<load.size();++i)load[i]=r[offset[K()]+i];
            weak.push_back(J::Object({{"variant",J::Number(variant)},{"valid","true"},{"weak_source_residual",J::Numbers(load)},
                {"mass_inverted_source_rate",J::Numbers(projected.rate)},{"mass_solve_defect",J::Number(Maximum(projected.freeMassResidual))},
                {"samples",J::List(samples)}}));
        }
        catch(const std::exception &error){weak.push_back(J::Object({{"variant",J::Number(variant)},{"valid","false"},{"reason",J::Quote(error.what())},
            {"witness",m_solveInvalidContext.empty()?"null":m_solveInvalidContext}}));}
    }
    const auto trace=DissipationTrace(s,time,data);std::vector<std::string> matching;
    for(const auto &[id,baseline]:trace[0])matching.push_back(J::Object({{"global_row",J::Number(offset[Epsilon()]+id)},
        {"retained",J::Number(s[offset[Epsilon()]+id])},{"baseline",J::Number(baseline)},
        {"policy_correction",J::Number(trace[1].at(id))},{"refreshed_target",J::Number(baseline+trace[1].at(id))},
        {"baseline_only_restriction",J::Boolean(m_baselineTraceIds.count(id))}}));
    Require(hash==RANSWallR3::Hash(s),"source readback modified selected coefficients");
    return J::Object({{"schema","1"},{"wall_policy",J::Quote(WallPolicyName())},{"time",J::Number(time)},
        {"state_sha256",J::Quote(hash)},{"matching_trace",J::List(matching)},
        {"fixed_external_coefficients",J::Constraints(data.fixedBoundary)},{"retained_state_unchanged","true"},{"production_governing_solves","0"},
        {"normalization_evaluations",J::Number(normalizations)},{"diagnostic_mass_inverse_applications",J::Number(massApplications)},
        {"native_mass_solve_defect",J::Number(Maximum(rate.freeMassResidual))},
        {"component_load_sum_defect",J::Number(componentLoadDefect)},{"component_mass_rate_sum_defect",J::Number(componentRateDefect)},{"points",J::List(points)},
        {"ray_budgets",J::List(budgets)},{"fixed_native_quadrature_source_comparisons",J::List(weak)},
        {"scope",J::Quote("same retained fields/time; chi counterfactual source only, not a lifted coupled state or production step")}});
}
}
