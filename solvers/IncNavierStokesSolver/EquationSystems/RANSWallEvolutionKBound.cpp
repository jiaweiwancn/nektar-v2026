// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallKBoundProjection.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace Nektar
{
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;
std::string RANSWallEvolutionPatch::PimpleKBoundSignature() const
{
    if(!PimpleKBoundEnabled())return ""; // Preserve Off checkpoint signatures.
    const double K0=m_parameters.model.Q()*m_robustUstar*m_robustUstar;
    std::ostringstream out;out<<std::setprecision(17)<<":k_bound="<<m_pimpleKBound
        <<":version=MassQR-v1:metric=physical_consistent_Mff:constraints=static_preK_rays_moving_witnesses_cubic_ground"
        <<":floor_relative="<<m_pimpleKFloorRelative<<":K0="<<K0<<":floor="<<K0*m_pimpleKFloorRelative
        <<":arithmetic_reserve="<<4096*std::numeric_limits<double>::epsilon()*K0;
    return out.str();
}

RANSWallEvolutionPatch::KBoundResult RANSWallEvolutionPatch::ApplyPimpleKBound(
    const Vector &raw,const Vector &preK,double time,const RANSWallEvolutionData &data) const
{
    Require(raw.size()==Size() && preK.size()==raw.size(),"k-bound state shape differs");
    KBoundResult result;result.state=raw;
    if(!PimpleKBoundEnabled())
    {result.success=true;result.report=J::Object({{"mode",J::Quote("Off")},{"calls","0"},{"activations","0"}});return result;}
    const double K0=m_parameters.model.Q()*m_robustUstar*m_robustUstar;
    const double floor=K0*m_pimpleKFloorRelative,reserve=4096*std::numeric_limits<double>::epsilon()*K0;
    const double target=floor+reserve,threshold=floor+reserve/2;
    struct Location {RANSWallGlobalPoint row;std::string inventory;};
    std::vector<Location> points,witnesses;
    std::vector<int> rows,index(size[K()],-1);
    std::unique_ptr<RANSWallKBound::Projection> projection;
    double rawMinimum=1e300,correctedMinimum=1e300,primal=0;
    Location rawWitness,correctedWitness;
    bool rawGroundOK=false,correctedGroundOK=false;
    const auto rawGround=WallTraceJSON(raw,rawGroundOK);
    std::string correctedGround=rawGround;
    Vector change(size[K()]);
    auto evaluate=[&](const Vector &s,Location &where)
    {
        double minimum=1e300;
        for(const auto &p:points)
        {const double value=p.row.Value(s,offset[K()]);if(value<minimum){minimum=value;where=p;}}
        return minimum;
    };
    try
    {
        Require(raw.size()==Size() && preK.size()==raw.size(),"k-bound state shape differs");
        InitialiseScalarLocations();
        const auto external=ExternalConstraints(time,data);
        for(int i=0;i<size[K()];++i)if(!external.count(offset[K()]+i))
        {index[i]=rows.size();rows.push_back(offset[K()]+i);}
        for(const auto &p:m_scalarLocations)if(p.block==K())points.push_back({p.row,p.inventory});
        // preK is admissible. Preserve its moving normalization-ray inventory
        // before touching the raw candidate or evaluating its normalization.
        const auto before=ScalarDomain(preK,time,data,true);
        Require(before.admissible,"k-bound pre-k state is inadmissible: "+before.failure);
        for(const auto &p:before.observations)if(p.block==K() && p.inventory.find("normalization_ray")!=std::string::npos)
            points.push_back({Point(K(),p.coordinate,p.side),p.inventory+":pre_k"});
        rawMinimum=evaluate(raw,rawWitness);
        for(int separation=0;separation<500;++separation)
        {
            correctedMinimum=evaluate(result.state,correctedWitness);
            auto ground=WallTraceR4(result.state);
            const auto worst=std::min_element(ground.begin(),ground.end(),[](const auto &a,const auto &b){return a.value-a.error<b.value-b.error;});
            Location cut=correctedWitness;double violation=threshold-correctedMinimum,lower=target;
            if(threshold-double(worst->value-worst->error)>violation)
            {violation=threshold-double(worst->value-worst->error);cut={worst->row,"continuous_ground_extremum"};lower+=2*double(worst->error);}
            if(violation<=0)
            {
                // Safe order: static/samplers first, every newly moving k ray
                // before normalization. Never turn a non-k failure into a cut.
                const auto domain=ScalarDomain(result.state,time,data,true);
                if(domain.admissible)
                {
                    bool movingViolation=false;
                    for(const auto &p:domain.observations)if(p.block==K() && p.value<threshold)
                    {cut={Point(K(),p.coordinate,p.side),p.inventory+":moving"};movingViolation=true;break;}
                    if(!movingViolation){result.success=true;break;}
                }
                else if(domain.failure=="negative_inner_k" || domain.failure=="nonpositive_K_or_outer_k")
                    cut={Point(K(),domain.limitingPoint),domain.limitingInventory+":moving"};
                else {result.failure="k_bound_non_k_domain_exit: "+domain.phase+": "+domain.failure;break;}
                points.push_back(cut);
            }
            result.activated=true;
            Vector row(rows.size());
            for(int j=0;j<cut.row.ids.size();++j)
            {const int col=index[cut.row.ids[j]];if(col>=0)row[col]+=cut.row.phi[j];}
            if(Maximum(row)==0.){result.failure="incompatible_prescribed_k_floor";break;}
            if(!projection)projection=std::make_unique<RANSWallKBound::Projection>(PicardKMass(rows));
            witnesses.push_back(cut);
            if(!projection->Add(row,lower-cut.row.Value(raw,offset[K()]),witnesses.size()-1))
            {result.failure=projection->failure;break;}
            // Every subproblem is relative to the original k*, never to the
            // preceding attempted repair. Non-k coefficients are never written.
            result.state=raw;
            for(int i=0;i<rows.size();++i)result.state[rows[i]]+=projection->delta[i];
        }
        if(!result.success && result.failure.empty())result.failure="k_bound_separation_limit";
    }
    catch(const std::exception &error){result.failure=error.what();}
    if(!points.empty())correctedMinimum=evaluate(result.state,correctedWitness);
    correctedGround=WallTraceJSON(result.state,correctedGroundOK);
    for(int i=0;i<size[K()];++i)change[i]=result.state[offset[K()]+i]-raw[offset[K()]+i];
    const auto load=Multiply(m_mass[K()],change);
    const double massNorm=std::sqrt(std::max(0.,RANSWallKBound::Dot(change,load)));
    Vector step(size[K()]);for(int i=0;i<size[K()];++i)step[i]=raw[offset[K()]+i]-preK[offset[K()]+i];
    const auto stepLoad=Multiply(m_mass[K()],step);
    const double stepNorm=std::sqrt(std::max(0.,RANSWallKBound::Dot(step,stepLoad)));
    result.integral=RANSWallKBound::Dot(m_constants[K()],load);
    Vector integrals(5);double nearHeight=1e300;
    for(const auto &el:fields[K()]->Elements())if(el.lower[2]==0)nearHeight=std::min(nearHeight,el.upper[2]);
    for(const auto &q:Quadrature(K()))
    {
        const double x=q.Value(change),v=q.weight*x;
        integrals[q.coordinate[2]<m_parameters.geometry.height?0:1]+=v;
        if(q.coordinate[2]<=nearHeight)integrals[2]+=v;
        integrals[x>=0?3:4]+=v;
    }
    std::vector<std::string> active;
    if(projection)
    {
        result.iterations=projection->iterations;result.massFactorApplications=projection->massFactorApplications;result.kktSolves=projection->kktSolves;
        for(int i=0;i<projection->active.size();++i)
        {
            const auto &c=projection->cuts[projection->active[i]];const auto &w=witnesses[c.witness];
            const double value=w.row.Value(result.state,offset[K()]);primal=std::max(primal,floor-value);
            active.push_back(J::Object({{"coordinate",J::Numbers(w.row.coordinate)},{"geometry",J::Number(w.row.geometry)},
                {"inventory",J::Quote(w.inventory)},{"raw",J::Number(w.row.Value(raw,offset[K()]))},{"corrected",J::Number(value)},
                {"floor",J::Number(floor)},{"normalized_dual_multiplier",J::Number(projection->lambda[i])}}));
        }
        if(result.success && projection->Stationarity()>1e-10*std::max(K0,RANSWallKBound::Norm(projection->y)))
        {result.success=false;result.failure="k_bound_stationarity_failure";}
    }
    auto location=[&](double value,const Location &w){return J::Object({{"value",J::Number(value)},
        {"coordinate",J::Numbers(w.row.coordinate)},{"geometry",J::Number(w.row.geometry)},{"inventory",J::Quote(w.inventory)}});};
    result.report=J::Object({{"mode",J::Quote(m_pimpleKBound)},{"algorithm",J::Quote("MassQR-v1")},
        {"calls","1"},{"activations",J::Number(result.activated?1:0)},{"success",J::Boolean(result.success)},{"failure",J::Quote(result.failure)},
        {"floor_relative",J::Number(m_pimpleKFloorRelative)},{"K0",J::Number(K0)},{"floor",J::Number(floor)},{"arithmetic_reserve",J::Number(reserve)},
        {"raw_minimum",location(rawMinimum,rawWitness)},{"corrected_minimum",location(correctedMinimum,correctedWitness)},
        {"raw_ground",rawGround},{"corrected_ground",correctedGround},
        {"raw_sha256",J::Quote(RANSWallR3::Hash(raw))},{"corrected_sha256",J::Quote(RANSWallR3::Hash(result.state))},
        {"mass_norm_correction",J::Number(massNorm)},{"mass_norm_raw_step",J::Number(stepNorm)},
        {"correction_to_step_ratio",J::Number(stepNorm>0?massNorm/stepNorm:0.)},{"integral_total",J::Number(result.integral)},
        {"integrals_inner_outer_near_positive_negative",J::Numbers(integrals)},{"near_ground_height",J::Number(nearHeight)},
        {"active_constraints",J::List(active)},{"active_set_iterations",J::Number(result.iterations)},
        {"mass_factorizations",J::Number(projection?1:0)},{"mass_factor_applications",J::Number(result.massFactorApplications)},
        {"mass_matrix_applications","2"},{"kkt_solves",J::Number(result.kktSolves)},{"governing_equation_solves","0"},
        {"primal_floor_defect",J::Number(std::max(0.,primal))},{"stationarity_mass_dual",J::Number(projection?projection->Stationarity():0.)},
        {"complementarity",J::Number(projection?projection->Complementarity():0.)},
        {"scope",J::Quote("declared k inventories and continuous cubic ground; no volume-wide polynomial certificate")}});
    return result;
}
}
