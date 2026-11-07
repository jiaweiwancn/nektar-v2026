// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallKBoundProjection.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>

namespace Nektar
{
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;
bool RANSWallEvolutionPatch::PimpleKBoundWorkParameter(const std::string &name)
{
    // Only resource/telemetry controls. Floor, metric, inventory, physical
    // parameters and actual algorithm choices remain in their native signature.
    return name=="RANSWALLPIMPLEKBOUNDMAXADD" || name=="RANSWALLPIMPLEKBOUNDMAXTOTAL" ||
        name=="RANSWALLPIMPLEKBOUNDMAXSEPARATIONS" || name=="RANSWALLPIMPLEKBOUNDMAXSECONDS";
}
std::string RANSWallEvolutionPatch::PimpleKBoundControls() const
{
    return J::Object({{"schema","1"},{"per_add_completed_work",J::Number(m_kBoundMaxAdd)},
        {"total_completed_work",J::Number(m_kBoundMaxTotal)},{"separation_passes",J::Number(m_kBoundMaxSeparations)},
        {"wall_seconds",J::Number(m_kBoundMaxSeconds)},{"trace",J::Quote(m_kBoundTrace)},
        {"trace_file",J::Quote(m_kBoundTraceFile)},{"dump_file",J::Quote(m_kBoundDumpFile)},
        {"mathematics",J::Quote("unchanged normalized mass-whitened two-pass-MGS working-set algorithm; resource limits are not physical parameters")}});
}

std::string RANSWallEvolutionPatch::PimpleKBoundSignature() const
{
    if(!PimpleKBoundEnabled())return ""; // Preserve Off checkpoint signatures.
    const double K0=m_parameters.model.Q()*m_robustUstar*m_robustUstar;
    std::ostringstream out;out<<std::setprecision(17)<<":k_bound="<<m_pimpleKBound
        <<(ConstitutiveKInventory()?":version=MassQR-ConstitutiveV2:metric=physical_consistent_Mff:constraints=shared_constitutive_all_profiles_contact_rays_retained_cuts_cubic_ground":
            ":version=MassQR-v1:metric=physical_consistent_Mff:constraints=static_preK_rays_moving_witnesses_cubic_ground")
        <<":floor_relative="<<m_pimpleKFloorRelative<<":K0="<<K0<<":floor="<<K0*m_pimpleKFloorRelative
        <<":arithmetic_reserve="<<4096*std::numeric_limits<double>::epsilon()*K0;
    return out.str();
}

RANSWallEvolutionPatch::KBoundResult RANSWallEvolutionPatch::ApplyPimpleKBound(
    const Vector &raw,const Vector &preK,double time,const RANSWallEvolutionData &data,bool projectionOnly) const
{
    Require(raw.size()==Size() && preK.size()==raw.size(),"k-bound state shape differs");
    KBoundResult result;result.state=raw;
    if(!PimpleKBoundEnabled())
    {result.success=true;result.report=J::Object({{"mode",J::Quote("Off")},{"calls","0"},{"activations","0"}});return result;}
    using Clock=std::chrono::steady_clock;
    const auto started=Clock::now();auto seconds=[&](){return std::chrono::duration<double>(Clock::now()-started).count();};
    double separationSeconds=0,inventorySeconds=0;
    const double K0=m_parameters.model.Q()*m_robustUstar*m_robustUstar;
    const double floor=K0*m_pimpleKFloorRelative,reserve=4096*std::numeric_limits<double>::epsilon()*K0;
    const double target=floor+reserve,threshold=floor+reserve/2;
    struct Location {RANSWallGlobalPoint row;std::string inventory;};
    std::vector<Location> points,witnesses;
    std::vector<RANSWallCoordinate> retainedAnchors;
    std::vector<int> rows,index(size[K()],-1);
    std::unique_ptr<RANSWallKBound::Projection> projection;WallMatrix projectionMass;
    std::ofstream traceStream;if(m_kBoundTraceFile!="None")
    {traceStream.open(m_kBoundTraceFile,std::ios::app);Require(traceStream.good(),"cannot open k-bound progress trace");}
    std::vector<std::string> traceRecords;int traceSequence=0;
    auto record=[&](const std::string &value)
    {
        ++traceSequence;
        if(m_kBoundTrace!="Off")
        {if(traceRecords.size()<64)traceRecords.push_back(value);else traceRecords[16+(traceSequence-65)%48]=value;}
        if(traceStream && m_kBoundTrace=="Full"){traceStream<<value<<'\n';if(traceSequence%64==0)traceStream.flush();}
    };
    record(J::Object({{"kind",J::Quote("correction_start")},{"time",J::Number(time)},
        {"raw_sha256",J::Quote(RANSWallR3::Hash(raw))},{"controls",PimpleKBoundControls()}}));
    double rawMinimum=1e300,correctedMinimum=1e300,primal=0;
    Location rawWitness,correctedWitness;
    bool rawGroundOK=false,correctedGroundOK=false;
    const auto rawGround=WallTraceJSON(raw,rawGroundOK);
    std::string correctedGround=rawGround;
    Vector change(size[K()]);int separations=0,inventoryRebuilds=0;
    std::map<std::string,int> familyCounts,discovered;

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
        Require(!projectionOnly || ConstitutiveKInventory(),"projection-only legacy replay is not a production-entry waiver");
        const auto beforeStart=Clock::now();const auto before=ScalarDomain(preK,time,data,true,ConstitutiveKInventory());
        inventorySeconds+=std::chrono::duration<double>(Clock::now()-beforeStart).count();
        if(!projectionOnly)Require(before.admissible,"k-bound pre-k state is inadmissible: "+before.failure);
        for(const auto &p:before.observations)if(p.block==K() && p.inventory.find("normalization_ray")!=std::string::npos)
            points.push_back({Point(K(),p.coordinate,p.side),p.inventory+":pre_k"});
        rawMinimum=evaluate(raw,rawWitness);
        for(int separation=0;separation<m_kBoundMaxSeparations;++separation)
        {
            if(seconds()>=m_kBoundMaxSeconds){result.failure="k_bound_correction_wall_time_limit";break;}
            ++separations;const auto separationStart=Clock::now();
            correctedMinimum=evaluate(result.state,correctedWitness);
            auto ground=WallTraceR4(result.state);
            const auto worst=std::min_element(ground.begin(),ground.end(),[](const auto &a,const auto &b){return a.value-a.error<b.value-b.error;});
            Location cut=correctedWitness;double violation=threshold-correctedMinimum,lower=target;
            if(threshold-double(worst->value-worst->error)>violation)
            {violation=threshold-double(worst->value-worst->error);cut={worst->row,"continuous_ground_extremum"};lower+=2*double(worst->error);}
            separationSeconds+=std::chrono::duration<double>(Clock::now()-separationStart).count();
            if(violation<=0)
            {
                // Safe order: static/samplers first, every newly moving k ray
                // before normalization. Never turn a non-k failure into a cut.
                const auto inventoryStart=Clock::now();
                const auto domain=ScalarDomain(result.state,time,data,true,ConstitutiveKInventory(),retainedAnchors);
                inventorySeconds+=std::chrono::duration<double>(Clock::now()-inventoryStart).count();
                ++inventoryRebuilds;familyCounts=domain.counts;
                if(domain.admissible)
                {
                    bool movingViolation=false;
                    for(const auto &p:domain.observations)if(p.block==K() && p.value<threshold)
                    {cut={Point(K(),p.coordinate,p.side),p.inventory+":moving"};movingViolation=true;break;}
                    if(!movingViolation)
                    {
                        for(const auto &p:domain.observations)if(p.block==K() && p.value<correctedMinimum)
                        {correctedMinimum=p.value;correctedWitness={Point(K(),p.coordinate,p.side),p.inventory};}
                        if(ConstitutiveKInventory())
                        {
                            const auto fullStart=Clock::now();const auto full=ScalarDomain(result.state,time,data,false,false,retainedAnchors);
                            inventorySeconds+=std::chrono::duration<double>(Clock::now()-fullStart).count();
                            if(!full.admissible){result.failure="k_bound_final_domain_exit: "+full.phase+": "+full.failure;break;}
                        }
                        result.success=true;break;
                    }
                }
                else if(domain.failure=="negative_inner_k" || domain.failure=="nonpositive_K_or_outer_k")
                    cut={Point(K(),domain.limitingPoint),domain.limitingInventory+":moving"};
                else {result.failure="k_bound_non_k_domain_exit: "+domain.phase+": "+domain.failure;break;}
                if(!ConstitutiveKInventory())points.push_back(cut);
            }
            if(ConstitutiveKInventory())points.push_back(cut); // retain inactive physical cuts as well
            ++discovered[cut.inventory];
            result.activated=true;
            Vector row(rows.size());
            for(int j=0;j<cut.row.ids.size();++j)
            {const int col=index[cut.row.ids[j]];if(col>=0)row[col]+=cut.row.phi[j];}
            if(Maximum(row)==0.){result.failure="incompatible_prescribed_k_floor";break;}
            if(!projection)
            {
                projectionMass=PicardKMass(rows);
                RANSWallKBound::Controls controls{m_kBoundMaxAdd,m_kBoundMaxTotal,std::max(.000001,m_kBoundMaxSeconds-seconds())};
                auto event=[&](const RANSWallKBound::Event &e)
                {
                    auto identity=[&](int id)
                    {
                        if(id<0 || id>=projection->cuts.size())return std::string("null");
                        const auto &c=projection->cuts[id];const auto &w=witnesses[c.witness];
                        return J::Object({{"cut",J::Number(id)},{"witness",J::Number(c.witness)},
                            {"row_identity",J::Quote(c.rowIdentity)},{"constraint_identity",J::Quote(c.constraintIdentity)},
                            {"family",J::Quote(w.inventory)},{"coordinate",J::Numbers(w.row.coordinate)}});
                    };
                    record(J::Object({{"kind",J::Quote(e.kind)},{"sequence",J::Number(traceSequence+1)},
                        {"add",J::Number(e.add)},{"per_add_work",J::Number(e.localWork)},{"total_work",J::Number(e.totalWork)},
                        {"entering",identity(e.entering)},{"leaving",identity(e.leaving)},{"multiplier_step",J::Number(e.step)},
                        {"active_size",J::Number(e.activeSize)},{"active_fingerprint",J::Quote(e.activeFingerprint)},
                        {"qr_min_diagonal",J::Number(e.minimumQR)},{"qr_max_diagonal",J::Number(e.maximumQR)},
                        {"last_completed_primal_objective",J::Number(e.objective)},
                        {"working_dual_violation",J::Number(e.dualViolation)},{"seconds",J::Number(seconds())}}));
                };
                projection=std::make_unique<RANSWallKBound::Projection>(projectionMass,controls,event);
            }
            witnesses.push_back(cut);retainedAnchors.push_back(cut.row.coordinate);
            record(J::Object({{"kind",J::Quote("separation")},{"sequence",J::Number(traceSequence+1)},
                {"pass",J::Number(separations)},{"inventory_rebuilds",J::Number(inventoryRebuilds)},
                {"floor_violation_known_inventory",J::Number(std::max(0.,floor-correctedMinimum))},
                {"ground_minimum",J::Number(double(worst->value))},{"entering_family",J::Quote(cut.inventory)},
                {"entering_coordinate",J::Numbers(cut.row.coordinate)},
                {"last_completed_primal_objective",J::Number(projection->Objective())},{"seconds",J::Number(seconds())}}));
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
    if(!points.empty())
    {Location fixed;const double minimum=evaluate(result.state,fixed);if(!result.success || minimum<correctedMinimum){correctedMinimum=minimum;correctedWitness=fixed;}}
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
        for(int i=0;i<projection->solvedActive.size();++i)
        {
            const auto &c=projection->cuts[projection->solvedActive[i]];const auto &w=witnesses[c.witness];
            const double value=w.row.Value(result.state,offset[K()]);primal=std::max(primal,floor-value);
            active.push_back(J::Object({{"coordinate",J::Numbers(w.row.coordinate)},{"geometry",J::Number(w.row.geometry)},
                {"inventory",J::Quote(w.inventory)},{"raw",J::Number(w.row.Value(raw,offset[K()]))},{"corrected",J::Number(value)},
                {"floor",J::Number(floor)},{"normalized_dual_multiplier",J::Number(projection->solvedLambda[i])}}));
        }
        if(result.success && projection->Stationarity()>1e-10*std::max(K0,RANSWallKBound::Norm(projection->y)))
        {result.success=false;result.failure="k_bound_stationarity_failure";}
    }
    // Incomplete working-set data is never paired with the previous primal.
    // Assess the actual returned state independently; an iteration limit is
    // not an infeasibility certificate and active-row feasibility is insufficient.
    const auto finalStart=Clock::now();ScalarDomainResult finalDomain;std::string finalDomainError;
    try{finalDomain=ScalarDomain(result.state,time,data,true,ConstitutiveKInventory(),retainedAnchors);}
    catch(const std::exception &error){finalDomainError=error.what();}
    inventorySeconds+=std::chrono::duration<double>(Clock::now()-finalStart).count();
    double groundFloorDefect=0;for(const auto &g:WallTraceR4(result.state))groundFloorDefect=std::max(groundFloorDefect,floor-double(g.value-g.error));
    const double fullPrimal=std::max({0.,floor-correctedMinimum,primal,groundFloorDefect});
    if(result.success && projection)
    {
        const double scale=std::max(K0*K0,projection->Objective());
        if(projection->DualViolation()>1e-12 || projection->Complementarity()>1e-10*scale ||
            std::abs(projection->Objective()-projection->DualObjective())>1e-9*scale)
        {result.success=false;result.failure="k_bound_completed_qp_optimality_check";}
    }
    std::string dump="null";
    if(projection && m_kBoundDumpFile!="None")
    {
        const auto massFile=m_kBoundDumpFile+".Mff.bin",cutsFile=m_kBoundDumpFile+".Cwhite.bin";
        std::ofstream mf(massFile,std::ios::binary),cf(cutsFile,std::ios::binary);
        for(const auto &r:projectionMass)mf.write(reinterpret_cast<const char*>(r.data()),r.size()*sizeof(double));
        std::vector<std::string> pool;
        for(const auto &c:projection->cuts)
        {
            cf.write(reinterpret_cast<const char*>(c.row.data()),c.row.size()*sizeof(double));const auto &w=witnesses[c.witness];
            pool.push_back(J::Object({{"rhs",J::Number(c.rhs)},{"witness",J::Number(c.witness)},
                {"row_identity",J::Quote(c.rowIdentity)},{"family",J::Quote(w.inventory)},
                {"coordinate",J::Numbers(w.row.coordinate)},{"physical_ids",J::Numbers(w.row.ids)},{"physical_phi",J::Numbers(w.row.phi)}}));
        }
        mf.close();cf.close();Require(mf.good() && cf.good(),"cannot save projection mass/cut pool");
        dump=J::Object({{"mass_file",J::Quote(massFile)},{"mass_sha256",J::Quote(RANSWallR3::HashFile(massFile))},
            {"mass_definition",J::Quote("PicardKMass on exact free k rows; physical consistent mass, NOT A_k")},
            {"free_rows",J::Numbers(rows)},{"cuts_file",J::Quote(cutsFile)},{"cuts_sha256",J::Quote(RANSWallR3::HashFile(cutsFile))},
            {"pool",J::List(pool)},{"solved_active",J::Numbers(projection->solvedActive)},
            {"solved_lambda",J::Numbers(projection->solvedLambda)},{"whitened_primal",J::Numbers(projection->y)},
            {"working_active",J::Numbers(projection->active)},{"working_lambda",J::Numbers(projection->lambda)},
            {"working_complete",J::Boolean(projection->workComplete)}});
        std::ofstream output(m_kBoundDumpFile);output<<dump<<'\n';Require(output.good(),"cannot write projection replay dump");
    }
    record(J::Object({{"kind",J::Quote("correction_end")},{"success",J::Boolean(result.success)},
        {"failure",J::Quote(result.failure)},{"seconds",J::Number(seconds())}}));if(traceStream)traceStream.flush();
    auto location=[&](double value,const Location &w){return J::Object({{"value",J::Number(value)},
        {"coordinate",J::Numbers(w.row.coordinate)},{"geometry",J::Number(w.row.geometry)},{"inventory",J::Quote(w.inventory)}});};
    std::vector<std::pair<std::string,std::string>> counts,cuts;
    for(auto [name,count]:familyCounts)counts.push_back({name,J::Number(count)});
    for(auto [name,count]:discovered)cuts.push_back({name,J::Number(count)});
    result.report=J::Object({{"mode",J::Quote(m_pimpleKBound)},{"algorithm",J::Quote(ConstitutiveKInventory()?"MassQR-ConstitutiveV2":"MassQR-v1")},
        {"inventory_policy",J::Quote(m_pimpleKInventory)},{"projection_only",J::Boolean(projectionOnly)},
        {"work_controls",PimpleKBoundControls()},{"telemetry_schema","1"},
        {"separation_passes",J::Number(separations)},{"inventory_rebuilds",J::Number(inventoryRebuilds)},
        {"inventory_counts",J::Object(counts)},{"discovered_constraint_families",J::Object(cuts)},
        {"calls","1"},{"activations",J::Number(result.activated?1:0)},{"success",J::Boolean(result.success)},{"failure",J::Quote(result.failure)},
        {"floor_relative",J::Number(m_pimpleKFloorRelative)},{"K0",J::Number(K0)},{"floor",J::Number(floor)},{"arithmetic_reserve",J::Number(reserve)},
        {"raw_minimum_scope",J::Quote("static and available pre-k ray locations; no inadmissible raw normalization evaluation")},
        {"raw_minimum",location(rawMinimum,rawWitness)},{"corrected_minimum",location(correctedMinimum,correctedWitness)},
        {"raw_ground",rawGround},{"corrected_ground",correctedGround},
        {"raw_sha256",J::Quote(RANSWallR3::Hash(raw))},{"corrected_sha256",J::Quote(RANSWallR3::Hash(result.state))},
        {"mass_norm_correction",J::Number(massNorm)},{"mass_norm_raw_step",J::Number(stepNorm)},
        {"correction_to_step_ratio_defined",J::Boolean(stepNorm>0.)},{"correction_to_step_ratio",J::Number(stepNorm>0?massNorm/stepNorm:0.)},{"integral_total",J::Number(result.integral)},
        {"integrals_inner_outer_near_positive_negative",J::Numbers(integrals)},{"near_ground_height",J::Number(nearHeight)},
        {"active_constraints",J::List(active)},{"active_set_iterations",J::Number(result.iterations)},
        {"mass_factorizations",J::Number(projection?1:0)},{"mass_factor_applications",J::Number(result.massFactorApplications)},
        {"mass_matrix_applications","2"},{"kkt_solves",J::Number(result.kktSolves)},{"governing_equation_solves","0"},
        {"primal_floor_defect",J::Number(fullPrimal)},{"active_row_floor_defect",J::Number(std::max(0.,primal))},
        {"current_inventory_passed",J::Boolean(finalDomain.admissible && finalDomainError.empty())},
        {"current_inventory",ScalarDomainJSON(finalDomain)},{"current_inventory_exception",J::Quote(finalDomainError)},
        {"ground_passed",J::Boolean(correctedGroundOK)},{"ground_floor_defect",J::Number(groundFloorDefect)},{"finite_pool_normalized_violation",J::Number(projection?projection->PoolViolation():0)},
        {"completed_primal_dual_coherent","true"},{"unfinished_working_set",J::Boolean(projection && !projection->workComplete)},
        {"last_completed_add",J::Number(projection?projection->solvedAdd:0)},{"adds",J::Number(projection?projection->adds:0)},
        {"maximum_work_one_add",J::Number(projection?projection->maximumAddIterations:0)},
        {"work_guard_checks",J::Number(projection?projection->workChecks:0)},
        {"rejected_total_limit_checks",J::Number(projection?projection->totalLimitChecks:0)},
        {"rejected_add_limit_checks",J::Number(projection?projection->addLimitChecks:0)},
        {"rejected_time_limit_checks",J::Number(projection?projection->timeLimitChecks:0)},
        {"legacy_guard_counter",J::Number(projection?projection->iterations+projection->totalLimitChecks:0)},
        {"unique_evaluation_rows",J::Number(projection?projection->rowVisits.size():0)},
        {"repeated_rows",J::Number(projection?projection->repeatedRows:0)},
        {"repeated_identical_constraints",J::Number(projection?projection->repeatedConstraints:0)},
        {"rediscovered_identical_constraints",J::Number(projection?projection->rediscoveredConstraints:0)},
        {"repeated_working_sets",J::Number(projection?projection->repeatedWorkingSets:0)},
        {"dependent_entries",J::Number(projection?projection->dependentEntries:0)},
        {"maximum_active_size",J::Number(projection?projection->maximumActive:0)},
        {"minimum_observed_qr_residual",J::Number(projection?projection->minimumObservedQR:0)},
        {"last_completed_qr_diagonal_range",J::Numbers(Vector{projection?projection->lastQRMin:0,projection?projection->lastQRMax:0})},
        {"primal_objective",J::Number(projection?projection->Objective():0)},
        {"dual_objective",J::Number(projection?projection->DualObjective():0)},
        {"primal_dual_gap_working_set",J::Number(projection?projection->Objective()-projection->DualObjective():0)},
        {"dual_feasibility_defect",J::Number(projection?projection->DualViolation():0)},
        {"seconds_total",J::Number(seconds())},{"seconds_algebra",J::Number(projection?projection->algebraSeconds:0)},
        {"seconds_factorization",J::Number(projection?projection->factorizationSeconds:0)},
        {"seconds_separation",J::Number(separationSeconds)},{"seconds_inventory",J::Number(inventorySeconds)},
        {"trace_events",J::Number(traceSequence)},{"bounded_trace",J::List(traceRecords)},
        {"full_diagnostic_dump",J::Quote(m_kBoundDumpFile)},{"stationarity_mass_dual",J::Number(projection?projection->Stationarity():0.)},
        {"complementarity",J::Number(projection?projection->Complementarity():0.)},
        {"scope",J::Quote("declared k inventories and continuous cubic ground; no volume-wide polynomial certificate")}});
    return result;
}
std::string RANSWallEvolutionPatch::PimpleKRateBudget(const Vector &q,const Vector &history,
    double time,const RANSWallEvolutionData &data) const
{
    Require(ScalarDomain(q,time,data).admissible,"k-rate budget requires an admissible selected state");
    auto resolution=ResidualResolutionComponents(q,history,0.,time,data,false);
    // Components deliberately do not populate the independent native load.
    // Supply it explicitly before the consistent-mass rate reconstruction.
    resolution.native["spatial_unmasked"]=Residual(q,history,0.,time,data,false,false);
    const auto ground=WallTraceR4(q);
    const auto worst=std::min_element(ground.begin(),ground.end(),[](const auto &a,const auto &b){return a.value<b.value;});
    const auto rate=RANSWallResolutionMassRate(*this,resolution,K(),worst->row.coordinate);
    const auto fixed=Point(K(),{1.60,0.,0.});
    std::vector<std::pair<std::string,std::string>> terms;
    Vector sum(2);
    for(const auto &[name,value]:rate.componentRates)
    {
        const Vector values{worst->row.Value(value),fixed.Value(value)};
        sum[0]+=values[0];sum[1]+=values[1];terms.push_back({name,J::Numbers(values)});
    }
    const Vector total{rate.contactRate,fixed.Value(rate.rate)};
    for(int i=0;i<2;++i)sum[i]-=total[i];
    const auto constraints=Constraints(q,time,data);double trace=0.;
    for(int id:m_traceIds)trace=std::max(trace,std::abs(q[offset[Epsilon()]+id]-constraints.at(offset[Epsilon()]+id)));
    return J::Object({{"success","true"},{"state_time",J::Number(time)},{"state_sha256",J::Quote(RANSWallR3::Hash(q))},
        {"ground_point",J::Numbers(worst->row.coordinate)},{"ground_geometry",J::Number(worst->row.geometry)},
        {"ground_value",J::Number(worst->value)},{"fixed_witness",J::Numbers(Vector{1.6,0,0})},{"fixed_geometry",J::Number(fixed.geometry)},
        {"rates_ground_and_fixed",J::Numbers(total)},{"component_rates",J::Object(terms)},
        {"component_sum_minus_native_rate",J::Numbers(sum)},{"mass_solve_residual",J::Number(Maximum(rate.freeMassResidual))},
        {"diagnostic_mass_rate_inverse_applications",J::Number(rate.componentRates.size()+1)},
        {"epsilon_trace_defect",J::Number(trace)},
        {"scope",J::Quote("instantaneous spatial k rate with zero essential rates; net destruction once; not a full-sweep h-to-zero consistency claim")}});
}

}
