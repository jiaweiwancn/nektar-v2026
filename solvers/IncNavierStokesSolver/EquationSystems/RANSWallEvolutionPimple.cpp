///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <filesystem>
#include <fstream>
#include <cstring>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>

namespace Nektar
{
using namespace RANSWallAlgebra;
namespace J = RANSWallSolveJson;

bool RANSRetainAttemptPacket(const std::string &policy,bool active,bool verification,long long step,int interval)
{return policy=="EveryAttempt" || active || verification || step<=1 || step%interval==0;}

// A temporal coupling operation, deliberately independent of every nonlinear
// convergence/search driver. All updates below have fraction exactly one.
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::SolvePimpleSingleSweep(
    Vector s, const Vector &history, double alpha, double time,
    const RANSWallEvolutionData &data, int &iterations) const
{
    RANSWallTiming::Scope timing("pimple.sweep");
    Require(s.size()==Size() && history.size()==s.size() &&
        std::isfinite(alpha) && alpha>0 && std::isfinite(time),"invalid PIMPLE stage input");
    const auto start=std::chrono::steady_clock::now();
    iterations=0;m_pendingPacket={};m_lastPacketPath.clear();m_lastPacketHash.clear();m_packetOutputFailure.clear();
    std::vector<std::pair<std::string,std::string>> record;
    int momentum=0,pressure=0,epsilon=0,k=0,massApplications=0;
    int momentumAttempts=0,pressureAttempts=0,epsilonAttempts=0,kAttempts=0;
    int boundCalls=0,boundActivations=0,boundIterations=0,boundMass=0,boundKKT=0;
    auto put=[&](const std::string &key,const std::string &value){record.emplace_back(key,value);};
    auto publish=[&](bool completed,const std::string &failure)
    {
        put("time",J::Number(time));put("alpha",J::Number(alpha));
        put("governing_solves_complete",J::Boolean(momentum==1 && pressure==1 && epsilon==1 && k==1));
        put("momentum_predictor_solve_attempts",J::Number(momentumAttempts));
        put("pressure_system_solve_attempts",J::Number(pressureAttempts));
        put("epsilon_transport_solve_attempts",J::Number(epsilonAttempts));
        put("k_transport_solve_attempts",J::Number(kAttempts));
        put("momentum_predictor_solve_completed",J::Boolean(momentum==1));
        put("pressure_system_solve_completed",J::Boolean(pressure==1));
        put("epsilon_transport_solve_completed",J::Boolean(epsilon==1));
        put("k_transport_solve_completed",J::Boolean(k==1));
        put("scheduled_complete",J::Boolean(completed));put("failure",J::Quote(failure));
        put("outer_sweeps",J::Number(iterations));
        put("momentum_predictor_solves",J::Number(momentum));
        put("pressure_system_solves",J::Number(pressure));
        put("epsilon_transport_solves",J::Number(epsilon));put("k_transport_solves",J::Number(k));
        put("mass_inverse_applications",J::Number(massApplications));
        put("additional_governing_equation_solves","0");
        put("additional_correction_or_restoration_solves",J::Number(boundKKT));
        put("bound_calls",J::Number(boundCalls));put("bound_activations",J::Number(boundActivations));
        put("bound_active_set_iterations",J::Number(boundIterations));put("bound_mass_factor_applications",J::Number(boundMass));
        put("bound_kkt_solves",J::Number(boundKKT));
        put("seconds",J::Number(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()));
        m_pimpleRecords.push_back(J::Object(record));
        std::cout<<"EP-AWF PIMPLE "<<m_pimpleRecords.back()<<std::endl;
    };
    auto check=[&](const Vector &q,const char *phase)
    {
        RANSWallTiming::Scope checkTiming(std::string("pimple.check.")+phase);
        Require(std::all_of(q.begin(),q.end(),[](double x){return std::isfinite(x);}),"PIMPLE nonfinite state");
        const auto domain=ScalarDomain(q,time,data);
        bool ground=false;const auto wall=WallTraceJSON(q,ground);
        put(std::string(phase)+"_domain",ScalarDomainJSON(domain));
        put(std::string(phase)+"_ground",wall);
        Require(domain.admissible && ground,std::string("PIMPLE ")+phase+" scalar/continuous-ground domain exit: "+domain.failure);
    };
    try
    {
        s=ExternalLift(s,time,data);m_lastFeasibleState=s;
        auto external=ExternalConstraints(time,data);
        for(int id:m_traceIds)external.erase(offset[Epsilon()]+id);
        put("history_sha256",J::Quote(RANSWallR3::Hash(history)));
        put("incoming_sha256",J::Quote(RANSWallR3::Hash(s)));
        check(s,"input");
        iterations=1;
        const auto base=s;
        RANSWallTiming::Scope flowTiming("flow.total");
        const auto flow=AssemblePicardBlock(0,s,history,alpha,time,data,external);
        const int nv=std::lower_bound(flow.rows.begin(),flow.rows.end(),DifferentialSize())-flow.rows.begin();
        const int np=flow.rows.size()-nv;
        WallMatrix a(nv,Vector(nv)),m(nv,Vector(nv)),g(nv,Vector(np)),d(np,Vector(nv));
        for(int i=0;i<nv;++i)
        {
            std::copy_n(flow.matrix[i].begin(),nv,a[i].begin());
            std::copy_n(flow.matrix[i].begin()+nv,np,g[i].begin());
            const int bi=flow.rows[i]/size[0]; // velocity components share the same space
            for(int j=0;j<nv;++j)
                if(flow.rows[j]>=offset[bi] && flow.rows[j]<offset[bi+1])
                    m[i][j]=m_mass[bi][flow.rows[i]-offset[bi]][flow.rows[j]-offset[bi]];
        }
        for(int i=0;i<np;++i)std::copy_n(flow.matrix[nv+i].begin(),nv,d[i].begin());
        Vector rhs(nv);
        for(int i=0;i<nv;++i)rhs[i]=-flow.residual[i];
        ++momentumAttempts;
        RANSWallTiming::Scope momentumFactorTiming("momentum.factorize");
        const Factor momentumFactor(a);momentumFactorTiming.Stop();
        RANSWallTiming::Scope massFactorTiming("pressure.mass_factorize");
        bool sameMass=m_velocityMassKey.size()==m.size();
        if(sameMass)for(int i=0;i<m.size();++i)if(m_velocityMassKey[i].size()!=m[i].size() || std::memcmp(m_velocityMassKey[i].data(),m[i].data(),m[i].size()*sizeof(double))!=0){sameMass=false;break;}
        const bool massHit=m_velocityMassReuse && m_velocityMassFactor && m_velocityMassRows==std::vector<int>(flow.rows.begin(),flow.rows.begin()+nv) && sameMass;
        std::shared_ptr<Factor> massOwner;
        if(massHit)massOwner=m_velocityMassFactor;
        else
        {
            massOwner=std::make_shared<Factor>(m);
            if(m_velocityMassReuse){m_velocityMassFactor=massOwner;m_velocityMassKey=m;m_velocityMassRows.assign(flow.rows.begin(),flow.rows.begin()+nv);}
        }
        const Factor &mass=*massOwner;massFactorTiming.Stop();
        put("velocity_mass_factor_reused",J::Boolean(massHit));
        put("velocity_mass_reuse_mode",J::Quote(m_velocityMassReuse?"ExactV1":"Off"));
        RANSWallTiming::Scope predictorTiming("momentum.solve_validate");
        const auto predictor=momentumFactor.Apply(rhs);
        auto predictorDefect=Multiply(a,predictor);
        for(int i=0;i<nv;++i)predictorDefect[i]-=rhs[i];
        Require(Maximum(predictorDefect)<=1e-10*std::max(1.,Maximum(rhs)),"PIMPLE predictor linear residual");
        ++momentum;predictorTiming.Stop();
        RANSWallTiming::Scope pressureBuildTiming("pressure.build");
        // B = (alpha M_free)^-1. This is never a full momentum inverse.
        WallMatrix bg(np,Vector(nv)),schur(np,Vector(np));
        for(int j=0;j<np;++j)
        {
            for(int i=0;i<nv;++i)rhs[i]=g[i][j];
            ++massApplications;bg[j]=mass.Apply(rhs);
            for(auto &x:bg[j])x/=alpha;
            const auto column=Multiply(d,bg[j]);
            for(int i=0;i<np;++i)schur[i][j]=column[i]-flow.matrix[nv+i][nv+j];
        }
        auto pressureRhs=Multiply(d,predictor);
        for(int i=0;i<np;++i)pressureRhs[i]+=flow.residual[nv+i];
        pressureBuildTiming.Stop();
        RANSWallTiming::Scope pressureSolveTiming("pressure.solve_validate");
        ++pressureAttempts;const auto dz=Factor(schur).Apply(pressureRhs);
        auto pressureDefect=Multiply(schur,dz);
        for(int i=0;i<np;++i)pressureDefect[i]-=pressureRhs[i];
        Require(Maximum(pressureDefect)<=1e-10*std::max(1.,Maximum(pressureRhs)),"PIMPLE pressure linear residual");
        ++pressure;pressureSolveTiming.Stop();
        RANSWallTiming::Scope correctionTiming("flow.correction_validate");
        Vector bGz(nv),Gz=Multiply(g,dz),du=predictor;
        for(int i=0;i<nv;++i)for(int j=0;j<np;++j)bGz[i]+=bg[j][i]*dz[j];
        for(int i=0;i<nv;++i){du[i]-=bGz[i];s[flow.rows[i]]+=du[i];}
        for(int j=0;j<np;++j)s[flow.rows[nv+j]]+=dz[j];
        auto step=du;step.insert(step.end(),dz.begin(),dz.end());
        auto defect=Multiply(flow.matrix,step);
        for(int i=0;i<defect.size();++i)defect[i]+=flow.residual[i];
        auto abGz=Multiply(a,bGz);Vector identity(nv),constraints(np);
        for(int i=0;i<nv;++i)identity[i]=defect[i]-(Gz[i]-abGz[i]);
        for(int i=0;i<np;++i)constraints[i]=defect[nv+i];
        const auto native=Residual(s,history,alpha,time,data);
        Vector regional(3);for(int b=InnerPressure();b<=Gauge();++b)
            for(int i=0;i<(b==Gauge()?1:size[b]);++i)
                regional[b-InnerPressure()]=std::max(regional[b-InnerPressure()],std::abs(native[offset[b]+i]));
        double boundary=0;for(auto [row,value]:external)boundary=std::max(boundary,std::abs(s[row]-value));
        Require(Maximum(identity)<=1e-9*std::max(1.,Maximum(Gz)) && Maximum(constraints)<1e-10 && Maximum(regional)<1e-10 && boundary<1e-13,
            "PIMPLE pressure constraint, boundary or momentum splitting identity failure");
        put("pressure_inverse",J::Quote("alpha^-1 consistent free-velocity mass inverse"));
        put("velocity_size",J::Number(nv));put("pressure_gauge_size",J::Number(np));
        put("predictor_linear_residual",J::Number(Maximum(predictorDefect)));
        put("pressure_linear_residual",J::Number(Maximum(pressureDefect)));
        put("corrected_continuity_gauge",J::Numbers(regional));put("boundary_defect",J::Number(boundary));
        put("momentum_splitting_identity_defect",J::Number(Maximum(identity)));
        put("frozen_momentum_defect",J::Number(Maximum(Vector(defect.begin(),defect.begin()+nv))));
        auto pressureDirection=Vector(Size());
        for(int j=0;j<np;++j)pressureDirection[flow.rows[nv+j]]=dz[j];
        const auto action=AffinePressureAction(pressureDirection,external);
        Vector actionError(nv);for(int i=0;i<nv;++i)actionError[i]=action[flow.rows[i]]-Gz[i];
        put("native_full_G_action_defect",J::Number(Maximum(actionError)));
        Require(Maximum(actionError)<1e-10,"PIMPLE omitted native pressure action");
        // Verification is opt-in and its extra solves are outside the trajectory counts.
        if(m_pimpleVerify)
        {
            const auto exact=PicardFlowCorrection(base,history,alpha,time,data,external);
            Vector difference=step;for(int i=0;i<step.size();++i)difference[i]-=exact.step[i];
            put("verification_exact_flow_coefficient_difference",J::Number(Maximum(difference)));
            put("verification_extra_momentum_inverse_applications",J::Number(np+1+4));
            put("verification_extra_pressure_solves","1");
            std::vector<std::string> actions;
            // Fixed native pressure action, frozen L, varying alpha tests the
            // small-dt approximation without changing physical trajectory.
            auto representative=Gz;
            if(Maximum(representative)==0.)for(int i=0;i<nv;++i)representative[i]=g[i][0];
            const auto inverseMass=mass.Apply(representative);
            for(double factor:{1.,4.,16.,64.})
            {
                auto shifted=a;for(int i=0;i<nv;++i)for(int j=0;j<nv;++j)shifted[i][j]+=(factor-1.)*alpha*m[i][j];
                const auto exactAction=Factor(shifted).Apply(representative);
                Vector error(nv);for(int i=0;i<nv;++i)error[i]=exactAction[i]-inverseMass[i]/(alpha*factor);
                auto norm=[&](const Vector &v){const auto mv=Multiply(m,v);return std::sqrt(std::max(0.,std::inner_product(v.begin(),v.end(),mv.begin(),0.)));};
                actions.push_back(J::Numbers(Vector{1/(alpha*factor),norm(error),norm(error)/std::max(1e-300,norm(exactAction))}));
            }
            put("verification_extra_mass_inverse_applications","1");
            put("verification_pressure_action_dt_absolute_relative",J::List(actions));
        }
        check(s,"flow");correctionTiming.Stop();flowTiming.Stop();
        const auto flowState=s;put("post_flow_sha256",J::Quote(RANSWallR3::Hash(s)));
        auto budget=[&](const char *phase,const Vector &q)
        {
            if(m_pimpleBudgetTime<0 || std::abs(time-m_pimpleBudgetTime)>.25/alpha)return;
            try
            {
                const double evaluationTime=std::string(phase)=="last_accepted"?time-1/alpha:time;
                put(std::string("k_budget_")+phase,PimpleKRateBudget(q,history,evaluationTime,data));
            }
            catch(const std::exception &error)
            {put(std::string("k_budget_")+phase,J::Object({{"success","false"},{"failure",J::Quote(error.what())}}));}
        };
        auto scalar=[&](int b)
        {
            RANSWallTiming::Scope scalarTiming(b==Epsilon()?"epsilon.total":"k.total");
            const auto block=AssemblePicardBlock(b,s,history,alpha,time,data,external);
            // Same-state consistency verifies actual signed source/rate and
            // physical mass/history. No transport solve occurs in this check.
            const auto original=Residual(s,history,alpha,time,data);
            Vector agreement(block.rows.size());
            for(int i=0;i<block.rows.size();++i)agreement[i]=block.residual[i]-original[block.rows[i]];
            Require(Maximum(agreement)<1e-9*std::max(1.,Maximum(original)),"PIMPLE scalar source/storage mismatch");
            auto rhs=block.residual;for(auto &x:rhs)x=-x;
            if(b==Epsilon())++epsilonAttempts;else ++kAttempts;
            RANSWallTiming::Scope scalarSolveTiming(b==Epsilon()?"epsilon.solve_validate":"k.solve_validate");
            const auto delta=Factor(block.matrix).Apply(rhs);
            auto linear=Multiply(block.matrix,delta);
            for(int i=0;i<linear.size();++i)linear[i]-=rhs[i];
            Require(Maximum(linear)<1e-10*std::max(1.,Maximum(rhs)),"PIMPLE scalar linear residual");
            if(b==Epsilon())++epsilon;else ++k;scalarSolveTiming.Stop();
            double traceColumns=0,traceMass=0,traceStorage=0,imposed=0;
            for(int i=0;i<block.rows.size();++i)
            {
                const int row=block.rows[i];
                const bool trace=b==Epsilon() && m_traceIds.count(row-offset[b]);
                if(trace)imposed=std::max(imposed,std::abs(block.residual[i]+delta[i]));
                else if(b==Epsilon())for(int j=0;j<block.rows.size();++j)if(m_traceIds.count(block.rows[j]-offset[b]))
                {
                    traceColumns=std::max(traceColumns,std::abs(block.matrix[i][j]));
                    const double massEntry=m_mass[b][row-offset[b]][block.rows[j]-offset[b]];
                    traceMass=std::max(traceMass,std::abs(alpha*massEntry));
                    traceStorage=std::max(traceStorage,std::abs(alpha*massEntry*(s[block.rows[j]]-history[block.rows[j]]+delta[j])));
                }
            }
            const auto preScalar=s;
            if(b==K())put("pre_k_sha256",J::Quote(RANSWallR3::Hash(s)));
            if(b==K()) {budget("last_accepted",base);budget("post_flow",flowState);budget("pre_k",preScalar);}
            for(int i=0;i<block.rows.size();++i)s[block.rows[i]]+=delta[i];
            const std::string name=b==Epsilon()?"epsilon":"k";
            put(name+"_frozen_linear_residual",J::Number(Maximum(linear)));
            put(name+"_same_state_source_storage_defect",J::Number(Maximum(agreement)));
            if(b==Epsilon())
            {
                put("imposed_trace_defect",J::Number(imposed));
                put("retained_trace_column_max",J::Number(traceColumns));put("trace_mass_column_max",J::Number(traceMass));
                put("trace_physical_storage_contribution_max",J::Number(traceStorage));
            }
            if(b==K())
            {
                const auto raw=s;put("raw_k_sha256",J::Quote(RANSWallR3::Hash(raw)));KBoundResult bounded;bounded.state=raw;bounded.success=true;
                if(PimpleKBoundEnabled())
                {
                    ++boundCalls;bounded=ApplyPimpleKBound(raw,preScalar,time,data);s=bounded.state;
                    boundActivations+=bounded.activated?1:0;boundIterations+=bounded.iterations;
                    boundMass+=bounded.massFactorApplications;boundKKT+=bounded.kktSolves;
                    put("k_bound",bounded.report);
                    RANSWallTiming::Scope boundAccountingTiming("k.bound_accounting");
                    Vector dk(block.rows.size());for(int i=0;i<dk.size();++i)dk[i]=s[block.rows[i]]-raw[block.rows[i]];
                    const auto added=Multiply(block.matrix,dk),storage=Multiply(PicardKMass(block.rows),dk);
                    auto corrected=linear;for(int i=0;i<corrected.size();++i)corrected[i]+=added[i];
                    Vector physicalStorage=storage;for(auto &x:physicalStorage)x*=alpha;
                    put("k_corrected_frozen_linear_residual",J::Number(Maximum(corrected)));
                    put("k_added_frozen_residual",J::Number(Maximum(added)));
                    put("k_added_storage_residual",J::Number(Maximum(physicalStorage)));
                    put("bound_residual_accounting_mass_applications","1");
                }
                put("corrected_k_sha256",J::Quote(RANSWallR3::Hash(s)));
                Require(RANSWallR3::Hash(Vector(preScalar.begin()+offset[Epsilon()],preScalar.begin()+offset[InnerPressure()]))==
                    RANSWallR3::Hash(Vector(s.begin()+offset[Epsilon()],s.begin()+offset[InnerPressure()])),"post-k operation changed solved epsilon");
                put("solved_epsilon_retained_bitwise","true");
                if(m_pimpleSnapshots)
                {
                    m_pendingPacket=[=](){
                    RANSWallTiming::Scope snapshotTiming("io.attempt_snapshot");
                    std::filesystem::create_directories(m_pimpleSnapshotDirectory);
                    const auto stem=(std::filesystem::path(m_pimpleSnapshotDirectory)/
                        ("k-t"+J::Number(time)+"-"+RANSWallR3::Hash(history).substr(0,12))).string();
                    const auto matrix=stem+".matrix.bin";
                    std::ofstream mat(matrix,std::ios::binary);
                    for(const auto &row:block.matrix)mat.write(reinterpret_cast<const char*>(row.data()),row.size()*sizeof(double));
                    mat.close();Require(mat.good(),"could not save frozen k matrix");
                    const auto packet=stem+".json";std::ofstream out(packet);
                    out<<J::Object({{"schema","1"},{"time",J::Number(time)},{"alpha",J::Number(alpha)},
                        {"committed",J::Numbers(base)},{"post_flow",J::Numbers(flowState)},{"pre_k",J::Numbers(preScalar)},
                        {"raw",J::Numbers(raw)},{"corrected",J::Numbers(s)},{"history",J::Numbers(history)},
                        {"fixed_reference",J::Numbers(m_reference)},{"fixed_external_coefficients",J::Constraints(data.fixedBoundary)},
                        {"reference_constraints",J::Constraints(m_referenceConstraints)},{"projection_work_controls",PimpleKBoundControls()},
                        {"transport_rows",J::Numbers(block.rows)},{"transport_residual",J::Numbers(block.residual)},
                        {"raw_increment",J::Numbers(delta)},{"raw_linear_residual",J::Numbers(linear)},
                        {"matrix_file",J::Quote(std::filesystem::absolute(matrix).string())},
                        {"matrix_sha256",J::Quote(RANSWallR3::HashFile(matrix))},
                        {"matrix_encoding",J::Quote("row-major native IEEE754 binary64")},
                        {"k_bound",PimpleKBoundEnabled()?bounded.report:"null"}})<<'\n';
                    out.close();Require(out.good(),"could not save raw/corrected k snapshot");
                    m_lastPacketPath=std::filesystem::absolute(packet).string();m_lastPacketHash=RANSWallR3::HashFile(packet);

                    };
                    const long long global=m_nativeAttempt;
                    const bool selected=RANSRetainAttemptPacket(m_packetPolicy,bounded.activated,m_pimpleVerify,global,m_packetInterval);
                    put("attempt_packet_policy",J::Quote(m_packetPolicy));
                    put("attempt_packet_selected",J::Boolean(selected));
                    if(selected || !bounded.success)
                    {FinishAttemptPacket(true);put("k_snapshot",J::Quote(m_lastPacketPath));put("k_snapshot_sha256",J::Quote(m_lastPacketHash));}
                }
                Require(bounded.success,"PIMPLE k corrector failed: "+bounded.failure);
            }
            check(s,name.c_str());
        };
        scalar(Epsilon()); // corrected flow, old k and epsilon; target evaluated once in assembly
        scalar(K());       // corrected flow and freshly solved epsilon; inner rate, not rate bound
        RANSWallTiming::Scope finalTiming("pimple.final_refresh_diagnostics");
        Refresh(s,time,&data);
        const auto target=Constraints(s,time,data);double refreshed=0;
        for(int id:m_traceIds)refreshed=std::max(refreshed,std::abs(s[offset[Epsilon()]+id]-target.at(offset[Epsilon()]+id)));
        put("refreshed_trace_defect",J::Number(refreshed));
        const auto norms=Norms(Residual(s,history,alpha,time,data),target,m_robustUstar);
        put("refreshed_original_block_norms",J::Numbers(norms));
        m_lastFeasibleState=s;
        finalTiming.Stop();publish(true,"");return s;
    }
    catch(const std::exception &error)
    {
        FlushRejectedPacket();
        try { if(m_pimpleSnapshots && m_lastPacketPath.empty() && m_packetOutputFailure.empty())
        {
            std::filesystem::create_directories(m_pimpleSnapshotDirectory);
            m_lastPacketPath=(std::filesystem::absolute(m_pimpleSnapshotDirectory)/("partial-t"+J::Number(time)+"-"+RANSWallR3::Hash(history).substr(0,12)+".json")).string();
            std::ofstream out(m_lastPacketPath);out<<J::Object({{"schema","1"},{"partial_attempt","true"},
                {"candidate",J::Numbers(s)},{"history",J::Numbers(history)},{"fixed_reference",J::Numbers(m_reference)},
                {"time",J::Number(time)},{"alpha",J::Number(alpha)},{"failure",J::Quote(error.what())},
                {"completed_U_p_epsilon_k",J::Numbers(Vector{double(momentum),double(pressure),double(epsilon),double(k)})},
                {"matrix_recovery",J::Quote("no completed post-k packet; frozen configuration/checkpoint/native callback history permit deterministic isolated stage replay")}})<<'\n';
            out.close();Require(out.good(),"could not preserve partial failed attempt");m_lastPacketHash=RANSWallR3::HashFile(m_lastPacketPath);
        }
        } catch(const std::exception &output){m_packetOutputFailure=output.what();}
        m_lastFeasibleState=s;publish(false,error.what());throw;
    }
}

void RANSWallEvolutionPatch::FinishAttemptPacket(bool rejected) const
{
    auto writer=std::move(m_pendingPacket);m_pendingPacket={};
    if(rejected && writer)
    {
        if(m_packetFault=="OnFlush")throw std::runtime_error("injected attempt-packet output failure");
        writer();
    }
}
void RANSWallEvolutionPatch::FlushRejectedPacket() const noexcept
{
    try { FinishAttemptPacket(true); }
    catch(const std::exception &error){m_packetOutputFailure=error.what();}
    catch(...){m_packetOutputFailure="nonstandard packet output exception";}
}
std::string RANSWallEvolutionPatch::AttemptPacketStatus() const
{
    return J::Object({{"policy",J::Quote(m_packetPolicy)},{"interval",J::Number(m_packetInterval)},
        {"native_attempt",J::Number(m_nativeAttempt)},{"pending_writer",J::Boolean(bool(m_pendingPacket))},
        {"output_failure",J::Quote(m_packetOutputFailure)},
        {"last_written_packet",m_lastPacketPath.empty()?"null":J::Quote(m_lastPacketPath)},
        {"sha256",m_lastPacketHash.empty()?"null":J::Quote(m_lastPacketHash)},
        {"crash_recovery",J::Quote("accepted native checkpoint/history authoritative; unselected in-flight packet is not SIGKILL-durable; frozen runtime/config plus checkpoint replay deterministically reconstructs its frozen A_k")}});
}

RANSWallEvolutionPatch::Diagnostic RANSWallEvolutionPatch::PhysicalSamples(
    const Vector &s,double time,const RANSWallEvolutionData &data) const
{
    RANSWallTiming::Scope timing("native.physical_samples");
    Diagnostic result;
    for(int b=0;b<Gauge();++b)
    {
        const auto name=std::to_string(b);
        for(const auto &q:Quadrature(b))
        {result["weights_"+name].push_back(q.weight);result["values_"+name].push_back(q.Value(s,offset[b]));}
    }
    for(const auto &q:SurfaceView(0,2,0.,-1,8))
    {
        const auto p=q.test.coordinate;const auto profile=StateProfile(s,p,time,data);
        SolverUtils::RANSCappedDestructionValues normal;
        if(!AcceptedAWF())normal=WallNormalization(s,p,*profile);
        SolverUtils::RANSTensor grad{};for(int a=0;a<m_dimension;++a)grad[PhysicalAxis(a)]=q.test.Gradient(s,offset[a]);
        const auto stress=SolverUtils::RANSWallInnerMeanStress(InnerStress(q.test.Value(s,offset[K()]),
            profile->Original().EddyViscosity(0.),grad),m_parameters.geometry.molecularNu,0.);
        double shear=stress[0][2],measure=q.weight;
        if(Composite())
        {const auto g=GeometryAt(p);shear=0.;for(int a=0;a<3;++a)for(int b=0;b<3;++b)shear+=g.tangent1[a]*stress[a][b]*g.normal[b];measure*=std::hypot(q.area[0],q.area[2]);}
        result["wall_weights"].push_back(measure);
        result["wall_stress"].push_back(shear);
        if(!AcceptedAWF()){result["wall_B"].push_back(normal.b);result["wall_lambda"].push_back(normal.lambda);}
        result["wall_K"].push_back(profile->K());
    }
    return result;
}
} // namespace Nektar
