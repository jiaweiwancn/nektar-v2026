///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>

namespace Nektar
{
using namespace RANSWallAlgebra;
namespace J = RANSWallSolveJson;

// A temporal coupling operation, deliberately independent of every nonlinear
// convergence/search driver. All updates below have fraction exactly one.
RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::SolvePimpleSingleSweep(
    Vector s, const Vector &history, double alpha, double time,
    const RANSWallEvolutionData &data, int &iterations) const
{
    Require(s.size()==Size() && history.size()==s.size() &&
        std::isfinite(alpha) && alpha>0 && std::isfinite(time),"invalid PIMPLE stage input");
    const auto start=std::chrono::steady_clock::now();
    iterations=0;
    std::vector<std::pair<std::string,std::string>> record;
    int momentum=0,pressure=0,epsilon=0,k=0,massApplications=0;
    auto put=[&](const std::string &key,const std::string &value){record.emplace_back(key,value);};
    auto publish=[&](bool completed,const std::string &failure)
    {
        put("time",J::Number(time));put("alpha",J::Number(alpha));
        put("scheduled_complete",J::Boolean(completed));put("failure",J::Quote(failure));
        put("outer_sweeps",J::Number(iterations));
        put("momentum_predictor_solves",J::Number(momentum));
        put("pressure_system_solves",J::Number(pressure));
        put("epsilon_transport_solves",J::Number(epsilon));put("k_transport_solves",J::Number(k));
        put("mass_inverse_applications",J::Number(massApplications));
        put("additional_correction_or_restoration_solves","0");
        put("seconds",J::Number(std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()));
        m_pimpleRecords.push_back(J::Object(record));
        std::cout<<"EP-AWF PIMPLE "<<m_pimpleRecords.back()<<std::endl;
    };
    auto check=[&](const Vector &q,const char *phase)
    {
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
        const Factor momentumFactor(a),mass(m);
        ++momentum;const auto predictor=momentumFactor.Apply(rhs);
        auto predictorDefect=Multiply(a,predictor);
        for(int i=0;i<nv;++i)predictorDefect[i]-=rhs[i];
        Require(Maximum(predictorDefect)<=1e-10*std::max(1.,Maximum(rhs)),"PIMPLE predictor linear residual");
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
        ++pressure;const auto dz=Factor(schur).Apply(pressureRhs);
        auto pressureDefect=Multiply(schur,dz);
        for(int i=0;i<np;++i)pressureDefect[i]-=pressureRhs[i];
        Require(Maximum(pressureDefect)<=1e-10*std::max(1.,Maximum(pressureRhs)),"PIMPLE pressure linear residual");
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
        check(s,"flow");
        auto scalar=[&](int b)
        {
            const auto block=AssemblePicardBlock(b,s,history,alpha,time,data,external);
            // Same-state consistency verifies actual signed source/rate and
            // physical mass/history. No transport solve occurs in this check.
            const auto original=Residual(s,history,alpha,time,data);
            Vector agreement(block.rows.size());
            for(int i=0;i<block.rows.size();++i)agreement[i]=block.residual[i]-original[block.rows[i]];
            Require(Maximum(agreement)<1e-9*std::max(1.,Maximum(original)),"PIMPLE scalar source/storage mismatch");
            auto rhs=block.residual;for(auto &x:rhs)x=-x;
            if(b==Epsilon())++epsilon;else ++k;
            const auto delta=Factor(block.matrix).Apply(rhs);
            auto linear=Multiply(block.matrix,delta);
            for(int i=0;i<linear.size();++i)linear[i]-=rhs[i];
            Require(Maximum(linear)<1e-10*std::max(1.,Maximum(rhs)),"PIMPLE scalar linear residual");
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
            check(s,name.c_str());
        };
        scalar(Epsilon()); // corrected flow, old k and epsilon; target evaluated once in assembly
        scalar(K());       // corrected flow and freshly solved epsilon; inner rate, not rate bound
        Refresh(s,time,&data);
        const auto target=Constraints(s,time,data);double refreshed=0;
        for(int id:m_traceIds)refreshed=std::max(refreshed,std::abs(s[offset[Epsilon()]+id]-target.at(offset[Epsilon()]+id)));
        put("refreshed_trace_defect",J::Number(refreshed));
        const auto norms=Norms(Residual(s,history,alpha,time,data),target,m_robustUstar);
        put("refreshed_original_block_norms",J::Numbers(norms));
        m_lastFeasibleState=s;
        publish(true,"");return s;
    }
    catch(const std::exception &error){m_lastFeasibleState=s;publish(false,error.what());throw;}
}

RANSWallEvolutionPatch::Diagnostic RANSWallEvolutionPatch::PhysicalSamples(
    const Vector &s,double time,const RANSWallEvolutionData &data) const
{
    Diagnostic result;
    for(int b=0;b<Gauge();++b)
    {
        const auto name=std::to_string(b);
        for(const auto &q:Quadrature(b))
        {result["weights_"+name].push_back(q.weight);result["values_"+name].push_back(q.Value(s,offset[b]));}
    }
    for(const auto &q:Surface(0,2,0.,-1,8))
    {
        const auto p=q.test.coordinate;const auto profile=StateProfile(s,p,time,data);
        const auto normal=WallNormalization(s,p,*profile);
        SolverUtils::RANSTensor grad{};for(int a=0;a<m_dimension;++a)grad[PhysicalAxis(a)]=q.test.Gradient(s,offset[a]);
        const auto stress=SolverUtils::RANSWallInnerMeanStress(InnerStress(q.test.Value(s,offset[K()]),
            profile->Original().EddyViscosity(0.),grad),m_parameters.geometry.molecularNu,0.);
        result["wall_weights"].push_back(q.weight);
        result["wall_stress"].push_back(stress[0][2]);result["wall_B"].push_back(normal.b);
        result["wall_lambda"].push_back(normal.lambda);result["wall_K"].push_back(profile->K());
    }
    return result;
}
} // namespace Nektar
