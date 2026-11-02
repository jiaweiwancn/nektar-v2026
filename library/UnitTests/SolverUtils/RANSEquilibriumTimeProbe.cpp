///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
using Engine = Nektar::RANSWallTimeIntegration;
using Vector = Engine::Vector;
using Snapshot = Engine::Snapshot;
const Vector initial{1.,1.};
const std::array<double,2> decay{2.,.7};
int checks=0;
std::map<std::string,double> errors;
struct Rejection { std::string name,message; };
std::vector<Rejection> rejections;

void Check(bool valid,const std::string &message)
{
    ++checks;
    if(!valid)throw std::runtime_error("M4 native time probe: "+message);
}
double Near(double actual,double expected,const std::string &name,double tolerance=1e-12)
{
    Check(std::isfinite(actual) && std::isfinite(expected),"nonfinite "+name);
    const double error=std::abs(actual-expected)/std::max(1.,std::abs(expected));
    errors[name]=std::max(errors[name],error);
    if(error>tolerance)
    {
        std::ostringstream diagnostic;diagnostic<<std::setprecision(17)<<name<<": actual="<<actual
            <<", expected="<<expected<<", scaled error="<<error<<", tolerance="<<tolerance;
        Check(false,diagnostic.str());
    }
    else Check(true,name);
    return error;
}
double Near(const Vector &actual,const Vector &expected,const std::string &name,double tolerance=1e-12)
{
    Check(actual.size()==expected.size(),name+" shape mismatch");
    double maximum=0;
    for(std::size_t i=0;i<actual.size();++i)maximum=std::max(maximum,Near(actual[i],expected[i],name,tolerance));
    return maximum;
}
template<class Function> void Reject(const std::string &name,Function function)
{
    bool rejected=false;std::string message;
    try{function();}catch(const std::exception &error){rejected=true;message=error.what();}
    Check(rejected,"invalid input accepted: "+name);
    rejections.push_back({name,message});
}
Vector Forcing(double time)
{
    return {.5+std::sin(time),.2+.3*std::cos(time)};
}
Vector Exact(double time)
{
    // y_p=c/a+(a*sin(t)-cos(t))/(a^2+1) for g=c+sin(t).
    // y_p=c/a+b*(a*cos(t)+sin(t))/(a^2+1) for g=c+b*cos(t).
    const double p0=.25+(2*std::sin(time)-std::cos(time))/5,
        p1=.2/.7+.3*(.7*std::cos(time)+std::sin(time))/1.49;
    return {p0+(initial[0]-(.25-1./5))*std::exp(-2*time),
            p1+(initial[1]-(.2/.7+.21/1.49))*std::exp(-.7*time)};
}
struct Stage
{
    std::size_t step;
    double time,lambda;
    Vector history,state;
};
struct ExplicitEvent
{
    double time;
    Vector state,rhs;
};
struct Trial
{
    std::string name;
    std::size_t order=0,nextStep=0;
    double dt=0;
    bool complete=false;
    Vector finalState,exactState;
    double rmsError=0,maxStageError=0;
    std::vector<Stage> stages;
    std::vector<ExplicitEvent> explicitEvents;
    std::vector<Snapshot> snapshots;
    std::vector<Vector> recurrence;
};
std::unique_ptr<Engine> MakeEngine(Trial &trial,std::size_t order,double dt)
{
    trial.order=order;trial.dt=dt;
    return std::make_unique<Engine>(order,dt,
        [&trial](const Vector &history,double time,double lambda)
        {
            Check(history.size()==decay.size(),"native implicit input shape");
            Vector state(history.size());
            for(std::size_t i=0;i<state.size();++i)state[i]=history[i]/(1+decay[i]*lambda);
            trial.stages.push_back({trial.nextStep,time,lambda,history,state});
            return state;
        },
        [&trial](const Vector &state,double time)
        {
            Check(state.size()==decay.size(),"native explicit input shape");
            const auto rhs=Forcing(time);trial.explicitEvents.push_back({time,state,rhs});return rhs;
        });
}
Vector ExpectedHistory(std::size_t order,std::size_t step,double dt,
                       const Vector &current,const Vector &previous)
{
    const auto g=Forcing(step*dt),old=Forcing((step?step-1:0)*dt);
    Vector history(2);
    for(std::size_t i=0;i<history.size();++i)
        history[i]=order==1 || step==0 ? current[i]+dt*g[i] :
            4*current[i]/3-previous[i]/3+(2*dt/3)*(2*g[i]-old[i]);
    return history;
}
Vector Advance(Engine &engine,Trial &trial,std::size_t step,
               const Vector &current,const Vector &previous)
{
    trial.nextStep=step;const auto before=trial.stages.size();
    const auto actual=engine.Advance();
    Check(trial.stages.size()==before+1,"native IMEX did not execute exactly one coupled stage");
    const auto &stage=trial.stages.back();
    const auto history=ExpectedHistory(trial.order,step,trial.dt,current,previous);
    const double lambda=(trial.order==1 || step==0)?trial.dt:2*trial.dt/3;
    Vector expected(2);for(std::size_t i=0;i<expected.size();++i)expected[i]=history[i]/(1+decay[i]*lambda);
    trial.maxStageError=std::max({trial.maxStageError,
        Near(stage.history,history,"native_stage_combined_history"),
        Near(stage.time,(step+1)*trial.dt,"native_stage_time"),
        Near(stage.lambda,lambda,"native_stage_lambda"),
        Near(stage.state,expected,"native_implicit_solution"),
        Near(actual,expected,"native_timestep_output")});
    Near(engine.Time(),(step+1)*trial.dt,"native_wrapper_time");
    Check(engine.Steps()==step+1,"native global step count mismatch");
    return expected;
}
void HistoryNonzero(const Snapshot &snapshot)
{
    Check(snapshot.history.size()==2*snapshot.order,"native history slot count");
    double magnitude=0;
    for(std::size_t slot=snapshot.order;slot<snapshot.history.size();++slot)
        for(double value:snapshot.history[slot])magnitude=std::max(magnitude,std::abs(value));
    Check(magnitude>1e-3,"native explicit history is zero");
}
double Compare(const Snapshot &actual,const Snapshot &expected,const std::string &name)
{
    Check(actual.order==expected.order && actual.steps==expected.steps,name+" order/step mismatch");
    double error=std::max({Near(actual.time,expected.time,name+"_time"),
                           Near(actual.dt,expected.dt,name+"_dt"),
                           Near(actual.state,expected.state,name+"_state"),
                           Near(actual.times,expected.times,name+"_time_labels")});
    Check(actual.history.size()==expected.history.size(),name+" history shape");
    for(std::size_t slot=0;slot<actual.history.size();++slot)
        error=std::max(error,Near(actual.history[slot],expected.history[slot],name+"_every_native_slot"));
    HistoryNonzero(actual);return error;
}
void InitialHistory(const Snapshot &snapshot)
{
    Near(snapshot.state,initial,"initial_native_state");
    Near(snapshot.history[0],initial,"initial_native_solution_slot");
    const auto forcing=Forcing(0);Vector derivative(2);
    for(std::size_t i=0;i<derivative.size();++i)derivative[i]=snapshot.dt*forcing[i];
    Near(snapshot.history[snapshot.order],derivative,"initial_native_scaled_explicit_slot");
    if(snapshot.order==2)
    {
        Near(snapshot.history[1],Vector(2,0.),"unused_startup_solution_slot");
        Near(snapshot.history[3],Vector(2,0.),"unused_startup_explicit_slot");
    }
}
void Run(Trial &trial,std::size_t order,double dt)
{
    trial.name="order"+std::to_string(order)+"-dt"+std::to_string(dt);
    auto engine=MakeEngine(trial,order,dt);engine->Initialise(initial,0);
    trial.snapshots.push_back(engine->Capture());InitialHistory(trial.snapshots.back());
    trial.recurrence.push_back(initial);
    const std::size_t steps=static_cast<std::size_t>(std::llround(1./dt));
    Near(steps*dt,1.,"declared_final_time");
    for(std::size_t step=0;step<steps;++step)
    {
        const auto expected=Advance(*engine,trial,step,trial.recurrence[step],
            trial.recurrence[step?step-1:0]);
        trial.recurrence.push_back(expected);trial.snapshots.push_back(engine->Capture());
        HistoryNonzero(trial.snapshots.back());
    }
    trial.finalState=engine->DifferentialState();trial.exactState=Exact(1.);
    for(std::size_t i=0;i<trial.finalState.size();++i)
        trial.rmsError+=(trial.finalState[i]-trial.exactState[i])*(trial.finalState[i]-trial.exactState[i]);
    trial.rmsError=std::sqrt(trial.rmsError/trial.finalState.size());
    Check(trial.rmsError>0 && std::isfinite(trial.rmsError),"invalid temporal error");
    trial.complete=true;
}
struct Restart
{
    std::size_t order=0,split=0;
    double maxDifference=0;
    Snapshot splitState,finalState;
    Trial trial;
};
Restart Resume(const Trial &baseline,std::size_t split)
{
    Restart result;result.order=baseline.order;result.split=split;result.splitState=baseline.snapshots.at(split);
    result.trial.name="restart-order"+std::to_string(result.order)+"-split"+std::to_string(split);
    result.trial.nextStep=split;
    auto engine=MakeEngine(result.trial,result.order,baseline.dt);engine->Restore(result.splitState);
    result.maxDifference=Compare(engine->Capture(),result.splitState,"restored_split");
    for(std::size_t step=split;step<10;++step)
        Advance(*engine,result.trial,step,baseline.recurrence[step],baseline.recurrence[step?step-1:0]);
    result.finalState=engine->Capture();
    result.maxDifference=std::max(result.maxDifference,Compare(result.finalState,baseline.snapshots.back(),"restart_uninterrupted"));
    result.trial.finalState=engine->DifferentialState();result.trial.exactState=Exact(1.);result.trial.complete=true;
    return result;
}
std::filesystem::path Reserve(const std::filesystem::path &base)
{
    for(std::size_t generation=0;generation<1000000;++generation)
    {
        const std::filesystem::path candidate=base.string()+(generation?"."+std::to_string(generation):"");
        std::error_code error;
        if(std::filesystem::create_directory(candidate,error))return candidate;
        Check(!error || error==std::errc::file_exists,"cannot reserve probe checkpoint directory");
    }
    throw std::runtime_error("checkpoint directory generation limit");
}
std::string ReadBytes(const std::filesystem::path &path)
{
    std::ifstream stream(path,std::ios::binary);Check(static_cast<bool>(stream),"cannot open checkpoint control input");
    return {std::istreambuf_iterator<char>(stream),std::istreambuf_iterator<char>()};
}
void WriteBytes(const std::filesystem::path &path,const std::string &bytes)
{
    std::ofstream stream(path,std::ios::binary|std::ios::trunc);
    stream.write(bytes.data(),static_cast<std::streamsize>(bytes.size()));stream.close();
    Check(static_cast<bool>(stream),"cannot write checkpoint control");
}
Vector Coupled(const Vector &state)
{
    auto result=state;result.insert(result.end(),{-.17,.23,0.});return result;
}
struct CheckpointResult
{
    std::size_t order;
    std::string first,second;
    double continuationDifference=0;
    Snapshot saved,restored;
    Trial continuation;
};
CheckpointResult Checkpoints(const Trial &baseline,const std::filesystem::path &root)
{
    const auto directory=Reserve(root/("order"+std::to_string(baseline.order)));
    Engine::Checkpoint saved;
    saved.signature="M4 split ODE: a=2,0.7; g=0.5+sin(t),0.2+0.3*cos(t); D2,pw,po,gauge; eta=1";
    saved.integration=baseline.snapshots.at(3);saved.coupledState=Coupled(saved.integration.state);
    saved.fixedReference=Coupled(initial);saved.fixedDefect={.03,-.02,0.,0.,0.};saved.eta=1;
    CheckpointResult result;result.order=baseline.order;result.saved=saved.integration;
    const auto base=(directory/"state").string();result.first=Engine::WriteCheckpoint(base,saved);
    const auto payload=std::filesystem::path(result.first)/"State.bin";const auto original=ReadBytes(payload);
    result.second=Engine::WriteCheckpoint(base,saved);
    Check(result.first!=result.second,"checkpoint generation was overwritten");
    Check(ReadBytes(payload)==original,"publishing a new generation changed the first payload");
    const auto read=[&](const std::string &path,const std::string &signature,std::size_t order,double dt,
                        std::size_t coupled,std::size_t differential)
    {
        return Engine::ReadCheckpoint(path,signature,order,dt,coupled,differential);
    };
    const auto loaded=read(result.first,saved.signature,baseline.order,baseline.dt,5,2);
    result.restored=loaded.integration;
    Compare(loaded.integration,saved.integration,"checkpoint_native_history");
    Near(loaded.coupledState,saved.coupledState,"checkpoint_coupled_pressure_gauge");
    Near(loaded.fixedReference,saved.fixedReference,"checkpoint_fixed_reference");
    Near(loaded.fixedDefect,saved.fixedDefect,"checkpoint_fixed_free_defect");
    Near(loaded.eta,saved.eta,"checkpoint_eta");Check(loaded.signature==saved.signature,"checkpoint signature changed");
    result.continuation.name="checkpoint-continuation-order"+std::to_string(baseline.order);
    result.continuation.nextStep=3;
    auto engine=MakeEngine(result.continuation,baseline.order,baseline.dt);engine->Restore(loaded.integration);
    for(std::size_t step=3;step<10;++step)
        Advance(*engine,result.continuation,step,baseline.recurrence[step],baseline.recurrence[step-1]);
    result.continuationDifference=Compare(engine->Capture(),baseline.snapshots.back(),"checkpoint_continuation");
    result.continuation.finalState=engine->DifferentialState();result.continuation.complete=true;
    Reject("checkpoint_signature",[&]{read(result.first,saved.signature+"different",baseline.order,baseline.dt,5,2);});
    Reject("checkpoint_order",[&]{read(result.first,saved.signature,baseline.order==1?2:1,baseline.dt,5,2);});
    Reject("checkpoint_dt",[&]{read(result.first,saved.signature,baseline.order,baseline.dt*1.25,5,2);});
    Reject("checkpoint_coupled_shape",[&]{read(result.first,saved.signature,baseline.order,baseline.dt,6,2);});
    Reject("checkpoint_differential_shape",[&]{read(result.first,saved.signature,baseline.order,baseline.dt,5,3);});
    auto corrupted=original;Check(corrupted.size()>32,"checkpoint payload unexpectedly short");
    corrupted[corrupted.size()/2]^=1;
    const auto corruptDirectory=Reserve(directory/"corrupt");WriteBytes(corruptDirectory/"State.bin",corrupted);
    Reject("checkpoint_corrupt",[&]{read(corruptDirectory.string(),saved.signature,baseline.order,baseline.dt,5,2);});
    const auto truncatedDirectory=Reserve(directory/"truncated");WriteBytes(truncatedDirectory/"State.bin",original.substr(0,16));
    Reject("checkpoint_truncated",[&]{read(truncatedDirectory.string(),saved.signature,baseline.order,baseline.dt,5,2);});
    const auto incompleteDirectory=Reserve(directory/"incomplete");WriteBytes(incompleteDirectory/"State.bin.tmp",original);
    Reject("checkpoint_incomplete_publication",[&]{read(incompleteDirectory.string(),saved.signature,baseline.order,baseline.dt,5,2);});
    auto incomplete=saved.integration;incomplete.history.pop_back();
    Reject("snapshot_incomplete_history",[&]{engine->Restore(incomplete);});
    auto inconsistent=saved.integration;inconsistent.state[0]+=.1;
    Reject("snapshot_current_slot_mismatch",[&]{engine->Restore(inconsistent);});
    auto missingAnchor=saved;missingAnchor.fixedDefect.clear();
    Reject("checkpoint_missing_fixed_anchor",[&]{Engine::WriteCheckpoint((directory/"invalid-anchor").string(),missingAnchor);});
    Check(ReadBytes(payload)==original,"bad controls modified the positive checkpoint payload");
    return result;
}
void String(std::ostream &out,const std::string &value)
{
    out<<'"';
    for(unsigned char c:value)
    {
        if(c=='"')out<<"\\\"";else if(c=='\\')out<<"\\\\";
        else if(c=='\n')out<<"\\n";else if(c=='\r')out<<"\\r";else if(c=='\t')out<<"\\t";
        else if(c<32)out<<"\\u"<<std::hex<<std::setfill('0')<<std::setw(4)<<static_cast<int>(c)<<std::dec<<std::setfill(' ');
        else out<<static_cast<char>(c);
    }
    out<<'"';
}
void Values(std::ostream &out,const Vector &value)
{
    out<<'[';for(std::size_t i=0;i<value.size();++i){if(i)out<<',';out<<value[i];}out<<']';
}
void NativeHistory(std::ostream &out,const Snapshot &snapshot)
{
    out<<"{\"order\":"<<snapshot.order<<",\"steps\":"<<snapshot.steps
       <<",\"time\":"<<snapshot.time<<",\"dt\":"<<snapshot.dt<<",\"state\":";Values(out,snapshot.state);
    out<<",\"times\":";Values(out,snapshot.times);out<<",\"slots\":[";
    for(std::size_t slot=0;slot<snapshot.history.size();++slot){if(slot)out<<',';Values(out,snapshot.history[slot]);}
    out<<"]}";
}
void TrialReport(std::ostream &out,const Trial &trial,bool refinement)
{
    out<<"{\"name\":";String(out,trial.name);out<<",\"order\":"<<trial.order<<",\"dt\":"<<trial.dt
       <<",\"complete\":"<<(trial.complete?"true":"false")<<",\"max_stage_error\":"<<trial.maxStageError
       <<",\"final\":";Values(out,trial.finalState);
    if(refinement)
    {
        out<<",\"exact\":";Values(out,trial.exactState);out<<",\"rms_error\":"<<trial.rmsError;
        if(!trial.snapshots.empty()){out<<",\"initial_history\":";NativeHistory(out,trial.snapshots.front());}
    }
    out<<",\"implicit_stages\":[";
    for(std::size_t k=0;k<trial.stages.size();++k)
    {
        if(k){out<<',';}const auto &stage=trial.stages[k];
        out<<"{\"step\":"<<stage.step<<",\"time\":"<<stage.time<<",\"lambda\":"<<stage.lambda
           <<",\"history\":";Values(out,stage.history);out<<",\"state\":";Values(out,stage.state);out<<'}';
    }
    out<<"],\"explicit_evaluations\":[";
    for(std::size_t k=0;k<trial.explicitEvents.size();++k)
    {
        if(k){out<<',';}const auto &event=trial.explicitEvents[k];out<<"{\"time\":"<<event.time<<",\"state\":";
        Values(out,event.state);out<<",\"rhs\":";Values(out,event.rhs);out<<'}';
    }
    out<<"]}";
}
}

int main(int argc,char **argv)
{
    std::filesystem::path output,checkpointRoot;
    for(int argument=1;argument<argc;++argument)
    {
        const std::string option=argv[argument];
        if((option=="--output" || option=="--checkpoint-root") && argument+1<argc)
        {
            if(option=="--output")output=argv[++argument];else checkpointRoot=argv[++argument];
        }
        else
        {
            std::cerr<<"Usage: RANSEquilibriumTimeProbe --output REPORT.json [--checkpoint-root DIR]\n";return 2;
        }
    }
    if(output.empty()){std::cerr<<"--output is required\n";return 2;}
    std::vector<Trial> trials;std::vector<Restart> restarts;std::vector<CheckpointResult> checkpoints;
    std::map<std::size_t,Vector> rates;std::string error;bool passed=false;
    try
    {
        Near(Exact(0),initial,"analytic_initial_condition");
        if(checkpointRoot.empty())checkpointRoot=std::filesystem::temp_directory_path()/"nektar-wall-time";
        std::filesystem::create_directories(checkpointRoot);checkpointRoot=Reserve(checkpointRoot/"time-probe");
        for(std::size_t order:{1u,2u})
        {
            const std::size_t base=trials.size();
            for(double dt:{.1,.05,.025}){trials.emplace_back();Run(trials.back(),order,dt);}
            for(std::size_t level=0;level<2;++level)
            {
                const double measured=std::log(trials[base+level].rmsError/trials[base+level+1].rmsError)/std::log(2.);
                rates[order].push_back(measured);
                Check(order==1?(measured>.85 && measured<1.15):(measured>1.8 && measured<2.2),
                      "incorrect native IMEX temporal order");
            }
            for(std::size_t split:{0u,1u,2u,3u})restarts.push_back(Resume(trials[base],split));
            checkpoints.push_back(Checkpoints(trials[base],checkpointRoot));
        }
        passed=true;
    }
    catch(const std::exception &failure){error=failure.what();std::cerr<<error<<'\n';}
    try
    {
        if(!output.parent_path().empty())std::filesystem::create_directories(output.parent_path());
        std::ofstream out(output);out<<std::setprecision(17);
        out<<"{\"schema\":1,\"verification_status\":\""<<(passed?"passed":"failed")<<"\",\"checks\":"<<checks
           <<",\"error\":";String(out,error);out<<",\"checkpoint_root\":";String(out,checkpointRoot.string());
        out<<",\"initial\":";Values(out,initial);out<<",\"decay\":[2,0.7],\"final_time\":1,\"refinement\":[";
        for(std::size_t k=0;k<trials.size();++k){if(k)out<<',';TrialReport(out,trials[k],true);}
        out<<"],\"temporal_rates\":{";std::size_t entry=0;
        for(const auto &rate:rates){if(entry++)out<<',';String(out,std::to_string(rate.first));out<<':';Values(out,rate.second);}
        out<<"},\"restarts\":[";
        for(std::size_t k=0;k<restarts.size();++k)
        {
            if(k){out<<',';}const auto &restart=restarts[k];out<<"{\"order\":"<<restart.order<<",\"split\":"<<restart.split
               <<",\"max_difference\":"<<restart.maxDifference<<",\"split_history\":";NativeHistory(out,restart.splitState);
            out<<",\"final_history\":";NativeHistory(out,restart.finalState);out<<",\"trial\":";
            TrialReport(out,restart.trial,false);out<<'}';
        }
        out<<"],\"checkpoints\":[";
        for(std::size_t k=0;k<checkpoints.size();++k)
        {
            if(k){out<<',';}const auto &checkpoint=checkpoints[k];out<<"{\"order\":"<<checkpoint.order<<",\"first\":";
            String(out,checkpoint.first);out<<",\"second\":";String(out,checkpoint.second);
            out<<",\"continuation_difference\":"<<checkpoint.continuationDifference<<",\"saved_history\":";
            NativeHistory(out,checkpoint.saved);out<<",\"loaded_history\":";NativeHistory(out,checkpoint.restored);
            out<<",\"trial\":";TrialReport(out,checkpoint.continuation,false);out<<'}';
        }
        out<<"],\"scaled_errors\":{";entry=0;
        for(const auto &item:errors){if(entry++)out<<',';String(out,item.first);out<<':'<<item.second;}
        out<<"},\"rejected_controls\":[";
        for(std::size_t k=0;k<rejections.size();++k)
        {
            if(k){out<<',';}out<<"{\"name\":";String(out,rejections[k].name);out<<",\"message\":";
            String(out,rejections[k].message);out<<'}';
        }
        out<<"]}\n";out.close();
        if(!out)throw std::runtime_error("cannot write time-probe report");
    }
    catch(const std::exception &failure){std::cerr<<failure.what()<<'\n';return 1;}
    if(passed)std::cout<<"Native wall IMEX verified: "<<checks<<" checks, "<<rejections.size()<<" rejected controls.\n";
    return passed?0:1;
}
