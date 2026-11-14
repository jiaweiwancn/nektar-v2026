// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_TIMING_H
#define NEKTAR_RANS_WALL_TIMING_H
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <map>
#include <vector>
#include <unistd.h>

namespace Nektar::RANSWallTiming
{
// Opt-in serial-pilot telemetry, outside every physical/session signature.
// No clock or I/O calls occur in disabled scopes after the first env lookup.
// Inclusive and exclusive intervals must not be added to each other.
class Scope;
struct State
{
    bool enabled=false;std::ofstream out;std::vector<Scope*> stack;
    long next=0,step=-1;double time=-1;
    State()
    {
        const char *path=std::getenv("NEKTAR_RANS_TIMING_FILE");
        if(!path || !*path)return;
        out.open(path,std::ios::app);
        if(!out)throw std::runtime_error("cannot open NEKTAR_RANS_TIMING_FILE");
        enabled=true;out<<std::setprecision(17);
    }
};
inline State &Get(){static State state;return state;}
class Context
{
    State &s;long previousStep;double previousTime;
public:
    Context(long step,double time):s(Get()),previousStep(s.step),previousTime(s.time){s.step=step;s.time=time;}
    ~Context(){s.step=previousStep;s.time=previousTime;}
};
class Scope
{
    using Clock=std::chrono::steady_clock;
    State *s=nullptr;Scope *parent=nullptr;std::string name;
    long id=0,step=-1;double time=-1,childrenWall=0,childrenCpu=0;
    Clock::time_point begin;std::clock_t cpu;
public:
    explicit Scope(const std::string &label)
    {
        auto &state=Get();if(!state.enabled)return;
        s=&state;name=label;id=++s->next;step=s->step;time=s->time;
        parent=s->stack.empty()?nullptr:s->stack.back();s->stack.push_back(this);
        begin=Clock::now();cpu=std::clock();
    }
    Scope(const Scope&)=delete;Scope &operator=(const Scope&)=delete;
    long Id()const{return id;}
    const std::string &Name()const{return name;}
    ~Scope(){Stop();}
    void Stop() noexcept
    {
        if(!s)return;
        const double wall=std::chrono::duration<double>(Clock::now()-begin).count();
        const double processCpu=double(std::clock()-cpu)/CLOCKS_PER_SEC;
        if(s->stack.empty() || s->stack.back()!=this)std::terminate();
        s->stack.pop_back();
        if(parent){parent->childrenWall+=wall;parent->childrenCpu+=processCpu;}
        s->out<<"{\"schema\":1,\"pid\":"<<getpid()<<",\"id\":"<<id<<",\"parent_id\":"<<(parent?parent->id:0)
            <<",\"step\":"<<step<<",\"time\":"<<time<<",\"name\":"<<std::quoted(name)
            <<",\"wall_seconds\":"<<wall<<",\"process_cpu_seconds\":"<<processCpu
            <<",\"exclusive_wall_seconds\":"<<std::max(0.,wall-childrenWall)
            <<",\"exclusive_process_cpu_seconds\":"<<std::max(0.,processCpu-childrenCpu)<<"}\n";
        s=nullptr;
    }
};
// Per-point audit buckets accumulate clocks in memory and serialize once.
// Buckets are disjoint wall intervals; parent Scope retains the unbucketed work.
class Accumulator
{
    std::map<std::string,std::pair<long,double>> values;
public:
    class Part
    {
        Accumulator *owner;std::string name;std::chrono::steady_clock::time_point began;
    public:
        Part(Accumulator &a,const std::string &n):owner(Get().enabled?&a:nullptr),name(n)
        {if(owner)began=std::chrono::steady_clock::now();}
        void Stop(){if(!owner)return;auto &v=owner->values[name];++v.first;v.second+=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();owner=nullptr;}
        ~Part(){Stop();}
    };
    ~Accumulator()
    {
        auto &s=Get();if(!s.enabled)return;
        for(const auto &[name,v]:values)s.out<<"{\"schema\":1,\"pid\":"<<getpid()<<",\"id\":"<<++s.next<<",\"parent_id\":0,\"step\":"<<s.step
            <<",\"time\":"<<s.time<<",\"name\":"<<std::quoted(name)<<",\"calls\":"<<v.first<<",\"wall_seconds\":"<<v.second
            <<",\"exclusive_wall_seconds\":"<<v.second<<",\"process_cpu_seconds\":0,\"exclusive_process_cpu_seconds\":0,\"clock_scope\":\"accumulated disjoint wall only; CPU not measured\"}\n";
    }
};
}
#endif
