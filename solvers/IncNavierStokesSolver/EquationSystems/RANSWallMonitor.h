// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_MONITOR_H
#define NEKTAR_RANS_WALL_MONITOR_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <boost/property_tree/json_parser.hpp>
#include <chrono>
#include <csignal>
#include <list>
#include <set>
#include <thread>

namespace Nektar::RANSWallMonitor
{
namespace J=RANSWallSolveJson;
namespace D=RANSWallR3;
using Tree=boost::property_tree::ptree;
inline Tree Read(const std::string &p){Tree t;boost::property_tree::read_json(p,t);return t;}
inline std::string Bits(double x)
{std::uint64_t u;std::memcpy(&u,&x,sizeof u);std::ostringstream s;s<<std::hex<<std::setfill('0')<<std::setw(16)<<u;return s.str();}
inline void Atomic(const std::string &path,const std::string &text,bool immutable=true)
{
    const std::filesystem::path p(path);std::filesystem::create_directories(p.parent_path());
    D::Check(!immutable || !std::filesystem::exists(p),"immutable monitor output exists: "+path);
    const auto tmp=path+".tmp-"+std::to_string(getpid());
    const int fd=::open(tmp.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);D::Check(fd>=0,"monitor temporary open failed");
    try{D::WriteAll(fd,text.data(),text.size());D::Check(::fsync(fd)==0,"monitor fsync failed");::close(fd);}
    catch(...){::close(fd);throw;}
    if(immutable){D::Check(::link(tmp.c_str(),path.c_str())==0,"immutable monitor publish failed");::unlink(tmp.c_str());}
    else D::Check(::rename(tmp.c_str(),path.c_str())==0,"monitor rename failed");
    D::SyncDirectory(p);
}
inline void Append(const std::string &path,const std::string &text)
{
    const int fd=::open(path.c_str(),O_WRONLY|O_APPEND|O_CREAT,0600);D::Check(fd>=0,"monitor stream open failed");
    try{D::WriteAll(fd,text.data(),text.size());D::WriteAll(fd,"\n",1);D::Check(::fsync(fd)==0,"monitor stream fsync failed");::close(fd);}
    catch(...){::close(fd);throw;}
}
// State.bin has already been closed/fsynced by the caller. D::SyncDirectory
// syncs its argument's parent: persist the checkpoint contents first, then
// the newly created checkpoint directory entry in its containing parent.
template<class Sync> inline void CheckpointDirectories(const std::string &checkpoint,Sync sync)
{
    auto directory=std::filesystem::path(checkpoint);
    if(directory.filename().empty())directory=directory.parent_path();
    sync(directory/"State.bin");sync(directory);
}
inline void CheckpointDirectories(const std::string &checkpoint)
{CheckpointDirectories(checkpoint,[](const std::filesystem::path &p){D::SyncDirectory(p);});}
inline void SameObject(const Tree &a,const Tree &b)
{
    D::Check(a.size()==b.size(),"monitor binding shape differs");std::set<std::string> seen;
    for(const auto &v:a)
    {
        D::Check(seen.insert(v.first).second && v.second.empty(),"duplicate/nonleaf monitor binding");
        D::Check(b.count(v.first)==1 && b.get<std::string>(v.first)==v.second.get_value<std::string>(),"monitor binding differs: "+v.first);
    }
}
inline void CheckPin(const Tree &t)
{
    const auto p=t.get<std::string>("path");D::Check(std::filesystem::file_size(p)==t.get<std::uint64_t>("bytes") &&
        D::HashFile(p)==t.get<std::string>("sha256"),"monitor pinned bytes differ");
}
inline std::string Wait(const std::string &directory,const std::string &stem,const std::string &readyHash,
    const std::string &binding,double seconds,int observerPid)
{
    const auto begin=std::chrono::steady_clock::now();const auto expected=Read(directory+"/"+stem+".ready.json").get_child("binding");
    D::Check(expected.size()==17,"incomplete native monitor binding");
    Tree declared;std::istringstream in(binding);boost::property_tree::read_json(in,declared);SameObject(expected,declared);
    const auto ack=directory+"/"+stem+".ack.json";
    while(true)
    {
        D::Check(::kill(observerPid,0)==0,"observer process exited before qualification");
        for(const auto &p:std::filesystem::directory_iterator(directory))
        {
            // Only the final atomic name is a published ACK. During link /
            // unlink publication both it and .ack.json.tmp-<pid> can exist;
            // before publication the temporary alone must never authorize or
            // reject a step. Consumed ACKs are separate durable evidence.
            const auto name=p.path().filename().string();const std::string suffix=".ack.json";
            if(name.size()>=suffix.size() && name.compare(name.size()-suffix.size(),suffix.size(),suffix)==0)
                D::Check(p.path()==std::filesystem::path(ack),"stale/duplicate observer ACK");
        }
        if(std::filesystem::exists(ack))
        {
            const auto a=Read(ack);const std::set<std::string> keys{"schema","status","ready_sha256","binding","result","qualified_cursor","continuation","inventory"};
            D::Check(a.size()==keys.size(),"unexpected/duplicate observer ACK fields");
            for(const auto &key:keys)D::Check(a.count(key)==1,"missing/duplicate observer ACK field");
            D::Check(a.count("schema")==1 && a.count("binding")==1 && a.count("status")==1 &&
                a.get<std::string>("schema")=="WJObserverAckV1" && a.get<std::string>("ready_sha256")==readyHash,"malformed/mismatched observer ACK");
            SameObject(expected,a.get_child("binding"));CheckPin(a.get_child("result"));
            CheckPin(a.get_child("inventory"));
            const auto result=Read(a.get<std::string>("result.path"));SameObject(expected,result.get_child("binding"));
            D::Check(result.get<std::string>("inventory_sha256")==a.get<std::string>("inventory.sha256"),"result actual inventory differs");
            D::Check(a.get<std::string>("status")=="PASS" && result.get<std::string>("status")=="PASS","independent observer qualification failed");
            CheckPin(a.get_child("qualified_cursor"));const auto cursor=Read(a.get<std::string>("qualified_cursor.path"));
            SameObject(expected,cursor.get_child("binding"));D::Check(cursor.get<std::string>("ready_sha256")==readyHash,"qualified cursor ready differs");
            D::Check(cursor.get<std::string>("inventory_sha256")==a.get<std::string>("inventory.sha256"),"cursor actual inventory differs");
            const auto next=a.get<std::string>("continuation");D::Check(next=="PROCEED" || next=="TERMINAL" || next=="RESOURCE_LIMIT","invalid continuation decision");
            D::Check(::rename(ack.c_str(),(directory+"/"+stem+".ack.consumed.json").c_str())==0,"cannot consume observer ACK");
            D::SyncDirectory(ack);return next;
        }
        D::Check(std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count()<seconds,"observer ACK timeout");
        std::this_thread::sleep_for(std::chrono::milliseconds(25));
    }
}

/** Execution-only bounded catalog. Eviction cannot suppress a call/definition. */
class FactorCatalog
{
    struct Entry {std::string key,id;};
    std::list<Entry> entries;std::size_t bytes=0;
    std::string directory,signature;std::uint64_t calls=0,definitions=0,evictions=0;
public:
    bool Enabled()const{return !directory.empty();}
    const std::string &Directory()const{return directory;}
    void Configure(const std::string &p,const std::string &s)
    {D::Check(!Enabled(),"factor catalog configured twice");directory=p;signature=s;std::filesystem::create_directories(p);}
    void Record(const std::string &pattern,const std::string &metadata)
    {
        const std::string key=signature+pattern;auto it=std::find_if(entries.begin(),entries.end(),[&](const Entry &e){return e.key==key;});std::string id;
        if(it!=entries.end()){id=it->id;entries.splice(entries.begin(),entries,it);}
        else
        {
            id="factor-"+std::to_string(++definitions);const auto path=directory+"/"+id+".pattern.bin";Atomic(path,pattern);
            Atomic(directory+"/"+id+".json",J::Object({{"schema",J::Quote("WJFactorDefinitionV1")},
                {"id",J::Quote(id)},{"scope_sha256",J::Quote(signature)},{"exact_pattern",D::Pin(path)},
                {"pattern_encoding",J::Quote("opaque exact structural key; native matrix metadata specifies its encoding; coefficients intentionally excluded")},{"metadata",metadata}}));
            // Huge patterns are still streamed, but are never resident catalog entries.
            while(!entries.empty() && (entries.size()>=16 || bytes+key.size()>32*1024*1024))
            {bytes-=entries.back().key.size();entries.pop_back();++evictions;}
            if(key.size()<=32*1024*1024){entries.push_front({key,id});bytes+=key.size();}
        }
        Append(directory+"/calls.jsonl",J::Object({{"call",J::Number(++calls)},{"definition",J::Quote(id)},
            {"native_step",J::Number(RANSWallTiming::Get().step)},{"native_time_bits",J::Quote(Bits(RANSWallTiming::Get().time))},
            {"timing_scope_id",J::Number(RANSWallTiming::Get().stack.empty()?0:RANSWallTiming::Get().stack.back()->Id())},
            {"resident_definitions",J::Number(entries.size())},{"resident_key_bytes",J::Number(bytes)},{"evictions",J::Number(evictions)}}));
    }
    std::string Summary()const{return J::Object({{"enabled",J::Boolean(Enabled())},{"calls",J::Number(calls)},
        {"definitions",J::Number(definitions)},{"resident_definitions",J::Number(entries.size())},{"resident_key_bytes",J::Number(bytes)},
        {"evictions",J::Number(evictions)},{"maximum_definitions","16"},{"maximum_key_bytes","33554432"}});}
};
inline FactorCatalog &Factors(){static FactorCatalog c;return c;}
}
#endif
