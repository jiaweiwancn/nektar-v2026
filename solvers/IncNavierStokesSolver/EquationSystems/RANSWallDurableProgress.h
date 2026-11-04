///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_DURABLE_PROGRESS_H
#define NEKTAR_RANS_WALL_DURABLE_PROGRESS_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <array>
#include <cerrno>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <fcntl.h>
#include <unistd.h>

namespace Nektar::RANSWallR3
{
inline void Check(bool value,const std::string &message)
{if(!value)throw std::runtime_error("R3 durable progress: "+message);}

/** Small streaming SHA256 for scientific payload identities, not encryption. */
class Sha256
{
public:
    void Update(const void *data,std::size_t size)
    {
        const auto *bytes=static_cast<const unsigned char*>(data);m_bytes+=size;
        while(size)
        {
            const auto take=std::min(size,64-m_used);
            std::memcpy(m_block.data()+m_used,bytes,take);m_used+=take;bytes+=take;size-=take;
            if(m_used==64){Transform();m_used=0;}
        }
    }
    std::string Finish() const
    {
        auto value=*this;const std::uint64_t bits=m_bytes*8;
        const unsigned char first=0x80,zero=0;value.Update(&first,1);
        while(value.m_used!=56)value.Update(&zero,1);
        std::array<unsigned char,8> length{};
        for(int i=0;i<8;++i)length[7-i]=static_cast<unsigned char>(bits>>(8*i));
        value.Update(length.data(),length.size());std::ostringstream out;out<<std::hex<<std::setfill('0');
        for(auto word:value.m_state)out<<std::setw(8)<<word;
        return out.str();
    }
private:
    static std::uint32_t Rotate(std::uint32_t x,int n){return (x>>n)|(x<<(32-n));}
    void Transform()
    {
        static constexpr std::array<std::uint32_t,64> constants{
            0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
            0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
            0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
            0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
            0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
            0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
            0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
            0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
        std::array<std::uint32_t,64> words{};
        for(int i=0;i<16;++i)for(int b=0;b<4;++b)words[i]=(words[i]<<8)|m_block[4*i+b];
        for(int i=16;i<64;++i)
        {
            const auto a=Rotate(words[i-15],7)^Rotate(words[i-15],18)^(words[i-15]>>3);
            const auto b=Rotate(words[i-2],17)^Rotate(words[i-2],19)^(words[i-2]>>10);
            words[i]=words[i-16]+a+words[i-7]+b;
        }
        auto a=m_state[0],b=m_state[1],c=m_state[2],d=m_state[3];
        auto e=m_state[4],f=m_state[5],g=m_state[6],h=m_state[7];
        for(int i=0;i<64;++i)
        {
            const auto s1=Rotate(e,6)^Rotate(e,11)^Rotate(e,25);
            const auto t1=h+s1+((e&f)^((~e)&g))+constants[i]+words[i];
            const auto s0=Rotate(a,2)^Rotate(a,13)^Rotate(a,22);
            const auto t2=s0+((a&b)^(a&c)^(b&c));
            h=g;g=f;f=e;e=d+t1;d=c;c=b;b=a;a=t1+t2;
        }
        m_state[0]+=a;m_state[1]+=b;m_state[2]+=c;m_state[3]+=d;
        m_state[4]+=e;m_state[5]+=f;m_state[6]+=g;m_state[7]+=h;
    }
    std::array<std::uint32_t,8> m_state{0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,
        0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    std::array<unsigned char,64> m_block{};
    std::size_t m_used=0;std::uint64_t m_bytes=0;
};
inline std::string Hash(const void *data,std::size_t bytes)
{Sha256 value;value.Update(data,bytes);return value.Finish();}
inline std::string Hash(const std::string &value){return Hash(value.data(),value.size());}
inline std::string Hash(const std::vector<double> &value){return Hash(value.data(),value.size()*sizeof(double));}
inline std::string HashFile(const std::string &path)
{
    std::ifstream in(path,std::ios::binary);Check(in.good(),"cannot hash "+path);
    Sha256 hash;std::array<char,65536> bytes{};
    while(in){in.read(bytes.data(),bytes.size());hash.Update(bytes.data(),in.gcount());}
    Check(in.eof(),"hash read failed "+path);return hash.Finish();
}
inline void WriteAll(int fd,const void *data,std::size_t bytes)
{
    const auto *at=static_cast<const char*>(data);
    while(bytes)
    {
        const auto done=::write(fd,at,bytes);
        if(done<0 && errno==EINTR)continue;
        Check(done>0,"payload/event write failed");at+=done;bytes-=done;
    }
}
inline void SyncDirectory(const std::filesystem::path &file)
{
    const auto parent=file.parent_path();const int fd=::open(parent.c_str(),O_RDONLY|O_DIRECTORY);
    Check(fd>=0,"cannot open payload directory");const int result=::fsync(fd);::close(fd);
    Check(result==0,"directory fsync failed");
}
inline std::string Pin(const std::string &path)
{
    const auto absolute=std::filesystem::absolute(path).string();
    return RANSWallSolveJson::Object({{"path",RANSWallSolveJson::Quote(absolute)},
        {"bytes",RANSWallSolveJson::Number(std::filesystem::file_size(absolute))},
        {"sha256",RANSWallSolveJson::Quote(HashFile(absolute))}});
}

/** R3-only immutable payloads, flushed event stream and atomic last snapshot. */
class DurableProgress
{
public:
    DurableProgress(const std::string &path,int solve,double alpha,double time,const std::vector<double> &history,
        const std::string &action):m_path(std::filesystem::absolute(path).string()),m_action(action),m_solve(solve),m_alpha(alpha),m_time(time)
    {
        std::uint16_t endian=1;
        Check(sizeof(double)==8 && std::numeric_limits<double>::is_iec559 &&
            *reinterpret_cast<unsigned char*>(&endian)==1,"requires little-endian IEEE754 binary64");
        std::filesystem::create_directories(std::filesystem::path(m_path).parent_path());
        m_fd=::open(m_path.c_str(),solve==1?O_WRONLY|O_CREAT|O_EXCL:O_WRONLY|O_APPEND,0644);
        Check(m_fd>=0,"refusing an existing first-solve progress stream or missing later stream");
        m_history=VectorPayload("history",history);
        Event("solve_started",{{"history",m_history},{"history_sha256",RANSWallSolveJson::Quote(Hash(history))},
            {"history_advanced","false"},{"payload_encoding",RANSWallSolveJson::Quote("little-endian IEEE754 binary64")}});
    }
    ~DurableProgress(){if(m_fd>=0)::close(m_fd);}
    DurableProgress(const DurableProgress&)=delete;
    DurableProgress &operator=(const DurableProgress&)=delete;

    std::string VectorPayload(const std::string &name,const std::vector<double> &value)
    {
        const auto path=Base()+"."+name+".bin";
        Atomic(path,[&](int fd){WriteAll(fd,value.data(),value.size()*sizeof(double));},false);
        return Pin(path);
    }
    std::string MatrixPayload(int iteration,const std::vector<std::vector<double>> &matrix)
    {
        const auto path=Base()+".J"+std::to_string(iteration)+".bin";
        Atomic(path,[&](int fd){for(const auto &row:matrix)WriteAll(fd,row.data(),row.size()*sizeof(double));},false);
        return Pin(path);
    }
    std::string Snapshot(int iteration,const std::vector<double> &state,const std::vector<double> &residual,
        const std::map<int,double> &constraints,const std::vector<double> &norms,const std::string &domain,
        const std::string &profiling,const std::string &arrival="null")
    {
        namespace J=RANSWallSolveJson;
        const auto prefix="iterate"+std::to_string(iteration);std::vector<double> pairs;
        for(auto [row,value]:constraints){pairs.push_back(row);pairs.push_back(value);}
        const auto statePin=VectorPayload(prefix+".state",state);
        const auto residualPin=VectorPayload(prefix+".residual",residual);
        const auto constraintsPin=VectorPayload(prefix+".constraints",pairs);
        const auto value=J::Object({{"schema","1"},{"policy",J::Quote("ScalarFeasibleR3:v1")},
            {"action_policy",J::Quote(m_action)},{"solve_number",J::Number(m_solve)},
            {"iteration",J::Number(iteration)},{"alpha",J::Number(m_alpha)},{"time",J::Number(m_time)},
            {"state",statePin},{"history",m_history},{"residual",residualPin},{"constraints",constraintsPin},
            {"constraint_pairs",J::Number(constraints.size())},{"scaled_block_norms",J::Numbers(norms)},
            {"scalar_domain",domain},{"profiling",profiling},{"arrival_update",arrival},
            {"is_complete_internal_iterate","true"},{"is_physical_time_sample","false"},{"history_advanced","false"}});
        const auto file=Base()+"."+prefix+".snapshot.json";Text(file,value,false);
        Text(m_path+".last.json",value,true);
        Event("iterate_completed",{{"iteration",J::Number(iteration)},{"snapshot",Pin(file)}});
        return Pin(file);
    }
    void Event(const std::string &kind,std::vector<std::pair<std::string,std::string>> fields={})
    {
        namespace J=RANSWallSolveJson;
        fields.insert(fields.begin(),{{"schema","1"},{"policy",J::Quote("ScalarFeasibleR3:v1")},
            {"action_policy",J::Quote(m_action)},{"solve_number",J::Number(m_solve)},
            {"sequence",J::Number(++m_sequence)},{"kind",J::Quote(kind)},
            {"alpha",J::Number(m_alpha)},{"time",J::Number(m_time)}});
        const auto value=J::Object(fields)+"\n";WriteAll(m_fd,value.data(),value.size());
        Check(::fsync(m_fd)==0,"event fsync failed");
    }
    void JacobianStatus(int iteration,const std::string &status,const std::string &matrix="null",
        const std::string &details="null")
    {
        namespace J=RANSWallSolveJson;
        const auto value=J::Object({{"iteration",J::Number(iteration)},{"status",J::Quote(status)},
            {"action_policy",J::Quote(m_action)},{"matrix",matrix},{"details",details}});
        Text(Base()+".J"+std::to_string(iteration)+".status.json",value,true);
        Event(status=="complete"?"jacobian_completed":"jacobian_started_incomplete",
            {{"iteration",J::Number(iteration)},{"matrix",matrix},{"details",details}});
    }
private:
    std::string Base() const{return m_path+".solve"+std::to_string(m_solve);}
    template<class Write> void Atomic(const std::string &path,Write write,bool replace)
    {
        Check(replace || !std::filesystem::exists(path),"refusing existing immutable payload "+path);
        const auto temporary=path+".tmp."+std::to_string(::getpid())+"."+std::to_string(++m_temporary);
        const int fd=::open(temporary.c_str(),O_WRONLY|O_CREAT|O_EXCL,0644);Check(fd>=0,"cannot reserve payload temporary");
        try{write(fd);Check(::fsync(fd)==0,"payload fsync failed");}
        catch(...){::close(fd);throw;}
        Check(::close(fd)==0,"payload close failed");std::filesystem::rename(temporary,path);SyncDirectory(path);
    }
    void Text(const std::string &path,const std::string &value,bool replace)
    {const auto bytes=value+"\n";Atomic(path,[&](int fd){WriteAll(fd,bytes.data(),bytes.size());},replace);}
    std::string m_path,m_action,m_history;int m_fd=-1,m_solve=0,m_sequence=0,m_temporary=0;
    double m_alpha,m_time;
};
} // namespace Nektar::RANSWallR3
#endif
