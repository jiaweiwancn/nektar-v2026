// SPDX-License-Identifier: MIT
#ifndef WJ_H26_ADVANCE_DRIVER_H
#define WJ_H26_ADVANCE_DRIVER_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <SolverUtils/Turbulence/AnalyticalWallFunction.h>
#include <numeric>
namespace H26
{
namespace Durable=Nektar::RANSWallR3;
struct MMS
{
    double xl=0,zl=0,Lx=1,Lz=1,nu=6.423173803526447e-5,muRef=.01,stress=.001,zv=0;
    std::string solution="polynomial",model="time",profile="linear";
    Geometry geo;
    std::array<double,4> Shape(double a)const
    {
        if(solution=="trigonometric")return {std::pow(std::sin(M_PI*a),2),M_PI*std::sin(2*M_PI*a),2*M_PI*M_PI*std::cos(2*M_PI*a),-4*std::pow(M_PI,3)*std::sin(2*M_PI*a)};
        return {a*a*(1-a),2*a-3*a*a,2-6*a,-6};
    }
    double Modulation(double t)const{return model=="time"?1+.1*std::sin(M_PI*t):1.;}
    std::array<double,3> Turbulent(double x,double z)const
    {
        if(profile=="linear")return {.01*(1+(x-xl)/Lx+(z-zl)/Lz),.01/Lx,.01/Lz};
        double r=z-geo.Ground(x),s=std::clamp((r/geo.H-.4)/.35,0.,1.),chi=s*s*s*(10+s*(-15+6*s)),dchi=(s>0&&s<1)?30*s*s*(s-1)*(s-1)/(.35*geo.H):0;
        double slope=.3819178635467331*.051,inner=slope*std::max(r-zv,0.),outer=slope*std::max(r,0.);
        double value=(1-chi)*inner+chi*outer,dr=(1-chi)*(r>zv?slope:0.)+chi*(r>0?slope:0.)+dchi*(outer-inner);
        return {value,-geo.Slope(x)*dr,dr};
    }
    std::array<double,3> Mu(double x,double z,double t)const
    {
        if(model=="nu")return {nu,0,0};if(model=="constant")return {muRef,0,0};auto v=Turbulent(x,z);double a=Modulation(t);return {nu+a*v[0],a*v[1],a*v[2]};
    }
    std::array<double,3> Value(double x,double z,double t)const
    {
        double a=std::exp(-t);auto f=Shape((x-xl)/Lx),g=Shape((z-zl)/Lz);return {a*f[0]*g[1],-a*Lz/Lx*f[1]*g[0],.01*std::exp(-2*t)*(std::pow((x-xl)/Lx,2)+std::pow((z-zl)/Lz,2))};
    }
    Tensor Gradient(double x,double z,double t)const
    {
        double a=std::exp(-t);auto f=Shape((x-xl)/Lx),g=Shape((z-zl)/Lz);Tensor G{};G[0][0]=a*f[1]*g[1]/Lx;G[0][2]=a*f[0]*g[2]/Lz;G[2][0]=-a*Lz/(Lx*Lx)*f[2]*g[0];G[2][2]=-G[0][0];return G;
    }
    Tensor Remainder(double x,double z,double t)const
    {
        Tensor T{};double X=(x-xl)/Lx,Z=(z-zl)/Lz,a=stress*std::exp(-t);T[0][0]=a*X*Z;T[2][2]=-T[0][0];T[0][2]=T[2][0]=a*(X+Z)/2;return T;
    }
    Vector Force(double x,double z,double t)const
    {
        auto u=Value(x,z,t);auto G=Gradient(x,z,t);auto mu=Mu(x,z,t);auto f=Shape((x-xl)/Lx),g=Shape((z-zl)/Lz);double a=std::exp(-t),b=.01*std::exp(-2*t),X=(x-xl)/Lx,Z=(z-zl)/Lz;
        Vector lap{a*(f[2]*g[1]/(Lx*Lx)+f[0]*g[3]/(Lz*Lz)),0,-a*Lz/Lx*(f[3]*g[0]/(Lx*Lx)+f[1]*g[2]/(Lz*Lz))};
        Vector result{};result[0]=-u[0]+u[0]*G[0][0]+u[1]*G[0][2]+2*b*X/Lx-mu[0]*lap[0]-2*G[0][0]*mu[1]-(G[0][2]+G[2][0])*mu[2]-stress*a*(Z/Lx+.5/Lz);
        result[2]=-u[1]+u[0]*G[2][0]+u[1]*G[2][2]+2*b*Z/Lz-mu[0]*lap[2]-(G[0][2]+G[2][0])*mu[1]-2*G[2][2]*mu[2]-stress*a*(.5/Lx-X/Lz);return result;
    }
    Vector PressureGradient(double x,double z,double t)const{return {.02*std::exp(-2*t)*(x-xl)/(Lx*Lx),0,.02*std::exp(-2*t)*(z-zl)/(Lz*Lz)};}
    std::array<Function,3> Functions(double t,double pressureMean=0)const
    {
        return {[=](double x,double z){return Value(x,z,t)[0];},[=](double x,double z){return Value(x,z,t)[1];},[=](double x,double z){return Value(x,z,t)[2]-pressureMean;}};
    }
};
struct History
{
    std::uint64_t step=0;double time=0,dt=0;
    Vec current,previous,convection,previousConvection;
};
std::string StateHash(const Vec&v){return Durable::Hash(v.data(),v.size()*sizeof(double));}
void WriteHistory(const fs::path&path,const std::string&signature,const History&h)
{
    if(fs::exists(path))throw std::runtime_error("immutable H26 checkpoint already exists");std::ostringstream s(std::ios::binary);s<<"H26_NATIVE_BDF2_V1\n"<<signature<<'\n';s.write(reinterpret_cast<const char*>(&h.step),sizeof(h.step));s.write(reinterpret_cast<const char*>(&h.time),sizeof(double));s.write(reinterpret_cast<const char*>(&h.dt),sizeof(double));
    for(const auto*v:{&h.current,&h.previous,&h.convection,&h.previousConvection}){std::uint64_t n=v->size();s.write(reinterpret_cast<const char*>(&n),sizeof(n));s.write(reinterpret_cast<const char*>(v->data()),n*sizeof(double));}std::string payload=s.str(),digest=Durable::Hash(payload);payload+=digest;
    fs::path tmp=path.string()+".tmp";int fd=::open(tmp.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);Durable::Check(fd>=0,"checkpoint create");Durable::WriteAll(fd,payload.data(),payload.size());Durable::Check(::fsync(fd)==0,"checkpoint sync");::close(fd);fs::rename(tmp,path);Durable::SyncDirectory(path);
}
History ReadHistory(const fs::path&path,const std::string&signature,int size,double dt)
{
    std::ifstream file(path,std::ios::binary);std::string data((std::istreambuf_iterator<char>(file)),{});Durable::Check(data.size()>64,"short checkpoint");std::string payload=data.substr(0,data.size()-64);Durable::Check(Durable::Hash(payload)==data.substr(data.size()-64),"checkpoint checksum");std::istringstream in(payload,std::ios::binary);std::string magic,saved;std::getline(in,magic);std::getline(in,saved);Durable::Check(magic=="H26_NATIVE_BDF2_V1"&&saved==signature,"checkpoint scientific identity");History h;in.read(reinterpret_cast<char*>(&h.step),sizeof(h.step));in.read(reinterpret_cast<char*>(&h.time),sizeof(double));in.read(reinterpret_cast<char*>(&h.dt),sizeof(double));Durable::Check(h.dt==dt&&std::abs(h.time-h.step*dt)<1e-13,"checkpoint time/history");
    for(auto*v:{&h.current,&h.previous,&h.convection,&h.previousConvection}){std::uint64_t n;in.read(reinterpret_cast<char*>(&n),sizeof(n));Durable::Check(n==std::uint64_t(size),"checkpoint state size");v->resize(n);in.read(reinterpret_cast<char*>(v->data()),n*sizeof(double));}Durable::Check(in.good(),"checkpoint read failed");return h;
}
void Advance(Flow&flow,const pt::ptree&req,std::ostream&out,const fs::path&xml,const fs::path&output,double constructionSeconds)
{
    double start=Clock();MMS m;m.geo=flow.Mapping();m.xl=flow.Elements().front().xl;m.zl=m.geo.d0;m.Lx=flow.Elements().back().xr-m.xl;m.Lz=flow.Elements().back().rr;
    m.model=req.get<std::string>("model","time");m.profile=req.get<std::string>("profile","linear");m.solution=req.get<std::string>("solution","polynomial");m.stress=req.get<double>("stress",.001);m.muRef=req.get<double>("mu_ref",.01);
    SolverUtils::AWFParameters wall;wall.cMu=572./6561.;wall.cLength=2.380402945263298;SolverUtils::AnalyticalWallFunction original({m.nu,m.geo.H,m.geo.h,m.geo.H/2},.008809015145515192,wall);m.zv=original.ViscosityInterceptHeight();
    int steps=req.get<int>("steps",20),warm=req.get<int>("warmup",0),stride=req.get<int>("checkpoint_stride",200),summaryStride=req.get<int>("summary_stride",20);double dt=req.get<double>("dt",.005),awf=req.get<double>("awf",1),tolerance=req.get<double>("tolerance",1e-11);
    std::ostringstream rule;rule<<std::setprecision(17)<<"H26-v2|"<<Durable::HashFile(xml.string())<<'|'<<dt<<'|'<<m.model<<'|'<<m.profile<<'|'<<m.solution<<'|'<<m.stress<<'|'<<m.muRef<<'|'<<awf<<'|'<<m.geo.curve<<'|'<<m.geo.d0;std::string signature=Durable::Hash(rule.str());
    double restoreSeconds=0;History h;h.dt=dt;
    auto restart=req.get<std::string>("restart","");if(!restart.empty()){double t=Clock();h=ReadHistory(restart,signature,2*flow.VelocitySize()+flow.PressureSize(),dt);restoreSeconds=Clock()-t;}
    else{h.current=Exact(flow,m.Functions(0));h.previous=h.current;h.convection=flow.Convection(h.current);h.previousConvection=h.convection;}
    Sparse knu=flow.Strain([&](double,double){return m.model=="constant"?m.muRef:m.nu;}),kt=flow.Strain([&](double x,double z){return m.model=="nu"||m.model=="constant"?0.:m.Turbulent(x,z)[0];});
    Sparse kReference=knu+kt,reference=flow.System(1.5/dt,kReference,awf),base=flow.System(1.5/dt,knu,awf),delta=reference-base,first=flow.System(1/dt,knu,awf)-base;double factorSeconds=flow.PreparePreconditioner(reference);
    // Exact time-separable analytic spatial source data: sampled continuous
    // derivatives, never A_h times a manufactured discrete state.
    MMS molecular=m;if(m.model!="constant")molecular.model="nu";
    auto quadratic=[&](double x,double z){auto u=m.Value(x,z,0);auto G=m.Gradient(x,z,0);auto p=m.PressureGradient(x,z,0);return Vector{u[0]*G[0][0]+u[1]*G[0][2]+p[0],0,u[0]*G[2][0]+u[1]*G[2][2]+p[2]};};
    Vec sourceQuadratic=flow.Load({[&](double x,double z){return quadratic(x,z)[0];},[&](double x,double z){return quadratic(x,z)[2];}});
    Vec sourceLinear=flow.Load({[&](double x,double z){return molecular.Force(x,z,0)[0]-quadratic(x,z)[0];},[&](double x,double z){return molecular.Force(x,z,0)[2]-quadratic(x,z)[2];}});
    Vec sourceTurbulent=flow.Load({[&](double x,double z){return m.Force(x,z,0)[0]-molecular.Force(x,z,0)[0];},[&](double x,double z){return m.Force(x,z,0)[2]-molecular.Force(x,z,0)[2];}});
    Vec sourceStress=flow.StressLoad([&](double x,double z){return m.Remainder(x,z,0);}),sourceAWF=flow.AWFLoad([&](double x,double z){auto g=m.PressureGradient(x,z,0);double hp=m.geo.Slope(x),tangent=(g[0]+hp*g[2])/(1+hp*hp);return Vector{-awf*tangent,0,-awf*hp*tangent};});
    Vec spatialExact=Exact(flow,m.Functions(0));double startup=constructionSeconds+Clock()-start;
    std::ofstream log(output.string()+".steps.jsonl");log<<std::setprecision(17);std::vector<fs::path>owned;std::vector<double>ordinary,checkpointTimes,means;std::vector<int>iterations;double maxMomentum=0,maxContinuity=0,minimumFeedback=1e100,minimumCorrection=1e100,minMod=1e100,maxMod=-1e100;
    std::array<double,6>cost{};int totalRebuilds=1;std::uint64_t beginStep=h.step;double beginTime=h.time;
    struct Record{int step,iterations,rebuilds;double time,wall,source,assembly,solve,setup,diagnostic,io,momentum,continuity,feedback,correction,modulation;};std::vector<Record>records;
    for(int n=0;n<steps;++n)
    {
        double began=Clock();int step=int(h.step)+1;double time=step*dt,mod=m.Modulation(time),alpha=(step==1?1.:1.5)/dt;minMod=std::min(minMod,mod);maxMod=std::max(maxMod,mod);
        Sparse K=knu+mod*kt;Vec hist,conv;if(step==1){hist=h.current/dt;conv=h.convection;}else{hist=(2*h.current-.5*h.previous)/dt;conv=2*h.convection-h.previousConvection;}
        double exponential=std::exp(-time);Vec exact=spatialExact;exact.head(2*flow.VelocitySize())*=exponential;exact.tail(flow.PressureSize())*=exponential*exponential;Vec lift=flow.Lift(exact);
        Vec load=flow.MassAction(hist)-conv+exponential*(sourceLinear+sourceStress+mod*sourceTurbulent)+exponential*exponential*(sourceQuadratic+sourceAWF);
        if(!flow.Gauge())load+=flow.TractionLoad(3,[&](double x,double z,const Vector&normal){auto G=m.Gradient(x,z,time);auto T=m.Remainder(x,z,time);auto values=m.Value(x,z,time);double mu=m.Mu(x,z,time)[0];Vector traction{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)traction[i]+=(mu*(G[i][j]+G[j][i])+T[i][j]-(i==j?values[2]:0.))*normal[j];return traction;});
        double sourceSeconds=Clock()-began,t=Clock();Sparse A=base+mod*delta;if(step==1)A+=first;Vec rhs=flow.LiftedRhs(alpha,K,awf,load,lift),x=flow.Restrict(h.current);double assemblySeconds=Clock()-t;
        auto info=flow.Solve(A,rhs,x,0.,tolerance);Vec next=flow.Expand(x,lift);t=Clock();Vec nextConvection=flow.Convection(next);double convectionSeconds=Clock()-t;
        t=Clock();Vec pressure=next.tail(flow.PressureSize()),oldPressure=h.current.tail(flow.PressureSize());Vec feedback=Vec::Zero(next.size());feedback.head(2*flow.VelocitySize())=awf*(flow.AWF()*pressure);double feedbackNorm=flow.VelocityDual(feedback)/(.01/m.Lx),correction=flow.PressureNorm(pressure-oldPressure)/.01;minimumFeedback=std::min(minimumFeedback,feedbackNorm);minimumCorrection=std::min(minimumCorrection,correction);maxMomentum=std::max(maxMomentum,info.momentum);maxContinuity=std::max(maxContinuity,info.continuity);double diagnostic=Clock()-t+info.residualSeconds;
        h.previous=h.current;h.current=next;h.previousConvection=h.convection;h.convection=nextConvection;h.step=step;h.time=time;
        t=Clock();if(step%summaryStride==0)log<<"{\"step\":"<<step<<",\"time\":"<<time<<",\"momentum\":"<<info.momentum<<",\"continuity\":"<<info.continuity<<",\"iterations\":"<<info.iterations<<",\"modulation\":"<<mod<<"}\n"<<std::flush;double io=Clock()-t;
        double wall=Clock()-began;if(n>=warm){ordinary.push_back(wall);iterations.push_back(info.iterations);cost[0]+=sourceSeconds;cost[1]+=assemblySeconds+info.solveSeconds+convectionSeconds;cost[2]+=info.factorSeconds;cost[3]+=diagnostic;cost[4]+=io;cost[5]+=wall-sourceSeconds-assemblySeconds-info.solveSeconds-convectionSeconds-info.factorSeconds-diagnostic-io;}
        records.push_back({step,info.iterations,info.rebuilds,time,wall,sourceSeconds,assemblySeconds,info.solveSeconds+convectionSeconds,info.factorSeconds,diagnostic,io,info.momentum,info.continuity,feedbackNorm,correction,mod});totalRebuilds+=info.rebuilds;
        if(stride>0&&step%stride==0){double t=Clock();fs::path p=output.string()+".checkpoint-"+std::to_string(step)+".bin";WriteHistory(p,signature,h);owned.push_back(p);while(owned.size()>2){fs::remove(owned.front());owned.erase(owned.begin());}checkpointTimes.push_back(Clock()-t);}
    }
    double endStart=Clock();fs::path final=output.string()+".final.bin";WriteHistory(final,signature,h);auto errors=flow.Error(h.current,m.Functions(h.time));if(flow.Gauge()){Vec exact=Exact(flow,m.Functions(h.time));double mean=0;Vec raw=flow.Project(m.Functions(h.time)[2],2),one=flow.Project([](double,double){return 1.;},2);mean=flow.PressureMean().dot(raw)/flow.PressureMean().dot(one);errors=flow.Error(h.current,m.Functions(h.time,mean));}
    double finalSeconds=Clock()-endStart;out<<"{\"mode\":\"advance\",\"scope\":\"prescribed scalar/stress native BDF2 flow; verification-only continuous MMS forcing\",\"signature\":\""<<signature<<"\",\"start_step\":"<<beginStep<<",\"start_time\":"<<beginTime<<",\"end_step\":"<<h.step<<",\"time\":"<<h.time<<",\"dt\":"<<dt<<",\"startup_seconds\":"<<startup<<",\"initial_preconditioner_seconds\":"<<factorSeconds<<",\"restart_seconds\":"<<restoreSeconds<<",\"final_output_and_assessment_seconds\":"<<finalSeconds<<",\"final_checkpoint\":\""<<final.string()<<"\",\"state_sha256\":\""<<StateHash(h.current)<<"\",\"previous_sha256\":\""<<StateHash(h.previous)<<"\",\"max_momentum_residual\":"<<maxMomentum<<",\"max_continuity_residual\":"<<maxContinuity<<",\"minimum_AWF_feedback\":"<<minimumFeedback<<",\"minimum_pressure_correction\":"<<minimumCorrection<<",\"coefficient_modulation_range\":["<<minMod<<','<<maxMod<<"],\"preconditioner_rebuilds\":"<<totalRebuilds<<",\"L2_errors\":["<<std::hypot(errors[0],errors[1])<<','<<errors[2]<<"],\"strong_divergence\":"<<flow.StrongDivergence(h.current)<<",\"ordinary_samples\": [";
    for(int i=0;i<ordinary.size();++i){if(i)out<<',';out<<ordinary[i];}out<<"],\"checkpoint_seconds\":[";for(int i=0;i<checkpointTimes.size();++i){if(i)out<<',';out<<checkpointTimes[i];}out<<"],\"disjoint_seconds\":[";for(int i=0;i<6;++i){if(i)out<<',';out<<cost[i];}out<<"],\"records\":[";
    for(int i=0;i<records.size();++i){if(i)out<<',';const auto&r=records[i];out<<"{\"step\":"<<r.step<<",\"time\":"<<r.time<<",\"seconds\":"<<r.wall<<",\"source\":"<<r.source<<",\"assembly\":"<<r.assembly<<",\"solve\":"<<r.solve<<",\"setup\":"<<r.setup<<",\"diagnostic\":"<<r.diagnostic<<",\"IO\":"<<r.io<<",\"iterations\":"<<r.iterations<<",\"rebuilds\":"<<r.rebuilds<<",\"momentum\":"<<r.momentum<<",\"continuity\":"<<r.continuity<<",\"AWF_feedback\":"<<r.feedback<<",\"pressure_correction\":"<<r.correction<<",\"modulation\":"<<r.modulation<<'}';}out<<"]}";
}
}
#endif
