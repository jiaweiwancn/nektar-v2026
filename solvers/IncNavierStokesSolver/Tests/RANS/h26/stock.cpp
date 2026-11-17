// SPDX-License-Identifier: MIT
// Actual stock VCSWeakPressure and native IMEX(2); instrumentation only.
#include "common.h"
#include "advance.h"
#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionSchemeWeakPressure.h>
#include <SolverUtils/Forcing/Forcing.h>
class PrescribedForce:public SolverUtils::Forcing
{
public:
    PrescribedForce(const LibUtilities::SessionReaderSharedPtr&s,const std::weak_ptr<SolverUtils::EquationSystem>&eq,const std::vector<MultiRegions::ContFieldSharedPtr>&fields,H26::MMS m,double dt):Forcing(s,eq),m_dt(dt)
    {
        int n=fields[0]->GetTotPoints();Array<OneD,NekDouble>x(n),z(n),dummy(n);fields[0]->GetCoords(x,z,dummy);m_NumVariable=2;m_Forcing=Array<OneD,Array<OneD,NekDouble>>(2);for(int c=0;c<2;++c){m_Forcing[c]=Array<OneD,NekDouble>(n);linear[c]=Vec(n);quadratic[c]=Vec(n);}
        for(int q=0;q<n;++q){auto u=m.Value(x[q],z[q],0);auto G=m.Gradient(x[q],z[q],0);auto p=m.PressureGradient(x[q],z[q],0);auto f=m.Force(x[q],z[q],0);for(int c=0;c<2;++c){int i=c?2:0;quadratic[c][q]=u[0]*G[i][0]+u[1]*G[i][2]+p[i];linear[c][q]=f[i]-quadratic[c][q];}}
    }
    double seconds=0;int calls=0;
private:
    // De-extrapolation makes the unmodified stock EX2 forcing combination
    // equal the prescribed new-time body force, including BE startup.
    double Weight(double rate,double t)const{double rho=std::exp(-rate*m_dt),n=std::round(t/m_dt),a=rho*rho/(2*rho-1);return a*std::pow(rho,n)+(rho-a)*std::pow(.5,n);}
    void v_InitObject(const Array<OneD,MultiRegions::ExpListSharedPtr>&,const unsigned int&,const TiXmlElement*)override{}
    void v_Apply(const Array<OneD,MultiRegions::ExpListSharedPtr>&,const Array<OneD,Array<OneD,NekDouble>>&,Array<OneD,Array<OneD,NekDouble>>&out,const NekDouble&t)override
    {double start=Clock(),a=Weight(1,t),b=Weight(2,t);for(int c=0;c<2;++c)for(int q=0;q<linear[c].size();++q){m_Forcing[c][q]=a*linear[c][q]+b*quadratic[c][q];out[c][q]+=m_Forcing[c][q];}seconds+=Clock()-start;++calls;}
    double m_dt;std::array<Vec,2>linear,quadratic;
};
class Stock:public VCSWeakPressure
{
public:
    Stock(const LibUtilities::SessionReaderSharedPtr&s,const SpatialDomains::MeshGraphSharedPtr&g):UnsteadySystem(s,g),VCSWeakPressure(s,g){}
    void Before(int n){v_PreIntegrate(n);}void After(int n){v_PostIntegrate(n);}
    void RegisterForce(const SolverUtils::ForcingSharedPtr&force){AddForcing(force);m_extrapolation->SetForcing(m_forcing);}
};
Vec ReadState(const std::vector<MultiRegions::ContFieldSharedPtr>&fields,Flow&flow)
{
    Vec result(2*flow.VelocitySize()+flow.PressureSize());for(int c=0;c<3;++c){auto map=fields[c]->GetLocalToGlobalMap();Array<OneD,NekDouble>global(map->GetNumGlobalCoeffs());map->LocalToGlobal(fields[c]->GetCoeffs(),global,true);for(int i=0;i<global.size();++i)result[(c<2?c*flow.VelocitySize():2*flow.VelocitySize())+i]=global[i];}
    if(flow.Gauge()){const Vec&one=flow.PressureOne();result.tail(flow.PressureSize())-=one*(flow.PressureMean().dot(result.tail(flow.PressureSize()))/flow.PressureMean().dot(one));}return result;
}
int main(int argc,char**argv)
{
    try
    {
        double began=Clock();if(argc!=4||fs::exists(argv[3]))throw std::runtime_error("stock usage or existing result");pt::ptree req;pt::read_json(argv[2],req);std::string xml=argv[1],opt=xml+".opt",program="h26-stock",option="--use-opt-file";{std::ofstream f(opt);f<<"<NEKTAR><COLLECTIONS DEFAULT=\"NoCollection\"/></NEKTAR>\n";}char*args[]{program.data(),option.data(),opt.data(),xml.data(),nullptr};auto session=LibUtilities::SessionReader::CreateInstance(4,args);auto graph=SpatialDomains::MeshGraphIO::Read(session);auto stock=std::make_shared<Stock>(session,graph);stock->SolverUtils::EquationSystem::InitObject();stock->DoInitialise(false);
        std::vector<MultiRegions::ContFieldSharedPtr>fields;for(const auto&field:stock->UpdateFields())fields.push_back(std::dynamic_pointer_cast<MultiRegions::ContField>(field));Geometry geometry;geometry.d0=req.get<double>("d0",.0875);Flow flow(fields,geometry);H26::MMS m;m.geo=geometry;m.xl=flow.Elements().front().xl;m.Lx=flow.Elements().back().xr-m.xl;m.zl=geometry.d0;m.Lz=flow.Elements().back().rr;m.model=req.get<std::string>("model","nu");m.muRef=req.get<double>("mu_ref",.01);m.stress=0;double dt=req.get<double>("dt",.005);stock->SetTimeStep(dt);
        auto forcing=std::make_shared<PrescribedForce>(session,stock,fields,m,dt);stock->RegisterForce(forcing);Vec initial=Exact(flow,m.Functions(0));
        for(int c=0;c<3;++c){Array<OneD,NekDouble>global(c<2?flow.VelocitySize():flow.PressureSize());for(int i=0;i<global.size();++i)global[i]=initial[(c<2?c*flow.VelocitySize():2*flow.VelocitySize())+i];fields[c]->GetLocalToGlobalMap()->GlobalToLocal(global,fields[c]->UpdateCoeffs());fields[c]->BwdTrans(fields[c]->GetCoeffs(),fields[c]->UpdatePhys());fields[c]->SetPhysState(false);}
        Array<OneD,Array<OneD,NekDouble>>values(2);for(int c=0;c<2;++c)values[c]=fields[c]->UpdatePhys();auto scheme=stock->GetTimeIntegrationScheme();scheme->InitializeScheme(dt,values,0,stock->GetTimeIntegrationSchemeOperators());
        Sparse K=flow.Strain([&](double,double){return m.model=="nu"?m.nu:m.muRef;}),A=flow.System(1.5/dt,K,0),Afirst=flow.System(1/dt,K,0);H26::History history;history.dt=dt;history.current=history.previous=initial;history.convection=history.previousConvection=flow.Convection(initial);Vec force0=flow.Load({[&](double x,double z){return m.Force(x,z,0)[0];},[&](double x,double z){return m.Force(x,z,0)[2];}});(void)force0;
        Vec forceLinear=flow.Load({[&](double x,double z){auto u=m.Value(x,z,0);auto G=m.Gradient(x,z,0);return m.Force(x,z,0)[0]-u[0]*G[0][0]-u[1]*G[0][2]-m.PressureGradient(x,z,0)[0];},[&](double x,double z){auto u=m.Value(x,z,0);auto G=m.Gradient(x,z,0);return m.Force(x,z,0)[2]-u[0]*G[2][0]-u[1]*G[2][2]-m.PressureGradient(x,z,0)[2];}}),forceQuadratic=force0-forceLinear;
        double startup=Clock()-began;int steps=req.get<int>("steps",420),warm=req.get<int>("warmup",20),stride=req.get<int>("checkpoint_stride",200);std::vector<double>samples,writes;std::array<double,6>cost{};std::ofstream log(std::string(argv[3])+".steps.jsonl");std::vector<fs::path>owned;double maxMomentum=0,maxContinuity=0;std::vector<std::array<double,7>>records;std::string signature=H26::Durable::Hash(std::string("H26-stock|")+H26::Durable::HashFile(xml)+std::to_string(dt));
        if(const char*operation=std::getenv("H26_OPERATION")){std::ifstream maps("/proc/self/maps");std::ofstream snapshot(std::string(operation)+".native-maps");snapshot<<maps.rdbuf();}
        for(int n=0;n<steps;++n)
        {
            double start=Clock(),previousSource=forcing->seconds;stock->SetTime(n*dt);stock->Before(n);values=scheme->TimeIntegrate(n,dt);stock->SetTime((n+1)*dt);for(int c=0;c<2;++c){fields[c]->SetPhys(values[c]);values[c]=fields[c]->UpdatePhys();fields[c]->FwdTransLocalElmt(values[c],fields[c]->UpdateCoeffs());fields[c]->SetPhysState(false);}stock->After(n);double equation=Clock()-start,source=forcing->seconds-previousSource,t=Clock();
            // Compare the actual split scheme to the same full-strain target;
            // these physical residuals are reported, never labelled algebraic
            // stock solve errors or forced to vanish by a hidden body load.
            Vec next=ReadState(fields,flow),hist,conv;if(n==0){hist=history.current/dt;conv=history.convection;}else{hist=(2*history.current-.5*history.previous)/dt;conv=2*history.convection-history.previousConvection;}double time=(n+1)*dt,a=std::exp(-time);Vec lift=flow.Lift(next),load=flow.MassAction(hist)-conv+a*forceLinear+a*a*forceQuadratic,rhs=flow.LiftedRhs((n?1.5:1)/dt,K,0,load,lift);auto residual=flow.Residual(n?A:Afirst,rhs,flow.Restrict(next));maxMomentum=std::max(maxMomentum,residual.momentum);maxContinuity=std::max(maxContinuity,residual.continuity);history.previous=history.current;history.current=next;history.previousConvection=history.convection;history.convection=flow.Convection(next);history.step=n+1;history.time=time;double diagnostic=Clock()-t;t=Clock();if((n+1)%20==0)log<<"{\"step\":"<<n+1<<",\"time\":"<<time<<",\"momentum_target_difference\":"<<residual.momentum<<",\"weak_continuity\":"<<residual.continuity<<"}\n"<<std::flush;double io=Clock()-t,wall=Clock()-start;
            if(n>=warm){samples.push_back(wall);cost[0]+=source;cost[1]+=equation-source;cost[3]+=diagnostic;cost[4]+=io;cost[5]+=wall-equation-diagnostic-io;}records.push_back({wall,source,equation-source,diagnostic,io,residual.momentum,residual.continuity});
            if(stride>0&&(n+1)%stride==0){t=Clock();fs::path checkpoint=std::string(argv[3])+".comparison-state-"+std::to_string(n+1)+".bin";H26::WriteHistory(checkpoint,signature,history);owned.push_back(checkpoint);while(owned.size()>2){fs::remove(owned.front());owned.erase(owned.begin());}writes.push_back(Clock()-t);}
        }
        double t=Clock();H26::WriteHistory(std::string(argv[3])+".final.bin",signature,history);Vec raw=flow.Project(m.Functions(history.time)[2],2),one=flow.Project([](double,double){return 1.;},2);double mean=flow.PressureMean().dot(raw)/flow.PressureMean().dot(one);auto errors=flow.Error(history.current,m.Functions(history.time,mean));double final=Clock()-t;std::ofstream out(argv[3]);out<<std::setprecision(17)<<"{\"mode\":\"actual_stock_VCSWeakPressure_IMEX2\",\"scope\":\"stock no-RANS, de-extrapolated prescribed body force; comparison checkpoint excludes stock pressure-BC/GLM histories and is not a qualified stock restart\",\"startup_seconds\":"<<startup<<",\"final_output_and_assessment_seconds\":"<<final<<",\"time\":"<<history.time<<",\"dt\":"<<dt<<",\"L2_errors\":["<<std::hypot(errors[0],errors[1])<<','<<errors[2]<<"],\"max_full_strain_target_difference\":"<<maxMomentum<<",\"max_weak_continuity\":"<<maxContinuity<<",\"forcing_calls\":"<<forcing->calls<<",\"ordinary_samples\":[";for(int i=0;i<samples.size();++i){if(i)out<<',';out<<samples[i];}out<<"],\"checkpoint_seconds\":[";for(int i=0;i<writes.size();++i){if(i)out<<',';out<<writes[i];}out<<"],\"disjoint_seconds\":[";for(int i=0;i<6;++i){if(i)out<<',';out<<cost[i];}out<<"],\"records\":[";for(int i=0;i<records.size();++i){if(i)out<<',';out<<'[';for(int j=0;j<7;++j){if(j)out<<',';out<<records[i][j];}out<<']';}out<<"]}\n";session->Finalise();return 0;
    }catch(const std::exception&e){std::cerr<<"H26 stock failure: "<<e.what()<<'\n';return 1;}
}
