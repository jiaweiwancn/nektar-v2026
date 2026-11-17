#include "common.h"
#include "advance.h"
#include "checks.h"
void Algebra(Flow&flow,const pt::ptree&request,std::ostream&o)
{
    double alpha=request.get<double>("alpha",300),awf=request.get<double>("awf",1.);
    auto mu=[](double x,double z){return .02*(1+x+z);};Sparse k=flow.Strain(mu),a=flow.System(alpha,k,awf);
    Vec target=Exact(flow,Polynomial(0.)),lift=flow.Lift(target),load=flow.ApplyFull(alpha,k,awf,target),b=flow.LiftedRhs(alpha,k,awf,load,lift),x=Vec::Zero(flow.Size());auto info=flow.Solve(a,b,x);
    Vec solved=flow.Expand(x,lift);double error=(solved-target).norm()/target.norm();Sparse k2=1.2*k,a2=flow.System(alpha,k2,awf);Vec b2=flow.LiftedRhs(alpha,k2,awf,load,lift);auto changed=flow.Solve(a2,b2,x,0.);Vec solved2=flow.Expand(x,lift);
    o<<"{\"mode\":\"algebra_only_not_continuous_MMS\",\"velocity_size\":"<<flow.VelocitySize()<<",\"pressure_size\":"<<flow.PressureSize()<<",\"dirichlet_size\":"<<flow.DirichletSize()<<",\"gauge\":"<<(flow.Gauge()?"true":"false")<<",\"alpha\":"<<alpha<<",\"awf\":"<<awf<<",\"error\":"<<error<<",\"solve\":";InfoJSON(o,info);o<<",\"changed_mu_retained_preconditioner\":";InfoJSON(o,changed);
    o<<",\"target\":";VectorJSON(o,target);o<<",\"lift\":";VectorJSON(o,lift);o<<",\"solution\":";VectorJSON(o,solved);o<<",\"rhs\":";VectorJSON(o,b);
    o<<",\"changed_rhs\":";VectorJSON(o,b2);o<<",\"changed_solution\":";VectorJSON(o,solved2);
    Vec knownK=flow.Project([](double x,double z){return .008809015145515192*(1+.1*x+.05*z);},2);o<<",\"known_k_coefficients\":";VectorJSON(o,knownK);o<<",\"known_mechanical_pressure_k_load\":";VectorJSON(o,flow.KnownKLoad(knownK));
    Vec traction=flow.TractionLoad(3,[&](double x,double z,const Vector&n){double ux=(2*x-3*x*x)*(2*z-3*z*z),uz=x*x*(1-x)*(2-6*z),vx=-(2-6*x)*z*z*(1-z),p=x*x+z*z;double tx=mu(x,z)*2*ux+.001*x*z-p,tz=-mu(x,z)*2*ux-.001*x*z-p,off=mu(x,z)*(uz+vx)+.0005*(x+z);return Vector{tx*n[0]+off*n[2],0,off*n[0]+tz*n[2]};});o<<",\"top_total_traction\":";VectorJSON(o,traction);
    for(const auto&entry:std::vector<std::pair<std::string,const Sparse*>>{{"mass",&flow.Mass()},{"pressure_mass",&flow.PressureMass()},{"D",&flow.Divergence()},{"G",&flow.Gradient()},{"C",&flow.AWF()},{"K",&k},{"system",&a}}){o<<",\""<<entry.first<<"\":";MatrixJSON(o,*entry.second);}
    o<<",\"mean\":";VectorJSON(o,flow.PressureMean());o<<",\"strong_divergence\":"<<flow.StrongDivergence(solved)<<",\"elements\":[";bool comma=false;
    for(const auto&e:flow.Elements()){if(comma)o<<',';comma=true;o<<"{\"geometry\":"<<e.geometry<<",\"bounds\":["<<e.xl<<','<<e.xr<<','<<e.rl<<','<<e.rr<<"],\"u\":[";for(int j=0;j<16;++j){if(j)o<<',';o<<e.u[j];}o<<"],\"p\":[";for(int j=0;j<9;++j){if(j)o<<',';o<<e.p[j];}o<<"],\"us\":[";for(int j=0;j<16;++j){if(j)o<<',';o<<e.us[j];}o<<"],\"ps\":[";for(int j=0;j<9;++j){if(j)o<<',';o<<e.ps[j];}o<<"],\"x\":";VectorJSON(o,e.x);o<<",\"z\":";VectorJSON(o,e.z);o<<",\"weights\":";VectorJSON(o,e.w);o<<'}';}
    o<<"]}";
}
int main(int argc,char**argv)
{
    try
    {
        if(argc!=4)throw std::runtime_error("usage: WJNativeFlowDriver case.xml request.json result.json");
        double began=Clock();if(fs::exists(argv[3]))throw std::runtime_error("immutable H26 result already exists");pt::ptree request;pt::read_json(argv[2],request);std::string xml=argv[1],opt=xml+".opt";{std::ofstream f(opt);f<<"<NEKTAR><COLLECTIONS DEFAULT=\"NoCollection\"/></NEKTAR>\n";}
        std::string program="h26-native-flow",option="--use-opt-file";char*args[]{program.data(),option.data(),opt.data(),xml.data(),nullptr};auto session=LibUtilities::SessionReader::CreateInstance(4,args);auto graph=SpatialDomains::MeshGraphIO::Read(session);
        std::vector<MultiRegions::ContFieldSharedPtr>fields;for(std::string name:{"u","v","p"})fields.push_back(MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,name,true,false));
        Geometry geometry;geometry.d0=request.get<double>("d0",0);geometry.curve=request.get<double>("curve",0);geometry.H=request.get<double>("H",.15);geometry.h=request.get<double>("h",.10854503542535461);
        Flow flow(fields,geometry);double construction=Clock()-began;if(const char*operation=std::getenv("H26_OPERATION")){std::ifstream maps("/proc/self/maps");std::ofstream snapshot(std::string(operation)+".native-maps");snapshot<<maps.rdbuf();}
        std::ofstream out(argv[3]);out<<std::setprecision(17);auto mode=request.get<std::string>("mode","algebra");if(mode=="algebra")Algebra(flow,request,out);else if(mode=="advance")H26::Advance(flow,request,out,xml,argv[3],construction);else if(mode=="adapter")H26::AdapterChecks(out);else if(mode=="nonlinear")H26::Nonlinear(flow,out);else if(mode=="formula")H26::Formula(out);else if(mode=="wall_trace")H26::WallTrace(flow,out);else throw std::runtime_error("unknown H26 mode");out<<'\n';out.close();session->Finalise();return 0;
    }
    catch(const std::exception&e){std::cerr<<"H26 failure: "<<e.what()<<'\n';return 1;}
}
