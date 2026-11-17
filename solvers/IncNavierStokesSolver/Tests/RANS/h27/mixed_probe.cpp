// SPDX-License-Identifier: MIT
#include "../h26/common.h"
#include <IncNavierStokesSolver/EquationSystems/WJNativeMixedFlow.h>

template<class V>void ListJSON(std::ostream&out,const V&v)
{out<<'[';for(int i=0;i<int(v.size());++i){if(i)out<<',';out<<v[i];}out<<']';}
int main(int argc,char**argv)
{
    try
    {
        if(argc!=4)throw std::runtime_error("usage: WJH27MixedProbe case.xml request.json result.json");
        if(fs::exists(argv[3]))throw std::runtime_error("immutable mixed result already exists");
        pt::ptree request;pt::read_json(argv[2],request);std::string xml=argv[1],opt=xml+".opt";
        {std::ofstream f(opt);f<<"<NEKTAR><COLLECTIONS DEFAULT=\"NoCollection\"/></NEKTAR>\n";}
        std::string program="h27-mixed-foundation",option="--use-opt-file";char*args[]{program.data(),option.data(),opt.data(),xml.data(),nullptr};
        auto session=LibUtilities::SessionReader::CreateInstance(4,args);auto graph=SpatialDomains::MeshGraphIO::Read(session);
        std::vector<MultiRegions::ContFieldSharedPtr>fields;
        for(std::string name:{"u","v","p"})fields.push_back(MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,name,true,false));
        auto scalar=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,"k",true,false);
        Geometry geometry;geometry.H=.15;geometry.h=.10854503542535461;MixedFlow flow(fields,geometry);flow.SetTransportedKSpace(scalar);
        const auto&core=flow.Core();int n=flow.VelocitySize(),p=flow.PressureSize();
        std::array<Function,3>functions{[](double x,double z){return x*x*(1-x)*(2*z-3*z*z);},[](double x,double z){return-(2*x-3*x*x)*z*z*(1-z);},[](double x,double z){return .01*(x*x+z*z-2./3);}};
        Vec target(2*n+p);for(int c=0;c<2;++c)target.segment(c*n,n)=flow.Project(functions[c],c);target.tail(p)=flow.Project(functions[2],2);
        auto mu=[](double x,double z){return .02*(1+x+z);};const double alpha=2,awf=1;
        Sparse k=core.Strain(mu),a=flow.System(alpha,k,awf);Vec lift=flow.Lift(target),load=core.ApplyFull(alpha,k,awf,target),rhs=flow.LiftedRhs(alpha,k,awf,load,lift),x=Vec::Zero(flow.Size());auto info=flow.Solve(a,rhs,x);Vec solution=flow.Expand(x,lift);
        Function kf=[](double x,double z){return .008809015145515192*(1+.1*x+.05*z+.3*x*x*x*(1+.2*z*z));};
        Vec kg=flow.ProjectTransportedK(kf),known=flow.KnownScalarKLoad(kg);Vec oldP=flow.Project(kf,2),oldLoad=core.KnownKLoad(oldP).head(2*n);
        std::ofstream out(argv[3]);out<<std::setprecision(17)<<"{\"foundation_only\":true,\"full_coupled_gate_pass\":false,\"case\":\""<<request.get<std::string>("case")<<"\",\"velocity_size\":"<<n<<",\"pressure_size\":"<<p<<",\"alpha\":2,\"awf\":1,\"H\":"<<geometry.H<<",\"h\":"<<geometry.h<<",\"curve\":0,\"gauge\":"<<(flow.Gauge()?"true":"false")<<",\"constant_pressure_action_ratio\":"<<flow.ConstantPressureAction()<<",\"gauge_multiplier\":"<<(flow.Gauge()?x[flow.Size()-1]:0)<<",\"free_velocity_ids\":[";ListJSON(out,flow.FreeVelocityIds()[0]);out<<',';ListJSON(out,flow.FreeVelocityIds()[1]);out<<"],\"v_to_u\":";ListJSON(out,flow.NativeVToCanonical());out<<",\"v_to_u_sign\":";ListJSON(out,flow.NativeVSigns());
        for(const auto&item:std::vector<std::pair<std::string,const Vec*>>{{"target",&target},{"lift",&lift},{"solution",&solution},{"rhs",&rhs},{"known_scalar_coefficients",&kg},{"known_scalar_load",&known},{"oldQ2k_load",&oldLoad}}){out<<",\""<<item.first<<"\":";VectorJSON(out,*item.second);}
        out<<",\"nativeVcoefftarget\":";VectorJSON(out,flow.NativeVelocity(target,1));out<<",\"mean\":";VectorJSON(out,flow.PressureMean());out<<",\"pressure_one\":";VectorJSON(out,flow.PressureOne());out<<",\"solve\":";InfoJSON(out,info);
        for(const auto&item:std::vector<std::pair<std::string,const Sparse*>>{{"mass",&core.Mass()},{"pressure_mass",&core.PressureMass()},{"D",&core.Divergence()},{"G",&core.Gradient()},{"C",&core.AWF()},{"K",&k},{"system",&a},{"known_scalar_action",&flow.KnownScalarAction()}}){out<<",\""<<item.first<<"\":";MatrixJSON(out,*item.second);}
        out<<",\"elements\":[";auto vm=fields[1]->GetLocalToGlobalMap(),km=scalar->GetLocalToGlobalMap();
        for(int i=0;i<int(flow.Elements().size());++i){if(i)out<<',';const auto&e=flow.Elements()[i];out<<"{\"geometry\":"<<e.geometry<<",\"bounds\":["<<e.xl<<','<<e.xr<<','<<e.rl<<','<<e.rr<<"],\"u\":";ListJSON(out,e.u);out<<",\"us\":";ListJSON(out,e.us);out<<",\"p\":";ListJSON(out,e.p);out<<",\"ps\":";ListJSON(out,e.ps);
            for(auto item:std::vector<std::pair<std::string,MultiRegions::ContFieldSharedPtr>>{{"v_native",fields[1]},{"k",scalar}}){auto map=item.second->GetLocalToGlobalMap();std::vector<int>ids;std::vector<double>signs;for(int j=0;j<16;++j){int l=item.second->GetCoeff_Offset(i)+j;ids.push_back(map->GetLocalToGlobalMap(l));signs.push_back(map->GetLocalToGlobalSign(l));}out<<",\""<<item.first<<"\":";ListJSON(out,ids);out<<",\""<<(item.first=="k"?"ks":"vs_native")<<"\":";ListJSON(out,signs);}out<<'}';}
        out<<"]}\n";out.close();session->Finalise();return 0;
    }
    catch(const std::exception&e){std::cerr<<"H27 mixed failure: "<<e.what()<<'\n';return 1;}
}
