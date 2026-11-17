// SPDX-License-Identifier: MIT
#include "../h26/common.h"
#include <IncNavierStokesSolver/EquationSystems/WJNativeMixedFlow.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeScalar.h>
#include <cmath>

template<class V>void ListScalarJSON(std::ostream&o,const V&v)
{o<<'[';for(int i=0;i<int(v.size());++i){if(i)o<<',';o<<v[i];}o<<']';}
void ScalarInfoJSON(std::ostream&o,const ScalarSolveInfo&i)
{o<<"{\"iterations\":"<<i.iterations<<",\"rebuilds\":"<<i.rebuilds<<",\"residual\":"<<i.residual<<",\"factor_seconds\":"<<i.factorSeconds<<",\"solve_seconds\":"<<i.solveSeconds<<'}';}
int main(int argc,char**argv)
{
    try
    {
        if(argc!=4)throw std::runtime_error("usage: WJH27ScalarProbe case.xml full|outer result.json");
        std::string mode=argv[2];if(mode!="full"&&mode!="outer")throw std::runtime_error("invalid scalar subset");
        if(fs::exists(argv[3]))throw std::runtime_error("immutable scalar result exists");
        std::string xml=argv[1],opt=xml+".opt",program="h27-scalar-foundation",option="--use-opt-file";
        {std::ofstream f(opt);f<<"<NEKTAR><COLLECTIONS DEFAULT=\"NoCollection\"/></NEKTAR>\n";}
        char*args[]{program.data(),option.data(),opt.data(),xml.data(),nullptr};auto session=LibUtilities::SessionReader::CreateInstance(4,args);auto graph=SpatialDomains::MeshGraphIO::Read(session);
        std::vector<MultiRegions::ContFieldSharedPtr>fields;for(std::string name:{"u","v","p"})fields.push_back(MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,name,true,false));
        auto fullK=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,"k",true,false);
        Geometry geo;geo.d0=.0875;geo.H=.15;geo.h=.10854503542535461;MixedFlow mixed(fields,geo);const auto&flow=mixed.Core();
        MultiRegions::ContFieldSharedPtr field=fullK;
        if(mode=="outer")
        {
            std::vector<unsigned int>ids;for(int i=0;i<fullK->GetExpSize();++i)if(flow.Elements()[i].rl>=geo.H-1e-12)ids.push_back(i);
            auto subset=MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(*fullK,ids,true);field=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,*subset);
        }
        const double kin=.008809015145515192,cl=2.3804029452632984,A=std::pow(kin,1.5)/cl;
        Function scale=mode=="full"?Function([=](double,double){return kin;}):Function([=](double,double z){return A/(z-geo.d0);});
        Scalar scalar(field,flow,scale,mode=="full"?std::vector<int>{0}:std::vector<int>{0,2});RuleFactory factory(flow,8,12);
        auto cutoff=[](const Element&,double x){return .03*(1+.03*x);};Rules rules=factory.Cutoff(cutoff);
        Function qf=[=](double x,double z){double r=z-geo.d0;return .01*(1+.2*x+.1*r+.3*x*x*x*(1+.1*r*r));};
        Vec target=scalar.Project(qf),lift=scalar.Lift(target),one=scalar.Project([](double,double){return 1.;});
        Vec velocity=Vec::Zero(2*flow.VelocitySize()+flow.PressureSize());velocity.head(flow.VelocitySize())=mixed.Project([=](double x,double z){return .2+.1*x+.03*(z-geo.d0);},0);velocity.segment(flow.VelocitySize(),flow.VelocitySize())=mixed.Project([=](double,double z){return .04+.02*(z-geo.d0);},1);
        std::vector<Vec>d,react,source;
        for(const auto&pointer:rules){const auto&r=*pointer;Vec a(r.w.size()),b(a.size()),f(a.size());for(int i=0;i<a.size();++i){double z=r.z[i]-geo.d0,c=cutoff(flow.Elements()[0],r.x[i]);a[i]=.003*(1+r.x[i]+2*z)*std::max(z,c);b[i]=.03*(1+.1*r.x[i]);f[i]=1+r.x[i]*r.x[i]+z*z*z+(1+.02*r.x[i])/std::max(z,c);}d.push_back(a);react.push_back(b);source.push_back(f);}
        Sparse k=scalar.Diffusion(rules,d,react),system=scalar.System(2.,k);Vec load=2*(scalar.Mass()*target)+k*target,rhs=scalar.LiftedRhs(2,k,load,lift),free=Vec::Zero(scalar.FreeSize());auto info=scalar.Solve(system,rhs,free);Vec solution=scalar.Expand(free,lift);
        Sparse k2=1.2*k,system2=scalar.System(2,k2);Vec rhs2=scalar.LiftedRhs(2,k2,load,lift);auto info2=scalar.Solve(system2,rhs2,free,.2);Vec solution2=scalar.Expand(free,lift);
        Vec adv=scalar.Advection(rules,velocity,target),src=scalar.Source(rules,source),flux=scalar.NaturalFlux([=](double x,double z,const Vector&){return .002*(1+x+.5*(z-geo.d0));},factory);
        Sparse transport=scalar.TransportMatrix(rules,velocity);
        Vec diffusionAction=scalar.DiffusionAction(rules,d,target,react);
        if((diffusionAction-k*target).cwiseAbs().maxCoeff()>1e-12)throw std::runtime_error("direct scalar action differs on full signed rows");
        bool incompatible=false;try{scalar.EssentialLift([=](double x,double z,const Vector&n){return qf(x,z)+(mode=="outer"&&n[0]<-.5?1e-5:0);},factory);}catch(const std::exception&){incompatible=true;}
        Vec gamma=Vec::Zero(scalar.Size());if(mode=="outer")gamma=scalar.EssentialLift([=](double x,double z,const Vector&n){return n[2]<-.5?A*std::pow(1+.02*x,3)/geo.H:A/(z-geo.d0);},factory);
        std::ofstream o(argv[3]);o<<std::setprecision(17)<<"{\"foundation_only\":true,\"full_coupled_gate_pass\":false,\"case\":\""<<mode<<"\",\"d0\":"<<geo.d0<<",\"H\":"<<geo.H<<",\"h\":"<<geo.h<<",\"alpha\":2,\"Kin\":"<<kin<<",\"Ain\":"<<A<<",\"velocity_size\":"<<flow.VelocitySize()<<",\"scalar_size\":"<<scalar.Size()<<",\"native_dirichlet_count\":"<<field->GetLocalToGlobalMap()->GetNumGlobalDirBndCoeffs()<<",\"essential_count\":"<<scalar.DirichletSize()<<",\"corner_conflict_rejected\":"<<(incompatible?"true":"false")<<",\"free_ids\":";ListScalarJSON(o,scalar.FreeIds());
        for(const auto&v:std::vector<std::pair<std::string,const Vec*>>{{"target",&target},{"lift",&lift},{"one",&one},{"solution",&solution},{"rhs",&rhs},{"changed_solution",&solution2},{"changed_rhs",&rhs2},{"velocity",&velocity},{"advection",&adv},{"source_load",&src},{"natural_flux_rhs",&flux},{"dependent_gamma_lift",&gamma},{"diffusion_action",&diffusionAction}}){o<<",\""<<v.first<<"\":";VectorJSON(o,*v.second);}
        for(const auto&m:std::vector<std::pair<std::string,const Sparse*>>{{"mass",&scalar.Mass()},{"K",&k},{"system",&system},{"transport",&transport}}){o<<",\""<<m.first<<"\":";MatrixJSON(o,*m.second);}
        o<<",\"solve\":";ScalarInfoJSON(o,info);o<<",\"changed_solve\":";ScalarInfoJSON(o,info2);o<<",\"flow_elements\":[";
        for(int i=0;i<int(flow.Elements().size());++i){if(i)o<<',';const auto&e=flow.Elements()[i];o<<"{\"geometry\":"<<e.geometry<<",\"bounds\":["<<e.xl<<','<<e.xr<<','<<e.rl<<','<<e.rr<<"],\"u\":";ListScalarJSON(o,e.u);o<<",\"us\":";ListScalarJSON(o,e.us);o<<'}';}
        o<<"],\"scalar_elements\":[";auto map=field->GetLocalToGlobalMap();
        for(int i=0;i<field->GetExpSize();++i){if(i)o<<',';int geometry=field->GetExp(i)->GetGeom()->GetGlobalID();const Element*e=nullptr;for(const auto&v:flow.Elements())if(v.geometry==geometry)e=&v;if(!e)throw std::runtime_error("scalar dump owner missing");std::vector<int>ids;std::vector<double>signs;for(int j=0;j<16;++j){int l=field->GetCoeff_Offset(i)+j;ids.push_back(map->GetLocalToGlobalMap(l));signs.push_back(map->GetLocalToGlobalSign(l));}o<<"{\"geometry\":"<<geometry<<",\"bounds\":["<<e->xl<<','<<e->xr<<','<<e->rl<<','<<e->rr<<"],\"s\":";ListScalarJSON(o,ids);o<<",\"ss\":";ListScalarJSON(o,signs);o<<'}';}
        o<<"]}\n";o.close();session->Finalise();return 0;
    }
    catch(const std::exception&e){std::cerr<<"H27 scalar failure: "<<e.what()<<'\n';return 1;}
}
