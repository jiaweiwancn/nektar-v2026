// SPDX-License-Identifier: MIT
#include "../h26/common.h"
#include <IncNavierStokesSolver/EquationSystems/WJNativePerturbation.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeMixedFlow.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeScalar.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <fcntl.h>
#include <iterator>
#include <map>
#include <set>
#include <sstream>
#include <unistd.h>

namespace
{
void Require(bool ok,const std::string&why){if(!ok)throw std::runtime_error(why);}
std::string Bytes(const std::string&path)
{std::ifstream f(path,std::ios::binary);Require(bool(f),"cannot read "+path);std::string b((std::istreambuf_iterator<char>(f)),{});Require(!f.bad(),"read failed "+path);return b;}
void Create(const std::string&path,const std::string&bytes)
{
    int fd=::open(path.c_str(),O_WRONLY|O_CREAT|O_EXCL,0600);Require(fd>=0,"create-only output failed");
    std::size_t n=0;while(n<bytes.size()){auto count=::write(fd,bytes.data()+n,bytes.size()-n);if(count<=0){::close(fd);throw std::runtime_error("output write failed");}n+=count;}
    const int sync=::fsync(fd),close=::close(fd);Require(sync==0&&close==0,"output flush/close failed");
}
template<class T>void List(std::ostream&o,const T&v){o<<'[';bool first=true;for(auto x:v){if(!first)o<<',';first=false;o<<x;}o<<']';}
std::vector<double> Numbers(const pt::ptree&t,const std::string&key)
{std::vector<double>v;for(auto&x:t.get_child(key))v.push_back(x.second.get_value<double>());return v;}
bool Close(double a,double b){return std::abs(a-b)<1e-12;}
std::vector<double> Edges(std::vector<double>v)
{std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end(),Close),v.end());return v;}
int Index(const std::vector<double>&v,double x)
{int answer=-1;for(int i=0;i<int(v.size());++i)if(Close(v[i],x)){Require(answer<0,"ambiguous tensor edge");answer=i;}Require(answer>=0,"missing tensor edge");return answer;}
std::string ShapeName(PerturbationShape s)
{return s==PerturbationShape::PrimaryQuartile?"PrimaryQuartile":"SecondaryLowerSkew";}
PerturbationShape Shape(const std::string&s)
{if(s=="PrimaryQuartile")return PerturbationShape::PrimaryQuartile;if(s=="SecondaryLowerSkew")return PerturbationShape::SecondaryLowerSkew;throw std::runtime_error("unknown perturbation shape");}
void LayoutJSON(std::ostream&o,const PerturbationLayout&l)
{o<<"{\"shape\":\""<<ShapeName(l.shape)<<"\",\"x_knots\":";List(o,l.xKnots);o<<",\"r_knots\":";List(o,l.rKnots);o<<",\"H\":"<<l.H<<",\"epsilon_peak_radius\":"<<l.epsilonPeakRadius<<'}';}
PerturbationLayout ReadLayout(const pt::ptree&t)
{PerturbationLayout l;l.shape=Shape(t.get<std::string>("shape"));auto x=Numbers(t,"x_knots"),r=Numbers(t,"r_knots");Require(x.size()==5&&r.size()==5,"wrong knot count");std::copy(x.begin(),x.end(),l.xKnots.begin());std::copy(r.begin(),r.end(),l.rKnots.begin());l.H=t.get<double>("H");l.epsilonPeakRadius=t.get<double>("epsilon_peak_radius");return l;}
struct Map{std::vector<int>ids;std::vector<double>signs;};
std::map<int,Map> Maps(const MultiRegions::ContFieldSharedPtr&f)
{std::map<int,Map>out;auto m=f->GetLocalToGlobalMap();for(int i=0;i<f->GetExpSize();++i){Map cell;Require(f->GetExp(i)->GetNcoeffs()==16,"probe requires scalar Q3");for(int j=0;j<16;++j){const int l=f->GetCoeff_Offset(i)+j;cell.ids.push_back(m->GetLocalToGlobalMap(l));cell.signs.push_back(m->GetLocalToGlobalSign(l));}Require(out.emplace(f->GetExp(i)->GetGeom()->GetGlobalID(),cell).second,"duplicate scalar geometry");}return out;}
void ElementsJSON(std::ostream&o,const Flow&flow,const std::map<int,Map>&k,const std::map<int,Map>&epsilon)
{
    o<<'[';bool first=true;for(const auto&e:flow.Elements()){
        if(!first)o<<',';first=false;o<<"{\"geometry\":"<<e.geometry<<",\"bounds\":["<<e.xl<<','<<e.xr<<','<<e.rl<<','<<e.rr<<"],\"u\":";List(o,e.u);o<<",\"us\":";List(o,e.us);o<<",\"p\":";List(o,e.p);o<<",\"ps\":";List(o,e.ps);
        auto emit=[&](const std::string&name,const std::string&sign,const std::map<int,Map>&map){o<<",\""<<name<<"\":";auto found=map.find(e.geometry);if(found==map.end())o<<"[]";else List(o,found->second.ids);o<<",\""<<sign<<"\":";if(found==map.end())o<<"[]";else List(o,found->second.signs);};
        emit("k","ks",k);emit("epsilon","epsilon_signs",epsilon);o<<'}';
    }o<<']';
}
}

int main(int argc,char**argv)
{
    try{
        Require(argc==5,"usage: WJH27PerturbationProbe case.xml reference.json request.json new-output.json");
        Require(!fs::exists(argv[4]),"immutable output exists");pt::ptree request;pt::read_json(argv[3],request);
        const std::string mode=request.get<std::string>("mode"),source=request.get<std::string>("source_identity");
        Require(mode=="layout"||mode=="project","invalid perturbation probe mode");
        Require(source.size()==64&&source.find_first_not_of("0123456789abcdef")==std::string::npos,"source identity malformed");
        Reference reference(argv[2],request.get<std::string>("reference_sha256"));const auto&p=reference.Parameters();
        Geometry geometry;geometry.d0=p.d0;geometry.H=p.H;geometry.h=p.roughness;
        std::string xml=argv[1],opt=xml+".opt",program="h27-perturbation-foundation",option="--use-opt-file";
        Require(fs::exists(opt),"read-only native fixture needs its existing .opt file");
        char*args[]{program.data(),option.data(),opt.data(),xml.data(),nullptr};
        auto session=LibUtilities::SessionReader::CreateInstance(4,args);auto graph=SpatialDomains::MeshGraphIO::Read(session);
        std::vector<MultiRegions::ContFieldSharedPtr>fields;for(std::string name:{"u","v","p"})fields.push_back(MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,name,true,false));
        auto kField=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,"k",true,false);
        MixedFlow mixed(fields,geometry);const auto&flow=mixed.Core();std::vector<double>xs,rs;
        for(const auto&e:flow.Elements()){xs.push_back(e.xl);xs.push_back(e.xr);rs.push_back(e.rl);rs.push_back(e.rr);}
        xs=Edges(xs);rs=Edges(rs);std::set<std::pair<int,int>>seen;
        for(const auto&e:flow.Elements()){
            int a=Index(xs,e.xl),b=Index(xs,e.xr),c=Index(rs,e.rl),d=Index(rs,e.rr);
            Require(b==a+1&&d==c+1&&seen.emplace(a,c).second,"not a common native tensor mesh");
        }
        Require(seen.size()==(xs.size()-1)*(rs.size()-1),"incomplete native tensor inventory");
        if(mode=="layout")Require(flow.Elements().size()==700,"initial layout must use retained 700-cell mesh");
        std::vector<unsigned int>outer;for(int i=0;i<kField->GetExpSize();++i){int id=kField->GetExp(i)->GetGeom()->GetGlobalID();auto found=std::find_if(flow.Elements().begin(),flow.Elements().end(),[&](const Element&e){return e.geometry==id;});Require(found!=flow.Elements().end(),"missing scalar geometry");if(found->rl>=p.H-1e-12)outer.push_back(i);}
        auto subset=MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(*kField,outer,true);
        auto eField=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,*subset);
        auto km=Maps(kField),em=Maps(eField);
        Scalar k(kField,flow,[&](double,double){return p.KIn;},{0});
        const double Ain=reference.Amplitude(p.KIn);
        Scalar epsilon(eField,flow,[&](double,double z){return Ain/(z-p.d0);},{0,2});
        RuleFactory factory(flow,12,16);
        std::ostringstream out;out<<std::setprecision(17)<<"{\"schema\":\"WJH27PerturbationNativeProbeV1\",\"status\":\"RECORDED\",\"mode\":\""<<mode<<"\",\"full_coupled_gate_pass\":false,\"developed_reference_qualified\":false,\"tensor_inventory_verified\":true,\"source_identity\":\""<<source<<"\",\"reference_identity\":\""<<reference.Identity()<<"\",\"mesh_xml_sha256\":\""<<RANSWallR3::Hash(Bytes(xml))<<"\",\"request_sha256\":\""<<RANSWallR3::Hash(Bytes(argv[3]))<<"\",\"helper_identity\":\""<<Perturbation::Identity<<"\",\"d0\":"<<p.d0<<",\"H\":"<<p.H<<",\"Kin\":"<<p.KIn<<",\"Ain\":"<<Ain<<",\"c_length\":"<<p.cLength<<",\"domain\":["<<xs.front()<<','<<xs.back()<<','<<rs.front()<<','<<rs.back()<<"],\"sizes\":["<<mixed.VelocitySize()<<','<<mixed.PressureSize()<<','<<k.Size()<<','<<epsilon.Size()<<"],\"x_edges\":";List(out,xs);out<<",\"r_edges\":";List(out,rs);
        if(mode=="layout"){
            out<<",\"layouts\":{";bool first=true;for(auto shape:{PerturbationShape::PrimaryQuartile,PerturbationShape::SecondaryLowerSkew}){
                auto layout=Perturbation::SelectLayout(xs,rs,p.H,shape);Perturbation seed(reference,layout,xs,rs);
                if(!first)out<<',';first=false;out<<'"'<<ShapeName(shape)<<"\":";LayoutJSON(out,layout);
            }out<<'}';
        }else{
            const std::string evidence=Bytes(request.get<std::string>("certificate"));const std::string identity=RANSWallR3::Hash(evidence);
            Require(identity==request.get<std::string>("certificate_sha256"),"immutable certificate SHA differs");
            pt::ptree certificate;std::istringstream stream(evidence);pt::read_json(stream,certificate);
            Require(certificate.get<std::string>("schema")=="WJH27PerturbationPeakCertificateV1"&&certificate.get<std::string>("status")=="PASS_MATHEMATICAL_SHAPE_ONLY","unqualified mathematical certificate");
            Require(certificate.get<std::string>("reference_identity")==reference.Identity()&&certificate.get<std::string>("source_identity")==source,"certificate source/reference binding differs");
            auto layout=ReadLayout(certificate.get_child("normalization.layout"));PerturbationSigns signs;
            signs.velocity=request.get<int>("velocity_sign",1);signs.k=request.get<int>("k_sign",1);signs.epsilon=request.get<int>("epsilon_sign",1);
            Perturbation seed(reference,layout,xs,rs,signs);PerturbationNormalization normalization;
            normalization.layout=layout;normalization.evidenceSha256=identity;
            normalization.unitPeakLower=certificate.get<double>("normalization.unitPeakLower");normalization.unitPeakEstimate=certificate.get<double>("normalization.unitPeakEstimate");normalization.unitPeakUpper=certificate.get<double>("normalization.unitPeakUpper");seed.SetNormalization(normalization);
            bool twice=false;try{seed.SetNormalization(normalization);}catch(const std::exception&){twice=true;}Require(twice,"normalization was mutable");
            const int n=mixed.VelocitySize();Vec velocity=Vec::Zero(2*n+mixed.PressureSize());
            for(int axis=0;axis<2;++axis)velocity.segment(axis*n,n)=mixed.Project([&](double x,double z){return seed.Velocity(x,z)[axis==0?0:2];},axis);
            Vec dk=k.Project([&](double x,double z){return seed.KIncrement(x,z).value;});
            Vec de=epsilon.Project([&](double x,double z){return seed.IndependentEpsilonIncrement(x,z).value;});
            Vec kBase=k.Project([&](double,double){return p.KIn;}),kNew=kBase+dk;
            auto traceK=[&](double x){const Element*owner=nullptr;for(const auto&e:flow.Elements())if(Close(e.rl,p.H)&&x>=e.xl-1e-12&&x<=e.xr+1e-12){owner=&e;break;}Require(owner,"missing outer H trace owner");Vec xx(1),zz(1);xx[0]=x;zz[0]=p.d0+p.H;return k.Evaluate(factory.Points(*owner,xx,zz),kNew)[0];};
            Vec gammaBase=epsilon.EssentialLift([&](double,double z,const Vector&){return Ain/(z-p.d0);},factory);
            Vec gammaNew=epsilon.EssentialLift([&](double x,double z,const Vector&normal){return normal[2]<-.5?reference.Amplitude(traceK(x))/p.H:Ain/(z-p.d0);},factory);
            Vec gammaIncrement=gammaNew-gammaBase,totalE=de+gammaIncrement;
            const double peakK=traceK(0.);
            bool cornerConflict=false;try{epsilon.EssentialLift([&](double,double z,const Vector&normal){return Ain/(z-p.d0)+(normal[0]<-.5?1e-5*Ain/p.H:0.);},factory);}catch(const std::exception&){cornerConflict=true;}Require(cornerConflict,"conflicting Gamma corner was accepted");
            Vec divergence=flow.Divergence()*velocity.head(2*n);Eigen::SimplicialLDLT<Sparse>pressureRiesz;pressureRiesz.compute(flow.PressureMass());Require(pressureRiesz.info()==Eigen::Success,"pressure Riesz failed");
            const double divRms=std::sqrt(std::max(0.,divergence.dot(pressureRiesz.solve(divergence))/seed.DomainArea()));
            Vec velocityLift=mixed.Lift(velocity),kLift=k.Lift(dk),epsilonLift=epsilon.Lift(de);
            out<<",\"certificate_sha256\":\""<<identity<<"\",\"base_mesh_xml_sha256\":\""<<certificate.get<std::string>("mesh_xml_sha256")<<"\",\"layout\":";LayoutJSON(out,layout);
            out<<",\"normalization\":{\"unitPeakLower\":"<<normalization.unitPeakLower<<",\"unitPeakEstimate\":"<<normalization.unitPeakEstimate<<",\"unitPeakUpper\":"<<normalization.unitPeakUpper<<",\"streamfunction_coefficient\":"<<seed.StreamfunctionCoefficient()<<"},\"signs\":["<<signs.velocity<<','<<signs.k<<','<<signs.epsilon<<"],\"support_area\":"<<seed.SupportArea()<<",\"outer_support_area\":"<<seed.OuterSupportArea()<<",\"support_intervals\":["<<seed.StreamwiseIntervals()<<','<<seed.NormalIntervals()<<','<<seed.OuterNormalIntervals()<<"],\"normalization_load_twice_rejected\":true,\"Gamma_corner_conflict_rejected\":true,\"synthetic_Gamma_base\":\"Kbase=Kin only; not an evolved physical reference\",\"native_divergence_dual_RMS\":"<<divRms<<",\"free_velocity_ids\":[";List(out,mixed.FreeVelocityIds()[0]);out<<',';List(out,mixed.FreeVelocityIds()[1]);out<<"],\"k_free_ids\":";List(out,k.FreeIds());out<<",\"epsilon_free_ids\":";List(out,epsilon.FreeIds());
            out<<",\"synthetic_interface_peak\":{\"Kbase\":"<<p.KIn<<",\"Knew\":"<<peakK<<",\"cutoff_base\":"<<reference.Cutoff(p.KIn)<<",\"cutoff_new\":"<<reference.Cutoff(peakK)<<",\"Gamma_base\":"<<Ain/p.H<<",\"Gamma_new_target\":"<<reference.Amplitude(peakK)/p.H<<'}';
            for(const auto&item:std::vector<std::pair<std::string,const Vec*>>{{"velocity",&velocity},{"k_increment",&dk},{"epsilon_independent_increment",&de},{"synthetic_k_base",&kBase},{"synthetic_k_new",&kNew},{"gamma_base_lift",&gammaBase},{"gamma_new_lift",&gammaNew},{"gamma_increment_lift",&gammaIncrement},{"epsilon_total_increment",&totalE},{"divergence_action",&divergence},{"velocity_essential_lift",&velocityLift},{"k_essential_lift",&kLift},{"epsilon_independent_essential_lift",&epsilonLift}}){out<<",\""<<item.first<<"\":";VectorJSON(out,*item.second);}
            out<<",\"helper_samples\":[";bool first=true;
            for(int i=0;i<4;++i)for(int j=0;j<4;++j)for(double fraction:{.21132486540518713,.78867513459481287}){
                const double x=layout.xKnots[i]+fraction*(layout.xKnots[i+1]-layout.xKnots[i]),r=layout.rKnots[j]+fraction*(layout.rKnots[j+1]-layout.rKnots[j]);const double z=p.d0+r;
                const auto u=seed.UnitVelocity(x,z);const auto g=seed.UnitGradient(x,z);const auto a=seed.KIncrement(x,z),b=seed.IndependentEpsilonIncrement(x,z);
                if(!first)out<<',';first=false;out<<"{\"x\":"<<x<<",\"r\":"<<r<<",\"z\":"<<z<<",\"unit_velocity\":";List(out,u);out<<",\"unit_gradient\":[";for(int c=0;c<3;++c){if(c)out<<',';List(out,g[c]);}out<<"],\"k_jet\":["<<a.value<<','<<a.dx<<','<<a.dr<<"],\"epsilon_jet\":["<<b.value<<','<<b.dx<<','<<b.dr<<"]}";
            }out<<']';
        }
        out<<",\"elements\":";ElementsJSON(out,flow,km,em);out<<"}\n";Create(argv[4],out.str());session->Finalise();
        std::cout<<"Recorded perturbation "<<mode<<" on "<<flow.Elements().size()<<" native elements\n";
    }catch(const std::exception&e){std::cerr<<"H27 perturbation probe: "<<e.what()<<'\n';return 2;}
    return 0;
}
