// SPDX-License-Identifier: MIT
// Native local volume rows for supplied original signed owners. No PDE
// trajectory, initial-history import or application forcing is performed.
#include "../h26/common.h"
#include <IncNavierStokesSolver/EquationSystems/WJNativeMixedFlow.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeRules.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeReference.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeCutoffPanels.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeWallV2.h>
void Require(bool ok,const char*message){if(!ok)throw std::runtime_error(message);}
std::vector<double> Numbers(const pt::ptree&tree,const std::string&key)
{std::vector<double>v;for(const auto&e:tree.get_child(key)){double x=e.second.get_value<double>();Require(std::isfinite(x),"nonfinite probe input");v.push_back(x);}return v;}
struct Trace{int owner;double left,right;std::array<double,4>c;};
int main(int argc,char**argv)
{
    try{
        Require(argc==5&&!fs::exists(argv[4]),"usage: create-only cutoff_panel_probe case.xml table.json request.json output.json");
        pt::ptree request;pt::read_json(argv[3],request);
        std::string xml=argv[1],opt=xml+".opt",program="h27-cutoff-panels",option="--use-opt-file";
        char*args[]{program.data(),option.data(),opt.data(),xml.data(),nullptr};
        auto session=LibUtilities::SessionReader::CreateInstance(4,args);auto graph=SpatialDomains::MeshGraphIO::Read(session);
        std::vector<MultiRegions::ContFieldSharedPtr>fields;
        for(std::string name:{"u","v","p"})fields.push_back(MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,name,true,false));
        Geometry geometry;geometry.H=.15;geometry.h=.10854503542535461;geometry.d0=.0875;
        MixedFlow flow(fields,geometry);RuleFactory factory(flow.Core(),8,12);Reference reference(argv[2]);
        const auto&p=reference.Parameters();std::vector<Trace>traces;
        for(const auto&item:request.get_child("traces")){
            const auto&j=item.second;Trace t;t.owner=j.get<int>("owner");t.left=j.get<double>("left");t.right=j.get<double>("right");
            auto c=Numbers(j,"coefficients");Require(c.size()==4&&t.right>t.left,"invalid supplied trace");std::copy(c.begin(),c.end(),t.c.begin());traces.push_back(t);
        }
        auto trace=[&](const Element&e)->const Trace&{
            const Trace*found=nullptr;
            for(const auto&t:traces)if(std::abs(t.left-e.xl)<1e-10&&std::abs(t.right-e.xr)<1e-10){Require(!found,"ambiguous supplied column owner");found=&t;}
            Require(found,"missing supplied column owner");return *found;
        };
        auto K=[&](const Element&e,double x){const auto&t=trace(e);double s=(x-t.left)/(t.right-t.left);return(1-s)*t.c[0]+s*t.c[1]+s*(1-s)*(t.c[2]+(4*s-2)*t.c[3]);};
        const auto range=reference.CutoffRadiusRange();
        CutoffPartition partition{range,[&](const Element&e){
            const auto&t=trace(e);std::vector<double>levels;
            for(double face:{e.rl,e.rr})if(face>=range[0]&&face<=range[1])levels.push_back(reference.KAtCutoff(face));
            auto roots=CubicCutoffPanels(t.c,levels);std::vector<double>x;
            const double eps=64*std::numeric_limits<double>::epsilon()*std::max({std::abs(e.xl),std::abs(e.xr),e.xr-e.xl});
            for(double s:roots){double point=t.left+(t.right-t.left)*s;if(point>e.xl+eps&&point<e.xr-eps)x.push_back(point);}return x;
        }};
        auto rules=factory.Cutoff([&](const Element&e,double x){return reference.Cutoff(K(e,x));},&partition);
        std::ofstream out(argv[4]);out<<std::setprecision(17)<<"{\"schema\":\"WJH27NativeCutoffPanelRowsV1\",\"source_identity\":\""<<WJ_H27_SOURCE_ID
            <<"\",\"policy\":\"CubicCutoffFacePanelsV1\",\"cells\":[";bool comma=false;
        for(const auto&item:request.get_child("cells")){
            const auto&j=item.second;int id=j.get<int>("owner");const Element*element=nullptr;const NativeRule*rule=nullptr;
            for(int i=0;i<int(flow.Elements().size());++i)if(flow.Elements()[i].geometry==id){element=&flow.Elements()[i];rule=rules[i].get();}
            Require(element&&rule,"requested original owner is absent");const auto&e=*element;const auto&r=*rule;
            auto bounds=Numbers(j,"bounds"),u=Numbers(j,"signed_u"),v=Numbers(j,"signed_v"),k=Numbers(j,"signed_k"),signs=Numbers(j,"k_signs");
            Require(bounds.size()==4&&u.size()==16&&v.size()==16&&k.size()==16&&signs.size()==16,"invalid local signed inventory");
            for(int i=0;i<4;++i)Require(std::abs(bounds[i]-std::array<double,4>{e.xl,e.xr,e.rl,e.rr}[i])<1e-12,"native/source owner bounds differ");
            Vec cu(16),cv(16),ck(16);Mat b=r.B,dx=r.Dx,dz=r.Dz;
            for(int i=0;i<16;++i){cu[i]=u[i];cv[i]=v[i];ck[i]=k[i];b.col(i)/=e.us[i];dx.col(i)/=e.us[i];dz.col(i)/=e.us[i];}
            Vec U=b*cu,W=b*cv,kval=b*ck,Ux=dx*cu,Uz=dz*cu,Wx=dx*cv,Wz=dz*cv,kx=dx*ck,kz=dz*ck;
            Vec residual=Vec::Zero(16),sink=Vec::Zero(16);double area=0;
            for(int q=0;q<r.w.size();++q){
                const double radius=r.z[q]-geometry.d0,interfaceK=K(e,r.x[q]),E=reference.Ew(interfaceK,radius),aux=reference.Aux(interfaceK,radius);
                Tensor G{};G[0][0]=Ux[q];G[0][2]=Uz[q];G[2][0]=Wx[q];G[2][2]=Wz[q];
                auto wall=AdaptAmended(kval[q],E,G,{0,0,1},p.nu,aux,radius,p.H,ConstitutiveBackend::FastFull3DV1);
                const double D=p.nu+p.Cd*kval[q]*kval[q]/(p.sigmaK*E);Require(r.w[q]>0&&D>0,"nonpositive current quadrature/diffusion");area+=r.w[q];
                for(int i=0;i<16;++i){const double phi=b(q,i)*signs[i],px=dx(q,i)*signs[i],pz=dz(q,i)*signs[i];
                    residual[i]+=r.w[q]*(D*(px*kx[q]+pz*kz[q])-(px*U[q]+pz*W[q])*kval[q]-phi*(wall.production-E));sink[i]+=r.w[q]*phi*E;}
            }
            if(comma)out<<',';comma=true;out<<"{\"owner\":"<<id<<",\"points\":"<<r.w.size()<<",\"area\":"<<area<<",\"k_volume\":";VectorJSON(out,residual);out<<",\"Ew_load\":";VectorJSON(out,sink);out<<'}';
        }
        const auto slope=reference.LogCutoffDerivativeBounds();out<<"],\"cutoff_range\":["<<range[0]<<','<<range[1]<<"],\"log_slope_bounds\":["<<slope[0]<<','<<slope[1]
            <<"],\"reference_identity\":\""<<reference.Identity()<<"\",\"inverse_cases\":[";
        for(int i=0;i<=20;++i){double radius=i==0?range[0]:i==20?range[1]:range[0]+(range[1]-range[0])*i/20;
            double value=reference.KAtCutoff(radius);if(i)out<<',';out<<"{\"radius\":"<<radius<<",\"K\":"<<value<<",\"cutoff\":"<<reference.Cutoff(value)<<",\"derivative\":"<<reference.CutoffDerivative(value)<<'}';}
        bool rejected=false;try{reference.KAtCutoff(std::nextafter(range[0],-std::numeric_limits<double>::infinity()));}catch(const std::invalid_argument&){rejected=true;}
        Require(rejected,"out-of-range inverse was clamped");
        out<<"],\"out_of_range_inverse_rejected\":true,\"full_H27_qualification\":false,\"TASK3_READY\":false}\n";
        out.close();Require(bool(out),"native cutoff report write failed");session->Finalise();return 0;
    }catch(const std::exception&e){std::cerr<<"Native cutoff panel diagnostic: "<<e.what()<<'\n';return 1;}
}
