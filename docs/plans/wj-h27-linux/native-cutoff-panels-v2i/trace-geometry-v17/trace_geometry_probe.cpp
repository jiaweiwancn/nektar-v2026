// SPDX-License-Identifier: MIT
// Diagnostic only: compare physical-coordinate H points and an exact native
// horizontal face with the independently written Modified_A edge polynomial.
#include "common.h"
#include <IncNavierStokesSolver/EquationSystems/WJNativeMixedFlow.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeRules.h>

std::array<double,2> Edge(const std::array<double,4>&c,double t,double width)
{
    const double bubble=t*(1-t),linear=c[2]+(4*t-2)*c[3];
    return{(1-t)*c[0]+t*c[1]+bubble*linear,
        (-c[0]+c[1]+(1-2*t)*linear+4*bubble*c[3])/width};
}
int main(int argc,char**argv)
{
    try{
        if(argc!=3||fs::exists(argv[2]))throw std::runtime_error("usage: create-only WJH27TraceGeometryProbe case.xml result.json");
        std::string xml=argv[1],opt=xml+".opt",program="h27-trace-geometry",option="--use-opt-file";
        char*args[]{program.data(),option.data(),opt.data(),xml.data(),nullptr};
        auto session=LibUtilities::SessionReader::CreateInstance(4,args);auto graph=SpatialDomains::MeshGraphIO::Read(session);
        std::vector<MultiRegions::ContFieldSharedPtr>fields;
        for(std::string name:{"u","v","p"})fields.push_back(MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,name,true,false));
        Geometry geometry;geometry.H=.15;geometry.h=.10854503542535461;geometry.d0=.0875;
        MixedFlow flow(fields,geometry);RuleFactory rules(flow.Core(),8,12);
        Vec synthetic(flow.VelocitySize());for(int j=0;j<synthetic.size();++j)synthetic[j]=std::cos(.7*(j+1));
        std::ofstream out(argv[2]);if(!out)throw std::runtime_error("cannot create trace geometry report");
        out<<std::setprecision(17)<<"{\"schema\":\"WJH27TraceGeometryDiagnosticV1\",\"diagnostic_only\":true,\"acceptance_limit\":2e-13,\"rows\":[";
        bool comma=false;double pointValue=0,pointGradient=0,faceValue=0,faceGradient=0;
        for(const auto&e:flow.Elements())if(std::abs(e.rl-geometry.H)<1e-12){
            Vec local(16);std::array<double,4>coeff;
            for(int j=0;j<16;++j)local[j]=synthetic[e.u[j]];
            for(int j=0;j<4;++j)coeff[j]=synthetic[e.u[j]]*e.us[j];
            auto observe=[&](const char*kind,const NativeRule&r,int i,double requestedX){
                const double width=e.xr-e.xl,t=(requestedX-e.xl)/width;auto exact=Edge(coeff,t,width);
                const double value=(r.B.row(i)*local)[0],gradient=(r.Dx.row(i)*local)[0];
                const double ev=std::abs(value-exact[0]),eg=std::abs(gradient-exact[1])*width;
                if(kind[0]=='p'){pointValue=std::max(pointValue,ev);pointGradient=std::max(pointGradient,eg);}
                else{faceValue=std::max(faceValue,ev);faceGradient=std::max(faceGradient,eg);}
                if(comma)out<<',';comma=true;
                out<<"{\"kind\":\""<<kind<<"\",\"owner\":"<<e.geometry<<",\"x\":"<<r.x[i]<<",\"requested_x\":"<<requestedX<<",\"z\":"<<r.z[i]
                   <<",\"bounds\":["<<e.xl<<','<<e.xr<<','<<e.rl<<','<<e.rr<<"],\"t\":"<<t
                   <<",\"eta_at_requested_H\":"<<2*(geometry.d0+geometry.H-geometry.d0-(e.rl+e.rr)/2)/(e.rr-e.rl)
                   <<",\"native_value\":"<<value<<",\"native_gradient\":"<<gradient<<",\"edge_value\":"<<exact[0]
                   <<",\"edge_gradient\":"<<exact[1]<<",\"value_error\":"<<ev<<",\"gradient_error_times_width\":"<<eg<<'}';
            };
            for(double t:{0.,.135,.555,.8,1.}){
                const double x=e.xl+t*(e.xr-e.xl);auto rule=rules.Points(e,Vec::Constant(1,x),Vec::Constant(1,geometry.d0+geometry.H));observe("physical_H_point",rule,0,x);
            }
            auto face=rules.Face(e,2,16);for(int i=0;i<face.x.size();++i)observe("native_lower_face",face,i,face.x[i]);
        }
        out<<"],\"physical_point_value_MAX\":"<<pointValue<<",\"physical_point_gradient_scaled_MAX\":"<<pointGradient
           <<",\"native_face_value_MAX\":"<<faceValue<<",\"native_face_gradient_scaled_MAX\":"<<faceGradient
           <<",\"physical_point_gate_pass\":"<<(std::max(pointValue,pointGradient)<=2e-13?"true":"false")
           <<",\"native_face_gate_pass\":"<<(std::max(faceValue,faceGradient)<=2e-13?"true":"false")
           <<",\"full_H27_qualification\":false,\"TASK3_READY\":false}\n";
        out.close();if(!out)throw std::runtime_error("trace geometry report write failed");session->Finalise();return 0;
    }catch(const std::exception&e){std::cerr<<"H27 trace geometry diagnostic: "<<e.what()<<'\n';return 1;}
}
