// SPDX-License-Identifier: MIT
#include "../h26/common.h"
#include <IncNavierStokesSolver/EquationSystems/WJNativeMixedFlow.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeFlowTangent.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <cmath>
#include <limits>

namespace J=Nektar::RANSWallSolveJson;
namespace
{
constexpr double Kin=.008809015145515192,Ein=.0023155240547992535,Nu=6.423173803526447e-5;
void Check(bool value,const char*message){if(!value)throw std::runtime_error(message);}
double Field(int component,double x,double r,bool direction)
{
    if(direction)return component==0?.02+.01*x-.013*r+.007*x*x*r-.004*x*r*r:
                                     -.017+.009*x+.011*r+.005*x*r+.002*x*x*r*r;
    return component==0?.03+.07*x+.11*r+.03*x*r+.017*x*x*r:
                        .01-.015*x+.023*r+.006*x*r+.009*x*r*r;
}
std::array<double,3> Frozen(double x,double r,double H)
{return {Kin*(1+.02*x+.03*r),Ein*(1+.01*x+.07*r),r<H?.0003+.001*r+.00001*x:0.};}
Vec Local(const Vec&v,const NativeRule&r,int offset)
{Vec out(16);for(int j=0;j<16;++j)out[j]=v[offset+r.u[j]];return out;}
std::vector<Tensor> Gradients(const Vec&v,const NativeRule&r,int n)
{
    std::vector<Tensor>out(r.w.size());
    for(int c=0;c<2;++c){Vec local=Local(v,r,c*n),dx=r.Dx*local,dz=r.Dz*local;
        for(int q=0;q<r.w.size();++q){out[q][c?2:0][0]=dx[q];out[q][c?2:0][2]=dz[q];}}
    return out;
}
Vec StressAction(const Rules&rules,int n,const Vec&velocity,const Geometry&geo,const Vec*direction=nullptr)
{
    Vec out=Vec::Zero(2*n);
    for(const auto&pointer:rules){const auto&r=*pointer;auto gradient=Gradients(velocity,r,n);
        std::vector<Tensor>dg;if(direction)dg=Gradients(*direction,r,n);
        std::array<Vec,4>stress;for(auto&v:stress)v=Vec(r.w.size());
        for(int q=0;q<r.w.size();++q){const double radius=r.z[q]-geo.d0;const auto f=Frozen(r.x[q],radius,geo.H);
            const Tensor value=direction?ApplyAmendedDirection(f[0],f[1],gradient[q],{0,0,1},Nu,f[2],radius,geo.H,
                0,0,dg[q],0).V:AdaptAmended(f[0],f[1],gradient[q],{0,0,1},Nu,f[2],radius,geo.H).V;
            for(int i=0;i<2;++i)for(int j=0;j<2;++j)stress[2*i+j][q]=r.w[q]*value[i?2:0][j?2:0];}
        for(int i=0;i<2;++i){Vec local=r.Dx.transpose()*stress[2*i]+r.Dz.transpose()*stress[2*i+1];
            for(int j=0;j<16;++j)out[i*n+r.u[j]]+=local[j];}}
    return out;
}
std::vector<OwnedPlanarVelocityTangents> Samples(const Rules&rules,int n,const Vec&velocity,
    const Geometry&geo,ConstitutiveBackend backend)
{
    std::vector<OwnedPlanarVelocityTangents>out;
    for(const auto&pointer:rules){const auto&r=*pointer;OwnedPlanarVelocityTangents cell;cell.geometry=r.geometry;
        const auto gradient=Gradients(velocity,r,n);
        for(int q=0;q<r.w.size();++q){const double radius=r.z[q]-geo.d0;const auto f=Frozen(r.x[q],radius,geo.H);
            cell.points.push_back(EvaluatePlanarVelocityTangent(f[0],f[1],gradient[q],Nu,f[2],radius,geo.H,backend));}
        out.push_back(std::move(cell));}
    return out;
}
std::string TangentJSON(const FrozenPlanarVelocityTangent&t)
{std::vector<std::string>out;for(const auto&v:t.dVdG)out.push_back(J::Rows(v));return J::List(out);}
}

int main(int argc,char**argv)
{
    try
    {
        Check(argc==3,"usage: WJH27FlowTangentProbe case.xml new-output.json");
        const std::string xml=argv[1],output=argv[2],opt=output+".opt";
        Check(!fs::exists(output)&&!fs::exists(opt),"immutable tangent probe output/options already exist");
        {std::ofstream f(opt);f<<"<NEKTAR><COLLECTIONS DEFAULT=\"NoCollection\"/></NEKTAR>\n";Check(bool(f),"cannot create tangent probe options");}
        std::string program="h27-flow-tangent",option="--use-opt-file",input=xml,options=opt;
        char*args[]{program.data(),option.data(),options.data(),input.data(),nullptr};
        auto session=LibUtilities::SessionReader::CreateInstance(4,args);auto graph=SpatialDomains::MeshGraphIO::Read(session);
        std::vector<MultiRegions::ContFieldSharedPtr>fields;
        for(std::string name:{"u","v","p"})fields.push_back(MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,name,true,false));
        Geometry geo;geo.d0=.0875;geo.H=.15;geo.h=.10854503542535461;
        MixedFlow mixed(fields,geo);const auto&flow=mixed.Core();const int n=flow.VelocitySize();
        Check(n<=512,"compact tangent probe is restricted to at most 512 velocity coefficients");
        RuleFactory factory(flow,6,8);const auto rules=factory.Cutoff([](const Element&,double x){return .025*(1+.01*x);});
        Vec velocity(2*n),direction(2*n);
        for(int c=0;c<2;++c){velocity.segment(c*n,n)=mixed.Project([=](double x,double z){return Field(c,x,z-geo.d0,false);},c);
            direction.segment(c*n,n)=mixed.Project([=](double x,double z){return Field(c,x,z-geo.d0,true);},c);}
        auto historical=Samples(rules,n,velocity,geo,ConstitutiveBackend::Historical);
        const auto fast=Samples(rules,n,velocity,geo,ConstitutiveBackend::FastFull3DV1);
        const Sparse matrix=AssemblePlanarVelocityTangent(rules,n,historical),fastMatrix=AssemblePlanarVelocityTangent(rules,n,fast);
        const Vec assembled=matrix*direction,direct=StressAction(rules,n,velocity,geo,&direction);
        const double actionScale=std::max(direct.cwiseAbs().maxCoeff(),1e-30);
        Check((assembled-direct).cwiseAbs().maxCoeff()/actionScale<2e-12,"four-column tangent assembly differs from combined exact direction");
        Check((matrix-fastMatrix).norm()/std::max(matrix.norm(),1e-30)<2e-11,"tangent backend parity failed");
        // A declared algebraic change of signed basis, derived from the actual
        // native rules, exercises negative signs even on an all-positive mesh.
        Vec sign(n),signedDirection=direction;
        for(int j=0;j<n;++j){sign[j]=j%3==1?-1.:1.;signedDirection[j]*=sign[j];signedDirection[n+j]*=sign[j];}
        Rules signedRules;
        for(const auto&rp:rules){auto r=std::make_shared<NativeRule>(*rp);for(int j=0;j<16;++j){const double s=sign[r->u[j]];
            r->us[j]*=s;r->B.col(j)*=s;r->Dx.col(j)*=s;r->Dz.col(j)*=s;}signedRules.push_back(r);}
        const Sparse signedMatrix=AssemblePlanarVelocityTangent(signedRules,n,historical);
        const Vec signedAction=signedMatrix*signedDirection;
        Vec expectedSigned=assembled;for(int j=0;j<2*n;++j)expectedSigned[j]*=sign[j%n];
        Check((signedAction-expectedSigned).cwiseAbs().maxCoeff()/actionScale<2e-12,"signed-basis tangent covariance failed");
        std::vector<std::string>fd;
        for(double step:{2e-5,1e-5}){Vec value=(StressAction(rules,n,velocity+step*direction,geo)-
            StressAction(rules,n,velocity-step*direction,geo))/(2*step);
            // VectorJSON inherits its stream precision. The nested stream must
            // preserve doubles just as the main output stream does below.
            std::ostringstream encoded;encoded<<std::setprecision(std::numeric_limits<double>::max_digits10);
            VectorJSON(encoded,value);
            std::ostringstream legacy;VectorJSON(legacy,value);
            std::vector<std::string>raw;for(int j=0;j<value.size();++j)raw.push_back(J::Number(value[j]));
            const double difference=(value-assembled).cwiseAbs().maxCoeff();
            const double scale=std::max(assembled.cwiseAbs().maxCoeff(),1e-30);
            fd.push_back(J::Object({{"step",J::Number(step)},{"action",encoded.str()},
                {"native_scalar_serialized_action",J::List(raw)},
                {"native_in_memory_max_abs_error",J::Number(difference)},
                {"native_in_memory_relative_error",J::Number(difference/scale)},
                {"legacy_precision6_action",legacy.str()}}));}
        std::vector<std::string>points;
        for(double ratio:{0.,.4,.400001,.55,.749999,.75,1.,2.})
        {
            const double radius=ratio*geo.H;const auto f=Frozen(.37,radius,geo.H);
            Tensor g{};g[0][0]=.11;g[0][2]=.8;g[2][0]=-.13;g[2][2]=-.07;
            auto tangent=EvaluatePlanarVelocityTangent(f[0],f[1],g,Nu,f[2],radius,geo.H);
            auto wall=AdaptAmended(f[0],f[1],g,{0,0,1},Nu,f[2],radius,geo.H);
            points.push_back(J::Object({{"name",J::Quote("weak_div_"+J::Number(ratio))},{"k",J::Number(f[0])},
                {"epsilon",J::Number(f[1])},{"nu",J::Number(Nu)},{"auxiliary",J::Number(f[2])},{"r",J::Number(radius)},
                {"H",J::Number(geo.H)},{"G",J::Rows(g)},{"mu",J::Number(wall.mu)},{"dVdG",TangentJSON(tangent)}}));
        }
        Tensor negative{};negative[0][2]=10*Ein/Kin;
        const auto negativeTangent=EvaluatePlanarVelocityTangent(Kin,Ein,negative,1e-8,0,geo.H,geo.H);
        Check(negativeTangent.dVdG[1][0][2]<0,"scientifically negative Newton shear tangent was modified/rejected");
        points.push_back(J::Object({{"name",J::Quote("negative_shear_unit_control")},{"k",J::Number(Kin)},
            {"epsilon",J::Number(Ein)},{"nu",J::Number(1e-8)},{"auxiliary",J::Number(0)},{"r",J::Number(geo.H)},
            {"H",J::Number(geo.H)},{"G",J::Rows(negative)},
            {"mu",J::Number(AdaptAmended(Kin,Ein,negative,{0,0,1},1e-8,0,geo.H,geo.H).mu)},
            {"dVdG",TangentJSON(negativeTangent)}}));
        std::vector<std::string>rejections;
        auto reject=[&](const std::string&name,const auto&operation){try{operation();}catch(const std::exception&e){
            rejections.push_back(J::Object({{"name",J::Quote(name)},{"reason",J::Quote(e.what())}}));return;}
            throw std::runtime_error("invalid tangent input accepted: "+name);};
        reject("owner_mismatch",[&]{auto s=historical;++s[0].geometry;AssemblePlanarVelocityTangent(rules,n,s);});
        reject("point_count",[&]{auto s=historical;s[0].points.pop_back();AssemblePlanarVelocityTangent(rules,n,s);});
        reject("nonfinite_tangent",[&]{auto s=historical;s[0].points[0].dVdG[0][0][0]=std::numeric_limits<double>::quiet_NaN();AssemblePlanarVelocityTangent(rules,n,s);});
        reject("negative_weight",[&]{auto r=std::make_shared<NativeRule>(*rules[0]);r->w[0]=-r->w[0];Rules bad=rules;bad[0]=r;AssemblePlanarVelocityTangent(bad,n,historical);});
        reject("duplicate_global_id",[&]{auto r=std::make_shared<NativeRule>(*rules[0]);r->u[1]=r->u[0];Rules bad=rules;bad[0]=r;AssemblePlanarVelocityTangent(bad,n,historical);});
        reject("spanwise_gradient",[&]{auto g=negative;g[1][1]=.01;EvaluatePlanarVelocityTangent(Kin,Ein,g,Nu,0,geo.H,geo.H);});
        reject("nonpositive_k",[&]{EvaluatePlanarVelocityTangent(-Kin,Ein,negative,Nu,0,geo.H,geo.H);});
        std::ofstream out(output);out<<std::setprecision(17)<<"{\"schema\":\"WJH27FrozenPlanarFlowTangentProbeV1\",\"status\":\"RECORDED\",\"velocity_size\":"<<n
            <<",\"d0\":"<<geo.d0<<",\"H\":"<<geo.H<<",\"nu\":"<<Nu<<",\"Kin\":"<<Kin<<",\"Ein\":"<<Ein
            <<",\"scope\":\"Frozen-scalar stress Jacobian only; no full coupled Newton convergence or performance claim\",\"velocity\":";VectorJSON(out,velocity);
        out<<",\"direction\":";VectorJSON(out,direction);out<<",\"assembled_action\":";VectorJSON(out,assembled);
        out<<",\"combined_direction_action\":";VectorJSON(out,direct);out<<",\"matrix\":";MatrixJSON(out,matrix);
        out<<",\"basis_sign_transform\":";VectorJSON(out,sign);out<<",\"sign_transformed_action\":";VectorJSON(out,signedAction);
        out<<",\"fast_matrix\":";MatrixJSON(out,fastMatrix);out<<",\"finite_difference_serialization\":\"MaxDigits10V2\",\"finite_differences\":"<<J::List(fd)<<",\"points\":"<<J::List(points)<<",\"rejections\":"<<J::List(rejections)<<",\"rules\":[";
        for(std::size_t i=0;i<rules.size();++i){if(i)out<<',';const auto&r=*rules[i];const Element*owner=nullptr;
            for(const auto&e:flow.Elements())if(e.geometry==r.geometry)owner=&e;Check(owner!=nullptr,"rule owner missing in probe");
            out<<"{\"geometry\":"<<r.geometry<<",\"bounds\":"<<J::Numbers(std::array<double,4>{owner->xl,owner->xr,owner->rl,owner->rr})
                <<",\"u\":"<<J::Numbers(r.u)<<",\"us\":"<<J::Numbers(r.us)<<",\"x\":";VectorJSON(out,r.x);
            out<<",\"z\":";VectorJSON(out,r.z);out<<",\"w\":";VectorJSON(out,r.w);out<<'}';}
        out<<"]}\n";out.close();Check(!out.fail(),"tangent probe output failed");session->Finalise();
        std::cout<<"Recorded native owned-rule frozen velocity tangent and controls\n";
    }
    catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 2;}
    return 0;
}
