// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeConstitutive.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>

using namespace Nektar::WJNative;
using namespace Nektar::SolverUtils;
namespace J=Nektar::RANSWallSolveJson;
namespace
{
struct State
{
    std::string name;double k,e,dk,de;RANSTensor g{},dg{};WJParameters p;
};
std::string Closure(const RANSClosureValues&v)
{
    return J::Object({{"cMu",J::Number(v.cMu)},{"eddyViscosity",J::Number(v.eddyViscosity)},
        {"scalarViscosity",J::Number(v.scalarViscosity)},{"strainInvariant",J::Number(v.strainInvariant)},
        {"rotationInvariant",J::Number(v.rotationInvariant)},{"divergence",J::Number(v.divergence)},
        {"production",J::Number(v.production)},{"productionDev",J::Number(v.productionDev)},
        {"kSource",J::Number(v.kSource)},{"epsilonProduction",J::Number(v.epsilonProduction)},
        {"epsilonDestruction",J::Number(v.epsilonDestruction)},{"epsilonSource",J::Number(v.epsilonSource)},
        {"momentumDiffusion",J::Number(v.momentumDiffusion)},{"kDiffusion",J::Number(v.kDiffusion)},
        {"epsilonDiffusion",J::Number(v.epsilonDiffusion)},{"yap",J::Number(v.yap)},
        {"R",J::Rows(v.reynoldsStress)},{"anisotropy",J::Rows(v.anisotropy)},
        {"deviatoricStress",J::Rows(v.deviatoricStress)},{"stressRemainder",J::Rows(v.stressRemainder)},
        {"quadraticStress",J::Rows(v.quadraticStress)},{"cubicStress",J::Rows(v.cubicStress)},
        {"root",J::Number(v.wj.n)},{"s",J::Number(v.wj.s)},{"w",J::Number(v.wj.w)},
        {"i4",J::Number(v.wj.i4)},{"i5",J::Number(v.wj.i5)},{"Q",J::Number(v.wj.q)},
        {"rootDerivative",J::Number(v.wj.rootDerivative)},{"rootResidual",J::Number(v.wj.rootResidual)},
        {"beta",J::Numbers(v.wj.beta)}});
}
std::string Status(const FastWJStatus&s)
{
    return J::Object({{"fallback",J::Quote(FastWJFallbackName(s.fallback))},
        {"root_iterations",J::Number(s.rootIterations)},
        {"normalized_root_residual",J::Number(s.normalizedRootResidual)},
        {"relative_root_condition",J::Number(s.relativeRootCondition)}});
}
double Clock(){return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count();}
}
int main(int argc,char**argv)
{
    try
    {
        if(argc!=2&&argc!=3)throw std::invalid_argument("usage: constitutive_probe new-output.json [benchmark-calls=20000]");
        if(std::filesystem::exists(argv[1]))throw std::invalid_argument("immutable output exists");
        const int calls=argc==3?std::stoi(argv[2]):20000;
        if(calls<=0)throw std::invalid_argument("positive benchmark calls required");
        std::vector<State>states;const double Kin=.008809015145515192,Ein=.0023155240547992535;
        auto append=[&](std::string name,double k,double e,const RANSTensor&g,int variant=0)
        {
            State s{name,k,e,.13*k,-.07*e,g};
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)s.dg[i][j]=.2*(e/k)*std::sin(1.+3*i+j);
            if(variant==1){s.p.cd*=.5;s.p.sigmaK=1.2;s.p.sigmaEpsilon=1.7;}
            if(variant==2){s.p.cd*=2;s.p.cEpsilon1=1.5;s.p.cEpsilon2=1.9;}
            states.push_back(s);
        };
        for(double factor:{.9,1.,1.1})for(double q:{0.,1e-10,.1,1.,3.386780140528718,5.538461538461538,10.,30.,100.})
        {RANSTensor g{};g[0][2]=q*Ein/Kin;append("parallel_"+J::Number(factor)+"_"+J::Number(q),Kin*factor,Ein*std::pow(factor,1.5),g);}
        for(double ratio:{.4,.5,.6,.75,1.})
        {RANSTensor g{};g[0][2]=(3.4+.5*(1-ratio))*Ein/(Kin*ratio);g[0][0]=.0001;g[2][2]=-.00009;
         append("collar_weak_div_"+J::Number(ratio),Kin*(1+.2*(1-ratio)),Ein/ratio,g);}
        RANSTensor rotation{};rotation[0][1]=1.7;rotation[1][0]=-1.7;rotation[1][2]=.3;rotation[2][1]=-.3;
        append("zero_strain_3D_rotation",Kin,Ein,rotation);
        RANSTensor dilation{};for(int i=0;i<3;++i)dilation[i][i]=.1;
        append("isotropic_raw_divergence",Kin,Ein,dilation);
        std::mt19937_64 random(279031);std::uniform_real_distribution<double>uniform(-1.,1.);
        for(double magnitude:{.01,.5,3.,20.})for(int j=0;j<12;++j)
        {RANSTensor g{};for(auto&row:g)for(double&v:row)v=magnitude*(Ein/Kin)*uniform(random);
         append("full3D_"+J::Number(magnitude)+"_"+std::to_string(j),Kin,Ein,g,j%3);}
        const std::size_t regularCount=states.size();
        RANSTensor small{};small[0][2]=1e-10;append("fallback_small_k_e",1e-130,1e-140,small);
        RANSTensor large{};large[0][2]=1e10;append("fallback_large_k_e",1e100,1e110,large);
        RANSTensor huge{};huge[0][2]=1e40;append("fallback_large_normalized_tensor",1.,1e-10,huge);
        std::vector<std::string>records;int valueFallbacks=0,directionFallbacks=0;
        for(const auto&s:states)
        {
            const auto old=EvaluateWJClosure(s.k,s.e,s.g,s.p,WJGradientPolicy::DevStrainRawWorkV1,0);
            const auto oldDirection=ApplyWJClosureDirection(s.k,s.e,s.g,s.dk,s.de,s.dg,s.p,WJGradientPolicy::DevStrainRawWorkV1,0);
            FastWJStatus valueStatus,directionStatus;
            const auto value=EvaluateFastWJClosure(s.k,s.e,s.g,s.p,&valueStatus);
            const auto direction=ApplyFastWJClosureDirection(s.k,s.e,s.g,s.dk,s.de,s.dg,s.p,&directionStatus);
            const auto pair=EvaluateFastWJValueAndDirection(s.k,s.e,s.g,s.dk,s.de,s.dg,s.p);
            valueFallbacks+=valueStatus.UsedFallback();directionFallbacks+=directionStatus.UsedFallback();
            records.push_back(J::Object({{"name",J::Quote(s.name)},{"k",J::Number(s.k)},{"epsilon",J::Number(s.e)},
                {"dk",J::Number(s.dk)},{"de",J::Number(s.de)},{"G",J::Rows(s.g)},{"dG",J::Rows(s.dg)},
                {"parameters",J::Numbers(std::array<double,5>{s.p.cd,s.p.cEpsilon1,s.p.cEpsilon2,s.p.sigmaK,s.p.sigmaEpsilon})},
                {"old",Closure(old)},{"fast",Closure(value)},{"old_direction",Closure(oldDirection)},
                {"fast_direction",Closure(direction)},{"paired_value",Closure(pair.value)},
                {"paired_direction",Closure(pair.direction)},{"value_status",Status(valueStatus)},
                {"direction_status",Status(directionStatus)},{"paired_status",Status(pair.status)}}));
        }
        int rejected=0;auto reject=[&](const auto&operation)
        {try{operation();}catch(const std::exception&){++rejected;return;}throw std::runtime_error("invalid input accepted");};
        reject([&]{EvaluateFastWJClosure(-Kin,Ein,{});});
        reject([&]{EvaluateFastWJClosure(Kin,0,{});});
        reject([&]{auto g=rotation;g[0][0]=std::numeric_limits<double>::infinity();EvaluateFastWJClosure(Kin,Ein,g);});
        reject([&]{ApplyFastWJClosureDirection(Kin,Ein,rotation,std::numeric_limits<double>::quiet_NaN(),0,{});});
        reject([&]{WJParameters p;p.c1=1.9;EvaluateFastWJClosure(Kin,Ein,rotation,p);});
        volatile double sink=0;
        for(int j=0;j<100;++j){const auto&s=states[j%regularCount];sink+=EvaluateFastWJClosure(s.k,s.e,s.g,s.p).production;}
        auto benchmark=[&](int kind)
        {
            int fallbacks=0;const double start=Clock();
            for(int j=0;j<calls;++j)
            {
                const auto&s=states[j%regularCount];RANSClosureValues value;FastWJStatus status;
                if(kind==0)value=EvaluateWJClosure(s.k,s.e,s.g,s.p,WJGradientPolicy::DevStrainRawWorkV1,0);
                if(kind==1)value=EvaluateFastWJClosure(s.k,s.e,s.g,s.p,&status);
                if(kind==2)value=ApplyWJClosureDirection(s.k,s.e,s.g,s.dk,s.de,s.dg,s.p,WJGradientPolicy::DevStrainRawWorkV1,0);
                if(kind==3)value=ApplyFastWJClosureDirection(s.k,s.e,s.g,s.dk,s.de,s.dg,s.p,&status);
                fallbacks+=status.UsedFallback();sink+=value.production+value.reynoldsStress[0][0]+value.eddyViscosity;
            }
            return J::Object({{"seconds",J::Number(Clock()-start)},{"calls",J::Number(calls)},
                              {"fallbacks",J::Number(fallbacks)}});
        };
        const auto oldValueTime=benchmark(0),fastValueTime=benchmark(1),oldDirectionTime=benchmark(2),fastDirectionTime=benchmark(3);
        const auto result=J::Object({{"schema",J::Quote("WJH27FastConstitutiveProbeV1")},{"status",J::Quote("RECORDED")},
            {"scope",J::Quote("Pure 3D WJ kernel equivalence and kernel-only microbenchmark; no coupled cost qualification")},
            {"value_fallbacks",J::Number(valueFallbacks)},{"direction_fallbacks",J::Number(directionFallbacks)},
            {"invalid_inputs_rejected",J::Number(rejected)},{"states",J::List(records)},
            {"benchmark",J::Object({{"old_value",oldValueTime},{"fast_value",fastValueTime},
                {"old_direction",oldDirectionTime},{"fast_direction",fastDirectionTime},{"volatile_checksum",J::Number(sink)}})}});
        std::ofstream output(argv[1]);output<<result<<'\n';if(!output.good())throw std::runtime_error("probe write failed");
        std::cout<<"Recorded "<<states.size()<<" full WJ states and kernel-only timings\n";
    }
    catch(const std::exception&error){std::cerr<<error.what()<<'\n';return 2;}
    return 0;
}
