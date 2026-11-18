// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeReference.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>

namespace J=Nektar::RANSWallSolveJson;
using Nektar::WJNative::Reference;

int main(int argc,char **argv)
{
    try
    {
        if(argc!=3&&argc!=4)throw std::invalid_argument("usage: reference_probe qualified-table.json new-output.json [policy]");
        if(std::filesystem::exists(argv[2]))throw std::invalid_argument("immutable probe output exists");
        const auto policy=Reference::ParsePolicy(argc==4?argv[3]:"Historical");
        Reference reference(argv[1],policy);const auto &p=reference.Parameters();
        const bool exact=policy==Nektar::WJNative::ReferencePolicy::ExactInnerFiniteNuV1;
        std::vector<std::string> dynamic,inlet,rejections;
        for(int i=0;i<=16;++i)
        {
            const double K=i==0?p.KMinimum:i==16?p.KMaximum:p.KIn*(.9+.2*i/16);
            const double cutoff=reference.Cutoff(K);
            std::vector<double> radii{0.,.1*p.H,.4*p.H,.6*p.H,.75*p.H,p.H};
            if(exact)for(double r:{cutoff,std::nextafter(cutoff,p.H),p.roughness})radii.push_back(r);
            for(double r:radii)
            {
                const auto a=reference.Auxiliary(K,r),e=reference.DependentEpsilon(K,r);
                dynamic.push_back(J::Object({{"K",J::Number(K)},{"r",J::Number(r)},
                    {"A",J::Number(reference.Amplitude(K))},{"dA",J::Number(reference.AmplitudeDerivative(K))},
                    {"ddA",J::Number(reference.AmplitudeSecondDerivative(K))},
                    {"c",J::Number(cutoff)},{"dc",J::Number(reference.CutoffDerivative(K))},
                    {"K_at_c",J::Number(reference.KAtCutoff(cutoff))},
                    {"aux",J::Numbers(std::array<double,3>{a.value,a.dK,a.dr})},
                    {"epsilon",J::Numbers(std::array<double,3>{e.value,e.dK,e.dr})}}));
            }
        }
        std::vector<double> radii;
        for(int i=0;i<=96;++i)radii.push_back(p.H*i/96);
        for(int i=1;i<=64;++i)radii.push_back(i==64?p.topRadius:p.H*std::exp(std::log(p.topRadius/p.H)*i/64));
        const double oldCutoff=reference.PrescribedInletCutoff();
        for(double r:{oldCutoff,std::nextafter(oldCutoff,0.),std::nextafter(oldCutoff,p.H),
                      std::nextafter(p.H,0.),p.H*(1-1e-10),p.H*1e-12})radii.push_back(r);
        if(exact)for(double r:{.4*p.H,.75*p.H,p.roughness,std::nextafter(p.H,p.topRadius)})radii.push_back(r);
        for(double r:radii)
        {
            const auto traction=reference.OutletTraction(r);
            inlet.push_back(J::Object({{"r",J::Number(r)},{"U",J::Number(reference.InletU(r))},
                {"k",J::Number(reference.InletK(r))},{"U_r",J::Number(reference.InletUGradient(r))},
                {"k_r",J::Number(reference.InletKGradient(r))},
                {"epsilon",r>=p.H?J::Number(reference.InletE(r)):"null"},
                {"epsilon_r",r>=p.H?J::Number(reference.InletEGradient(r)):"null"},
                {"Pi",exact?J::Number(reference.InletModifiedPressure(r)):"null"},
                {"traction",J::Numbers(traction)}}));
        }
        auto rejected=[&](const char *name,const auto &operation)
        {
            bool didReject=false;std::string message;
            try{operation();}catch(const std::exception &e){didReject=true;message=e.what();}
            if(!didReject)throw std::runtime_error(std::string("missing rejection: ")+name);
            rejections.push_back(J::Object({{"case",J::Quote(name)},{"reason",J::Quote(message)}}));
        };
        rejected("K below interval",[&]{reference.Cutoff(std::nextafter(p.KMinimum,0.));});
        rejected("K above interval",[&]{reference.Cutoff(std::nextafter(p.KMaximum,std::numeric_limits<double>::infinity()));});
        rejected("nonfinite K",[&]{reference.Amplitude(std::numeric_limits<double>::quiet_NaN());});
        rejected("negative radius",[&]{reference.Auxiliary(p.KIn,-1e-15);});
        rejected("above top",[&]{reference.OutletTraction(std::nextafter(p.topRadius,std::numeric_limits<double>::infinity()));});
        rejected("transported inner inlet epsilon",[&]{reference.InletE(.5*p.H);});
        rejected("dependent outer epsilon",[&]{reference.DependentEpsilon(p.KIn,2*p.H);});
        rejected("wrong table identity",[&]{Reference other(argv[1],std::string(64,'0'));});
        const auto result=J::Object({{"schema",J::Quote("WJH27NativeReferenceProbeV1")},
            {"status",J::Quote("RECORDED")},{"table_sha256",J::Quote(reference.Identity())},
            {"historical_inlet_cutoff",J::Number(reference.PrescribedInletCutoff())},
            {"reference_policy",J::Quote(reference.PolicyName())},
            {"top_tangential_traction",J::Number(reference.TopTangentialTraction())},
            {"top_epsilon_outward_flux",J::Number(reference.TopEpsilonOutwardFlux())},
            {"dynamic",J::List(dynamic)},{"inlet",J::List(inlet)},{"rejections",J::List(rejections)}});
        std::ofstream output(argv[2]);output<<result<<'\n';
        if(!output.good())throw std::runtime_error("probe output write failed");
        std::cout<<"native reference probe recorded\n";
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 2;}
    return 0;
}
