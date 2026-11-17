// SPDX-License-Identifier: MIT
// Continuous verification only: no mesh, projection, time step or application.
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupledMMS.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeWallV2.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <fcntl.h>
#include <unistd.h>
#include <cerrno>
#include <cstring>
#include <iostream>

using namespace Nektar::WJNative;
namespace J=Nektar::RANSWallSolveJson;

std::string Jet(const MMSJet &v)
{
    return J::Numbers(std::array<double,7>{v.value,v.dx,v.dr,v.dt,v.dxx,v.dxr,v.drr});
}

void WriteExclusive(const char *path,const std::string &text)
{
    int fd=::open(path,O_WRONLY|O_CREAT|O_EXCL,0600);
    if(fd<0)throw std::runtime_error(std::string("immutable MMS probe output: ")+std::strerror(errno));
    std::size_t offset=0;
    while(offset<text.size())
    {
        const ssize_t count=::write(fd,text.data()+offset,text.size()-offset);
        if(count<0&&errno==EINTR)continue;
        if(count<=0){::close(fd);throw std::runtime_error("MMS probe output write failed");}
        offset+=std::size_t(count);
    }
    if(::close(fd)!=0)throw std::runtime_error("MMS probe output close failed");
}

int main(int argc,char **argv)
{
    try
    {
        if(argc<3||(argc-3)%2!=0)throw std::invalid_argument("usage: WJH27CoupledMMSProbe qualified-table.json new-output.json [--body-amplitude VALUE] [--variant Historical|PolynomialScalarsV1]");
        double bodyAmplitude=0;
        MMSVariant variant=MMSVariant::Historical;
        bool seenBody=false,seenVariant=false;
        for(int argument=3;argument<argc;argument+=2)
        {
            const std::string option=argv[argument],value=argv[argument+1];
            if(option=="--body-amplitude"&&!seenBody)
            {
                std::size_t consumed=0;bodyAmplitude=std::stod(value,&consumed);
                if(consumed!=value.size())throw std::invalid_argument("invalid MMS body amplitude");
                seenBody=true;
            }
            else if(option=="--variant"&&!seenVariant)
            {
                if(value=="Historical")variant=MMSVariant::Historical;
                else if(value=="PolynomialScalarsV1")variant=MMSVariant::PolynomialScalarsV1;
                else throw std::invalid_argument("unknown MMS variant");
                seenVariant=true;
            }
            else throw std::invalid_argument("unknown or repeated MMS probe option");
        }
        Reference reference(argv[1]);const auto &p=reference.Parameters();
        Geometry geometry;geometry.d0=p.d0;geometry.H=p.H;geometry.h=p.roughness;
        CoupledMMS mms(reference,geometry,.5,bodyAmplitude,variant);
        std::vector<std::string> states,boundaries;
        const Vector nx{1,0,0},nr{0,0,1};
        auto point=[&](const std::string &label,double x,double r,double time)
        {
            const double z=geometry.d0+r;
            const auto f=mms.Fields(x,z,time);
            const auto K=mms.InterfaceK(x,time);
            const auto force=mms.Forcing(x,z,time);
            const auto auxiliary=reference.Auxiliary(K.value,r);
            const auto stress=AdaptAmended(f.k.value,force.destruction,f.gradient,nr,
                                           p.nu,auxiliary.value,r,p.H);
            const auto epsilon=force.outer?Jet(mms.OuterEpsilon(x,z,time)):"null";
            states.push_back(J::Object({
                {"label",J::Quote(label)},{"x",J::Number(x)},{"r",J::Number(r)},
                {"z",J::Number(z)},{"time",J::Number(time)},{"outer",J::Boolean(force.outer)},
                {"u",Jet(f.u)},{"w",Jet(f.w)},{"Pi",Jet(f.Pi)},{"k",Jet(f.k)},
                {"K",Jet(K)},{"epsilon",epsilon},
                {"gradient",J::Rows(f.gradient)},{"gradient_x",J::Rows(f.gradientX)},
                {"gradient_r",J::Rows(f.gradientR)},
                {"A",J::Number(reference.Amplitude(K.value))},
                {"cutoff",J::Number(reference.Cutoff(K.value))},
                {"auxiliary",J::Number(auxiliary.value)},
                {"R",J::Rows(stress.R)},{"V",J::Rows(stress.V)},
                {"mu",J::Number(stress.mu)},{"chi",J::Number(stress.chi)},
                {"production",J::Number(force.production)},
                {"destruction",J::Number(force.destruction)},
                {"diffusion_k",J::Number(force.diffusionK)},
                {"diffusion_e",J::Number(force.diffusionE)},
                {"momentum_forcing",J::Numbers(force.momentum)},
                {"k_forcing",J::Number(force.k)},{"epsilon_forcing",J::Number(force.epsilon)},
                {"awf",J::Numbers(force.awf)},{"awf_at_H",J::Numbers(mms.AWFAtH(x,time))},
                {"physical_body",J::Numbers(force.body)},
                {"physical_body_at_H",J::Numbers(mms.BodyAtH(x,time))},
                {"k_outward_flux_x_r",J::Numbers(std::array<double,2>{
                    mms.KOutwardFlux(x,z,time,nx),mms.KOutwardFlux(x,z,time,nr)})},
                {"epsilon_outward_flux_x_r",force.outer?J::Numbers(std::array<double,2>{
                    mms.EpsilonOutwardFlux(x,z,time,nx),mms.EpsilonOutwardFlux(x,z,time,nr)}):"null"}}));
        };
        auto boundary=[&](const char *label,double x,double r,double time,const Vector &normal)
        {
            const double z=geometry.d0+r;
            boundaries.push_back(J::Object({{"label",J::Quote(label)},
                {"x",J::Number(x)},{"r",J::Number(r)},{"time",J::Number(time)},
                {"normal",J::Numbers(normal)},
                {"k_flux",J::Number(mms.KOutwardFlux(x,z,time,normal))},
                {"epsilon_flux",r>=p.H?J::Number(mms.EpsilonOutwardFlux(x,z,time,normal)):"null"}}));
        };
        for(double time:{0.,.13,.41})
        {
            for(double x:{.137,.61})
            {
                const double c=reference.Cutoff(mms.InterfaceK(x,time).value);
                for(const auto &[label,r]:std::vector<std::pair<std::string,double>>{
                    {"below_cutoff",.5*c},{"near_cutoff_below",.95*c},
                    {"near_cutoff_above",1.05*c},{"chi_zero",.3*p.H},
                    {"collar_lower",.4*p.H},{"collar_middle",.5*p.H},
                    {"roughness_below",p.roughness-1e-4},{"roughness_above",p.roughness+1e-4},
                    {"collar_upper",.75*p.H},{"H_minus",p.H-1e-4},
                    {"H_trace",p.H},{"H_plus",p.H+1e-4},
                    {"outer_middle",.3},{"outer_high",.47}})
                    point(label,x,r,time);
                boundary("ground",x,0,time,Vector{0,0,-1});
                boundary("top",x,.5,time,nr);
                boundary("epsilon_H",x,p.H,time,Vector{0,0,-1});
            }
            boundary("left",0,.3,time,Vector{-1,0,0});
            boundary("right",1,.3,time,nx);
        }
        const auto result=J::Object({
            {"schema",J::Quote("WJH27ContinuousCoupledMMSProbeV1")},
            {"status",J::Quote("RECORDED")},{"full_coupled_gate_pass","false"},
            {"admissibility_qualified","false"},
            {"formula_identity",J::Quote(CoupledMMS::Identity)},
            {"mms_variant",J::Quote(mms.VariantName())},
            {"body_formula_identity",J::Quote(CoupledMMS::BodyIdentity)},
            {"body_amplitude",J::Number(mms.BodyAmplitude())},
            {"model_identity",J::Quote(WallModelV2)},
            {"table_sha256",J::Quote(reference.Identity())},
            {"jet_order",J::List({J::Quote("value"),J::Quote("dx"),J::Quote("dr"),
                J::Quote("dt"),J::Quote("dxx"),J::Quote("dxr"),J::Quote("drr")})},
            {"source_formula",J::Object({
                {"s",J::Quote("sqrt(Kin)*(1+.012*sin(.6*t)+.008*x*cos(.8*t))")},
                {"k",J::Quote(variant==MMSVariant::Historical?
                    "s^2*(1+.005*(1+.2*sin(.7*t))*r^2*(r-H)^2/(H^2*R^2))":
                    "s^2*(1+.005*(1+.2*sin(.7*t))*r*(r-H)/(H*R))")},
                {"epsilon_outer",J::Quote(variant==MMSVariant::Historical?
                    "s^3/(cl*r)*(1+.01*sin(.5*t)*(r-H)^2/R^2)":
                    "s^3/(cl*H)*(1+.01*sin(.5*t)*(r-H)^2/R^2)")},
                {"psi",J::Quote(".02*(1+.1*cos(.9*t))*r^2/2+.05*(1+.2*sin(.7*t))*x^2*(1-x)*r^2*(R-r)")},
                {"Pi",J::Quote(".0002*(1+.3*sin(.8*t))*(x^2+r^2/R^2-2/3)")},
                {"momentum",J::Quote("Ut+U.grad(U)-div(V)+grad(Pi)-body-1_(r<=h)*(2*uH*uxH+Pi_xH-2*k_xH/3-bxH)*ex")},
                {"physical_body",J::Quote("body_amplitude*(2e-4*(1+x+.5*r)*cos(.3*t),0,3e-4*(1+r)*sin(.2*t))")},
                {"physical_body_at_H",J::Quote("physical_body(x,H,t); prescribed independently of manufactured momentum load")}})},
            {"parameters",J::Object({{"Kin",J::Number(p.KIn)},{"H",J::Number(p.H)},
                {"h",J::Number(p.roughness)},{"d0",J::Number(p.d0)},{"R",J::Number(.5)},
                {"body_amplitude",J::Number(mms.BodyAmplitude())},
                {"nu",J::Number(p.nu)},{"Cd",J::Number(p.Cd)},
                {"sigmaK",J::Number(p.sigmaK)},{"sigmaE",J::Number(p.sigmaE)},
                {"cl",J::Number(p.cLength)},{"qe",J::Number(p.qe)},
                {"kappa",J::Number(p.kappa)},{"tau0",J::Number(p.tau0)}})},
            {"states",J::List(states)},{"boundary_fluxes",J::List(boundaries)}});
        WriteExclusive(argv[2],result+"\n");
        std::cout<<"continuous MMS probe recorded "<<states.size()<<" states\n";
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 2;}
    return 0;
}
