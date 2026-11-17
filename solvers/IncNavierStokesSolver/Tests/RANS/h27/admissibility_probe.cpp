// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeAdmissibility.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>

using namespace Nektar::WJNative;
namespace J=Nektar::RANSWallSolveJson;
namespace
{
Tensor Direction(double theta)
{
    const Vector n{std::cos(theta),0,std::sin(theta)},l{-std::sin(theta),0,std::cos(theta)};
    Tensor h{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)h[i][j]=l[i]*n[j];return h;
}
double Contract(const Tensor &a,const Tensor &b)
{ double v=0;for(int i=0;i<3;++i)for(int j=0;j<3;++j)v+=a[i][j]*b[i][j];return v; }
struct State{std::string name;double r;Tensor g;double nu=6.423173803526447e-5;};
}
int main(int argc,char **argv)
{
    try
    {
        if(argc!=2)throw std::invalid_argument("usage: admissibility_probe new-output.json");
        if(std::filesystem::exists(argv[1]))throw std::invalid_argument("immutable output exists");
        const double k=.008809015145515192,e=.0023155240547992535,nu=6.423173803526447e-5,H=.15;
        const double pi=std::acos(-1.);std::vector<State> states;
        const std::vector<std::array<double,4>> gradients{{0,0,0,0},{0,.9,0,0},{0,25,0,0},
            {.3,.7,-.2,-.3},{.4,.2,.1,.1},{-.3,-.7,.2,.3},{1e-12,.9,-1e-12,1e-12}};
        for(double ratio:{0.,.4,.400001,.5,.6,.749999,.75,1.,2.})
            for(std::size_t j=0;j<gradients.size();++j)
            {
                const auto &a=gradients[j];Tensor g{};g[0][0]=a[0];g[0][2]=a[1];g[2][0]=a[2];g[2][2]=a[3];
                states.push_back({"r"+std::to_string(ratio)+"_g"+std::to_string(j),H*ratio,g});
            }
        Tensor base{};base[0][2]=25;
        const double baseAngle=CheckPlanarState(k,e,base,nu,.003,H,H).principalAngle;
        for(double offset:{0.,1e-12,1e-7})
        {
            const double alpha=pi/2-baseAngle+offset,c=std::cos(alpha),s=std::sin(alpha);
            const Tensor rotation{{{c,0,-s},{0,1,0},{s,0,c}}};Tensor rotated{};
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int l=0;l<3;++l)for(int m=0;m<3;++m)
                rotated[i][j]+=rotation[i][l]*base[l][m]*rotation[j][m];
            states.push_back({"stationary_phi_near_pi_"+J::Number(offset),H,rotated});
        }
        // Constitutive unit-test control, not an application parameter change:
        // q=10 is past the simple-shear stress maximum. With this small nu,
        // the admissible rank-one vertical shear tangent is strictly negative.
        Tensor negative{};negative[0][2]=10*e/k;
        states.push_back({"negative_principal_unit_control",H,negative,1e-8});
        std::vector<std::string> records;int negativeCovariance=0,negativePrincipal=0;
        for(const auto &state:states)
        {
            const double auxiliary=.00005+.02*state.r,stateNu=state.nu;
            const auto wall=AdaptAmended(k,e,state.g,{0,0,1},stateNu,auxiliary,state.r,H);
            const auto result=CheckPlanarState(k,e,state.g,stateNu,auxiliary,state.r,H);
            const auto covarianceOnly=CheckPlanarState(k,e,state.g,stateNu,auxiliary,state.r,H,false);
            if(covarianceOnly.principalEvaluated || covarianceOnly.covarianceOverK!=result.covarianceOverK)
                throw std::runtime_error("covariance-only API contract differs");
            negativeCovariance+=result.covarianceOverK<0;negativePrincipal+=result.principalOverMu<0;
            std::vector<double> grid;constexpr int angles=192;
            for(int j=0;j<angles;++j)
            {
                const auto h=Direction(pi*j/angles);
                const auto a=ApplyAmendedDirection(k,e,state.g,{0,0,1},stateNu,auxiliary,state.r,H,0,0,h,0);
                grid.push_back(Contract(h,a.V)/wall.mu);
            }
            const auto minimumDirection=Direction(result.principalAngle);
            const auto minimumAction=ApplyAmendedDirection(k,e,state.g,{0,0,1},stateNu,auxiliary,state.r,H,0,0,minimumDirection,0);
            Tensor dg{};dg[0][0]=.13*e/k;dg[0][2]=-.27*e/k;dg[2][0]=.19*e/k;dg[2][2]=-.07*e/k;
            const double dk=.01*k,de=-.02*e,da=.03*auxiliary,dr=state.r==0?0:.01*H;
            const auto action=ApplyAmendedDirection(k,e,state.g,{0,0,1},stateNu,auxiliary,state.r,H,dk,de,dg,da,dr);
            records.push_back(J::Object({{"name",J::Quote(state.name)},{"k",J::Number(k)},{"epsilon",J::Number(e)},
                {"nu",J::Number(stateNu)},{"auxiliary",J::Number(auxiliary)},{"r",J::Number(state.r)},{"H",J::Number(H)},
                {"gradient",J::Rows(state.g)},{"R",J::Rows(wall.R)},{"mu",J::Number(wall.mu)},{"chi",J::Number(wall.chi)},
                {"covariance_over_k",J::Number(result.covarianceOverK)},
                {"principal_over_mu",J::Number(result.principalOverMu)},{"principal_raw_over_mu",J::Number(result.principalRawOverMu)},
                {"principal_reserve_over_mu",J::Number(result.principalReserveOverMu)},
                {"theta",J::Number(result.principalAngle)},{"stationary_residual",J::Number(result.stationaryResidual)},
                {"stationary_candidates",J::Number(result.stationaryCandidates)},
                {"native_at_minimum",J::Number(Contract(minimumDirection,minimumAction.V)/wall.mu)},
                {"grid",J::Numbers(grid)},{"direction",J::Object({{"dk",J::Number(dk)},{"de",J::Number(de)},
                    {"da",J::Number(da)},{"dr",J::Number(dr)},{"dg",J::Rows(dg)},{"R",J::Rows(action.R)},
                    {"V",J::Rows(action.V)},{"mu",J::Number(action.mu)},{"production",J::Number(action.production)}})}}));
        }
        Tensor unsupported{};unsupported[1][0]=.1;bool rejected=false;
        try{CheckPlanarState(k,e,unsupported,nu,.0001,0,H);}catch(const std::exception&){rejected=true;}
        if(!rejected)throw std::runtime_error("unsupported 3D gradient accepted");
        const std::string result=J::Object({{"schema",J::Quote("WJH27PlanarAdmissibilityProbeV1")},
            {"scope",J::Quote("Current planar frozen-scalar principal and full covariance; no 3D principal claim")},
            {"status",J::Quote("RECORDED")},{"negative_covariance_count",J::Number(negativeCovariance)},
            {"negative_principal_count",J::Number(negativePrincipal)},{"unsupported_spanwise_rejected","true"},
            {"states",J::List(records)}});
        std::ofstream output(argv[1]);output<<result<<'\n';if(!output.good())throw std::runtime_error("write failed");
        std::cout<<"planar admissibility states recorded: "<<records.size()<<'\n';
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 2;}
    return 0;
}
