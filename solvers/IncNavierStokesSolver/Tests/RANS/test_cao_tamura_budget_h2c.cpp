// SPDX-License-Identifier: MIT
#include <SolverUtils/Turbulence/RANSAnalyticalNeutralInflow.h>
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <SolverUtils/Turbulence/RANSTwoRegionHill.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <fstream>
#include <iostream>
using namespace Nektar::SolverUtils;
namespace J=Nektar::RANSWallSolveJson;
void Check(bool yes,const char *message){if(!yes)throw std::runtime_error(message);}
int main(int argc,char **argv)
{
    try
    {
        Check(argc==3,"usage: H2cChecks profiles.csv output.json");
        const auto model=RANSNeutralReference::CLS();const auto wall=CalibratedAWFParameters(model);
        const double nu=.04/10500.,us=.051/.6867027875750893,inf=us/.051,K=model.Q()*us*us,kap=model.Kappa();
        struct Row{double z,U,w;};std::vector<Row> rows;std::ifstream in(argv[1]);std::string line;std::getline(in,line);
        while(std::getline(in,line))
        {std::vector<std::string> c;std::istringstream ss(line);std::string a;while(std::getline(ss,a,','))c.push_back(a);
         if(c[0]=="f3a_R" && std::stod(c[16])>0)rows.push_back({.04*std::stod(c[12]),inf*std::stod(c[14]),std::stod(c[16])});}
        std::vector<std::string> fits;double z0=0,z0Sensitivity=0;
        for(double minimum:{.01,.015})
        {
            long double weight=0,intercept=0;int count=0;
            for(const auto &r:rows)if(r.z>=minimum && r.z<=.051){weight+=r.w;intercept+=r.w*(r.U-us/kap*std::log(r.z));++count;}
            intercept/=weight;const double rough=std::exp(-intercept/(us/kap));if(minimum==.01)z0=rough;else z0Sensitivity=rough;
            long double error=0;std::vector<std::string> points;
            for(const auto &r:rows)if(r.z>=minimum && r.z<=.051)
            {const double e=us/kap*std::log(r.z/rough)-r.U;error+=r.w*e*e;
             points.push_back(J::Numbers(std::vector<double>{r.z,r.U,r.w,e,e/inf}));}
            Check(count==(minimum==.01?7:6),"unexpected frozen mean-fit membership");
            fits.push_back(J::Object({{"minimum_z",J::Number(minimum)},{"maximum_z","0.051"},{"count",J::Number(count)},
                {"z0_eff",J::Number(rough)},{"RMSE_over_Uinf",J::Number(std::sqrt(error/weight)/inf)},{"points_z_U_weight_residual_residualOverUinf",J::List(points)}}));
        }
        std::vector<std::string> heights;
        for(double H:{.0055,.006,.0065})
        {
            const double target=us/kap*std::log(H/z0);double left=H*1e-8,right=H*(1-1e-9);
            auto f=[&](double h){return us*us*AnalyticalWallFunction({nu,H,h,H/2},K,wall).Integrals().j0-target;};
            const double fl=f(left),fr=f(right);Check(fl*fr<0,"roughness root is not bracketed");
            int iterations=0;for(;iterations<90;++iterations){const double m=(left+right)/2;if(f(m)*fl>0)left=m;else right=m;if(right-left<1e-16)break;}
            const double h=(left+right)/2;
            double sl=H*1e-8,sr=H*(1-1e-9);const double target2=us/kap*std::log(H/z0Sensitivity);
            auto f2=[&](double r){return us*us*AnalyticalWallFunction({nu,H,r,H/2},K,wall).Integrals().j0-target2;};
            const double sf=f2(sl);Check(sf*f2(sr)<0,"sensitivity roughness root unbracketed");
            for(int j=0;j<90;++j){const double m=(sl+sr)/2;if(sf*f2(m)>0)sl=m;else sr=m;}
            AWFProfileOptions o;o.formulation=AWFProfileFormulation::CLSCalibratedBudgetAWF;o.quadratureOrder=32;
            EquilibriumAWFProfiles p({nu,H,h,H,H},K,model,{},o,wall);const auto &a=p.Original();
            const double EC=p.Dissipation(H);long double P=0,E=0;double minK=1e100;
            for(const auto &q:p.Quadrature()){P+=q.weight*p.ReferenceProduction(q.r);E+=q.weight*p.Dissipation(q.r);minK=std::min(minK,p.ReferenceK(q.r));}
            std::vector<std::string> convergence;
            for(int order:{16,32,64})
            {auto options=o;options.quadratureOrder=order;EquilibriumAWFProfiles v({nu,H,h,H,H},K,model,{},options,wall);
             long double pr=0,di=0;for(const auto &q:v.Quadrature()){pr+=q.weight*v.ReferenceProduction(q.r);di+=q.weight*v.Dissipation(q.r);}
             convergence.push_back(J::Numbers(std::vector<double>{double(order),v.Reference().cutoff,double(pr),double(di),v.ReferenceK(0)}));}
            std::vector<std::string> curved;
            for(double curvature:{-19.73920880217871,0.,19.73920880217871})for(double ratio:{.8,1.,1.2})
            {
                EquilibriumAWFProfiles v({nu,H,h,H,H},ratio*K,model,{curvature,0,.1},o,wall);
                curved.push_back(J::Numbers(std::vector<double>{curvature,ratio,v.Reference().cutoff,v.Reference().margin,
                    v.Reference().budgetResidual,v.ReferenceK(0),v.Dissipation(H)}));
            }
            const RANSAnalyticalNeutralInflow ref(model,us,H,z0,1,{-.2,0,0},{1,0,0},{0,0,1});
            double neutralError=0,strainTrace=0,symmetry=0,workError=0;
            for(double z:{H,.04,.1,.8})
            {
                const auto v=ref.AtHeight(z);RANSTensor g{};g[0][2]=v.dU;
                const auto c=EvaluateRANSClosure(RANSModel::CLSHighRe,v.k,v.epsilon,g,z);
                RANSRobustVariant reg;reg.beta=.4;reg.gamma=v.gamma;reg.referenceStrain=v.strain;
                const auto r=EvaluateRANSRobustClosure(RANSModel::CLSHighRe,v.k,v.epsilon,g,z,{},0,reg);
                neutralError=std::max(neutralError,std::abs(c.eddyViscosity-v.nuT));
                for(int i=0;i<3;++i){strainTrace+=v.strain[i][i];for(int j=0;j<3;++j){symmetry=std::max(symmetry,std::abs(v.strain[i][j]-v.strain[j][i]));neutralError=std::max(neutralError,std::abs(c.reynoldsStress[i][j]-r.reynoldsStress[i][j]));}}
                double w=0;for(int i=0;i<3;++i)for(int j=0;j<3;++j)w-=r.reynoldsStress[i][j]*g[i][j];workError=std::max(workError,std::abs(w-r.production));
            }
            Check(neutralError<1e-12 && workError<1e-12 && strainTrace==0 && symmetry==0,"neutral analytical provider/law consistency");
            Check(minK>0 && std::abs(double(P-E))<1e-11 && std::abs(f(h))<1e-11,"inner reference budget/root");
            RANSTwoRegionHill geometry(.04,.1,H,.8,-.2,.6);
            heights.push_back(J::Object({{"H",J::Number(H)},{"h_s",J::Number(h)},{"sensitivity_h_s",J::Number((sl+sr)/2)},
                {"quadrature_columns",J::Quote("order,cutoff,production,dissipation,reference_ground_k")},{"quadrature",J::List(convergence)},{"root_bracket",J::Numbers(std::vector<double>{H*1e-8,H*(1-1e-9),fl,fr})},
                {"root_residual",J::Number(f(h))},{"root_iterations",J::Number(iterations)},
                {"U_H",J::Number(target)},{"K0",J::Number(K)},{"E_C_H",J::Number(EC)},
                {"cutoff",J::Number(p.Reference().cutoff)},{"budget_P_E",J::Numbers(std::vector<double>{double(P),double(E)})},
                {"min_inner_reference_k",J::Number(minK)},{"k_ground",J::Number(p.ReferenceK(0))},
                {"viscosity_intercept",J::Number(a.ViscosityInterceptHeight())},
                {"H_plus",J::Number(us*H/nu)},{"displaced_H_plus_diagnostic",J::Number(us*(H-.0035)/nu)},
                {"ratios_H_z0_b_hs",J::Numbers(std::vector<double>{H/z0,H/.005,H/h})},
                {"molecular_fraction",J::Number(nu/a.TotalViscosity(H))},{"clearance_to_7mm",J::Number(.007-H)},
                {"reach_ratio",J::Number(H*geometry.CurvatureBound())},
                {"inner_outer_shear",J::Numbers(std::vector<double>{a.EddyViscosity(H)*us*us/a.TotalViscosity(H),us*us})},
                {"matching_k_flux",J::Number(p.ReferenceKFlux(H))},
                {"curved_columns",J::Quote("curvature,K/K0,cutoff,margin,budgetResidual,k_ground,E_C_H")},{"curved",J::List(curved)},
                {"provider_error",J::Number(neutralError)},{"stress_work_error",J::Number(workError)}}));
        }
        std::ofstream out(argv[2]);out<<J::Object({{"q_star",J::Number(model.Q())},{"kappa_CLS",J::Number(kap)},
            {"cMu_wall",J::Number(wall.cMu)},{"cLength_wall",J::Number(wall.cLength)},
            {"u_tau",J::Number(us)},{"nu",J::Number(nu)},{"Uinf",J::Number(inf)},
            {"fit",J::List(fits)},{"heights",J::List(heights)},
            {"implied_U_Hhill",J::Number(us/kap*std::log(.04/z0))},{"implied_Re_H",J::Number(.04*us/kap*std::log(.04/z0)/nu)},
            {"physical_steps","0"},{"status",J::Quote("passed")}})<<'\n';Check(out.good(),"output failed");return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
