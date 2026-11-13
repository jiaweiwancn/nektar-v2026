// SPDX-License-Identifier: MIT
#include <SolverUtils/Turbulence/RANSAnalyticalNeutralInflow.h>
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <SolverUtils/Turbulence/RANSRobustVariant.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <boost/property_tree/json_parser.hpp>
#include <fstream>
#include <iostream>
#include <memory>
using namespace Nektar::SolverUtils;
namespace J=Nektar::RANSWallSolveJson;
void check(bool a,const char *s){if(!a)throw std::runtime_error(s);}
int main(int argc,char **argv)
{
 try {
  check(argc==3,"usage: H2gChecks retained-H2c-preflight.json output.json");
  boost::property_tree::ptree retained;boost::property_tree::read_json(argv[1],retained);
  const double us=.07426793792419734,nu=3.8095238095238094e-6,H=.0065,z0=.00012222264395414748;
  const auto cls=RANSNeutralReference::CLS();const double kap=cls.Kappa();
  RANSParameters c;c.sigmaEpsilon=kap*kap/((c.cEpsilon2-c.cEpsilon1)*std::sqrt(c.cMu));c.yapCoefficient=0;
  const auto stdm=RANSNeutralReference::Standard(c);const double U=us/kap*std::log(H/z0);
  check(std::abs(stdm.Kappa()-kap)<2e-16,"standard equilibrium kappa");
  std::vector<std::string> models;std::vector<std::shared_ptr<EquilibriumAWFProfiles>> profiles;
  double rootResidual=0;
  for(const auto &m:{cls,stdm}) {
   const bool standard=m.Model()==RANSModel::StandardKEpsilon;const double K=m.Q()*us*us;auto wall=CalibratedAWFParameters(m);
   const double initialLow=H*1e-8,initialHigh=H*(1-1e-9);double lo=initialLow,hi=initialHigh;
   auto f=[&](double h){return us*us*AnalyticalWallFunction({nu,H,h,H/2},K,wall).Integrals().j0-U;};const double fl=f(lo),fh=f(hi);check(fl*fh<0,"roughness root bracket");
   for(int i=0;i<100;++i){double mid=(lo+hi)/2;if(mid==lo||mid==hi)break;if(fl*f(mid)>0)lo=mid;else hi=mid;}
   const double h=(lo+hi)/2;rootResidual=std::max(rootResidual,std::abs(f(h)));
   if(!standard)check(std::abs(h-.0022472751006336437)<1e-15,"retained CLS roughness differs");
   AWFProfileOptions o;o.formulation=standard?AWFProfileFormulation::StandardCalibratedBudgetAWF:AWFProfileFormulation::CLSCalibratedBudgetAWF;o.quadratureOrder=32;
   auto p=std::make_shared<EquilibriumAWFProfiles>(AWFProfileGeometry{nu,H,h,H,H},K,m,AWFNormalMetric{},o,wall);profiles.push_back(p);
   const auto &a=p->Original();const auto j=a.Integrals();const double budget=std::pow(us,4)*(j.j0-nu*j.j2),A=std::pow(K,1.5)/wall.cLength,cut=H*std::exp(1-budget/A);
   check(std::abs(cut-p->Reference().cutoff)<1e-13,"independent cutoff root differs");
   double P=0,E=0,minK=1e100;for(const auto &q:p->Quadrature()){P+=q.weight*p->ReferenceProduction(q.r);E+=q.weight*p->Dissipation(q.r);minK=std::min(minK,p->ReferenceK(q.r));}
   check(std::abs(P-E)<1e-12 && minK>0,"reference budget or positivity");
   double neutral=0,work=0,epsBalance=0,inactive=0;
   for(double z:{H,.04,.1,.8}){
    const auto v=m.AtHeight(z,us,H,U);RANSTensor g{};g[0][2]=v.velocityGradient;
    const auto b=EvaluateRANSClosure(m.Model(),v.k,v.epsilon,g,z,m.Parameters());
    neutral=std::max(neutral,std::abs(b.eddyViscosity-v.eddyViscosity));neutral=std::max(neutral,std::abs(b.production-v.epsilon));
    long double contraction=0;for(int i=0;i<3;++i)for(int j=0;j<3;++j){contraction-=b.reynoldsStress[i][j]*g[i][j];if(standard)inactive=std::max({inactive,std::abs(b.quadraticStress[i][j]),std::abs(b.cubicStress[i][j]),std::abs(b.yap)});}
    work=std::max(work,std::abs(double(contraction)-b.production));
    // d/dz[(nu_t/sigma_e) epsilon'] plus the complete source, including CLS Yap.
    epsBalance=std::max(epsBalance,std::abs(us*us*us*us/(m.Parameters().sigmaEpsilon*z*z)+b.epsilonSource));
   }
   check(neutral<1e-12 && work<1e-12 && epsBalance<1e-10 && inactive==0,"neutral closure/source balance");
   std::vector<std::string> curved;for(double curvature:{-19.73920880217871,0.,19.73920880217871})for(double ratio:{.8,1.,1.2}){
    EquilibriumAWFProfiles v({nu,H,h,H,H},ratio*K,m,{curvature,0,.1},o,wall);curved.push_back(J::Numbers(std::vector<double>{curvature,ratio,v.Reference().cutoff,v.Reference().budgetResidual,v.ReferenceK(0)}));}
   models.push_back(J::Object({{"name",J::Quote(standard?"StandardKEpsilon":"CLSHighRe")},{"q",J::Number(m.Q())},{"kappa",J::Number(m.Kappa())},{"K0",J::Number(K)},{"sigmaEpsilon",J::Number(m.Parameters().sigmaEpsilon)},
    {"cMu_wall",J::Number(wall.cMu)},{"cLength_wall",J::Number(wall.cLength)},{"h_s",J::Number(h)},{"root_bracket",J::Numbers(std::vector<double>{initialLow,initialHigh,fl,fh})},{"root_residual",J::Number(f(h))},{"viscosity_intercept",J::Number(a.ViscosityInterceptHeight())},
    {"cutoff",J::Number(p->Reference().cutoff)},{"budget_P_E",J::Numbers(std::vector<double>{P,E})},{"P_budget_independent",J::Number(budget)},{"A_epsilon",J::Number(A)},{"independent_cutoff",J::Number(cut)},{"E_C_H",J::Number(p->Dissipation(H))},{"j0_j1_j2",J::Numbers(std::vector<double>{j.j0,j.j1,j.j2})},{"J1_over_J0",J::Number(j.j1/j.j0)},
    {"min_inner_reference_k",J::Number(minK)},{"k_ground",J::Number(p->ReferenceK(0))},{"neutral_error",J::Number(neutral)},{"stress_work_error",J::Number(work)},{"epsilon_balance_error",J::Number(epsBalance)},{"inactive_terms_max",J::Number(inactive)},{"curved_K_screen",J::List(curved)}}));
  }
  double uerr=0,nterr=0,eerr=0,kgerr=0,koffset=0;const double dK=(stdm.Q()-cls.Q())*us*us;
  for(int i=0;i<=160;++i){double z=i==160?H:H*i/160;const auto &a=*profiles[0],&b=*profiles[1];uerr=std::max(uerr,std::abs(a.Original().Velocity(z,{U,0},{0,0})[0]-b.Original().Velocity(z,{U,0},{0,0})[0]));nterr=std::max(nterr,std::abs(a.Original().EddyViscosity(z)-b.Original().EddyViscosity(z)));eerr=std::max(eerr,std::abs(a.Dissipation(z)-b.Dissipation(z)));kgerr=std::max(kgerr,std::abs(a.ReferenceKGradient(z)-b.ReferenceKGradient(z)));koffset=std::max(koffset,std::abs((b.ReferenceK(z)-a.ReferenceK(z))-dK));}
  check(uerr<1e-11 && nterr<1e-14 && eerr<1e-10 && kgerr<1e-8 && koffset<1e-10,"paired inner flat identities");
  std::ofstream out(argv[2]);out<<J::Object({{"version",J::Quote("H2gNativeCalibrationV1")},{"passed","true"},{"u_tau",J::Number(us)},{"nu",J::Number(nu)},{"H",J::Number(H)},{"z0_eff",J::Number(z0)},{"kappa",J::Number(kap)},{"U_H",J::Number(U)},{"models",J::List(models)},{"inner_identity_errors_U_nu_epsilon_kgradient_koffset",J::Numbers(std::vector<double>{uerr,nterr,eerr,kgerr,koffset})},{"delta_K0",J::Number(dK)},{"governing_solves","0"}})<<std::endl;
 }catch(const std::exception &e){std::cerr<<e.what()<<std::endl;return 1;}
}
