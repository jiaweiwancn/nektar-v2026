// SPDX-License-Identifier: MIT
#include <SolverUtils/Turbulence/RANSModels.h>
#include <SolverUtils/Turbulence/RANSTwoRegionHill.h>
#include <SolverUtils/Turbulence/RANSRobustVariant.h>
#include <SolverUtils/Turbulence/RANSWallInnerStress.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <boost/property_tree/json_parser.hpp>
#include <fstream>
#include <map>
#include <iostream>
using namespace Nektar::SolverUtils;
namespace J=Nektar::RANSWallSolveJson;
int main(int argc,char **argv)
{
 try{
  if(argc!=3)throw std::runtime_error("usage: H2gConstitutiveProbe input.json output.json");
  boost::property_tree::ptree root;boost::property_tree::read_json(argv[1],root);std::vector<std::string> output;
  RANSTwoRegionHill geometry(.04,.1,.0065,.8,-.2,.6);std::map<std::pair<double,double>,double> distances;
  const double us=root.get<double>("u_tau"),kap=root.get<double>("kappa");RANSParameters standard;standard.sigmaEpsilon=kap*kap/((standard.cEpsilon2-standard.cEpsilon1)*std::sqrt(standard.cMu));standard.yapCoefficient=0;
  for(const auto &entry:root.get_child("tuples"))
  {
   const auto &v=entry.second;const double k=v.get<double>("k"),e=v.get<double>("epsilon"),z=v.get<double>("physical_z"),nt=v.get<double>("inner_nu_t",0);const bool inner=v.get<bool>("inner",false);const int law=v.get<int>("law"),branch=v.get<int>("branch",0);RANSTensor g{};int i=0;for(auto &row:v.get_child("gradient")){int j=0;for(auto &x:row.second)g[i][j++]=x.second.get_value<double>();++i;}
   if(law<0 || law>5)throw std::runtime_error("unknown frozen law");const auto model=law==5?RANSModel::StandardKEpsilon:RANSModel::CLSHighRe;const auto par=law==5?standard:RANSParameters{};
   const double x=v.get<double>("physical_x",0);const auto key=std::make_pair(x,z);double wall=z;
   if(!v.get<bool>("flat_reference",false)){if(!distances.count(key))distances[key]=geometry.PhysicalWallDistance({x,0,z}).distance;wall=distances.at(key);}
   const auto c=branch==0?EvaluateRANSClosure(model,k,e,g,wall,par):EvaluateRANSClosureLinearization(model,k,e,g,wall,par,0,branch==1);RANSTensor Q=c.quadraticStress,C=c.cubicStress,R=c.reynoldsStress,reg{};
   for(int a=0;a<3;++a)for(int b=0;b<3;++b)
   {
    if(law==0 || law==1){R[a][b]-=.6*C[a][b];C[a][b]*=.4;}
    if(law==3 || law==4){R[a][b]-=Q[a][b]+C[a][b];Q[a][b]=C[a][b]=0;}
    if(law==0 || law==4){const double sr=(a==0&&b==2)||(a==2&&b==0)?us/(kap*z):0;reg[a][b]=-kap*us*z*(g[a][b]+g[b][a]-sr);R[a][b]+=reg[a][b];}
   }
   if(inner){R=EvaluateRANSWallInnerStress(k,nt,g,RANSWallInnerStressModel::Linear).reynoldsStress;Q={};C={};reg={};}
   long double P=0;for(int a=0;a<3;++a)for(int b=0;b<3;++b)P-=R[a][b]*g[a][b];
   output.push_back(J::Object({{"physical_wall_distance",J::Number(wall)},{"incoming_height",J::Number(z)},{"tensor",J::Rows(R)},{"quadratic",J::Rows(Q)},{"cubic",J::Rows(C)},{"regularization",J::Rows(reg)},{"Cmu",J::Number(c.cMu)},{"nu_t",J::Number(inner?nt:c.eddyViscosity)},{"production",J::Number(P)},{"epsilon_production",J::Number(par.cEpsilon1*e/k*P)},{"Yap",J::Number(c.yap)},{"strain",J::Number(c.strainInvariant)},{"rotation",J::Number(c.rotationInvariant)}}));
  }
  std::ofstream out(argv[2]);out<<J::Object({{"tuples",J::List(output)},{"governing_solves","0"}})<<std::endl;
 }catch(const std::exception &e){std::cerr<<e.what()<<std::endl;return 1;}
}
