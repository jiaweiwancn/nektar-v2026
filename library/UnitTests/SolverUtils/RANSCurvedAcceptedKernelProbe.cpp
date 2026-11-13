// SPDX-License-Identifier: MIT
#include <SolverUtils/Turbulence/EquilibriumAWFProfiles.h>
#include <SolverUtils/Turbulence/RANSWallNormalMap.h>
#include <boost/property_tree/json_parser.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
using namespace Nektar::SolverUtils;
int main(int argc,char **argv)
{
    try
    {
        if(argc!=3)throw std::invalid_argument("usage: kernel input.json output.json");
        boost::property_tree::ptree p;boost::property_tree::read_json(argv[1],p);
        const double nu=p.get<double>("nu",1.5e-5),H=p.get<double>("H",.5),hr=p.get<double>("hr",.05),K=p.get<double>("K");
        const auto m=RANSResolvedModel::Resolve(ParseRANSModel(p.get<std::string>("model")),RANSKappaPolicy::Natural,.41,RANSWallFormulation::AcceptedBudgetAWF,nu);
        AWFProfileOptions o;o.formulation=AWFProfileFormulation::AcceptedBudgetAWF;o.quadratureOrder=p.get<int>("order",64);
        if(p.get<bool>("verification",true))o.curvedCapability=AWFCurvedCapability::AcceptedVerificationV1;
        const AWFProfileGeometry g{nu,H,hr,p.get<double>("sample",H),p.get<double>("tail",H)};
        const AWFNormalMetric metric{p.get<double>("gamma",0),0};
        EquilibriumAWFProfiles f(g,K,m,metric,o,CalibratedAWFParameters(m));
        std::ofstream out(argv[2]);out<<std::setprecision(17)<<"{\"status\":\"pass\",\"cutoff\":"<<f.Reference().cutoff<<",\"budget\":"<<f.Reference().productionBudget<<",\"minimum_budget\":"<<f.Reference().minimumBudget<<",\"budget_residual\":"<<f.Reference().budgetResidual<<",\"samples\":[";
        int i=0;for(const auto &r:p.get_child("r")){const double z=r.second.get_value<double>();if(i++)out<<',';out<<'['<<z<<','<<metric.Jacobian(z)<<','<<f.ReferenceProduction(z)<<','<<f.Dissipation(z)<<','<<f.ReferenceK(z)<<','<<f.ReferenceKGradient(z)<<','<<f.ReferenceKFlux(z)<<','<<f.ReferenceKDiffusivity(z)<<']';}
        out<<"]}\n";return 0;
    }
    catch(const std::exception &e){if(argc==3){std::ofstream out(argv[2]);out<<"{\"status\":\"rejected\",\"reason\":"<<std::quoted(e.what())<<"}\n";}std::cerr<<e.what()<<'\n';return 1;}
}
