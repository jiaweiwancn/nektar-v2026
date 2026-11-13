// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_MODEL_SPEC_H
#define NEKTAR_RANS_WALL_MODEL_SPEC_H
#include <IncNavierStokesSolver/EquationSystems/RANSWallPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <cmath>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <boost/property_tree/json_parser.hpp>
namespace Nektar
{
/** One effective native model, shared by evolution and read-only probes.
 * The existing CLS default and signature remain unchanged. Accepted standard
 * cases opt in to this calibration; an ordinary default standard is distinct. */
inline void LoadRANSWallModel(const LibUtilities::SessionReaderSharedPtr &s,
                             RANSWallPatchParameters &p)
{
    using namespace SolverUtils;
    using RANSWallAlgebra::Require;
    std::string name,calibration,wall;
    s->LoadSolverInfo("TurbulenceModel",name,"CLSHighRe");
    s->LoadSolverInfo("RANSWallModelCalibration",calibration,"");
    s->LoadSolverInfo("RANSWallFormulation",wall,"EPAWF");
    Require(name=="CLSHighRe" || name=="StandardKEpsilon","unknown canonical turbulence model");
    RANSParameters c;
    if(name=="StandardKEpsilon" && calibration=="StandardNeutralCalibratedV1")
    {
        c.sigmaEpsilon=std::pow(RANSNeutralReference::CLS().Kappa(),2)/((c.cEpsilon2-c.cEpsilon1)*std::sqrt(c.cMu));
        c.yapCoefficient=0;
        p.modelCalibration=calibration;
    }
    else Require(calibration.empty(),"unsupported turbulence calibration/model pair");
    for(const auto &v:std::vector<std::pair<std::string,double>>{
        {"RANSCmu",c.cMu},{"RANSCepsilon1",c.cEpsilon1},{"RANSCepsilon2",c.cEpsilon2},
        {"RANSSigmaK",c.sigmaK},{"RANSSigmaEpsilon",c.sigmaEpsilon},{"RANSYapCoefficient",c.yapCoefficient}})
        if(s->DefinesParameter(v.first))Require(std::abs(s->GetParameter(v.first)-v.second)<=1e-13*std::max(1.,std::abs(v.second)),"inconsistent effective calibrated coefficient: "+v.first);
    p.model=name=="CLSHighRe"?RANSNeutralReference::CLS(c):RANSNeutralReference::Standard(c);
    Require(wall!="StandardCalibratedBudgetAWF" || p.modelCalibration=="StandardNeutralCalibratedV1","standard AWF needs explicit calibrated model");
    Require(wall!="CLSCalibratedBudgetAWF" || name=="CLSHighRe","CLS AWF cannot be relabelled standard");
    if(calibration=="StandardNeutralCalibratedV1")
    {
        s->LoadSolverInfo("RANSWallModelManifestSHA256",p.modelManifest,"");
        Require(p.modelManifest.size()==64,"standard case needs model/calibration manifest identity");
    }
    double us;s->LoadParameter("RANSWallFrictionVelocity",us,.4);
    if(calibration=="StandardNeutralCalibratedV1")
    {
        std::string file;s->LoadSolverInfo("RANSWallModelManifestFile",file,"");
        Require(!file.empty() && RANSWallR3::HashFile(file)==p.modelManifest,"model/calibration manifest content mismatch");
        boost::property_tree::ptree m;boost::property_tree::read_json(file,m);
        Require(m.get<std::string>("version")==calibration && m.get<std::string>("model")==name,"model/calibration manifest interpretation mismatch");
        const auto awf=CalibratedAWFParameters(p.model);
        for(const auto &v:std::vector<std::pair<std::string,double>>{
            {"coefficients.Cmu",c.cMu},{"coefficients.Cepsilon1",c.cEpsilon1},{"coefficients.Cepsilon2",c.cEpsilon2},
            {"coefficients.sigmaK",c.sigmaK},{"coefficients.sigmaEpsilon",c.sigmaEpsilon},{"coefficients.YapCoefficient",0},
            {"q",p.model.Q()},{"kappa",p.model.Kappa()},{"K0",p.model.Q()*us*us},{"u_tau",us},
            {"wall.cMu",awf.cMu},{"wall.cLength",awf.cLength},{"physical_terms.Gamma",0},{"physical_terms.gamma0",0}})
            Require(std::abs(m.get<double>(v.first)-v.second)<=1e-12*std::max(1.,std::abs(v.second)),"native/model manifest coefficient mismatch: "+v.first);
        double H,rough,visc,z0;s->LoadParameter("RANSWallHeight",H,0);s->LoadParameter("RANSWallRoughness",rough,0);s->LoadParameter("Kinvis",visc,0);s->LoadParameter("RANSWallNeutralZ0",z0,0);
        Require(m.get<double>("H")==H && m.get<double>("wall.h_s")==rough && m.get<double>("nu")==visc && m.get<double>("z0_eff")==z0,"model manifest physical calibration mismatch");
    }

    for(const auto &v:std::vector<std::pair<std::string,double>>{{"Qstar",p.model.Q()},{"Kappa",p.model.Kappa()},{"K0",p.model.Q()*us*us}})
        if(s->DefinesParameter(v.first))Require(std::abs(s->GetParameter(v.first)-v.second)<=1e-12*std::max(1.,std::abs(v.second)),"native coefficient/reference expression mismatch: "+v.first);
}
inline std::string RANSWallModelJSON(const RANSWallPatchParameters &p,double us,
                                   const std::string &variant,double beta,double gamma0,double scale)
{
    namespace J=RANSWallSolveJson;using namespace SolverUtils;
    const auto &c=p.model.Parameters();const bool cls=p.model.Model()==RANSModel::CLSHighRe;
    return J::Object({{"version",J::Quote("CanonicalWallModelV1")},{"name",J::Quote(cls?"CLSHighRe":"StandardKEpsilon")},
        {"calibration",J::Quote(p.modelCalibration.empty()?(cls?"DeployedCLSV1":"DefaultStandardV1"):p.modelCalibration)},
        {"manifest_sha256",J::Quote(p.modelManifest)},{"coefficients",J::Object({{"Cmu",J::Number(c.cMu)},{"Cepsilon1",J::Number(c.cEpsilon1)},{"Cepsilon2",J::Number(c.cEpsilon2)},{"sigmaK",J::Number(c.sigmaK)},{"sigmaEpsilon",J::Number(c.sigmaEpsilon)},{"yapCoefficient",J::Number(c.yapCoefficient)},{"yapLength",J::Number(c.yapLength)}})},
        {"q",J::Number(p.model.Q())},{"kappa",J::Number(p.model.Kappa())},{"K0",J::Number(p.model.Q()*us*us)},
        {"quadratic_active",J::Boolean(cls)},{"cubic_active",J::Boolean(cls && beta!=0)},{"Yap_active",J::Boolean(cls && c.yapCoefficient>0)},
        {"variant",J::Quote(variant)},{"beta",J::Number(beta)},{"gamma0",J::Number(gamma0)},{"Gamma",J::Number(scale)},
        {"physical_regularization_active",J::Boolean(cls && (gamma0>0 || scale>0))}});
}
}
#endif
