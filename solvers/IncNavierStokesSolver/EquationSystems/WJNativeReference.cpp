// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeReference.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <boost/property_tree/json_parser.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iterator>
#include <limits>
#include <stdexcept>

namespace Nektar::WJNative
{
namespace
{
using Tree = boost::property_tree::ptree;
void Require(bool condition, const std::string &message)
{
    if (!condition) throw std::invalid_argument("H27 reference: "+message);
}
void Finite(double v, const std::string &name)
{ Require(std::isfinite(v), "nonfinite "+name); }
double Number(const Tree &tree, const std::string &path)
{ double v=tree.get<double>(path); Finite(v,path); return v; }
void Same(double a, double b, const std::string &name)
{
    Finite(a,name); Finite(b,name);
    Require(std::abs(a-b)<=1e-13*std::max(std::abs(a),std::abs(b)),name+" differs");
}
std::vector<double> Numbers(const Tree &tree, const std::string &path)
{
    std::vector<double> values;
    for (const auto &entry:tree.get_child(path))
    {
        Require(entry.first.empty(),path+" is not an array");
        double value=entry.second.get_value<double>(); Finite(value,path);
        values.push_back(value);
    }
    Require(!values.empty(),path+" is empty"); return values;
}
void HashSyntax(const std::string &hash)
{
    Require(hash.size()==64 && hash.find_first_not_of("0123456789abcdef")==std::string::npos,
            "expected SHA256 identity is malformed");
}
long double HistoricalQ(long double r,long double A,long double cutoff,
                       long double d,long double a,long double nu)
{
    const long double z=a*r/d;
    const long double sink=r<=cutoff?r/cutoff:1+std::log(r/cutoff);
    // Algebraically identical to log1p(z)+nu/(d+a*r)-nu/d, but avoids
    // subtracting nearly equal molecular terms at the ground.
    const long double production=std::log1p(z)-(nu/d)*z/(1+z);
    return A*(sink-production);
}
}

double Reference::Polynomial::Evaluate(double argument) const
{
    Finite(argument,"polynomial argument");
    Require(argument>=lower && argument<=upper,"polynomial argument outside qualified interval");
    const double x=(2*argument-lower-upper)/(upper-lower);
    double b1=0,b2=0;
    for (std::size_t j=coefficients.size()-1;j>0;--j)
    {
        const double b=2*x*b1-b2+coefficients[j]; b2=b1; b1=b;
    }
    const double value=x*b1-b2+coefficients[0]; Finite(value,"Clenshaw result"); return value;
}

Reference::Polynomial Reference::Polynomial::Differentiate() const
{
    const std::size_t n=coefficients.size()-1;
    Require(n>=1,"cannot differentiate constant inlet panel");
    std::vector<double> d(n+2,0.);
    for (std::size_t j=n;j>0;--j) d[j-1]=d[j+1]+2*j*coefficients[j];
    d[0]*=.5;
    d.resize(n);
    for (double &v:d) v*=2/(upper-lower);
    return {lower,upper,std::move(d)};
}

ReferencePolicy Reference::ParsePolicy(const std::string &name)
{
    if (name=="Historical") return ReferencePolicy::Historical;
    if (name=="ExactInnerFiniteNuV1") return ReferencePolicy::ExactInnerFiniteNuV1;
    throw std::invalid_argument("H27 reference: unknown WJReferencePolicy "+name);
}

const char *Reference::PolicyName() const
{
    return m_policy==ReferencePolicy::ExactInnerFiniteNuV1?"ExactInnerFiniteNuV1":"Historical";
}

Reference::Reference(const std::string &json, ReferencePolicy policy)
    : Reference(json, policy==ReferencePolicy::ExactInnerFiniteNuV1 ?
                QualifiedExactInnerSha256 : QualifiedTableSha256, InnerKTransport::ScalarV2, policy)
{
}

Reference::Reference(const std::string &json, const std::string &expectedSha256,
                     InnerKTransport transport)
    : Reference(json,expectedSha256,transport,ReferencePolicy::Historical)
{
}

Reference::Reference(const std::string &json, const std::string &expectedSha256,
                     InnerKTransport transport, ReferencePolicy policy)
    : m_policy(policy)
{
    const bool exact=policy==ReferencePolicy::ExactInnerFiniteNuV1;
    Require(exact || policy==ReferencePolicy::Historical,"unsupported reference policy");
    Require(!exact || (expectedSha256==QualifiedExactInnerSha256 &&
                      transport==InnerKTransport::ScalarV2),"exact-inner identity or transport differs");
    HashSyntax(expectedSha256);
    std::ifstream input(json,std::ios::binary);
    Require(input.good(),"cannot read table "+json);
    const std::string bytes((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());
    Require(!input.bad() && !bytes.empty() && bytes.size()<1024*1024,"invalid table size/read");
    m_identity=RANSWallR3::Hash(bytes);
    Require(m_identity==expectedSha256,"immutable table SHA256 differs");
    m_boundaryIdentity=m_identity;
    Tree table; std::istringstream stream(bytes); boost::property_tree::read_json(stream,table);
    Require(table.get<std::string>("schema")==
        (exact?"WJH27ExactInnerFiniteNuTablesV1":"WJH27V2ReferenceTablesV1"),"schema differs");
    Require(table.get<std::string>("status")=="PASS","table is not qualified PASS");
    Require(table.get<std::string>("model")=="IsotropicLowerWallFullWJCollarV2","model identity differs");
    Require(table.get<std::string>("coefficient_convention")==
        "Ascending Chebyshev T_j coefficients; x=(2*argument-lo-hi)/(hi-lo)","coefficient convention differs");
    Require(table.get<std::string>("parameters.inner_diffusion")=="scalar","inner diffusion differs");
    if (exact)
    {
        Require(table.get<std::string>("policy")==PolicyName() &&
                table.get<std::string>("transport")==TransportName(),"exact-inner policy/transport binding differs");
        // The pinned byte identity authenticates the component contents as
        // well as these independently generated canonical component digests.
        for (const auto &[path, hash]:std::array<std::array<const char *,2>,7>{{
            {{"bindings.calibration","e3ade3777175e06286e1391bb7d4436dba305a46068fa27751f539eec9bfcb60"}},
            {{"bindings.cutoff_table","97101fe763cd4b0dc6bb2fffabc3d24c53b85fc34733b97efcd2755fc31cbc41"}},
            {{"bindings.inlet","e13609d631071d1372931444ccd110729a4789b339ac06634a6d8e969de33d6a"}},
            {{"bindings.loads","8a849654f638a1942cf6dab04f51c0f2f95223053bbfdb8b96f64bb43bae6f0d"}},
            {{"bindings.geometry","67bcb96b0d8696aeda4286aa79b6fc193e7239f8e9e131473bfcba5dd34c4dfe"}},
            {{"bindings.units","b86156fa2c62f755705607bcbc1d1de949373b42d212609d2b97720007466805"}},
            {{"bindings.sources","ca5457f7f440982afcb13eba504d54a725de55da386f2584e3bee0af9fd5080e"}}
        }})
            Require(table.get<std::string>(path)==hash,std::string(path)+" differs");
    }
    else
    {
        Require(table.get<std::string>("sources.generator.sha256")==
            "65664ef400cfd90f5bb4d83e24b08e0f9ac52d1f6f4385b87012360accd526e0","generator identity differs");
        Require(table.get<std::string>("sources.reference_oracle.sha256")==
            "f5c97e5c1f5426df4ac363b8d8dd99b84d636514a0cb899dec1219ceda4b35e9","oracle identity differs");
        Require(table.get<std::string>("sources.amended_model_contract.sha256")==
            "727b59a125429f8ecd1cb865be627de59b489cdd20773589a74961bdd9ae1922","amended model identity differs");
        Require(table.get<std::string>("inlet_binding.calibration.sha256")==
            "00007be74f8689b09d6d295b21e458041918e21e2e5a7a081e2a8d0ccc75c454","historical calibration identity differs");
    }

    auto &p=m_p;
    p.KIn=Number(table,"parameters.K"); p.H=Number(table,"parameters.H");
    p.roughness=Number(table,"parameters.h"); p.nu=Number(table,"parameters.nu");
    p.Cd=Number(table,"parameters.cd"); p.sigmaK=Number(table,"parameters.sigma_k");
    p.sigmaE=Number(table,"guard.sigma_e"); p.qe=Number(table,"constants.q_e");
    p.Ne=Number(table,"constants.N_e"); p.kappa=Number(table,"constants.kappa");
    p.cLength=Number(table,"constants.c_length"); p.d0=Number(table,"inlet_binding.d0");
    p.topRadius=Number(table,"inlet_binding.top_model_radius"); p.uTau=Number(table,"inlet_binding.u_tau");
    Same(p.KIn,.008809015145515192,"fixed inlet K"); Same(p.H,.15,"fixed H");
    Same(p.roughness,exact?.11541061001735275:.10854503542535461,"fixed roughness"); Same(p.nu,6.423173803526447e-5,"physical nu");
    Same(p.Cd,572./6561.,"independent Cd"); Same(p.sigmaK,1.,"sigma_k"); Same(p.sigmaE,1.3,"sigma_e");
    Same(p.qe,1/std::sqrt(572./6561.),"q_e"); Same(p.Ne,81./20.,"N_e");
    Same(p.kappa,std::sqrt(1.3*(1.82-1.44)*std::sqrt(572./6561.)),"kappa");
    Same(p.cLength,p.kappa*std::pow(p.qe,1.5),"c_length"); Same(p.d0,.0875,"physical displacement");
    Same(p.topRadius,20.-p.d0,"top model radius"); Same(p.uTau,.051,"fixed friction velocity");
    Same(Number(table,"parameters.endpoint_flux_over_A"),0.,"reference endpoint flux");
    Same(Number(table,"guard.reference_endpoint_flux_over_A"),0.,"reference endpoint flux guard");
    Same(Number(table,"inlet_binding.K_in"),p.KIn,"bound inlet K");
    Same(Number(table,"inlet_binding.roughness"),p.roughness,"bound roughness");
    const auto interval=Numbers(table,"K_interval"), factors=Numbers(table,"K_factor_interval");
    Require(interval.size()==2 && factors.size()==2,"K interval shape");
    Same(factors[0],.9,"lower K factor"); Same(factors[1],1.1,"upper K factor");
    p.KMinimum=interval[0]; p.KMaximum=interval[1];
    Same(p.KMinimum,.9*p.KIn,"lower K bound"); Same(p.KMaximum,1.1*p.KIn,"upper K bound");

    auto polynomial=[&](const Tree &node,const std::string &values)
    {
        const auto domain=Numbers(node,"domain"); Require(domain.size()==2 && domain[1]>domain[0],"polynomial domain");
        Polynomial result{domain[0],domain[1],Numbers(node,values)};
        Require(result.coefficients.size()==node.get<std::size_t>("degree")+1,"coefficient count differs");
        return result;
    };
    const auto &cutoff=table.get_child("cutoff_table");
    Require(cutoff.get<std::string>("argument")=="K" && cutoff.get<std::string>("value")=="log(c/H)","cutoff argument/value differs");
    m_logCutoff=polynomial(cutoff,"coefficients");
    Same(m_logCutoff.lower,p.KMinimum,"cutoff lower domain"); Same(m_logCutoff.upper,p.KMaximum,"cutoff upper domain");
    m_logCutoffDerivative={p.KMinimum,p.KMaximum,Numbers(cutoff,"logcutoff_dK_coefficients")};
    const auto &c=m_logCutoff.coefficients; const std::size_t n=c.size()-1;
    Require(n>=1 && m_logCutoffDerivative.coefficients.size()==n,"cutoff derivative degree differs");
    std::vector<double> derivative(n+2,0.);
    for (std::size_t j=n;j>0;--j) derivative[j-1]=derivative[j+1]+2*j*c[j];
    derivative[0]*=.5;
    for (std::size_t j=0;j<n;++j)
        Require(std::abs(derivative[j]*2/(p.KMaximum-p.KMinimum)-m_logCutoffDerivative.coefficients[j])
                    <=1e-10*(1+std::abs(m_logCutoffDerivative.coefficients[j])),"analytic derivative coefficients differ");
    // |T_j|<=1 on the complete admitted interval. Enclose the derivative
    // derived from the value polynomial itself, including arithmetic reserve.
    const long double scale=2.L/(p.KMaximum-p.KMinimum),constant=derivative[0]*scale;
    long double tail=0;
    for(std::size_t j=1;j<n;++j)tail+=std::abs(derivative[j]*scale);
    const long double reserve=32*(n+1)*std::numeric_limits<double>::epsilon()*(std::abs(constant)+tail);
    m_logCutoffSlopeBounds={std::nextafter(double(constant-tail-reserve),-std::numeric_limits<double>::infinity()),
        std::nextafter(double(constant+tail+reserve),std::numeric_limits<double>::infinity())};
    Require(std::isfinite(m_logCutoffSlopeBounds[0])&&m_logCutoffSlopeBounds[1]<0,
        "cutoff table has no certified strictly decreasing inverse");
    Require(Number(cutoff,"maximum_budget_or_ground_flux_error")<=1e-10 &&
            Number(cutoff,"maximum_K_in_times_derivative_error")<=1e-8,"cutoff qualification tolerance failed");
    Require(cutoff.get<int>("heldout_K_count")>=33,"missing complete-IVP validation coverage");
    for (const auto &entry:cutoff.get_child("validations"))
    {
        const auto &v=entry.second;
        Require(std::abs(Number(v,"ground_flux_over_A"))<=1e-10 &&
                std::abs(Number(v,"independent_budget_over_A"))<=1e-10 &&
                Number(v,"min_k_over_K")>0 && Number(v,"minimum_shear")>0 &&
                Number(v,"minimum_shear_tangent")>0,"reference validation witness is inadmissible");
    }
    p.tau0=p.KIn/p.qe; Same(p.tau0,p.uTau*p.uTau,"fixed total traction");
    if (exact)
    {
        p.z0=Number(table,"inlet_binding.z0e");
        Same(p.z0,.005676270321216798,"fixed outer aerodynamic roughness");
        Same(Number(table,"geometry.H"),p.H,"matching radius");
        Same(Number(table,"geometry.d0"),p.d0,"geometry displacement");
        Same(Number(table,"geometry.rT"),p.topRadius,"geometry top");
        Same(Number(table,"geometry.force_support"),p.roughness,"force support");
        Same(Number(table,"units.L0_m"),.04,"reference length");
        Same(Number(table,"units.H_m"),p.H*.04,"physical matching height");
        Same(Number(table,"units.interface_elevation_m"),(p.H+p.d0)*.04,"physical matching elevation");
        const auto &metrics=cutoff.get_child("maxima");
        Require(Number(metrics,"relative_cutoff_error")<=1e-10 &&
                Number(metrics,"K_times_log_derivative_error")<=1e-8 &&
                Number(metrics,"epsilon_K_derivative_fixed_radius_error")<=1e-8 &&
                Number(metrics,"inverse_K_error_over_Kin")<=1e-9 &&
                Number(metrics,"represented_inverse_residual")<=256*std::numeric_limits<double>::epsilon()*p.H,
                "exact-inner cutoff controls failed");
        m_inletA=Amplitude(p.KIn);
        Same(m_inletA,Number(table,"inlet_binding.A_in"),"outer inlet amplitude");
        Same(m_inletA/p.H,Number(table,"inlet_binding.epsilon_in_H"),"outer epsilon inlet/H corner");
        m_inletCutoff=Cutoff(p.KIn);
        Same(m_inletCutoff,Number(table,"inlet.cutoff"),"prescribed inlet cutoff");
        const auto edges=Numbers(table,"geometry.material_edges");
        Require(edges.size()==6,"exact-inner material edge count");
        const std::array<double,6> required{0.,m_inletCutoff,.4*p.H,.75*p.H,p.roughness,p.H};
        Require(std::is_sorted(required.begin(),required.end()),"exact-inner material order differs");
        for (std::size_t j=0;j<edges.size();++j) Same(edges[j],required[j],"material edge");
        m_referenceShear=p.tau0+p.nu*p.qe*m_inletA/(p.H*p.KIn);
        Same(m_referenceShear,Number(table,"calibration.reference_shear"),"reference total shear");
        Same(m_referenceShear,Number(table,"loads.interface_r"),"interface total shear");
        Same(TopTangentialTraction(),Number(table,"loads.top_x"),"top tangent");
        Same(TopEpsilonOutwardFlux(),Number(table,"loads.top_e_out"),"top epsilon flux");
        Same(Number(table,"loads.interface_x"),-p.KIn/p.Ne,"interface axial traction");
        for (const auto &entry:table.get_child("inlet.scaled_errors"))
            Require(Number(table,"inlet.scaled_errors."+entry.first)<=1e-8,"inlet qualification failed");
        double end=0.;
        for (const auto &entry:table.get_child("inlet.panels"))
        {
            const auto &node=entry.second;
            const auto domain=Numbers(node,"domain");
            Require(domain.size()==2 && domain[0]==end && domain[1]>end && domain[1]<=p.H,
                    "exact-inner inlet panel partition differs");
            end=domain[1];
            auto read=[&](const char *field)
            {
                auto coefficients=Numbers(node,field);
                Require(coefficients.size()>=2 && coefficients.size()<=64,"invalid inlet panel degree");
                return Polynomial{domain[0],end,std::move(coefficients)};
            };
            InletPanel panel{read("U"),read("k"),read("Pi"),{},{}};
            panel.du=panel.u.Differentiate();
            panel.dk=panel.k.Differentiate();
            if (!m_exactInnerInlet.empty())
            {
                const auto &previous=m_exactInnerInlet.back();
                Require(std::abs(previous.u.Evaluate(domain[0])-panel.u.Evaluate(domain[0]))<=1e-10 &&
                        std::abs(previous.k.Evaluate(domain[0])-panel.k.Evaluate(domain[0]))<=1e-10*p.KIn,
                        "inlet panel field discontinuity");
            }
            m_exactInnerInlet.push_back(std::move(panel));
        }
        Require(end==p.H,"exact-inner inlet partition does not reach H");
        for (double edge:edges)
            Require(edge==0. || std::any_of(m_exactInnerInlet.begin(),m_exactInnerInlet.end(),
                [edge](const InletPanel &panel){return panel.u.upper==edge;}),
                "material edge is absent from inlet partition");
        Require(std::abs(InletU(0.))<=1e-10 &&
                std::abs(InletU(p.H)-p.uTau/p.kappa*std::log(p.H/p.z0))<=1e-10 &&
                std::abs(InnerPanel(p.H).k.Evaluate(p.H)-p.KIn)<=1e-10*p.KIn &&
                std::abs(InnerPanel(p.H).dk.Evaluate(p.H))*p.H<=1e-10*p.KIn &&
                std::abs(InletUGradient(p.H)-p.uTau/(p.kappa*p.H))*p.H<=1e-10 &&
                std::abs(InletModifiedPressure(p.H))<=1e-10*p.tau0,
                "exact-inner inlet endpoint continuity failed");
        for (double K:{p.KMinimum,p.KIn,p.KMaximum}) { Cutoff(K); Auxiliary(K,0.); }
        return;
    }
    const auto &outlet=table.get_child("outlet_traction");
    Same(Number(outlet,"lower_constant.tx"),p.KIn/p.Ne,"lower outlet tx");
    Same(Number(outlet,"lower_constant.tz"),p.tau0,"lower outlet tz");
    int index=0;
    for (const auto &entry:outlet.get_child("segments"))
    {
        Require(index<3,"too many outlet segments"); const auto &segment=entry.second;
        Require(segment.get<std::string>("argument")== (index==2?"log(r/H)":"r"),"outlet coordinate differs");
        m_outlet[index]=polynomial(segment,"tx_coefficients");
        Same(Number(segment,"tz_constant"),p.tau0,"outlet shear traction");
        Require(Number(segment,"maximum_tx_error_over_tau")<=1e-10,"outlet interpolation qualification failed");
        ++index;
    }
    Require(index==3,"missing outlet segments");
    Same(m_outlet[0].lower,.4*p.H,"outlet lower interface"); Same(m_outlet[0].upper,.75*p.H,"outlet collar start");
    Same(m_outlet[1].lower,.75*p.H,"outlet collar start"); Same(m_outlet[1].upper,p.H,"outlet H");
    Same(m_outlet[2].lower,0.,"outer log origin"); Same(m_outlet[2].upper,std::log(p.topRadius/p.H),"outer log top");
    for (const auto &entry:outlet.get_child("continuity"))
        Require(std::isfinite(entry.second.get_value<double>()) && entry.second.get_value<double>()<=1e-10,
                "outlet endpoint continuity failed");
    Require(std::abs(m_outlet[0].Evaluate(.4*p.H)-p.KIn/p.Ne)<=1e-10*p.tau0 &&
            std::abs(m_outlet[0].Evaluate(.75*p.H)-m_outlet[1].Evaluate(.75*p.H))<=1e-10*p.tau0 &&
            std::abs(m_outlet[1].Evaluate(p.H)-m_outlet[2].Evaluate(0.))<=1e-10*p.tau0,
            "native outlet endpoint continuity failed");

    m_inletA=Amplitude(p.KIn); m_inletSlope=p.kappa*p.uTau;
    m_inletD=p.nu+Aux(p.KIn,0.);
    const double z=m_inletSlope*p.H/m_inletD;
    const double budget=std::log1p(z)+p.nu/(m_inletD+m_inletSlope*p.H)-p.nu/m_inletD;
    // This OLD cutoff appears ONLY in the fixed historical inlet k profile.
    // The ordinary v2 dependent epsilon uses its qualified cutoff table.
    m_inletCutoff=p.H*std::exp(1-budget);
    Require(m_inletCutoff>0 && m_inletCutoff<p.H,"historical inlet cutoff is inadmissible");
    p.z0=p.H/(1+z); Same(p.z0,.005676270321216798,"fixed fitted aerodynamic roughness");
    Same(m_inletA,Number(table,"inlet_binding.A_in"),"inlet amplitude");
    Same(m_inletA/p.H,Number(table,"inlet_binding.epsilon_in_H"),"epsilon inlet/H corner");
    for (double K:{p.KMinimum,p.KIn,p.KMaximum}) { Cutoff(K); Auxiliary(K,0.); }
    Require(transport==InnerKTransport::ScalarV2||transport==InnerKTransport::AuxiliaryHistoricalFixtureV1,
            "unknown inner-k transport policy");
    if (transport==InnerKTransport::AuxiliaryHistoricalFixtureV1) BindFixtureTransport();
}

void Reference::CheckK(double K) const
{ Finite(K,"K"); Require(K>=m_p.KMinimum && K<=m_p.KMaximum,"K outside qualified [.9,1.1] inlet interval"); }
void Reference::CheckRadius(double r) const
{ Finite(r,"radius"); Require(r>=0 && r<=m_p.topRadius,"radius outside fixed physical domain"); }
double Reference::Cutoff(double K) const
{
    CheckK(K);
    const double logCutoff=m_transport==InnerKTransport::ScalarV2 ? m_logCutoff.Evaluate(K) : HistoricalLogCutoff(K)[0];
    const double c=m_p.H*std::exp(logCutoff);
    Require(std::isfinite(c) && c>0 && c<.4*m_p.H,"cutoff outside selected 0<c<.4H branch"); return c;
}
double Reference::CutoffDerivative(double K) const
{
    const double slope=m_transport==InnerKTransport::ScalarV2 ? m_logCutoffDerivative.Evaluate(K) : HistoricalLogCutoff(K)[1];
    const double v=Cutoff(K)*slope;
    Finite(v,"cutoff derivative");return v;
}

const char *Reference::TransportName() const
{
    return m_transport==InnerKTransport::ScalarV2 ? "ScalarCdV2" : "AuxiliaryHistoricalInnerKFixtureV1";
}

double Reference::KDiffusivity(double k, double epsilon, double auxiliary, bool outer) const
{
    if (!outer && m_transport==InnerKTransport::AuxiliaryHistoricalFixtureV1)
        return m_p.nu+auxiliary;
    return m_p.nu+m_p.Cd*k*k/(m_p.sigmaK*epsilon);
}

std::array<double,2> Reference::HistoricalLogCutoff(double K) const
{
    const auto ground=Auxiliary(K,0.);
    const double a=ground.dr, nu=m_p.nu, d=nu+ground.value, B=d+a*m_p.H;
    const double dK=ground.dK, BK=dK+a*m_p.H/(2*K);
    const double J=std::log1p(a*m_p.H/d)+nu/B-nu/d;
    const double JK=BK*(B-nu)/(B*B)-dK*(d-nu)/(d*d);
    Require(std::isfinite(J)&&std::isfinite(JK)&&JK>0,"invalid historical cutoff slope");
    return {1-J,-JK};
}

void Reference::BindFixtureTransport()
{
    Require(m_boundaryIdentity==QualifiedTableSha256,
            "the inner-k fixture requires the exact retained v2 boundary table");
    const long double nu=m_p.nu, H=m_p.H;
    const long double multiplier=10.7L*m_p.kappa/std::sqrt(static_cast<long double>(m_p.qe));
    auto roughness=[&](long double K){return m_p.roughness*std::sqrt(K)/(70*nu);};
    const long double X0=roughness(m_p.KMinimum), X1=roughness(m_p.KMaximum);
    auto m=[](long double X){return 1-.79L*std::pow(X,-.28L);};
    Require(multiplier>1 && X0>1 && m(X0)>std::max(0.L,.5L-.4L*std::pow(X0,.7L)) &&
            .28L*std::log(X1)<2, "historical roughness branch is not separated across the interval");
    auto bounds=[&](long double K,long double X){
        const long double a=m_p.kappa*std::sqrt(K/m_p.qe), logX=std::log(X), q=m(X)*logX;
        const long double d=nu*(1+multiplier*std::expm1(q));
        const long double factor=multiplier*std::exp(q)*nu/d;
        const long double dq=(m(X)+.2212L*std::pow(X,-.28L)*logX)/2;
        return std::array<long double,5>{a,d,d+a*H,factor,dq};
    };
    const auto lo=bounds(m_p.KMinimum,X0), hi=bounds(m_p.KMaximum,X1);
    // s = K d'(K)/d(K) is a product of a decreasing and an increasing
    // positive factor. Bound each factor by its endpoint. Set p=aH/B.
    // Then K J' = (1-nu/B) p (1/2-s) + (nu/d-nu/B) s.
    const long double s0=hi[3]*lo[4], s1=lo[3]*hi[4];
    Require(s0>0 && s1<.5L && hi[1]<lo[2], "historical cutoff monotonicity bound failed");
    const long double j0=(1-nu/lo[2])*(lo[0]*H/hi[2])*(.5L-s1)+(nu/hi[1]-nu/lo[2])*s0;
    const long double j1=(1-nu/hi[2])*(.5L-s0)+(nu/lo[1]-nu/hi[2])*s1;
    constexpr double reserve=1e-9;
    m_logCutoffSlopeBounds={double(-j1/m_p.KMinimum)-reserve, double(-j0/m_p.KMaximum)+reserve};
    Require(m_logCutoffSlopeBounds[1]<0 && m_logCutoffSlopeBounds[0]<=m_logCutoffSlopeBounds[1],
            "historical cutoff slope enclosure is not strictly decreasing");
    m_transport=InnerKTransport::AuxiliaryHistoricalFixtureV1;
    m_identity=RANSWallR3::Hash(std::string("WJH27|")+TransportName()+"|retained-v2-boundaries|"+m_boundaryIdentity);
    const auto range=CutoffRadiusRange();
    Require(0<range[0] && range[0]<range[1] && range[1]<.4*m_p.H,
            "historical cutoff leaves the admitted lower-wall domain");
}
std::array<double,2> Reference::CutoffRadiusRange() const
{ return {Cutoff(m_p.KMaximum),Cutoff(m_p.KMinimum)}; }
double Reference::KAtCutoff(double radius) const
{
    Finite(radius,"cutoff face radius");const auto range=CutoffRadiusRange();
    Require(radius>=range[0]&&radius<=range[1],"cutoff face radius outside qualified inverse interval");
    if(radius==range[0])return m_p.KMaximum;if(radius==range[1])return m_p.KMinimum;
    double lower=m_p.KMinimum,upper=m_p.KMaximum;
    for(int iteration=0;iteration<80;++iteration){
        const double middle=(lower+upper)/2;if(middle==lower||middle==upper)break;
        if(Cutoff(middle)>radius)lower=middle;else upper=middle;
    }
    const double K=(lower+upper)/2;
    Require(std::abs(Cutoff(K)-radius)<=256*std::numeric_limits<double>::epsilon()*m_p.H,
        "cutoff face inverse failed its current residual");return K;
}
double Reference::Amplitude(double K) const
{ CheckK(K); return K*std::sqrt(K)/m_p.cLength; }
double Reference::AmplitudeDerivative(double K) const
{ CheckK(K); return 1.5*std::sqrt(K)/m_p.cLength; }
double Reference::AmplitudeSecondDerivative(double K) const
{ CheckK(K); return .75/(std::sqrt(K)*m_p.cLength); }

ReferenceValue Reference::Auxiliary(double K,double r) const
{
    CheckK(K); CheckRadius(r);
    const double root=std::sqrt(K),x=m_p.roughness*root/(70*m_p.nu),logx=std::log(x);
    const double power1=std::pow(x,.7),power2=std::pow(x,-.28);
    const double m1=.5-.4*power1,m2=1-.79*power2;
    Require(m1!=m2,"roughness correlation derivative at branch intersection");
    const double exponent=m1>m2?m1:m2;
    const double exponentK=(m1>m2?-.28*power1:.2212*power2)/(2*K);
    const double b=-10.7*std::expm1(exponent*logx);
    const double bK=-10.7*std::exp(exponent*logx)*(exponentK*logx+exponent/(2*K));
    const double zv=m_p.nu*b/root,zvK=m_p.nu/root*(bK-b/(2*K));
    Require(std::isfinite(zv) && zv<0,"outside selected negative-intercept roughness branch");
    const double a=m_p.kappa*std::sqrt(K/m_p.qe);
    ReferenceValue value{a*(r-zv),a*(r-zv)/(2*K)-a*zvK,a};
    Require(std::isfinite(value.value) && value.value>0 && std::isfinite(value.dK),"invalid auxiliary viscosity");
    return value;
}

ReferenceValue Reference::DependentEpsilon(double K,double r) const
{
    CheckRadius(r); Require(r<=m_p.H,"dependent wall epsilon requested outside inner domain");
    const double A=Amplitude(K),Ap=AmplitudeDerivative(K);
    if (m_transport==InnerKTransport::AuxiliaryHistoricalFixtureV1)
    {
        const auto law=HistoricalLogCutoff(K);
        const double c=m_p.H*std::exp(law[0]);
        Require(c>0 && c<.4*m_p.H,"cutoff leaves the admitted lower-wall domain");
        if (r<=c) return {A/c,(Ap-A*law[1])/c,0.};
        return {A/r,Ap/r,-A/(r*r)};
    }
    const double c=Cutoff(K);
    // At the kink use the plateau-side derivative; MMS derivative evaluation
    // must avoid r=c. This prescribed max law is not a positivity clamp.
    if (r<=c) return {A/c,Ap/c-A*CutoffDerivative(K)/(c*c),0.};
    return {A/r,Ap/r,-A/(r*r)};
}

double Reference::InletQ(double r) const
{
    return double(HistoricalQ(r,m_inletA,m_inletCutoff,m_inletD,m_inletSlope,m_p.nu));
}
const Reference::InletPanel &Reference::InnerPanel(double r) const
{
    auto found=std::lower_bound(m_exactInnerInlet.begin(),m_exactInnerInlet.end(),r,
        [](const InletPanel &panel,double radius){return panel.u.upper<radius;});
    Require(found!=m_exactInnerInlet.end(),"exact-inner panel outside qualified radius");
    return *found;
}

double Reference::InletModifiedPressure(double r) const
{
    CheckRadius(r);
    Require(m_policy==ReferencePolicy::ExactInnerFiniteNuV1,"pressure panel requires exact-inner policy");
    return r<=m_p.H?InnerPanel(r).pressure.Evaluate(r):0.;
}

double Reference::InletU(double r) const
{
    CheckRadius(r);
    if (r<=m_p.H && m_policy==ReferencePolicy::ExactInnerFiniteNuV1)
        return InnerPanel(r).u.Evaluate(r);
    return r<=m_p.H?m_p.tau0/m_inletSlope*std::log1p(m_inletSlope*r/m_inletD)
                     :m_p.uTau/m_p.kappa*std::log(r/m_p.z0);
}
double Reference::InletK(double r) const
{
    CheckRadius(r); if (r>=m_p.H) return m_p.KIn;
    if (m_policy==ReferencePolicy::ExactInnerFiniteNuV1) return InnerPanel(r).k.Evaluate(r);
    using Real=long double;
    struct Panel { Real a,b,budget; };
    const Real budget=2e-12L*m_p.KIn;
    std::vector<Panel> pending;
    if(r<m_inletCutoff)
    {
        const Real leftBudget=budget*(Real(m_inletCutoff)-r)/(Real(m_p.H)-r);
        pending.push_back({r,m_inletCutoff,leftBudget});
        pending.push_back({m_inletCutoff,m_p.H,budget-leftBudget});
    }
    else pending.push_back({r,m_p.H,budget});
    Real correction=0,totalError=0;
    while(!pending.empty())
    {
        const Panel panel=pending.back();pending.pop_back();
        const Real half=(panel.b-panel.a)/2,mean=panel.a+half;
        // Installed Boost's adaptive routine scales the estimate/L1 but not
        // its returned local error. Own the mapping and subdivision here:
        // integrate on [-1,1], with the physical Jacobian INSIDE f and no
        // Boost recursion. Both estimates and errors now have k units.
        const auto mapped=[&](Real t)
        {
            const Real x=mean+half*t;
            return half*HistoricalQ(x,m_inletA,m_inletCutoff,m_inletD,m_inletSlope,m_p.nu)
                       /(Real(m_inletD)+Real(m_inletSlope)*x);
        };
        Real error=0,l1=0,coarseError=0;
        const Real value=boost::math::quadrature::gauss_kronrod<Real,61>::integrate(
            mapped,Real(-1),Real(1),0,Real(2e-13L),&error,&l1);
        const Real coarse=boost::math::quadrature::gauss_kronrod<Real,31>::integrate(
            mapped,Real(-1),Real(1),0,Real(2e-13L),&coarseError);
        error=std::max(error,std::abs(value-coarse));
        auto failure=[&](const char *reason)
        {
            std::ostringstream message;message<<std::setprecision(21)<<reason
                <<"; requested_r="<<r<<" panel_a="<<panel.a<<" panel_b="<<panel.b
                <<" value="<<value<<" error="<<error<<" L1="<<l1
                <<" GK31_error="<<coarseError<<" panel_budget="<<panel.budget;
            throw std::runtime_error(message.str());
        };
        if(!std::isfinite(value)||!std::isfinite(coarse)||!std::isfinite(error)||
           !std::isfinite(coarseError)||!std::isfinite(l1)||error<0||coarseError<0||l1<0)
            failure("nonfinite/invalid historical inlet k quadrature");
        if(error<=panel.budget){correction+=value;totalError+=error;continue;}
        if(!(mean>panel.a&&mean<panel.b))
            failure("historical inlet k quadrature exhausted representable interval progress");
        pending.push_back({mean,panel.b,panel.budget/2});
        pending.push_back({panel.a,mean,panel.budget/2});
    }
    Require(totalError<=budget,"historical inlet k quadrature exceeded global absolute error guard");
    const double k=double(Real(m_p.KIn)-correction);
    Require(std::isfinite(k) && k>0,"historical inlet k is nonpositive"); return k;
}
double Reference::InletE(double r) const
{ CheckRadius(r); Require(r>=m_p.H,"no transported inner inlet epsilon"); return m_inletA/r; }
double Reference::InletUGradient(double r) const
{
    CheckRadius(r);
    if (r<=m_p.H && m_policy==ReferencePolicy::ExactInnerFiniteNuV1)
        return InnerPanel(r).du.Evaluate(r);
    return r<=m_p.H?m_p.tau0/(m_inletD+m_inletSlope*r):m_p.uTau/(m_p.kappa*r);
}
double Reference::InletKGradient(double r) const
{
    CheckRadius(r);
    if (r>=m_p.H) return 0.;
    if (m_policy==ReferencePolicy::ExactInnerFiniteNuV1) return InnerPanel(r).dk.Evaluate(r);
    return InletQ(r)/(m_inletD+m_inletSlope*r);
}
double Reference::InletEGradient(double r) const
{ return -InletE(r)/r; }
std::array<double,2> Reference::OutletTraction(double r) const
{
    CheckRadius(r);
    if (m_policy==ReferencePolicy::ExactInnerFiniteNuV1)
    {
        const auto &p=m_p;
        const double tz=r<=p.H?m_referenceShear:p.tau0+p.nu*p.uTau/(p.kappa*r);
        return {-p.KIn/p.Ne-2*InletModifiedPressure(r),tz};
    }
    double tx=m_p.KIn/m_p.Ne;
    if (r>.4*m_p.H && r<=.75*m_p.H) tx=m_outlet[0].Evaluate(r);
    else if (r>.75*m_p.H && r<=m_p.H) tx=m_outlet[1].Evaluate(r);
    else if (r>m_p.H) tx=m_outlet[2].Evaluate(r==m_p.topRadius?m_outlet[2].upper:std::log(r/m_p.H));
    return {tx,m_p.tau0};
}
double Reference::TopTangentialTraction() const
{
    return OutletTraction(m_p.topRadius)[1];
}
double Reference::TopEpsilonOutwardFlux() const
{ return m_p.nu*m_inletA/(m_p.topRadius*m_p.topRadius)+m_p.Cd*m_p.KIn*m_p.KIn/(m_p.sigmaE*m_p.topRadius); }
}
