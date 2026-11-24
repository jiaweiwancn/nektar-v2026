// SPDX-License-Identifier: MIT
#include "../h26/common.h"
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeAdmissibility.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeMixedFlow.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeReference.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeRules.h>

namespace J=Nektar::RANSWallSolveJson;

namespace
{
Tensor Rotate(const Tensor &rotation, const Tensor &value)
{
    Tensor result{};
    for (int i=0;i<3;++i) for(int j=0;j<3;++j)
        for (int k=0;k<3;++k) for(int l=0;l<3;++l)
            result[i][j]+=rotation[i][k]*value[k][l]*rotation[j][l];
    return result;
}

std::string RotationObservations()
{
    const double k=.008809015145515192,E=.0023155240547992535;
    const double nu=6.423173803526447e-5,H=.15,aux=.0015;
    Tensor gradient{};
    gradient[0][0]=.023;gradient[0][2]=.58;
    gradient[2][0]=-.12;gradient[2][2]=.05;
    std::vector<std::string> rows;
    for (double r : {.03,.078,H,.25}) for(double angle : {0.,.31,1.06,-1.7})
    {
        const double c=std::cos(angle),s=std::sin(angle);
        const Tensor Q{{{c,0,-s},{0,1,0},{s,0,c}}};
        const auto G=Rotate(Q,gradient);
        const Vector n{-s,0,c};
        const auto wall=AdaptAmended(k,E,G,n,nu,aux,r,H,ConstitutiveBackend::FastFull3DV1);
        const auto guard=CheckPlanarState(k,E,G,nu,aux,r,H,true,&wall,
                                         ConstitutiveBackend::FastFull3DV1,nullptr,n);
        Vector traction{};
        for (int i=0;i<3;++i) for(int j=0;j<3;++j)
            traction[i]+=(wall.V[i][j]-(i==j?.001:0.))*n[j];
        rows.push_back(J::Object({{"radius",J::Number(r)},{"angle",J::Number(angle)},
            {"R",J::Rows(wall.R)},{"V",J::Rows(wall.V)},{"traction",J::Numbers(traction)},
            {"production",J::Number(wall.production)},{"covariance",J::Number(guard.covarianceOverK)},
            {"principal",J::Number(guard.principalOverMu)},
            {"principal_raw",J::Number(guard.principalRawOverMu)}}));
    }
    return J::List(rows);
}
}

int main(int argc,char **argv)
{
    try
    {
        if(argc!=4 || fs::exists(argv[3]))
            throw std::runtime_error("usage: WJH27HillGeometryProbe session.xml table.json new-output.json");
        std::string xml=argv[1],opt=xml+".opt",program="wj-hill-geometry",option="--use-opt-file";
        char *args[]{program.data(),option.data(),opt.data(),xml.data(),nullptr};
        auto session=LibUtilities::SessionReader::CreateInstance(4,args);
        auto graph=SpatialDomains::MeshGraphIO::Read(session);
        std::vector<MultiRegions::ContFieldSharedPtr> fields;
        for (const std::string name : {"u","v","p"})
            fields.push_back(MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
                session,graph,name,true,false));
        const Reference ref(argv[2],ReferencePolicy::ExactInnerFiniteNuV1);
        Geometry geometry;
        geometry.H=ref.Parameters().H;
        geometry.d0=ref.Parameters().d0;
        geometry.h=ref.Parameters().roughness;
        auto &map=geometry.normalMap.emplace();
        map.kind=SolverUtils::RANSWallNormalMap::Kind::CosineHill;
        map.compositeTunnel=true;
        map.tunnelTop=20.;
        map.displacement=geometry.d0;
        map.height=geometry.H;
        map.lengthX=2.5;
        session->LoadParameter("WJHillHeight",map.amplitude,1.);

        const MixedFlow flow(fields,geometry);
        const RuleFactory rules(flow.Core(),6,8);
        const Vec first=flow.Project([&](double a,double z) {
            const double r=z-geometry.d0;
            return .1+.23*a+.07*a*r-.12*r*r;
        },0);
        const Vec second=flow.Project([&](double a,double z) {
            const double r=z-geometry.d0;
            return -.04+.11*r+.027*a*a;
        },1);
        std::vector<std::string> observations,faces,cutoff;
        for (const auto &e : flow.Elements())
        {
            const std::array<double,4> bounds{e.xl,e.xr,e.rl,e.rr};
            const bool outer=e.rl>=geometry.H-1e-12;
            if (std::abs(e.rl-geometry.H)<1e-12 || std::abs(e.rr-geometry.H)<1e-12)
            {
                const int side=outer?2:3;
                const auto face=rules.Face(e,side,16);
                faces.push_back(J::Object({{"bounds",J::Numbers(bounds)},{"side",J::Number(side)},
                    {"measure",J::Number(face.w.sum())},{"unit_normal",J::Numbers(face.Outward(3))}}));
            }
            // Both feet, the crest, the matching and outer regions use the
            // actual native basis and signed coefficients at off-grid points.
            if ((e.rr-e.rl<.0075 && std::abs(e.rr-geometry.H)>1e-12) ||
                e.xl < -5. || e.xr > 5.)
                continue;
            Vec locals(16), localw(16);
            for(int k=0;k<16;++k)
            {
                locals[k]=first[e.u[k]];
                localw[k]=second[e.u[k]];
            }
            for(double fraction : {.23,.61})
            {
                const double a=e.xl+fraction*(e.xr-e.xl);
                const double rho=e.rl+(.15+fraction)*(e.rr-e.rl);
                const auto M=geometry.At(a,rho,e.xl,e.xr,outer);
                const auto rule=rules.Points(e,Vec::Constant(1,a),Vec::Constant(1,rho+geometry.d0));
                observations.push_back(J::Object({{"bounds",J::Numbers(bounds)},{"a",J::Number(a)},
                    {"rho",J::Number(rho)},{"det",J::Number(M.determinant)},{"position",J::Numbers(M.point)},
                    {"inverse",J::Rows(M.inverse)},{"normal",J::Numbers(M.normal)},
                    {"grad",J::Numbers(std::array<double,2>{rule.Dx.row(0).dot(locals),
                                                          rule.Dz.row(0).dot(locals)})},
                    {"divergence",J::Number(rule.Dx.row(0).dot(locals)+rule.Dz.row(0).dot(localw))}}));
            }
        }
        const double Kin=ref.Parameters().KIn;
        for (int i=1;i<40;++i)
        {
            const double k=Kin*(.9+i*.005);
            const auto below=ref.DependentEpsilon(k,.01);
            const auto above=ref.DependentEpsilon(k,.04);
            cutoff.push_back(J::Object({{"K",J::Number(k)},{"c",J::Number(ref.Cutoff(k))},
                {"dc",J::Number(ref.CutoffDerivative(k))},{"EwK_below",J::Number(below.dK)},
                {"EwK_above",J::Number(above.dK)},{"traceK",J::Number(ref.AmplitudeDerivative(k)/geometry.H)}}));
        }
        std::ofstream output(argv[3]);
        output<<J::Object({{"height",J::Number(map.amplitude)},{"points",J::List(observations)},
                     {"faces",J::List(faces)},{"rotation",RotationObservations()},
                     {"cutoff_derivatives",J::List(cutoff)}})<<'\n';
        session->Finalise();
        return output ? 0 : 1;
    }
    catch(const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 2;
    }
}
