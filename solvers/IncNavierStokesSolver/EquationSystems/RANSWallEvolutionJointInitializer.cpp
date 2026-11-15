///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Diagnostic-only algebra; compiled into the resolution probe, not production.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <fstream>
#include <numeric>
#include <chrono>
#include <set>
namespace Nektar
{
using namespace SolverUtils;
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;
std::string RANSWallEvolutionPatch::JointInitializerAlgebra(const Vector &cold,const Vector &start,const WallMatrix &Z,const std::string &path) const
{
    Require(CalibratedStandard() && m_seedM && Z.size()==41,"H2i explicit retained 41-direction algebra");
    const int nv=Z.size(),n=m_seedFree.size();const auto begun=std::chrono::steady_clock::now();
    auto binary=[&](const std::string &name,const WallMatrix &a){std::ofstream f(path+"/"+name+".f64",std::ios::binary);for(const auto &v:a)f.write(reinterpret_cast<const char*>(v.data()),v.size()*sizeof(double));f.close();Require(f.good(),"joint algebra binary write");return J::Object({{"file",J::Quote(name+".f64")},{"shape",J::Numbers(std::vector<int>{int(a.size()),int(a.empty()?0:a[0].size())})},{"dtype",J::Quote("<f8")},{"sha256",J::Quote(RANSWallR3::HashFile(path+"/"+name+".f64"))}});};
    WallMatrix free(nv+1,Vector(n));std::set<int> allowed(m_seedFree.begin(),m_seedFree.end());
    for(int i=0;i<n;++i)free[0][i]=start[m_seedFree[i]]-cold[m_seedFree[i]];
    Vector qError(nv),dError(nv);double fixed=0;
    for(int j=0;j<nv;++j){Require(Z[j].size()==Size(),"H2i velocity dimension");for(int i=0;i<n;++i)free[j+1][i]=Z[j][m_seedFree[i]];for(int i=0;i<Size();++i)if(!allowed.count(i))fixed=std::max(fixed,std::abs(Z[j][i]));for(double x:Multiply(m_seedQ,free[j+1]))qError[j]=std::max(qError[j],std::abs(x));}
    Require(fixed==0 && *std::max_element(qError.begin(),qError.end())<1e-9,"H2i direction violates essential or independent continuity rows");
    WallMatrix mass(nv+1,Vector(nv+1));for(int j=0;j<=nv;++j){const auto v=Multiply(m_seedMass,free[j]);for(int i=0;i<=nv;++i)mass[i][j]=std::inner_product(free[i].begin(),free[i].end(),v.begin(),0.);}
    // Original regional continuity rows, including all dependent/shared rows.
    WallMatrix D(size[InnerPressure()]+size[OuterPressure()],Vector(n));std::vector<int> index(Size(),-1);for(int i=0;i<n;++i)index[m_seedFree[i]]=i;
    for(const auto *set:{&m_bodyInner,&m_outer})for(const auto &entry:*set)
    {const auto &q=entry.v;const auto &p=entry.p;const int off=q.coordinate[2]<m_parameters.geometry.height?0:size[InnerPressure()];for(int i=0;i<p.ids.size();++i)for(int a=0;a<2;++a)for(int j=0;j<q.ids.size();++j)if(index[offset[a]+q.ids[j]]>=0)D[off+p.ids[i]][index[offset[a]+q.ids[j]]]+=entry.weight*p.phi[i]*q.derivative[PhysicalAxis(a)][j];}
    for(int j=0;j<nv;++j)for(double x:Multiply(D,free[j+1]))dError[j]=std::max(dError[j],std::abs(x));
    Require(*std::max_element(dError.begin(),dError.end())<1e-10,"H2i direction violates a full regional continuity row");
    std::vector<std::pair<std::string,std::string>> files{{"velocity_mass",binary("velocity_mass",mass)},{"mass_free",binary("mass_free",m_seedMass)},{"continuity_Q",binary("continuity_Q",m_seedQ)},{"continuity_D",binary("continuity_D",D)}};
    for(int order:{16,20,24})
    {
        WallMatrix div(nv+1,Vector(nv+1));
        for(const auto &q:Quadrature(0,order))
        {Vector v(nv+1);for(int j=0;j<=nv;++j){const auto &s=j?Z[j-1]:start;v[j]=q.Gradient(s,offset[0])[0]+q.Gradient(s,offset[1])[2];}for(int i=0;i<=nv;++i)for(int j=0;j<=nv;++j)div[i][j]+=q.weight*v[i]*v[j];}
        files.push_back({"divergence_GL"+std::to_string(order),binary("divergence_GL"+std::to_string(order),div)});
    }
    return J::Object({{"version",J::Quote("JointNativeInitializerAlgebraV1")},{"velocity_directions",J::Number(nv)},{"files",J::Object(files)},{"free_velocity_rows",J::Numbers(m_seedFree)},{"independent_continuity_rank",J::Number(m_seedQ.size())},{"all_regional_continuity_rows",J::Number(D.size())},{"independent_row_errors",J::Numbers(qError)},{"all_regional_row_errors",J::Numbers(dError)},{"fixed_direction_max",J::Number(fixed)},{"B_sha256",J::Quote(RANSWallR3::Hash(start))},{"cold_sha256",J::Quote(RANSWallR3::Hash(cold))},{"seconds",J::Number(std::chrono::duration<double>(std::chrono::steady_clock::now()-begun).count())},{"governing_solves","0"}});
}

std::string RANSWallEvolutionPatch::JointInitializerAccessibility(const Vector &s,const WallMatrix &Z,const std::vector<InitializerWitness> &witnesses,const RANSWallEvolutionData &data) const
{
    Require(CalibratedStandard() && m_seedM && Z.size()==41 && witnesses.size()<=6,"H2i bounded native accessibility");
    const int n=m_seedFree.size();WallMatrix z(41,Vector(n)),G(41,Vector(41));for(int j=0;j<41;++j)for(int i=0;i<n;++i)z[j][i]=Z[j][m_seedFree[i]];
    for(int j=0;j<41;++j){const auto m=Multiply(m_seedMass,z[j]);for(int k=0;k<41;++k)G[k][j]=std::inner_product(z[k].begin(),z[k].end(),m.begin(),0.);}
    const Factor gram(G);std::vector<std::string> out;
    for(const auto &w:witnesses)
    {
        const auto row=GeometryPoint(0,w.geometry,w.chart);const bool inner=fields[0]->Elements()[row.element].upper[2]<=m_parameters.geometry.height;const double k=GeometryPoint(K(),w.geometry,w.chart).Value(s,offset[K()]);const auto prof=StateProfile(s,w.chart,0,data);const double e=inner?prof->Dissipation(Distance(w.chart)):GeometryPoint(Epsilon(),w.geometry,w.chart).Value(s,offset[Epsilon()]);const double nt=inner?prof->Original().EddyViscosity(Distance(w.chart)):.09*k*k/e;
        Vector g(n);for(int a=0;a<2;++a)for(int j=0;j<row.ids.size();++j){const int global=offset[a]+row.ids[j];const auto it=std::lower_bound(m_seedFree.begin(),m_seedFree.end(),global);if(it==m_seedFree.end() || *it!=global)continue;Vector unit(Size());unit[global]=1.;const auto d=row.Gradient(unit,offset[a]);for(int b=0;b<3;++b)g[it-m_seedFree.begin()]+=-2*nt*w.eigenvector[PhysicalAxis(a)]*w.eigenvector[b]*d[b];}
        auto response=m_seedM->Apply(g);const auto mu=m_seedSchur->Apply(Multiply(m_seedQ,response));for(int j=0;j<mu.size();++j)for(int i=0;i<n;++i)response[i]-=m_seedInverse[j][i]*mu[j];const auto mr=Multiply(m_seedMass,response);const double norm=std::inner_product(response.begin(),response.end(),mr.begin(),0.);Vector overlap(41);for(int j=0;j<41;++j)overlap[j]=std::inner_product(z[j].begin(),z[j].end(),mr.begin(),0.);const auto a=gram.Apply(overlap);double error=0;for(double x:Multiply(m_seedQ,response))error=std::max(error,std::abs(x));
        out.push_back(J::Object({{"chart",J::Numbers(w.chart)},{"physical",J::Numbers(Physical(w.chart))},{"geometry",J::Number(w.geometry)},{"element",J::Number(row.element)},{"inner",J::Boolean(inner)},{"k",J::Number(k)},{"epsilon",J::Number(e)},{"eigenvector",J::Numbers(w.eigenvector)},{"Riesz_mass_squared",J::Number(norm)},{"eta_squared",norm>1e-24?J::Number(std::inner_product(a.begin(),a.end(),overlap.begin(),0.)/norm):"null"},{"Q_response_max",J::Number(error)},{"span_derivative",J::Numbers(overlap)}}));
    }
    return J::Object({{"version",J::Quote("JointNativeAccessibilityV1")},{"state_sha256",J::Quote(RANSWallR3::Hash(s))},{"witnesses",J::List(out)},{"scope",J::Quote("linear essential/weak-continuity accessibility; not progress under active covariance/cap/quadratic constraints")},{"governing_solves","0"}});
}

std::string RANSWallEvolutionPatch::FlatEnergyProjectionDiagnostic(const Vector &base,const RANSWallEvolutionData &data) const
{
    Require(!Mapped() && AnalyticalNeutralInflow() && !OpenBoundary(),"H2i H1 projection requires compatible flat verification case");
    auto s=base;const auto constraints=Constraints(base,0,data);const double H=m_parameters.geometry.height;const auto profile=Profile(InflowAtHeight(H).k,{lower[0],0,0},1);std::vector<std::string> results;
    for(int b:{0,K()})
    {
        const auto n=size[b];WallMatrix A=m_mass[b];Vector rhs(n);std::map<int,double> fixed;
        for(const auto &[i,v]:constraints)if(i>=offset[b] && i<offset[b]+n)fixed[i-offset[b]]=base[i];
        // A diagnostic matching trace constraint preserves endpoint K/epsilon.
        if(b==K())
        {const auto boundary=fields[b]->Boundary([](const auto&){return 0.;},[&](int axis,int side,const auto &el){return axis==2 && ((side && el.upper[2]==H)||(!side && el.lower[2]==H));});for(auto [i,v]:boundary)fixed[i]=base[offset[b]+i];}
        for(const auto &q:Quadrature(b,16))
        {
            const double z=q.coordinate[2];const bool inner=z<H;const auto ref=InflowAtHeight(std::max(z,H));const double value=b==0?(inner?profile->Original().Velocity(z,{InflowAtHeight(H).U,0},{0,0})[0]:ref.U):(inner?profile->ReferenceK(z):ref.k);const double derivative=b==0?(inner?profile->Original().VelocityGradient(z,{InflowAtHeight(H).U,0},{0,0})[0]:ref.dU):(inner?profile->ReferenceKGradient(z):0.);
            for(int i=0;i<q.ids.size();++i){rhs[q.ids[i]]+=q.weight*(q.phi[i]*value+H*H*q.derivative[2][i]*derivative);for(int j=0;j<q.ids.size();++j)for(int d:{0,2})A[q.ids[i]][q.ids[j]]+=q.weight*H*H*q.derivative[d][i]*q.derivative[d][j];}
        }
        std::vector<int> free;for(int i=0;i<n;++i)if(!fixed.count(i))free.push_back(i);WallMatrix F(free.size(),Vector(free.size()));Vector r(free.size());for(int i=0;i<free.size();++i){r[i]=rhs[free[i]];for(auto [j,v]:fixed)r[i]-=A[free[i]][j]*v;for(int j=0;j<free.size();++j)F[i][j]=A[free[i]][free[j]];}
        const auto x=Factor(F).Apply(r);double residual=0;for(int i=0;i<free.size();++i){s[offset[b]+free[i]]=x[i];double v=-r[i];for(int j=0;j<free.size();++j)v+=F[i][j]*x[j];residual=std::max(residual,std::abs(v));}
        results.push_back(J::Object({{"block",J::Number(b)},{"free_rows",J::Number(free.size())},{"fixed_rows",J::Number(fixed.size())},{"setup_linear_residual_max",J::Number(residual)}}));
    }
    const auto domain=ScalarDomain(s,0,data,true);return J::Object({{"version",J::Quote("SameSpaceFlatH1ProjectionV1")},{"base_sha256",J::Quote(RANSWallR3::Hash(base))},{"state_sha256",J::Quote(RANSWallR3::Hash(s))},{"state",J::Numbers(s)},{"H1_length",J::Number(H)},{"assembly_order","16"},{"fields_changed",J::Quote("free U and k only; retained prescribed and k interface traces; W/epsilon/pressures unchanged")},{"setup",J::List(results)},{"scalar_domain",ScalarDomainJSON(domain)},{"baseline",FlatSpatialDiagnostic(base,data)},{"projected",domain.admissible?FlatSpatialDiagnostic(s,data):"null"},{"governing_solves","0"},{"physical_steps","0"},{"adopted","false"}});
}
}
