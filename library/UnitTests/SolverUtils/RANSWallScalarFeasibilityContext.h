///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Probe-only fixed-velocity context. No production initialization or evolution.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_SCALAR_FEASIBILITY_CONTEXT_H
#define NEKTAR_RANS_SCALAR_FEASIBILITY_CONTEXT_H
namespace Nektar::RANSH2h
{
using Patch=RANSWallEvolutionPatch;
using Vector=RANSWallAlgebra::Vector;
using RANSWallAlgebra::Require;
using SolverUtils::RANSTensor;
namespace J=RANSWallSolveJson;

inline void FixedB(const Patch &p,const Vector &s,const Vector &b,const Vector &cold)
{
    Require(s.size()==b.size() && b.size()==cold.size(),"H2h state dimensions");
    const auto rows=p.InitializerRows(cold);std::set<int> allowed(rows.k.begin(),rows.k.end());
    allowed.insert(rows.epsilon.begin(),rows.epsilon.end());
    for(int i=0;i<s.size();++i)if(!allowed.count(i))
        Require(std::memcmp(&s[i],&b[i],sizeof(double))==0,"H2h changed archived B non-scalar bits");
}
inline double EvalError(const RANSWallGlobalPoint &p,const Vector &s,int off)
{
    const double centre=s[off+p.ids[0]]/p.constant[0];long double v=std::abs(centre);
    for(int j=0;j<p.ids.size();++j)v+=std::abs(p.phi[j])*(std::abs(s[off+p.ids[j]])+std::abs(centre*p.constant[j]));
    return 64*std::numeric_limits<double>::epsilon()*double(v);
}
inline std::pair<double,double> RowBounds(const RANSWallGlobalPoint &p)
{
    double sum=0,norm=0,amp=1/std::abs(p.constant[0]);
    for(int j=0;j<p.ids.size();++j){sum+=p.phi[j]*p.constant[j];norm+=std::abs(p.phi[j]);amp+=std::abs(p.phi[j])*(1+std::abs(p.constant[j]/p.constant[0]));}
    norm+=std::abs((1-sum)/p.constant[0]);return {norm,amp};
}
inline std::string ProfileHash(const Patch::ScalarDomainResult &d)
{
    Vector v;for(const auto &p:d.pendingProfiles)
    {
        v.insert(v.end(),p.anchor.begin(),p.anchor.end());const auto &f=*p.profile;
        const auto &r=f.Reference();v.insert(v.end(),{r.cutoff,r.productionBudget,r.minimumBudget});
        for(double x:{0.,r.cutoff,f.Geometry().height})
            v.insert(v.end(),{f.ReferenceK(x),f.Dissipation(x),f.Original().EddyViscosity(x)});
    }
    return RANSWallR3::Hash(v);
}
inline std::string Export(const Patch &p,const Vector &b,const Vector &cold,
    const RANSWallEvolutionData &data,const boost::property_tree::ptree &packet)
{
    const auto started=std::chrono::steady_clock::now();
    Require(p.CalibratedStandard(),"H2h context requires calibrated standard model");
    FixedB(p,b,b,cold);p.ValidateOuterInitializer(b,cold,data,true);
    const auto path=packet.get<std::string>("context_directory"),file=packet.get<std::string>("scalar_basis_file");
    Require(RANSWallR3::HashFile(file)==packet.get<std::string>("scalar_basis_sha256"),"H2h saved scalar basis checksum");
    Patch::WallMatrix z(16,Vector(b.size()));std::ifstream in(file,std::ios::binary);
    for(auto &r:z)in.read(reinterpret_cast<char*>(r.data()),r.size()*sizeof(double));
    Require(in.good() && in.peek()==std::char_traits<char>::eof(),"H2h scalar basis dimensions");
    const auto rows=p.InitializerRows(cold);std::set<int> kr(rows.k.begin(),rows.k.end()),er(rows.epsilon.begin(),rows.epsilon.end());
    for(int a=0;a<16;++a)for(int j=0;j<b.size();++j)
        if(!(a<8?kr:er).count(j))Require(z[a][j]==0.,"H2h scalar basis changes protected row");
    const auto d=p.ScalarDomain(b,0,data,true),cd=p.ScalarDomain(cold,0,data,true);
    Require(d.admissible && cd.admissible && d.observations.size()==cd.observations.size(),"H2h scalar domain identity");
    const auto profileHash=ProfileHash(d);Require(profileHash==ProfileHash(cd),"H2h archived B profile differs from cold");
    Vector previous;for(const auto &entry:packet.get_child("B_scalar_coefficients"))previous.push_back(entry.second.get_value<double>());
    Require(previous.size()==16,"H2h previous scalar coordinates");
    std::vector<Vector> trials;for(double t:{.01,.03,.07})
    {
        Vector s=b;for(int i:rows.k)for(int a=0;a<16;++a)s[i]-=t*previous[a]*z[a][i];
        for(int i:rows.epsilon)for(int a=0;a<16;++a)s[i]-=t*previous[a]*z[a][i];
        FixedB(p,s,b,cold);p.ValidateOuterInitializer(s,cold,data,true);
        Require(ProfileHash(p.ScalarDomain(s,0,data,true))==profileHash,"H2h legal scalar increment changed profiles");trials.push_back(std::move(s));
    }
    std::ofstream out(path+"/points.f64",std::ios::binary),members(path+"/members.i32",std::ios::binary);
    std::map<std::string,int> families;std::vector<std::string> names;int n=0,native=0;double maxError=0,maxGradientChange=0,maxEigenError=0;
    auto emit=[&](const RANSWallGlobalPoint &v,bool inner,bool dependent,bool covariance,double dependentValue,
        int family,int radial,int tangential,int block,int original,double weight,int set)
    {
        const auto kp=p.GeometryPoint(p.K(),v.geometry,v.coordinate);const auto ep=inner?kp:p.GeometryPoint(p.Epsilon(),v.geometry,v.coordinate);
        const double k=kp.Value(b,p.offset[p.K()]),kc=kp.Value(cold,p.offset[p.K()]);
        const auto profile=inner?p.StateProfile(b,v.coordinate,0,data):nullptr;
        const double e=inner?profile->Dissipation(v.coordinate[2]):(dependent?dependentValue:ep.Value(b,p.offset[p.Epsilon()]));
        const double ec=inner || dependent?e:ep.Value(cold,p.offset[p.Epsilon()]);
        RANSTensor g{};for(int a=0;a<2;++a)g[p.PhysicalAxis(a)]=v.Gradient(b,p.offset[a]);
        Vector r(v.coordinate.begin(),v.coordinate.end());const auto phys=p.Physical(v.coordinate);r.insert(r.end(),phys.begin(),phys.end());
        for(const auto &row:g)r.insert(r.end(),row.begin(),row.end());
        r.insert(r.end(),{k,e,kc,ec,double(inner),double(covariance),double(dependent),weight,inner?profile->Original().EddyViscosity(v.coordinate[2]):0.});
        for(int a=0;a<8;++a)r.push_back(inner?0.:kp.Value(z[a],p.offset[p.K()]));
        for(int a=8;a<16;++a)r.push_back(inner || dependent?0.:ep.Value(z[a],p.offset[p.Epsilon()]));
        r.push_back(EvalError(kp,b,p.offset[p.K()]));r.push_back(inner || dependent?0.:EvalError(ep,b,p.offset[p.Epsilon()]));
        for(int a=0;a<16;++a)r.push_back(inner || (a>=8 && dependent)?0.:EvalError(a<8?kp:ep,z[a],p.offset[a<8?p.K():p.Epsilon()]));
        const auto kb=RowBounds(kp),eb=RowBounds(ep);r.insert(r.end(),{kb.first,eb.first,kb.second,eb.second});
        Require(r.size()==62,"H2h scalar context width");out.write(reinterpret_cast<const char*>(r.data()),r.size()*sizeof(double));
        const std::int32_t m[]{n++,set,original,family,v.geometry,v.element,radial,tangential,block,int(dependent)};members.write(reinterpret_cast<const char*>(m),sizeof(m));
        for(int t=0;t<trials.size();++t)
        {
            const double scale=std::vector<double>{.01,.03,.07}[t];double ka=k,ea=e;
            for(int a=0;a<8;++a)ka-=scale*previous[a]*r[24+a];
            for(int a=8;a<16;++a)ea-=scale*previous[a]*r[24+a];
            maxError=std::max(maxError,std::abs(ka-kp.Value(trials[t],p.offset[p.K()]))/std::max(1e-12,std::abs(kc)));
            if(!inner && !dependent)maxError=std::max(maxError,std::abs(ea-ep.Value(trials[t],p.offset[p.Epsilon()]))/std::max(1e-12,std::abs(ec)));
            for(int a=0;a<2;++a){const auto tg=v.Gradient(trials[t],p.offset[a]);for(int j=0;j<3;++j)maxGradientChange=std::max(maxGradientChange,std::abs(tg[j]-g[p.PhysicalAxis(a)][j]));}
        }
        if(!inner && covariance)
        {
            RANSTensor S{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)S[i][j]=g[i][j]+g[j][i];
            const auto ev=J::Eigenvalues(p.OuterClosure(k,e,g,v.coordinate).reynoldsStress),sv=J::Eigenvalues(S);
            maxEigenError=std::max(maxEigenError,std::abs(ev[0]-((2./3)*k-.09*k*k/e*sv[2])));
        }
    };
    int index=0;for(const auto &q:d.observations)
    {
        if(q.block==p.K() || q.block==p.Epsilon() || q.dependentTrace)
        {
            if(!families.count(q.inventory)){families[q.inventory]=names.size();names.push_back(J::Quote(q.inventory));}
            const auto o=p.ResolveObservation(q);emit(*o.velocity,o.inner,q.dependentTrace,q.block==p.K() || q.dependentTrace,q.value,families[q.inventory],o.radialSide,o.tangentialSide,q.block,index,q.weight,0);++native;
        }
        ++index;
    }
    int qi=0;for(const auto &q:p.Quadrature(0,20))
    {const bool inner=p.fields[0]->Elements()[q.element].upper[2]<=p.Parameters().geometry.height;emit(q,inner,false,true,0,-1,inner?-1:1,0,p.K(),qi++,q.weight,1);}
    out.close();members.close();Require(out.good() && members.good(),"H2h context output");
    Require(maxError<1e-10 && maxGradientChange==0. && maxEigenError<1e-11,"H2h direct-native affine/eigenvalue verification failed");
    auto fileJSON=[&](const std::string &file,int width,const std::string &type){return J::Object({{"file",J::Quote(file)},{"shape",J::Numbers(std::vector<int>{n,width})},{"dtype",J::Quote(type)},{"sha256",J::Quote(RANSWallR3::HashFile(path+"/"+file))}});};
    return J::Object({{"version",J::Quote("FixedBScalarContextV1")},{"state_sha256",J::Quote(RANSWallR3::Hash(b))},{"base_sha256",J::Quote(RANSWallR3::Hash(cold))},{"basis_sha256",J::Quote(packet.get<std::string>("scalar_basis_sha256"))},{"model_specification",p.ModelMetadata()},
        {"files",J::Object({{"points",fileJSON("points.f64",62,"<f8")},{"members",fileJSON("members.i32",10,"<i4")}})},
        {"families",J::List(names)},{"native_memberships",J::Number(native)},{"GL20_points",J::Number(qi)},{"profile_sha256",J::Quote(profileHash)},
        {"verification",J::Object({{"direct_legal_trials","3"},{"relative_affine_error",J::Number(maxError)},{"gradient_change",J::Number(maxGradientChange)},{"native_tensor_eigenvalue_error",J::Number(maxEigenError)},{"frozen_non_scalar_bits","true"},{"protected_profile_invariant","true"}})},
        {"scalar_rows",J::Object({{"k",J::Numbers(rows.k)},{"epsilon",J::Numbers(rows.epsilon)}})},
        {"seconds",J::Number(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count())},{"governing_solves","0"},{"physical_commits","0"}});
}
}
#endif
