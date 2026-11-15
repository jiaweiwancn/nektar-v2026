///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// H2i full retained directions, scalar-bearing B; no governing solves.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_JOINT_FEASIBILITY_CONTEXT_H
#define NEKTAR_RANS_JOINT_FEASIBILITY_CONTEXT_H
namespace Nektar::RANSH2i
{
using Patch=RANSWallEvolutionPatch;using Vector=RANSWallAlgebra::Vector;
using RANSWallAlgebra::Require;namespace J=RANSWallSolveJson;
inline Patch::WallMatrix ReadBasis(const boost::property_tree::ptree &packet,int count,int size)
{
    const auto file=packet.get<std::string>("basis_file");Require(RANSWallR3::HashFile(file)==packet.get<std::string>("basis_sha256"),"H2i saved basis checksum");Patch::WallMatrix z(count,Vector(size));std::ifstream in(file,std::ios::binary);for(auto &v:z)in.read(reinterpret_cast<char*>(v.data()),v.size()*sizeof(double));Require(in.good() && in.peek()==std::char_traits<char>::eof(),"H2i saved basis dimensions");return z;
}
inline std::string Export(const Patch &p,const Vector &B,const Vector &cold,const RANSWallEvolutionData &data,const boost::property_tree::ptree &packet)
{
    Require(p.CalibratedStandard(),"H2i requires calibrated standard model");const auto started=std::chrono::steady_clock::now();p.ValidateOuterInitializer(B,cold,data,true);
    const bool extended=packet.get<std::string>("operation")=="h2j_context";const int nv=extended?packet.get<int>("velocity_direction_count"):41;Require(nv>=41 && nv<=47,"joint velocity count");const int totalDirections=nv+16,ko=23+6*nv,eo=ko+8,width=eo+8;const auto z=ReadBasis(packet,totalDirections,B.size());const auto rows=p.InitializerRows(cold);const std::set<int> vrows(rows.velocity.begin(),rows.velocity.end()),krows(rows.k.begin(),rows.k.end()),erows(rows.epsilon.begin(),rows.epsilon.end());
    for(int a=0;a<totalDirections;++a)for(int i=0;i<B.size();++i)if(!(a<nv?vrows:(a<nv+8?krows:erows)).count(i))Require(z[a][i]==0,"H2i direction changes protected row");
    const auto path=packet.get<std::string>("context_directory");Patch::WallMatrix velocity(z.begin(),z.begin()+nv);const auto algebra=p.JointInitializerAlgebra(cold,B,velocity,path);
    const auto domain=p.ScalarDomain(B,0,data,true),cd=p.ScalarDomain(cold,0,data,true);const auto ph=RANSH2h::ProfileHash(domain);Require(ph==RANSH2h::ProfileHash(cd),"H2i B profiles differ from immutable cold");
    Vector bv,bs;for(const auto &v:packet.get_child("B_scalar_coefficients"))bs.push_back(v.second.get_value<double>());Require(bs.size()==16,"H2i scalar start dimension");
    std::vector<Vector> trial(3,Vector(nv));if(extended){int t=0;for(const auto &entry:packet.get_child("velocity_verification")){Require(t<3,"joint verification count");int j=0;for(const auto &v:entry.second){Require(j<nv,"joint verification width");trial[t][j++]=v.second.get_value<double>();}Require(j==nv,"joint trial dimension");++t;}Require(t==3,"three joint trials required");}else for(int t=0;t<3;++t)for(int j=0;j<nv;++j)trial[t][j]=(t+1)*1e-12*std::sin((j+1.)*(t+1.));
    std::vector<Vector> tests;for(int t=0;t<3;++t)
    {auto s=B;for(int i:rows.velocity)for(int j=0;j<nv;++j)s[i]+=trial[t][j]*z[j][i];for(int b:{0,1})for(int i:b==0?rows.k:rows.epsilon)for(int j=0;j<16;++j)s[i]-=.01*(t+1)*bs[j]*z[nv+j][i];p.ValidateOuterInitializer(s,cold,data,true);Require(ph==RANSH2h::ProfileHash(p.ScalarDomain(s,0,data,true)),"joint trial changed protected profile");tests.push_back(std::move(s));}
    std::ofstream points(path+"/points.f64",std::ios::binary),members(path+"/members.i32",std::ios::binary),memberWeights(path+"/member_weights.f64",std::ios::binary);std::map<std::string,int> seen,families;std::vector<std::string> names;int n=0,total=0;double gradError=0,scalarError=0,tensorError=0,trialMin=1e300,trialMax=-1e300;Vector trialCap(3),trialMass(3),trialDiv(3);
    auto emit=[&](const RANSWallGlobalPoint &v,bool inner,bool dependent,bool covariance,double dep,int radial,int tangent,int block,const std::string &family,int set,int original,double weight)
    {
        std::ostringstream key;key<<v.geometry<<':'<<v.element<<':'<<radial<<':'<<tangent<<':'<<block<<':'<<dependent<<':'<<RANSWallR3::Hash(v.coordinate.data(),3*sizeof(double))<<':'<<set;
        int id;auto found=seen.find(key.str());if(found!=seen.end())id=found->second;
        else
        {
            id=n++;seen[key.str()]=id;const auto kr=p.GeometryPoint(p.K(),v.geometry,v.coordinate);const auto ep=inner?kr:p.GeometryPoint(p.Epsilon(),v.geometry,v.coordinate);Vector r(v.coordinate.begin(),v.coordinate.end());const auto phys=p.Physical(v.coordinate);r.insert(r.end(),phys.begin(),phys.end());SolverUtils::RANSTensor g{};
            for(int a=0;a<2;++a){g[p.PhysicalAxis(a)]=v.Gradient(B,p.offset[a]);r.push_back(g[p.PhysicalAxis(a)][0]);r.push_back(g[p.PhysicalAxis(a)][2]);}for(int a=0;a<2;++a)r.push_back(v.Value(B,p.offset[a]));
            const double k=kr.Value(B,p.offset[p.K()]),kc=kr.Value(cold,p.offset[p.K()]);const auto profile=inner?p.StateProfile(B,v.coordinate,0,data):nullptr;const double e=inner?profile->Dissipation(v.coordinate[2]):(dependent?dep:ep.Value(B,p.offset[p.Epsilon()])),ec=inner||dependent?e:ep.Value(cold,p.offset[p.Epsilon()]);const double nt=inner?profile->Original().EddyViscosity(v.coordinate[2]):0;
            r.insert(r.end(),{k,kc,e,ec,nt,double(inner),double(covariance),double(dependent),weight});for(int a=0;a<2;++a)r.push_back(v.Value(cold,p.offset[a]));
            for(int j=0;j<nv;++j){for(int a=0;a<2;++a){const auto d=v.Gradient(z[j],p.offset[a]);r.push_back(d[0]);r.push_back(d[2]);}for(int a=0;a<2;++a)r.push_back(v.Value(z[j],p.offset[a]));}
            for(int j=nv;j<nv+8;++j){r.push_back(inner?0.:kr.Value(z[j],p.offset[p.K()]));}
            for(int j=nv+8;j<totalDirections;++j){r.push_back(inner||dependent?0.:ep.Value(z[j],p.offset[p.Epsilon()]));}
            Require(r.size()==width,"joint context record width");points.write(reinterpret_cast<const char*>(r.data()),r.size()*sizeof(double));
            for(int t=0;t<3;++t)
            {
                for(int a=0;a<2;++a){if(extended)trialCap[t]=std::max(trialCap[t],std::abs(v.Value(tests[t],p.offset[a])-v.Value(cold,p.offset[a])));auto d=g[p.PhysicalAxis(a)];for(int j=0;j<nv;++j)for(int b:{0,2})d[b]+=trial[t][j]*r[23+j*6+a*2+(b==2)];const auto direct=v.Gradient(tests[t],p.offset[a]);for(int b:{0,2})gradError=std::max(gradError,std::abs(d[b]-direct[b])/std::max(1.,std::abs(direct[b])));}
                double kt=k,et=e;for(int j=0;j<8;++j){kt-=.01*(t+1)*bs[j]*r[ko+j];et-=.01*(t+1)*bs[j+8]*r[eo+j];}const double actualK=kr.Value(tests[t],p.offset[p.K()]),actualE=inner||dependent?et:ep.Value(tests[t],p.offset[p.Epsilon()]);scalarError=std::max({scalarError,std::abs(kt-actualK)/kc,std::abs(et-actualE)/ec});if(!inner){trialMin=std::min({trialMin,actualK/kc,actualE/ec});trialMax=std::max({trialMax,actualK/kc,actualE/ec});}
            }
            if(!inner && covariance){SolverUtils::RANSTensor S{};for(int a=0;a<3;++a)for(int b=0;b<3;++b)S[a][b]=g[a][b]+g[b][a];const auto sv=J::Eigenvalues(S),ev=J::Eigenvalues(p.OuterClosure(k,e,g,v.coordinate).reynoldsStress);tensorError=std::max(tensorError,std::abs(ev[0]-((2./3)*k-.09*k*k/e*sv[2])));}
        }
        int f=-1;if(set==0){if(!families.count(family)){families[family]=names.size();names.push_back(J::Quote(family));}f=families[family];}
        const std::int32_t m[]{total++,id,set,original,f,v.geometry,v.element,radial,tangent,block,int(dependent),int(covariance)};members.write(reinterpret_cast<const char*>(m),sizeof(m));memberWeights.write(reinterpret_cast<const char*>(&weight),sizeof(double));
    };
    int ix=0;for(const auto &o:domain.observations){if(o.block==p.K() || o.block==p.Epsilon() || o.dependentTrace){const auto v=p.ResolveObservation(o);emit(*v.velocity,v.inner,o.dependentTrace,o.block==p.K()||o.dependentTrace,o.value,v.radialSide,v.tangentialSide,o.block,o.inventory,0,ix,o.weight);}++ix;}
    const int nativeMembers=total;for(int order:{20,24}){int i=0;for(const auto &q:p.Quadrature(0,order)){const bool inner=p.fields[0]->Elements()[q.element].upper[2]<=p.Parameters().geometry.height;emit(q,inner,false,true,0,0,0,p.K(),"GL"+std::to_string(order),order,i++,q.weight);}}
    points.close();members.close();memberWeights.close();Require(points.good() && members.good() && memberWeights.good(),"joint context output");Require(gradError<1e-10 && scalarError<1e-10 && tensorError<1e-11 && trialMin>=.5-1e-12 && trialMax<=2+1e-12,"H2i direct native response checks");
    if(extended)for(const auto &q:p.Quadrature(0,16))for(int t=0;t<3;++t){double div=0;for(int a=0;a<2;++a){const double d=q.Value(tests[t],p.offset[a])-q.Value(cold,p.offset[a]);trialMass[t]+=q.weight*d*d;div+=q.Gradient(tests[t],p.offset[a])[p.PhysicalAxis(a)];}trialDiv[t]+=q.weight*div*div;}
    if(extended)for(int t=0;t<3;++t)Require(std::sqrt(trialMass[t])<=.12+1e-10 && std::sqrt(trialDiv[t])<=.06030809942199109+1e-12 && trialCap[t]<=.5+1e-10,"joint meaningful trial violates physical caps");
    auto file=[&](const std::string &name,int rows,int width,const std::string &type){return J::Object({{"file",J::Quote(name)},{"shape",J::Numbers(std::vector<int>{rows,width})},{"dtype",J::Quote(type)},{"sha256",J::Quote(RANSWallR3::HashFile(path+"/"+name))}});};
    std::vector<std::pair<std::string,std::string>> blocks;auto block=[&](std::string name,int offset,int width,std::vector<int> shape){blocks.push_back({name,J::Object({{"offset",J::Number(offset)},{"width",J::Number(width)},{"shape",J::Numbers(shape)}})});};block("chart",0,3,{3});block("physical",3,3,{3});block("gradient_B",6,4,{2,2});block("velocity_B",10,2,{2});for(auto [name,off]:std::vector<std::pair<std::string,int>>{{"k_B",12},{"k_cold",13},{"epsilon_B",14},{"epsilon_cold",15},{"inner_nu",16},{"inner",17},{"covariance",18},{"dependent",19},{"weight",20}})block(name,off,1,{1});block("velocity_cold",21,2,{2});block("velocity_basis",23,6*nv,{nv,6});block("k_basis",ko,8,{8});block("epsilon_basis",eo,8,{8});
    return J::Object({{"version",J::Quote(extended?"FixedKJointContextV2":"FixedKJointContextV1")},{"native_coefficient_offsets",J::Numbers(p.offset)},{"native_coefficient_sizes",J::Numbers(p.size)},{"direction_counts",J::Object({{"velocity",J::Number(nv)},{"k_anchor","8"},{"epsilon","8"}})},{"observation_blocks",J::Object(blocks)},{"files",J::Object({{"points",file("points.f64",n,width,"<f8")},{"members",file("members.i32",total,12,"<i4")},{"member_weights",file("member_weights.f64",total,1,"<f8")}})},{"algebra",algebra},{"families",J::List(names)},{"native_memberships",J::Number(nativeMembers)},{"all_memberships",J::Number(total)},{"unique_evaluations",J::Number(n)},{"rows",J::Object({{"velocity",J::Numbers(rows.velocity)},{"k",J::Numbers(rows.k)},{"epsilon",J::Numbers(rows.epsilon)}})},{"model_specification",p.ModelMetadata()},{"B_sha256",J::Quote(RANSWallR3::Hash(B))},{"cold_sha256",J::Quote(RANSWallR3::Hash(cold))},{"basis_sha256",J::Quote(packet.get<std::string>("basis_sha256"))},{"profile_sha256",J::Quote(ph)},
      {"verification",J::Object({{"gradient_relative_error",J::Number(gradError)},{"scalar_relative_error",J::Number(scalarError)},{"native_tensor_error",J::Number(tensorError)},{"three_nonzero_trials_scalar_ratios",J::Numbers(Vector{trialMin,trialMax})},{"protected_rows_and_profiles","true"},{"meaningful_trial_mass_squared",J::Numbers(trialMass)},{"meaningful_trial_divergence_squared",J::Numbers(trialDiv)},{"meaningful_trial_sampled_caps",J::Numbers(trialCap)}})},{"seconds",J::Number(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count())},{"governing_solves","0"}});
}
}
#endif
