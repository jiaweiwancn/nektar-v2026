// SPDX-License-Identifier: MIT
// Small native coefficient fixture, never a full-size CFD initialization.
std::string FlatInletChecks(RANSWallEvolutionPatch &patch,
    const LibUtilities::SessionReaderSharedPtr &serial,const SpatialDomains::MeshGraphSharedPtr &graph,
    const LibUtilities::CommSharedPtr &comm,const std::string &output)
{
    Require(patch.Size()<600 && patch.InitialVelocityMode()=="FlatInletTraceV1","initializer fixture must be small and explicitly enabled");
    auto bits=[](const Vector &a,const Vector &b){return a.size()==b.size() && (a.empty() || std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0);};
    Vector pre(patch.Size());const double top=patch.upper[2];
    auto u=[&](const RANSWallCoordinate &q){const double t=q[2]/top;return .3*t+.1*t*t+.02*q[0]*t/(patch.upper[0]-patch.lower[0]);};
    for(int b=0;b<patch.Gauge();++b)
    {
        const auto c=patch.fields[b]->Project([&](const RANSWallCoordinate &q){return b==0?u(q):b==1?0.:.01+.001*b+.0001*q[2]/top;});
        std::copy(c.begin(),c.end(),pre.begin()+patch.offset[b]);
    }
    const auto essential=patch.fields[0]->Boundary(u,[&](int axis,int side,const RANSWallGlobalField::Element &el){return (axis==0 && !side && el.lower[0]==patch.lower[0]) || (axis==2 && !side && el.lower[2]==0);});
    for(auto [id,value]:essential)pre[id]=value;
    // Include a genuine signed-zero essential normal coefficient.
    const auto records=patch.fields[0]->InletTraceElements();const int zero=records.front().ids.front();pre[patch.offset[1]+zero]=-0.;
    const auto frozenReference=patch.Reference(),frozenDefect=patch.Defect();
    const auto post=patch.ApplyInitialVelocity(pre),again=patch.ApplyInitialVelocity(post);
    Require(bits(post,again),"FlatInletTraceV1 is not bitwise idempotent");
    Require(bits(Vector(pre.begin()+patch.offset[2],pre.end()),Vector(post.begin()+patch.offset[2],post.end())) &&
        bits(frozenReference,patch.Reference()) && bits(frozenDefect,patch.Defect()),"initializer changed nonvelocity/reference data");
    Require(std::signbit(post[patch.offset[1]+zero]),"initializer lost prescribed signed zero");
    for(auto [id,value]:essential){(void)value;Require(bits(Vector{pre[id]},Vector{post[id]}),"initializer changed native essential trace bits");}
    double divergence=0,dx=0;std::vector<std::string> moments;
    for(int side=0;side<2;++side)
    {
        const int b=patch.InnerPressure()+side;Vector load(patch.size[b]);
        for(const auto &q:patch.Quadrature(b,16))
        {
            const auto gu=patch.Point(0,q.coordinate).Gradient(post,patch.offset[0]),gw=patch.Point(1,q.coordinate).Gradient(post,patch.offset[1]);
            const double div=gu[0]+gw[2];divergence=std::max(divergence,std::abs(div));dx=std::max(dx,std::abs(gu[0]));
            for(int j=0;j<q.ids.size();++j)load[q.ids[j]]+=q.weight*q.phi[j]*div;
        }
        Require(Maximum(load)<2e-15,"initializer synthetic Q1 continuity moments");moments.push_back(J::Numbers(load));
    }
    Require(divergence<1e-12 && dx<1e-12,"initializer physical divergence");
    auto expect=[&](auto action,const std::string &needle){bool failed=false;try{action();}catch(const std::exception &e){failed=std::string(e.what()).find(needle)!=std::string::npos;}Require(failed,"missing intended rejection: "+needle);};
    auto bad=records;bad.front().element.lower[2]+=.00001;
    expect([&]{RANSWallGlobalField::FlatInletTraceCoefficients(Vector(pre.begin(),pre.begin()+patch.size[0]),bad,patch.lower[0]);},"FlatInletTraceV1");
    bad=records;bad.front().element.degree+=1;
    expect([&]{RANSWallGlobalField::FlatInletTraceCoefficients(Vector(pre.begin(),pre.begin()+patch.size[0]),bad,patch.lower[0]);},"FlatInletTraceV1");
    bad=records;bad.front().ids[1]=bad.front().ids[0];
    expect([&]{RANSWallGlobalField::FlatInletTraceCoefficients(Vector(pre.begin(),pre.begin()+patch.size[0]),bad,patch.lower[0]);},"FlatInletTraceV1");
    // Change one native global orientation consistently and require equivariance.
    bad=records;const int flipped=records.front().ids[4];Vector v(pre.begin(),pre.begin()+patch.size[0]);v[flipped]*=-1;
    for(auto &el:bad)for(int j=0;j<el.ids.size();++j)if(el.ids[j]==flipped)el.signs[j]*=-1;
    auto oriented=RANSWallGlobalField::FlatInletTraceCoefficients(v,bad,patch.lower[0]);oriented[flipped]*=-1;
    Require(bits(oriented,Vector(post.begin(),post.begin()+patch.size[0])),"initializer native sign equivariance");
    auto nonzero=pre;nonzero[patch.offset[1]+zero]=1e-8;
    expect([&]{patch.ApplyInitialVelocity(nonzero);},"nonzero prescribed normal");
    auto map=patch.Mapping();map.slopeX=.1;expect([&]{patch.ValidateFlatInletTraceMap(map);},"horizontal zero-hill");
    map=patch.Mapping();map.kind=RANSWallNormalMap::Kind::CosineHill;map.amplitude=.1;expect([&]{patch.ValidateFlatInletTraceMap(map);},"horizontal zero-hill");
    map=patch.Mapping();map.heightSlopeX=.01;expect([&]{patch.ValidateFlatInletTraceMap(map);},"horizontal zero-hill");
    double amplitude=.01;serial->SetParameter("RANSWallAcceptedKAmplitude",amplitude);
    expect([&]{RANSWallEvolutionPatch p(serial,graph,comm,patch.Parameters(),patch.Mapping());},"forbids initial perturbations");
    amplitude=0.;serial->SetParameter("RANSWallAcceptedKAmplitude",amplitude);
    const auto active=RANSWallPhysicsSignature(patch,serial,0.,false,.051,{1.,0.},false);
    serial->SetSolverInfo("RANSWallInitialVelocity","Off");
    RANSWallEvolutionPatch off(serial,graph,comm,patch.Parameters(),patch.Mapping());
    const auto legacy=RANSWallPhysicsSignature(off,serial,0.,false,.051,{1.,0.},false);
    auto stripped=active;const std::string tag=":initial_velocity=FlatInletTraceV1:lifted_native_inlet_coefficients";
    const auto where=stripped.find(tag);Require(where!=std::string::npos,"initializer signature tag absent");stripped.erase(where,tag.size());
    Require(stripped==legacy && off.InitialVelocityMode()=="Off" && off.InitialVelocityCalls()==0,"default-off signature semantics changed");
    int callbacks=0;
    RANSWallTimeIntegration time(1,.01,[&](const Vector &,double,double)->Vector{++callbacks;throw std::runtime_error("fixture must not advance");},[](const Vector &x,double){return Vector(x.size());});
    time.Initialise(Vector(post.begin(),post.begin()+patch.DifferentialSize()),0.);
    RANSWallTimeIntegration::Checkpoint cp;cp.signature=active;cp.integration=time.Capture();cp.coupledState=post;cp.fixedReference=pre;cp.fixedDefect=Vector(pre.size());
    const auto path=RANSWallTimeIntegration::WriteCheckpoint(output+".small-checkpoint",cp);
    const auto read=RANSWallTimeIntegration::ReadCheckpoint(path,active,1,.01,patch.Size(),patch.DifferentialSize());
    Require(bits(read.coupledState,post) && bits(read.fixedReference,pre) && callbacks==0,"small initializer checkpoint bits/callbacks differ");
    expect([&]{RANSWallTimeIntegration::ReadCheckpoint(path,legacy,1,.01,patch.Size(),patch.DifferentialSize());},"signature");
    int negatives=0;for(const auto &el:records)for(double sign:el.signs)negatives+=sign<0;
    return J::Object({{"status",J::Quote("PASS")},{"scope",J::Quote("small native synthetic coefficient fixture, no CFD CP0/governing solve")},
        {"state_size",J::Number(patch.Size())},{"native_negative_signs",J::Number(negatives)},
        {"physical_divergence_max",J::Number(divergence)},{"all_regional_Q1_moments",J::List(moments)},
        {"idempotent","true"},{"essential_and_nonvelocity_bits","true"},{"map_partition_ID_normal_perturbation_rejections","true"},
        {"default_off_and_signature_rejection","true"},{"synthetic_checkpoint",J::Quote(path)},
        {"governing_attempts","0"},{"full_size_CP0s","0"},{"new_initialization_solves","0"}});
}
