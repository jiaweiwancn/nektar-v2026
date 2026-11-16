// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_H25_SEED_H
#define NEKTAR_RANS_WALL_H25_SEED_H

// Explicit H25 transient producer. Ordinary Restore and StepZeroIMEX1 retain
// their old strict contracts; this is never an evolved-state dt conversion.
double H25Spline(const std::vector<double> &t,double x,int i,int p)
{
    if(p==0)return t[i]<=x && x<t[i+1]?1.:0.;
    return (x-t[i])/(t[i+p]-t[i])*H25Spline(t,x,i,p-1)+
        (t[i+p+1]-x)/(t[i+p+1]-t[i+1])*H25Spline(t,x,i+1,p-1);
}
double H25SplineDerivative(const std::vector<double> &t,double x)
{return 3/(t[3]-t[0])*H25Spline(t,x,0,2)-3/(t[4]-t[1])*H25Spline(t,x,1,2);}

std::string H25Seed(RANSWallEvolutionPatch &patch,const LibUtilities::SessionReaderSharedPtr &serial,
    double us,const Tree &packet,RANSWallEvolutionData &data)
{
    Require(patch.Size()==24688 && patch.MappedBlock() && patch.InitialVelocityMode()=="FlatInletTraceV1",
        "H25 seed requires the unchanged full700-cell compatible control");
    double dt;serial->LoadParameter("TimeStep",dt);double sourceDt=packet.get<double>("source_dt");
    Require(dt>0 && sourceDt==.0025,"H25 seed source is the immutable original fine CP0");
    serial->SetParameter("TimeStep",sourceDt);OpenImportedState imported;
    try{imported=ImportOpenState(patch,serial,us,packet,data);}catch(...){serial->SetParameter("TimeStep",dt);throw;}
    serial->SetParameter("TimeStep",dt);const auto &source=imported.cp;
    Require(source.integration.steps==0 && RANSWallMonitor::Bits(source.integration.time)==RANSWallMonitor::Bits(0.) &&
        source.integration.order==1 && source.eta==0 && RANSWallR3::Hash(source.coupledState)==
        "631caa32a59ff43253f1fcac4effd955eb68b8529d9c28f29e615cf14ca0f077","H25 seed is strictly original IMEX1 spatial state0");
    auto state=source.coupledState;const bool perturbed=packet.get<bool>("perturbed");
    std::vector<double> xt,zt;double amplitude=0.;
    if(perturbed)
    {
        std::set<double> x,z;
        for(const auto &e:patch.fields[0]->Elements())
        {x.insert(e.lower[0]);x.insert(e.upper[0]);z.insert(patch.Physical(e.lower)[2]);z.insert(patch.Physical(e.upper)[2]);}
        std::vector<double> xs(x.begin(),x.end()),zs(z.begin(),z.end());
        Require(xs.size()==15 && zs.size()==51,"H25 transient knot inventory differs");
        const auto h=std::find(zs.begin(),zs.end(),patch.Physical({0.,0.,.15})[2]);
        Require(h!=zs.end() && h-zs.begin()>=2 && zs.end()-h>=3,"H25 matching knot missing");
        xt.assign(xs.begin()+4,xs.begin()+9);zt.assign(h-2,h+3);
        const double fn=H25Spline(xt,xt[2],0,3),gn=H25Spline(zt,*h,0,3);
        const double gd=H25SplineDerivative(zt,(zt[1]+zt[2])/2);
        Require(fn>0 && gn>0 && gd!=0,"degenerate H25 compact cubic seed");
        amplitude=.001/(fn*gd);const double K=81/(2*std::sqrt(143.))*us*us;
        std::vector<RANSWallEvolutionFunction> functions(patch.fields.size(),[](const RANSWallCoordinate&,double){return 0.;});
        functions[0]=[=](const RANSWallCoordinate &p,double){return amplitude*H25Spline(xt,p[0],0,3)*H25SplineDerivative(zt,p[2]);};
        functions[1]=[=](const RANSWallCoordinate &p,double){return -amplitude*H25SplineDerivative(xt,p[0])*H25Spline(zt,p[2],0,3);};
        functions[2]=[=](const RANSWallCoordinate &p,double){return .01*K*H25Spline(xt,p[0],0,3)/fn*H25Spline(zt,p[2],0,3)/gn;};
        const double d0=patch.Mapping().displacement,kappa=patch.Parameters().model.Kappa();
        functions[3]=[=](const RANSWallCoordinate &p,double){return .01*us*us*us/(kappa*std::max(.15,p[2]-d0))*
            H25Spline(xt,p[0],0,3)/fn*H25Spline(zt,p[2],0,3)/gn;};
        const auto increment=patch.Project(functions,0.);
        for(int j=0;j<patch.DifferentialSize();++j)if(increment[j]!=0.)state[j]+=increment[j];
        // This is declared new initial-data construction. Native current-K
        // boundary/trace lifting happens once BEFORE native history production.
        state=patch.Lift(state,0.,data);
        Require(RANSWallR3::Hash(Vector(state.begin()+patch.offset[patch.InnerPressure()],state.end()))==
            RANSWallR3::Hash(Vector(source.coupledState.begin()+patch.offset[patch.InnerPressure()],source.coupledState.end())),
            "transient seed changed stored regional pressure");
        for(const auto &[row,target]:patch.Constraints(source.coupledState,0.,data))
            if(row<patch.offset[patch.Epsilon()])Require(RANSWallR3::Hash(Vector{state[row]})==
                RANSWallR3::Hash(Vector{source.coupledState[row]}),"transient seed changed prescribed velocity/k bits");
    }
    const auto signature=RANSWallPhysicsSignature(patch,serial,0.,false,us,{1.,0.},false);
    std::size_t explicitCalls=0;
    RANSWallTimeIntegration native(1,dt,[](const Vector&,double,double)->Vector{throw std::runtime_error("H25 seed attempted governing solve");},
        [&](const Vector &q,double t){++explicitCalls;auto full=state;std::copy(q.begin(),q.end(),full.begin());return patch.ExplicitBody(full,t,data);});
    native.Initialise(Vector(state.begin(),state.begin()+patch.DifferentialSize()),0.);
    RANSWallTimeIntegration::Checkpoint cp;cp.signature=signature;cp.integration=native.Capture();cp.coupledState=state;
    cp.fixedReference=source.fixedReference;cp.fixedDefect=source.fixedDefect;cp.eta=source.eta;
    const auto scalar=patch.ScalarDomain(state,0.,data,true);Require(scalar.admissible,"H25 seed scalar/current-K domain failed");
    patch.ValidateOutlet(state,"H25 native transient seed");
    const auto audit=patch.CompleteRestrictedAudit(state,state,0.,0.,data);
    Require(audit.at("wj_hard_admissible")[0]==1,"H25 seed hard audit failed");
    const auto bound=patch.ApplyPimpleKBound(state,state,0.,data);
    Require(bound.success && !bound.activated && RANSWallR3::Hash(bound.state)==RANSWallR3::Hash(state),"H25 seed bound changed state");
    const auto checkpoint=RANSWallTimeIntegration::WriteCheckpoint(packet.get<std::string>("checkpoint_base"),cp);
    const int fd=::open((checkpoint+"/State.bin").c_str(),O_RDONLY);Require(fd>=0 && ::fsync(fd)==0,"H25 seed fsync failed");::close(fd);
    RANSWallMonitor::CheckpointDirectories(checkpoint);
    std::vector<std::string> histories,fields;
    for(const auto &v:cp.integration.history)histories.push_back(J::Numbers(v));
    for(int b=0;b<patch.Gauge();++b)fields.push_back(Field(patch,b));
    std::vector<std::pair<std::string,std::string>> a;for(const auto &[n,v]:audit)a.push_back({n,J::Numbers(v)});
    Require(RANSWallR3::HashFile(imported.cpPath+"/State.bin")==imported.checkpointHash &&
        RANSWallR3::HashFile(imported.reportPath)==imported.reportHash,"H25 seed modified original bytes");
    return J::Object({{"schema",J::Quote("WJH25NativeSeedV1")},{"status",J::Quote("PASS")},{"failure",J::Quote("")},
        {"steps","0"},{"time","0"},{"dt",J::Number(dt)},{"signature",J::Quote(signature)},
        {"state",J::Numbers(state)},{"fixed_reference",J::Numbers(cp.fixedReference)},{"fixed_defect",J::Numbers(cp.fixedDefect)},
        {"eta","0"},{"native_history",J::List(histories)},{"native_history_times",J::Numbers(cp.integration.times)},
        {"fields",J::List(fields)},{"offsets",J::Numbers(patch.offset)},{"map_parameters",J::Numbers(MapParameters(patch.Mapping()))},
        {"parameters",J::Numbers(Read(imported.native,"parameters"))},{"geometry_signature",J::Quote(patch.Mapping().Signature())},
        {"quadrature_orders",J::Numbers(Vector{double(patch.Parameters().tangentialOrder),double(patch.Parameters().normalOrder)})},
        {"accepted_states",J::List({J::Numbers(state)})},{"wj_audit_packets",J::List({patch.WJAuditPacket(audit,state,state,0.,0.,"startup")})},
        {"native_hard_audit",J::Object(a)},{"native_full_audit_status",J::Quote("PASS")},{"initial_k_bound_check",bound.report},
        {"last_stage_history","[]"},{"last_stage_alpha","0"},{"last_stage_time","0"},{"checkpoint",J::Quote(checkpoint)},
        {"source_checkpoint_sha256",J::Quote(imported.checkpointHash)},{"source_report_sha256",J::Quote(imported.reportHash)},
        {"common_original_qh0_sha256",J::Quote(RANSWallR3::Hash(source.coupledState))},{"seed_state_sha256",J::Quote(RANSWallR3::Hash(state))},
        {"perturbed",J::Boolean(perturbed)},{"streamfunction_x_knots",J::Numbers(xt)},{"streamfunction_z_knots",J::Numbers(zt)},
        {"streamfunction_amplitude",J::Number(amplitude)},{"scalar_fraction",J::Number(perturbed?.01:0.)},
        {"velocity_scale",J::Number(perturbed?.001:0.)},{"native_history_initializations","1"},{"explicit_initialization_calls",J::Number(explicitCalls)},
        {"spatial_work",J::Numbers(patch.InitializationWork())},{"initial_velocity_applications",J::Number(patch.InitialVelocityCalls())},{"governing_solves","0"}});
}
#endif
