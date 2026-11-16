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

OpenImportedState H25ImportFine(const Tree &packet,const LibUtilities::CommSharedPtr &comm)
{
    // Verify the immutable fine source in ITS context, never under target
    // geometry/quadrature. No ordinary importer or restart rule is relaxed.
    const auto &pin=packet.get_child("source_xml");RANSWallMonitor::CheckPin(pin);
    Require(pin.get<std::string>("sha256")=="515caab634d2d4391bb77d3a18be4e4ee375449e11ff2aea62d94275c393e51d",
        "H25 target seed source XML is not the immutable original fine CP0");
    char program[]="h25-fine-source";char *args[]{program,nullptr};
    std::vector<std::string> files{pin.get<std::string>("path")};
    auto serial=LibUtilities::SessionReader::CreateInstance(1,args,files,LibUtilities::CommSerial::create(1,args));serial->InitSession();
    serial->SetSolverInfo("GlobalSysSoln","DirectFull");
    RANSWallPatchParameters p;
    serial->LoadParameter("RANSWallChi",p.chi,1.);serial->LoadParameter("RANSWallHeight",p.geometry.height,2.);
    serial->LoadParameter("RANSWallRoughness",p.geometry.roughness,.1);serial->LoadParameter("Kinvis",p.geometry.molecularNu,1.5e-5);
    serial->LoadParameter("RANSWallSampleHeight",p.geometry.sampleHeight,1.);serial->LoadParameter("RANSWallTailStart",p.geometry.tailStart,1.);
    serial->LoadParameter("RANSWallTangentialOrder",p.tangentialOrder,8);serial->LoadParameter("RANSWallNormalOrder",p.normalOrder,16);
    LoadWJNativeConfiguration(serial,p);
    auto graph=SpatialDomains::MeshGraphIO::Read(serial);
    RANSWallEvolutionPatch fine(serial,graph,comm,p,LoadRANSWallNormalMap(serial,2,p.geometry.height));
    Require(fine.Size()==24688 && fine.MappedBlock() && fine.InitialVelocityMode()=="FlatInletTraceV1" &&
        p.tangentialOrder==12 && p.normalOrder==16,"H25 source compatible fine context differs");
    RANSWallEvolutionData data(2);fine.LoadOpenBoundaryData(serial,data);
    double us;serial->LoadParameter("RANSWallFrictionVelocity",us);
    return ImportOpenState(fine,serial,us,packet,data);
}

Vector H25TargetReference(RANSWallEvolutionPatch &patch,const LibUtilities::SessionReaderSharedPtr &serial,
    RANSWallEvolutionData &data)
{
    // Same physical InitialConditions Functions and reference/defect ordering
    // as native v_DoInitialise. Freeze reference BEFORE FlatInletTraceV1.
    auto function=[&](const std::string &family,const std::string &variable,bool required)
        ->RANSWallEvolutionFunction
    {
        Require(!required || serial->DefinesFunction(family,variable),"H25 missing target physical Function");
        if(!serial->DefinesFunction(family,variable))return [](const RANSWallCoordinate&,double){return 0.;};
        Require(serial->GetFunctionType(family,variable)==LibUtilities::eFunctionTypeExpression,"H25 target Functions must be expressions");
        const auto equation=serial->GetFunction(family,variable);
        return [equation,&patch](const RANSWallCoordinate &chart,double t)
        {const auto p=patch.Physical(chart);return equation->Evaluate(p[0],p[2],0.,t);};
    };
    std::vector<RANSWallEvolutionFunction> functions;
    for(int b=0;b<patch.Gauge();++b)
    {
        const std::string variable=b<2?serial->GetVariable(b):b==patch.K()?"k":b==patch.Epsilon()?"epsilon":"p";
        auto inner=function("InitialConditionsInner",variable,b!=patch.Epsilon() && b!=patch.OuterPressure());
        auto outer=function("InitialConditionsOuter",variable,b!=patch.InnerPressure());
        if(b==patch.InnerPressure())functions.push_back(inner);
        else if(b==patch.OuterPressure() || b==patch.Epsilon())functions.push_back(outer);
        else functions.push_back([inner,outer,H=patch.Parameters().geometry.height](const RANSWallCoordinate &p,double t)
            {return p[2]<=H?inner(p,t):outer(p,t);});
    }
    auto reference=patch.Lift(patch.Project(functions,0.),0.,data);
    patch.SetReference(reference,data,0.,0.);
    return patch.ApplyInitialVelocity(reference);
}

std::string H25Seed(RANSWallEvolutionPatch &patch,const LibUtilities::SessionReaderSharedPtr &serial,
    double us,const Tree &packet,RANSWallEvolutionData &data)
{
    Require((patch.Size()==24688 || patch.Size()==12445) && patch.MappedBlock() && patch.InitialVelocityMode()=="FlatInletTraceV1",
        "H25 seed requires the qualified700/350-cell compatible control");
    double dt;serial->LoadParameter("TimeStep",dt);double sourceDt=packet.get<double>("source_dt");
    Require(dt>0 && sourceDt==.0025,"H25 seed source is the immutable original fine CP0");
    const bool target=packet.get<std::string>("operation")=="h25_target_seed";
    OpenImportedState imported;
    if(target)imported=H25ImportFine(packet,patch.Comm());
    else
    {
        Require(patch.Size()==24688,"fine dt-only seed cannot import into a different mesh");
        RANSWallMonitor::CheckPin(packet.get_child("source_xml"));
        serial->SetParameter("TimeStep",sourceDt);
        try{imported=ImportOpenState(patch,serial,us,packet,data);}catch(...){serial->SetParameter("TimeStep",dt);throw;}
        serial->SetParameter("TimeStep",dt);
    }
    const auto &source=imported.cp;
    Require(source.integration.steps==0 && RANSWallMonitor::Bits(source.integration.time)==RANSWallMonitor::Bits(0.) &&
        source.integration.order==1 && source.eta==0 && RANSWallR3::Hash(source.coupledState)==
        "631caa32a59ff43253f1fcac4effd955eb68b8529d9c28f29e615cf14ca0f077","H25 seed is strictly original IMEX1 spatial state0");
    auto state=target?H25TargetReference(patch,serial,data):source.coupledState;
    const auto baseState=state;
    const auto fixedReference=target?patch.Reference():source.fixedReference;
    const auto fixedDefect=target?patch.Defect():source.fixedDefect;
    const bool perturbed=packet.get<bool>("perturbed");
    std::vector<double> xt,zt;double amplitude=0.;
    if(perturbed)
    {
        std::set<double> x,z;
        for(const auto &e:patch.fields[0]->Elements())
        {x.insert(e.lower[0]);x.insert(e.upper[0]);z.insert(patch.Physical(e.lower)[2]);z.insert(patch.Physical(e.upper)[2]);}
        std::vector<double> xs(x.begin(),x.end()),zs(z.begin(),z.end());
        Require(xs.size()==15 && (zs.size()==51 || zs.size()==26),"H25 transient knot inventory differs");
        const auto h=std::find(zs.begin(),zs.end(),patch.Physical({0.,0.,.15})[2]);
        Require(h!=zs.end() && h-zs.begin()>=2 && zs.end()-h>=3,"H25 matching knot missing");
        xt.assign(xs.begin()+4,xs.begin()+9);const int stride=zs.size()==51?2:1;
        Require(h-zs.begin()>=2*stride && zs.end()-h>2*stride,"H25 parent-knot neighborhood missing");
        for(int j=-2;j<=2;++j)zt.push_back(*(h+j*stride));
        const double fn=H25Spline(xt,xt[2],0,3),gn=H25Spline(zt,*h,0,3);
        const double gd=H25SplineDerivative(zt,(zt[1]+zt[2])/2);
        Require(fn>0 && gn>0 && gd!=0,"degenerate H25 compact cubic seed");
        amplitude=.001/(fn*gd);const double K=81/(2*std::sqrt(143.))*us*us;
        std::vector<RANSWallEvolutionFunction> functions(patch.fields.size(),[](const RANSWallCoordinate&,double){return 0.;});
        // Project consumes CHART coordinates. The compact physical function
        // and its physical derivatives must be evaluated after this map.
        functions[0]=[=,&patch](const RANSWallCoordinate &chart,double){const auto p=patch.Physical(chart);
            return amplitude*H25Spline(xt,p[0],0,3)*H25SplineDerivative(zt,p[2]);};
        functions[1]=[=,&patch](const RANSWallCoordinate &chart,double){const auto p=patch.Physical(chart);
            return -amplitude*H25SplineDerivative(xt,p[0])*H25Spline(zt,p[2],0,3);};
        functions[2]=[=,&patch](const RANSWallCoordinate &chart,double){const auto p=patch.Physical(chart);
            return .01*K*H25Spline(xt,p[0],0,3)/fn*H25Spline(zt,p[2],0,3)/gn;};
        const double d0=patch.Mapping().displacement,kappa=patch.Parameters().model.Kappa();
        const double physicalH=patch.Physical({0.,0.,patch.Parameters().geometry.height})[2]-d0;
        functions[3]=[=,&patch](const RANSWallCoordinate &chart,double){const auto p=patch.Physical(chart);
            return .01*us*us*us/(kappa*std::max(physicalH,p[2]-d0))*
            H25Spline(xt,p[0],0,3)/fn*H25Spline(zt,p[2],0,3)/gn;};
        const auto increment=patch.Project(functions,0.);
        for(int j=0;j<patch.DifferentialSize();++j)if(increment[j]!=0.)state[j]+=increment[j];
        // This is declared new initial-data construction. Native current-K
        // boundary/trace lifting happens once BEFORE native history production.
        state=patch.Lift(state,0.,data);
        Require(RANSWallR3::Hash(Vector(state.begin()+patch.offset[patch.InnerPressure()],state.end()))==
            RANSWallR3::Hash(Vector(baseState.begin()+patch.offset[patch.InnerPressure()],baseState.end())),
            "transient seed changed stored regional pressure");
        for(const auto &[row,target]:patch.Constraints(baseState,0.,data))
            if(row<patch.offset[patch.Epsilon()])Require(RANSWallR3::Hash(Vector{state[row]})==
                RANSWallR3::Hash(Vector{baseState[row]}),"transient seed changed prescribed velocity/k bits");
    }
    const auto signature=RANSWallPhysicsSignature(patch,serial,0.,false,us,{1.,0.},false);
    std::size_t explicitCalls=0;
    RANSWallTimeIntegration native(1,dt,[](const Vector&,double,double)->Vector{throw std::runtime_error("H25 seed attempted governing solve");},
        [&](const Vector &q,double t){++explicitCalls;auto full=state;std::copy(q.begin(),q.end(),full.begin());return patch.ExplicitBody(full,t,data);});
    native.Initialise(Vector(state.begin(),state.begin()+patch.DifferentialSize()),0.);
    RANSWallTimeIntegration::Checkpoint cp;cp.signature=signature;cp.integration=native.Capture();cp.coupledState=state;
    cp.fixedReference=fixedReference;cp.fixedDefect=fixedDefect;cp.eta=source.eta;
    if(signature!=source.signature)
    {
        bool rejected=false;try{RANSWallTimeIntegration::ReadCheckpoint(imported.cpPath,signature,1,dt,patch.Size(),patch.DifferentialSize());}
        catch(const std::exception &){rejected=true;}
        Require(rejected,"H25 ordinary wrong-signature restart refusal was lost");
    }
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
    return J::Object({{"schema",J::Quote(target?"WJH25TargetNativeSeedV2":"WJH25NativeSeedV2")},{"status",J::Quote("PASS")},{"failure",J::Quote("")},
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
        {"source_signature",J::Quote(source.signature)},{"target_discretization",J::Boolean(target)},
        {"base_state",J::Numbers(baseState)},{"base_state_sha256",J::Quote(RANSWallR3::Hash(baseState))},
        {"ordinary_wrong_signature_restart_rejected",J::Boolean(signature!=source.signature)},
        {"seed_coordinates",J::Quote("physical; chart mapped before evaluation; epsilon scale A/max(physical_H,r)")},
        {"common_original_qh0_sha256",J::Quote(RANSWallR3::Hash(source.coupledState))},{"seed_state_sha256",J::Quote(RANSWallR3::Hash(state))},
        {"perturbed",J::Boolean(perturbed)},{"streamfunction_x_knots",J::Numbers(xt)},{"streamfunction_z_knots",J::Numbers(zt)},
        {"streamfunction_amplitude",J::Number(amplitude)},{"scalar_fraction",J::Number(perturbed?.01:0.)},
        {"velocity_scale",J::Number(perturbed?.001:0.)},{"native_history_initializations","1"},{"explicit_initialization_calls",J::Number(explicitCalls)},
        {"spatial_work",J::Numbers(patch.InitializationWork())},{"initial_velocity_applications",J::Number(patch.InitialVelocityCalls())},{"governing_solves","0"},
        {"construction_scope",J::Quote(target?"strict fine source validation; target physical projection/lift/reference/defect/FlatInletTrace; compact increment/lift; one target native history":"strict fine source; compact increment/lift; one own-dt native history")}});
}
#endif
