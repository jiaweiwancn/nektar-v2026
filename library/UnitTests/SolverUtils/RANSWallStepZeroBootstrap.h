// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_STEP_ZERO_BOOTSTRAP_H
#define NEKTAR_RANS_WALL_STEP_ZERO_BOOTSTRAP_H

// A separate producer, never an ordinary restart or a checkpoint-header edit.
// The caller strictly imports the source under its OWN signature and dt first.
RANSWallTimeIntegration::Checkpoint StepZeroIMEX1(
    const RANSWallTimeIntegration::Checkpoint &source,const std::string &signature,
    double dt,const RANSWallTimeIntegration::ExplicitCallback &rhs,
    std::size_t &explicitCalls)
{
    Require(source.integration.order==1 && source.integration.steps==0 &&
        RANSWallMonitor::Bits(source.integration.time)==RANSWallMonitor::Bits(0.),
        "step-zero bootstrap requires genuine IMEX1 absolute step/time zero");
    Require(signature!=source.signature && dt>0 && dt!=source.integration.dt,
        "step-zero bootstrap requires its own dt-bearing identity");
    RANSWallTimeIntegration native(1,dt,
        [](const Vector &,double,double)->Vector
        {throw std::runtime_error("step-zero bootstrap attempted a governing implicit solve");},
        [&](const Vector &q,double t){++explicitCalls;return rhs(q,t);});
    native.Initialise(source.integration.state,0.);
    RANSWallTimeIntegration::Checkpoint coarse;coarse.signature=signature;coarse.integration=native.Capture();
    coarse.coupledState=source.coupledState;coarse.fixedReference=source.fixedReference;
    coarse.fixedDefect=source.fixedDefect;coarse.eta=source.eta;
    Require(RANSWallR3::Hash(coarse.integration.state)==RANSWallR3::Hash(source.integration.state) &&
        coarse.integration.steps==0 && coarse.integration.time==0 && coarse.integration.dt==dt,
        "native step-zero history construction changed its physical primary state");
    return coarse;
}

std::string StepZeroBootstrap(RANSWallEvolutionPatch &patch,
    const LibUtilities::SessionReaderSharedPtr &serial,double us,
    const Tree &packet,RANSWallEvolutionData &data,bool fixture)
{
    double dt;serial->LoadParameter("TimeStep",dt);
    Require(patch.WJ() && patch.OpenBoundary() && patch.InitialVelocityMode()=="FlatInletTraceV1",
        "step-zero helper requires the compatible unforced WJ open control");
    const auto coarseSignature=RANSWallPhysicsSignature(patch,serial,0.,false,us,{1.,0.},false);
    auto sourceDt=packet.get<double>("source_dt");
    RANSWallMonitor::CheckPin(packet.get_child("source_xml"));
    Require(sourceDt>0 && dt!=sourceDt,"bootstrap source dt must differ");
    // The signature still includes TIMESTEP. No generic validation is relaxed.
    // All current physical/mesh/Function data must also satisfy the original
    // source signature; the external XML allowlist separately pins source bytes.
    serial->SetParameter("TimeStep",sourceDt);
    OpenImportedState imported;
    try{imported=ImportOpenState(patch,serial,us,packet,data);}
    catch(...){serial->SetParameter("TimeStep",dt);throw;}
    serial->SetParameter("TimeStep",dt);
    const auto &source=imported.cp;Require(source.eta==0.,"step-zero bootstrap requires unchanged eta0");
    std::size_t explicitCalls=0;
    auto rhs=[&](const Vector &q,double t)
    {auto full=source.coupledState;std::copy(q.begin(),q.end(),full.begin());return patch.ExplicitBody(full,t,data);};
    const auto coarse=StepZeroIMEX1(source,coarseSignature,dt,rhs,explicitCalls);
    auto reject=[&](auto call){bool bad=false;try{call();}catch(const std::exception &){bad=true;}
        Require(bad,"step-zero/wrong-dt refusal was lost");};
    reject([&]{RANSWallTimeIntegration::ReadCheckpoint(imported.cpPath,coarseSignature,1,dt,patch.Size(),patch.DifferentialSize());});
    if(fixture)
    {
        Require(patch.Size()<=3000,"bootstrap component fixture is small");
        auto bad=source;bad.integration.steps=1;
        reject([&]{std::size_t n=0;StepZeroIMEX1(bad,coarseSignature,dt,rhs,n);});
        bad=source;bad.integration.time=sourceDt;
        reject([&]{std::size_t n=0;StepZeroIMEX1(bad,coarseSignature,dt,rhs,n);});
        bad=source;bad.integration.order=2;
        reject([&]{std::size_t n=0;StepZeroIMEX1(bad,coarseSignature,dt,rhs,n);});
        RANSWallTimeIntegration restored(1,dt,[](const Vector&,double,double)->Vector{throw std::runtime_error("unexpected solve");},rhs);
        reject([&]{restored.Restore(source.integration);});
    }
    Require(patch.InitializationWork()==Vector({0.,0.}) && patch.InitialVelocityCalls()==0,
        "bootstrap constructed new spatial fields");
    const auto audit=fixture?RANSWallEvolutionPatch::Diagnostic{}:
        patch.CompleteRestrictedAudit(coarse.coupledState,coarse.coupledState,0.,0.,data);
    std::string bound="{}";
    if(!fixture)
    {
        const auto k=patch.ApplyPimpleKBound(coarse.coupledState,coarse.coupledState,0.,data);
        bound=k.report;Require(k.success && !k.activated && RANSWallR3::Hash(k.state)==RANSWallR3::Hash(coarse.coupledState),"coarse CP0 bound changed state");
        Require(patch.ScalarDomain(coarse.coupledState,0.,data,true).admissible,"coarse CP0 scalar domain failed");
        patch.ValidateOutlet(coarse.coupledState,"step-zero bootstrap");
        Require(audit.at("wj_hard_admissible")[0]==1,"coarse CP0 native WJ hard audit failed");
    }
    const std::string checkpoint=fixture?"":RANSWallTimeIntegration::WriteCheckpoint(packet.get<std::string>("checkpoint_base"),coarse);
    if(!fixture){const int fd=::open((checkpoint+"/State.bin").c_str(),O_RDONLY);Require(fd>=0 && ::fsync(fd)==0,"bootstrap checkpoint fsync failed");::close(fd);RANSWallMonitor::CheckpointDirectories(checkpoint);}
    std::vector<std::string> history;for(const auto &h:coarse.integration.history)history.push_back(J::Numbers(h));
    std::vector<std::pair<std::string,std::string>> a;for(const auto &[n,v]:audit)a.emplace_back(n,J::Numbers(v));
    const auto packetText=fixture?"{}":patch.WJAuditPacket(audit,coarse.coupledState,coarse.coupledState,0.,0.,"startup");
    const std::string result=J::Object({{"schema",J::Quote("WJStepZeroIMEX1BootstrapV1")},{"status",J::Quote("PASS")},
        {"failure",J::Quote("")},{"steps","0"},{"time","0"},{"dt",J::Number(dt)},{"signature",J::Quote(coarseSignature)},
        {"state",J::Numbers(coarse.coupledState)},{"fixed_reference",J::Numbers(coarse.fixedReference)},{"fixed_defect",J::Numbers(coarse.fixedDefect)},
        {"eta",J::Number(coarse.eta)},{"native_history",J::List(history)},{"native_history_times",J::Numbers(coarse.integration.times)},
        {"accepted_states",J::List({J::Numbers(coarse.coupledState)})},{"wj_audit_packets",J::List({packetText})},
        {"native_hard_audit",J::Object(a)},{"initial_k_bound_check",bound},{"native_full_audit_status",J::Quote(fixture?"NOT_RUN":"PASS")},
        {"state_sha256",J::Quote(RANSWallR3::Hash(coarse.coupledState))},{"checkpoint",J::Quote(checkpoint)},
        {"source_report_sha256",J::Quote(imported.reportHash)},{"source_checkpoint_sha256",J::Quote(imported.checkpointHash)},
        {"source_signature",J::Quote(source.signature)},{"source_dt",J::Number(sourceDt)},
        {"native_history_initializations","1"},{"explicit_initialization_calls",J::Number(explicitCalls)},
        {"initialization_work",J::Numbers(patch.InitializationWork())},{"initial_velocity_calls","0"},{"governing_solves","0"},
        {"ordinary_wrong_dt_restart_rejected","true"},{"fixture",J::Boolean(fixture)}});
    Require(RANSWallR3::HashFile(imported.cpPath+"/State.bin")==imported.checkpointHash && RANSWallR3::HashFile(imported.reportPath)==imported.reportHash,
        "bootstrap changed immutable source bytes");
    return result;
}
#endif
