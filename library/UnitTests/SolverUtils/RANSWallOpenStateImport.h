// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_OPEN_STATE_IMPORT_H
#define NEKTAR_RANS_WALL_OPEN_STATE_IMPORT_H
// Shared strict importer extracted from open_checks. No initializer, lift,
// history synthesis or permissive batch-specific checkpoint path.
struct OpenImportedState
{
    std::string reportPath,reportHash,cpPath,checkpointHash;
    Tree native;RANSWallTimeIntegration::Checkpoint cp;Vector history;
    double alpha,dt;std::size_t producerSteps;std::vector<std::string> checkpointHistory;
};
OpenImportedState ImportOpenState(RANSWallEvolutionPatch &patch,
    const LibUtilities::SessionReaderSharedPtr &serial,double us,
    const Tree &packet,RANSWallEvolutionData &data)
{
            const auto reportPath=packet.get<std::string>("native_report");
            const auto reportHash=RANSWallR3::HashFile(reportPath);
            Require(reportHash==packet.get<std::string>("native_report_sha256"),"open producing report SHA256 differs");
            Tree native;boost::property_tree::read_json(reportPath,native);
            Require(native.get<std::string>("failure").empty() || packet.get<bool>("allow_failed_producer",false),
                "open readback of an accepted predecessor in a failed report requires an explicit diagnostic role");
            const auto signature=RANSWallPhysicsSignature(patch,serial,0.,false,us,{1.,0.},false);
            Require(signature==native.get<std::string>("signature"),"open native physical signature differs");
            const auto cpPath=packet.get<std::string>("native_checkpoint");
            const auto checkpointHash=RANSWallR3::HashFile(cpPath+"/State.bin");
            Require(checkpointHash==packet.get<std::string>("checkpoint_sha256"),"open checkpoint bytes differ");
            double dt;serial->LoadParameter("TimeStep",dt);
            const auto cp=RANSWallTimeIntegration::ReadCheckpoint(cpPath,signature,1,dt,patch.Size(),patch.DifferentialSize());
            const auto &state=cp.coupledState;const double time=cp.integration.time;
            if(native.get<std::string>("retention_mode","")=="BoundedV1")
            {
                Require(native.get<std::string>("schema")=="WJCompletedStepV1","bounded strict read requires a genuine completed producer");
                const auto &side=packet.get_child("completed_sidecar"),&cursor=packet.get_child("qualified_cursor");
                RANSWallMonitor::CheckPin(side);RANSWallMonitor::CheckPin(cursor);
                const auto q=RANSWallMonitor::Read(cursor.get<std::string>("path"));const auto &b=q.get_child("binding");
                Require(q.get<std::string>("schema")=="WJQualificationCommitV1","strict bounded cursor was not committed after valid ACK");
                RANSWallMonitor::CheckPin(q.get_child("consumed_ack"));RANSWallMonitor::CheckPin(q.get_child("prepared_cursor"));
                const auto consumed=RANSWallMonitor::Read(q.get<std::string>("consumed_ack.path"));
                RANSWallMonitor::SameObject(b,consumed.get_child("binding"));Require(consumed.get<std::string>("status")=="PASS","strict cursor consumed a failed ACK");
                Require(b.get<std::string>("checkpoint_sha256")==checkpointHash && b.get<std::string>("report_sha256")==reportHash &&
                    b.get<std::string>("callback_sidecar_sha256")==side.get<std::string>("sha256") &&
                    b.get<std::string>("state_sha256")==RANSWallR3::Hash(state) && b.get<std::string>("time_bits")==RANSWallMonitor::Bits(time),
                    "bounded strict qualified cursor/producer differs");
                const auto callback=RANSWallMonitor::Read(side.get<std::string>("path"));
                Require(callback.get<std::string>("schema")=="WJCompletedCallbackV1" && callback.get<std::string>("state_sha256")==RANSWallR3::Hash(state) &&
                    callback.get<std::string>("time_bits")==RANSWallMonitor::Bits(time) && callback.get<double>("alpha_mass")==1/dt &&
                    RANSWallR3::Hash(Read(callback,"callback_history"))==b.get<std::string>("callback_history_sha256"),"bounded strict genuine callback differs");
            }
            const auto producerSteps=native.get<std::size_t>("steps");
            const auto &accepted=native.get_child("accepted_states"),&packets=native.get_child("wj_audit_packets");
            Require(!accepted.empty() && accepted.size()==packets.size() && producerSteps+1>=accepted.size(),"open accepted lineage shape differs");
            const auto first=producerSteps+1-accepted.size();
            Require(first<=cp.integration.steps && cp.integration.steps<=producerSteps,"open checkpoint outside producer accepted lineage");
            auto ai=accepted.begin(),pi=packets.begin();std::advance(ai,cp.integration.steps-first);std::advance(pi,cp.integration.steps-first);
            Vector acceptedState;for(const auto &v:ai->second)acceptedState.push_back(v.second.get_value<double>());
            Require(RANSWallR3::Hash(state)==RANSWallR3::Hash(acceptedState) && time==pi->second.get<double>("time"),"open selected state/time producer differs");
            Require(RANSWallR3::Hash(cp.fixedReference)==RANSWallR3::Hash(Read(native,"fixed_reference")),"open fixed reference differs");
            Require(RANSWallR3::Hash(Vector(state.begin(),state.begin()+patch.DifferentialSize()))==RANSWallR3::Hash(cp.integration.state),
                "open native integrator primary state differs from coupled fields");
            std::vector<std::string> checkpointHistory;
            for(const auto &slot:cp.integration.history)checkpointHistory.push_back(J::Quote(RANSWallR3::Hash(slot)));
            if(cp.integration.steps==producerSteps)
            {
                Require(time==native.get<double>("time") && RANSWallR3::Hash(state)==RANSWallR3::Hash(Read(native,"state")),"open terminal state differs");
                std::size_t slot=0;for(const auto &v:native.get_child("native_history"))
                {Vector h;for(const auto &x:v.second)h.push_back(x.second.get_value<double>());
                    Require(slot<cp.integration.history.size() && RANSWallR3::Hash(h)==RANSWallR3::Hash(cp.integration.history[slot++]),"open native GLM history differs");}
                Require(slot==cp.integration.history.size() && Read(native,"native_history_times")==cp.integration.times,"open native history time labels differ");
            }
            patch.RestoreReference(cp.fixedReference,cp.fixedDefect,cp.eta,data,time);
            Vector history=state;double alpha=0;
            if(cp.integration.steps>0)
            {
                const auto &audit=pi->second;history=Read(audit,"callback_history");alpha=audit.get<double>("alpha_mass");
                Require(history.size()==state.size() && alpha==1/dt && native.get<int>("implicit_callback_count")>0 &&
                    audit.get<std::string>("state_sha256")==RANSWallR3::Hash(state) &&
                    audit.get<std::string>("stage_history_sha256")==RANSWallR3::Hash(history),"open accepted audit/callback binding differs");
                if(cp.integration.steps==producerSteps)Require(native.get<double>("last_stage_time")==time &&
                    native.get<double>("last_stage_alpha")==alpha && RANSWallR3::Hash(Read(native,"last_stage_history"))==RANSWallR3::Hash(history),
                    "open completed callback must come from its producing report, not a zero-step restart");
            }
    Require(RANSWallR3::Hash(cp.fixedDefect)==RANSWallR3::Hash(Read(native,"fixed_defect")),"open fixed defect differs");
    return {reportPath,reportHash,cpPath,checkpointHash,std::move(native),cp,std::move(history),alpha,dt,producerSteps,std::move(checkpointHistory)};
}
#endif
