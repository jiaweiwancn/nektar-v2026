// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_COMPLETE_AUDIT_BATCH_H
#define NEKTAR_RANS_WALL_COMPLETE_AUDIT_BATCH_H
#include <chrono>
std::string ImmutableRowChecks(RANSWallEvolutionPatch &patch,
    const LibUtilities::SessionReaderSharedPtr &serial,double us,const Tree &packet,RANSWallEvolutionData &data)
{
    Require(patch.Size()<=3000,"immutable row fixture is small/read-only");
    auto imported=ImportOpenState(patch,serial,us,packet,data);
    const auto checks=patch.EvaluationReuseChecks(imported.cp.coupledState,imported.history,
        imported.cp.integration.time,data,true);
    // Native block/space and separate patch/map tests. No reference projection,
    // initializer or governing solve is performed in these row-only patches.
    const auto el=*std::find_if(patch.fields[0]->Elements().begin(),patch.fields[0]->Elements().end(),
        [&](const auto &e){return e.lower[2]>=patch.Parameters().geometry.height;});
    RANSWallCoordinate q{(el.lower[0]+el.upper[0])/2,0,(el.lower[2]+el.upper[2])/2};
    const auto original=patch.GeometryRow(0,el.geometry,q),pressure=patch.GeometryRow(patch.OuterPressure(),el.geometry,q);
    Require(original->phi.size()!=pressure->phi.size(),"different native spaces aliased");
    const auto graph=SpatialDomains::MeshGraphIO::Read(serial);
    {
        RANSWallEvolutionPatch other(serial,graph,serial->GetComm(),patch.Parameters(),patch.Mapping());
        const auto row=other.GeometryRow(0,el.geometry,q);
        Require(row.get()!=original.get() && row->ids==original->ids && row->phi==original->phi && row->constant==original->constant && row->derivative==original->derivative,
            "different native patches aliased or changed identical rows");
    }
    {
        auto map=patch.Mapping();map.displacement+=.001;
        RANSWallEvolutionPatch other(serial,graph,serial->GetComm(),patch.Parameters(),map);
        const auto row=other.GeometryRow(0,el.geometry,q);
        Require(row.get()!=original.get() && row->phi==original->phi && row->derivative!=original->derivative,
            "synthetic different physical map inherited a transformed row");
    }
    Require(patch.InitializationWork()==Vector({0.,0.}) && patch.InitialVelocityCalls()==0,"read-only row fixture initialized fields");
    return J::Object({{"status",J::Quote("PASS")},{"checks",checks},{"separate_patch_space_map","true"},
        {"scope",J::Quote("synthetic/read-only native rows; altered displacement tests identity only, never a case/height candidate")},
        {"full_audit_calls",J::Number(patch.CompleteWJCalls())},{"governing_solves","0"},{"physical_commits","0"}});
}

std::string CompleteAuditBatch(RANSWallEvolutionPatch &patch,
    const LibUtilities::SessionReaderSharedPtr &serial,double us,const Tree &packet,
    RANSWallEvolutionData &data,const std::string &output)
{
    Require(patch.Size()==24688 && patch.CompleteWJCalls()==0,"complete target batch scope/call count");
    const auto &requests=packet.get_child("states");Require(requests.size()==6,"complete batch has exactly six saved roles");
    const std::array<std::size_t,6>sequence{0,1,2,3,4,0};std::vector<std::string>results;
    const auto start=std::chrono::steady_clock::now();const double maximum=packet.get<double>("maximum_seconds");
    Require(maximum>0 && maximum<=900,"bounded complete batch deadline");
    std::string referenceHash,defectHash,qh0Hash;
    const auto progress=[&](const std::string &status,const std::string &reason="")
    {RANSWallMonitor::Atomic(output+".progress.json",J::Object({{"schema",J::Quote("WJH17AuditProgressV1")},
        {"status",J::Quote(status)},{"failure",J::Quote(reason)},{"completed_evaluations",J::Number(results.size())},
        {"full_audit_calls",J::Number(patch.CompleteWJCalls())},{"governing_solves","0"}}),false);};
    progress("STARTED");
    try
    {
        std::size_t n=0;
        for(const auto &request:requests)
        {
            const auto began=std::chrono::steady_clock::now();const auto beforeCalls=patch.CompleteWJCalls();
            data.Touch();auto imported=ImportOpenState(patch,serial,us,request.second,data);
            Require(patch.CompleteWJCalls()==beforeCalls,"strict import unexpectedly evaluated complete WJ audit");
            const auto &cp=imported.cp;Require(cp.integration.steps==sequence[n],"wrong absolute step in complete audit sequence");
            const auto stateHash=RANSWallR3::Hash(cp.coupledState);
            if(!n)
            {
                qh0Hash=stateHash;referenceHash=RANSWallR3::Hash(cp.fixedReference);defectHash=RANSWallR3::Hash(cp.fixedDefect);
                Require(qh0Hash=="631caa32a59ff43253f1fcac4effd955eb68b8529d9c28f29e615cf14ca0f077","original target qh0 differs");
            }
            Require(referenceHash==RANSWallR3::Hash(cp.fixedReference) && defectHash==RANSWallR3::Hash(cp.fixedDefect),"fixed reference/defect changed across saved states");
            if(!sequence[n])Require(imported.alpha==0 && stateHash==qh0Hash,"initial state/history role differs");
            else Require(imported.alpha==400,"completed audit alpha differs");
            const auto before=patch.ImmutableRowsJSON();const auto importedAt=std::chrono::steady_clock::now();
            const auto audit=patch.CompleteAuditCommitment(cp.coupledState,imported.history,imported.alpha,cp.integration.time,data);
            const auto afterAudit=std::chrono::steady_clock::now();
            const auto constraints=patch.Constraints(cp.coupledState,cp.integration.time,data);
            std::vector<double>rows,values;for(auto [r,v]:constraints){rows.push_back(r);values.push_back(v);}
            const auto spatial=patch.Residual(cp.coupledState,cp.coupledState,0.,cp.integration.time,data,false);
            const auto stage=imported.alpha>0?patch.Residual(cp.coupledState,imported.history,imported.alpha,cp.integration.time,data,false):Vector{};
            const auto exterior=patch.ExteriorResidual(cp.coupledState,cp.integration.time,data);
            const auto boundary=patch.OpenBoundaryDiagnostic(cp.coupledState,imported.history,imported.alpha,cp.integration.time,data);
            const auto roughness=patch.RoughnessBudget(cp.coupledState,cp.integration.time,data);
            const auto energy=patch.EnergyBudget(cp.coupledState,imported.history,imported.alpha,cp.integration.time,data);
            Require(patch.CompleteWJCalls()==beforeCalls+1 && patch.InitializationWork()==Vector({0.,0.}) && patch.InitialVelocityCalls()==0,
                "saved audit performed extra audit/initialization work");
            Require(imported.reportHash==RANSWallR3::HashFile(imported.reportPath) && imported.checkpointHash==RANSWallR3::HashFile(imported.cpPath+"/State.bin"),"saved audit modified imported bytes");
            const auto elapsed=std::chrono::duration<double>(std::chrono::steady_clock::now()-began).count();
            std::vector<std::string> nativeHistory;
            for(const auto &slot:cp.integration.history)nativeHistory.push_back(J::Quote(RANSWallR3::Hash(slot)));
            results.push_back(J::Object({{"step",J::Number(sequence[n])},{"time",J::Number(cp.integration.time)},
                {"alpha",J::Number(imported.alpha)},{"initial_history_role",J::Quote(sequence[n]?"genuine_completed_callback":"initial_alpha0_no_completed_callback")},
                {"state_sha256",J::Quote(stateHash)},{"history_sha256",J::Quote(RANSWallR3::Hash(imported.history))},
                {"producer_sha256",J::Quote(imported.reportHash)},{"checkpoint_file_sha256",J::Quote(imported.checkpointHash)},
                {"integrator_primary_sha256",J::Quote(RANSWallR3::Hash(cp.integration.state))},{"native_history_sha256",J::List(nativeHistory)},
                {"native_history_times",J::Numbers(cp.integration.times)},{"reference_sha256",J::Quote(referenceHash)},
                {"defect_sha256",J::Quote(defectHash)},{"qh0_sha256",J::Quote(qh0Hash)},{"complete",audit},
                {"constraint_rows_sha256",J::Quote(RANSWallR3::Hash(rows))},{"constraint_values_sha256",J::Quote(RANSWallR3::Hash(values))},
                {"spatial_residual_sha256",J::Quote(RANSWallR3::Hash(spatial))},{"completed_residual_sha256",J::Quote(RANSWallR3::Hash(stage))},
                {"exterior_load_sha256",J::Quote(RANSWallR3::Hash(exterior))},{"boundary_exact_native_json_sha256",J::Quote(RANSWallR3::Hash(boundary))},
                {"boundary_comparison",J::Quote("exact native max_digits10 serialization, plus direct IEEE exterior-vector hash; no rounded observer stream")},
                {"roughness_budget_sha256",J::Quote(RANSWallR3::Hash(roughness))},{"energy_budget_sha256",J::Quote(RANSWallR3::Hash(energy))},
                {"import_seconds",J::Number(std::chrono::duration<double>(importedAt-began).count())},
                {"audit_seconds",J::Number(std::chrono::duration<double>(afterAudit-importedAt).count())},
                {"operation_seconds",J::Number(elapsed)},{"rows_before",before},{"rows_after",patch.ImmutableRowsJSON()}}));
            ++n;progress("RUNNING");
            const double spent=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
            if(n==1)Require(spent+1.25*elapsed*5<maximum,"RESOURCE_LIMIT complete batch first-pilot remaining projection");
            Require(spent<maximum,"RESOURCE_LIMIT complete batch deadline");
        }
        Require(patch.CompleteWJCalls()==6,"complete batch total differs");progress("COMPLETE");
    }
    catch(const std::exception &e){progress("FAIL",e.what());throw;}
    return J::Object({{"schema",J::Quote("WJH17CompleteBatchV1")},{"status",J::Quote("PASS")},
        {"results",J::List(results)},{"full_audit_calls",J::Number(patch.CompleteWJCalls())},
        {"immutable_rows",patch.ImmutableRowsJSON()},{"dynamic_reuse",patch.EvaluationReuseJSON()},
        {"dynamic_scope",J::Quote("no outer evaluation scope; every imported state/data generation recomputes dynamic results; rows alone persist")},
        {"governing_solves","0"},{"physical_commits","0"}});
}
#endif
