// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallMonitor.h>
#include <iostream>
namespace M=Nektar::RANSWallMonitor;
namespace J=Nektar::RANSWallSolveJson;
namespace D=Nektar::RANSWallR3;
int main(int argc,char **argv)
{
    try
    {
        D::Check(argc==2,"monitor checks require an isolated output path");const std::string root=argv[1];
        D::Check(!std::filesystem::exists(root),"fixture output exists");std::filesystem::create_directories(root);
        M::FactorCatalog c;c.Configure(root+"/catalog","fixed-mesh-signs-source");
        for(int i=0;i<1000;++i)c.Record("changed-numerical-zero-pattern-"+std::to_string(i),"{}");
        for(int i=0;i<100;++i)c.Record("changed-numerical-zero-pattern-999","{}");
        auto summary=c.Summary();M::Tree s;std::istringstream input(summary);boost::property_tree::read_json(input,s);
        D::Check(s.get<int>("resident_definitions")==16 && s.get<int>("evictions")==984 && s.get<int>("definitions")==1000 && s.get<int>("calls")==1100,"factor retention is unbounded or suppressed a definition/call");
        std::vector<std::pair<std::string,std::string>> fields;
        for(const auto &name:{"nonce","absolute_step","time_bits","alpha_bits","state_sha256","callback_history_sha256","checkpoint_sha256",
            "callback_sidecar_sha256","report_sha256","original_qh0_sha256","fixed_reference_sha256","fixed_defect_sha256","physical_signature_sha256",
            "structure_sha256","source_sha256","observer_sha256","inventory_rule_sha256"})fields.push_back({name,J::Quote(name)});
        const auto binding=J::Object(fields);
        auto check=[&](const std::string &name,const std::string &fault,bool pass)
        {
            const auto dir=root+"/"+name;std::filesystem::create_directories(dir);const auto ready=dir+"/step-1.ready.json";
            M::Atomic(ready,J::Object({{"binding",binding}}));const auto readyHash=D::HashFile(ready);
            M::Atomic(dir+"/inventory.json","{}");const auto inventoryHash=D::HashFile(dir+"/inventory.json");
            M::Atomic(dir+"/result.json",J::Object({{"status",J::Quote("PASS")},{"binding",binding},{"inventory_sha256",J::Quote(inventoryHash)}}));
            M::Atomic(dir+"/cursor.json",J::Object({{"binding",binding},{"ready_sha256",J::Quote(readyHash)},{"inventory_sha256",J::Quote(inventoryHash)}}));
            std::string actual=binding;
            if(fault=="binding")actual.replace(actual.find("\"nonce\":\"nonce\""),15,"\"nonce\":\"stale\"");
            if(fault!="timeout" && fault!="temporary-only")M::Atomic(dir+"/step-1.ack.json",J::Object({{"schema",J::Quote("WJObserverAckV1")},{"status",J::Quote(fault=="FAIL"?"FAIL":"PASS")},
                {"ready_sha256",J::Quote(readyHash)},{"binding",actual},{"result",D::Pin(dir+"/result.json")},{"qualified_cursor",D::Pin(dir+"/cursor.json")},
                {"inventory",D::Pin(dir+"/inventory.json")},{"continuation",J::Quote("TERMINAL")}}));
            if(fault=="duplicate")M::Atomic(dir+"/step-0.ack.json","{}");
            if(fault=="publishing")D::Check(::link((dir+"/step-1.ack.json").c_str(),(dir+"/step-1.ack.json.tmp-publishing").c_str())==0,"cannot represent atomic dual-name window");
            if(fault=="temporary-only")M::Atomic(dir+"/step-1.ack.json.tmp-publishing","unpublished/incomplete");
            bool success=false;try{success=M::Wait(dir,"step-1",readyHash,binding,.06,getpid())=="TERMINAL";}catch(const std::exception&){}
            D::Check(success==pass,"incorrect protocol decision: "+name);
        };
        check("pass","",true);check("atomic-dual-name","publishing",true);check("failed","FAIL",false);check("stale","binding",false);check("duplicate","duplicate",false);check("timeout","timeout",false);check("temporary-only","temporary-only",false);
        bool overwritten=false;try{M::Atomic(root+"/pass/inventory.json","replacement");overwritten=true;}catch(const std::exception&){}
        D::Check(!overwritten,"immutable output overwritten");std::cout<<"bounded catalog and protocol checks PASS "<<summary<<'\n';return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
