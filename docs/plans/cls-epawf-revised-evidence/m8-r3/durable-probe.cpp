// SPDX-License-Identifier: MIT
// Metadata/I/O fixture only; no constitutive or PDE calculation.
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <iostream>

int main(int argc,char **argv)
{
    namespace D=Nektar::RANSWallR3;namespace J=Nektar::RANSWallSolveJson;
    try
    {
        D::Check(argc==3,"usage: durable-probe PROGRESS complete|hold");
        D::Check(D::Hash(std::string{})=="e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855","empty SHA256 vector");
        D::Check(D::Hash(std::string("abc"))=="ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad","abc SHA256 vector");
        D::Check(D::Hash(std::string(1000000,'a'))=="cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0","long SHA256 vector");
        const std::vector<double> history{1,2,3},state{1.1,2.1,3.1},residual{.1,.2,.3},norms{.3};
        D::DurableProgress progress(argv[1],1,400.,.0025,history,"VerifiedR3Action:v1");
        progress.Event("protocol",{{"tolerance","2e-12"},{"maximum_Newton_updates","12"},
            {"benchmark_only","false"},{"no_auxiliary_alpha_sequence","true"}});
        progress.Snapshot(0,state,residual,{{0,1.1}},norms,"{\"admissible\":true}","{}");
        progress.JacobianStatus(0,"incomplete");
        if(std::string(argv[2])=="hold")for(;;)::pause();
        const auto matrix=progress.MatrixPayload(0,{{1,0,0},{0,1,0},{0,0,1}});
        progress.JacobianStatus(0,"complete",matrix,"{}");
        progress.Event("solve_terminal",{{"converged","false"},
            {"terminal",J::Quote("durability_fixture_not_a_PDE_solve")}});
        return 0;
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 1;}
}
