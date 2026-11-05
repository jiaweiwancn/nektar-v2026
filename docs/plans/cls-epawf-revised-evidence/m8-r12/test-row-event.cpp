// R12 metadata regression; no constitutive or PDE calculation.
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
int main(int argc,char **argv)
{
    namespace D=Nektar::RANSWallR3;namespace J=Nektar::RANSWallSolveJson;
    D::Check(argc==2,"test-row-event PROGRESS");
    D::DurableProgress progress(argv[1],1,400.,.0025,{1.},"VerifiedR3Action:v1");
    progress.Event("r12_original_cut",{{"id","0"},{"family",J::Quote("continuous_wall_reserve")},
        {"identity",J::Quote(D::Hash(std::vector<double>{1.,2.}))},
        {"row",progress.VectorPayload("row",{1.,2.})}});
}
