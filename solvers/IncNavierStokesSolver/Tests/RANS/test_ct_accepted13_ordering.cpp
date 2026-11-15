// SPDX-License-Identifier: MIT
// Read-only use of the actual guarded native factor and refinement code.
#include <IncNavierStokesSolver/EquationSystems/RANSWallSerialSparse.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <fstream>
#include <iostream>
#include <cstdint>
#include <cstring>
using namespace Nektar::RANSWallAlgebra;
namespace J=Nektar::RANSWallSolveJson;
int main(int argc,char **argv)
{
    if(argc!=5 && argc!=6)return 1;
    try
    {
        std::ifstream in(argv[1],std::ios::binary);char tag[8];in.read(tag,8);
        Require(std::memcmp(tag,"RANSCSR1",8)==0,"wrong canonical CSR encoding");
        std::uint64_t rows,cols,nnz;in.read(reinterpret_cast<char*>(&rows),8);in.read(reinterpret_cast<char*>(&cols),8);in.read(reinterpret_cast<char*>(&nnz),8);
        Require(rows>0 && rows==cols && rows<=70000 && nnz<=20000000,"CSR capacity/shape");
        std::vector<std::uint64_t> ptr(rows+1);std::vector<std::uint32_t> ids(nnz);Vector val(nnz);
        in.read(reinterpret_cast<char*>(ptr.data()),ptr.size()*8);in.read(reinterpret_cast<char*>(ids.data()),ids.size()*4);
        in.read(reinterpret_cast<char*>(val.data()),val.size()*8);Require(in.good(),"truncated CSR");
        Require(ptr[0]==0 && ptr.back()==nnz,"invalid CSR pointers");
        SparseMatrix matrix(rows,cols);
        for(int i=0;i<rows;++i)for(auto j=ptr[i];j<ptr[i+1];++j)
        {Require(ids[j]<cols && val[j]!=0.,"invalid/explicit-zero CSR entry");matrix.Entry(i,ids[j])=val[j];}
        matrix.Finalize();Require(matrix.NonZeros()==nnz,"CSR canonical entries changed");
        Vector rhs(rows);std::ifstream r(argv[2]);for(auto &v:rhs)r>>v;Require(!r.fail(),"truncated RHS");
        const std::string name=argv[3];Require(name=="COLAMD" || name=="AMD" || name=="Natural" || name=="ScotchND0V1" || (name=="ScotchVectorND0V1" || name=="ScotchVectorNDKeepV1"),"unknown predictor ordering");
        const auto ordering=name=="ScotchVectorNDKeepV1"?SparsePredictorOrdering::ScotchVectorNDKeepV1:name=="ScotchVectorND0V1"?SparsePredictorOrdering::ScotchVectorND0V1:name=="ScotchND0V1"?SparsePredictorOrdering::ScotchND0V1:name=="AMD"?SparsePredictorOrdering::AMD:name=="Natural"?SparsePredictorOrdering::Natural:SparsePredictorOrdering::COLAMD;
        std::vector<int> groups;if((name=="ScotchVectorND0V1" || name=="ScotchVectorNDKeepV1"))
        {Require(argc==6,"missing canonical vector groups");groups.resize(rows);std::ifstream g(argv[5]);for(auto &v:groups)g>>v;Require(!g.fail(),"truncated canonical vector groups");}
        std::vector<int> permutation;
        if(name=="ScotchND0V1" || (name=="ScotchVectorND0V1" || name=="ScotchVectorNDKeepV1"))permutation=SparseScotchPermutation(matrix,groups);
        const SparseFactor factor(matrix,false,ordering,groups);const auto x=factor.Apply(rhs);
        auto residual=matrix.Apply(x);for(int i=0;i<rows;++i)residual[i]-=rhs[i];
        // Tiny signed-cancellation right-hand side uses the same factor and
        // unchanged native absolute/relative true-residual criterion.
        auto tiny=rhs;for(auto &v:tiny)v*=1e-12;
        const auto tx=factor.Apply(tiny);auto tr=matrix.Apply(tx);for(int i=0;i<rows;++i)tr[i]-=tiny[i];
        Vector cancelling(rows);for(int i=0;i<rows;++i)cancelling[i]=(i%2?1.:-1.)*1e-14/(1+i%7);
        const auto crhs=matrix.Apply(cancelling),cx=factor.Apply(crhs);auto cres=matrix.Apply(cx);
        for(int i=0;i<rows;++i)cres[i]-=crhs[i];
        std::ofstream out(argv[4]);out<<J::Object({{"status",J::Quote("pass")},{"ordering",J::Quote(name)},
            {"rows",J::Number(rows)},{"nnz",J::Number(nnz)},{"factor_nnz",J::Number(factor.FactorNonZeros())},
            {"rhs_maximum",J::Number(Maximum(rhs))},{"true_residual",J::Number(Maximum(residual))},
            {"tiny_rhs_maximum",J::Number(Maximum(tiny))},{"tiny_true_residual",J::Number(Maximum(tr))},
            {"cancelling_rhs_maximum",J::Number(Maximum(crhs))},{"cancelling_true_residual",J::Number(Maximum(cres))},
            {"old_to_new",J::Numbers(permutation)},
            {"solution",J::Numbers(x)},{"tiny_solution",J::Numbers(tx)}})<<'\n';return 0;
    }
    catch(const std::exception &e)
    {std::ofstream out(argv[4]);out<<J::Object({{"status",J::Quote("blocked_or_failed")},{"failure",J::Quote(e.what())}})<<'\n';std::cerr<<e.what()<<'\n';return 2;}
}
