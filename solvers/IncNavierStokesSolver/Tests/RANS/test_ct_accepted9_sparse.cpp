// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallSerialSparse.h>
#include <cmath>
#include <iostream>
using namespace Nektar::RANSWallAlgebra;
int main()
{
    try
    {
        SparseMatrix m(3,3);m.Entry(0,0)=4;m.Entry(0,1)=1;m.Entry(1,0)=1;m.Entry(1,1)=3;m.Entry(1,2)=1;m.Entry(2,1)=1;m.Entry(2,2)=2;m.Finalize();
        const Vector expected{2,-1,3},rhs{7,2,5};const auto x=SparseFactor(m,true).Apply(rhs);
        for(int i=0;i<3;++i)Require(std::abs(x[i]-expected[i])<2e-13,"consistent mass inverse differs from independent solution");
        const auto copy=m;Require(copy.NonZeros()==7 && copy.Get(0,2)==0,"CSR canonical storage");
        const auto restricted=m.Restrict({2,0});const auto v=restricted.Apply({2,3});Require(v==Vector({4,12}),"signed-ID restriction/action ordering");
        SparseMatrix a(3,3);a.Entry(0,0)=2;a.Entry(0,1)=-1;a.Entry(1,1)=3;a.Entry(1,2)=1;a.Entry(2,0)=1;a.Entry(2,2)=2;a.Finalize();
        const auto y=SparseFactor(a).Apply({2,-1,4.5});for(int i=0;i<3;++i)Require(std::abs(y[i]-Vector({.5,-1,2})[i])<2e-13,"nonsymmetric sparse LU solution");
        const auto p=SparseScotchPermutation(a);Require(p==SparseScotchPermutation(a),"Scotch ordering is nondeterministic");
        const auto sy=SparseFactor(a,false,SparsePredictorOrdering::ScotchND0V1).Apply({2,-1,4.5});
        for(int i=0;i<3;++i)Require(std::abs(sy[i]-Vector({.5,-1,2})[i])<2e-13,"Scotch canonical solve differs from independent solution");
        SparseMatrix vector(48,48);Vector known(48);std::vector<int> groups(48);
        for(int i=0;i<48;++i)
        {
            known[i]=(i%2?-1.:1.)/(1+i);groups[i]=i%24;vector.Entry(i,i)=5.;
            if(i%24)vector.Entry(i,i-1)=-1.;
            if(i%24<23)vector.Entry(i,i+1)=-.7;
            vector.Entry(i,(i+24)%48)=.125;
        }
        vector.Finalize();const auto vp=SparseScotchPermutation(vector,groups);
        Require(vp==SparseScotchPermutation(vector,groups) && vp!=groups,"vector permutation control has no rearrangement");
        const auto vx=SparseFactor(vector,false,SparsePredictorOrdering::ScotchVectorNDKeepV1,groups).Apply(vector.Apply(known));
        for(int i=0;i<48;++i)Require(std::abs(vx[i]-known[i])<2e-13,"vector ND solve differs from independent manufactured values");
        bool bad=false;try{m.Entry(0,0)=0.;}catch(const std::exception&){bad=true;}Require(bad,"finalized matrix must be immutable");
        bad=false;try{SparseMatrix invalid(-1,2);}catch(const std::exception&){bad=true;}Require(bad,"invalid sparse shape accepted");
        bad=false;try{SparseMatrix singular(2,2);singular.Entry(0,0)=1;singular.Finalize();SparseFactor invalid(singular);}catch(const std::exception&){bad=true;}Require(bad,"singular sparse factor accepted");
        Require(std::abs(rhs[0]/m.Get(0,0)-x[0])>.1,"mass-lumping negative control is nondiscriminating");
        m.Write("ct9-algebra.csr");std::cout<<"CT9 CSR, direct mass/LU, restriction and negative controls pass\n";return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
