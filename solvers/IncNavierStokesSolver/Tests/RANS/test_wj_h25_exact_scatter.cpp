// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallMappedAlgebra.h>
#include <array>
#include <cstring>
#include <iostream>

int main()
{
    using namespace Nektar::RANSWallAlgebra;
    try
    {
        for(bool sparse:{false,true})
        {
            MappedMatrix ordinary(37,37,sparse),scatter(37,37,sparse);
            std::vector<double*> slots;
            for(int i=0;i<37;++i)for(int j=0;j<37;++j)slots.push_back(scatter.AssemblyEntry(i,j));
            const std::array<double,9> values{{0.,-0.,1.,-1.,1e-300,-1e-300,1e6,1e-10,-1e6}};
            for(int q=0;q<111;++q)for(int i=0;i<37;++i)for(int j=0;j<37;++j)
            {
                const double x=values[(q+3*i+7*j)%values.size()];ordinary.Add(i,j,x);
                MappedMatrix::AddAssemblyEntry(slots[i*37+j],x);
            }
            for(int i=0;i<37;++i)for(int j=0;j<37;++j)
            {const double a=ordinary.Get(i,j),b=scatter.Get(i,j);Require(std::memcmp(&a,&b,sizeof a)==0,"exact scatter coefficient mismatch");}
            Require(ordinary.Pattern()==scatter.Pattern() && ordinary.Nonzeros()==scatter.Nonzeros(),"scatter zero pattern differs");
            Vector rhs(37);for(int i=0;i<37;++i)rhs[i]=i*.03125;
            Require(ordinary.Apply(rhs)==scatter.Apply(rhs),"scatter action differs");
            slots.clear();ordinary.ClearRow(3);scatter.ClearRow(3);ordinary.Set(3,3,1.);scatter.Set(3,3,1.);
            Require(ordinary.Apply(rhs)==scatter.Apply(rhs),"post-scatter row replacement differs");
        }
        std::cout<<"{\"status\":\"PASS\",\"coefficient_action_pattern_bitwise\":true,\"dense_sparse\":true,\"matrix_instances\":2,\"production_speedup\":\"NOT_MEASURED\"}\n";
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
