// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallMappedAlgebra.h>
#include <StdRegions/StdQuadExp.h>
#include <LibUtilities/LinearAlgebra/NekLinSys.hpp>
#include <chrono>
#include <cmath>
#include <iostream>

int main()
{
    using namespace Nektar;
    using namespace Nektar::RANSWallAlgebra;
    using Clock=std::chrono::steady_clock;
    try
    {
        const LibUtilities::BasisKey x(LibUtilities::eModified_A,4,
            LibUtilities::PointsKey(12,LibUtilities::eGaussGaussLegendre));
        const LibUtilities::BasisKey z(LibUtilities::eModified_A,4,
            LibUtilities::PointsKey(16,LibUtilities::eGaussGaussLegendre));
        StdRegions::StdQuadExp expansion(x,z);
        StdRegions::StdMatrixKey key(StdRegions::eMass,expansion.DetShapeType(),expansion);
        const auto full=expansion.GetStdMatrix(key);
        // One affine physical mass, identical coefficient bits/RHS/free rows.
        // x-left essential modes removed; no pressure/scalar equivalence claim.
        std::vector<int> free;for(int j=0;j<16;++j)if(j%4!=0)free.push_back(j);
        const int n=free.size();const double jacobian=2.5*.001/4;
        auto native=MemoryManager<DNekMat>::AllocateSharedPtr(n,n,0.);
        MappedMatrix mapped(n,n,true);mapped.ConsistentMass();
        for(int i=0;i<n;++i)for(int j=0;j<n;++j)
        {const double v=jacobian*(*full)(free[i],free[j]);(*native)(i,j)=v;mapped.Set(i,j,v);}
        const auto nativeStart=Clock::now();DNekLinSys ordinary(native);
        const double nativeSetup=std::chrono::duration<double>(Clock::now()-nativeStart).count();
        const auto mappedStart=Clock::now();MappedBandFactor custom(mapped,true);
        const double mappedSetup=std::chrono::duration<double>(Clock::now()-mappedStart).count();
        Vector truth(n);for(int i=0;i<n;++i)truth[i]=std::sin(.7*i)+.2;
        const auto rhs=mapped.Apply(truth);DNekVec b(n,0.),a(n,0.);
        for(int i=0;i<n;++i)b[i]=rhs[i];ordinary.Solve(b,a);const auto c=custom.Apply(rhs);
        double error=0,residual=0;auto back=mapped.Apply(c);
        for(int i=0;i<n;++i){error=std::max(error,std::abs(a[i]-c[i]));residual=std::max(residual,std::abs(back[i]-rhs[i]));}
        Require(error<=2e-11 && residual<=2e-11,"same frozen mass/native backend parity failed");
        double checksum=0;constexpr int count=2000;
        const auto t0=Clock::now();for(int k=0;k<count;++k){b[k%n]+=.00001;ordinary.Solve(b,a);checksum+=a[k%n];}
        const double nativeSeconds=std::chrono::duration<double>(Clock::now()-t0).count();
        Vector changing=rhs;const auto t1=Clock::now();for(int k=0;k<count;++k){changing[k%n]+=.00001;auto out=custom.Apply(changing);checksum+=out[k%n];}
        const double mappedSeconds=std::chrono::duration<double>(Clock::now()-t1).count();
        std::cout.precision(17);
        std::cout<<"{\"status\":\"PASS\",\"free_rows\":"<<n<<",\"degree\":3,\"quadrature\":[12,16],\"physical_jacobian\":"<<jacobian
            <<",\"max_solution_difference\":"<<error<<",\"max_mapped_residual\":"<<residual
            <<",\"native_factor_seconds\":"<<nativeSetup<<",\"mapped_factor_seconds\":"<<mappedSetup
            <<",\"native_solve_seconds\":"<<nativeSeconds<<",\"mapped_solve_seconds\":"<<mappedSeconds<<",\"applications_each\":"<<count
            <<",\"checksum\":"<<checksum<<",\"scope\":\"identical constrained frozen local mass only; not full stock/WJ equations or a production trajectory\"}\n";
        return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
