// SPDX-License-Identifier: MIT
// Read-only exact flat-geometry primitive experiment; no governing linkage.
#include <SolverUtils/Turbulence/RANSWallNormalMap.h>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>
using namespace Nektar::SolverUtils;
namespace {thread_local std::string error;}
extern "C" const char *ct16_flat_metric_error(){return error.c_str();}
extern "C" int ct16_flat_metric_run(std::size_t n,const double *coordinates,double height,
    const double *bounds,int repeats,int policy,double *out)
{
    try
    {
        if(n>65536 || repeats<1 || repeats>32)throw std::invalid_argument("prototype bound");
        RANSWallNormalMap map;map.height=height;RANSVector lower{bounds[0],0,bounds[1]},upper{bounds[2],0,bounds[3]};
        map.Validate(lower,upper,.1);if(!map.Identity())throw std::invalid_argument("not identity flat geometry");
        const auto seed=map.At(lower);
        for(int repeat=0;repeat<repeats;++repeat)for(std::size_t i=0;i<n;++i)
        {
            const RANSVector p{coordinates[3*i],coordinates[3*i+1],coordinates[3*i+2]};
            if(policy)
            {
                if(!std::all_of(p.begin(),p.end(),[](double x){return std::isfinite(x);}) || p[1]!=0 || p[0]<lower[0] || p[0]>upper[0] || p[2]<lower[2] || p[2]>upper[2])
                    throw std::invalid_argument("outside certified exact flat chart");
                out[3*i]=seed.height;out[3*i+1]=seed.scale;out[3*i+2]=0.;
            }
            else {out[3*i]=map.At(p).height;out[3*i+1]=map.At(p).scale;out[3*i+2]=!map.Identity();}
        }
        return 0;
    }
    catch(const std::exception &e){error=e.what();return 1;}
}
