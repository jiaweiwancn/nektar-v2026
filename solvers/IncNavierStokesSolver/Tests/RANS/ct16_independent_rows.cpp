// SPDX-License-Identifier: MIT
// Independent Python audit rows: preserve the same BLAS calls and their order.
#include <cblas.h>
#include <cstdint>
#include <vector>

extern "C" int ct16_row_abi() { return 1; }
extern "C" int ct16_gather_blas_rows(std::int64_t n,int width,
    const std::int64_t *ids,const double *phi,const double *derivative,
    const double *state,std::int64_t stateSize,std::int64_t offset,
    double *values,double *gradients)
{
    if(n<0 || n>262144 || width<1 || width>128 || stateSize<0 || offset<0 || offset>stateSize)return 1;
    // Validate all IDs before writing any result. The Python owner supplies
    // strongly held, contiguous, immutable arrays with checked full shapes.
    for(std::int64_t i=0;i<n*width;++i)
        if(ids[i]<0 || ids[i]>=stateSize-offset)return 2;
    try
    {
        std::vector<double> current(width);
        for(std::int64_t i=0;i<n;++i)
        {
            for(int j=0;j<width;++j)current[j]=state[offset+ids[i*width+j]];
            values[i]=cblas_ddot(width,phi+i*width,1,current.data(),1);
            cblas_dgemv(CblasRowMajor,CblasNoTrans,3,width,1.,
                derivative+i*3*width,width,current.data(),1,0.,gradients+i*3,1);
        }
    }
    catch(...){return 3;}
    return 0;
}
