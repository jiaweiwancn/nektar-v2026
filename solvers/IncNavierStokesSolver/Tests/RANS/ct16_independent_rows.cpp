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

// A scope-owned immutable row catalogue. The Python owner keeps every input
// allocation alive and read-only. No state pointer or computed value is stored.
#include <array>
#include <algorithm>
#include <cmath>
#include <memory>
namespace
{
struct RowPlan
{
    std::int64_t n,stateSize;
    int width[3];
    const std::int64_t *ids[3],*pressureOffsets;
    const double *phi[3],*d[3];
    const unsigned char *epsilonPresent;
    std::array<std::int64_t,7> offsets;
};
bool ValidPlanRow(const RowPlan &p,std::int64_t i)
{
    for(int field=0;field<3;++field)
    {
        if(field==1 && !p.epsilonPresent[i])continue;
        const auto first=field==0?p.offsets[0]:field==1?p.offsets[3]:p.pressureOffsets[i];
        if(field==2 && first!=p.offsets[4] && first!=p.offsets[5])return false;
        const auto limit=field==0?std::min({p.offsets[1]-p.offsets[0],p.offsets[2]-p.offsets[1],p.offsets[3]-p.offsets[2]}):
            field==1?p.offsets[4]-p.offsets[3]:first==p.offsets[4]?p.offsets[5]-first:p.offsets[6]-first;
        const int w=p.width[field];
        for(int j=0;j<w;++j)
        {
            if(p.ids[field][i*w+j]<0 || p.ids[field][i*w+j]>=limit ||
                !std::isfinite(p.phi[field][i*w+j]))return false;
            for(int a=0;a<3;++a)if(!std::isfinite(p.d[field][(i*3+a)*w+j]))return false;
        }
    }
    return p.epsilonPresent[i]<=1;
}
}
extern "C" int ct16_plan_abi() { return 1; }
extern "C" void *ct16_plan_create(std::int64_t n,int vWidth,int eWidth,int pWidth,
    const std::int64_t *vIds,const double *vPhi,const double *vD,
    const std::int64_t *eIds,const double *ePhi,const double *eD,const unsigned char *epsilonPresent,
    const std::int64_t *pIds,const double *pPhi,const double *pD,const std::int64_t *pOffsets,
    const std::int64_t *offsets,std::int64_t stateSize)
{
    if(n<1 || n>262144 || stateSize<1 || stateSize>70000 ||
        vWidth<1 || vWidth>128 || eWidth<1 || eWidth>128 || pWidth<1 || pWidth>128)return nullptr;
    if(offsets[0]!=0 || offsets[6]>stateSize)return nullptr;
    for(int b=0;b<6;++b)if(offsets[b]<0 || offsets[b]>=offsets[b+1])return nullptr;
    try
    {
        auto p=std::make_unique<RowPlan>();p->n=n;p->stateSize=stateSize;
        p->width[0]=vWidth;p->width[1]=eWidth;p->width[2]=pWidth;
        p->ids[0]=vIds;p->ids[1]=eIds;p->ids[2]=pIds;
        p->phi[0]=vPhi;p->phi[1]=ePhi;p->phi[2]=pPhi;
        p->d[0]=vD;p->d[1]=eD;p->d[2]=pD;
        p->epsilonPresent=epsilonPresent;p->pressureOffsets=pOffsets;
        std::copy_n(offsets,7,p->offsets.begin());
        for(std::int64_t i=0;i<n;++i)if(!ValidPlanRow(*p,i))return nullptr;
        return p.release();
    }
    catch(...){return nullptr;}
}
extern "C" void ct16_plan_free(void *handle) { delete static_cast<RowPlan*>(handle); }
extern "C" int ct16_plan_point(const void *handle,std::int64_t i,
    const double *state,std::int64_t stateSize,double *output)
{
    const auto p=static_cast<const RowPlan*>(handle);
    if(!p || stateSize!=p->stateSize || i<0 || i>=p->n || !ValidPlanRow(*p,i))return 1;
    for(int b=0;b<5;++b)
    {
        if(b==3 && !p->epsilonPresent[i]){std::fill_n(output+4*b,4,0.);continue;}
        const int field=b<3?0:b==3?1:2,w=p->width[field];
        const auto offset=b<4?p->offsets[b]:p->pressureOffsets[i];
        const int code=ct16_gather_blas_rows(1,w,p->ids[field]+i*w,p->phi[field]+i*w,
            p->d[field]+i*3*w,state,stateSize,offset,output+4*b,output+4*b+1);
        if(code)return code;
    }
    return 0;
}
extern "C" int ct16_plan_grid(const void *handle,std::int64_t count,
    const double *state,std::int64_t stateSize,double *values,double *gradients)
{
    const auto p=static_cast<const RowPlan*>(handle);
    if(!p || count<0 || count>p->n || stateSize!=p->stateSize)return 1;
    for(std::int64_t i=0;i<count;++i)if(!ValidPlanRow(*p,i))return 2;
    std::array<double,20> point;
    for(std::int64_t i=0;i<count;++i)
    {
        const int code=ct16_plan_point(handle,i,state,stateSize,point.data());if(code)return code;
        for(int b=0;b<5;++b)
        {values[5*i+b]=point[4*b];std::copy_n(point.data()+4*b+1,3,gradients+(5*i+b)*3);}
    }
    return 0;
}

// Preserve the original scalar NumPy norm's nine-entry BLAS reduction. Each
// matrix norm is fresh; this helper retains no data or physical classification.
extern "C" int ct16_norm9(std::int64_t n,const double *matrices,double *norms)
{
    if(n<0 || n>262144)return 1;
    for(std::int64_t i=0;i<n;++i)
        norms[i]=std::sqrt(cblas_ddot(9,matrices+9*i,1,matrices+9*i,1));
    return 0;
}
