// SPDX-License-Identifier: MIT
// A bounded read-only prototype; not linked into the governing runtime.
#include <IncNavierStokesSolver/EquationSystems/RANSWallSharedVector.h>
#include <array>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
using Nektar::RANSWallSharedVector;
namespace
{
struct Row
{
    RANSWallSharedVector<int> ids;
    RANSWallSharedVector<double> phi,constant;
    std::array<RANSWallSharedVector<double>,3> derivative;
    int offset=0;
};
using State=std::vector<double>;
void require(bool value,const char *message){if(!value)throw std::runtime_error(message);}
// Literal same-order current point Value then Gradient arithmetic.
__attribute__((noinline)) std::array<double,4> original(const Row &r,const State &state)
{
    const auto &ids=r.ids;const auto &constant=r.constant;const auto &phi=r.phi;const auto &derivative=r.derivative;const int offset=r.offset;
    require(offset>=0 && !ids.empty() && constant.size()==ids.size() && phi.size()==ids.size() && constant[0]!=0,"invalid value row");
    const double valueCentre=state.at(offset+ids[0])/constant[0];long double value=valueCentre;
    for(std::size_t j=0;j<ids.size();++j)value+=(long double)phi[j]*(state.at(offset+ids[j])-valueCentre*constant[j]);
    require(offset>=0 && !ids.empty() && constant.size()==ids.size() && constant[0]!=0,"invalid gradient row");
    const double centre=state.at(offset+ids[0])/constant[0];std::array<long double,3> gradient{};
    for(int d=0;d<3;++d)
    {
        require(derivative[d].size()==ids.size(),"invalid derivative row");
        for(std::size_t j=0;j<ids.size();++j)gradient[d]+=(long double)derivative[d][j]*(state.at(offset+ids[j])-centre*constant[j]);
    }
    return {static_cast<double>(value),static_cast<double>(gradient[0]),static_cast<double>(gradient[1]),static_cast<double>(gradient[2])};
}
__attribute__((noinline)) std::array<double,4> packed(const Row &r,const State &state)
{
    const auto n=r.ids.size();if(n>128)return original(r,state);
    require(r.offset>=0 && n && r.constant.size()==n && r.phi.size()==n && r.constant[0]!=0,"invalid value row");
    for(int d=0;d<3;++d)require(r.derivative[d].size()==n,"invalid derivative row");
    const auto ids=r.ids.data();const auto constant=r.constant.data();const auto phi=r.phi.data();
    const double centre=state.at(r.offset+ids[0])/constant[0];std::array<double,128> delta;
    for(std::size_t j=0;j<n;++j)delta[j]=state.at(r.offset+ids[j])-centre*constant[j];
    long double value=centre;for(std::size_t j=0;j<n;++j)value+=(long double)phi[j]*delta[j];
    std::array<long double,3> gradient{};
    for(int d=0;d<3;++d){const auto row=r.derivative[d].data();for(std::size_t j=0;j<n;++j)gradient[d]+=(long double)row[j]*delta[j];}
    return {static_cast<double>(value),static_cast<double>(gradient[0]),static_cast<double>(gradient[1]),static_cast<double>(gradient[2])};
}
thread_local std::string error;
}
extern "C" const char *ct16_contraction_error(){return error.c_str();}
extern "C" void *ct16_contraction_make(std::size_t n,const std::size_t *ptr,const int *ids,const double *phi,const double *constant,
    const double *dx,const double *dy,const double *dz,const int *offsets)
{
    try
    {
        require(n<=65536,"prototype row bound");auto rows=std::make_unique<std::vector<Row>>(n);
        for(std::size_t i=0;i<n;++i)
        {
            require(ptr[i+1]>ptr[i] && ptr[i+1]-ptr[i]<=128,"prototype row shape bound");auto &r=rows->at(i);const auto a=ptr[i],b=ptr[i+1];
            r.ids=std::vector<int>(ids+a,ids+b);r.phi=std::vector<double>(phi+a,phi+b);r.constant=std::vector<double>(constant+a,constant+b);r.offset=offsets[i];
            for(int d=0;d<3;++d)r.derivative[d]=std::vector<double>((d==0?dx:d==1?dy:dz)+a,(d==0?dx:d==1?dy:dz)+b);
            r.ids.Share();r.phi.Share();r.constant.Share();for(auto &d:r.derivative)d.Share();
        }
        return rows.release();
    }
    catch(const std::exception &e){error=e.what();return nullptr;}
}
extern "C" int ct16_contraction_run(void *handle,std::size_t size,const double *data,int repeats,int policy,double *out)
{
    try
    {
        require(handle && size<=70000 && repeats>0 && repeats<=32,"prototype runtime bound");
        const State state(data,data+size);const auto &rows=*static_cast<std::vector<Row>*>(handle);
        for(int repeat=0;repeat<repeats;++repeat)for(std::size_t i=0;i<rows.size();++i)
        {
            const auto result=policy?packed(rows[i],state):original(rows[i],state);
            std::memcpy(out+4*i,result.data(),4*sizeof(double));
        }
        return 0;
    }
    catch(const std::exception &e){error=e.what();return 1;}
}
extern "C" void ct16_contraction_free(void *handle){delete static_cast<std::vector<Row>*>(handle);}
