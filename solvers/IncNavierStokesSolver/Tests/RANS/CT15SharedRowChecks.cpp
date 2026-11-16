// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallSharedVector.h>
#include <cstring>
#include <iostream>
#include <stdexcept>
using Nektar::RANSWallSharedVector;
int main()
{
    const auto check=[](bool ok){if(!ok)throw std::runtime_error("shared row ownership changed");};
    RANSWallSharedVector<double> ordinary(std::vector<double>{-0.,1.,3.});
    auto deep=ordinary;const auto &a=ordinary,&b=deep;
    check(a.data()!=b.data() && std::memcmp(a.data(),b.data(),3*sizeof(double))==0);
    ordinary.Share();auto frozen=ordinary;const auto &f=frozen;
    check(a.data()==f.data());frozen[1]=4.;check(a[1]==1. && f[1]==4. && a.data()!=f.data());
    auto writable=ordinary;for(auto &x:writable)x+=1.;check(a[1]==1. && writable[1]==2.);
    auto assigned=writable;assigned=ordinary;check(static_cast<const decltype(assigned)&>(assigned).data()==a.data());
    assigned.assign(4,2.);check(a.size()==3 && assigned.size()==4);
    auto moved=std::move(frozen);check(moved[1]==4. && a[1]==1.);
    auto zeros=RANSWallSharedVector<double>(25,0.);zeros.Share();auto z=zeros;
    z[0]=-0.;check(std::memcmp(static_cast<const decltype(zeros)&>(zeros).data(),static_cast<const decltype(z)&>(z).data(),sizeof(double))!=0);
    const std::vector<double> saved=a;check(std::memcmp(saved.data(),a.data(),3*sizeof(double))==0);
    std::cout<<"shared row deep-copy, frozen-read, mutation-detach, assignment, move and signed-zero checks pass\n";
}
