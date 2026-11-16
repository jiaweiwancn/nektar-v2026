// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_SHARED_VECTOR_H
#define NEKTAR_RANS_WALL_SHARED_VECTOR_H
#include <memory>
#include <utility>
#include <vector>
namespace Nektar
{
// Ordinary rows retain ordinary vector storage/copy behavior. Share() freezes
// a completed immutable row; subsequent copies share it, and every mutable
// access first detaches. There is no state-dependent numerical payload.
template<class T> class RANSWallSharedVector
{
    std::vector<T> values;
    std::shared_ptr<const std::vector<T>> shared;
    const std::vector<T> &Read() const {return shared?*shared:values;}
    std::vector<T> &Write()
    {if(shared){values=*shared;shared.reset();}return values;}
public:
    RANSWallSharedVector()=default;
    RANSWallSharedVector(std::vector<T> v):values(std::move(v)){}
    RANSWallSharedVector(std::size_t n,const T &v=T{}):values(n,v){}
    RANSWallSharedVector(const RANSWallSharedVector &v)
    {if(v.shared)shared=v.shared;else values=v.values;}
    RANSWallSharedVector(RANSWallSharedVector &&)=default;
    RANSWallSharedVector &operator=(const RANSWallSharedVector &v)
    {if(this!=&v){if(v.shared){values.clear();shared=v.shared;}else{values=v.values;shared.reset();}}return *this;}
    RANSWallSharedVector &operator=(RANSWallSharedVector &&)=default;
    RANSWallSharedVector &operator=(std::vector<T> v){values=std::move(v);shared.reset();return *this;}
    void Share(){if(!shared)shared=std::make_shared<const std::vector<T>>(std::move(values));}
    operator const std::vector<T>&() const {return Read();}
    auto begin() const {return Read().cbegin();}
    auto end() const {return Read().cend();}
    auto begin(){return Write().begin();}
    auto end(){return Write().end();}
    auto size() const {return Read().size();}
    bool empty() const {return Read().empty();}
    const T &operator[](std::size_t i) const {return Read()[i];}
    T &operator[](std::size_t i){return Write()[i];}
    const T &at(std::size_t i) const {return Read().at(i);}
    T &at(std::size_t i){return Write().at(i);}
    const T *data() const {return Read().data();}
    T *data(){return Write().data();}
    const T &front() const {return Read().front();}
    const T &back() const {return Read().back();}
    void push_back(T v){Write().push_back(std::move(v));}
    void reserve(std::size_t n){Write().reserve(n);}
    void resize(std::size_t n){Write().resize(n);}
    void assign(std::size_t n,const T &v){Write().assign(n,v);}
    void clear(){values.clear();shared.reset();}
    friend bool operator==(const RANSWallSharedVector &a,const RANSWallSharedVector &b){return a.Read()==b.Read();}
};
}
#endif
