// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_SHARED_INVENTORY_H
#define NEKTAR_RANS_WALL_SHARED_INVENTORY_H
#include <memory>
#include <vector>
namespace Nektar
{
// Immutable reads share row/observation storage. The few construction and
// append sites detach before writing; copies retain order and membership.
template<class T> class RANSWallSharedInventory
{
    std::shared_ptr<std::vector<T>> values=std::make_shared<std::vector<T>>();
    void Detach(){if(!values.unique())values=std::make_shared<std::vector<T>>(*values);}
public:
    auto begin() const {return values->cbegin();}
    auto end() const {return values->cend();}
    std::size_t size() const {return values->size();}
    bool empty() const {return values->empty();}
    const T &front() const {return values->front();}
    const T &operator[](std::size_t i) const {return (*values)[i];}
    void push_back(T value){Detach();values->push_back(std::move(value));}
    void clear(){values=std::make_shared<std::vector<T>>();}
    void reserve(std::size_t n){Detach();values->reserve(n);}
    void shrink_to_fit(){Detach();values->shrink_to_fit();}
};
}
#endif
