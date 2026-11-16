// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_TENSOR_LOOKUP_H
#define NEKTAR_RANS_WALL_TENSOR_LOOKUP_H
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <utility>
#include <vector>

namespace Nektar
{
// Candidate index only: callers still apply their original containment, side
// and geometry-ID predicates. No coordinates, rows or state values are cached.
// A nonconforming/overlapping grid uses the existing complete scan instead.
class RANSWallTensorLookup
{
    using Interval=std::pair<double,double>;
    std::array<std::vector<Interval>,2> intervals;
    std::vector<int> cells;
    bool ready=false;
    static std::vector<int> Containing(const std::vector<Interval> &v,double x,double tolerance)
    {
        std::vector<int> result;
        if(!std::isfinite(x))return result;
        // These are exactly the caller's rounded lower-tol / upper+tol
        // comparisons, not the inequivalent x+tol / x-tol rearrangement.
        auto i=std::lower_bound(v.begin(),v.end(),x,[&](const Interval &a,double b){return a.second+tolerance<b;});
        for(;i!=v.end() && x>=i->first-tolerance;++i)
            if(x<=i->second+tolerance)result.push_back(int(i-v.begin()));
        return result;
    }
public:
    template<class Elements> void Build(const Elements &elements,int dimension)
    {
        ready=false;cells.clear();for(auto &v:intervals)v.clear();
        if(dimension!=2 || elements.empty())return;
        for(const auto &e:elements)for(int d=0;d<2;++d)
        {
            const int a=d?2:0;
            if(!std::isfinite(e.lower[a]) || !std::isfinite(e.upper[a]) || !(e.lower[a]<e.upper[a]))return;
            intervals[d].emplace_back(e.lower[a],e.upper[a]);
        }
        for(auto &v:intervals)
        {
            std::sort(v.begin(),v.end());v.erase(std::unique(v.begin(),v.end()),v.end());
            for(std::size_t j=1;j<v.size();++j)if(v[j-1].second>v[j].first)return;
        }
        // Storage is bounded by the actual conforming tensor topology.
        if(intervals[0].size()*intervals[1].size()!=elements.size())return;
        cells.assign(elements.size(),-1);
        for(std::size_t e=0;e<elements.size();++e)
        {
            const auto &v=elements[e];
            const int x=std::lower_bound(intervals[0].begin(),intervals[0].end(),Interval{v.lower[0],v.upper[0]})-intervals[0].begin();
            const int z=std::lower_bound(intervals[1].begin(),intervals[1].end(),Interval{v.lower[2],v.upper[2]})-intervals[1].begin();
            auto &cell=cells[x*intervals[1].size()+z];if(cell>=0)return;cell=int(e);
        }
        ready=std::all_of(cells.begin(),cells.end(),[](int e){return e>=0;});
    }
    bool Ready()const{return ready;}
    std::vector<int> Candidates(const std::array<double,3> &p,double tolerance)const
    {
        std::vector<int> result;
        if(!ready)return result;
        const auto x=Containing(intervals[0],p[0],tolerance),z=Containing(intervals[1],p[2],tolerance);
        for(int i:x)for(int j:z)result.push_back(cells[i*intervals[1].size()+j]);
        std::sort(result.begin(),result.end());return result;
    }
};
}
#endif
