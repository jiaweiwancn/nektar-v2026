// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_IMMUTABLE_ROWS_H
#define NEKTAR_RANS_WALL_IMMUTABLE_ROWS_H
#include <algorithm>
#include <array>
#include <cstdint>
#include <list>
#include <map>
#include <memory>
#include <stdexcept>
#include <type_traits>

namespace Nektar
{
/** Serial patch-generation-owned exact rows. Accounting is a conservative
 * allocation upper bound: vector capacities +32 per buffer, object allocation
 * +32, shared control/deleter +96, map/list/key allocation +384 per entry, and
 * fixed catalog/control overhead1024. No bucket array or unbounded key trace.
 * Evicted shared rows remain in liveRowBytes until the last consumer releases
 * them; an exhausted retained budget bypasses storage instead of raising caps.
 */
template<class Row> class RANSWallImmutableRows
{
public:
    using Key=std::array<std::uint64_t,7>;
    struct Statistics
    {
        std::uint64_t requests=0,hits=0,misses=0,constructs=0,evictions=0,bypasses=0,
            oversized=0,retainedBypasses=0,residentEntries=0,residentRowBytes=0,
            liveRowBytes=0,peakLiveBytes=0,peakEntries=0,maximumRequestGap=0,
            minimumRowBytes=~std::uint64_t(0),maximumRowBytes=0;
        std::array<std::uint64_t,32> requestGapHistogram{};
    };
    static constexpr std::uint64_t fixedBytes=1024,entryBytes=384;
    RANSWallImmutableRows(std::size_t entries,std::uint64_t bytes):m_entries(entries),m_bytes(bytes)
    {
        if(!entries || entries>524288 || bytes<fixedBytes || bytes>512ULL*1024*1024)
            throw std::runtime_error("invalid immutable row cache limits");
        m_stats->peakLiveBytes=fixedBytes;
    }
    ~RANSWallImmutableRows(){while(!m_lru.empty())Evict();}
    RANSWallImmutableRows(const RANSWallImmutableRows&)=delete;
    static std::uint64_t RowBytes(const Row &r)
    {
        std::uint64_t n=sizeof(Row)+32+96;
        const auto add=[&](const auto &v){n+=v.capacity()*sizeof(typename std::decay_t<decltype(v)>::value_type)+32;};
        add(r.ids);add(r.phi);add(r.constant);for(const auto &v:r.derivative)add(v);return n;
    }
    template<class Factory> std::shared_ptr<const Row> Get(const Key &key,Factory factory)
    {
        const auto tick=++m_stats->requests;auto it=m_rows.find(key);
        if(it!=m_rows.end())
        {
            ++m_stats->hits;auto gap=tick-it->second.tick;m_stats->maximumRequestGap=std::max(m_stats->maximumRequestGap,gap);
            unsigned bin=0;for(auto v=gap;v>1 && bin<31;v>>=1)++bin;++m_stats->requestGapHistogram[bin];
            it->second.tick=tick;m_lru.splice(m_lru.begin(),m_lru,it->second.position);return it->second.row;
        }
        ++m_stats->misses;++m_stats->constructs;Row row=factory();const auto bytes=RowBytes(row);
        m_stats->minimumRowBytes=std::min(m_stats->minimumRowBytes,bytes);m_stats->maximumRowBytes=std::max(m_stats->maximumRowBytes,bytes);
        if(bytes+entryBytes+fixedBytes>m_bytes)
        {++m_stats->oversized;++m_stats->bypasses;return std::make_shared<const Row>(std::move(row));}
        while(!m_lru.empty() && (m_rows.size()>=m_entries || LiveBytes()+bytes+entryBytes>m_bytes))Evict();
        if(LiveBytes()+bytes+entryBytes>m_bytes)
        {++m_stats->retainedBypasses;++m_stats->bypasses;return std::make_shared<const Row>(std::move(row));}
        auto stats=m_stats;
        const auto raw=new Row(std::move(row));m_stats->liveRowBytes+=bytes;
        auto result=std::shared_ptr<const Row>(raw,[stats,bytes](const Row *p){delete p;stats->liveRowBytes-=bytes;});
        m_lru.push_front(key);
        try{m_rows.emplace(key,Entry{result,m_lru.begin(),tick,bytes});}
        catch(...){m_lru.pop_front();throw;}
        ++m_stats->residentEntries;m_stats->residentRowBytes+=bytes;
        m_stats->peakEntries=std::max(m_stats->peakEntries,m_stats->residentEntries);
        m_stats->peakLiveBytes=std::max(m_stats->peakLiveBytes,LiveBytes());return result;
    }
    std::uint64_t LiveBytes()const{return fixedBytes+m_stats->liveRowBytes+entryBytes*m_stats->residentEntries;}
    std::uint64_t ResidentBytes()const{return fixedBytes+m_stats->residentRowBytes+entryBytes*m_stats->residentEntries;}
    std::uint64_t RetainedBytes()const{return m_stats->liveRowBytes-m_stats->residentRowBytes;}
    std::shared_ptr<const Statistics> Stats()const{return m_stats;}
    std::size_t EntryLimit()const{return m_entries;}
    std::uint64_t ByteLimit()const{return m_bytes;}
private:
    struct Entry{std::shared_ptr<const Row> row;typename std::list<Key>::iterator position;std::uint64_t tick,bytes;};
    void Evict()
    {
        auto it=m_rows.find(m_lru.back());m_stats->residentRowBytes-=it->second.bytes;--m_stats->residentEntries;
        m_rows.erase(it);m_lru.pop_back();++m_stats->evictions;
    }
    const std::size_t m_entries;const std::uint64_t m_bytes;
    std::shared_ptr<Statistics> m_stats=std::make_shared<Statistics>();
    std::map<Key,Entry> m_rows;std::list<Key> m_lru;
};
}
#endif
