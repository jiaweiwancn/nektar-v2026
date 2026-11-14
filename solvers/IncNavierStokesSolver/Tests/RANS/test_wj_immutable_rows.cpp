// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallImmutableRows.h>
#include <vector>
#include <cmath>
#include <cstring>
#include <iostream>
struct Row
{
    std::vector<int> ids;std::vector<double> phi,constant;std::array<std::vector<double>,3>derivative;
};
using C=Nektar::RANSWallImmutableRows<Row>;
static void Check(bool p,const char *s){if(!p)throw std::runtime_error(s);}
static std::uint64_t Bits(double x){std::uint64_t b;std::memcpy(&b,&x,sizeof b);return b;}
static Row Make(double value=1.){Row r;r.ids={0,1,2,3};r.phi={value,2,3,4};r.constant={1,1,0,0};for(auto &v:r.derivative)v={-1,1,2,3};return r;}
int main()
{
    try
    {
        const auto footprint=C::RowBytes(Make());C::Key a{},b{};b[0]=1;
        C c(2,C::fixedBytes+2*(footprint+C::entryBytes));int builds=0;
        auto factory=[&](){++builds;return Make();};auto held=c.Get(a,factory);auto other=c.Get(b,factory);
        Check(c.Get(a,factory)==held && builds==2,"exact hit reconstructed a row");
        C::Key key{};key[0]=2;auto transient=c.Get(key,factory);
        Check(c.Stats()->retainedBypasses==1 && c.Stats()->residentEntries==0 && c.RetainedBytes()==2*footprint,"held evicted rows escaped budget");
        Check(held->phi==Make().phi && transient->phi==Make().phi,"eviction invalidated consumer");
        other.reset();auto admitted=c.Get(key,factory);Check(c.Stats()->residentEntries==1 && c.LiveBytes()<=c.ByteLimit(),"retained release did not permit bounded admission");
        for(int i=0;i<10000;++i){key[0]=i+3;auto r=c.Get(key,factory);Check(r->derivative==Make().derivative && c.LiveBytes()<=c.ByteLimit(),"bounded reconstruction differs");}
        C tiny(1,C::fixedBytes+1);auto bypass=tiny.Get(a,factory);Check(tiny.Stats()->oversized==1 && tiny.Stats()->residentEntries==0,"oversize storage bypass failed");
        C exact(8,65536);std::vector<C::Key>keys;
        for(double x:{0.,-0.,std::nextafter(0.,1.),.15,std::nextafter(.15,1.)}){auto k=a;k[4]=Bits(x);keys.push_back(k);exact.Get(k,factory);}
        auto lower=a,upper=a;lower[1]=17;lower[2]=2;lower[3]=0;upper[1]=18;upper[2]=3;upper[3]=2;
        exact.Get(lower,factory);exact.Get(upper,factory);Check(exact.Stats()->residentEntries==7,"coordinate bits or owner sides merged");
        C patchB(8,65536);auto changed=patchB.Get(keys[0],[]{return Make(7.);});Check(changed->phi!=exact.Get(keys[0],factory)->phi,"different patch/space aliased");
        std::shared_ptr<const Row>survivor;std::shared_ptr<const C::Statistics>stats;
        {C lifetime(1,65536);survivor=lifetime.Get(a,factory);stats=lifetime.Stats();}
        Check(survivor->phi==Make().phi && stats->residentEntries==0 && stats->liveRowBytes==footprint,"catalog destruction lost retained row accounting");
        survivor.reset();Check(stats->liveRowBytes==0,"last shared release was not charged out");
        std::cout<<"exact row catalog PASS; accounted fixture row bytes="<<footprint<<" evictions="<<c.Stats()->evictions<<" peak="<<c.Stats()->peakLiveBytes<<'\n';return 0;
    }catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
