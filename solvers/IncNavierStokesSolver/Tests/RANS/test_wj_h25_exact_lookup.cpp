// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallTensorLookup.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
struct Element { int geometry;std::array<double,3> lower,upper; };
int Owner(const std::vector<Element> &e,const std::array<double,3> &p,int side,
          const Nektar::RANSWallTensorLookup *index)
{
    for(double tolerance:{0.,1e-11})
    {
        std::vector<int> candidates;
        if(index)candidates=index->Candidates(p,tolerance);
        else for(int i=0;i<int(e.size());++i)candidates.push_back(i);
        int selected=-1;
        for(int i:candidates)
        {
            bool contains=true;
            for(int a:{0,2})contains=contains && p[a]>=e[i].lower[a]-tolerance && p[a]<=e[i].upper[a]+tolerance;
            const double radial=tolerance?1e-12:0.;
            if(!contains || (side>0 && p[2]>=e[i].upper[2]-radial) || (side<0 && p[2]<=e[i].lower[2]+radial))continue;
            if(selected<0 || e[i].geometry<e[selected].geometry)selected=i;
        }
        if(selected>=0)return selected;
    }
    return -1;
}
void Require(bool value,const char *why){if(!value)throw std::runtime_error(why);}
}
int main()
{
    try
    {
        std::vector<Element> e;
        const std::vector<double> xs{-10.,-7.5,-0.,.25,1.,15.},zs{0.,.001,.1,.15,.150000000001,1.,20.};
        for(int x=0;x<int(xs.size())-1;++x)for(int z=0;z<int(zs.size())-1;++z)
            e.push_back({700-17*int(e.size()),{xs[x],0.,zs[z]},{xs[x+1],0.,zs[z+1]}});
        std::reverse(e.begin(),e.end());
        Nektar::RANSWallTensorLookup index;index.Build(e,2);Require(index.Ready(),"conforming index missing");
        std::size_t queries=0;
        std::vector<double> x=xs,z=zs;
        for(double q:xs)for(double d:{-1e-11,-1e-12,1e-12,1e-11})x.push_back(q+d);
        for(double q:zs)for(double d:{-1e-11,-1e-12,1e-12,1e-11})z.push_back(q+d);
        for(double q:xs){x.push_back(std::nextafter(q,-INFINITY));x.push_back(std::nextafter(q,INFINITY));}
        for(double q:zs){z.push_back(std::nextafter(q,-INFINITY));z.push_back(std::nextafter(q,INFINITY));}
        x.push_back(+0.);z.push_back(-0.);
        for(double a:x)for(double b:z)for(int side:{-1,0,1})
        {const std::array<double,3> p{a,0.,b};Require(Owner(e,p,side,&index)==Owner(e,p,side,nullptr),"face/ULP/side owner differs");++queries;}
        for(int i=0;i<10000;++i)
        {const std::array<double,3> p{-10.+25.*((i*2377)%10007)/10007.,0.,20.*((i*541)%10009)/10009.};
         for(int side:{-1,0,1}){Require(Owner(e,p,side,&index)==Owner(e,p,side,nullptr),"moving coordinate owner differs");++queries;}}
        for(double bad:{NAN,INFINITY,-INFINITY})for(int a:{0,2})
        {std::array<double,3> p{.2,0.,.11};p[a]=bad;Require(Owner(e,p,0,&index)==Owner(e,p,0,nullptr),"invalid coordinate failure differs");}
        auto overlap=e;overlap[0].lower[2]=0.;index.Build(overlap,2);Require(!index.Ready(),"overlapping partition did not fall back");
        index.Build(e,3);Require(!index.Ready(),"unsupported3D index accepted");
        e.pop_back();index.Build(e,2);Require(!index.Ready(),"missing cell did not fall back");
        std::cout<<"{\"status\":\"PASS\",\"exact_owner_comparisons\":"<<queries
            <<",\"covers\":\"physical faces, radial sides, native-ID ordering, ULP, signed zero, moving points, rejection/fallback\"}\n";
        return 0;
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
}
