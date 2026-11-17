// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeInitialFields.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>

namespace J=Nektar::RANSWallSolveJson;
using namespace Nektar::WJNative;
namespace
{
void Check(bool condition,const std::string&reason){if(!condition)throw std::runtime_error(reason);}
}
int main(int argc,char**argv)
{
    try{
        if(argc!=3)throw std::invalid_argument("usage: initial_fields_probe packet.json new-output.json");
        if(std::filesystem::exists(argv[2]))throw std::invalid_argument("immutable probe output exists");
        InitialFields fields(argv[1]);const auto&m=fields.Metadata();const auto&g=fields.GeometryData();
        fields.ValidateBindings(m.sourceIdentity,m.referenceIdentity,m.model,m.manufactured,m.mmsVariant,g);
        std::vector<std::string>rows,rejections;
        for(const auto&e:fields.Elements()){
            const std::array<std::array<double,2>,7>points{{{0,0},{1,0},{0,1},{1,1},{.17,.63},{.79,.26},{.5,.5}}};
            for(auto point:points){const double x=e.xl+(e.xr-e.xl)*point[0],r=e.rl+(e.rr-e.rl)*point[1],z=r+g.d0;
                std::vector<std::string>owned,automatic;
                for(unsigned f=0;f<5;++f){
                    owned.push_back(f==3&&e.epsilon.empty()?"null":J::Number(fields.Evaluate(f,x,z,e.geometry)));
                    const double tol=64*std::numeric_limits<double>::epsilon()*std::max({g.H,std::abs(r),std::abs(g.H)});
                    automatic.push_back(f==3&&r<g.H-tol?"null":J::Number(fields.Evaluate(f,x,z)));
                }
                rows.push_back(J::Object({{"owner",J::Number(e.geometry)},{"x",J::Number(x)},{"physical_z",J::Number(z)},
                    {"owned",J::List(owned)},{"automatic",J::List(automatic)}}));
            }
        }
        auto reject=[&](const std::string&name,const auto&operation){
            try{operation();}catch(const std::invalid_argument&error){rejections.push_back(J::Object({{"name",J::Quote(name)},{"reason",J::Quote(error.what())}}));return;}
            throw std::runtime_error("invalid initial-field operation was accepted: "+name);
        };
        const double x=(g.left+g.right)/2;
        reject("inner_epsilon",[&]{fields.Evaluate(3,x,g.d0+g.H/2);});
        reject("outside_domain",[&]{fields.Evaluate(0,g.right+g.H,g.d0+g.H);});
        reject("invalid_field",[&]{fields.Evaluate(5,x,g.d0+g.H);});
        reject("absent_owner",[&]{fields.Evaluate(0,x,g.d0+g.H,std::numeric_limits<int>::max());});
        reject("nonfinite_coordinate",[&]{fields.Evaluate(0,std::numeric_limits<double>::quiet_NaN(),g.d0+g.H);});
        bool innerOwner=false;
        for(const auto&e:fields.Elements())if(e.epsilon.empty()&&std::abs(e.rr-g.H)<=64*std::numeric_limits<double>::epsilon()*g.H){
            reject("explicit_inner_owner_at_H",[&]{fields.Evaluate(3,(e.xl+e.xr)/2,g.d0+g.H,e.geometry);});innerOwner=true;break;}
        Check(innerOwner,"probe packet lacks an inner H owner");
        Check(fields.Elements().size()>1,"probe needs multiple native owners");
        const auto&second=fields.Elements()[1];
        reject("wrong_explicit_owner",[&]{fields.Evaluate(0,(second.xl+second.xr)/2,g.d0+(second.rl+second.rr)/2,fields.Elements()[0].geometry);});
        reject("source_binding",[&]{fields.ValidateBindings(m.sourceIdentity+"x",m.referenceIdentity,m.model,m.manufactured,m.mmsVariant,g);});
        reject("reference_binding",[&]{fields.ValidateBindings(m.sourceIdentity,m.referenceIdentity+"x",m.model,m.manufactured,m.mmsVariant,g);});
        reject("model_binding",[&]{fields.ValidateBindings(m.sourceIdentity,m.referenceIdentity,m.model+"x",m.manufactured,m.mmsVariant,g);});
        reject("manufactured_binding",[&]{fields.ValidateBindings(m.sourceIdentity,m.referenceIdentity,m.model,!m.manufactured,m.mmsVariant,g);});
        reject("variant_binding",[&]{fields.ValidateBindings(m.sourceIdentity,m.referenceIdentity,m.model,m.manufactured,m.mmsVariant+"x",g);});
        reject("geometry_binding",[&]{auto wrong=g;wrong.top+=g.H;fields.ValidateBindings(m.sourceIdentity,m.referenceIdentity,m.model,m.manufactured,m.mmsVariant,wrong);});
        Check(fields.SameInventory(fields.Elements(),fields.Sizes()),"identical inventory failed");
        auto reordered=fields.Elements();std::reverse(reordered.begin(),reordered.end());
        Check(fields.SameInventory(reordered,fields.Sizes()),"element ordering changed inventory equality");
        auto changed=fields.Elements();changed[0].us[0]*=-1;
        Check(!fields.SameInventory(changed,fields.Sizes()),"changed sign accepted as same inventory");
        changed=fields.Elements();std::swap(changed[0].u[0],changed[0].u[1]);
        Check(!fields.SameInventory(changed,fields.Sizes()),"changed map accepted as same inventory");
        changed=fields.Elements();changed[0].xr+=1e-6*g.H;
        Check(!fields.SameInventory(changed,fields.Sizes()),"changed bounds accepted as same inventory");
        auto sizes=fields.Sizes();++sizes[2];Check(!fields.SameInventory(fields.Elements(),sizes),"changed size accepted as same inventory");
        const std::string output=J::Object({{"schema",J::Quote("WJH27InitialFieldsProbeV1")},{"status",J::Quote("RECORDED")},
            {"packet_sha256",J::Quote(m.packetSha256)},{"state_sha256",J::Quote(m.stateSha256)},
            {"role",J::Quote(m.role==InitialFieldRole::QualifiedReference?"QualifiedReference":"Candidate")},
            {"qualified_reference",J::Boolean(m.qualifiedReference)},{"sizes",J::Numbers(fields.Sizes())},
            {"same_inventory_controls",J::Quote("identity/reordering accepted; signs/maps/bounds/sizes changes rejected")},
            {"samples",J::List(rows)},{"rejections",J::List(rejections)}});
        std::ofstream file(argv[2]);file<<output<<'\n';file.flush();Check(bool(file),"initial fields probe output write failed");file.close();Check(!file.fail(),"initial fields probe output close failed");
        std::cout<<"Recorded "<<rows.size()<<" owner/value samples and "<<rejections.size()<<" rejection controls\n";
    }catch(const std::exception&error){std::cerr<<error.what()<<'\n';return 2;}
    return 0;
}
