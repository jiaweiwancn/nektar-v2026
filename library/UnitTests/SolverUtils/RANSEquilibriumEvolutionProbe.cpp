///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <LibUtilities/Communication/CommSerial.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>

using namespace Nektar;
using Vector=RANSWallAlgebra::Vector;
using WallMatrix=RANSWallAlgebra::DenseMatrix;
using RANSWallAlgebra::Require;
namespace
{
void Write(std::ostream &o,const Vector &v){o<<'[';for(int i=0;i<v.size();++i){if(i)o<<',';o<<v[i];}o<<']';}
void Rows(std::ostream &o,const WallMatrix &v){o<<'[';for(int i=0;i<v.size();++i){if(i)o<<',';Write(o,v[i]);}o<<']';}
void Describe(std::ostream &o,const std::shared_ptr<RANSWallGlobalField> &field)
{
    o<<"{\"size\":"<<field->Size()<<",\"mesh\":[";const auto map=field->Field()->GetLocalToGlobalMap();
    for(int e=0;e<field->Elements().size();++e)
    {
        if(e)o<<',';
        const auto &g=field->Elements()[e];const int nc=field->Field()->GetExp(e)->GetNcoeffs(),start=field->Field()->GetCoeff_Offset(e);
        o<<"{\"geometry\":"<<g.geometry<<",\"degree\":"<<g.degree<<",\"lower\":";Write(o,Vector(g.lower.begin(),g.lower.end()));
        o<<",\"upper\":";Write(o,Vector(g.upper.begin(),g.upper.end()));o<<",\"ids\":[";
        for(int j=0;j<nc;++j){if(j)o<<',';o<<map->GetLocalToGlobalMap(start+j);}o<<"],\"sign\":[";
        for(int j=0;j<nc;++j){if(j)o<<',';o<<map->GetLocalToGlobalSign(start+j);}o<<"]}";
    }
    o<<"]}";
}
void Generic(std::ostream &out,const LibUtilities::SessionReaderSharedPtr &serial,const SpatialDomains::MeshGraphSharedPtr &graph,
    const LibUtilities::CommSharedPtr &comm,RANSWallPatchParameters parameters,double chi)
{
    parameters.chi=chi;RANSWallEvolutionPatch patch(serial,graph,comm,parameters);const int d=patch.Dimension();
    const double K=parameters.model.Q()*.16,H=parameters.geometry.height,Z=patch.upper[2];
    using F=RANSWallEvolutionFunction;std::vector<F> f;
    f.push_back([](const RANSWallCoordinate &p,double){return (.1+.03*p[0]+.02*p[1])*p[2];});
    if(d==3)f.push_back([](const RANSWallCoordinate &p,double){return .04*p[2]+.01*p[0]*p[1]*p[2];});
    f.push_back([Z](const RANSWallCoordinate &p,double){return .02*p[2]+.01*(p[0]+p[1])*p[2]*(Z-p[2]);});
    f.push_back([K,rK=parameters.geometry.sampleHeight](const RANSWallCoordinate &p,double){return K+.02*p[0]+.03*p[1]-.02*(p[2]-rK);});
    f.push_back([K,H,chi,q=parameters.model.Q(),kappa=parameters.model.Kappa()](const RANSWallCoordinate &p,double){
        const double sampled=K+.02*p[0]+.03*p[1];return (1-chi)*std::pow(sampled,1.5)/(2.55*H)+chi*std::pow(sampled/q,1.5)/(kappa*H);});
    f.push_back([](const RANSWallCoordinate &p,double){return -.01*p[0]+.015*p[1]+.02*p[2];});
    f.push_back([](const RANSWallCoordinate &p,double){return .04*p[0]+.03*p[1]+.05*p[2];});
    RANSWallEvolutionData data(d);for(int b=0;b<patch.InnerPressure();++b)data.boundary[b]=f[b];
    data.innerBody[0]=[](const RANSWallCoordinate &p,double){return .01+.02*p[0];};
    data.outerBody[0]=[](const RANSWallCoordinate &p,double){return .04*p[0]+.01*p[2];};
    if(d==3)
    {
        data.innerBody[1]=[](const RANSWallCoordinate &p,double){return .015-.01*p[1]+.02*p[2];};
        data.outerBody[1]=[](const RANSWallCoordinate &p,double){return .02*p[1]-.01*p[2];};
    }
    data.innerBody[d-1]=[](const RANSWallCoordinate &p,double){return -.02+.03*p[2];};
    data.outerBody[d-1]=[](const RANSWallCoordinate &p,double){return -.03+.02*p[2];};
    data.innerKSource=[](const RANSWallCoordinate &p,double){return .03+.01*p[0]+.04*p[1]+.02*p[2];};
    data.outerKSource=[](const RANSWallCoordinate &p,double){return -.02+.01*p[0]+.02*p[1]+.03*p[2];};
    data.outerEpsilonSource=[](const RANSWallCoordinate &p,double){return .01*(1+p[0]+p[1]+p[2]);};
    auto state=patch.Project(f,.2),history=state;for(int j=0;j<patch.DifferentialSize();++j)history[j]*=.97;
    const auto constraints=patch.Constraints(state,.2,data);
    const auto lifted=patch.Lift(state,.2,data);
    const auto residual=patch.Residual(state,history,5.3,.2,data,false);
    out<<"{\"chi\":"<<chi<<",\"alpha\":5.3,\"time\":0.2,\"state\":";Write(out,state);
    out<<",\"history\":";Write(out,history);out<<",\"lifted_state\":";Write(out,lifted);out<<",\"residual\":";Write(out,residual);
    out<<",\"constraints\":[";int count=0;for(auto [row,value]:constraints){if(count++)out<<',';Write(out,{static_cast<double>(row),value});}
    out<<"],\"forcing\":[";count=0;
    for(double x:{.13,.47,.91})
    {
        RANSWallCoordinate p{patch.lower[0]+x*(patch.upper[0]-patch.lower[0]),d==3?.37*(patch.upper[1]-patch.lower[1])+patch.lower[1]:0.,H};
        const auto force=patch.Forcing(lifted,p,.2,data);if(count++)out<<',';Write(out,{p[0],p[1],patch.SampleK(lifted,p),force[0],force[1]});
    }
    out<<"],\"constraint_columns\":[";count=0;
    for(int block=0;block<patch.InnerPressure();++block)
    {
        const auto found=constraints.lower_bound(patch.offset[block]);
        Require(found!=constraints.end() && found->first<patch.offset[block]+patch.size[block],"missing constrained column control");
        const int column=found->first;const double h=2e-6;auto plus=state,minus=state;
        plus[column]+=h;minus[column]-=h;
        const auto rp=patch.Residual(plus,history,5.3,.2,data,false),rm=patch.Residual(minus,history,5.3,.2,data,false);
        if(count++)out<<',';
        Vector derivative(rp.size());for(int i=0;i<rp.size();++i)derivative[i]=(rp[i]-rm[i])/(2*h);
        out<<"{\"column\":"<<column<<",\"derivative\":";Write(out,derivative);out<<'}';
    }
    out<<"]}";
}
}
int main(int argc,char **argv)
{
    try
    {
        Require(argc==3,"usage: RANSEquilibriumEvolutionProbe SESSION.xml REPORT.json");
        auto world=LibUtilities::SessionReader::CreateInstance(2,argv);world->InitSession();
        char program[]="wall-global-probe";char *args[]{program,nullptr};auto filenames=world->GetFilenames();
        auto serial=LibUtilities::SessionReader::CreateInstance(1,args,filenames,LibUtilities::CommSerial::create(1,args));serial->InitSession();
        serial->SetSolverInfo("GlobalSysSoln","DirectFull");auto graph=SpatialDomains::MeshGraphIO::Read(serial);
        RANSWallPatchParameters parameters;serial->LoadParameter("RANSWallHeight",parameters.geometry.height,2.);
        serial->LoadParameter("RANSWallRoughness",parameters.geometry.roughness,.1);serial->LoadParameter("Kinvis",parameters.geometry.molecularNu,1.5e-5);
        serial->LoadParameter("RANSWallSampleHeight",parameters.geometry.sampleHeight,1.);serial->LoadParameter("RANSWallTailStart",parameters.geometry.tailStart,1.);
        serial->LoadParameter("RANSWallTangentialOrder",parameters.tangentialOrder,12);serial->LoadParameter("RANSWallNormalOrder",parameters.normalOrder,16);
        RANSWallEvolutionPatch patch(serial,graph,world->GetComm(),parameters);const int d=patch.Dimension();
        std::ostringstream bytes;bytes<<std::setprecision(17)<<"{\"schema\":1,\"dimension\":"<<d<<",\"ranks\":"<<world->GetComm()->GetSize()<<",\"offsets\":[";
        for(int b=0;b<patch.offset.size();++b){if(b)bytes<<',';bytes<<patch.offset[b];}bytes<<"],\"parameters\":";
        Write(bytes,{parameters.geometry.molecularNu,parameters.geometry.height,parameters.geometry.roughness,parameters.geometry.sampleHeight,
            parameters.geometry.tailStart,parameters.model.Q()*.16,parameters.model.Q(),parameters.model.Kappa()});
        bytes<<",\"fields\":[";int fieldCount=0;
        for(const auto &field:patch.fields)
        {
            if(fieldCount++)bytes<<',';
            Describe(bytes,field);
        }
        bytes<<"],\"native_evaluations\":[";fieldCount=0;
        for(const auto &field:patch.fields)
        {
            if(fieldCount++)bytes<<',';
            Vector values(field->Size());for(int i=0;i<values.size();++i)values[i]=.1*std::sin(.3*i)+.2*std::cos(.17*i);
            bytes<<"{\"coefficients\":";Write(bytes,values);bytes<<",\"samples\":[";int count=0;
            for(const auto &e:field->Elements())
            {
                RANSWallCoordinate p{e.lower[0]+.37*(e.upper[0]-e.lower[0]),d==3?e.lower[1]+.41*(e.upper[1]-e.lower[1]):0.,e.lower[2]+.43*(e.upper[2]-e.lower[2])};
                const auto row=field->Point(p);const auto gradient=row.Gradient(values);const auto native=field->NativeEvaluate(values,p);
                const Vector fast{row.Value(values),gradient[0],gradient[1],gradient[2]};
                for(int a=0;a<4;++a)Require(std::abs(native[a]-fast[a])<1e-10*std::max(1.,std::abs(native[a])),"fast/native tensor evaluation differs");
                if(count++)bytes<<',';
                Vector record(11);
                std::copy(p.begin(),p.end(),record.begin());
                std::copy(fast.begin(),fast.end(),record.begin()+3);
                std::copy(native.begin(),native.end(),record.begin()+7);Write(bytes,record);
            }
            bytes<<"]}";
        }
        RANSWallEvolutionData zero(d);const auto ref=patch.NeutralReference(.4,{1,0});zero.fixedBoundary=patch.FixedReferenceBoundary(ref);
        const auto constrained=patch.Constraints(ref,0,zero);
        std::vector<int> free;for(int b=0;b<d;++b)for(int i=0;i<patch.size[b];++i)if(!constrained.count(patch.offset[b]+i))free.push_back(patch.offset[b]+i);
        const int np=patch.size[patch.InnerPressure()]+patch.size[patch.OuterPressure()];WallMatrix B(np,Vector(free.size()));
        for(int side=0;side<2;++side)
        {
            const int block=patch.InnerPressure()+side,start=side?patch.size[patch.InnerPressure()]:0;
            const auto region=side?patch.outerVelocity:patch.innerVelocity;
            for(const auto &q:region->Quadrature(8))
            {
                const auto v=patch.fields[0]->Point(q.coordinate),p=patch.fields[block]->Point(q.coordinate);
                for(int j=0;j<free.size();++j)
                {
                    const int b=free[j]/patch.size[0],id=free[j]%patch.size[0],axis=patch.PhysicalAxis(b);
                    for(int i=0;i<v.ids.size();++i)if(v.ids[i]==id)for(int k=0;k<p.ids.size();++k)
                        B[start+p.ids[k]][j]+=q.weight*p.phi[k]*v.derivative[axis][i];
                }
            }
        }
        bytes<<"],\"pressure_divergence\":";Rows(bytes,B);bytes<<",\"free_velocity\":[";
        for(int i=0;i<free.size();++i){if(i)bytes<<',';bytes<<free[i];}bytes<<']';
        if(world->GetComm()->GetRank()==0)
        {std::ofstream layout(std::string(argv[2])+".layout.json");layout<<bytes.str()<<"}\n";}
        bytes<<",\"generic_cases\":[";
        int count=0;for(double chi:{0.,.37,1.}){if(count++)bytes<<',';Generic(bytes,serial,graph,world->GetComm(),parameters,chi);}bytes<<"]}\n";
        if(world->GetComm()->GetRank()==0){std::ofstream output(argv[2]);output<<bytes.str();Require(static_cast<bool>(output),"cannot write evolution probe report");}
        world->Finalise();return 0;
    }
    catch(const std::exception &error){std::cerr<<error.what()<<'\n';return 1;}
}
