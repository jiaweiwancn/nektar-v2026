///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <fstream>
#include <iomanip>
#include <iostream>

using namespace Nektar;
using namespace Nektar::RANSWallAlgebra;
void Write(std::ostream &out,const Vector &v)
{
    out << '[';
    for(int i=0;i<v.size();++i) { if(i)out<<',';out<<v[i]; }
    out << ']';
}
int main(int argc,char **argv)
{
    try
    {
        Require(argc==3,"Use SESSION.xml REPORT.json");
        const std::string output=argv[2];
        auto session=LibUtilities::SessionReader::CreateInstance(2,argv);
        auto graph=SpatialDomains::MeshGraphIO::Read(session);
        double H;session->LoadParameter("RANSWallHeight",H,2.);
        RANSWallFieldLayout layout(session,graph,H);
        std::ofstream out(output);out<<std::setprecision(17)<<"{\"schema\":1,\"layout\":[";
        int count=0;
        for(const auto &field:{layout.velocity,layout.innerVelocity,layout.outerVelocity,
            layout.k,layout.epsilon,layout.innerPressure,layout.outerPressure})
        {
            if(count++)out<<',';
            auto f=[](double x,double r){return 1+.2*x+.03*r*r;};
            const auto coefficients=field->Project(f);
            auto native=field->NativeLoad(f);Vector direct(field->Size());
            for(const auto &q:field->Quadrature())for(int j=0;j<q.ids.size();++j)
                direct[q.ids[j]]+=q.weight*q.phi[j]*f(q.x,q.r);
            for(int i=0;i<direct.size();++i)Require(std::abs(native[i]-direct[i])<1e-11,"native CG dual assembly");
            const auto &e=field->Elements()[0];const auto q=field->Point(.37*(e.left+e.right),e.bottom+.41*(e.top-e.bottom));
            const auto g=q.Gradient(coefficients);
            Require(std::abs(q.Value(coefficients)-f(q.x,q.r))<1e-10 &&
                std::abs(g[0]-.2)<1e-10 && std::abs(g[1]-.06*q.r)<1e-10,"native physical polynomial/derivative");
            field->Store(coefficients);
            out<<"{\"elements\":"<<field->Elements().size()<<",\"global_coefficients\":"<<field->Size()
               <<",\"local_coefficients\":"<<field->Field()->GetNcoeffs()<<",\"probe\":";
            Write(out,{q.x,q.r,q.Value(coefficients),g[0],g[1]});out<<",\"native_load\":";
            Write(out,native);out<<'}';
        }
        out<<"],\"stage_unknowns\":"<<layout.Size()<<"}\n";
        Require(static_cast<bool>(out),"cannot write layout report");
        std::cout<<"Native multiregion layout verified: "<<layout.Size()<<" stage unknowns.\n";
        session->Finalise();
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
    return 0;
}
