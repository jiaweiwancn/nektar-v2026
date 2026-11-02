///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallInterface.h>
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
void WriteMatrix(std::ostream &out,const RANSWallAlgebra::DenseMatrix &m)
{
    out<<'[';
    for(int i=0;i<m.size();++i){if(i)out<<',';Write(out,m[i]);}
    out<<']';
}
void GenericCase(std::ostream &out,const std::shared_ptr<RANSWallFieldLayout> &layout,RANSWallPatchParameters parameters,double chi)
{
    parameters.chi=chi;RANSWallPatch patch(layout,parameters);
    const double K=parameters.model.Q()*.16,H=parameters.geometry.height;
    SolverUtils::AWFProfileOptions options;options.requestedWeight=chi;
    SolverUtils::EquilibriumAWFProfiles profile(parameters.geometry,K,parameters.model,{},options);
    const double ep=profile.Dissipation(H),Z=layout->Top();
    RANSWallPatchData data;std::array<RANSWallTimeFunction,6> f;
    f[0]=[](double x,double r,double){return (.1+.03*x)*r;};
    f[1]=[Z](double x,double r,double){return .02*r+.01*x*r*(Z-r);};
    f[2]=[K,rK=parameters.geometry.sampleHeight](double,double r,double){return K+.01*(r-rK);};
    f[3]=[ep](double,double,double){return ep;};
    f[4]=[](double x,double r,double){return -.01*x+.02*r;};
    f[5]=[](double x,double r,double){return .04*x+.05*r;};
    for(int b=0;b<4;++b)data.boundary[b]=f[b];
    data.innerBody[0]=[](double x,double,double){return .01+.02*x;};
    data.innerBody[1]=[](double,double r,double){return -.02+.03*r;};
    data.outerBody[0]=[](double x,double r,double){return .04*x+.01*r;};
    data.outerBody[1]=[](double,double r,double){return -.03+.02*r;};
    data.innerKSource=[](double x,double r,double){return .03+.01*x+.02*r;};
    data.outerKSource=[](double x,double r,double){return -.02+.01*x+.03*r;};
    data.outerEpsilonSource=[](double x,double r,double){return .01*(1+x+r);};
    auto state=patch.Project(f,.2),history=state;
    for(int j=0;j<layout->offset[layout->InnerPressure];++j)history[j]*=.97;
    const auto residual=patch.Residual(state,history,5.3,.2,data);
    const auto constraints=patch.Constraints(state,.2,data);
    auto lifted=state;for(auto [row,value]:constraints)lifted[row]=value;
    out<<"{\"chi\":"<<chi<<",\"alpha\":5.3,\"time\":0.2,\"state\":";Write(out,state);
    out<<",\"lifted_state\":";Write(out,lifted);out<<",\"history\":";Write(out,history);
    out<<",\"residual\":";Write(out,residual);out<<",\"constraints\":[";
    int count=0;for(auto [row,value]:constraints){if(count++)out<<',';Write(out,{static_cast<double>(row),value});}
    out<<"],\"forcing_samples\":[";count=0;
    for(double fraction:{.13,.47,.91})
    {
        if(count++)out<<',';
        const double x=layout->Left()+fraction*(layout->Right()-layout->Left());
        Write(out,{x,patch.SampleK(lifted,x),patch.Forcing(lifted,x,.2,data)});
    }
    out<<"]}";
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
        auto layoutPtr=std::make_shared<RANSWallFieldLayout>(session,graph,H);
        auto &layout=*layoutPtr;
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
            const auto &e=field->Elements()[0];const auto q=field->Point(e.left+.37*(e.right-e.left),e.bottom+.41*(e.top-e.bottom));
            const auto g=q.Gradient(coefficients);
            Require(std::abs(q.Value(coefficients)-f(q.x,q.r))<1e-10 &&
                std::abs(g[0]-.2)<1e-10 && std::abs(g[1]-.06*q.r)<1e-10,"native physical polynomial/derivative");
            field->Store(coefficients);
            out<<"{\"elements\":"<<field->Elements().size()<<",\"global_coefficients\":"<<field->Size()
               <<",\"local_coefficients\":"<<field->Field()->GetNcoeffs()<<",\"probe\":";
            Write(out,{q.x,q.r,q.Value(coefficients),g[0],g[1]});out<<",\"native_load\":";
            Write(out,native);out<<",\"mesh\":[";
            const auto map=field->Field()->GetLocalToGlobalMap();
            for(int e=0;e<field->Elements().size();++e)
            {
                if(e)out<<',';
                const auto &el=field->Elements()[e];
                out<<"{\"degree\":"<<el.degree<<",\"geometry\":"<<el.geometry<<",\"bounds\":";
                Write(out,{el.left,el.right,el.bottom,el.top});out<<",\"ids\":[";
                for(int j=0;j<field->Field()->GetExp(e)->GetNcoeffs();++j)
                {if(j)out<<',';out<<map->GetLocalToGlobalMap(field->Field()->GetCoeff_Offset(e)+j);}
                out<<"],\"sign\":[";
                for(int j=0;j<field->Field()->GetExp(e)->GetNcoeffs();++j)
                {if(j)out<<',';out<<map->GetLocalToGlobalSign(field->Field()->GetCoeff_Offset(e)+j);}
                out<<"]}";
            }
            out<<"]}";
        }
        out<<"],\"stage_unknowns\":"<<layout.Size()<<",\"offsets\":[";
        for(int b=0;b<layout.Blocks;++b){if(b)out<<',';out<<layout.offset[b];}out<<']';
        RANSWallPatchParameters parameters;parameters.geometry.height=H;
        session->LoadParameter("Kinvis",parameters.geometry.molecularNu,1.5e-5);
        session->LoadParameter("RANSWallRoughness",parameters.geometry.roughness,.1);
        session->LoadParameter("RANSWallSampleHeight",parameters.geometry.sampleHeight,1.);
        session->LoadParameter("RANSWallTailStart",parameters.geometry.tailStart,1.);
        RANSWallPatch patch(layoutPtr,parameters);
        const double k=parameters.model.Q()*.16;
        SolverUtils::EquilibriumAWFProfiles profile(parameters.geometry,k,parameters.model);
        const double ep=profile.Dissipation(H);
        RANSWallPatchData data;
        data.boundary[0]=[](double,double r,double){return .1*r;};
        data.boundary[2]=[&](double,double,double){return k;};
        data.boundary[3]=[&](double,double,double){return ep;};
        std::array<RANSWallTimeFunction,6> functions;
        functions[0]=data.boundary[0];functions[1]=data.boundary[1];functions[2]=data.boundary[2];
        functions[3]=data.boundary[3];functions[4]=functions[5]=[](double,double,double){return 0.;};
        const auto state=patch.Project(functions,0);
        const auto residual=patch.Residual(state,state,5,0,data);
        out<<",\"state\":";Write(out,state);out<<",\"residual\":";Write(out,residual);
        out<<",\"residual_norms\":";Write(out,patch.ResidualNorms(residual,patch.Constraints(state,0,data)));
        RANSWallInterface interface(layout);
        out<<",\"trace_constraint\":";WriteMatrix(out,interface.TraceConstraint());
        const auto loads=interface.Loads([](double x,double){return .17+.03*x;});
        out<<",\"shared_loads\":[";Write(out,loads[0]);out<<',';Write(out,loads[1]);out<<']';
        const auto merged=interface.Merge(loads);Require(Maximum(merged)<1e-13,"opposite common-flux dual loads");
        out<<",\"merged_shared_load\":";Write(out,merged);
        const auto pressure=interface.ReferencePressure(parameters.model,k);
        Require(pressure.velocityMaximum<1e-11 && std::abs(pressure.globalMean)<1e-12 &&
            pressure.tractionMismatch<1e-11,"native pressure/gauge/full-traction coupling");
        out<<",\"pressure\":{\"divergence_matrix\":";WriteMatrix(out,pressure.divergence);
        out<<",\"system\":";WriteMatrix(out,pressure.system);out<<",\"rhs\":";Write(out,pressure.rhs);
        out<<",\"solution\":";Write(out,pressure.solution);out<<",\"metrics\":";
        Write(out,{pressure.velocityMaximum,pressure.globalMean,pressure.pressureJump,pressure.tractionMismatch});out<<'}';
        out<<",\"parameters\":";
        Write(out,{parameters.geometry.molecularNu,H,parameters.geometry.roughness,parameters.geometry.sampleHeight,
            parameters.geometry.tailStart,k,parameters.model.Q(),parameters.model.Kappa()});
        out<<",\"generic_cases\":[";int generic=0;
        for(double chi:{0.,.37,1.}){if(generic++)out<<',';GenericCase(out,layoutPtr,parameters,chi);}out<<']';
        out<<"}\n";
        Require(static_cast<bool>(out),"cannot write layout report");
        std::cout<<"Native multiregion layout verified: "<<layout.Size()<<" stage unknowns.\n";
        session->Finalise();
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
    return 0;
}
