///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallInterface.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <fstream>
#include <cmath>
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
void FixtureSources(std::ostream &out,const LibUtilities::SessionReaderSharedPtr &session,double left,double right,double H,double Z)
{
    if(!session->DefinesFunction("ScalarSourceOuter","eps_tilde"))return;
    out<<",\"fixture_source_probes\":[";
    for(int j=0;j<3;++j)
    {
        if(j)out<<',';
        const double x=(left+right)/2,r=H+(j+1)*(Z-H)/4;
        Write(out,{x,r,session->GetFunction("ScalarSourceOuter","k")->Evaluate(x,r),
            session->GetFunction("ScalarSourceOuter","eps_tilde")->Evaluate(x,r)});
    }
    out<<']';
}
void GenericCase(std::ostream &out,const std::shared_ptr<RANSWallFieldLayout> &layout,RANSWallPatchParameters parameters,double chi)
{
    parameters.chi=chi;RANSWallPatch patch(layout,parameters);
    const double K=parameters.model.Q()*.16,H=parameters.geometry.height;
    const double Z=layout->Top();
    RANSWallPatchData data;std::array<RANSWallTimeFunction,6> f;
    f[0]=[](double x,double r,double){return (.1+.03*x)*r;};
    f[1]=[Z](double x,double r,double){return .02*r+.01*x*r*(Z-r);};
    f[2]=[K,rK=parameters.geometry.sampleHeight](double x,double r,double){return K+.02*x-.02*(r-rK);};
    f[3]=[K,H,chi,q=parameters.model.Q(),kappa=parameters.model.Kappa()](double x,double,double){
        const double sampled=K+.02*x;
        return (1-chi)*std::pow(sampled,1.5)/(2.55*H)+chi*std::pow(sampled/q,1.5)/(kappa*H);};
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
    out<<']';
    auto potential=f;
    potential[0]=potential[1]=[](double,double,double){return 0.;};
    potential[5]=[K,rK=parameters.geometry.sampleHeight](double x,double r,double){
        return .31*x-.47*r+2./3*(K+.02*x-.02*(r-rK));};
    potential[4]=[](double x,double r,double){return 2*x+.1*r;};
    auto potentialData=data;potentialData.outerBody[0]=[](double,double,double){return .31;};
    potentialData.outerBody[1]=[](double,double,double){return -.47;};
    const auto potentialState=patch.Project(potential,.2);
    out<<",\"potential_cancellation\":[";int potentialCount=0;
    for(double fraction:{.13,.47,.91})
    {
        if(potentialCount++)out<<',';
        const double x=layout->Left()+fraction*(layout->Right()-layout->Left());
        const double value=patch.Forcing(potentialState,x,.2,potentialData);
        Require(std::abs(value)<1e-10,"outer pressure/body-force potential did not cancel");Write(out,{x,value});
    }
    out<<']';
    if(chi==.37)
    {
        const auto jacobian=patch.Jacobian(state,history,5.3,.2,data);
        const auto inner=patch.CondensableVelocity(constraints);
        RANSWallCondensation reduction(jacobian,residual,inner);
        auto rhs=residual;for(double &x:rhs)x=-x;
        const auto full=RANSWallAlgebra::Solve(jacobian,rhs),condensed=reduction.Solve();
        Vector difference(full.size());for(int i=0;i<full.size();++i)difference[i]=full[i]-condensed[i];
        Require(Maximum(difference)<1e-9*std::max(1.,Maximum(full)),"full/Schur Newton step mismatch");
        Vector zero(state.size());RANSWallCondensation action(jacobian,zero,inner);
        Vector direction(action.Retained().size());
        for(int j=0;j<direction.size();++j)
            if(!constraints.count(action.Retained()[j]) && action.Retained()[j]!=layout->offset[layout->Gauge])
                direction[j]=std::sin(.7*(j+1))/std::sqrt(direction.size());
        double mean=0;
        for(int j=0;j<direction.size();++j)
            for(int b=0;b<2;++b)
            {
                const int index=action.Retained()[j]-layout->offset[layout->InnerPressure+b];
                if(index>=0 && index<patch.PressureMeanLoads()[b].size())mean+=patch.PressureMeanLoads()[b][index]*direction[j];
            }
        for(int b=0;b<2;++b)
        {
            const auto constant=layout->Field(static_cast<RANSWallFieldLayout::Block>(layout->InnerPressure+b)).Project([](double,double){return 1.;});
            for(int j=0;j<direction.size();++j)
            {
                const int index=action.Retained()[j]-layout->offset[layout->InnerPressure+b];
                if(index>=0 && index<constant.size())direction[j]-=mean*constant[index];
            }
        }
        const auto lift=action.Recover(direction),fullAction=Multiply(jacobian,lift),reducedAction=Multiply(action.Matrix(),direction);
        out<<",\"condensation\":{\"jacobian\":";WriteMatrix(out,jacobian);
        out<<",\"inner\":[";for(int j=0;j<inner.size();++j){if(j)out<<',';out<<inner[j];}out<<']';
        out<<",\"retained\":[";for(int j=0;j<action.Retained().size();++j){if(j)out<<',';out<<action.Retained()[j];}out<<']';
        out<<",\"schur\":";WriteMatrix(out,reduction.Matrix());out<<",\"reduced_residual\":";Write(out,reduction.Residual());
        out<<",\"full_step\":";Write(out,full);out<<",\"condensed_step\":";Write(out,condensed);
        out<<",\"retained_direction\":";Write(out,direction);out<<",\"lifted_direction\":";Write(out,lift);
        out<<",\"full_action\":";Write(out,fullAction);out<<",\"reduced_action\":";Write(out,reducedAction);
        out<<",\"finite_difference_actions\":[";int count=0;
        for(double h:{1e-3,3e-4,1e-4})
        {
            if(count++)out<<',';
            auto plus=state,minus=state;for(int j=0;j<state.size();++j){plus[j]+=h*lift[j];minus[j]-=h*lift[j];}
            const auto a=patch.Residual(plus,history,5.3,.2,data),b=patch.Residual(minus,history,5.3,.2,data);
            Vector value(a.size());for(int j=0;j<value.size();++j)value[j]=(a[j]-b[j])/(2*h);
            out<<"{\"step\":"<<h<<",\"action\":";Write(out,value);out<<'}';
        }
        out<<"]}";
    }
    out<<'}';
}
int main(int argc,char **argv)
{
    try
    {
        Require(argc==3 || argc==4,"Use SESSION.xml REPORT.json [parallel|layout|sources]");
        const std::string mode=argc==4?argv[3]:"all";
        Require(mode=="all" || mode=="parallel" || mode=="layout" || mode=="sources","unknown native probe mode");
        const bool parallelOnly=mode=="parallel",layoutOnly=mode=="layout" || mode=="sources",skipLayoutChecks=parallelOnly || layoutOnly;
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
            const auto coefficients=skipLayoutChecks?Vector(field->Size()):field->Project(f);
            auto native=skipLayoutChecks?Vector(field->Size()):field->NativeLoad(f);Vector direct(field->Size());
            if(!skipLayoutChecks)
            {
                for(const auto &q:field->Quadrature())for(int j=0;j<q.ids.size();++j)
                    direct[q.ids[j]]+=q.weight*q.phi[j]*f(q.x,q.r);
                for(int i=0;i<direct.size();++i)Require(std::abs(native[i]-direct[i])<1e-11,"native CG dual assembly");
            }
            const auto &e=field->Elements()[0];const auto q=field->Point(e.left+.37*(e.right-e.left),e.bottom+.41*(e.top-e.bottom));
            const auto g=q.Gradient(coefficients);
            if(!skipLayoutChecks)Require(std::abs(q.Value(coefficients)-f(q.x,q.r))<1e-10 &&
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
            out<<']';
            if(!skipLayoutChecks)
            {
                out<<",\"nonpolynomial_projection\":";
                Write(out,field->Project([](double x,double r){return std::sin(.4*x+.2*r)+.01*r*r*r;}));
            }
            out<<'}';
        }
        out<<"],\"stage_unknowns\":"<<layout.Size()<<",\"offsets\":[";
        for(int b=0;b<layout.Blocks;++b){if(b)out<<',';out<<layout.offset[b];}out<<']';
        RANSWallPatchParameters parameters;parameters.geometry.height=H;
        session->LoadParameter("Kinvis",parameters.geometry.molecularNu,1.5e-5);
        session->LoadParameter("RANSWallRoughness",parameters.geometry.roughness,.1);
        session->LoadParameter("RANSWallSampleHeight",parameters.geometry.sampleHeight,1.);
        session->LoadParameter("RANSWallTailStart",parameters.geometry.tailStart,1.);
        session->LoadParameter("RANSWallNormalOrder",parameters.normalOrder,16);
        session->LoadParameter("RANSWallTangentialOrder",parameters.tangentialOrder,8);
        const double k=parameters.model.Q()*.16;
        if(layoutOnly)
        {
            out<<",\"layout_only\":true,\"parameters\":";
            Write(out,{parameters.geometry.molecularNu,H,parameters.geometry.roughness,parameters.geometry.sampleHeight,
                parameters.geometry.tailStart,k,parameters.model.Q(),parameters.model.Kappa()});
            if(mode=="sources")FixtureSources(out,session,layout.Left(),layout.Right(),H,layout.Top());
            out<<"}\n";Require(static_cast<bool>(out),"cannot write native layout metadata");session->Finalise();return 0;
        }
        RANSWallPatch patch(layoutPtr,parameters,!parallelOnly);
        SolverUtils::EquilibriumAWFProfiles profile(parameters.geometry,k,parameters.model);
        const double ep=profile.Dissipation(H);
        out<<",\"parallel_only\":"<<(parallelOnly?"true":"false")<<",\"parallel_cases\":[";
        int parallelCount=0;
        for(double force:{-.02,0.,.03})
        {
            if(parallelCount++)out<<',';
            const auto reduction=patch.ParallelReduction(k,profile.Reference().matchingVelocity,force);
            out<<"{\"radial_unknowns\":"<<reduction.radialUnknowns<<",\"metrics\":";Write(out,reduction.metrics);
            out<<",\"coefficients\":";Write(out,reduction.coefficients);out<<",\"samples\":";WriteMatrix(out,reduction.samples);out<<'}';
        }
        out<<']';
        if(parallelOnly)
        {
            out<<",\"parameters\":";Write(out,{parameters.geometry.molecularNu,H,parameters.geometry.roughness,
                parameters.geometry.sampleHeight,parameters.geometry.tailStart,k,parameters.model.Q(),parameters.model.Kappa()});
            out<<"}\n";Require(static_cast<bool>(out),"cannot write parallel report");session->Finalise();return 0;
        }
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
        const auto topReference=parameters.model.AtHeight(layout.Top(),std::sqrt(k/parameters.model.Q()),H,profile.Reference().matchingVelocity);
        out<<",\"top_reference\":";
        Write(out,{layout.Top(),topReference.velocity,topReference.k,topReference.epsilon,topReference.eddyViscosity,
            topReference.velocityGradient,topReference.epsilonGradient,topReference.epsilonGradientFlux});
        out<<",\"top_reference_trace_samples\":[";int topCount=0;
        for(const auto &item:std::vector<std::pair<std::shared_ptr<RANSWallNativeField>,double>>{
            {layout.velocity,topReference.velocity},{layout.k,topReference.k},{layout.epsilon,topReference.epsilon}})
        {
            Vector trace(item.first->Size());
            const auto values=item.first->Boundary([&](double,double){return item.second;},
                [&](int edge,const RANSWallNativeField::Element &e){return edge==2 && std::abs(e.top-layout.Top())<1e-11;});
            for(auto [id,value]:values)trace[id]=value;
            if(topCount++)out<<',';
            Vector samples;
            for(double fraction:{.13,.47,.91})samples.push_back(item.first->Point(layout.Left()+fraction*(layout.Right()-layout.Left()),layout.Top(),-1).Value(trace));
            Write(out,samples);
        }
        out<<']';
        Vector ineffectiveGroundLoad(layout.velocity->Size());
        const auto noSlip=layout.velocity->Boundary([](double,double){return 0.;},
            [](int edge,const RANSWallNativeField::Element &e){return edge==0 && std::abs(e.bottom)<1e-11;});
        for(const auto &p:interface.Points())
        {
            const auto ground=layout.velocity->Point(p.x,0,1);
            for(int j=0;j<ground.ids.size();++j)if(!noSlip.count(ground.ids[j]))
                ineffectiveGroundLoad[ground.ids[j]]+=p.weight*ground.phi[j]*profile.Reference().traction;
        }
        out<<",\"strong_Dirichlet_ground_load_control\":";Write(out,ineffectiveGroundLoad);
        Require(Maximum(ineffectiveGroundLoad)<1e-14,"no-slip free tests have a nonzero ground trace");
        out<<",\"parameters\":";
        Write(out,{parameters.geometry.molecularNu,H,parameters.geometry.roughness,parameters.geometry.sampleHeight,
            parameters.geometry.tailStart,k,parameters.model.Q(),parameters.model.Kappa()});
        out<<",\"generic_cases\":[";int generic=0;
        for(double chi:{0.,.37,1.}){if(generic++)out<<',';GenericCase(out,layoutPtr,parameters,chi);}out<<']';
        FixtureSources(out,session,layout.Left(),layout.Right(),H,layout.Top());
        out<<"}\n";
        Require(static_cast<bool>(out),"cannot write layout report");
        std::cout<<"Native multiregion layout verified: "<<layout.Size()<<" stage unknowns.\n";
        session->Finalise();
    }
    catch(const std::exception &e){std::cerr<<e.what()<<'\n';return 1;}
    return 0;
}
