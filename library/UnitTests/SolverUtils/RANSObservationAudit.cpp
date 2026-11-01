///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSTerrainWall.h>
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <MultiRegions/ContField.h>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <sstream>

using namespace Nektar;
using namespace Nektar::SolverUtils;
namespace
{
using Fields=RANSWallBoundary::Fields;
using Coordinate=RANSPhysicalPointSampler::Coordinate;
const std::array<std::string,6> stressNames{{"rans_R_xx","rans_R_yy","rans_R_zz","rans_R_xy","rans_R_xz","rans_R_yz"}};
double Dot(const RANSVector &a,const RANSVector &b)
{ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
std::vector<std::string> Split(const std::string &line)
{
    std::vector<std::string> result;std::istringstream stream(line);std::string value;
    while (std::getline(stream,value,',')) result.push_back(value);
    return result;
}
RANSVector Physical(const Fields &values,int q,int dim)
{ return {values[0][q],dim==3 ? values[1][q] : 0.,values[dim-1][q]}; }
Coordinate Native(const RANSVector &p,int dim)
{ return dim==2 ? Coordinate{p[0],p[2],0.} : p; }
}

int main(int argc,char **argv)
{
    if (argc != 5)
    {
        std::cerr << "Usage: RANSObservationAudit session.xml field.fld points.csv output.csv\n"
                     "Use --manufactured instead of field.fld for explicit expression tests.\n";
        return 2;
    }
    try
    {
        char scotch[]="--use-scotch";
        char *arguments[]{argv[0],argv[1],scotch,nullptr};
        auto session=LibUtilities::SessionReader::CreateInstance(3,arguments);
        auto comm=session->GetComm();
        auto graph=SpatialDomains::MeshGraphIO::Read(session);
        const int dim=graph->GetMeshDimension(),nf=dim+2;
        ASSERTL0((dim==2 || dim==3) && graph->GetSpaceDimension()==dim,"Observation supports full 2D/3D meshes");
        Array<OneD,MultiRegions::ExpListSharedPtr> fields(nf+1);
        for (int f=0;f<=nf;++f)
            fields[f]=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,session->GetVariable(f),true,false);
        auto field=fields[0];const int nq=field->GetTotPoints();
        const bool manufactured=std::string(argv[2])=="--manufactured";
        std::map<std::string,Array<OneD,NekDouble>> imported;
        LibUtilities::FieldMetaDataMap metadata;
        Array<OneD,NekDouble> x(nq),y(nq),z(nq);field->GetCoords(x,y,z);
        if (manufactured)
        {
            metadata["Time"]="0";
            for (const auto &name : session->GetVariables())
            {
                Array<OneD,NekDouble> raw(nq),values(nq);
                session->GetFunction("InitialConditions",name)->Evaluate(x,y,z,0.,raw);
                field->FwdTransLocalElmt(raw,field->UpdateCoeffs());
                field->BwdTrans(field->GetCoeffs(),values);imported[name]=values;
            }
        }
        else
        {
            std::vector<LibUtilities::FieldDefinitionsSharedPtr> definitions;
            std::vector<std::vector<NekDouble>> data;
            Array<OneD,int> ids(field->GetExpSize());std::set<int> expected;
            for (int e=0;e<ids.size();++e) { ids[e]=field->GetExp(e)->GetGeom()->GetGlobalID();expected.insert(ids[e]); }
            LibUtilities::FieldIO::CreateForFile(session,argv[2])->Import(argv[2],definitions,data,metadata,ids);
            ASSERTL0(!definitions.empty() && metadata.count("Time"),"Missing native field/time");
            ASSERTL0(metadata.count("RANSObservationConfiguration") &&
                     metadata.at("RANSObservationConfiguration")==RANSObservationConfiguration(session),
                     "Observation configuration differs from checkpoint (or field predates observation binding)");
            for (auto name : definitions.front()->m_fields)
            {
                Vmath::Zero(field->GetNcoeffs(),field->UpdateCoeffs(),1);std::set<int> seen;
                for (int b=0;b<definitions.size();++b)
                {
                    ASSERTL0(definitions[b]->m_numHomogeneousDir==0,"Homogeneous extensions are not supported");
                    for (int id : definitions[b]->m_elementIDs)
                        if (expected.count(id)) ASSERTL0(seen.insert(id).second,"Duplicate native field element");
                    field->ExtractDataToCoeffs(definitions[b],data[b],name,field->UpdateCoeffs());
                }
                ASSERTL0(seen==expected,"Incomplete native field coverage");
                Array<OneD,NekDouble> values(nq);field->BwdTrans(field->GetCoeffs(),values);imported[name]=values;
            }
        }
        Fields state(nf),body(nf),pressureForce(dim),convection(dim);
        for (int f=0;f<nf;++f)
        { state[f]=imported.at(session->GetVariable(f));body[f]=Array<OneD,NekDouble>(nq,0.); }
        RANSWallBoundary::Gradients gradient(dim);
        Array<OneD,NekDouble> mechanical(nq);
        for (int q=0;q<nq;++q) mechanical[q]=imported.at("p")[q]-2*state[dim][q]/3;
        for (int i=0;i<dim;++i)
        {
            gradient[i]=Fields(dim);pressureForce[i]=Array<OneD,NekDouble>(nq);
            convection[i]=Array<OneD,NekDouble>(nq,0.);
            field->PhysDeriv(MultiRegions::DirCartesianMap[i],mechanical,pressureForce[i]);
            for (int j=0;j<dim;++j)
            {
                gradient[i][j]=Array<OneD,NekDouble>(nq);
                field->PhysDeriv(MultiRegions::DirCartesianMap[j],state[i],gradient[i][j]);
                for (int q=0;q<nq;++q) convection[i][q]+=state[j][q]*gradient[i][j][q];
            }
        }
        RANSParameters parameters;
        session->LoadParameter("RANSCmu",parameters.cMu,.09);
        session->LoadParameter("RANSCepsilon1",parameters.cEpsilon1,1.44);
        session->LoadParameter("RANSCepsilon2",parameters.cEpsilon2,1.92);
        session->LoadParameter("RANSsigmaK",parameters.sigmaK,1.);
        session->LoadParameter("RANSsigmaEpsilon",parameters.sigmaEpsilon,1.3);
        session->LoadParameter("RANSYapCoefficient",parameters.yapCoefficient,.83);
        session->LoadParameter("RANSYapLength",parameters.yapLength,2.5);
        const bool cls=session->GetSolverInfo("TurbulenceModel")=="CLSHighRe";
        const auto reference=cls ? RANSNeutralReference::CLS(parameters) : RANSNeutralReference::Standard(parameters);
        double origin,nu,H,P;session->LoadParameter("RANSWallOrigin",origin,0.);
        session->LoadParameter("Kinvis",nu);session->LoadParameter("RANSStripHeight",H);
        session->LoadParameter("RANSStripPointP",P);
        const auto terrain=ReadRANSWallGeometry(session,origin);
        const bool log=session->GetSolverInfo("RANSWallModel")=="None";
        std::string advection,dealias;
        session->LoadSolverInfo("AdvectionType",advection,"Convective");
        session->LoadSolverInfo("SpectralHPDealiasing",dealias,"False");
        const bool auxiliarySupported=!log && !session->DefinesElement("NEKTAR/FORCING") &&
            advection=="Convective" && dealias=="False";
        std::unique_ptr<RANSTerrainWall> wall;
        if (auxiliarySupported)
        {
            wall=std::make_unique<RANSTerrainWall>(session,fields,dim,terrain,reference,true);
            wall->Update(state,state,gradient,pressureForce,convection,body,1.);
        }
        Fields full(nf+6);
        for (int f=0;f<nf;++f) full[f]=state[f];
        const std::array<std::array<int,2>,6> tensor{{{0,0},{1,1},{2,2},{0,1},{0,2},{1,2}}};
        for (int i=0;i<6;++i) full[nf+i]=Array<OneD,NekDouble>(nq);
        for (int q=0;q<nq;++q)
        {
            RANSTensor g{};
            for (int i=0;i<dim;++i) for (int j=0;j<dim;++j)
                g[i==1 && dim==2 ? 2 : i][j==1 && dim==2 ? 2 : j]=gradient[i][j][q];
            const auto c=EvaluateRANSClosure(reference.Model(),state[dim][q],state[dim+1][q],g,
                terrain.Closest(x[q],dim==2 ? y[q] : z[q]).distance,parameters);
            for (int i=0;i<6;++i) full[nf+i][q]=c.reynoldsStress[tensor[i][0]][tensor[i][1]];
        }
        for (int i=0;i<6;++i)
        {
            if (imported.count(stressNames[i])) full[nf+i]=imported.at(stressNames[i]);
            else if (manufactured && session->DefinesFunction("RANSProbeTestTensor",stressNames[i]))
                session->GetFunction("RANSProbeTestTensor",stressNames[i])->Evaluate(x,y,z,0.,full[nf+i]);
            else
            {
                field->FwdTransLocalElmt(full[nf+i],field->UpdateCoeffs());
                field->BwdTrans(field->GetCoeffs(),full[nf+i]);
            }
        }
        std::ifstream input(argv[3]);ASSERTL0(input.good(),"Cannot read observation points");
        std::string line;std::getline(input,line);ASSERTL0(line=="id,x,span,z","Expected points header id,x,span,z");
        std::vector<std::string> names;std::vector<Coordinate> points,matching;
        std::vector<RANSWallFrame> frames;std::set<std::string> unique;
        while (std::getline(input,line))
        {
            if (line.empty()) continue;
            auto words=Split(line);ASSERTL0(words.size()==4 && unique.insert(words[0]).second,"Invalid/duplicate observation ID");
            RANSVector p{};
            for (int i=0;i<3;++i)
            {
                size_t consumed=0;p[i]=std::stod(words[i+1],&consumed);
                ASSERTL0(consumed==words[i+1].size() && std::isfinite(p[i]),"Invalid observation coordinate");
            }
            names.push_back(words[0]);points.push_back(Native(p,dim));
            ASSERTL0(dim==3 || p[1]==0.,"2D observation points require span=0");
            frames.push_back(terrain.Closest(p[0],p[2]));
            matching.push_back(Native(terrain.Point(frames.back().x,H,p[1]),dim));
        }
        ASSERTL0(!points.empty(),"No observation points");
        auto sampler=RANSPhysicalPointSampler::Volume(field,points,dim);
        auto matchingSampler=RANSPhysicalPointSampler::Volume(field,matching,dim);
        const auto actual=sampler.Sample(full),atH=matchingSampler.Sample(state);
        Fields auxiliary;Array<OneD,int> wallIDs;
        if (wall) auxiliary=wall->ObserveAuxiliary(matching,wallIDs);
        std::ofstream out;if (comm->GetRank()==0) out.open(argv[4]);
        if (comm->GetRank()==0) out << std::setprecision(17)
            << "id,time,x,span,z,foot_x,eta,origin,kind,element,owner,inverse_residual,u,v,w,k,epsilon_model,R_xx,R_yy,R_zz,R_xy,R_xz,R_yz,aux_KP,aux_Ft,aux_Fspan,aux_dsUP,aux_dspanVP,aux_UP,aux_VP,aux_UA_P,aux_W_H_target,aux_H_tangent_defect,aux_H_normal_defect,aux_chi,aux_cutoff\n";
        AWFParameters awfParameters;
        session->LoadParameter("RANSAWFCmu",awfParameters.cMu,.09);
        session->LoadParameter("RANSAWFCLength",awfParameters.cLength,2.55);
        session->LoadParameter("RANSAWFSmoothIntercept",awfParameters.smoothIntercept,10.7);
        for (int q=0;q<points.size();++q)
        {
            std::string kind="outside_resolved";std::vector<std::string> values(24);
            auto value=[&](int index,double v) { ASSERTL0(std::isfinite(v),"Nonfinite observation");std::ostringstream s;s<<std::setprecision(17)<<v;values[index]=s.str(); };
            const auto &frame=frames[q];
            if (frame.distance>=H-1e-10 && sampler.Found(q))
            {
                kind="resolved_native";const auto u=Physical(actual,q,dim);
                for (int i=0;i<3;++i) value(i,u[i]);
                value(3,actual[dim][q]);value(4,actual[dim+1][q]);
                for (int i=0;i<6;++i) value(5+i,actual[nf+i][q]);
            }
            else if (frame.distance>=0 && frame.distance<H)
            {
                kind=log ? "unsupported_log_layer" : "unsupported_auxiliary_forcing";
                if (wall && wallIDs[q]>=0 && matchingSampler.Found(q))
                {
                    kind="auxiliary_velocity";
                    if (terrain.terrain==RANSTerrain::CosineHill &&
                        std::abs(std::abs(frame.x-terrain.centre)-terrain.hillHalfLength)<1e-10)
                        kind="unsupported_C1_foot";
                    else
                    {
                        const auto h=Physical(atH,q,dim);
                        const AWFVector uh{Dot(h,frame.tangent),h[1]},force{auxiliary[1][q],auxiliary[2][q]};
                        AnalyticalWallFunction local({nu,H,auxiliary[5][q],P},auxiliary[0][q],awfParameters);
                        const auto u=ReconstructAWFVelocity(local,frame,frame.distance,uh,force,auxiliary[3][q],auxiliary[4][q]);
                        const auto limit=ReconstructAWFVelocity(local,frame,H,uh,force,auxiliary[3][q],auxiliary[4][q]);
                        for (int i=0;i<3;++i) value(i,u[i]);
                        for (int i=0;i<5;++i) value(11+i,auxiliary[i][q]);
                        value(16,auxiliary[6][q]);value(17,auxiliary[7][q]);
                        value(18,local.Velocity(P,uh,force)[0]);value(19,auxiliary[8][q]);
                        value(20,Dot(limit,frame.tangent)-Dot(h,frame.tangent));
                        value(21,Dot(limit,frame.normal)-Dot(h,frame.normal));
                        value(22,auxiliary[9][q]);value(23,auxiliary[10][q]);
                    }
                }
            }
            if (comm->GetRank()==0)
            {
                const auto physical=dim==2 ? RANSVector{points[q][0],0.,points[q][1]} : points[q];
                out << names[q] << ',' << metadata.at("Time") << ',' << physical[0] << ',' << physical[1] << ',' << physical[2]
                    << ',' << frame.x << ',' << frame.distance << ',' << (manufactured ? "manufactured" :
                        metadata.count("RANSInitializationOnly") ? "prepared_initial" : "checkpoint") << ',' << kind
                    << ',' << sampler.Geometry(q) << ',' << sampler.Owner(q) << ',' << sampler.Residual(q);
                for (const auto &v : values) out << ',' << v;
                out << '\n';
            }
        }
        if (comm->GetRank()==0)
        {
            ASSERTL0(out.good(),"Cannot write observation output");
            std::ofstream provenance(std::string(argv[4])+".native.json");
            provenance << "{\"configuration\":\"" << (manufactured ? "manufactured" :
                metadata.at("RANSObservationConfiguration")) << "\"}\n";
            ASSERTL0(provenance.good(),"Cannot write observation configuration receipt");
        }
        session->Finalise();
    }
    catch (const std::exception &e) { std::cerr << e.what() << '\n';return 1; }
}
