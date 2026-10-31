///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <MultiRegions/ContField.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <SolverUtils/Turbulence/RANSWallGeometry.h>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>

using namespace Nektar;
int main(int argc, char **argv)
{
    if (argc != 4) { std::cerr << "Usage: RANSWallFieldAudit session.xml field.fld output.csv (or --initial output.csv / --terrain-initial output.fld)\n"; return 2; }
    try
    {
        char *arguments[]{argv[0],argv[1],nullptr};
        auto session = LibUtilities::SessionReader::CreateInstance(2,arguments);
        ASSERTL0(session->GetComm()->GetSize() == 1,"Import all partitions with one audit process");
        auto graph = SpatialDomains::MeshGraphIO::Read(session);
        auto field = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,graph,session->GetVariable(0),true,false);
        if (std::string(argv[2]) == "--terrain-initial")
        {
            SolverUtils::RANSWallGeometry terrain;
            std::string kind; session->LoadSolverInfo("RANSWallGeometry",kind,"Plane");
            terrain.terrain = kind == "CosineHill" ? SolverUtils::RANSTerrain::CosineHill : SolverUtils::RANSTerrain::Plane;
            session->LoadParameter("RANSTerrainHeight",terrain.hillHeight,.04);
            session->LoadParameter("RANSTerrainHalfLength",terrain.hillHalfLength,.1);
            session->LoadParameter("RANSTerrainCentre",terrain.centre,0);
            session->LoadParameter("RANSTerrainSlope",terrain.slope,0);
            session->LoadParameter("RANSWallOrigin",terrain.origin,0);
            terrain.Validate();
            double H,N,star,kap,matching,qstar,top;
            session->LoadParameter("RANSStripHeight",H);
            session->LoadParameter("RANSStripPointN",N);
            session->LoadParameter("Ustar",star); session->LoadParameter("Kappa",kap);
            session->LoadParameter("MatchingVelocity",matching);session->LoadParameter("Qstar",qstar);session->LoadParameter("Z",top);
            const int dim = graph->GetMeshDimension(),nq = field->GetTotPoints();
            Array<OneD,NekDouble> x(nq,0.0),y(nq,0.0),z(nq,0.0),raw(nq),projected(nq);
            field->GetCoords(x,y,z);
            std::vector<SolverUtils::RANSWallFrame> frames;
            for (int q = 0; q < nq; ++q) frames.push_back(terrain.Closest(x[q],dim == 2 ? y[q] : z[q]));
            auto primitive = [&](double eta) { return matching*eta+star/kap*eta*(std::log(eta/H)-1); };
            auto speed = [&](double eta) { return matching+star/kap*std::log(eta/H); };
            // Blend G(normal distance) through N into the smooth roof-compatible
            // streamfunction G(e), e=z-h(x)*((top-z)/(top-h(x)))^2, by 2N.
            // Near H/N the tangential log profile is independent of wall foot,
            // so its initial strip flux has zero tangential divergence. The
            // subsequent AWF condition still solves integrated continuity.
            auto base = [&](double px,double pz) {
                const auto surface = terrain.At(px);
                const double height = surface.z-terrain.origin;
                const double slope = surface.tangent[2]/surface.tangent[0],gap = top-height;
                const double e = pz-height*std::pow((top-pz)/gap,2);
                const double ex = -slope*std::pow(top-pz,2)*(top+height)/std::pow(gap,3);
                const double ez = 1+2*height*(top-pz)/(gap*gap);
                return std::array<double,3>{{e,speed(e)*ex,speed(e)*ez}};
            };
            auto definitions = field->GetFieldDefinitions();
            std::vector<std::vector<NekDouble>> data(definitions.size());
            std::ofstream out(std::string(argv[3])+".csv");
            out << std::setprecision(17) << "field,x,y,z,input,projected,foot_x,distance,normal_x,normal_z\n";
            for (const auto &name : session->GetVariables())
            {
                for (int q = 0; q < nq; ++q)
                {
                    const auto &f = frames[q];
                    const double eta = f.distance,vertical = (dim == 2 ? y[q] : z[q])-terrain.origin;
                    ASSERTL0(eta > 0 && vertical > 0,"Terrain start requires positive heights");
                    const double r = std::clamp((eta-N)/N,0.0,1.0);
                    const double chi = 1-10*r*r*r+15*std::pow(r,4)-6*std::pow(r,5);
                    const double derivative = -30*r*r*(1-r)*(1-r)/N;
                    const auto b = base(x[q],vertical);
                    const double a = chi*speed(eta)+derivative*(primitive(eta)-primitive(b[0]));
                    const double gx = (1-chi)*b[1]+a*f.normal[0];
                    const double gz = (1-chi)*b[2]+a*f.normal[2];
                    raw[q] = name == "u" ? gz :
                        name == (dim == 2 ? "v" : "w") ? -gx :
                        name == "k" ? qstar*star*star :
                        (name == "eps" || name == "eps_tilde") ? star*star*star/kap*(chi/eta+(1-chi)/b[0]) : 0;
                }
                field->FwdTransBndConstrained(raw,field->UpdateCoeffs());
                field->AvgAssemble();field->GlobalToLocal();field->BwdTrans(field->GetCoeffs(),projected);
                for (int q = 0; q < nq; ++q)
                    out << name << ',' << x[q] << ',' << y[q] << ',' << z[q] << ',' << raw[q] << ',' << projected[q]
                        << ',' << frames[q].x << ',' << frames[q].distance << ',' << frames[q].normal[0] << ',' << frames[q].normal[2] << '\n';
                for (int block = 0; block < definitions.size(); ++block)
                {
                    definitions[block]->m_fields.push_back(name);
                    field->AppendFieldData(definitions[block],data[block],field->UpdateCoeffs());
                }
            }
            LibUtilities::FieldMetaDataMap metadata{{"Time","0"},{"ChkFileNum","0"}};
            LibUtilities::FieldIO::CreateDefault(session)->Write(argv[3],definitions,data,metadata,false);
            session->Finalise();return 0;
        }
        if (std::string(argv[2]) == "--initial")
        {
            Array<OneD,NekDouble> x(field->GetTotPoints(),0.0),y(x.size(),0.0),z(x.size(),0.0),raw(x.size()),reconstructed(x.size());
            field->GetCoords(x,y,z);
            std::ofstream out(argv[3]); out << std::setprecision(17) << "field,x,y,z,input,projected\n";
            for (const auto &name : session->GetVariables())
            {
                session->GetFunction("InitialConditions",name)->Evaluate(x,y,z,0,raw);
                field->FwdTransLocalElmt(raw,field->UpdateCoeffs());
                field->AvgAssemble(); field->GlobalToLocal();
                field->BwdTrans(field->GetCoeffs(),reconstructed);
                double a = raw[0], b = raw[0], c = reconstructed[0], d = reconstructed[0];
                for (int q = 0; q < raw.size(); ++q)
                {
                    a = std::min(a,raw[q]); b = std::max(b,raw[q]);
                    c = std::min(c,reconstructed[q]); d = std::max(d,reconstructed[q]);
                    out << name << ',' << x[q] << ',' << y[q] << ',' << z[q] << ',' << raw[q] << ',' << reconstructed[q] << '\n';
                }
                std::cout << name << " input [" << a << ',' << b << "] projected [" << c << ',' << d << "]\n";
            }
            double jmin = std::numeric_limits<double>::max(), jmax = -jmin;
            for (int e = 0; e < field->GetExpSize(); ++e)
                for (double j : field->GetExp(e)->GetGeomFactors()->GetJac())
                { jmin = std::min(jmin,j); jmax = std::max(jmax,j); }
            std::cout << "Jacobian [" << jmin << ',' << jmax << "]\n";
            session->Finalise(); return 0;
        }
        std::vector<LibUtilities::FieldDefinitionsSharedPtr> definitions;
        std::vector<std::vector<NekDouble>> data;
        LibUtilities::FieldMetaDataMap metadata;
        LibUtilities::FieldIO::CreateForFile(session,argv[2])->Import(argv[2],definitions,data,metadata);
        ASSERTL0(!definitions.empty() && metadata.count("Time"),"Empty field or missing time");
        auto names = definitions.front()->m_fields;
        Array<OneD,Array<OneD,NekDouble>> values(names.size());
        std::set<int> expected;
        for (int e = 0; e < field->GetExpSize(); ++e) expected.insert(field->GetExp(e)->GetGeom()->GetGlobalID());
        for (int f = 0; f < names.size(); ++f)
        {
            std::set<int> seen;
            Vmath::Zero(field->GetNcoeffs(),field->UpdateCoeffs(),1);
            for (int b = 0; b < definitions.size(); ++b)
            {
                auto &d = definitions[b];
                ASSERTL0(d->m_fields == names && d->m_numHomogeneousDir == 0,"Inconsistent field definitions");
                for (int id : d->m_elementIDs) ASSERTL0(expected.count(id) && seen.insert(id).second,"Duplicate field element");
                field->ExtractDataToCoeffs(d,data[b],names[f],field->UpdateCoeffs());
            }
            ASSERTL0(seen == expected,"Incomplete MPI field coverage");
            values[f] = Array<OneD,NekDouble>(field->GetTotPoints());
            field->BwdTrans(field->GetCoeffs(),values[f]);
        }
        Array<OneD,NekDouble> x(field->GetTotPoints(),0.0),y(x.size(),0.0),z(x.size(),0.0);
        field->GetCoords(x,y,z);
        std::ofstream out(argv[3]); ASSERTL0(out.good(),"Cannot open wall audit output");
        out << "element,point,x,y,z,time";
        for (const auto &name : names) out << ',' << name;
        out << '\n' << std::setprecision(17);
        for (int e = 0; e < field->GetExpSize(); ++e)
            for (int q = 0; q < field->GetExp(e)->GetTotPoints(); ++q)
            {
                const int offset = field->GetPhys_Offset(e)+q;
                out << field->GetExp(e)->GetGeom()->GetGlobalID() << ',' << q << ',' << x[offset] << ',' << y[offset]
                    << ',' << z[offset] << ',' << metadata.at("Time");
                for (const auto &v : values) out << ',' << v[offset];
                out << '\n';
            }
        ASSERTL0(out.good(),"Cannot finish wall audit output");
        std::cout << "Exported " << field->GetTotPoints() << " physical points and " << names.size() << " fields\n";
        session->Finalise();
    }
    catch (const std::exception &e) { std::cerr << e.what() << '\n'; return 1; }
}
