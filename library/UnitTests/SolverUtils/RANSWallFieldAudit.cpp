///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <LibUtilities/BasicUtils/FieldIO.h>
#include <MultiRegions/ContField.h>
#include <SpatialDomains/MeshGraphIO.h>
#include <SolverUtils/Turbulence/RANSWallGeometry.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSession.h>
#include <LocalRegions/HexExp.h>
#include <StdRegions/StdQuadExp.h>
#include <StdRegions/StdSegExp.h>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <set>
#include <functional>

using namespace Nektar;
namespace
{
void ProjectReferenceFace(const std::shared_ptr<StdRegions::StdQuadExp> &face,
    const Array<OneD,NekDouble> &input,Array<OneD,NekDouble> &coefficients)
{
    const int nx = face->GetNumPoints(0),ny = face->GetNumPoints(1),nc = face->GetNcoeffs();
    Vmath::Zero(nc,coefficients,1);
    for (int edge = 0; edge < 4; ++edge)
    {
        auto segment = MemoryManager<StdRegions::StdSegExp>::AllocateSharedPtr(face->GetBasis(edge%2)->GetBasisKey());
        Array<OneD,NekDouble> values(segment->GetTotPoints()),trace(segment->GetNcoeffs());
        for (int q = 0; q < values.size(); ++q)
        {
            // TraceCoeffMap uses increasing tensor coordinates also on the
            // upper/left edges, as in LocalRegions::QuadExp. The native
            // StdQuadExp constrained transform instead reverses these edges.
            const int index = edge == 0 ? q : edge == 1 ? nx-1+nx*q :
                              edge == 2 ? nx*(ny-1)+q : nx*q;
            values[q] = input[index];
        }
        segment->FwdTransBndConstrained(values,trace);
        Array<OneD,unsigned int> map;Array<OneD,int> sign;
        face->GetTraceToElementMap(edge,map,sign);
        for (int q = 0; q < trace.size(); ++q) coefficients[map[q]] = sign[q]*trace[q];
    }
    Array<OneD,unsigned int> interior;face->GetInteriorMap(interior);
    if (interior.size())
    {
        Array<OneD,NekDouble> residual(face->GetTotPoints()),load(nc);
        face->BwdTrans(coefficients,residual);
        Vmath::Vsub(residual.size(),input,1,residual,1,residual,1);
        face->IProductWRTBase(residual,load);
        StdRegions::StdMatrixKey key(StdRegions::eMass,face->DetShapeType(),*face);
        const auto inverse = face->GetStdStaticCondMatrix(key)->GetBlock(1,1);
        for (int i = 0; i < interior.size(); ++i)
        {
            double value = 0;
            for (int j = 0; j < interior.size(); ++j) value += (*inverse)(i,j)*load[interior[j]];
            coefficients[interior[i]] = value;
        }
    }
}

void ProjectTerrainInitial(const MultiRegions::ContFieldSharedPtr &field,
                          const Array<OneD,NekDouble> &raw,int dimension)
{
    if (dimension == 2)
    {
        field->FwdTransBndConstrained(raw,field->UpdateCoeffs());
        return;
    }
    // HexExp has no FwdTransBndConstrained implementation. Constrain all
    // local face traces with native reference-quad projections, then solve
    // the reference interior mass block. Local tensor coordinates avoid any
    // global-face orientation ambiguity and preserve span-independent input.
    for (int e = 0; e < field->GetExpSize(); ++e)
    {
        auto hex = std::dynamic_pointer_cast<LocalRegions::HexExp>(field->GetExp(e));
        ASSERTL0(hex,"Terrain initial preparation requires quadrilaterals/hexahedra");
        const int nc = hex->GetNcoeffs(),np = hex->GetTotPoints();
        auto reference = MemoryManager<StdRegions::StdHexExp>::AllocateSharedPtr(
            hex->GetBasis(0)->GetBasisKey(),hex->GetBasis(1)->GetBasisKey(),hex->GetBasis(2)->GetBasisKey());
        Array<OneD,NekDouble> coefficients = field->UpdateCoeffs()+field->GetCoeff_Offset(e);
        Array<OneD,NekDouble> input = raw+field->GetPhys_Offset(e);
        Vmath::Zero(nc,coefficients,1);
        std::vector<int> seen(nc,0);
        for (int face = 0; face < hex->GetNtraces(); ++face)
        {
            auto first = hex->GetTraceBasisKey(face,0),second = hex->GetTraceBasisKey(face,1);
            auto surface = MemoryManager<StdRegions::StdQuadExp>::AllocateSharedPtr(first,second);
            Array<OneD,NekDouble> values(surface->GetTotPoints()),trace(surface->GetNcoeffs());
            Array<OneD,int> physicalMap;hex->GetTracePhysMap(face,physicalMap);
            ASSERTL0(physicalMap.size() == values.size(),"Initial face quadrature mismatch");
            for (int j = 0; j < values.size(); ++j) values[j] = input[physicalMap[j]];
            ProjectReferenceFace(surface,values,trace);
            Array<OneD,unsigned int> map;Array<OneD,int> sign;
            hex->GetTraceToElementMap(face,map,sign,StdRegions::eDir1FwdDir1_Dir2FwdDir2);
            for (int j = 0; j < trace.size(); ++j)
            {
                const double value = sign[j]*trace[j];
                if (seen[map[j]])
                    ASSERTL0(std::abs(coefficients[map[j]]-value) <= 1e-10*std::max(1.,std::abs(value)),
                             "Inconsistent initial face/edge traces: element "+std::to_string(e)+
                             " face "+std::to_string(face)+" previous "+std::to_string(seen[map[j]]-1)+
                             " mode "+std::to_string(map[j])+" old "+std::to_string(coefficients[map[j]])+
                             " new "+std::to_string(value));
                coefficients[map[j]] = value;seen[map[j]] = face+1;
            }
        }
        Array<OneD,unsigned int> interior;hex->GetInteriorMap(interior);
        if (interior.size())
        {
            Array<OneD,NekDouble> residual(np),load(nc);
            hex->BwdTrans(coefficients,residual);
            Vmath::Vsub(np,input,1,residual,1,residual,1);
            reference->IProductWRTBase(residual,load);
            StdRegions::StdMatrixKey key(StdRegions::eMass,reference->DetShapeType(),*reference);
            const auto inverse = reference->GetStdStaticCondMatrix(key)->GetBlock(1,1);
            for (int i = 0; i < interior.size(); ++i)
            {
                double value = 0;
                for (int j = 0; j < interior.size(); ++j) value += (*inverse)(i,j)*load[interior[j]];
                coefficients[interior[i]] = value;
            }
        }
    }
}
}
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
        if (std::string(argv[2]) == "--terrain-initial" || std::string(argv[2]) == "--experimental-initial")
        {
            const bool experimental=std::string(argv[2])=="--experimental-initial";
            ASSERTL0(experimental==session->DefinesSolverInfo("RANSExperimentalSignature"),
                     "Use --experimental-initial for experimental profiles; legacy log preparation is separate");
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
            double jacMin=std::numeric_limits<double>::max(),jacMax=0.;
            for (int e=0;e<field->GetExpSize();++e)
                for (double j : field->GetExp(e)->GetGeomFactors()->GetJac())
                {
                    ASSERTL0(std::isfinite(j) && j>0,"Initial geometry requires finite positive Jacobians");
                    jacMin=std::min(jacMin,j);jacMax=std::max(jacMax,j);
                }
            std::ofstream geometryReport(std::string(argv[3])+".geometry.json");
            geometryReport << std::setprecision(17) << "{\"jacobian_min\":" << jacMin
                           << ",\"jacobian_max\":" << jacMax << "}\n";
            ASSERTL0(geometryReport.good(),"Cannot write initial geometry audit");
            Array<OneD,NekDouble> x(nq,0.0),y(nq,0.0),z(nq,0.0),raw(nq),projected(nq);
            field->GetCoords(x,y,z);
            std::vector<SolverUtils::RANSWallFrame> frames;
            for (int q = 0; q < nq; ++q) frames.push_back(terrain.Closest(x[q],dim == 2 ? y[q] : z[q]));
            auto evaluate = [&](const std::string &name,double eta) {
                auto f=session->GetFunction("RANSExperimentalInlet",name);
                return dim==2 ? f->Evaluate(0.,eta,0.,0.) : f->Evaluate(0.,0.,eta,0.);
            };
            auto primitive = [&](double eta) { return matching*eta+star/kap*eta*(std::log(eta/H)-1); };
            auto speed = [&](double eta) { return experimental ? evaluate("u",eta) : matching+star/kap*std::log(eta/H); };
            auto difference = [&](double a,double b) {
                if (!experimental) return primitive(a)-primitive(b);
                // Only used in the N..2N initial streamfunction blending band.
                // Independent Gauss integration of the actual inlet expression
                // avoids inventing a logarithmic primitive for measured data.
                constexpr double gx[8]={.09501250983763744,.2816035507792589,.4580167776572274,.6178762444026438,
                    .755404408355003,.8656312023878318,.9445750230732326,.9894009349916499};
                constexpr double gw[8]={.1894506104550685,.1826034150449236,.1691565193950025,.1495959888165767,
                    .1246289712555339,.0951585116824928,.0622535239386479,.027152459411754095};
                if (a==b) return 0.;
                const double low=std::min(a,b),high=std::max(a,b);double sum=0.;
                ASSERTL0(low>0,"Nonpositive experimental initial-profile coordinate");
                std::vector<double> cuts{low,high};
                for (double c=low*1.4;c<high;c*=1.4) cuts.push_back(c);
                for (const std::string key : {"RANSExperimentBlendBegin","RANSExperimentBlendEnd"})
                {
                    double c;session->LoadParameter(key,c);
                    if (c>low && c<high) cuts.push_back(c);
                }
                std::sort(cuts.begin(),cuts.end());
                for (int i=1;i<cuts.size();++i)
                {
                    const double mid=.5*(cuts[i]+cuts[i-1]),half=.5*(cuts[i]-cuts[i-1]);
                    for (int j=0;j<8;++j) sum+=half*gw[j]*(speed(mid-half*gx[j])+speed(mid+half*gx[j]));
                }
                return a>b ? sum : -sum;
            };
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
                    const double a = chi*speed(eta)+(derivative==0 ? 0 : derivative*difference(eta,b[0]));
                    const double gx = (1-chi)*b[1]+a*f.normal[0];
                    const double gz = (1-chi)*b[2]+a*f.normal[2];
                    if (name=="u") raw[q]=gz;
                    else if (name==(dim==2 ? "v" : "w")) raw[q]=-gx;
                    else if (experimental && (name=="k" || name=="eps" || name=="eps_tilde" || name=="p"))
                        raw[q]=chi*evaluate(name,eta)+(1-chi)*evaluate(name,b[0]);
                    else raw[q]=name=="k" ? qstar*star*star :
                        (name=="eps" || name=="eps_tilde") ? star*star*star/kap*(chi/eta+(1-chi)/b[0]) : 0.;
                }
                ProjectTerrainInitial(field,raw,dim);
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
            if (experimental)
            {
                metadata["RANSObservationConfiguration"]=RANSObservationConfiguration(session);
                metadata["RANSInitializationOnly"]="1";
            }
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
