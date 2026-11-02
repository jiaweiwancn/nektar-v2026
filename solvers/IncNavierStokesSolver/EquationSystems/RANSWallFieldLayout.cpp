///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallFieldLayout.h>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <StdRegions/StdSegExp.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>

namespace Nektar
{
namespace RANSWallAlgebra
{
void Require(bool valid, const std::string &message)
{ if (!valid) throw std::runtime_error("EP-AWF patch: "+message); }
Vector Solve(const DenseMatrix &matrix, const Vector &rhs)
{
    const int n = rhs.size();
    Require(matrix.size() == rhs.size() && n > 0, "invalid dense system");
    Vector column(n*n), solution(rhs); std::vector<int> pivots(n); int info;
    for (int i = 0; i < n; ++i)
    {
        Require(matrix[i].size() == rhs.size(), "nonsquare dense system");
        for (int j = 0; j < n; ++j) column[i+n*j] = matrix[i][j];
    }
    Lapack::Dgetrf(n, n, column.data(), n, pivots.data(), info);
    Require(info == 0, "singular uncondensed/condensed matrix");
    Lapack::Dgetrs('N', n, 1, column.data(), n, pivots.data(), solution.data(), n, info);
    Require(info == 0 && std::all_of(solution.begin(), solution.end(),
        [](double x) { return std::isfinite(x); }), "invalid dense solution");
    return solution;
}
Vector Multiply(const DenseMatrix &matrix, const Vector &value)
{
    Vector result(matrix.size());
    for (std::size_t i = 0; i < matrix.size(); ++i)
    {
        Require(matrix[i].size() == value.size(), "invalid matrix action");
        long double sum = 0;
        for (std::size_t j = 0; j < value.size(); ++j) sum += (long double)matrix[i][j]*value[j];
        result[i] = sum;
    }
    return result;
}
double Maximum(const Vector &value)
{
    double result = 0;
    for (double x : value) { Require(std::isfinite(x), "nonfinite vector"); result = std::max(result, std::abs(x)); }
    return result;
}
}
using namespace RANSWallAlgebra;
namespace
{
std::array<double,6> PhysicalBounds(SpatialDomains::Geometry *geom)
{
    std::array<double,6> box{};
    geom->GetVertex(0)->GetCoords(box[0],box[1],box[2]);
    for(int d=0;d<3;++d)box[d+3]=box[d];
    for(int v=1;v<geom->GetNumVerts();++v)
    {
        double x,y,z;geom->GetVertex(v)->GetCoords(x,y,z);
        const std::array<double,3> point{x,y,z};
        for(int d=0;d<3;++d){box[d]=std::min(box[d],point[d]);box[d+3]=std::max(box[d+3],point[d]);}
    }
    return box;
}
}

double RANSWallPoint::Value(const Vector &field, int offset) const
{
    // Evaluate about a native constant. This preserves the partition of unity
    // without cancellation of O(1) values on thin elements.
    const double constant=field.at(offset+ids[vertexPositions[0]])*vertexSigns[0];
    long double result = constant;
    for (std::size_t j = 0; j < ids.size(); ++j)
    {
        double coefficient=field.at(offset+ids[j]);
        for(int v=0;v<4;++v)if(j==vertexPositions[v])coefficient-=constant*vertexSigns[v];
        result += (long double)phi[j]*coefficient;
    }
    return result;
}
std::array<double, 2> RANSWallPoint::Gradient(const Vector &field, int offset) const
{
    long double x = 0, r = 0;
    const double constant=field.at(offset+ids[vertexPositions[0]])*vertexSigns[0];
    for (std::size_t j = 0; j < ids.size(); ++j)
    {
        double coefficient=field.at(offset+ids[j]);
        for(int v=0;v<4;++v)if(j==vertexPositions[v])coefficient-=constant*vertexSigns[v];
        x += (long double)dx[j]*coefficient; r += (long double)dr[j]*coefficient;
    }
    return {static_cast<double>(x), static_cast<double>(r)};
}

RANSWallNativeField::RANSWallNativeField(const LibUtilities::SessionReaderSharedPtr &session,
    const MultiRegions::ExpListSharedPtr &full, const std::vector<unsigned int> &elements)
{
    Require(session->GetComm()->GetSize() == 1 && !elements.empty(), "M3 native backend requires a nonempty serial subdomain");
    auto subset = MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(*full, elements, true);
    m_field = MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session, *subset);
    for (int e = 0; e < m_field->GetExpSize(); ++e)
    {
        const auto exp = m_field->GetExp(e); const auto geom = exp->GetGeom();
        Require(geom->GetShapeType() == LibUtilities::eQuadrilateral && geom->GetCoordim() == 2,
                "M3 supports physical 2D quadrilaterals only");
        Require(exp->GetGeomFactors()->GetGtype()==SpatialDomains::eRegular,
                "M3 requires affine geometry; curved patches belong to M5");
        const auto box = PhysicalBounds(geom);
        Element el{box[0], box[3], box[1], box[4], geom->GetGlobalID(),
                   static_cast<int>(exp->GetBasisNumModes(0))-1, {}, {}};
        Require(el.degree >= 2 && exp->GetBasisNumModes(1) == el.degree+1 &&
            exp->GetBasisType(0) == LibUtilities::eModified_A &&
            exp->GetBasisType(1) == LibUtilities::eModified_A, "M3 requires equal-order Modified_A tensor bases");
        const std::array<std::array<double, 2>, 4> vertices{{{el.left,el.bottom},
            {el.right,el.bottom},{el.right,el.top},{el.left,el.top}}};
        for (int v = 0; v < 4; ++v)
        {
            double x,y,z; geom->GetVertex(v)->GetCoords(x,y,z);
            Require(std::abs(x-vertices[v][0])<1e-11 && std::abs(y-vertices[v][1])<1e-11,
                    "M3 requires canonical affine axis-aligned quads");
        }
        Array<OneD,NekDouble> coef(exp->GetNcoeffs(),0.), values(exp->GetTotPoints()),
            dx(exp->GetTotPoints()), dr(exp->GetTotPoints());
        for (int i = 0; i < exp->GetNcoeffs(); ++i)
        {
            coef[i] = 1; exp->BwdTrans(coef,values); exp->PhysDeriv(0,values,dx); exp->PhysDeriv(1,values,dr);
            el.derivativeX.emplace_back(dx.size(),dx.data()); el.derivativeR.emplace_back(dr.size(),dr.data());
            coef[i] = 0;
        }
        m_elements.push_back(std::move(el));
    }
}
int RANSWallNativeField::Size() const { return m_field->GetLocalToGlobalMap()->GetNumGlobalCoeffs(); }

RANSWallPoint RANSWallNativeField::ElementPoint(int element, double x, double r) const
{
    const auto &el = m_elements.at(element);
    Require(x >= el.left-1e-11 && x <= el.right+1e-11 && r >= el.bottom-1e-11 && r <= el.top+1e-11,
            "physical sampler outside its native element");
    const auto exp = m_field->GetExp(element); const auto map = m_field->GetLocalToGlobalMap();
    Array<OneD,NekDouble> local(2); local[0] = 2*(x-el.left)/(el.right-el.left)-1;
    local[1] = 2*(r-el.bottom)/(el.top-el.bottom)-1;
    RANSWallPoint row; row.x=x; row.r=r; row.element=element; row.geometry=el.geometry;
    for(int v=0;v<4;++v)
    {
        row.vertexPositions[v]=exp->GetVertexMap(v);
        row.vertexSigns[v]=map->GetLocalToGlobalSign(m_field->GetCoeff_Offset(element)+row.vertexPositions[v]);
    }
    for (int i = 0; i < exp->GetNcoeffs(); ++i)
    {
        const int pos = m_field->GetCoeff_Offset(element)+i;
        const double sign = map->GetLocalToGlobalSign(pos);
        row.ids.push_back(map->GetLocalToGlobalMap(pos));
        row.phi.push_back(sign*exp->PhysEvaluateBasis(local,i));
        row.dx.push_back(sign*exp->StdPhysEvaluate(local,el.derivativeX[i]));
        row.dr.push_back(sign*exp->StdPhysEvaluate(local,el.derivativeR[i]));
    }
    return row;
}
RANSWallPoint RANSWallNativeField::Point(double x, double r, int side) const
{
    int selected=-1;
    for (int e = 0; e < m_elements.size(); ++e)
    {
        const auto &el = m_elements[e];
        if (x < el.left-1e-11 || x > el.right+1e-11 || r < el.bottom-1e-11 || r > el.top+1e-11) continue;
        if ((side > 0 && r >= el.top-1e-12) || (side < 0 && r <= el.bottom+1e-12)) continue;
        if (selected < 0 || el.geometry < m_elements[selected].geometry) selected=e;
    }
    Require(selected >= 0, "physical sample not found on requested one-sided subdomain");
    return ElementPoint(selected,x,r);
}
std::vector<RANSWallPoint> RANSWallNativeField::Quadrature(int order) const
{
    std::vector<RANSWallPoint> result;
    for (int e = 0; e < m_elements.size(); ++e)
    {
        const auto &el = m_elements[e]; const auto exp = m_field->GetExp(e);
        const int n = order ? order : std::max(8,el.degree+3);
        auto rule = LibUtilities::PointsManager()[LibUtilities::PointsKey(n,LibUtilities::eGaussGaussLegendre)];
        const auto z=rule->GetZ(),w=rule->GetW();
        for (int j = 0; j < n; ++j) for (int i = 0; i < n; ++i)
        {
            auto row = ElementPoint(e,el.left+(z[i]+1)*(el.right-el.left)/2,
                                      el.bottom+(z[j]+1)*(el.top-el.bottom)/2);
            row.weight = w[i]*w[j]*(el.right-el.left)*(el.top-el.bottom)/4;
            result.push_back(std::move(row));
        }
    }
    return result;
}
RANSWallAlgebra::DenseMatrix RANSWallNativeField::Mass() const
{
    RANSWallAlgebra::DenseMatrix result(Size(),Vector(Size()));
    for (const auto &q : Quadrature())
        for (int i=0;i<q.ids.size();++i) for(int j=0;j<q.ids.size();++j)
            result[q.ids[i]][q.ids[j]] += q.weight*q.phi[i]*q.phi[j];
    return result;
}
Vector RANSWallNativeField::Project(const Function &f) const
{
    const auto quadrature=Quadrature();bool first=true,constant=true;double value=0;
    for (const auto &q : quadrature)
    {
        const double sampled=f(q.x,q.r);Require(std::isfinite(sampled),"nonfinite projection data");
        if(first){value=sampled;first=false;}else if(sampled!=value)constant=false;
    }
    if(constant)
    {
        Vector result(Size());const auto map=m_field->GetLocalToGlobalMap();
        for(int e=0;e<m_field->GetExpSize();++e)for(int v=0;v<4;++v)
        {
            const int pos=m_field->GetCoeff_Offset(e)+m_field->GetExp(e)->GetVertexMap(v);
            result[map->GetLocalToGlobalMap(pos)]=value*map->GetLocalToGlobalSign(pos);
        }
        return result;
    }
    // A local native transform supplies a good seed. Solve the physical L2
    // defect, rather than forming a large background load and inverting its
    // rounded mass matrix. The correction retains the same global L2 problem.
    Vector result(Size()),count(Size());const auto map=m_field->GetLocalToGlobalMap();
    for(int e=0;e<m_field->GetExpSize();++e)
    {
        const auto exp=m_field->GetExp(e);const int nq=exp->GetTotPoints();
        Array<OneD,NekDouble> x(nq),r(nq),z(nq),values(nq),coeff(exp->GetNcoeffs());
        exp->GetCoords(x,r,z);
        for(int q=0;q<nq;++q)values[q]=f(x[q],r[q]);
        exp->FwdTrans(values,coeff);
        for(int j=0;j<coeff.size();++j)
        {
            const int pos=m_field->GetCoeff_Offset(e)+j,id=map->GetLocalToGlobalMap(pos);
            result[id]+=map->GetLocalToGlobalSign(pos)*coeff[j];count[id]+=1;
        }
    }
    for(int j=0;j<result.size();++j){Require(count[j]>0,"unmapped projection coefficient");result[j]/=count[j];}
    const auto mass=Mass();
    for(int iteration=0;iteration<3;++iteration)
    {
        std::vector<long double> load(Size(),0);
        for(const auto &q:quadrature)
        {
            const double defect=f(q.x,q.r)-q.Value(result);
            for(int j=0;j<q.ids.size();++j)load[q.ids[j]]+=(long double)q.weight*q.phi[j]*defect;
        }
        Vector rhs(load.begin(),load.end());const auto correction=Solve(mass,rhs);
        for(int j=0;j<result.size();++j)result[j]+=correction[j];
    }
    return result;
}
std::map<int,double> RANSWallNativeField::Boundary(const Function &f,
    const std::function<bool(int,const Element &)> &select) const
{
    std::map<int,double> result;
    const auto map=m_field->GetLocalToGlobalMap();
    for(int e=0;e<m_elements.size();++e)
    {
        const auto &el=m_elements[e]; const int modes=el.degree+1;
        StdRegions::StdSegExp seg(LibUtilities::BasisKey(LibUtilities::eModified_A,
            modes,LibUtilities::PointsKey(std::max(8,modes+2),LibUtilities::eGaussLobattoLegendre)));
        const auto z=seg.GetBasis(0)->GetZ();
        Array<OneD,NekDouble> v(seg.GetTotPoints()),coeff(modes);
        for(int edge=0;edge<4;++edge) if(select(edge,el))
        {
            for(int q=0;q<v.size();++q)
            {
                double x=el.left+(z[q]+1)*(el.right-el.left)/2,
                       r=el.bottom+(z[q]+1)*(el.top-el.bottom)/2;
                if(edge==0) r=el.bottom;
                if(edge==1) x=el.right;
                if(edge==2) r=el.top;
                if(edge==3) x=el.left;
                v[q]=f(x,r);
            }
            seg.FwdTransBndConstrained(v,coeff);
            for(int j=0;j<modes;++j)
            {
                const int mode=edge==0?j:edge==1?1+modes*j:edge==2?j+modes:modes*j;
                const int pos=m_field->GetCoeff_Offset(e)+mode, id=map->GetLocalToGlobalMap(pos);
                const double value=map->GetLocalToGlobalSign(pos)*coeff[j];
                if(result.count(id)) Require(std::abs(result[id]-value)<1e-9*std::max(1.,std::abs(value)),"inconsistent corner/edge Dirichlet data");
                result[id]=value;
            }
        }
    }
    return result;
}
void RANSWallNativeField::Store(const Vector &state,int offset) const
{
    Array<OneD,NekDouble> global(Size());
    for(int i=0;i<Size();++i) global[i]=state.at(offset+i);
    m_field->GetLocalToGlobalMap()->GlobalToLocal(global,m_field->UpdateCoeffs());
    m_field->BwdTrans(m_field->GetCoeffs(),m_field->UpdatePhys()); m_field->SetPhysState(true);
}
Vector RANSWallNativeField::NativeLoad(const Function &f) const
{
    const int nq=m_field->GetTotPoints();
    Array<OneD,NekDouble> x(nq),r(nq),z(nq),values(nq),load(m_field->GetNcoeffs()),global(Size());
    m_field->GetCoords(x,r,z);
    for(int i=0;i<nq;++i) values[i]=f(x[i],r[i]);
    m_field->IProductWRTBase(values,load);m_field->GetLocalToGlobalMap()->Assemble(load,global);
    return {global.begin(),global.end()};
}

RANSWallFieldLayout::RANSWallFieldLayout(const LibUtilities::SessionReaderSharedPtr &session,
    const SpatialDomains::MeshGraphSharedPtr &graph,double height) : m_height(height)
{
    Require(height>0 && session->GetComm()->GetSize()==1,"M3 flat layout is serial with positive H");
    auto full=MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(session,graph,true,"u");
    auto pressure=MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(session,graph,true,"p");
    std::vector<unsigned int> all,inner,outer;
    m_left=std::numeric_limits<double>::infinity();m_right=-m_left;m_top=0;
    for(int e=0;e<full->GetExpSize();++e)
    {
        const auto box=PhysicalBounds(full->GetExp(e)->GetGeom());all.push_back(e);
        m_left=std::min(m_left,box[0]);m_right=std::max(m_right,box[3]);m_top=std::max(m_top,box[4]);
        Require(box[1]>=-1e-12 && !(box[1]<height-1e-12 && box[4]>height+1e-12),"element crosses ground or matching height");
        if(box[4]<=height+1e-12) inner.push_back(e); else outer.push_back(e);
        Require(pressure->GetExp(e)->GetGeom()->GetGlobalID()==full->GetExp(e)->GetGeom()->GetGlobalID(),"pressure/velocity element order mismatch");
    }
    velocity=std::make_shared<RANSWallNativeField>(session,full,all);
    innerVelocity=std::make_shared<RANSWallNativeField>(session,full,inner);
    outerVelocity=std::make_shared<RANSWallNativeField>(session,full,outer);
    k=std::make_shared<RANSWallNativeField>(session,full,all);
    epsilon=std::make_shared<RANSWallNativeField>(session,full,outer);
    innerPressure=std::make_shared<RANSWallNativeField>(session,pressure,inner);
    outerPressure=std::make_shared<RANSWallNativeField>(session,pressure,outer);
    for(int b=0;b<Gauge;++b) size[b]=Field(static_cast<Block>(b)).Size();
    size[Gauge]=1; for(int b=1;b<Blocks;++b)offset[b]=offset[b-1]+size[b-1];
}
const RANSWallNativeField &RANSWallFieldLayout::Field(Block b) const
{
    switch(b)
    {
        case U:case W:return *velocity;case K:return *k;case Epsilon:return *epsilon;
        case InnerPressure:return *innerPressure;case OuterPressure:return *outerPressure;
        default:throw std::invalid_argument("Gauge has no expansion");
    }
}
}
