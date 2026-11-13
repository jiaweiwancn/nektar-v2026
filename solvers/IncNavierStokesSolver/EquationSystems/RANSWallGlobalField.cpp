///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallGlobalField.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <StdRegions/StdSegExp.h>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <set>
#include <stdexcept>

namespace Nektar
{
namespace
{
using Vector = RANSWallAlgebra::Vector;
using WallMatrix = RANSWallAlgebra::DenseMatrix;
using Coordinate = RANSWallCoordinate;
void Require(bool valid, const std::string &message)
{
    if (!valid) throw std::runtime_error("Global wall field: " + message);
}
bool Equal(double a,double b)
{
    return std::abs(a-b) <= 1e-11;
}
Coordinate Physical(double x,double y,double z,int dimension)
{
    return dimension==2 ? Coordinate{x,0,y} : Coordinate{x,y,z};
}
double Sample(const RANSWallFunction &function,const Coordinate &coordinate)
{
    const double value=function(coordinate);
    Require(std::isfinite(value),"nonfinite projection/boundary data");
    return value;
}
int TensorSize(int modes,int dimension)
{
    return dimension==2 ? modes*modes : modes*modes*modes;
}
std::array<int,3> Indices(int index,int modes)
{
    return {index%modes,(index/modes)%modes,index/(modes*modes)};
}
}

struct RANSWallGlobalField::BasisData
{
    struct Values { Vector value,derivative,secondDerivative; };
    std::shared_ptr<StdRegions::StdSegExp> segment;
    std::vector<Array<OneD,NekDouble>> derivative;
    std::vector<Array<OneD,NekDouble>> secondDerivative;
    mutable std::map<double,std::shared_ptr<const Values>> cache;

    explicit BasisData(const LibUtilities::BasisKey &key)
        :segment(std::make_shared<StdRegions::StdSegExp>(key))
    {
        const int points=segment->GetTotPoints();
        const auto db=segment->GetBasis(0)->GetDbdata();
        for(int j=0;j<segment->GetNcoeffs();++j)
        {
            derivative.emplace_back(points,db.data()+j*points);
            secondDerivative.emplace_back(points,0.);
            segment->StdPhysDeriv(derivative.back(),secondDerivative.back());
        }
    }
    std::shared_ptr<const Values> Evaluate(double coordinate) const
    {
        coordinate=std::clamp(coordinate,-1.,1.);
        const auto found=cache.find(coordinate);
        if(found!=cache.end()) return found->second;
        auto result=std::make_shared<Values>();
        Array<OneD,NekDouble> local(1);local[0]=coordinate;
        for(int j=0;j<segment->GetNcoeffs();++j)
        {
            result->value.push_back(segment->PhysEvaluateBasis(local,j));
            result->derivative.push_back(segment->StdPhysEvaluate(local,derivative[j]));
            result->secondDerivative.push_back(segment->StdPhysEvaluate(local,secondDerivative[j]));
        }
        // Moving cutoffs can create many radial coordinates during a Jacobian.
        // Bound the cache; shared results stay valid through the current row.
        if(cache.size()>=8192) cache.clear();
        cache.emplace(coordinate,result);
        return result;
    }
};

double RANSWallGlobalPoint::Value(const Vector &state,int offset) const
{
    Require(offset>=0 && !ids.empty() && constant.size()==ids.size() &&
                phi.size()==ids.size() && constant[0]!=0,"invalid value row");
    const double centre=state.at(offset+ids[0])/constant[0];
    long double value=centre;
    for(std::size_t j=0;j<ids.size();++j)
        value+=(long double)phi[j]*(state.at(offset+ids[j])-centre*constant[j]);
    return value;
}
std::array<double,3> RANSWallGlobalPoint::Gradient(const Vector &state,int offset) const
{
    Require(offset>=0 && !ids.empty() && constant.size()==ids.size() &&
                constant[0]!=0,"invalid gradient row");
    const double centre=state.at(offset+ids[0])/constant[0];
    std::array<long double,3> value{};
    for(int d=0;d<3;++d)
    {
        Require(derivative[d].size()==ids.size(),"invalid derivative row");
        for(std::size_t j=0;j<ids.size();++j)
            value[d]+=(long double)derivative[d][j]*(state.at(offset+ids[j])-centre*constant[j]);
    }
    return {static_cast<double>(value[0]),static_cast<double>(value[1]),static_cast<double>(value[2])};
}

RANSWallGlobalField::RANSWallGlobalField(
    const LibUtilities::SessionReaderSharedPtr &session,
    const MultiRegions::ExpListSharedPtr &full,
    const std::vector<unsigned int> &elements,int dimension)
    :m_dimension(dimension),m_axes(dimension==2?std::array<int,3>{0,2,-1}:
                                                 std::array<int,3>{0,1,2})
{
    Require((dimension==2 || dimension==3) && session->GetComm()->GetSize()==1 &&
                !elements.empty(),"requires nonempty serial 2D/3D native field");
    auto subset=MemoryManager<MultiRegions::ExpList>::AllocateSharedPtr(*full,elements,true);
    m_field=MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(session,*subset);
    for(int e=0;e<m_field->GetExpSize();++e)
    {
        const auto exp=m_field->GetExp(e);const auto geometry=exp->GetGeom();
        Require(geometry->GetCoordim()==dimension &&
                    geometry->GetShapeType()==(dimension==2?LibUtilities::eQuadrilateral:LibUtilities::eHexahedron),
                "requires physical quads/hexahedra of the declared dimension");
        Require(exp->GetGeomFactors()->GetGtype()==SpatialDomains::eRegular,
                "requires affine Cartesian geometry");
        Element element; element.geometry=geometry->GetGlobalID();
        element.degree=exp->GetBasisNumModes(0)-1;
        Require(element.degree>=1,"Modified_A degree must be at least one");
        element.lower.fill(std::numeric_limits<double>::infinity());
        element.upper.fill(-std::numeric_limits<double>::infinity());
        for(int v=0;v<geometry->GetNumVerts();++v)
        {
            double x=0,y=0,z=0;geometry->GetVertex(v)->GetCoords(x,y,z);
            const auto point=Physical(x,y,z,dimension);
            for(int d=0;d<3;++d)
            {
                element.lower[d]=std::min(element.lower[d],point[d]);
                element.upper[d]=std::max(element.upper[d],point[d]);
            }
        }
        std::array<std::shared_ptr<BasisData>,3> bases{};
        for(int d=0;d<dimension;++d)
        {
            const auto basis=exp->GetBasis(d);
            Require(basis->GetBasisType()==LibUtilities::eModified_A &&
                        basis->GetNumModes()==element.degree+1 &&
                        basis->GetNumPoints()>=basis->GetNumModes(),
                    "requires equal-order Modified_A tensor bases");
            Require(element.upper[m_axes[d]]>element.lower[m_axes[d]],"nonpositive physical element extent");
            const auto key=std::make_tuple(basis->GetNumModes(),basis->GetNumPoints(),
                                           static_cast<int>(basis->GetPointsType()));
            auto found=m_basis.find(key);
            if(found==m_basis.end())
                found=m_basis.emplace(key,std::make_shared<BasisData>(basis->GetBasisKey())).first;
            bases[d]=found->second;
        }
        const int modes=element.degree+1;
        Require(exp->GetNcoeffs()==TensorSize(modes,dimension),"unexpected tensor coefficient count");
        for(int v=0;v<geometry->GetNumVerts();++v)
        {
            // Native geometric vertex ordering differs from tensor packing.
            const auto tensor=Indices(exp->GetVertexMap(v,false),modes);
            double x=0,y=0,z=0;geometry->GetVertex(v)->GetCoords(x,y,z);
            const auto point=Physical(x,y,z,dimension);
            for(int d=0;d<dimension;++d)
            {
                Require(tensor[d]==0 || tensor[d]==1,"unexpected native vertex packing");
                const double expected=tensor[d]?element.upper[m_axes[d]]:element.lower[m_axes[d]];
                Require(Equal(point[m_axes[d]],expected),"requires canonical axis-aligned element vertices");
            }
        }
        m_elements.push_back(element);m_elementBasis.push_back(std::move(bases));
    }
}
int RANSWallGlobalField::Size() const
{
    return m_field->GetLocalToGlobalMap()->GetNumGlobalCoeffs();
}
RANSWallGlobalPoint RANSWallGlobalField::ElementPoint(int e,const Coordinate &coordinate) const
{
    const auto &element=m_elements.at(e);
    Require(std::all_of(coordinate.begin(),coordinate.end(),[](double x){return std::isfinite(x);}),
            "nonfinite physical sample");
    if(m_dimension==2) Require(Equal(coordinate[1],0),"2D sample has nonzero physical span coordinate");
    std::array<std::shared_ptr<const BasisData::Values>,3> basis{};
    std::array<double,3> inverse{};
    for(int d=0;d<m_dimension;++d)
    {
        const int a=m_axes[d];
        Require(coordinate[a]>=element.lower[a]-1e-11 && coordinate[a]<=element.upper[a]+1e-11,
                "sample outside its native element");
        inverse[d]=2/(element.upper[a]-element.lower[a]);
        basis[d]=m_elementBasis[e][d]->Evaluate((coordinate[a]-element.lower[a])*inverse[d]-1);
    }
    const int modes=element.degree+1,count=TensorSize(modes,m_dimension);
    const auto map=m_field->GetLocalToGlobalMap();const int start=m_field->GetCoeff_Offset(e);
    RANSWallGlobalPoint result;result.coordinate=coordinate;result.geometry=element.geometry;result.element=e;
    result.ids.reserve(count);result.phi.reserve(count);result.constant.reserve(count);
    for(auto &derivative:result.derivative) derivative.reserve(count);
    for(int j=0;j<count;++j)
    {
        const auto index=Indices(j,modes);
        const double sign=map->GetLocalToGlobalSign(start+j);
        result.ids.push_back(map->GetLocalToGlobalMap(start+j));
        double value=sign,constant=sign;
        for(int d=0;d<m_dimension;++d)
        {
            value*=basis[d]->value[index[d]];
            if(index[d]>=2) constant=0;
        }
        result.phi.push_back(value);result.constant.push_back(constant);
        std::array<double,3> derivative{};
        for(int d=0;d<m_dimension;++d)
        {
            double entry=sign*inverse[d];
            for(int q=0;q<m_dimension;++q)
                entry*=q==d?basis[q]->derivative[index[q]]:basis[q]->value[index[q]];
            derivative[m_axes[d]]=entry;
        }
        for(int d=0;d<3;++d) result.derivative[d].push_back(derivative[d]);
    }
    return result;
}
std::array<RANSWallGlobalField::Vector,3> RANSWallGlobalField::RadialGradientDerivative(
    const RANSWallGlobalPoint &point) const
{
    const auto &element=m_elements.at(point.element);
    Require(point.geometry==element.geometry,"radial Hessian has mismatched native geometry");
    const int radial=m_dimension-1,modes=element.degree+1,count=TensorSize(modes,m_dimension);
    Require(point.ids.size()==static_cast<std::size_t>(count),"radial Hessian row shape differs");
    std::array<std::shared_ptr<const BasisData::Values>,3> basis{};
    std::array<double,3> inverse{};
    for(int d=0;d<m_dimension;++d)
    {
        const int axis=m_axes[d];inverse[d]=2/(element.upper[axis]-element.lower[axis]);
        basis[d]=m_elementBasis[point.element][d]->Evaluate(
            (point.coordinate[axis]-element.lower[axis])*inverse[d]-1);
    }
    std::array<Vector,3> result;for(auto &row:result)row.assign(count,0.);
    const auto map=m_field->GetLocalToGlobalMap();const int start=m_field->GetCoeff_Offset(point.element);
    for(int j=0;j<count;++j)
    {
        const auto index=Indices(j,modes);
        for(int d=0;d<m_dimension;++d)
        {
            double value=map->GetLocalToGlobalSign(start+j)*inverse[d]*inverse[radial];
            for(int a=0;a<m_dimension;++a)
                value*=a==radial && a==d?basis[a]->secondDerivative[index[a]]:
                    a==radial || a==d?basis[a]->derivative[index[a]]:basis[a]->value[index[a]];
            result[m_axes[d]][j]=value;
        }
    }
    return result;
}
RANSWallGlobalPoint RANSWallGlobalField::Point(const Coordinate &coordinate,int side) const
{
    return ElementPoint(LocateElement(coordinate,side),coordinate);
}
int RANSWallGlobalField::LocateElement(const Coordinate &coordinate,int side) const
{
    Require(side>=-1 && side<=1,"unknown radial sampling side");
    int selected=-1;
    for(std::size_t e=0;e<m_elements.size();++e)
    {
        const auto &element=m_elements[e];bool contains=true;
        for(int d=0;d<m_dimension;++d)
        {
            const int a=m_axes[d];
            contains=contains && coordinate[a]>=element.lower[a]-1e-11 && coordinate[a]<=element.upper[a]+1e-11;
        }
        if(!contains || (side>0 && coordinate[2]>=element.upper[2]-1e-12) ||
                        (side<0 && coordinate[2]<=element.lower[2]+1e-12)) continue;
        if(selected<0 || element.geometry<m_elements[selected].geometry) selected=e;
    }
    Require(selected>=0,"physical sample not found on requested radial side");
    return selected;
}
std::vector<RANSWallGlobalPoint> RANSWallGlobalField::Quadrature(int order) const
{
    Require(order>=0,"negative quadrature order");
    std::vector<RANSWallGlobalPoint> result;
    for(std::size_t e=0;e<m_elements.size();++e)
    {
        const auto &element=m_elements[e];const int n=order?order:std::max(8,element.degree+3);
        auto rule=LibUtilities::PointsManager()[LibUtilities::PointsKey(n,LibUtilities::eGaussGaussLegendre)];
        const auto z=rule->GetZ(),weights=rule->GetW();
        for(int k=0;k<(m_dimension==3?n:1);++k)for(int j=0;j<n;++j)for(int i=0;i<n;++i)
        {
            const std::array<int,3> index{i,j,k};Coordinate coordinate{};double weight=1;
            for(int d=0;d<m_dimension;++d)
            {
                const int a=m_axes[d];const double half=(element.upper[a]-element.lower[a])/2;
                coordinate[a]=element.lower[a]+(z[index[d]]+1)*half;
                weight*=weights[index[d]]*half;
            }
            auto point=ElementPoint(e,coordinate);point.weight=weight;result.push_back(std::move(point));
        }
    }
    return result;
}
RANSWallGlobalField::WallMatrix RANSWallGlobalField::Mass() const
{
    std::vector<std::vector<long double>> mass(Size(),std::vector<long double>(Size()));
    for(const auto &q:Quadrature())
        for(std::size_t i=0;i<q.ids.size();++i)for(std::size_t j=0;j<q.ids.size();++j)
            mass[q.ids[i]][q.ids[j]]+=(long double)q.weight*q.phi[i]*q.phi[j];
    WallMatrix result(Size(),Vector(Size()));
    for(int i=0;i<Size();++i)for(int j=0;j<Size();++j)result[i][j]=mass[i][j];
    return result;
}
RANSWallGlobalField::Vector RANSWallGlobalField::Project(const RANSWallFunction &function) const
{
    Require(static_cast<bool>(function),"missing projection function");
    const auto quadrature=Quadrature();bool first=true,isConstant=true;double constant=0;
    for(const auto &q:quadrature)
    {
        const double value=Sample(function,q.coordinate);
        if(first){first=false;constant=value;}else if(value!=constant)isConstant=false;
    }
    const auto map=m_field->GetLocalToGlobalMap();Vector result(Size());
    if(isConstant)
    {
        for(int e=0;e<m_field->GetExpSize();++e)
            for(int v=0;v<m_field->GetExp(e)->GetGeom()->GetNumVerts();++v)
            {
                const int pos=m_field->GetCoeff_Offset(e)+m_field->GetExp(e)->GetVertexMap(v,false);
                result[map->GetLocalToGlobalMap(pos)]=map->GetLocalToGlobalSign(pos)*constant;
            }
        return result;
    }
    Vector count(Size());
    for(int e=0;e<m_field->GetExpSize();++e)
    {
        const auto exp=m_field->GetExp(e);const int n=exp->GetTotPoints();
        Array<OneD,NekDouble> x(n,0.),y(n,0.),z(n,0.),values(n),coeff(exp->GetNcoeffs());
        exp->GetCoords(x,y,z);
        for(int q=0;q<n;++q)
        {
            auto point=Physical(x[q],y[q],z[q],m_dimension);
            // Native GetCoords reconstructs an affine canonical quad/hex by
            // floating-point basis evaluation. Its endpoint can be a few ulps
            // outside the exact element and hence outside a certified map.
            // Canonicalize only exterior roundoff, never an interior sample or
            // an actual out-of-element request. Physical geometry guards stay
            // strict; no field value or constitutive state is clipped.
            const auto &element=m_elements[e];
            for(int axis=0;axis<3;++axis)if(m_dimension==3 || axis!=1)
            {
                const double lo=element.lower[axis],hi=element.upper[axis];
                const double tol=16*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(lo),std::abs(hi)});
                if(point[axis]<lo && lo-point[axis]<=tol)point[axis]=lo;
                if(point[axis]>hi && point[axis]-hi<=tol)point[axis]=hi;
                Require(point[axis]>=lo && point[axis]<=hi,"native projection coordinate outside canonical affine element");
            }
            values[q]=Sample(function,point);
        }
        exp->FwdTrans(values,coeff);
        for(int j=0;j<exp->GetNcoeffs();++j)
        {
            const int pos=m_field->GetCoeff_Offset(e)+j,id=map->GetLocalToGlobalMap(pos);
            result[id]+=map->GetLocalToGlobalSign(pos)*coeff[j];count[id]+=1;
        }
    }
    for(int j=0;j<Size();++j){Require(count[j]>0,"unmapped projection coefficient");result[j]/=count[j];}
    const auto mass=Mass();
    for(int iteration=0;iteration<3;++iteration)
    {
        std::vector<long double> load(Size());
        for(const auto &q:quadrature)
        {
            const double defect=Sample(function,q.coordinate)-q.Value(result);
            for(std::size_t j=0;j<q.ids.size();++j)load[q.ids[j]]+=(long double)q.weight*q.phi[j]*defect;
        }
        const auto correction=RANSWallAlgebra::Solve(mass,Vector(load.begin(),load.end()));
        for(int j=0;j<Size();++j)result[j]+=correction[j];
    }
    return result;
}
std::map<int,double> RANSWallGlobalField::Boundary(
    const RANSWallFunction &function,const BoundaryPredicate &predicate) const
{
    Require(static_cast<bool>(function) && static_cast<bool>(predicate),"missing boundary function/predicate");
    std::map<int,double> result;const auto map=m_field->GetLocalToGlobalMap();
    for(std::size_t e=0;e<m_elements.size();++e)
    {
        const auto &element=m_elements[e];const int modes=element.degree+1,n=std::max(8,modes+2);
        const LibUtilities::BasisKey key(LibUtilities::eModified_A,modes,
            LibUtilities::PointsKey(n,LibUtilities::eGaussLobattoLegendre));
        StdRegions::StdSegExp segment(key);
        const auto z=segment.GetBasis(0)->GetZ();
        for(int fixed=0;fixed<m_dimension;++fixed)for(int side=0;side<2;++side)
        {
            if(!predicate(m_axes[fixed],side,element))continue;
            std::vector<int> free;for(int d=0;d<m_dimension;++d)if(d!=fixed)free.push_back(d);
            const int nq=m_dimension==2?n:n*n,nc=m_dimension==2?modes:modes*modes;
            Array<OneD,NekDouble> values(nq),coefficients(nc);
            for(int q=0;q<nq;++q)
            {
                const std::array<int,2> index{q%n,q/n};Coordinate coordinate{};
                coordinate[m_axes[fixed]]=side?element.upper[m_axes[fixed]]:element.lower[m_axes[fixed]];
                for(std::size_t d=0;d<free.size();++d)
                {
                    const int a=m_axes[free[d]];
                    coordinate[a]=element.lower[a]+(z[index[d]]+1)*(element.upper[a]-element.lower[a])/2;
                }
                values[q]=Sample(function,coordinate);
            }
            if(m_dimension==2)segment.FwdTransBndConstrained(values,coefficients);
            else
            {
                // The native StdQuadExp constrained transform reverses the
                // top/left physical edges although its trace map uses increasing
                // tensor order. Two commuting native segment projections give
                // consistent edges/corners and the same interior-bubble L2
                // equations on this affine tensor face.
                Array<OneD,NekDouble> temporary(modes*n),physical(n),modal(modes);
                for(int row=0;row<n;++row)
                {
                    for(int i=0;i<n;++i)physical[i]=values[i+n*row];
                    segment.FwdTransBndConstrained(physical,modal);
                    for(int i=0;i<modes;++i)temporary[i+modes*row]=modal[i];
                }
                for(int i=0;i<modes;++i)
                {
                    for(int row=0;row<n;++row)physical[row]=temporary[i+modes*row];
                    segment.FwdTransBndConstrained(physical,modal);
                    for(int j=0;j<modes;++j)coefficients[i+modes*j]=modal[j];
                }
            }
            for(int j=0;j<nc;++j)
            {
                const std::array<int,2> index{j%modes,j/modes};std::array<int,3> tensor{};tensor[fixed]=side;
                for(std::size_t d=0;d<free.size();++d)tensor[free[d]]=index[d];
                const int mode=tensor[0]+modes*tensor[1]+modes*modes*tensor[2];
                const int pos=m_field->GetCoeff_Offset(e)+mode,id=map->GetLocalToGlobalMap(pos);
                const double value=map->GetLocalToGlobalSign(pos)*coefficients[j];
                const auto found=result.find(id);
                Require(found==result.end() || std::abs(found->second-value)<1e-9*std::max(1.,std::abs(value)),
                        "inconsistent edge/corner boundary data, geometry="+std::to_string(element.geometry)+
                        ", axis="+std::to_string(m_axes[fixed])+", side="+std::to_string(side)+
                        ", global="+std::to_string(id)+", previous="+(found==result.end()?std::string("unset"):std::to_string(found->second))+
                        ", current="+std::to_string(value));
                result[id]=value;
            }
        }
    }
    return result;
}
void RANSWallGlobalField::Store(const Vector &state,int offset) const
{
    Require(offset>=0 && static_cast<std::size_t>(offset)+Size()<=state.size(),"invalid stored coefficient shape");
    Array<OneD,NekDouble> global(Size());for(int j=0;j<Size();++j)global[j]=state[offset+j];
    m_field->GetLocalToGlobalMap()->GlobalToLocal(global,m_field->UpdateCoeffs());
    m_field->BwdTrans(m_field->GetCoeffs(),m_field->UpdatePhys());m_field->SetPhysState(true);
}
std::array<double,4> RANSWallGlobalField::NativeEvaluate(
    const Vector &state,const Coordinate &coordinate,int offset,int side) const
{
    Require(offset>=0 && static_cast<std::size_t>(offset)+Size()<=state.size(),"invalid native probe coefficient shape");
    const auto point=Point(coordinate,side);const auto exp=m_field->GetExp(point.element);
    Array<OneD,NekDouble> coefficients(exp->GetNcoeffs()),values(exp->GetTotPoints()),
        derivative(exp->GetTotPoints()),local(m_dimension);
    for(int d=0;d<m_dimension;++d)
    {
        const int a=m_axes[d];const auto &element=m_elements[point.element];
        local[d]=2*(coordinate[a]-element.lower[a])/(element.upper[a]-element.lower[a])-1;
    }
    for(std::size_t j=0;j<point.ids.size();++j)
        coefficients[j]=m_field->GetLocalToGlobalMap()->GetLocalToGlobalSign(m_field->GetCoeff_Offset(point.element)+j)*
                        state.at(offset+point.ids[j]);
    exp->BwdTrans(coefficients,values);std::array<double,4> result{};
    result[0]=exp->StdPhysEvaluate(local,values);
    for(int d=0;d<m_dimension;++d)
    {
        exp->PhysDeriv(d,values,derivative);
        result[1+m_axes[d]]=exp->StdPhysEvaluate(local,derivative);
    }
    return result;
}

// This independent native transform deliberately retains the same coefficient
// ordering and sign map while requiring explicit element ownership.
std::array<double,4> RANSWallGlobalField::NativeEvaluateElement(
    const Vector &state,int element,const Coordinate &coordinate,int offset) const
{
    Require(offset>=0 && static_cast<std::size_t>(offset)+Size()<=state.size(),"invalid native probe coefficient shape");
    const auto point=ElementPoint(element,coordinate);const auto exp=m_field->GetExp(point.element);
    Array<OneD,NekDouble> coefficients(exp->GetNcoeffs()),values(exp->GetTotPoints()),
        derivative(exp->GetTotPoints()),local(m_dimension);
    for(int d=0;d<m_dimension;++d)
    {
        const int a=m_axes[d];const auto &element=m_elements[point.element];
        local[d]=2*(coordinate[a]-element.lower[a])/(element.upper[a]-element.lower[a])-1;
    }
    for(std::size_t j=0;j<point.ids.size();++j)
        coefficients[j]=m_field->GetLocalToGlobalMap()->GetLocalToGlobalSign(m_field->GetCoeff_Offset(point.element)+j)*
                        state.at(offset+point.ids[j]);
    exp->BwdTrans(coefficients,values);std::array<double,4> result{};
    result[0]=exp->StdPhysEvaluate(local,values);
    for(int d=0;d<m_dimension;++d)
    {
        exp->PhysDeriv(d,values,derivative);
        result[1+m_axes[d]]=exp->StdPhysEvaluate(local,derivative);
    }
    return result;
}

RANSWallGlobalField::Vector RANSWallGlobalField::RadialProject(
    const std::function<double(double)> &function,
    const std::vector<std::array<double,3>> &constraints,const std::vector<double> &physicalBreakpoints) const
{
    Require(static_cast<bool>(function),"missing radial projection function");
    using IntervalKey=std::pair<double,double>;
    struct Interval { double lower,upper,area=0;int degree;std::vector<int> ids;std::size_t element; };
    std::map<IntervalKey,Interval> intervals;
    for(std::size_t e=0;e<m_elements.size();++e)
    {
        const auto &element=m_elements[e];const IntervalKey key{element.lower[2],element.upper[2]};
        double area=element.upper[0]-element.lower[0];
        if(m_dimension==3)area*=element.upper[1]-element.lower[1];
        auto found=intervals.find(key);
        if(found==intervals.end())intervals.emplace(key,Interval{key.first,key.second,area,element.degree,{},e});
        else
        {
            Require(found->second.degree==element.degree,"nonuniform horizontal degree in radial restriction");
            found->second.area+=area;
        }
    }
    Require(!intervals.empty(),"empty radial partition");
    bool first=true;double previous=0;
    using Key=std::tuple<double,double,int>;std::map<Key,int> radial;
    auto number=[&](const Key &key)
    {
        auto found=radial.find(key);
        if(found!=radial.end())return found->second;
        const int id=radial.size();radial.emplace(key,id);return id;
    };
    for(auto &item:intervals)
    {
        auto &interval=item.second;
        Require(first || Equal(previous,interval.lower),"radial intervals overlap or do not form a conforming chain");
        first=false;previous=interval.upper;
        interval.ids.push_back(number({interval.lower,interval.lower,0}));
        interval.ids.push_back(number({interval.upper,interval.upper,0}));
        for(int j=2;j<=interval.degree;++j)interval.ids.push_back(number({interval.lower,interval.upper,j}));
    }
    const int n=radial.size();std::vector<int> injection(Size(),-1);Vector sign(Size());
    const auto map=m_field->GetLocalToGlobalMap();
    for(std::size_t e=0;e<m_elements.size();++e)
    {
        const auto &element=m_elements[e];const int modes=element.degree+1;
        const auto &interval=intervals.at({element.lower[2],element.upper[2]});
        for(int j=0;j<TensorSize(modes,m_dimension);++j)
        {
            const auto index=Indices(j,modes);bool active=true;
            for(int d=0;d<m_dimension-1;++d)active=active && index[d]<2;
            if(!active)continue;
            const int pos=m_field->GetCoeff_Offset(e)+j,id=map->GetLocalToGlobalMap(pos),
                      radialId=interval.ids[index[m_dimension-1]];
            const double orientation=map->GetLocalToGlobalSign(pos);
            Require(injection[id]<0 || (injection[id]==radialId && sign[id]==orientation),
                    "incompatible native radial injection");
            injection[id]=radialId;sign[id]=orientation;
        }
    }
    WallMatrix mass(n,Vector(n));std::vector<long double> load(n);bool constant=true,firstSample=true;double value=0;
    auto rule=LibUtilities::PointsManager()[LibUtilities::PointsKey(32,LibUtilities::eGaussGaussLegendre)];
    const auto gauss=rule->GetZ(),weights=rule->GetW();
    for(const auto &item:intervals)
    {
        const auto &interval=item.second;const auto basis=m_elementBasis[interval.element][m_dimension-1];
        std::set<double> panels{interval.lower,interval.upper};
        for(double r:physicalBreakpoints)
        {Require(std::isfinite(r),"nonfinite radial integration breakpoint");if(r>interval.lower && r<interval.upper)panels.insert(r);}
        double r=1e-4;int guard=0;
        while(r<interval.upper)
        {
            if(r>interval.lower)panels.insert(r);
            r*=1.5;Require(++guard<10000,"invalid radial quadrature range");
        }
        for(auto lower=panels.begin(),upper=std::next(lower);upper!=panels.end();++lower,++upper)
            for(int q=0;q<gauss.size();++q)
            {
                const double radius=*lower+(gauss[q]+1)*(*upper-*lower)/2,
                    weight=weights[q]*(*upper-*lower)/2*interval.area;
                const double sampled=function(radius);Require(std::isfinite(sampled),"nonfinite radial projection data");
                if(firstSample){value=sampled;firstSample=false;}else if(sampled!=value)constant=false;
                const auto row=basis->Evaluate(2*(radius-interval.lower)/(interval.upper-interval.lower)-1);
                for(std::size_t i=0;i<interval.ids.size();++i)
                {
                    const int a=interval.ids[i];load[a]+=(long double)weight*row->value[i]*sampled;
                    for(std::size_t j=0;j<interval.ids.size();++j)
                        mass[a][interval.ids[j]]+=weight*row->value[i]*row->value[j];
                }
            }
    }
    WallMatrix rows;Vector targets;
    for(const auto &constraint:constraints)
    {
        const double height=constraint[0],side=constraint[1];
        Require(std::isfinite(height) && std::isfinite(constraint[2]) && (side==0 || side==-1 || side==1),
                "invalid radial projection constraint");
        const Element *selected=nullptr;
        for(const auto &element:m_elements)
        {
            if(height<element.lower[2]-1e-11 || height>element.upper[2]+1e-11 ||
                (side<0 && height<=element.lower[2]+1e-12) ||
                (side>0 && height>=element.upper[2]-1e-12))continue;
            if(!selected || element.geometry<selected->geometry)selected=&element;
        }
        Require(selected,"radial constraint outside requested one-sided domain");
        Coordinate coordinate{};
        for(int d=0;d<m_dimension-1;++d)
        {
            const int a=m_axes[d];coordinate[a]=(selected->lower[a]+selected->upper[a])/2;
        }
        coordinate[2]=height;const auto point=Point(coordinate,static_cast<int>(side));
        Vector row(n);
        for(std::size_t j=0;j<point.ids.size();++j)
        {
            const int id=point.ids[j];
            if(injection[id]>=0)row[injection[id]]+=sign[id]*(side==0?point.phi[j]:point.derivative[2][j]);
        }
        bool duplicate=false;
        for(std::size_t k=0;k<rows.size();++k)
        {
            double error=0,scale=1;for(int j=0;j<n;++j)
            {
                error=std::max(error,std::abs(row[j]-rows[k][j]));scale=std::max(scale,std::abs(row[j]));
            }
            if(error<1e-13*scale)
            {
                Require(std::abs(targets[k]-constraint[2])<1e-11*std::max(1.,std::abs(constraint[2])),
                        "inconsistent duplicate radial constraints");
                duplicate=true;break;
            }
        }
        if(!duplicate){rows.push_back(std::move(row));targets.push_back(constraint[2]);}
    }
    Vector coefficients(n);
    if(constant)
    {
        for(const auto &key:radial)if(std::get<0>(key.first)==std::get<1>(key.first))coefficients[key.second]=value;
        const auto actual=RANSWallAlgebra::Multiply(rows,coefficients);
        for(std::size_t j=0;j<actual.size();++j)
            if(std::abs(actual[j]-targets[j])>1e-12*std::max(1.,std::abs(targets[j])))constant=false;
    }
    if(!constant)
    {
        Require(rows.size()<=static_cast<std::size_t>(n),"too many radial constraints");
        WallMatrix system(n+rows.size(),Vector(n+rows.size()));Vector rhs(n+rows.size());
        double volume=0;
        for(const auto &item:intervals)volume+=item.second.area*(item.second.upper-item.second.lower);
        Require(volume>0 && std::isfinite(volume),"invalid radial projection volume");
        for(int i=0;i<n;++i)
        {
            rhs[i]=load[i]/volume;
            for(int j=0;j<n;++j)system[i][j]=mass[i][j]/volume;
        }
        for(std::size_t k=0;k<rows.size();++k)
        {
            long double square=0;for(double entry:rows[k])square+=(long double)entry*entry;
            const double norm=std::sqrt(square);Require(norm>0,"empty radial constraint row");
            rhs[n+k]=targets[k]/norm;
            for(int j=0;j<n;++j)system[j][n+k]=system[n+k][j]=rows[k][j]/norm;
        }
        const auto solved=RANSWallAlgebra::Solve(system,rhs);
        std::copy(solved.begin(),solved.begin()+n,coefficients.begin());
    }
    const auto actual=RANSWallAlgebra::Multiply(rows,coefficients);
    for(std::size_t j=0;j<actual.size();++j)
        Require(std::abs(actual[j]-targets[j])<1e-9*std::max(1.,std::abs(targets[j])),
                "radial KKT constraint residual");
    Vector result(Size());
    for(int j=0;j<Size();++j)if(injection[j]>=0)result[j]=sign[j]*coefficients[injection[j]];
    return result;
}
std::vector<RANSWallGlobalField::InletTraceElement>
RANSWallGlobalField::InletTraceElements() const
{
    Require(m_dimension==2,"FlatInletTraceV1 requires 2D Modified_A tensor fields");
    std::vector<InletTraceElement> result;const auto map=m_field->GetLocalToGlobalMap();
    for(int e=0;e<m_elements.size();++e)
    {
        InletTraceElement item;item.element=m_elements[e];
        const int start=m_field->GetCoeff_Offset(e),n=m_field->GetExp(e)->GetNcoeffs();
        for(int j=0;j<n;++j){item.ids.push_back(map->GetLocalToGlobalMap(start+j));item.signs.push_back(map->GetLocalToGlobalSign(start+j));}
        result.push_back(std::move(item));
    }
    return result;
}
RANSWallGlobalField::Vector RANSWallGlobalField::FlatInletTraceCoefficients(
    const Vector &values,const std::vector<InletTraceElement> &elements,double inlet)
{
    Require(!values.empty() && !elements.empty() && std::isfinite(inlet),"empty FlatInletTraceV1 data");
    using Interval=std::pair<double,double>;
    using Role=std::tuple<double,double,int,double,double,int>;
    struct Radial {int degree=-1;Vector trace;std::vector<Interval> columns;};
    std::map<Interval,Radial> radial;
    std::map<Role,std::pair<int,double>> roles;std::map<int,Role> ids;
    for(const auto &item:elements)
    {
        const auto &el=item.element;const int m=el.degree+1;
        Require(m>=2 && el.lower[0]>=inlet && el.upper[0]>el.lower[0] && el.upper[2]>el.lower[2] &&
            el.lower[1]==0 && el.upper[1]==0 && item.ids.size()==std::size_t(m*m) && item.signs.size()==item.ids.size(),
            "inconsistent FlatInletTraceV1 interval/degree/shape");
        auto &r=radial[{el.lower[2],el.upper[2]}];
        Require(r.degree<0 || r.degree==el.degree,"inconsistent FlatInletTraceV1 horizontal degrees");r.degree=el.degree;
        r.columns.emplace_back(el.lower[0],el.upper[0]);
        for(int j=0;j<m;++j)for(int i=0;i<m;++i)
        {
            const int n=i+m*j,id=item.ids[n];const double sign=item.signs[n];
            Require(id>=0 && id<values.size() && std::isfinite(values[id]) && (sign==1 || sign==-1),"invalid FlatInletTraceV1 native ID/sign");
            const double x=i<2?(i?el.upper[0]:el.lower[0]):el.lower[0],
                         z=j<2?(j?el.upper[2]:el.lower[2]):el.lower[2];
            const Role role{x,i<2?x:el.upper[0],i<2?0:i,z,j<2?z:el.upper[2],j<2?0:j};
            const auto old=roles.find(role);Require(old==roles.end() || old->second==std::make_pair(id,sign),"inconsistent FlatInletTraceV1 shared native trace ID/sign");roles[role]={id,sign};
            const auto alias=ids.find(id);Require(alias==ids.end() || alias->second==role,"aliased FlatInletTraceV1 native ID");ids[id]=role;
        }
        if(el.lower[0]==inlet)
        {
            Require(r.trace.empty(),"duplicate FlatInletTraceV1 inlet interval");
            for(int j=0;j<m;++j)r.trace.push_back(item.signs[m*j]*values[item.ids[m*j]]);
        }
    }
    bool first=true;double previous=0;std::vector<Interval> columns;
    for(auto &[interval,r]:radial)
    {
        Require((first || previous==interval.first) && !r.trace.empty(),"nonconforming FlatInletTraceV1 normal partition/inlet");previous=interval.second;
        std::sort(r.columns.begin(),r.columns.end());double end=inlet;
        for(const auto &c:r.columns){Require(c.first==end,"nonconforming FlatInletTraceV1 horizontal partition");end=c.second;}
        Require(first || columns==r.columns,"inconsistent FlatInletTraceV1 columns across normal intervals");columns=r.columns;first=false;
    }
    Vector result(values.size(),std::numeric_limits<double>::quiet_NaN());
    for(const auto &item:elements)
    {
        const auto &el=item.element;const int m=el.degree+1;const auto &trace=radial.at({el.lower[2],el.upper[2]}).trace;
        for(int j=0;j<m;++j)for(int i=0;i<m;++i)
        {
            const int n=i+m*j,id=item.ids[n];const double v=item.signs[n]*(i<2?trace[j]:0.);
            Require(!std::isfinite(result[id]) || result[id]==v,"inconsistent FlatInletTraceV1 injected shared coefficient");result[id]=v;
        }
    }
    for(double v:result)Require(std::isfinite(v),"unassigned FlatInletTraceV1 coefficient");
    return result;
}
RANSWallGlobalField::Vector RANSWallGlobalField::FlatInletTrace(const Vector &values,double inlet) const
{
    Require(values.size()==std::size_t(Size()),"FlatInletTraceV1 coefficient shape differs");
    return FlatInletTraceCoefficients(values,InletTraceElements(),inlet);
}

}
