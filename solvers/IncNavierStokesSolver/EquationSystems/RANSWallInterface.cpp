///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallInterface.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <algorithm>
#include <cmath>
#include <set>
namespace Nektar
{
using namespace RANSWallAlgebra;
using WallMatrix=RANSWallAlgebra::DenseMatrix;
namespace { bool Equal(double a,double b){return std::abs(a-b)<1e-11;} }
RANSWallInterface::RANSWallInterface(const RANSWallFieldLayout &layout,int order):m_layout(layout)
{
    const auto &l=m_layout;
    auto trace=l.innerVelocity->Boundary([](double,double){return 0.;},
        [&](int edge,const RANSWallNativeField::Element &e){return edge==2 && Equal(e.top,l.Height());});
    for(const auto &item:trace)m_traceIds.push_back(item.first);
    m_traceConstraint=WallMatrix(trace.size(),Vector(l.innerVelocity->Size()+l.outerVelocity->Size()));
    auto rule=LibUtilities::PointsManager()[LibUtilities::PointsKey(order,LibUtilities::eGaussGaussLegendre)];
    const auto x=rule->GetZ(),w=rule->GetW();
    for(const auto &e:l.innerVelocity->Elements())if(Equal(e.top,l.Height()))
        for(int i=0;i<x.size();++i)
        {
            const double point=e.left+(x[i]+1)*(e.right-e.left)/2;
            auto full=l.velocity->Point(point,l.Height(),-1);full.weight=w[i]*(e.right-e.left)/2;m_points.push_back(full);
            const auto a=l.innerVelocity->Point(point,l.Height(),-1),b=l.outerVelocity->Point(point,l.Height(),1);
            for(int row=0;row<m_traceIds.size();++row)
            {
                double eta=0;for(int j=0;j<a.ids.size();++j)if(a.ids[j]==m_traceIds[row])eta=a.phi[j];
                for(int j=0;j<a.ids.size();++j)m_traceConstraint[row][a.ids[j]]+=full.weight*eta*a.phi[j];
                for(int j=0;j<b.ids.size();++j)m_traceConstraint[row][l.innerVelocity->Size()+b.ids[j]]-=full.weight*eta*b.phi[j];
            }
        }
    for(int side=0;side<2;++side)
    {
        const auto sub=side?l.outerVelocity:l.innerVelocity;
        auto &map=m_injection[side];map.resize(sub->Size(),{-1,0});
        const auto sm=sub->Field()->GetLocalToGlobalMap(),fm=l.velocity->Field()->GetLocalToGlobalMap();
        for(int e=0;e<sub->Elements().size();++e)
        {
            int full=-1;for(int f=0;f<l.velocity->Elements().size();++f)
                if(sub->Elements()[e].geometry==l.velocity->Elements()[f].geometry)full=f;
            Require(full>=0,"native interface injection element missing");
            for(int j=0;j<sub->Field()->GetExp(e)->GetNcoeffs();++j)
            {
                const int si=sub->Field()->GetCoeff_Offset(e)+j,fi=l.velocity->Field()->GetCoeff_Offset(full)+j;
                const int id=sm->GetLocalToGlobalMap(si),global=fm->GetLocalToGlobalMap(fi);
                const double sign=fm->GetLocalToGlobalSign(fi)/sm->GetLocalToGlobalSign(si);
                if(map[id].first>=0)Require(map[id]==std::make_pair(global,sign),"inconsistent native interface injection");
                map[id]={global,sign};
            }
        }
    }
}
std::array<RANSWallInterface::Vector,2> RANSWallInterface::Loads(const Function &flux) const
{
    std::array<Vector,2> result{Vector(m_layout.innerVelocity->Size()),Vector(m_layout.outerVelocity->Size())};
    for(const auto &p:m_points)
    {
        const double common=flux(p.x,p.r);
        for(int side=0;side<2;++side)
        {
            const auto sub=side?m_layout.outerVelocity:m_layout.innerVelocity;
            const auto q=sub->Point(p.x,p.r,side?1:-1);
            for(int j=0;j<q.ids.size();++j)result[side][q.ids[j]]+=(side?-1:1)*p.weight*q.phi[j]*common;
        }
    }
    return result;
}
RANSWallInterface::Vector RANSWallInterface::Merge(const std::array<Vector,2> &loads) const
{
    Vector result(m_layout.velocity->Size());
    for(int side=0;side<2;++side)
    {
        Require(loads[side].size()==m_injection[side].size(),"native interface load size mismatch");
        for(int j=0;j<loads[side].size();++j)result[m_injection[side][j].first]+=m_injection[side][j].second*loads[side][j];
    }
    return result;
}
RANSWallInterface::PressureTest RANSWallInterface::ReferencePressure(const SolverUtils::RANSNeutralReference &model,double k) const
{
    const auto &l=m_layout;
    auto boundary=l.velocity->Boundary([](double,double){return 0.;},[&](int edge,const RANSWallNativeField::Element &e){
        return (edge==0 && Equal(e.bottom,0)) || (edge==1 && Equal(e.right,l.Right())) ||
            (edge==2 && Equal(e.top,l.Top())) || (edge==3 && Equal(e.left,l.Left()));});
    PressureTest result;
    for(int component=0;component<2;++component)for(int i=0;i<l.velocity->Size();++i)
        if(!boundary.count(i))result.freeVelocity.push_back(component*l.velocity->Size()+i);
    const int nv=result.freeVelocity.size(),np=l.innerPressure->Size()+l.outerPressure->Size();
    result.divergence=WallMatrix(np,Vector(nv));result.system=WallMatrix(nv+np+1,Vector(nv+np+1));
    result.rhs=Vector(nv+np+1);Vector pressureMean(np);const auto mass=l.velocity->Mass();
    for(int i=0;i<nv;++i)for(int j=0;j<nv;++j)
        if(result.freeVelocity[i]/l.velocity->Size()==result.freeVelocity[j]/l.velocity->Size())
            result.system[i][j]=mass[result.freeVelocity[i]%l.velocity->Size()][result.freeVelocity[j]%l.velocity->Size()];
    const double t0=k/model.Q(),ustar=std::sqrt(t0),H=l.Height();
    SolverUtils::RANSTensor gradient{};gradient[0][2]=ustar/(model.Kappa()*H);
    const auto closure=SolverUtils::EvaluateRANSClosure(model.Model(),k,ustar*ustar*ustar/(model.Kappa()*H),gradient,H,model.Parameters(),0);
    std::array<SolverUtils::RANSTensor,2> viscous{};
    viscous[0][0][2]=viscous[0][2][0]=t0;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)viscous[1][i][j]=-closure.reynoldsStress[i][j]+(i==j?2*k/3:0);
    for(int side=0;side<2;++side)
    {
        const auto velocity=side?l.outerVelocity:l.innerVelocity;
        const auto pressure=side?l.outerPressure:l.innerPressure;
        const int offset=side?l.innerPressure->Size():0;
        for(const auto &raw:velocity->Quadrature())
        {
            const auto v=l.velocity->Point(raw.x,raw.r),p=pressure->Point(raw.x,raw.r);
            for(int a=0;a<p.ids.size();++a)
            {
                pressureMean[offset+p.ids[a]]+=raw.weight*p.phi[a]/((l.Right()-l.Left())*l.Top());
                for(int i=0;i<nv;++i)
                {
                    const int component=result.freeVelocity[i]/l.velocity->Size(),id=result.freeVelocity[i]%l.velocity->Size();
                    for(int b=0;b<v.ids.size();++b)if(v.ids[b]==id)
                        result.divergence[offset+p.ids[a]][i]+=raw.weight*p.phi[a]*(component?v.dr[b]:v.dx[b]);
                }
            }
            for(int i=0;i<nv;++i)
            {
                const int component=result.freeVelocity[i]/l.velocity->Size(),id=result.freeVelocity[i]%l.velocity->Size(),physical=component?2:0;
                for(int b=0;b<v.ids.size();++b)if(v.ids[b]==id)
                    result.rhs[i]-=raw.weight*(v.dx[b]*viscous[side][physical][0]+v.dr[b]*viscous[side][physical][2]);
            }
        }
    }
    for(int p=0;p<np;++p)
    {
        for(int v=0;v<nv;++v){result.system[v][nv+p]=-result.divergence[p][v];result.system[nv+p][v]=result.divergence[p][v];}
        result.system[nv+p][nv+np]=result.system[nv+np][nv+p]=pressureMean[p];
    }
    result.solution=RANSWallAlgebra::Solve(result.system,result.rhs);
    result.velocityMaximum=Maximum(Vector(result.solution.begin(),result.solution.begin()+nv));result.globalMean=0;
    for(int p=0;p<np;++p)result.globalMean+=pressureMean[p]*result.solution[nv+p];
    Vector inner(result.solution.begin()+nv,result.solution.begin()+nv+l.innerPressure->Size()),
        outer(result.solution.begin()+nv+l.innerPressure->Size(),result.solution.begin()+nv+np);
    const double x=(l.Left()+l.Right())/2;
    result.pressureJump=l.innerPressure->Point(x,H,-1).Value(inner)-l.outerPressure->Point(x,H,1).Value(outer);
    result.tractionMismatch=std::abs(-result.pressureJump+viscous[0][2][2]-viscous[1][2][2]);
    return result;
}
}
