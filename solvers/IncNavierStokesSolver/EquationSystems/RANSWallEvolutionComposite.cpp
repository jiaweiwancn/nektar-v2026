// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionResolution.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <cmath>
#include <cstring>
#include <numeric>

namespace Nektar
{
using namespace SolverUtils;
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;
namespace
{
double StreamCurvature(const RANSWallNormalMap::Sample &g)
{
    double k=0;for(int i=0;i<3;++i)for(int j=0;j<3;++j)k+=g.tangent1[i]*g.B[i][j]*g.tangent1[j];
    return k;
}
}

RANSWallCoordinate RANSWallEvolutionPatch::CompositeChartCoordinate(const RANSWallCoordinate &raw) const
{
    auto p=raw;
    if(p[1]!=0.)Require(false,"composite 2D chart has no spanwise coordinate");
    // Native GetCoords can place a fitted endpoint a few ulps outside its
    // declared affine box. Normalize only that geometric roundoff; never
    // alter a field value or admit a genuinely out-of-domain query.
    for(int a:{0,2})
    {
        const double tolerance=16*std::numeric_limits<double>::epsilon()*std::max({1.,std::abs(lower[a]),std::abs(upper[a])});
        if(!std::isfinite(p[a]) || p[a]<lower[a]-tolerance || p[a]>upper[a]+tolerance)
            Require(false,"composite native chart query outside declared geometry: axis="+std::to_string(a)+" value="+J::Number(p[a]));
        if(p[a]<lower[a]){++m_chartEndpointSnaps;m_chartEndpointOvershoot=std::max(m_chartEndpointOvershoot,lower[a]-p[a]);p[a]=lower[a];}
        if(p[a]>upper[a]){++m_chartEndpointSnaps;m_chartEndpointOvershoot=std::max(m_chartEndpointOvershoot,p[a]-upper[a]);p[a]=upper[a];}
    }
    return p;
}
RANSWallNormalMap::Sample RANSWallEvolutionPatch::GeometryAt(const RANSWallCoordinate &raw,int radialSide,int footSide) const
{
    if(!Composite())return m_mapping.At(raw,footSide);
    const auto p=CompositeChartCoordinate(raw);
    if(!m_composite)Require(false,"composite chart queried before geometry construction");
    // For H=0 amplitude the normal chart is exactly the identity, including
    // its floating arithmetic. No eta round trip is needed in this limit.
    if(m_mapping.amplitude==0)return m_mapping.At(p,footSide);
    const double H=m_parameters.geometry.height;
    if(p[2]<H || (std::abs(p[2]-H)<1e-13 && radialSide<=0))return m_composite->Inner(p[0],std::clamp(p[2],0.,H),footSide);
    const double scale=upper[2]-H;const auto v=m_composite->Outer(p[0],std::clamp((p[2]-H)/scale,0.,1.),footSide);
    auto q=m_composite->Inner(p[0],H,footSide);q.point=v.point;q.jacobian=v.jacobian;q.inverse=v.inverse;
    for(int a=0;a<3;++a){q.jacobian[a][2]/=scale;q.inverse[2][a]*=scale;}
    q.determinant=v.determinant/scale;q.distance=p[2];return q;
}
RANSVector RANSWallEvolutionPatch::GeometryCofactor(const RANSWallCoordinate &p,int axis,int side,int footSide,int radialSide) const
{
    if(!Composite())return m_mapping.Cofactor(p,axis,side,footSide);
    const auto q=GeometryAt(p,radialSide,footSide);RANSVector a{};
    for(int j=0;j<3;++j)a[j]=side*q.determinant*q.inverse[axis][j];
    return a;
}
std::string RANSWallEvolutionPatch::GeometrySignature() const
{return Composite()?"CompositeHillV1:native-outer=[Hw,Z]:eta=(zHat-Hw)/(Z-Hw):"+m_composite->Signature():m_mapping.Signature();}
RANSWallCoordinate RANSWallEvolutionPatch::ChartFromPhysical(const RANSWallCoordinate &p,int side) const
{
    if(!Composite())return Mapped()?m_mapping.Inverse(p):p;
    const auto q=m_composite->Locate(p,side);auto c=q.chart;
    if(!q.inner)c[2]=m_parameters.geometry.height+c[2]*(upper[2]-m_parameters.geometry.height);
    return c;
}
RANSWallGlobalPoint RANSWallEvolutionPatch::PhysicalPoint(int block,const RANSWallCoordinate &p,int side,int footSide) const
{
    const auto c=ChartFromPhysical(p,side);
    for(int e=0;e<fields[block]->Elements().size();++e)
    {
        const auto &el=fields[block]->Elements()[e];bool yes=true;
        for(int a:{0,2})yes=yes && c[a]>=el.lower[a]-1e-11 && c[a]<=el.upper[a]+1e-11;
        if(side<0 && std::abs(c[2]-m_parameters.geometry.height)<1e-11)yes=yes && el.upper[2]<=m_parameters.geometry.height+1e-11;
        if(side>0 && std::abs(c[2]-m_parameters.geometry.height)<1e-11)yes=yes && el.lower[2]>=m_parameters.geometry.height-1e-11;
        if(footSide && (std::abs(c[0]-m_mapping.lengthX)<1e-11 || std::abs(c[0]+m_mapping.lengthX)<1e-11))
            yes=yes && (footSide<0?el.upper[0]<=c[0]+1e-11:el.lower[0]>=c[0]-1e-11);
        if(yes)return Transform(fields[block]->ElementPoint(e,c));
    }
    throw std::runtime_error("physical observation has no compatible native element/side");
}
double RANSWallEvolutionPatch::OuterWallDistance(const RANSWallCoordinate &p) const
{
    if(!Composite())return Distance(p);
    if(p[2]<=m_parameters.geometry.height)return p[2]; // certified normal strip, including Gamma_H
    // Fixed inventory certificates belong to this immutable geometry. Only
    // genuinely additional/moving queries use bounded transient storage.
    auto &cache=m_staticWallCoordinates.count(p)?m_wallDistances:m_dynamicWallDistances;
    const auto old=cache.find(p);if(old!=cache.end()){++m_wallDistanceHits;return old->second.distance;}
    const auto start=std::chrono::steady_clock::now();++m_wallDistanceBuilds;
    const auto result=m_composite->PhysicalWallDistance(Physical(p),1e-10);
    m_wallDistanceSeconds+=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    Require(result.distance>0 && result.errorBound<=1.01e-10,"uncertified outer physical wall distance");
    if(&cache==&m_dynamicWallDistances && cache.size()>=4096)cache.clear();
    cache.emplace(p,result);return result.distance;
}

double RANSWallEvolutionPatch::SeedRadial(int block,double z,int derivative) const
{
    double integral=0.;const auto &panels=m_seedRadial[block];
    for(const auto &c:panels)
    {
        if(z<c[0])break;
        const long double h=c[1]-c[0],t=std::min(1.L,(z-c[0])/h);
        if(derivative<0){integral+=h*t*(c[2]+t*(c[3]/2+t*(c[4]/3+t*c[5]/4)));if(z<=c[1])return integral;}
        else if(z<=c[1]+1e-14)return derivative==2?(2*c[4]+6*t*c[5])/(h*h):derivative?(c[3]+t*(2*c[4]+3*t*c[5]))/h:c[2]+t*(c[3]+t*(c[4]+t*c[5]));
    }
    if(derivative<0)return integral;
    throw std::runtime_error("seed radial evaluation outside represented panels");
}

RANSVector RANSWallEvolutionPatch::CompositeSeedVelocity(const RANSWallCoordinate &raw,int footSide) const
{
    const auto p=CompositeChartCoordinate(raw);
    if(m_compositeFlatSeed.empty())Require(false,"composite seed target needs its explicit radial setup projection");
    const double H=m_parameters.geometry.height,B=AnalyticalNeutralInflow()?.2:.4,z=Physical(p)[2];
    auto value=[&](double r){return SeedRadial(0,r);};
    auto derivative=[&](double r){return SeedRadial(0,r,1);};
    const auto g=m_composite->Inner(p[0],std::min(p[2],H),footSide);
    if(m_compositeSeedPolicy=="PhysicalHeightInterfaceV2")
    {
        // A setup choice, not a wall-law/source change: interface speed may
        // follow the incoming profile at its PHYSICAL height. The complete
        // streamfunction (including its normal velocity) is then projected.
        const auto gh=m_composite->Inner(p[0],H,footSide);
        const double h=gh.point[2],hs=gh.jacobian[2][0],dx=gh.jacobian[0][0],nz=gh.normal[2],UH=value(H);
        const double fp=gh.tangent1[2]/gh.tangent1[0],q=gh.baseArea,kappa=StreamCurvature(gh);
        const double fpp=-kappa*q*q*q;
        const double fppp=std::abs(p[0])<m_mapping.lengthX?m_mapping.amplitude*std::pow(std::acos(-1.)/m_mapping.lengthX,3)/2*std::sin(std::acos(-1.)*p[0]/m_mapping.lengthX):0.;
        const double nzs=fp*kappa,nzss=-(fpp*fpp+fp*fppp)/(q*q*q)+3*fp*fp*fpp*fpp/std::pow(q,5);
        const double hss=fpp+H*nzss;
        const double a=value(h)/(UH*nz);
        const double as=derivative(h)*hs/(UH*nz)-value(h)*nzs/(UH*nz*nz);
        const double ass=(SeedRadial(0,h,2)*hs*hs+derivative(h)*hss)/(UH*nz)-2*derivative(h)*hs*nzs/(UH*nz*nz)-value(h)*nzss/(UH*nz*nz)+2*value(h)*nzs*nzs/(UH*nz*nz*nz);
        const double FH=SeedRadial(0,H,-1);
        if(p[2]<=H)
        {
            const double ps=as*SeedRadial(0,p[2],-1),pr=a*value(p[2]);
            return {ps*g.inverse[0][2]+pr*g.inverse[2][2],0,-ps*g.inverse[0][0]-pr*g.inverse[2][0]};
        }
        if(z>=B)return {value(z),0,0};
        const double J=1+H*kappa,kappas=-fppp/std::pow(q,3)+3*fp*fpp*fpp/std::pow(q,5),qs=fp*fpp/q;
        const double I=gh.inverse[0][2],Is=fpp/(q*q*J)-fp*(2*q*qs*J+q*q*H*kappas)/(q*q*q*q*J*J);
        const double delta=B-h,t=(z-h)/delta,ts=-hs*(1-t)/delta;
        const double C=a*FH-SeedRadial(0,h,-1),Cs=as*FH-value(h)*hs,D=as*FH*I,Ds=ass*FH*I+as*FH*Is;
        const double h00=1-3*t*t+2*t*t*t,h10=t-2*t*t+t*t*t,d00=-6*t+6*t*t,d10=1-4*t+3*t*t;
        const double u=value(z)+C*d00/delta+D*d10;
        const double ps=Cs*h00+C*d00*ts+(-hs*D+delta*Ds)*h10+delta*D*d10*ts;
        return {u,0,-ps/dx};
    }
    if(p[2]<=H)return {value(p[2])*g.tangent1[0],0,value(p[2])*g.tangent1[2]};
    if(z>=B)return {value(z),0,0};
    auto integral=[&](double top){return SeedRadial(0,top,-1);};
    const double h=g.point[2],hs=g.jacobian[2][0],dx=g.jacobian[0][0],delta=B-h,t=(z-h)/delta;
    const double U=value(h),UH=value(H),C=integral(H)-integral(h),D=UH*g.normal[2]-U;
    // n_z'= -slope*slope'/(1+slope^2)^(3/2), with signed curvature
    // kappa=-f''/(1+f'^2)^(3/2). Thus n_z'=slope*kappa.
    const double slope=g.tangent1[2]/g.tangent1[0],nz_s=slope*StreamCurvature(g);
    const double Cs=-U*hs,Ds=UH*nz_s-derivative(h)*hs,ts=-hs*(1-t)/delta;
    const double h00=1-3*t*t+2*t*t*t,h10=t-2*t*t+t*t*t,d00=-6*t+6*t*t,d10=1-4*t+3*t*t;
    const double u=value(z)+C*d00/delta+D*d10;
    const double psi_s=Cs*h00+C*d00*ts+(-hs*D+delta*Ds)*h10+delta*D*d10*ts;
    return {u,0,-psi_s/dx};
}

double RANSWallEvolutionPatch::CompositeSeedStreamfunction(const RANSWallCoordinate &p,int footSide) const
{
    const double H=m_parameters.geometry.height,B=AnalyticalNeutralInflow()?.2:.4;
    const auto gh=m_composite->Inner(p[0],H,footSide);const double h=gh.point[2],z=Physical(p)[2],nz=gh.normal[2],UH=SeedRadial(0,H);
    const bool adjusted=m_compositeSeedPolicy=="PhysicalHeightInterfaceV2";
    const double a=adjusted?SeedRadial(0,h)/(UH*nz):1.;
    if(p[2]<=H)return a*SeedRadial(0,p[2],-1);
    if(z>=B)return SeedRadial(0,z,-1);
    const double slope=gh.tangent1[2]/gh.tangent1[0],nzs=slope*StreamCurvature(gh);
    const double as=adjusted?SeedRadial(0,h,1)*gh.jacobian[2][0]/(UH*nz)-SeedRadial(0,h)*nzs/(UH*nz*nz):0.;
    const double C=a*SeedRadial(0,H,-1)-SeedRadial(0,h,-1);
    const double D=adjusted?as*SeedRadial(0,H,-1)*gh.inverse[0][2]:UH*nz-SeedRadial(0,h);
    const double width=B-h,t=(z-h)/width;
    return SeedRadial(0,z,-1)+C*(1-3*t*t+2*t*t*t)+width*D*(t-2*t*t+t*t*t);
}

std::string RANSWallEvolutionPatch::CompositeSeedDerivativeChecks() const
{
    Require(Composite() && AnalyticalNeutralInflow(),"signed seed checks require analytical composite setup");
    std::vector<std::string> rows;double maximum=0;
    const double H=m_parameters.geometry.height,L=m_mapping.lengthX;
    for(double s:{0.,.059553093889933906,-L,L})for(int side:{-1,1})
    {
        if(std::abs(s)!=L && side==-1)continue;
        const int foot=std::abs(s)==L?side:0;
        for(double r:{.6*H,H,.08})
        {
            const RANSWallCoordinate p{s,0,r};const auto metric=GeometryAt(p,r<=H?-1:1,foot);
            const auto velocity=CompositeSeedVelocity(p,foot);
            auto differentiated=[&](double step)
            {
                Vector d(3);
                for(int axis:{0,2})
                {
                    const double h=axis==0?step:step*.1;auto a=p,b=p;
                    int direction=0;if(axis==0 && foot)direction=foot;if(axis==2 && r==H)direction=-1;
                    if(direction){a[axis]+=direction*h;b[axis]+=direction*2*h;d[axis]=direction*(-3*CompositeSeedStreamfunction(p,foot)+4*CompositeSeedStreamfunction(a)-CompositeSeedStreamfunction(b))/(2*h);}
                    else{a[axis]+=h;b[axis]-=h;d[axis]=(CompositeSeedStreamfunction(a)-CompositeSeedStreamfunction(b))/(2*h);}
                }
                RANSVector v{};for(int j:{0,2}){v[0]+=d[j]*metric.inverse[j][2];v[2]-=d[j]*metric.inverse[j][0];}return v;
            };
            const auto first=differentiated(2e-6),second=differentiated(1e-6);double error=0,halving=0;
            for(int j:{0,2}){error=std::max(error,std::abs(velocity[j]-second[j]));halving=std::max(halving,std::abs(first[j]-second[j]));}
            auto targetDivergence=[&](double scale)
            {
                RANSTensor chart{};const auto v0=CompositeSeedVelocity(p,foot);
                for(int axis:{0,2})
                {
                    const double h=scale*(axis==0?1e-7:1e-8);auto a=p,b=p;int direction=axis==0?foot:(r==H?-1:0);
                    if(direction){a[axis]+=direction*h;b[axis]+=direction*2*h;const auto va=CompositeSeedVelocity(a),vb=CompositeSeedVelocity(b);for(int i:{0,2})chart[i][axis]=direction*(-3*v0[i]+4*va[i]-vb[i])/(2*h);}
                    else{a[axis]+=h;b[axis]-=h;const auto va=CompositeSeedVelocity(a),vb=CompositeSeedVelocity(b);for(int i:{0,2})chart[i][axis]=(va[i]-vb[i])/(2*h);}
                }
                double div=0;for(int i:{0,2})for(int j:{0,2})div+=chart[i][j]*metric.inverse[j][i];return div;
            };
            maximum=std::max(maximum,error);const auto ground=m_composite->Inner(s,0,foot);
            rows.push_back(J::Object({{"chart",J::Numbers(p)},{"physical",J::Numbers(metric.point)},{"tangential_side",J::Number(foot)},
                {"signed_stream_curvature",J::Number(StreamCurvature(ground))},{"sorted_principal",J::Numbers(Vector{ground.kappa1,ground.kappa2})},
                {"velocity",J::Numbers(velocity)},{"streamfunction_FD_velocity",J::Numbers(second)},
                {"velocity_error",J::Number(error)},{"step_halving",J::Number(halving)},{"target_divergence",J::Number(targetDivergence(.5))},{"divergence_step_halving",J::Number(std::abs(targetDivergence(1)-targetDivergence(.5)))}}));
        }
    }
    // Physical incoming reference differs from zHat at this nonflat point.
    const RANSWallCoordinate c{.04,0,.1};const auto ref=InflowAt(Physical(c));
    const double expected=m_robustUstar/(m_parameters.model.Kappa()*Physical(c)[2]);
    Require(std::abs(ref.dU/expected-1)<1e-14 && std::abs(Physical(c)[2]-c[2])>1e-3,"physical neutral reference diagnostic height");
    Require(maximum<2e-5,"signed seed velocity disagrees with independently differentiated represented streamfunction");
    return J::Object({{"passed","true"},{"maximum_velocity_error",J::Number(maximum)},{"samples",J::List(rows)},
        {"physical_reference_point",J::Numbers(Physical(c))},{"native_zHat",J::Number(c[2])},{"exact_reference_gradient_departure",J::Number(ref.dU/expected-1)},
        {"exact_reference_epsilon_ratio",J::Number(ref.epsilon/(std::pow(m_robustUstar,3)/(m_parameters.model.Kappa()*Physical(c)[2])))}});
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::CompositeSeed() const
{
    Require(Composite() && PrescribedInflow() && AcceptedAWF() && OpenBoundary(),"mapped seed requires the composite finite/open accepted-AWF contract");
    RANSWallTiming::Scope timing("startup.composite_seed");m_compositeFlatSeed=FlatFiniteSeed();
    const double H=m_parameters.geometry.height,B=AnalyticalNeutralInflow()?.2:.4;auto s=m_compositeFlatSeed;
    if(m_mapping.amplitude==0){m_compositeSeedInfo=J::Object({{"method",J::Quote("exact composite-flat radial seed")},{"governing_solves","0"}});return s;}
    for(int block:{0,K(),Epsilon()})
    {
        const int which=block==0?0:block==K()?1:2;m_seedRadial[which].clear();std::set<std::pair<double,double>> bounds;
        for(const auto &el:fields[block]->Elements())bounds.emplace(el.lower[2],el.upper[2]);
        for(auto [a,b]:bounds)
        {
            std::array<long double,4> v{};for(int j=0;j<4;++j)v[j]=Point(block,{lower[0],0,a+j*(b-a)/3}).Value(s,offset[block]);
            m_seedRadial[which].push_back({a,b,v[0],-5.5L*v[0]+9*v[1]-4.5L*v[2]+v[3],9*v[0]-22.5L*v[1]+18*v[2]-4.5L*v[3],-4.5L*v[0]+13.5L*v[1]-13.5L*v[2]+4.5L*v[3]});
        }
    }
    for(int b=0;b<2;++b)
    {
        auto values=PhysicalProject(b,[&](const RANSWallCoordinate &p,double){return CompositeSeedVelocity(p)[PhysicalAxis(b)];},0.);
        std::copy(values.begin(),values.end(),s.begin()+offset[b]);
    }
    // Scalar seed is the unchanged native inlet completion extended in the
    // chart, with its exact endpoint K and dependent epsilon. It is NOT used
    // as incoming reference height; gamma/S_ref always use physical height.
    // The scalar tail is constant above native .4 (physical .5 is a safe
    // whole-domain tail certificate). No interpolation clipping is performed.
    std::map<int,double> fixed;
    for(int b=0;b<InnerPressure();++b)
    {
        auto boundary=fields[b]->Boundary([](const auto&){return 0.;},[&](int axis,int side,const auto &el)
        {return (axis==0 && !side && el.lower[0]==lower[0]) || (axis==2 && !side && b<K() && el.lower[2]==0.) ||
            (AnalyticalNeutralInflow() && axis==2 && side && b==1 && el.upper[2]==upper[2]) ||
            (axis==2 && !side && b==Epsilon() && el.lower[2]==H);});
        for(auto [id,v]:boundary)fixed[offset[b]+id]=m_compositeFlatSeed[offset[b]+id];
    }
    // Keep the intended inner budget-reference k exactly. This is seed setup,
    // not a prescribed ground-k condition in physical evolution.
    for(const auto &q:Quadrature(K(),4))if(q.coordinate[2]<H)for(int id:q.ids)fixed[offset[K()]+id]=m_compositeFlatSeed[offset[K()]+id];
    // Only the historical finite completion has an exactly constant tail
    // whose coefficients may be frozen wholesale. A represented logarithm
    // above a curved chart is not exactly divergence free; leave its velocity
    // coefficients available to the physical-mass continuity projection.
    if(!AnalyticalNeutralInflow())for(int b=0;b<2;++b)for(int e=0;e<fields[b]->Elements().size();++e)
        if(fields[b]->Elements()[e].lower[2]>=B)for(int id:fields[b]->ElementPoint(e,fields[b]->Elements()[e].lower).ids)
            fixed[offset[b]+id]=AnalyticalNeutralInflow()?(b==0?s[offset[b]+id]:0.):(b==0?InflowAtHeight(B).U:0.)*m_constants[b][id];
    // Consistent physical L2 projection with the exact setup constraints.
    for(int b=0;b<InnerPressure();++b)
    {
        std::vector<int> rows;for(int j=0;j<size[b];++j)if(!fixed.count(offset[b]+j))rows.push_back(j);
        Vector rhs(rows.size());WallMatrix matrix(rows.size(),Vector(rows.size()));
        for(int i=0;i<rows.size();++i)
        {for(auto [id,v]:fixed)if(id>=offset[b] && id<offset[b]+size[b])rhs[i]-=m_mass[b][rows[i]][id-offset[b]]*(v-s[id]);
         for(int j=0;j<rows.size();++j)matrix[i][j]=m_mass[b][rows[i]][rows[j]];}
        if(!rows.empty()){const auto d=RANSWallAlgebra::Solve(matrix,rhs);for(int i=0;i<rows.size();++i)s[offset[b]+rows[i]]+=d[i];}
    }
    for(auto [id,v]:fixed)s[id]=v;
    std::vector<int> free;for(int b=0;b<2;++b)for(int j=0;j<size[b];++j)if(!fixed.count(offset[b]+j))free.push_back(offset[b]+j);
    std::vector<int> index(Size(),-1);for(int i=0;i<free.size();++i)index[free[i]]=i;
    WallMatrix D(size[InnerPressure()]+size[OuterPressure()],Vector(free.size()));Vector defect(D.size());
    std::vector<const VolumePoint*> continuity;for(const auto &q:m_bodyInner)continuity.push_back(&q);for(const auto &q:m_outer)continuity.push_back(&q);
    for(const auto *entry:continuity)
    {
        const auto &q=entry->v;const double weight=entry->weight;
        const int b=q.coordinate[2]<H?InnerPressure():OuterPressure();const auto p=entry->p;
        const int start=b==InnerPressure()?0:size[InnerPressure()];const double div=q.Gradient(s,offset[0])[0]+q.Gradient(s,offset[1])[2];
        for(int i=0;i<p.ids.size();++i)
        {const int row=start+p.ids[i];defect[row]+=weight*p.phi[i]*div;
         for(int a=0;a<2;++a)for(int j=0;j<q.ids.size();++j)if(index[offset[a]+q.ids[j]]>=0)
             D[row][index[offset[a]+q.ids[j]]]+=weight*p.phi[i]*q.derivative[PhysicalAxis(a)][j];}
    }
    WallMatrix Q;Vector beta;
    for(int i=0;i<D.size();++i)
    {
        auto v=D[i];double b=-defect[i],original=std::sqrt(std::inner_product(v.begin(),v.end(),v.begin(),0.));
        for(int pass=0;pass<2;++pass)for(int j=0;j<Q.size();++j)
        {double a=std::inner_product(v.begin(),v.end(),Q[j].begin(),0.);for(int k=0;k<v.size();++k)v[k]-=a*Q[j][k];b-=a*beta[j];}
        const double norm=std::sqrt(std::inner_product(v.begin(),v.end(),v.begin(),0.));
        if(norm>1e-11*std::max(1.,original)){for(double &x:v)x/=norm;Q.push_back(v);beta.push_back(b/norm);}
        else Require(std::abs(b)<1e-10,"mapped seed continuity constraints incompatible with fixed data");
    }
    WallMatrix mass(free.size(),Vector(free.size()));
    for(int i=0;i<free.size();++i)for(int j=0;j<free.size();++j)
        if(free[i]/size[0]==free[j]/size[0])mass[i][j]=m_mass[0][free[i]%size[0]][free[j]%size[0]];
    m_seedM=std::make_shared<Factor>(mass);const auto &M=*m_seedM;WallMatrix inverse,schur(Q.size(),Vector(Q.size()));
    for(const auto &row:Q)inverse.push_back(M.Apply(row));
    for(int i=0;i<Q.size();++i)for(int j=0;j<Q.size();++j)schur[i][j]=std::inner_product(Q[i].begin(),Q[i].end(),inverse[j].begin(),0.);
    m_seedSchur=std::make_shared<Factor>(schur);
    m_seedFree=free;m_seedQ=Q;m_seedInverse=inverse;m_seedMass=mass;
    const auto multiplier=m_seedSchur->Apply(beta);Vector delta(free.size());
    for(int j=0;j<Q.size();++j)for(int i=0;i<free.size();++i)delta[i]+=inverse[j][i]*multiplier[j];
    for(int i=0;i<free.size();++i)s[free[i]]+=delta[i];
    auto remaining=Multiply(D,delta);
    for(int i=0;i<D.size();++i)remaining[i]+=defect[i];
    Require(Maximum(remaining)<1e-10,"mapped seed native continuity projection failed");
    m_compositeSeedInfo=J::Object({{"initializer_policy",J::Quote(m_compositeSeedPolicy)},{"method",J::Quote(AnalyticalNeutralInflow()?"physical streamfunction lift blending to exact incoming neutral law above z=.2; no constant tail; physical-mass/native continuity setup":"physical streamfunction lift; unstrained tail above z=.4; physical consistent-mass constrained native projection")},
        {"setup_radial_projections","3"},{"setup_field_projections","2"},{"setup_continuity_mass_factors","1"},
        {"setup_physical_projection_mass_solves","6"},{"setup_fixed_constraint_mass_solves","4"},
        {"setup_continuity_schur_solves","1"},{"independent_continuity_rows",J::Number(Q.size())},
        {"native_weak_continuity_max",J::Number(Maximum(remaining))},{"velocity_projection_mass_norm",J::Number(std::sqrt(std::max(0.,std::inner_product(delta.begin(),delta.end(),Multiply(mass,delta).begin(),0.))))},
        {"governing_solves","0"},{"physical_steps","0"}});return s;
}

std::string RANSWallEvolutionPatch::CompositeFeasibilityData(const Vector &s,const RANSWallEvolutionData &data) const
{
    Require(Composite() && AnalyticalNeutralInflow() && m_seedM && m_seedSchur,"feasibility data needs the corrected native seed projection");
    const auto started=std::chrono::steady_clock::now();const double Z=upper[2],xmin=lower[0],span=upper[0]-xmin;
    std::vector<Vector> directions;std::vector<std::string> definitions;double constraintError=0;
    // 16 and 32 columns are the only declared search spaces. Smooth physical
    // streamfunctions vanish with their normal derivative at the floor and
    // inlet; roof psi=0 permits tangential velocity there. No z=.2 restriction.
    for(double cx:{-.075,-.025,.025,.075,-.13,-.045,.045,.18})for(double height:{.025,.07,.18,.5})
    {
        const double width=cx==.18?.2:.065;
        Vector d(Size());
        auto target=[&](const RANSWallCoordinate &c,int component)
        {
            const auto p=Physical(c);const double x=p[0],z=p[2];
            const double pi=std::acos(-1.),L=m_mapping.lengthX;
            const double f=std::abs(x)<L?m_mapping.amplitude*std::pow(std::cos(pi*x/(2*L)),2):0.;
            const double fp=std::abs(x)<L?-m_mapping.amplitude*pi/(2*L)*std::sin(pi*x/L):0.;
            const double h=z-f,t=(x-xmin)/span,gauss=std::exp(-std::pow((x-cx)/width,2));
            const double A=t*t*gauss,Ax=(2*t/span-2*t*t*(x-cx)/(width*width))*gauss;
            const double E=std::exp(-h/height),B=h*h*(Z-z)*E;
            const double Bz=(2*h*(Z-z)-h*h-h*h*(Z-z)/height)*E;
            const double Bx=-fp*(2*h-h*h/height)*(Z-z)*E;
            return component==0?A*Bz:-(Ax*B+A*Bx);
        };
        for(int b=0;b<2;++b){auto v=PhysicalProject(b,[&](const auto &p,double){return target(p,b);},0);std::copy(v.begin(),v.end(),d.begin()+offset[b]);}
        Vector v(m_seedFree.size());for(int i=0;i<v.size();++i)v[i]=d[m_seedFree[i]];
        Vector defect(m_seedQ.size());for(int j=0;j<defect.size();++j)defect[j]=std::inner_product(m_seedQ[j].begin(),m_seedQ[j].end(),v.begin(),0.);
        const auto lambda=m_seedSchur->Apply(defect);
        for(int j=0;j<lambda.size();++j)for(int i=0;i<v.size();++i)v[i]-=m_seedInverse[j][i]*lambda[j];
        const auto mv=Multiply(m_seedMass,v);const double norm=std::sqrt(std::inner_product(v.begin(),v.end(),mv.begin(),0.));
        Require(norm>1e-14,"degenerate projected seed basis");
        std::fill(d.begin(),d.end(),0.);for(int i=0;i<v.size();++i)d[m_seedFree[i]]=v[i]/norm;
        for(const auto &q:m_seedQ)constraintError=std::max(constraintError,std::abs(std::inner_product(q.begin(),q.end(),v.begin(),0.))/norm);
        directions.push_back(d);definitions.push_back(J::Numbers(Vector{cx,height,width,norm}));
    }
    Require(constraintError<1e-9,"homogeneous projected basis violates native continuity");
    const auto domain=ScalarDomain(s,0,data,true);Require(domain.admissible,"feasibility scalar context inadmissible");
    std::vector<std::string> observations;std::set<std::tuple<int,int,int,bool,RANSWallCoordinate>> seen;
    auto tensor=[](const RANSTensor &g){std::vector<std::string> v;for(const auto &r:g)v.push_back(J::Numbers(r));return J::List(v);};
    for(const auto &p:domain.observations)
    {
        if(p.block!=K() && !p.dependentTrace)continue;
        const auto origin=ResolveObservation(p);const auto &row=*origin.velocity;
        if(!seen.emplace(row.geometry,row.element,origin.radialSide,p.dependentTrace,p.coordinate).second)continue;
        RANSTensor G{};for(int a=0;a<2;++a)G[PhysicalAxis(a)]=row.Gradient(s,offset[a]);
        const double k=EvaluationRow(K(),GeometryElement(K(),row.geometry),p.coordinate)->Value(s,offset[K()]);
        const double r=Distance(p.coordinate);double e=0,nt=0;
        if(origin.inner)nt=StateProfile(s,p.coordinate,0,data)->Original().EddyViscosity(r);
        else e=p.dependentTrace?p.value:EvaluationRow(Epsilon(),GeometryElement(Epsilon(),row.geometry),p.coordinate)->Value(s,offset[Epsilon()]);
        std::vector<std::string> gradients,velocities;
        for(const auto &d:directions){RANSTensor g{};for(int a=0;a<2;++a)g[PhysicalAxis(a)]=row.Gradient(d,offset[a]);gradients.push_back(tensor(g));velocities.push_back(J::Numbers(Vector{row.Value(d,offset[0]),row.Value(d,offset[1])}));}
        observations.push_back(J::Object({{"chart",J::Numbers(p.coordinate)},{"physical",J::Numbers(Physical(p.coordinate))},{"inner",J::Boolean(origin.inner)},
            {"geometry",J::Number(row.geometry)},{"element",J::Number(row.element)},{"side",J::Number(origin.radialSide)},{"dependent_epsilon",J::Boolean(p.dependentTrace)},
            {"family",J::Quote(p.inventory)},{"gradient",tensor(G)},{"gradient_basis",J::List(gradients)},{"velocity_basis",J::List(velocities)},{"k",J::Number(k)},{"epsilon",J::Number(e)},{"inner_nu_t",J::Number(nt)}}));
    }
    std::vector<std::string> basis,gram;
    for(const auto &v:directions)basis.push_back(J::Numbers(v));
    for(const auto &v:directions)
    {Vector free(m_seedFree.size());for(int i=0;i<free.size();++i)free[i]=v[m_seedFree[i]];const auto mv=Multiply(m_seedMass,free);Vector row;
     for(const auto &w:directions){double dot=0;for(int i=0;i<free.size();++i)dot+=w[m_seedFree[i]]*mv[i];row.push_back(dot);}gram.push_back(J::Numbers(row));}
    return J::Object({{"version",J::Quote("VelocityNullspacePhaseIV1")},{"base",J::Numbers(s)},{"basis",J::List(basis)},
        {"basis_parameters_cx_height_width_original_massnorm",J::List(definitions)},{"mass_normalized",J::Boolean(true)},{"physical_mass_Gram",J::List(gram)},
        {"homogeneous_continuity_max",J::Number(constraintError)},{"free_velocity_rows",J::Numbers(m_seedFree)},{"velocity_size",J::Number(size[0])},
        {"observations",J::List(observations)},{"seconds",J::Number(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count())},
        {"setup_only","true"},{"governing_solves","0"}});
}

std::string RANSWallEvolutionPatch::FlatSpatialDiagnostic(const Vector &s,const RANSWallEvolutionData &data) const
{
    Require(!Mapped() && AnalyticalNeutralInflow() && !OpenBoundary(),"zero-step spatial comparison requires compatible affine neutral boundaries");
    const auto resolved=ResidualResolutionComponents(s,s,0,0,data,false,false,false,true);
    const auto analytical=ResidualResolutionComponents(s,s,0,0,data,false,false,true,true);
    const double H=m_parameters.geometry.height,K0=InflowAtHeight(H).k;
    const auto profile=Profile(K0,{lower[0],0,0},1);const auto &total=resolved.terms.at("spatial_sum");
    WallMatrix duals;std::vector<std::vector<int>> frees;Vector energies;
    for(int b=0;b<InnerPressure();++b)
    {
        std::vector<int> free;for(int i=0;i<size[b];++i)if(!resolved.constraints.count(offset[b]+i))free.push_back(i);
        WallMatrix M(free.size(),Vector(free.size()));Vector v(free.size());
        for(int i=0;i<free.size();++i){v[i]=total[offset[b]+free[i]];for(int j=0;j<free.size();++j)M[i][j]=m_mass[b][free[i]][free[j]];}
        auto d=Factor(M).Apply(v);energies.push_back(std::inner_product(v.begin(),v.end(),d.begin(),0.));duals.push_back(d);frees.push_back(free);
    }
    auto norm=[&](Vector v){return Norms(v,resolved.constraints,m_robustUstar);};
    std::vector<std::pair<std::string,std::string>> terms,analyticTerms,panels;Vector panelSum(Size());
    for(const auto &[name,v]:resolved.terms)
    {
        if(name.rfind("radial:",0)==0)
        {
            Vector projection;
            for(int b=0;b<InnerPressure();++b){double product=0;for(int j=0;j<frees[b].size();++j)product+=v[offset[b]+frees[b][j]]*duals[b][j];projection.push_back(energies[b]>0?product/energies[b]:0.);}
            panels.push_back({name,J::Object({{"scaled_dual_norm",J::Numbers(norm(v))},{"signed_dual_projection_on_total",J::Numbers(projection)}})});
            for(int i=0;i<Size();++i)panelSum[i]+=v[i];
        }
        else if(std::find(RANSWallResolutionSpatialTerms().begin(),RANSWallResolutionSpatialTerms().end(),name)!=RANSWallResolutionSpatialTerms().end())
            terms.push_back({name,J::Numbers(norm(v))});
    }
    for(const auto &name:RANSWallResolutionSpatialTerms())analyticTerms.push_back({name,J::Numbers(norm(analytical.terms.at(name)))});
    // Compatible verification boundaries have no open loads; only the gauge
    // lies outside volume panels. Test differential free rows, with all global
    // sharing and cross terms retained before contraction.
    double sumError=0;for(int b=0;b<InnerPressure();++b)for(int id:frees[b])sumError=std::max(sumError,std::abs(panelSum[offset[b]+id]-total[offset[b]+id]));
    Require(sumError<1e-10,"radial free-load decomposition failed");
    std::map<std::pair<double,double>,Vector> error;
    for(const auto &q:Quadrature(0,32))
    {
        const auto &el=fields[0]->Elements()[q.element];auto &v=error[{el.lower[2],el.upper[2]}];if(v.empty())v=Vector(8);
        const double z=q.coordinate[2];const bool inner=z<H;const auto ref=InflowAtHeight(std::max(z,H));
        const double U=inner?profile->Original().Velocity(z,{InflowAtHeight(H).U,0},{0,0})[0]:ref.U;
        const double dU=inner?profile->Original().VelocityGradient(z,{InflowAtHeight(H).U,0},{0,0})[0]:ref.dU;
        RANSTensor g{},gr{};for(int a=0;a<2;++a)g[PhysicalAxis(a)]=q.Gradient(s,offset[a]);gr[0][2]=dU;
        const double k=q.Value(s,offset[K()]),kr=inner?profile->ReferenceK(z):K0;
        const double e=inner?profile->Dissipation(z):Point(Epsilon(),q.coordinate).Value(s,offset[Epsilon()]);
        RANSTensor actualStress{},exactStress{};double actualProduction=0,exactProduction=0;
        if(inner){const auto a=InnerStress(k,profile->Original().EddyViscosity(z),g),b=InnerStress(kr,profile->Original().EddyViscosity(z),gr);actualStress=a.reynoldsStress;exactStress=b.reynoldsStress;actualProduction=a.production;exactProduction=b.production;}
        else{const auto a=OuterClosure(k,e,g,q.coordinate),b=OuterClosure(kr,ref.epsilon,gr,q.coordinate);actualStress=a.reynoldsStress;exactStress=b.reynoldsStress;actualProduction=a.production;exactProduction=b.production;}
        const Vector d{q.Value(s,offset[0])-U,g[0][2]-dU,k-kr,q.Gradient(s,offset[K()])[2]-(inner?profile->ReferenceKGradient(z):0.),inner?0:e-ref.epsilon,
            actualStress[0][2]-exactStress[0][2]-(inner?m_parameters.geometry.molecularNu*(g[0][2]+g[2][0]-dU):0.),actualStress[2][2]-exactStress[2][2]-(inner?2*m_parameters.geometry.molecularNu*g[2][2]:0.),actualProduction-exactProduction};
        for(int j=0;j<8;++j)v[j]+=q.weight*d[j]*d[j];
    }
    std::vector<std::string> errs,jumps;
    for(auto [bounds,v]:error){for(auto &x:v)x=std::sqrt(x);errs.push_back(J::Object({{"bounds",J::Numbers(Vector{bounds.first,bounds.second})},{"dimensional_errors",J::Numbers(v)}}));}
    for(auto [bounds,v]:error)if(bounds.second<upper[2])
    {
        const double z=bounds.second;const RANSWallCoordinate c{(lower[0]+upper[0])/2,0,z};
        std::vector<double> flux;
        for(int b:{0,K(),Epsilon()})
        {
            if(b==Epsilon() && z<=H){flux.push_back(0);continue;}
            const auto minus=Point(b,c,-1),plus=Point(b,c,1);flux.push_back(plus.Gradient(s,offset[b])[2]-minus.Gradient(s,offset[b])[2]);
        }
        std::vector<std::string> sided;
        for(int side:{-1,1})
        {
            const auto row=Point(0,c,side);const bool inner=fields[0]->Elements()[row.element].upper[2]<=H;
            RANSTensor g{};for(int a=0;a<2;++a)g[PhysicalAxis(a)]=row.Gradient(s,offset[a]);
            const auto kr=Point(K(),c,side);const double k=kr.Value(s,offset[K()]);const auto dk=kr.Gradient(s,offset[K()]);
            const double pi=Point(inner?InnerPressure():OuterPressure(),c,side).Value(s,offset[inner?InnerPressure():OuterPressure()]);
            RANSTensor stress{};double diffusion=0,epsilonFlux=0;
            if(inner){const double nt=profile->Original().EddyViscosity(z);stress=RANSWallInnerMeanStress(InnerStress(k,nt,g),m_parameters.geometry.molecularNu,pi);diffusion=m_parameters.geometry.molecularNu+nt/m_parameters.model.Parameters().sigmaK;}
            else{const auto er=Point(Epsilon(),c,side);const auto law=OuterClosure(k,er.Value(s,offset[Epsilon()]),g,c);diffusion=law.kDiffusion;
                epsilonFlux=-law.epsilonDiffusion*er.Gradient(s,offset[Epsilon()])[2];
                for(int a=0;a<3;++a)for(int b=0;b<3;++b)stress[a][b]=-law.reynoldsStress[a][b]+(a==b?2*k/3-pi:0);}
            Vector traction;for(int a=0;a<3;++a)traction.push_back(stress[a][2]);
            sided.push_back(J::Object({{"radial_side",J::Number(side)},{"geometry",J::Number(row.geometry)},{"inner",J::Boolean(inner)},
                {"selected_mean_traction_on_plus_z",J::Numbers(traction)},{"physical_k_diffusive_flux_plus_z",J::Number(-diffusion*dk[2])},
                {"physical_epsilon_diffusive_flux_plus_z",inner?"null":J::Number(epsilonFlux)}}));
        }
        jumps.push_back(J::Object({{"z",J::Number(z)},{"raw_derivative_jumps_U_k_epsilon",J::Numbers(flux)},{"one_sided_physical_fluxes",J::List(sided)}}));
    }
    auto analytic=analytical.terms.at("spatial_sum");
    return J::Object({{"operation",J::Quote("H2dFlatSpatialV1")},{"physical_steps","0"},{"governing_solves","0"},{"setup_projection_only","true"},
        {"size",J::Number(Size())},{"state",J::Numbers(s)},{"selected_spatial_norms",J::Numbers(norm(total))},{"analytical_free_load_norms",J::Numbers(norm(analytic))},
        {"physical_reference",FiniteInflowChecks(s)},{"radial_errors",J::List(errs)},{"error_columns",J::Quote("U,gradU,k,gradk,outer_epsilon,shear_stress,normal_stress,production")},
        {"radial_components",J::Object(panels)},{"global_components",J::Object(terms)},{"analytical_component_norms",J::Object(analyticTerms)},{"one_sided_derivative_jumps",J::List(jumps)},
        {"radial_global_free_sum_error",J::Number(sumError)},{"scope",J::Quote("zero-step representation diagnostic; shared global rows combined before full-mass contraction; signed projections retain cancellation; p-enriched diagnostic does not qualify production positivity")}});
}

std::string RANSWallEvolutionPatch::MovingFootChecks(const Vector &s,const RANSWallEvolutionData &data) const
{
    const double H=m_parameters.geometry.height,x=-m_mapping.lengthX;
    const auto wall=Profile(SampleK(s,{x,0,0}),{x,0,0},1,+1);
    const auto nearest=std::min_element(wall->Quadrature().begin(),wall->Quadrature().end(),
        [&](const auto &a,const auto &b){return std::abs(a.r-.6*H)<std::abs(b.r-.6*H);});
    const double r=nearest->r;
    const RANSWallCoordinate p{x,0,r};auto test=s;
    std::fill(test.begin()+offset[0],test.begin()+offset[0]+size[0],0.);
    const auto v=fields[1]->Project([&](const RANSWallCoordinate &c)
        {return 1e9*std::max(0.,c[0]-x)*c[2]*(H-c[2]);});
    std::copy(v.begin(),v.end(),test.begin()+offset[1]);
    std::vector<std::string> sides;double eigen[2];int i=0;
    for(int foot:{-1,1})
    {
        const auto row=RayObservationRow(0,p,foot);RANSTensor g{};
        for(int a=0;a<2;++a)g[PhysicalAxis(a)]=row.Gradient(test,offset[a]);
        const auto profile=Profile(SampleK(test,p),p,1,foot);
        const double k=row.Value(test,offset[K()]);
        const auto stress=InnerStress(k,profile->Original().EddyViscosity(r),g);
        eigen[i++]=J::Eigenvalues(stress.reynoldsStress)[0];
        sides.push_back(J::Object({{"tangential_side",J::Number(foot)},{"geometry",J::Number(row.geometry)},
            {"element",J::Number(row.element)},{"minimum_covariance",J::Number(eigen[i-1])}}));
    }
    Require(eigen[0]>0 && eigen[1]<0,"moving foot fixture did not separate the one-sided gradients");
    const auto domain=ScalarDomain(test,0,data,true);Require(domain.admissible,"foot fixture changed scalar admissibility");
    int minus=0,plus=0;bool offGrid=true;
    for(const auto &o:domain.observations)
        if(o.coordinate==p && o.inventory=="constitutive_foot_minus_budget_ray")minus=o.geometry+1;
        else if(o.coordinate==p && o.inventory=="constitutive_foot_plus_budget_ray")plus=o.geometry+1;
    for(const auto &o:ReconstructionLocations(K(),true))if(o.row.coordinate==p)offGrid=false;
    Require(minus && plus && minus!=plus && offGrid,"moving foot sides absent, aliased or on fixed grid");
    const auto before=RANSWallR3::Hash(test);const auto audit=CompleteRestrictedAudit(test,test,0,0,data);
    Require(before==RANSWallR3::Hash(test) && audit.at("accepted_awf_hard_admissible")[0]==0,"negative moving covariance escaped hard audit");
    return J::Object({{"passed","true"},{"point",J::Numbers(p)},{"sides",J::List(sides)},
        {"off_fixed_reconstruction_grid",J::Boolean(offGrid)},{"minus_plus_inventory_geometries",J::Numbers(Vector{double(minus-1),double(plus-1)})},
        {"audit",J::Object({{"accepted_awf_hard_admissible",J::Numbers(audit.at("accepted_awf_hard_admissible"))}})},
        {"physical_solves","0"}});
}

std::string RANSWallEvolutionPatch::CompositeChecks(const Vector &s,double time,const RANSWallEvolutionData &data) const
{
    Require(AcceptedAWF() && PrescribedInflow() && OpenBoundary(),"H2b operators require the declared finite/open budget AWF");
    const auto before=RANSWallR3::Hash(s);const double H=m_parameters.geometry.height;
    if(Composite())
    {
        auto invalid=upper;invalid[0]+=1e-6;bool rejected=false;
        try{(void)GeometryAt(invalid);}catch(const std::exception&){rejected=true;}
        Require(rejected,"composite geometry accepted a truly outside native query");
        auto endpoint=lower;endpoint[0]=std::nextafter(lower[0],-std::numeric_limits<double>::infinity());
        Require(GeometryAt(endpoint).point==GeometryAt(lower).point,"one-ulp native endpoint normalization differs");
    }
    const auto manufactured=[&](const RANSWallCoordinate &p){return 1.+.2*p[0]+.3*p[2]+.1*p[0]*p[2];};
    const auto v=fields[0]->Project(manufactured);
    double inverseError=0,gradientError=0,interfacePoint=0,interfaceArea=0,minJ=1e300,massError=0;
    int sample=0;Vector convergence;std::vector<std::string> integralChecks;
    for(int order:{8,16,32})
    {
        double volume=0,div=0,flux=0;Vector load(size[0]);
        for(const auto &q:Quadrature(0,order))
        {
            volume+=q.weight;const auto g=q.Gradient(v);div+=q.weight*(g[0]+.3*g[2]);
            for(int i=0;i<q.ids.size();++i)load[q.ids[i]]+=q.weight*q.phi[i]*q.Value(v);
            const auto &el=fields[0]->Elements()[q.element];int foot=0;
            for(double f:{-m_mapping.lengthX,m_mapping.lengthX})if(std::abs(q.coordinate[0]-f)<1e-12)foot=(el.lower[0]+el.upper[0])/2<f?-1:1;
            const auto map=GeometryAt(q.coordinate,el.upper[2]<=H?-1:1,foot);minJ=std::min(minJ,map.determinant);
            if(order==8 && ++sample%37==0)
            {
                const auto recovered=ChartFromPhysical(map.point,el.upper[2]<=H?-1:1);
                for(int a:{0,2})inverseError=std::max(inverseError,std::abs(recovered[a]-q.coordinate[a]));
                const double hs=1e-7,hr=1e-8;auto a=q.coordinate,b=q.coordinate,c=q.coordinate,d=q.coordinate;
                a[0]+=hs;b[0]-=hs;c[2]+=hr;d[2]-=hr;
                const auto pa=Physical(a),pb=Physical(b),pc=Physical(c),pd=Physical(d);
                const double A=(pa[0]-pb[0])/(2*hs),B=(pc[0]-pd[0])/(2*hr),C=(pa[2]-pb[2])/(2*hs),D=(pc[2]-pd[2])/(2*hr),det=A*D-B*C;
                const double fs=.2+.1*q.coordinate[2],fr=.3+.1*q.coordinate[0];
                gradientError=std::max({gradientError,std::abs(g[0]-(fs*D-fr*C)/det),std::abs(g[2]-(-fs*B+fr*A)/det)});
            }
        }
        for(int axis:{0,2})for(int side:{-1,1})for(const auto &f:SurfaceView(0,axis,side<0?lower[axis]:upper[axis],side,order))
            flux+=f.weight*f.test.Value(v)*(f.area[0]+.3*f.area[2]);
        const auto action=Multiply(m_mass[0],v);double error=0;for(int i=0;i<action.size();++i)error=std::max(error,std::abs(action[i]-load[i]));
        massError=error;convergence.push_back(volume);
        integralChecks.push_back(J::Object({{"order",J::Number(order)},{"volume",J::Number(volume)},
            {"divergence_integral",J::Number(div)},{"boundary_flux",J::Number(flux)},{"divergence_flux_error",J::Number(std::abs(div-flux))},
            {"physical_mass_action_error",J::Number(error)}}));
        if(order==32)Require(std::abs(div-flux)<2e-9,"curved native divergence/flux identity failed");
    }
    const double expected=(upper[0]-lower[0])*upper[2]-(Composite()?m_mapping.amplitude*m_mapping.lengthX:0.);
    Require(std::abs(convergence.back()-expected)<1e-9 && minJ>0 && gradientError<2e-6 && inverseError<2e-10 && massError<1e-9,
        "composite native map/gradient/volume/mass qualification failed");
    for(const auto &f:SurfaceView(0,2,H,1,16))
    {
        const auto p=f.test.coordinate;const auto lo=GeometryAt(p,-1),hi=GeometryAt(p,1);
        const auto a=GeometryCofactor(p,2,1,0,-1),b=GeometryCofactor(p,2,-1,0,1);
        for(int j=0;j<3;++j){interfacePoint=std::max(interfacePoint,std::abs(lo.point[j]-hi.point[j]));interfaceArea=std::max(interfaceArea,std::abs(a[j]+b[j]));}
    }
    Require(interfacePoint<1e-12 && interfaceArea<1e-12,"composite interface does not coincide with opposite area vectors");
    auto external=ExternalConstraints(time,data);
    // The governing epsilon block owns its dependent trace, including the
    // inlet/interface corner, exactly as in PimpleSingleSweep.
    for(int id:m_traceIds)external.erase(offset[Epsilon()]+id);
    const auto flow=AssemblePicardBlock(0,s,s,37.,time,data,external);
    Vector direction(Size());
    for(int b:{InnerPressure(),OuterPressure()})
    {const auto coeff=fields[b]->Project([&](const RANSWallCoordinate &p){return .3+.07*p[0]+.2*p[2];});for(int j=0;j<size[b];++j)direction[offset[b]+j]=coeff[j];}
    const auto action=AffinePressureAction(direction,external);double matrixG=0;
    for(int i=0;i<flow.rows.size();++i)if(flow.rows[i]<offset[K()])
    {double value=0;for(int j=0;j<flow.rows.size();++j)value+=flow.matrix[i][j]*direction[flow.rows[j]];matrixG=std::max(matrixG,std::abs(value-action[flow.rows[i]]));}
    auto plus=s,minus=s;const double step=1e-4;for(int j=0;j<Size();++j){plus[j]+=step*direction[j];minus[j]-=step*direction[j];}
    const auto rp=Residual(plus,s,37.,time,data,false),rm=Residual(minus,s,37.,time,data,false);double nativeG=0;
    for(int i=0;i<offset[K()];++i)nativeG=std::max(nativeG,std::abs((rp[i]-rm[i])/(2*step)-action[i]));
    Require(matrixG<1e-10 && nativeG<1e-8,"composite complete native pressure/roughness action failed");
    const auto parts=ResidualResolutionComponents(s,s,0,time,data);const auto native=Residual(s,s,0,time,data,false,false);double componentError=0;
    for(int i=0;i<Size();++i)componentError=std::max(componentError,std::abs(parts.terms.at("total_unmasked")[i]-native[i]));
    Require(componentError<1e-9,"mapped native residual decomposition differs");
    auto scalar=AssemblePicardBlock(K(),s,s,37.,time,data,external);auto epsilon=AssemblePicardBlock(Epsilon(),s,s,37.,time,data,external);
    const auto open=DevelopingFlux(s);double outletTraction=0;
    for(int j=0;j<size[0];++j)outletTraction+=m_constants[0][j]*open.traction[j];
    Require(std::abs(outletTraction-m_outletTotalPressure*upper[2])<1e-12,"nonzero outlet traction has incorrect physical measure");
    std::vector<std::string> boundary;
    for(double sign:{-1.,1.})
    {
        auto q=s;for(int j=0;j<size[0];++j){q[j]=sign*m_constants[0][j];q[offset[1]+j]=0;}
        const auto flux=DevelopingFlux(q);Vector errors;
        for(int b=0;b<InnerPressure();++b)
        {
            double sum=0,independent=0;for(int j=0;j<size[b];++j)sum+=m_constants[b][j]*flux.advection[offset[b]+j];
            for(const auto &f:SurfaceView(b,0,upper[0],1,32))independent+=f.weight*sign*f.test.Value(sign>0?q:m_backflowSeed,offset[b]);
            errors.push_back(std::abs(sum-independent));
        }
        Require(Maximum(errors)<1e-10,"mapped outflow/backflow native flux mismatch");
        boundary.push_back(J::Object({{"outward_sign",J::Number(sign)},{"errors",J::Numbers(errors)}}));
    }
    // Frozen full selected stress divergence theorem for every native test
    // coefficient. Constant tensor is a manufactured flux, not an RANS root.
    const auto &oq=m_outer.front();RANSTensor g{};
    for(int a=0;a<2;++a)g[PhysicalAxis(a)]=oq.v.Gradient(s,offset[a]);
    const double ok=oq.v.Value(s,offset[K()]);
    const auto law=OuterClosure(ok,oq.e.Value(s,offset[Epsilon()]),g,oq.v.coordinate);RANSTensor sigma{};
    for(int a=0;a<3;++a)for(int b=0;b<3;++b)sigma[a][b]=-law.reynoldsStress[a][b]+(a==b?2*ok/3-.013:0.);
    Vector tensorErrors;
    for(int b:{0,1,K(),Epsilon()})
    {
        Vector r(size[b]);const RANSVector flux=b<K()?sigma[PhysicalAxis(b)]:RANSVector{.07,0,.13};
        for(const auto &q:Quadrature(b,32))for(int j=0;j<q.ids.size();++j)r[q.ids[j]]+=q.weight*(q.derivative[0][j]*flux[0]+q.derivative[2][j]*flux[2]);
        for(int axis:{0,2})for(int side:{-1,1})
        {
            const double lower=b==Epsilon()?H:0.;const double position=axis==0?(side<0?this->lower[0]:upper[0]):(side<0?lower:upper[2]);
            for(const auto &f:SurfaceView(b,axis,position,side,32))for(int j=0;j<f.test.ids.size();++j)r[f.test.ids[j]]-=f.weight*f.test.phi[j]*(flux[0]*f.area[0]+flux[2]*f.area[2]);
        }
        tensorErrors.push_back(Maximum(r));
    }
    Require(Maximum(tensorErrors)<1e-9,"curved full-stress/scalar-flux native identity failed");
    // Projection and interface-epsilon are distinct diagnostic fixtures.
    auto raw=s;const auto ground=Point(K(),{upper[0],0,0});int id=-1;
    for(int j=0;j<ground.ids.size();++j)if(std::abs(ground.phi[j])>.5 && !external.count(offset[K()]+ground.ids[j])){id=j;break;}
    Require(id>=0,"mapped outlet-ground k unexpectedly prescribed");raw[offset[K()]+ground.ids[id]]-=(ground.Value(s,offset[K()])+1e-4)/ground.phi[id];
    const auto repair=ApplyPimpleKBound(raw,s,time,data,true);
    Require(repair.success && repair.activated && ScalarDomain(repair.state,time,data).admissible,"mapped physical-mass projection fixture failed");
    Vector delta(size[K()]);for(int j=0;j<size[K()];++j)delta[j]=repair.state[offset[K()]+j]-raw[offset[K()]+j];
    const auto load=Multiply(m_mass[K()],delta);
    const double objective=.5*std::inner_product(delta.begin(),delta.end(),load.begin(),0.);
    long double independentObjective=0.;for(const auto &q:Quadrature(K(),32))
        independentObjective+=(long double).5*q.weight*std::pow(q.Value(delta),2);
    Require(std::abs(objective-double(independentObjective))<1e-10*std::max(1e-30,objective),
        "mapped QP objective differs from independent physical integration");
    for(int j=0;j<Size();++j)if(j<offset[K()] || j>=offset[Epsilon()] || external.count(j))
        Require(std::memcmp(&s[j],&repair.state[j],sizeof(double))==0,"mapped correction changed essential/non-k data");
    auto bad=s;const auto ep=Point(Epsilon(),{0,0,H},1);id=-1;
    for(int j=0;j<ep.ids.size();++j)if(std::abs(ep.phi[j])>.5 && !external.count(offset[Epsilon()]+ep.ids[j])){id=j;break;}
    Require(id>=0,"crest interface epsilon fixture needs a free matching coefficient");bad[offset[Epsilon()]+ep.ids[id]]-=(ep.Value(s,offset[Epsilon()])+1e-5)/ep.phi[id];
    const auto rejected=ScalarDomain(bad,time,data);
    Require(!rejected.admissible && rejected.failure.find("epsilon")!=std::string::npos,"actual mapped H-plus epsilon was concealed by dependent trace");
    Refresh(s,time,&data);Require(before==RANSWallR3::Hash(s),"mapped operator diagnostics changed selected fields");
    Vector massAction=Multiply(m_mass[0],v),masked;for(int row:flow.rows)masked.push_back(action[row]);
    return J::Object({{"status",J::Quote("passed")},{"geometry",J::Quote(GeometrySignature())},{"inverse_max",J::Number(inverseError)},
        {"native_endpoint_roundoff_snaps",J::Number(m_chartEndpointSnaps)},{"native_endpoint_max_overshoot",J::Number(m_chartEndpointOvershoot)},
        {"gradient_fd_max",J::Number(gradientError)},{"minimum_J",J::Number(minJ)},{"integral_checks",J::List(integralChecks)},
        {"expected_volume",J::Number(expected)},{"interface_point_area_errors",J::Numbers(Vector{interfacePoint,interfaceArea})},
        {"full_native_G_matrix_fd_errors",J::Numbers(Vector{matrixG,nativeG})},{"native_component_error",J::Number(componentError)},
        {"frozen_selected_stress_and_scalar_flux_errors",J::Numbers(tensorErrors)},{"outflow_backflow",J::List(boundary)},
        {"outlet_datum_traction_integral",J::Number(outletTraction)},{"Pout",J::Number(m_outletTotalPressure)},
        {"projection_only",repair.report},{"interface_actual_epsilon_rejection",ScalarDomainJSON(rejected)},
        {"projection_physical_mass_objective",J::Number(objective)},
        {"projection_independent_GL32_objective",J::Number(independentObjective)},
        {"physical_mass_action",J::Numbers(massAction)},{"native_pressure_action",J::Numbers(action)},
        {"native_spatial_residual",J::Numbers(native)},{"scalar_frozen_load",J::Numbers(scalar.residual)},
        {"epsilon_frozen_load",J::Numbers(epsilon.residual)},{"scalar_free_rows",J::Numbers(Vector(scalar.rows.begin(),scalar.rows.end()))},
        {"epsilon_free_rows",J::Numbers(Vector(epsilon.rows.begin(),epsilon.rows.end()))},{"governing_solves","0"},{"physical_commits","0"},
        {"diagnostic_block_assemblies","3"},{"projection_fixture_calls","1"},{"selected_state_sha256",J::Quote(before)}});
}
RANSWallEvolutionPatch::OuterInitializerRows RANSWallEvolutionPatch::InitializerRows(const Vector &base) const
{
    Require(base.size()==Size(),"initializer base size");OuterInitializerRows out;
    const auto essential=FixedReferenceBoundary(base);std::set<int> fixedK,fixedE;
    const auto map=fields[K()]->Field()->GetLocalToGlobalMap();
    for(int e=0;e<fields[K()]->Elements().size();++e)if(fields[K()]->Elements()[e].upper[2]<=m_parameters.geometry.height+1e-13)
    {const int start=fields[K()]->Field()->GetCoeff_Offset(e),n=fields[K()]->Field()->GetExp(e)->GetNcoeffs();for(int j=0;j<n;++j)fixedK.insert(map->GetLocalToGlobalMap(start+j));}
    for(int id:m_traceIds)fixedE.insert(id);
    for(int b=0;b<InnerPressure();++b)for(int i=0;i<size[b];++i)
    {
        const int row=offset[b]+i;if(essential.count(row))continue;
        if(b<K())out.velocity.push_back(row);
        else if(b==K() && !fixedK.count(i))out.k.push_back(row);
        else if(b==Epsilon() && !fixedE.count(i))out.epsilon.push_back(row);
    }
    return out;
}

std::string RANSWallEvolutionPatch::ValidateOuterInitializer(const Vector &s,const Vector &base,const RANSWallEvolutionData &data,bool coupled) const
{
    Require(s.size()==Size() && base.size()==Size(),"coupled initializer size mismatch");
    const auto rows=InitializerRows(base);std::set<int> permitted(rows.velocity.begin(),rows.velocity.end());
    if(coupled){permitted.insert(rows.k.begin(),rows.k.end());permitted.insert(rows.epsilon.begin(),rows.epsilon.end());}
    int changed=0;for(int i=0;i<Size();++i)
    {Require(std::isfinite(s[i]),"nonfinite initializer coefficient");if(!permitted.count(i))Require(std::memcmp(&s[i],&base[i],sizeof(double))==0,"coupled initializer changed frozen coefficient bits at "+std::to_string(i));else if(std::memcmp(&s[i],&base[i],sizeof(double)))++changed;}
    double boundary=0,endpoint=0,trace=0;
    for(int b=0;b<InnerPressure();++b)for(const auto &q:SurfaceView(b,0,lower[0],-1,12))boundary=std::max(boundary,std::abs(q.test.Value(s,offset[b])-q.test.Value(base,offset[b])));
    for(const auto &q:SurfaceView(K(),2,m_parameters.geometry.height,-1,12))endpoint=std::max(endpoint,std::abs(q.test.Value(s,offset[K()])-q.test.Value(base,offset[K()])));
    for(const auto &q:SurfaceView(Epsilon(),2,m_parameters.geometry.height,-1,12))trace=std::max(trace,std::abs(q.test.Value(s,offset[Epsilon()])-q.test.Value(base,offset[Epsilon()])));
    Require(boundary==0 && endpoint==0 && trace==0,"coupled initializer changed inlet, endpoint or transported trace functional");
    const auto d0=ScalarDomain(base,0,data,true),d1=ScalarDomain(s,0,data,true);
    Require(d0.admissible && d1.admissible,"coupled initializer scalar domain rejected");
    Require(d0.observations.size()==d1.observations.size(),"fixed K initializer changed inventory cardinality");
    double lowK=1e300,highK=0,lowE=1e300,highE=0;int locations=0;
    for(int i=0;i<d0.observations.size();++i)
    {
        const auto &a=d0.observations[i],&b=d1.observations[i];
        Require(a.coordinate==b.coordinate && a.geometry==b.geometry && a.side==b.side && a.tangentialSide==b.tangentialSide && a.inventory==b.inventory && a.dependentTrace==b.dependentTrace,"fixed K initializer moved inventory/profile geometry");
        if(a.dependentTrace)Require(std::memcmp(&a.value,&b.value,sizeof(double))==0,"coupled initializer changed dependent epsilon trace");
        if(!a.dependentTrace && a.coordinate[2]>m_parameters.geometry.height && (a.block==K() || a.block==Epsilon()))
        {const double ratio=b.value/a.value;if(a.block==K()){lowK=std::min(lowK,ratio);highK=std::max(highK,ratio);}else{lowE=std::min(lowE,ratio);highE=std::max(highE,ratio);}++locations;}
    }
    Require(lowK>=.5-1e-12 && highK<=2+1e-12 && lowE>=.5-1e-12 && highE<=2+1e-12,"coupled initializer scalar envelope rejected");
    return J::Object({{"version",J::Quote("OuterScalarInitializerV1")},{"base_sha256",J::Quote(RANSWallR3::Hash(base))},{"state_sha256",J::Quote(RANSWallR3::Hash(s))},{"velocity_rows",J::Numbers(rows.velocity)},{"outer_k_rows",J::Numbers(rows.k)},{"outer_epsilon_rows",J::Numbers(rows.epsilon)},
        {"changed_permitted_rows",J::Number(changed)},{"frozen_bits_preserved","true"},{"boundary_from_immutable_base","true"},{"inlet_functional_difference",J::Number(boundary)},{"endpoint_difference",J::Number(endpoint)},{"trace_difference",J::Number(trace)},
        {"K_envelope",J::Numbers(Vector{lowK,highK})},{"epsilon_envelope",J::Numbers(Vector{lowE,highE})},{"fixed_inventory_locations",J::Number(locations)},{"moving_inventory_unchanged","true"}});
}

std::string RANSWallEvolutionPatch::CoupledInitializerContext(const Vector &s,const RANSWallEvolutionData &data,const WallMatrix &directions,const std::string &path) const
{
    Require(Composite() && AnalyticalNeutralInflow() && m_seedM && directions.size()==32,"H2e needs corrected seed and 32 saved native directions");
    const auto started=std::chrono::steady_clock::now();const auto allowed=InitializerRows(s);WallMatrix all=directions;double continuity=0,frozen=0;const int nv=32;
    for(const auto &d:directions)
    {
        Require(d.size()==Size(),"saved direction dimension");Vector f(m_seedFree.size());for(int i=0;i<f.size();++i)f[i]=d[m_seedFree[i]];
        for(const auto &q:m_seedQ)continuity=std::max(continuity,std::abs(std::inner_product(q.begin(),q.end(),f.begin(),0.)));
        const std::set<int> free(m_seedFree.begin(),m_seedFree.end());for(int i=0;i<Size();++i)if(!free.count(i))frozen=std::max(frozen,std::abs(d[i]));
    }
    Require(continuity<1e-9 && frozen==0,"saved direction violates frozen rows or continuity");
    const double H=m_parameters.geometry.height,span=upper[0]-lower[0];
    for(int b:{K(),Epsilon()})for(double cx:{-.075,-.025,.025,.18})for(double height:{.08,.35})
    {
        const double width=cx==.18?.2:.065;Vector d(Size());
        auto projected=PhysicalProject(b,[&](const RANSWallCoordinate &p,double)
        {if(p[2]<=H)return 0.;const double t=(p[0]-lower[0])/span,r=(p[2]-H)/height;return Point(b,p).Value(s,offset[b])*t*t*std::exp(-std::pow((p[0]-cx)/width,2))*r*r*std::exp(-r);},0.);
        const auto &rows=b==K()?allowed.k:allowed.epsilon;for(int i:rows)d[i]=projected[i-offset[b]];all.push_back(std::move(d));
    }
    auto binary=[&](const std::string &name,const WallMatrix &a)
    {std::ofstream out(path+"/"+name+".f64",std::ios::binary);Require(out.good(),"cannot create H2e binary context");for(const auto &r:a)out.write(reinterpret_cast<const char*>(r.data()),r.size()*sizeof(double));out.close();Require(out.good(),"H2e binary write failed");return J::Object({{"file",J::Quote(name+".f64")},{"shape",J::Numbers(std::vector<int>{int(a.size()),int(a.empty()?0:a[0].size())})},{"sha256",J::Quote(RANSWallR3::HashFile(path+"/"+name+".f64"))},{"dtype",J::Quote("<f8")}});};
    WallMatrix gram(nv,Vector(nv));for(int i=0;i<nv;++i){Vector f(m_seedFree.size());for(int j=0;j<f.size();++j)f[j]=all[i][m_seedFree[j]];auto mf=Multiply(m_seedMass,f);for(int k=0;k<nv;++k)for(int j=0;j<f.size();++j)gram[i][k]+=all[k][m_seedFree[j]]*mf[j];}
    WallMatrix div(nv+1,Vector(nv+1)),strain=div,roof=div,sg(16,Vector(16)),smooth=sg;double volume=0,outerVolume=0,maxDiv=0;
    for(const auto &q:Quadrature(0,16))
    {
        const bool inner=fields[0]->Elements()[q.element].upper[2]<=H;Vector d(nv+1);std::vector<RANSTensor> S(nv+1);
        for(int j=0;j<=nv;++j){const auto &v=j?all[j-1]:s;RANSTensor g{};for(int a=0;a<2;++a)g[PhysicalAxis(a)]=q.Gradient(v,offset[a]);d[j]=g[0][0]+g[2][2];for(int a=0;a<3;++a)for(int b=0;b<3;++b)S[j][a][b]=g[a][b]+g[b][a];}
        volume+=q.weight;maxDiv=std::max(maxDiv,std::abs(d[0]));
        for(int i=0;i<=nv;++i)for(int j=0;j<=nv;++j)div[i][j]+=q.weight*d[i]*d[j];
        if(!inner)
        {
            outerVolume+=q.weight;const auto ref=RegularizedReferenceAt(q.coordinate);for(int a=0;a<3;++a)for(int b=0;b<3;++b)S[0][a][b]-=ref.referenceStrain[a][b];
            for(int i=0;i<=nv;++i)for(int j=0;j<=nv;++j)for(int a=0;a<3;++a)for(int b=0;b<3;++b)strain[i][j]+=q.weight*ref.gamma*S[i][a][b]*S[j][a][b];
            for(int block:{K(),Epsilon()})
            {const auto row=GeometryPoint(block,q.geometry,q.coordinate);const double base=row.Value(s,offset[block]);const auto bg=row.Gradient(s,offset[block]);Vector v(8);std::vector<RANSVector> g(8);const int off=block==K()?0:8;
             for(int i=0;i<8;++i){v[i]=row.Value(all[nv+off+i],offset[block])/base;g[i]=row.Gradient(all[nv+off+i],offset[block]);for(int a=0;a<3;++a)g[i][a]=g[i][a]/base-v[i]*bg[a]/base;}
             for(int i=0;i<8;++i)for(int j=0;j<8;++j){sg[off+i][off+j]+=q.weight*v[i]*v[j];for(int a=0;a<3;++a)smooth[off+i][off+j]+=q.weight*.02*.02*g[i][a]*g[j][a];}}
        }
    }
    for(auto &r:sg)for(auto &v:r){v/=outerVolume;}
    for(auto &r:smooth)for(auto &v:r){v/=outerVolume;}
    for(const auto &f:SurfaceView(0,2,upper[2],1,16)){Vector d(nv+1);for(int j=0;j<=nv;++j)d[j]=f.test.Gradient(j?all[j-1]:s,offset[0])[0];double area=0;for(double v:f.area)area+=v*v;for(int i=0;i<=nv;++i)for(int j=0;j<=nv;++j)roof[i][j]+=f.weight*std::sqrt(area)*d[i]*d[j];}
    std::vector<std::pair<std::string,std::string>> files{{"base",binary("base",WallMatrix{s})},{"basis",binary("basis",all)},{"mass_gram",binary("mass_gram",gram)},{"divergence",binary("divergence",div)},{"strain",binary("strain",strain)},{"roof",binary("roof",roof)},{"scalar_gram",binary("scalar_gram",sg)},{"scalar_smoothness",binary("scalar_smoothness",smooth)}};
    // Every original membership is retained; only identical native evaluations
    // share numerical storage. The key retains both sided identities and the
    // actual/dependent epsilon context, including signed binary64 coordinates.
    const auto domain=ScalarDomain(s,0,data,true);Require(domain.admissible,"H2e base scalar gate");
    std::map<std::string,int> seen,families;std::vector<std::string> familyNames;std::ofstream points(path+"/observations.f64",std::ios::binary),membership(path+"/membership.i32",std::ios::binary);int n=0,total=0;const int width=18+nv*6+16;
    for(const auto &p:domain.observations)
    {
        if(p.block!=K() && p.block!=Epsilon() && !p.dependentTrace)continue;
        const auto origin=ResolveObservation(p);const auto &row=*origin.velocity;
        std::ostringstream key;key<<row.geometry<<':'<<row.element<<':'<<origin.radialSide<<':'<<origin.tangentialSide<<':'<<p.dependentTrace<<':'<<p.block<<':'<<RANSWallR3::Hash(p.coordinate.data(),3*sizeof(double))<<':'<<RANSWallR3::Hash(&p.value,p.dependentTrace?sizeof(double):0);
        auto it=seen.find(key.str());int id;
        if(it!=seen.end())id=it->second;
        else
        {
            id=n++;seen.emplace(key.str(),id);Vector r;for(double v:p.coordinate)r.push_back(v);for(double v:Physical(p.coordinate))r.push_back(v);
            RANSTensor g{};for(int a=0;a<2;++a)g[PhysicalAxis(a)]=row.Gradient(s,offset[a]);for(int a:{0,2})for(int b:{0,2})r.push_back(g[a][b]);for(int a=0;a<2;++a)r.push_back(row.Value(s,offset[a]));
            const auto kr=GeometryPoint(K(),row.geometry,p.coordinate);double e=0,nt=0;if(origin.inner)nt=StateProfile(s,p.coordinate,0,data)->Original().EddyViscosity(Distance(p.coordinate));else e=p.dependentTrace?p.value:GeometryPoint(Epsilon(),row.geometry,p.coordinate).Value(s,offset[Epsilon()]);
            r.insert(r.end(),{kr.Value(s,offset[K()]),e,nt,double(origin.inner),double(p.block==K() || p.dependentTrace),double(p.dependentTrace)});
            for(int j=0;j<nv;++j){for(int a=0;a<2;++a){const auto grad=row.Gradient(all[j],offset[a]);r.push_back(grad[0]);r.push_back(grad[2]);}for(int a=0;a<2;++a)r.push_back(row.Value(all[j],offset[a]));}
            for(int b:{K(),Epsilon()})for(int j=0;j<8;++j)r.push_back(origin.inner || (b==Epsilon() && p.dependentTrace)?0.:GeometryPoint(b,row.geometry,p.coordinate).Value(all[nv+(b==K()?0:8)+j],offset[b]));
            Require(r.size()==width,"H2e record width");points.write(reinterpret_cast<const char*>(r.data()),r.size()*sizeof(double));
        }
        int family;if(families.count(p.inventory))family=families[p.inventory];else{family=familyNames.size();families[p.inventory]=family;familyNames.push_back(J::Quote(p.inventory));}
        const std::int32_t m[]{total++,id,family,row.geometry,row.element,origin.radialSide,origin.tangentialSide,p.block,p.dependentTrace};membership.write(reinterpret_cast<const char*>(m),sizeof(m));
    }
    points.close();membership.close();Require(points.good() && membership.good(),"H2e observation output failed");
    files.push_back({"observations",J::Object({{"file",J::Quote("observations.f64")},{"shape",J::Numbers(std::vector<int>{n,width})},{"dtype",J::Quote("<f8")},{"sha256",J::Quote(RANSWallR3::HashFile(path+"/observations.f64"))}})});
    files.push_back({"membership",J::Object({{"file",J::Quote("membership.i32")},{"shape",J::Numbers(std::vector<int>{total,9})},{"dtype",J::Quote("<i4")},{"sha256",J::Quote(RANSWallR3::HashFile(path+"/membership.i32"))}})});
    return J::Object({{"version",J::Quote("NativeOuterScalarContextV1")},{"base_sha256",J::Quote(RANSWallR3::Hash(s))},{"files",J::Object(files)},{"families",J::List(familyNames)},
        {"rows",J::Object({{"velocity",J::Numbers(allowed.velocity)},{"k",J::Numbers(allowed.k)},{"epsilon",J::Numbers(allowed.epsilon)}})},{"homogeneous_continuity_max",J::Number(continuity)},{"frozen_direction_max",J::Number(frozen)},
        {"volume",J::Number(volume)},{"outer_volume",J::Number(outerVolume)},{"base_divergence_L2",J::Number(std::sqrt(div[0][0]))},{"base_maximum_quadrature_divergence",J::Number(maxDiv)},
        {"seconds",J::Number(std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count())},{"governing_solves","0"},{"physical_commits","0"}});
}

}
