// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
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
        else if(z<=c[1]+1e-14)return derivative?(c[3]+t*(2*c[4]+3*t*c[5]))/h:c[2]+t*(c[3]+t*(c[4]+t*c[5]));
    }
    if(derivative<0)return integral;
    throw std::runtime_error("seed radial evaluation outside represented panels");
}

RANSVector RANSWallEvolutionPatch::CompositeSeedVelocity(const RANSWallCoordinate &raw) const
{
    const auto p=CompositeChartCoordinate(raw);
    if(m_compositeFlatSeed.empty())Require(false,"composite seed target needs its explicit radial setup projection");
    const double H=m_parameters.geometry.height,B=.4,z=Physical(p)[2];
    auto value=[&](double r){return SeedRadial(0,r);};
    auto derivative=[&](double r){return SeedRadial(0,r,1);};
    const auto g=m_composite->Inner(p[0],std::min(p[2],H));
    if(p[2]<=H)return {value(p[2])*g.tangent1[0],0,value(p[2])*g.tangent1[2]};
    if(z>=B)return {value(z),0,0};
    auto integral=[&](double top){return SeedRadial(0,top,-1);};
    const double h=g.point[2],hs=g.jacobian[2][0],dx=g.jacobian[0][0],delta=B-h,t=(z-h)/delta;
    const double U=value(h),UH=value(H),C=integral(H)-integral(h),D=UH*g.normal[2]-U;
    // n_z'= -slope*slope'/(1+slope^2)^(3/2), with signed curvature
    // kappa=-f''/(1+f'^2)^(3/2). Thus n_z'=slope*kappa.
    const double slope=g.tangent1[2]/g.tangent1[0],nz_s=slope*g.kappa1;
    const double Cs=-U*hs,Ds=UH*nz_s-derivative(h)*hs,ts=-hs*(1-t)/delta;
    const double h00=1-3*t*t+2*t*t*t,h10=t-2*t*t+t*t*t,d00=-6*t+6*t*t,d10=1-4*t+3*t*t;
    const double u=value(z)+C*d00/delta+D*d10;
    const double psi_s=Cs*h00+C*d00*ts+(-hs*D+delta*Ds)*h10+delta*D*d10*ts;
    return {u,0,-psi_s/dx};
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::CompositeSeed() const
{
    Require(Composite() && m_finiteInflow && AcceptedAWF() && OpenBoundary(),"mapped seed requires the composite finite/open accepted-AWF contract");
    RANSWallTiming::Scope timing("startup.composite_seed");m_compositeFlatSeed=FlatFiniteSeed();
    const double H=m_parameters.geometry.height,B=.4;auto s=m_compositeFlatSeed;
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
            (axis==2 && !side && b==Epsilon() && el.lower[2]==H);});
        for(auto [id,v]:boundary)fixed[offset[b]+id]=m_compositeFlatSeed[offset[b]+id];
    }
    // Keep the intended inner budget-reference k exactly. This is seed setup,
    // not a prescribed ground-k condition in physical evolution.
    for(const auto &q:Quadrature(K(),4))if(q.coordinate[2]<H)for(int id:q.ids)fixed[offset[K()]+id]=m_compositeFlatSeed[offset[K()]+id];
    for(int b=0;b<2;++b)for(int e=0;e<fields[b]->Elements().size();++e)
        if(fields[b]->Elements()[e].lower[2]>=B)for(int id:fields[b]->ElementPoint(e,fields[b]->Elements()[e].lower).ids)
            fixed[offset[b]+id]=(b==0?FixedInflow().AtHeight(B).U:0.)*m_constants[b][id];
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
    const Factor M(mass);WallMatrix inverse,schur(Q.size(),Vector(Q.size()));
    for(const auto &row:Q)inverse.push_back(M.Apply(row));
    for(int i=0;i<Q.size();++i)for(int j=0;j<Q.size();++j)schur[i][j]=std::inner_product(Q[i].begin(),Q[i].end(),inverse[j].begin(),0.);
    const auto multiplier=RANSWallAlgebra::Solve(schur,beta);Vector delta(free.size());
    for(int j=0;j<Q.size();++j)for(int i=0;i<free.size();++i)delta[i]+=inverse[j][i]*multiplier[j];
    for(int i=0;i<free.size();++i)s[free[i]]+=delta[i];
    auto remaining=Multiply(D,delta);
    for(int i=0;i<D.size();++i)remaining[i]+=defect[i];
    Require(Maximum(remaining)<1e-10,"mapped seed native continuity projection failed");
    m_compositeSeedInfo=J::Object({{"method",J::Quote("physical streamfunction lift; unstrained tail above z=.4; physical consistent-mass constrained native projection")},
        {"setup_radial_projections","3"},{"setup_field_projections","2"},{"setup_continuity_mass_factors","1"},
        {"setup_physical_projection_mass_solves","6"},{"setup_fixed_constraint_mass_solves","4"},
        {"setup_continuity_schur_solves","1"},{"independent_continuity_rows",J::Number(Q.size())},
        {"native_weak_continuity_max",J::Number(Maximum(remaining))},{"velocity_projection_mass_norm",J::Number(std::sqrt(std::max(0.,std::inner_product(delta.begin(),delta.end(),Multiply(mass,delta).begin(),0.))))},
        {"governing_solves","0"},{"physical_steps","0"}});return s;
}

std::string RANSWallEvolutionPatch::CompositeChecks(const Vector &s,double time,const RANSWallEvolutionData &data) const
{
    Require(AcceptedAWF() && FiniteInflowReference() && OpenBoundary(),"H2b operators require the declared finite/open budget AWF");
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
}
