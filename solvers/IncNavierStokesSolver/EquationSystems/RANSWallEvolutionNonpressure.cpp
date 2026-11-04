///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>

namespace Nektar
{
using namespace SolverUtils;
namespace
{
using Vector=RANSWallAlgebra::Vector;
using RANSWallAlgebra::Require;
using Coordinate=RANSWallCoordinate;
using Sparse=std::map<int,long double>;
using Clock=std::chrono::steady_clock;

// A local scalar value/direction pair. This is private component algebra,
// not a replacement value model or a whole-residual differencing mechanism.
struct D
{
    long double v=0,d=0;
    D()=default;
    D(long double value,long double direction=0):v(value),d(direction){}
};
D operator+(D a,D b){return {a.v+b.v,a.d+b.d};}
D operator-(D a,D b){return {a.v-b.v,a.d-b.d};}
D operator-(D a){return {-a.v,-a.d};}
D operator*(D a,D b){return {a.v*b.v,a.d*b.v+a.v*b.d};}
D operator/(D a,D b){return {a.v/b.v,(a.d-a.v/b.v*b.d)/b.v};}
D Sqrt(D a)
{
    Require(a.v>=0,"R3 square root outside its base domain");
    const long double r=std::sqrt(a.v);
    return {r,r>0?a.d/(2*r):0}; // selected zero norm subgradient at the origin
}
D Exp(D a){const long double e=std::exp(a.v);return {e,e*a.d};}
D Expm1(D a){return {std::expm1(a.v),std::exp(a.v)*a.d};}
D Log(D a){return {std::log(a.v),a.d/a.v};}
D Pow(D a,long double p){const long double v=std::pow(a.v,p);return {v,p*std::pow(a.v,p-1)*a.d};}
D Max(D a,D b){return a.v>=b.v?a:b;}
D Min(D a,D b){return a.v<=b.v?a:b;}
D PositivePart(D a){return a.v>0?a:D{};}
double V(D a){return static_cast<double>(a.v);}
double T(D a){return static_cast<double>(a.d);}
using DT=std::array<std::array<D,3>,3>;
DT Product(const DT &a,const DT &b)
{
    DT result{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)for(int k=0;k<3;++k)
        result[i][j]=result[i][j]+a[i][k]*b[k][j];
    return result;
}
struct OuterAction
{
    RANSTensor stress{};double viscosity=0,kSource=0,epsilonSource=0;
};
OuterAction OuterDerivative(double k,double epsilon,const RANSTensor &g,double r,
    double dk,double de,const RANSTensor &dg,const RANSParameters &p,
    bool cls,bool strainBranch,const RANSRobustVariant &variant)
{
    D K{k,dk},E{epsilon,de},time=K/E,ss{},oo{},div{};DT G{},s{},o{};
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)G[i][j]={g[i][j],dg[i][j]};
    for(int i=0;i<3;++i)
    {
        div=div+G[i][i];
        for(int j=0;j<3;++j)
        {
            s[i][j]=time*(G[i][j]+G[j][i]);o[i][j]=time*(G[i][j]-G[j][i]);
            ss=ss+s[i][j]*s[i][j];oo=oo+o[i][j]*o[i][j];
        }
    }
    D c{p.cMu};
    if(cls)
    {
        const D q=Sqrt((strainBranch?ss:oo)/D{2});
        const D factor=q.v>8?D{1}:-Expm1(-D{.36L}*Exp(D{.75L}*q));
        if(q.v>1)
        {const D reciprocal=D{1}/q/Sqrt(q);c=D{.3L}*factor*reciprocal/(reciprocal+D{.35L});}
        else c=D{.3L}*factor/(D{1}+D{.35L}*q*Sqrt(q));
    }
    const D nt=c*K*time;DT quadratic{},cubic{},stress{};
    if(cls)
    {
        const auto s2=Product(s,s),o2=Product(o,o),os=Product(o,s),so=Product(s,o);
        const auto s2o=Product(s2,o),os2=Product(o,s2);
        const D c4=-D{10}*c*c,c6=-D{5}*c*c,c7=D{5}*c*c;
        for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        {
            const D delta{i==j?1.L:0.L};
            quadratic[i][j]=K*c*(-D{.1L}*(s2[i][j]-ss*delta/D{3})+
                D{.1L}*(os[i][j]-so[i][j])+D{.26L}*(-o2[i][j]-oo*delta/D{3}));
            cubic[i][j]=K*c*(c4*(s2o[i][j]-os2[i][j])+(c6*ss+c7*oo)*s[i][j]);
        }
    }
    D production{};OuterAction result;
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {
        // Form the turbulent mean stress directly, before adding reduced pi.
        // Isotropic k cancels here, but remains in -R:G for production.
        stress[i][j]=K*c*s[i][j]-quadratic[i][j]-D{variant.beta}*cubic[i][j]+
            D{variant.gamma}*(G[i][j]+G[j][i]-D{variant.referenceStrain[i][j]});
        production=production+(stress[i][j]-(i==j?D{2}*K/D{3}:D{}))*G[i][j];
        result.stress[i][j]=T(stress[i][j]);
    }
    D yap{};const D scale=E*E/K;
    if(cls)
    {
        const D ratio=K*Sqrt(K)/(D{p.yapLength*r}*E);
        // One fixed generalized derivative: inactive at the Yap knot.
        yap=scale*D{p.yapCoefficient}*PositivePart(ratio-D{1})*ratio*ratio;
    }
    result.viscosity=T(nt);result.kSource=T(production-E);
    result.epsilonSource=T(D{p.cEpsilon1}*E/K*production-D{p.cEpsilon2}*scale+yap);
    return result;
}

struct QD {D r,w;};
struct PreD {double nt=0,E=0,Dpre=0,ntK=0,EK=0,DK=0,Echi=0,Dchi=0;};
class ProfileDerivative
{
public:
    explicit ProfileDerivative(const EquilibriumAWFProfiles &profile,const RANSNeutralReference &model)
        :p(profile),Q(model.Q()),kappa(model.Kappa()),K(profile.K(),1)
    {
        Require(p.Metric().kappa1==0 && p.Metric().kappa2==0,"R3 profile differential requires the flat ray metric");
        const auto &g=p.Geometry();const auto &w=p.Original().Parameters();
        const D root=Sqrt(K);a=D{w.cMu}*D{w.cLength}*root;
        D intercept{w.smoothIntercept};
        if(g.roughness>0)
        {
            const D x=D{g.roughness}*root/(D{70}*D{g.molecularNu});
            const D exponent=Max(D{.5L}-D{.4L}*Pow(x,.7L),D{1}-D{.79L}*Pow(x,-.28L));
            intercept=-D{w.smoothIntercept}*Expm1(exponent*Log(x));
        }
        zv=D{g.molecularNu}*intercept/root;ze=D{2}*D{w.cLength}*D{g.molecularNu}/root;
        zv.v=p.Original().ViscosityInterceptHeight();ze.v=p.Original().DissipationCutoff();
        cutoff={p.Reference().cutoff,0};Gauss();
        if(p.Reference().weight>0)
        {
            const D partialK=Budget(K,cutoff),partialR=Budget(D{K.v},D{cutoff.v,1});
            Require(partialR.d!=0 && std::isfinite(partialR.d),"R3 implicit cutoff derivative is singular");
            cutoff.d=-partialK.d/partialR.d;
        }
        rule=Rule(K,p.Reference().weight>0?cutoff:D{});
        Require(rule.size()==p.Quadrature().size(),"R3 moving quadrature topology differs from the native profile");
        for(std::size_t i=0;i<rule.size();++i)
        {
            const auto &native=p.Quadrature()[i];
            Require(std::abs(V(rule[i].r)-native.r)<=2e-13*std::max(1.,native.r) &&
                std::abs(V(rule[i].w)-native.weight)<=2e-13*std::max(1.,native.weight),"R3 moving rule does not reproduce the native profile");
            rule[i].r.v=native.r;rule[i].w.v=native.weight;
        }
        D production{},minimum{};
        for(const auto &q:Rule(K,D{}))
        {production=production+q.w*P0(K,q.r);minimum=minimum+q.w*Corrected(K,q.r,D{g.tailStart});}
        const D margin=(production-minimum)/minimum;
        marginK=T(margin);reynoldsK=g.height/(2*g.molecularNu*std::sqrt(p.K()));
    }
    const std::vector<QD> &Quadrature() const{return rule;}
    double MarginK() const{return marginK;}
    double ReynoldsK() const{return reynoldsK;}
    double CutoffK() const{return T(cutoff);}
    PreD At(double radius,double radiusK=0) const
    {
        const D r{radius,radiusK},nt=Viscosity(K,r),original=Original(K,r),production=P0(K,r);
        const double chi=p.Reference().weight;D e=original,d=original;
        double ec=0;
        if(chi>0)
        {
            const D corrected=Corrected(K,r,cutoff);
            e=D{1-chi}*original+D{chi}*corrected;
            d=D{1-chi}*original+D{chi}*production;ec=V(corrected-original);
        }
        const auto native=p.AtProduction(radius,0);
        return {p.Original().EddyViscosity(radius),native.dissipation,native.dissipation-native.redistribution,
            T(nt),T(e),T(d),ec,V(production-original)};
    }
private:
    D Viscosity(D k,D r) const
    {
        const D z{zv.v,zv.d*k.d};
        return D{a.v,a.d*k.d}*PositivePart(r-z);
    }
    D Original(D k,D r) const
    {return k*Sqrt(k)/D{p.Original().Parameters().cLength}/Max(r,D{ze.v,ze.d*k.d});}
    D P0(D k,D r) const
    {
        const D nt=Viscosity(k,r),tau=k/D{Q},nu{p.Geometry().molecularNu};
        return nt*tau*tau/((nu+nt)*(nu+nt));
    }
    D Corrected(D k,D r,D rc) const
    {return Pow(k/D{Q},1.5L)/D{kappa}/Max(r,rc);}
    void Gauss()
    {
        const int n=p.Options().quadratureOrder;nodes.resize(n);weights.resize(n);
        for(int j=0;j<(n+1)/2;++j)
        {
            long double x=std::cos(std::acos(-1.L)*(j+.75L)/(n+.5L));
            for(int it=0;it<100;++it)
            {
                long double value=1,previous=0;
                for(int k=1;k<=n;++k){const long double next=((2*k-1)*x*value-(k-1)*previous)/k;previous=value;value=next;}
                const long double derivative=n*(x*value-previous)/(x*x-1),delta=value/derivative;x-=delta;
                if(std::abs(delta)<=4*std::numeric_limits<long double>::epsilon())break;
            }
            long double value=1,previous=0;
            for(int k=1;k<=n;++k){const long double next=((2*k-1)*x*value-(k-1)*previous)/k;previous=value;value=next;}
            const long double derivative=n*(x*value-previous)/(x*x-1);
            nodes[j]=static_cast<double>(-x);nodes[n-1-j]=static_cast<double>(x);
            weights[j]=weights[n-1-j]=static_cast<double>(2/((1-x*x)*derivative*derivative));
        }
    }
    std::vector<QD> Rule(D k,D rc) const
    {
        const auto &g=p.Geometry();const double H=g.height;std::vector<D> breaks{{0},{H}};
        auto add=[&](D value)
        {if(value.v>0 && value.v<H){value.v=static_cast<double>(value.v);if(value.v>0)breaks.push_back(value);}};
        for(double r:p.Options().physicalBreakpoints)add(D{r});
        for(D r:{D{g.roughness},D{zv.v,zv.d*k.d},D{ze.v,ze.d*k.d},D{g.tailStart},rc})add(r);
        const auto &w=p.Original().Parameters();D seed=Min(D{g.tailStart},D{g.molecularNu}/(D{w.cMu*w.cLength}*Sqrt(k)));
        if(rc.v>0)seed=Min(seed,rc);
        for(int j=0;j<4096 && seed.v<H;++j,seed=seed*D{2})add(seed);
        std::stable_sort(breaks.begin(),breaks.end(),[](D x,D y){return x.v<y.v;});
        breaks.erase(std::unique(breaks.begin(),breaks.end(),[](D x,D y){return x.v==y.v;}),breaks.end());
        std::vector<QD> result;
        for(std::size_t j=1;j<breaks.size();++j)
        {
            const D half=(breaks[j]-breaks[j-1])/D{2},mid=breaks[j-1]+half;
            for(int i=0;i<p.Options().quadratureOrder;++i)
            {
                auto r=mid+half*D{nodes[i]},weight=half*D{weights[i]};
                r.v=static_cast<double>(r.v);weight.v=static_cast<double>(weight.v);
                result.push_back({r,weight});
            }
        }
        return result;
    }
    D Budget(D k,D rc) const
    {D result{};for(const auto &q:Rule(k,rc))result=result+q.w*(Corrected(k,q.r,rc)-P0(k,q.r));return result;}
    const EquilibriumAWFProfiles &p;double Q,kappa,marginK=0,reynoldsK=0;
    D K,a,zv,ze,cutoff;
    std::vector<double> nodes,weights;std::vector<QD> rule;
};

struct Basis
{
    RANSWallGlobalPoint row;Vector value;std::array<Vector,3> gradient,hessian;
};
Basis LinearBasis(const RANSWallGlobalPoint &row,const RANSWallGlobalField &field,bool hessian)
{
    Basis result;result.row=row;result.value=row.phi;result.gradient=row.derivative;
    if(hessian)result.hessian=field.RadialGradientDerivative(row);
    for(auto &values:result.hessian)if(values.empty())values.assign(row.ids.size(),0.);
    const long double inverse=1.L/row.constant[0];long double sum=0;
    for(std::size_t j=0;j<row.ids.size();++j)sum+=static_cast<long double>(row.phi[j])*row.constant[j];
    result.value[0]+=static_cast<double>((1-sum)*inverse);
    for(int a=0;a<3;++a)
    {
        sum=0;for(std::size_t j=0;j<row.ids.size();++j)sum+=static_cast<long double>(row.derivative[a][j])*row.constant[j];
        result.gradient[a][0]-=static_cast<double>(sum*inverse);
    }
    return result;
}
RANSVector RadialGradient(const Basis &basis,const Vector &s,int offset)
{
    const auto &row=basis.row;const double centre=s[offset+row.ids[0]]/row.constant[0];RANSVector result{};
    for(int a=0;a<3;++a)
    {long double sum=0;for(std::size_t j=0;j<row.ids.size();++j)sum+=static_cast<long double>(basis.hessian[a][j])*(s[offset+row.ids[j]]-centre*row.constant[j]);result[a]=static_cast<double>(sum);}
    return result;
}
void Append(Sparse &target,const Sparse &source,long double factor)
{if(factor!=0)for(auto [i,value]:source)target[i]+=factor*value;}
Sparse ValueRow(const Basis &basis,int offset)
{Sparse result;for(std::size_t j=0;j<basis.row.ids.size();++j)result[offset+basis.row.ids[j]]+=basis.value[j];return result;}
struct LimitValue {double value=0,slope=0;};
LimitValue Limit(double ratio,const RANSCappedDestructionControls &controls)
{
    if(ratio>=controls.theta)return {1,0};
    if(!controls.smooth)return {ratio/controls.theta,1/controls.theta};
    if(ratio<=0)return {0,0};
    const double t=ratio/controls.theta,one=1-t;
    return {t*t*t*(10+t*(-15+6*t)),30*t*t*one*one/controls.theta};
}
struct Node
{
    Basis k;
    PreD pre;
    double radius=0,radiusK=0,weight=0,weightK=0,kValue=0,kRadial=0;
    LimitValue limiter;
};
struct RayLinear
{
    std::shared_ptr<EquilibriumAWFProfiles> profile;
    std::unique_ptr<ProfileDerivative> differential;
    RANSCappedDestructionValues normalization;
    Sparse sampled,chi,b,m,lambda;
    std::vector<Node> nodes;
};
struct FieldState
{
    Basis v,p,e;
    RANSVector u{},gk{},ge{},uz{},gkz{};
    RANSTensor G{},Gz{},stress{};
    double k=0,epsilon=0,pi=0,piz=0,dk=0,de=0,sk=0,se=0,weight=0;
};
struct Variation
{
    RANSVector u{},gk{},ge{};RANSTensor G{};
    double k=0,epsilon=0,pi=0,nt=0,preD=0,sampledK=0,radius=0,weight=0;
};
struct Loads
{
    std::map<int,long double> value;
    void Add(int offset,const Basis &test,double source,const RANSVector &flux,
        double dsource,const RANSVector &dflux,double weight,double dr,double dw)
    {
        const auto &q=test.row;
        for(std::size_t j=0;j<q.ids.size();++j)
        {
            long double base=static_cast<long double>(q.phi[j])*source;
            long double action=static_cast<long double>(q.phi[j])*dsource+
                static_cast<long double>(dr)*q.derivative[2][j]*source;
            for(int a=0;a<3;++a)
            {
                base+=static_cast<long double>(q.derivative[a][j])*flux[a];
                action+=static_cast<long double>(q.derivative[a][j])*dflux[a]+
                    static_cast<long double>(dr)*test.hessian[a][j]*flux[a];
            }
            value[offset+q.ids[j]]+=weight*action+dw*base;
        }
    }
};
} // namespace

RANSWallEvolutionPatch::WallMatrix RANSWallEvolutionPatch::AssembleNonpressureR3(
    const Vector &state,double time,const RANSWallEvolutionData &data,bool imposeConstraints,
    int region,std::vector<std::string> *metadata) const
{
    const auto started=Clock::now();
    Require(m_dimension==2 && !Mapped() && m_transitionEnabled && m_cappedDestruction && Size()<=3000,
        "R3 nonpressure actions require the declared flat2D Upstream+Capped dense patch");
    Require(m_comm->GetSize()==1,"R3 nonpressure dense action currently requires one rank");
    Require(state.size()==static_cast<std::size_t>(Size()) && region>=-1 && region<=1,
        "invalid R3 nonpressure state/region");
    const auto constraintStarted=Clock::now();
    const auto constraints=Constraints(state,time,data);const auto s=ExternalLift(state,time,data);
    const double constraintSeconds=std::chrono::duration<double>(Clock::now()-constraintStarted).count();
    std::vector<std::vector<long double>> jac(Size(),std::vector<long double>(DifferentialSize()));
    std::size_t localKernels=0,volumePoints=0,movingPoints=0,normalizations=0,profileCount=0,traceSamples=0;
    std::string outerBranches,yapBranches;
    auto scatter=[&](const Loads &load,int column,long double factor=1)
    {if(column<DifferentialSize() && factor!=0)for(auto [row,value]:load.value)jac[row][column]+=factor*value;};
    auto rank=[&](const Loads &load,const Sparse &row)
    {for(auto [column,value]:row)scatter(load,column,value);};
    auto basis=[&](int block,const Coordinate &p,int side=0,bool hessian=false)
    {return LinearBasis(Point(block,p,side),*fields[block],hessian);};

    struct ForceRows {RANSVector pressure{},convection{};std::array<Sparse,3> dpressure,dconvection,total;};
    std::map<Coordinate,ForceRows> forceRows;
    auto forcing=[&](const Coordinate &chart)->const ForceRows&
    {
        auto H=chart;H[2]=m_parameters.geometry.height;auto found=forceRows.find(H);
        if(found!=forceRows.end())return found->second;
        ForceRows f;const auto parts=ForcingParts(s,H,time,data);f.pressure=parts[0];f.convection=parts[1];
        const auto v=basis(0,H,1),k=basis(K(),H,1);RANSVector u{};RANSTensor G{};
        for(int a=0;a<m_dimension;++a){u[PhysicalAxis(a)]=v.row.Value(s,offset[a]);G[PhysicalAxis(a)]=v.row.Gradient(s,offset[a]);}
        for(std::size_t j=0;j<k.row.ids.size();++j)f.dpressure[0][offset[K()]+k.row.ids[j]]-=2*k.gradient[0][j]/3;
        for(int component=0;component<m_dimension;++component)for(std::size_t j=0;j<v.row.ids.size();++j)
        {
            // Flat tangent force is (u_x d_x)u_x+u_x div_Gamma u.
            // In2D the tangent projector has only x; wall-normal velocity is
            // retained by matching/sensor samples but not tangent convection.
            const int a=PhysicalAxis(component),column=offset[component]+v.row.ids[j];
            if(a==0)f.dconvection[0][column]+=2*(v.value[j]*G[0][0]+u[0]*v.gradient[0][j]);
        }
        for(int a=0;a<3;++a){f.total[a]=f.dpressure[a];Append(f.total[a],f.dconvection[a],1);}
        return forceRows.emplace(H,std::move(f)).first->second;
    };

    auto chiRow=[&](const Coordinate &chart,const ProfileDerivative &base)->Sparse
    {
        Sparse result;
        if(m_parameters.chi==0 || WallCoordinate(chart)>=m_transition.end)return result;
        const auto input=TransitionInput(s,chart,time,data);
        if(input.budgetMargin<=m_transition.state.marginOff || input.reynolds<=m_transition.state.reynoldsOff)return result;
        auto H=chart;H[2]=m_parameters.geometry.height;auto a=H,b=H,sample=H;
        a[2]=m_sensorA;b[2]=m_sensorB;sample[2]=m_parameters.geometry.sampleHeight;
        const auto vh=basis(0,H,1),va=basis(0,a),vb=basis(0,b),ka=basis(K(),a),ea=basis(Epsilon(),a),ks=basis(K(),sample);
        AWFTransitionPhysicalInput physical;physical.normal={0,0,1};physical.sampledK=SampleK(s,H);
        physical.height=m_parameters.geometry.height;physical.molecularNu=m_parameters.geometry.molecularNu;physical.Q=m_parameters.model.Q();
        for(int j=0;j<m_dimension;++j)
        {
            const int axis=PhysicalAxis(j);
            physical.matchingVelocity[axis]=vh.row.Value(s,offset[j]);physical.matchingReference[axis]=vh.row.Value(m_transitionReference,offset[j]);
            physical.velocityA[axis]=va.row.Value(s,offset[j]);physical.velocityAReference[axis]=va.row.Value(m_transitionReference,offset[j]);
            physical.velocityB[axis]=vb.row.Value(s,offset[j]);physical.velocityBReference[axis]=vb.row.Value(m_transitionReference,offset[j]);
        }
        physical.kA=ka.row.Value(s,offset[K()]);physical.epsilonA=ea.row.Value(s,offset[Epsilon()]);
        physical.epsilonAReference=ea.row.Value(m_transitionReference,offset[Epsilon()]);
        const auto &f=forcing(H);physical.pressure=f.pressure;physical.convection=f.convection;
        RANSWallEvolutionData zero(m_dimension);const auto reference=ForcingParts(m_transitionReference,H,time,zero);
        physical.pressureReference=reference[0];physical.convectionReference=reference[1];
        struct PhysicalDirection {AWFTransitionPhysicalInput input;};
        std::map<int,PhysicalDirection> directions;
        for(int j=0;j<m_dimension;++j)
        {
            const int axis=PhysicalAxis(j);
            for(std::size_t i=0;i<vh.row.ids.size();++i)directions[offset[j]+vh.row.ids[i]].input.matchingVelocity[axis]+=vh.value[i];
            for(std::size_t i=0;i<va.row.ids.size();++i)directions[offset[j]+va.row.ids[i]].input.velocityA[axis]+=va.value[i];
            for(std::size_t i=0;i<vb.row.ids.size();++i)directions[offset[j]+vb.row.ids[i]].input.velocityB[axis]+=vb.value[i];
        }
        for(std::size_t i=0;i<ka.row.ids.size();++i)directions[offset[K()]+ka.row.ids[i]].input.kA+=ka.value[i];
        for(std::size_t i=0;i<ks.row.ids.size();++i)directions[offset[K()]+ks.row.ids[i]].input.sampledK+=ks.value[i];
        for(std::size_t i=0;i<ea.row.ids.size();++i)directions[offset[Epsilon()]+ea.row.ids[i]].input.epsilonA+=ea.value[i];
        for(int axis=0;axis<3;++axis)
        {
            for(auto [column,value]:f.dpressure[axis])directions[column].input.pressure[axis]+=static_cast<double>(value);
            for(auto [column,value]:f.dconvection[axis])directions[column].input.convection[axis]+=static_cast<double>(value);
        }
        for(const auto &[column,direction]:directions)
        {
            AWFTransitionInput dinput;dinput.budgetMargin=base.MarginK()*direction.input.sampledK;
            dinput.reynolds=base.ReynoldsK()*direction.input.sampledK;
            dinput.sensors=EvaluateEquilibriumAWFTransitionSensors(physical,direction.input).derivative;
            const double value=EvaluateEquilibriumAWFTransition(input,m_transition,dinput).derivative.chi;
            if(value!=0)result[column]=value;
        }
        return result;
    };

    std::map<std::pair<Coordinate,bool>,std::shared_ptr<RayLinear>> cache;
    auto rayLinear=[&](const Coordinate &chart,bool baseline)->const RayLinear&
    {
        auto key=chart;key[2]=0;const auto cacheKey=std::make_pair(key,baseline);auto found=cache.find(cacheKey);
        if(found!=cache.end())return *found->second;
        auto ray=std::make_shared<RayLinear>();const double K0=SampleK(s,chart);
        auto sample=chart;sample[2]=m_parameters.geometry.sampleHeight;
        ray->sampled=ValueRow(basis(K(),sample),offset[K()]);
        const auto base=Profile(K0,chart,0);
        ProfileDerivative baseDerivative(*base,m_parameters.model);++profileCount;
        ray->chi=baseline?Sparse{}:chiRow(chart,baseDerivative);
        ray->profile=baseline?base:StateProfile(s,chart,time,data);
        ray->differential=std::make_unique<ProfileDerivative>(*ray->profile,m_parameters.model);++profileCount;
        RANSCappedDestructionInput input;input.sampledK=K0;
        const auto &rule=ray->profile->Quadrature();const auto &dRule=ray->differential->Quadrature();
        for(std::size_t i=0;i<rule.size();++i)
        {
            auto point=chart;point[2]=rule[i].r;Node n;n.radius=rule[i].r;n.radiusK=T(dRule[i].r);
            n.weight=rule[i].weight;n.weightK=T(dRule[i].w);n.k=basis(K(),point);
            n.kValue=n.k.row.Value(s,offset[K()]);n.kRadial=n.k.row.Gradient(s,offset[K()])[2];
            n.pre=ray->differential->At(n.radius,n.radiusK);n.limiter=Limit(n.kValue/K0,m_destruction);
            input.k.push_back(n.kValue);input.preE.push_back(n.pre.E);input.preD.push_back(n.pre.Dpre);input.metricWeights.push_back(n.weight);
            ray->nodes.push_back(std::move(n));
        }
        ray->normalization=EvaluateRANSCappedDestruction(input,m_destruction);++normalizations;
        const long double B=ray->normalization.budgetD,b=ray->normalization.b;long double bk=0,bc=0;
        for(const auto &n:ray->nodes)
        {
            const long double common=n.weight*n.pre.Dpre;
            bk+=(static_cast<long double>(n.weightK)*n.pre.Dpre+static_cast<long double>(n.weight)*n.pre.DK)*(n.limiter.value-b)+
                common*n.limiter.slope*(n.kRadial*n.radiusK/K0-n.kValue/(K0*K0));
            bc+=static_cast<long double>(n.weight)*n.pre.Dchi*(n.limiter.value-b);
            const long double factor=common*n.limiter.slope/(K0*B);
            for(std::size_t j=0;j<n.k.row.ids.size();++j)ray->b[offset[K()]+n.k.row.ids[j]]+=factor*n.k.value[j];
        }
        Append(ray->b,ray->sampled,bk/B);Append(ray->b,ray->chi,bc/B);
        if(ray->normalization.b>=m_destruction.delta)ray->m=ray->b;
        else Append(ray->lambda,ray->b,1/m_destruction.delta);
        return *cache.emplace(cacheKey,std::move(ray)).first->second;
    };

    auto stateAt=[&](const RANSWallGlobalPoint &v,const RANSWallGlobalPoint &p,
        const RANSWallGlobalPoint *epsilon,double weight,bool moving)->FieldState
    {
        FieldState q;q.v=LinearBasis(v,*fields[0],moving);q.p=LinearBasis(p,*fields[epsilon?OuterPressure():InnerPressure()],moving);
        if(epsilon)q.e=LinearBasis(*epsilon,*fields[Epsilon()],false);
        for(int a=0;a<m_dimension;++a)
        {
            const int axis=PhysicalAxis(a);q.u[axis]=v.Value(s,offset[a]);q.G[axis]=v.Gradient(s,offset[a]);
            q.uz[axis]=q.G[axis][2];if(moving)q.Gz[axis]=RadialGradient(q.v,s,offset[a]);
        }
        q.k=v.Value(s,offset[K()]);q.gk=v.Gradient(s,offset[K()]);if(moving)q.gkz=RadialGradient(q.v,s,offset[K()]);
        const int pressure=epsilon?OuterPressure():InnerPressure();q.pi=p.Value(s,offset[pressure]);q.piz=p.Gradient(s,offset[pressure])[2];
        if(epsilon){q.epsilon=epsilon->Value(s,offset[Epsilon()]);q.ge=epsilon->Gradient(s,offset[Epsilon()]);}
        q.weight=weight;return q;
    };

    if(region!=1)for(const auto &nativeRay:m_rays)
    {
        const auto &ray=rayLinear(nativeRay.coordinate,false);const auto &f=forcing(nativeRay.coordinate);
        const auto force=CartesianForcing(s,nativeRay.coordinate,time,data,false);
        Loads sampledAction,chiAction,normalizationAction;std::array<Loads,3> forceAction;
        for(const auto &n:ray.nodes)
        {
            ++volumePoints;if(n.radiusK!=0 || n.weightK!=0)++movingPoints;
            auto point=nativeRay.coordinate;point[2]=n.radius;
            const auto v=Point(0,point),pressure=Point(InnerPressure(),point);
            Require(v.ids==n.k.row.ids,"R3 k/velocity native maps differ");
            auto q=stateAt(v,pressure,nullptr,nativeRay.weight*n.weight,true);
            Require(data.innerKSource(point,time)==0,"R3 moving-source action requires the declared unforced inner scalar source");
            const auto inner=InnerStress(q.k,n.pre.nt,q.G);
            q.stress=RANSWallInnerMeanStress(inner,m_parameters.geometry.molecularNu,q.pi);
            q.dk=m_parameters.geometry.molecularNu+n.pre.nt/m_parameters.model.Parameters().sigmaK;
            const double m=ray.normalization.m,K0=ray.profile->K();
            const auto source=EvaluateRANSCappedDestructionPoint(q.k,K0,n.pre.E,n.pre.Dpre,
                ray.normalization.b,m,m_destruction);
            q.sk=inner.production-source.destruction;
            std::array<RANSWallInnerStressAction,6> partial{};RANSTensor zero{};
            partial[0]=RANSWallInnerStressDirectional(q.k,n.pre.nt,q.G,1,0,zero,m_innerStress);
            partial[1]=RANSWallInnerStressDirectional(q.k,n.pre.nt,q.G,0,1,zero,m_innerStress);
            for(int component=0;component<2;++component)for(int axis=0;axis<2;++axis)
            {
                auto h=zero;h[PhysicalAxis(component)][PhysicalAxis(axis)]=1;
                partial[2+2*component+axis]=RANSWallInnerStressDirectional(q.k,n.pre.nt,q.G,0,0,h,m_innerStress);
            }
            localKernels+=partial.size();
            auto evaluate=[&](const Variation &d)->Loads
            {
                RANSTensor stress{};double production=partial[0].production*d.k+partial[1].production*d.nt;
                for(int a=0;a<3;++a)for(int b=0;b<3;++b)
                    stress[a][b]=partial[0].turbulentStress[a][b]*d.k+partial[1].turbulentStress[a][b]*d.nt+
                        m_parameters.geometry.molecularNu*(d.G[a][b]+d.G[b][a])-(a==b?d.pi:0);
                for(int component=0;component<2;++component)for(int axis=0;axis<2;++axis)
                {
                    const auto &action=partial[2+2*component+axis];const double factor=d.G[PhysicalAxis(component)][PhysicalAxis(axis)];
                    production+=action.production*factor;
                    for(int a=0;a<3;++a)for(int b=0;b<3;++b)stress[a][b]+=action.turbulentStress[a][b]*factor;
                }
                const double destruction=(d.preD*n.limiter.value+n.pre.Dpre*n.limiter.slope*
                    (d.k/K0-q.k*d.sampledK/(K0*K0)))/m;
                const double sk=production-destruction;
                Loads load;
                for(int component=0;component<2;++component)
                {
                    const int a=PhysicalAxis(component);RANSVector flux{},change{};
                    for(int b=0;b<3;++b)
                    {flux[b]=q.stress[a][b]-q.u[a]*q.u[b];change[b]=stress[a][b]-d.u[a]*q.u[b]-q.u[a]*d.u[b];}
                    const double rough=n.radius<=m_parameters.geometry.roughness?force[a]:0;
                    load.Add(offset[component],q.v,-rough,flux,0,change,q.weight,d.radius,d.weight);
                }
                RANSVector kflux{},dkflux{};
                for(int a=0;a<3;++a)
                {
                    kflux[a]=q.dk*q.gk[a]-q.u[a]*q.k;
                    dkflux[a]=q.dk*d.gk[a]+d.nt/m_parameters.model.Parameters().sigmaK*q.gk[a]-d.u[a]*q.k-q.u[a]*d.k;
                }
                load.Add(offset[K()],q.v,-q.sk,kflux,-sk,dkflux,q.weight,d.radius,d.weight);
                const double div=q.G[0][0]+q.G[2][2],ddiv=d.G[0][0]+d.G[2][2];
                load.Add(offset[InnerPressure()],q.p,div,{},ddiv,{},q.weight,d.radius,d.weight);
                return load;
            };
            for(int component=0;component<2;++component)for(std::size_t j=0;j<q.v.row.ids.size();++j)
            {
                Variation d;const int a=PhysicalAxis(component);d.u[a]=q.v.value[j];
                for(int b=0;b<3;++b)d.G[a][b]=q.v.gradient[b][j];
                scatter(evaluate(d),offset[component]+q.v.row.ids[j]);
            }
            for(std::size_t j=0;j<q.v.row.ids.size();++j)
            {
                Variation d;d.k=q.v.value[j];for(int b=0;b<3;++b)d.gk[b]=q.v.gradient[b][j];
                scatter(evaluate(d),offset[K()]+q.v.row.ids[j]);
            }
            Variation kd;kd.sampledK=1;kd.radius=n.radiusK;kd.weight=nativeRay.weight*n.weightK;
            kd.k=q.gk[2]*n.radiusK;kd.pi=q.piz*n.radiusK;kd.nt=n.pre.ntK;kd.preD=n.pre.DK;
            for(int a=0;a<3;++a)
            {
                kd.u[a]=q.uz[a]*n.radiusK;kd.gk[a]=q.gkz[a]*n.radiusK;
                for(int b=0;b<3;++b)kd.G[a][b]=q.Gz[a][b]*n.radiusK;
            }
            const auto kload=evaluate(kd);for(auto [row,value]:kload.value)sampledAction.value[row]+=value;
            Variation cd;cd.preD=n.pre.Dchi;
            const auto cload=evaluate(cd);for(auto [row,value]:cload.value)chiAction.value[row]+=value;
            for(std::size_t j=0;j<q.v.row.ids.size();++j)
            {
                normalizationAction.value[offset[K()]+q.v.row.ids[j]]-=static_cast<long double>(q.weight)*q.v.row.phi[j]*source.destruction/m;
                if(n.radius<=m_parameters.geometry.roughness)for(int component=0;component<2;++component)
                    forceAction[PhysicalAxis(component)].value[offset[component]+q.v.row.ids[j]]-=static_cast<long double>(q.weight)*q.v.row.phi[j];
            }
        }
        rank(sampledAction,ray.sampled);rank(chiAction,ray.chi);rank(normalizationAction,ray.m);
        for(int a=0;a<3;++a)rank(forceAction[a],f.total[a]);
    }

    if(region!=0)for(const auto &native:m_outer)
    {
        ++volumePoints;auto q=stateAt(native.v,native.p,&native.e,native.weight,false);
        Require(data.outerKSource(q.v.row.coordinate,time)==0 && data.outerEpsilonSource(q.v.row.coordinate,time)==0,
            "R3 analytic outer action requires the declared unforced scalar sources");
        const auto closure=OuterClosure(q.k,q.epsilon,q.G,q.v.row.coordinate);
        q.dk=closure.kDiffusion;q.de=closure.epsilonDiffusion;q.sk=closure.kSource;q.se=closure.epsilonSource;
        for(int a=0;a<3;++a)for(int b=0;b<3;++b)
            q.stress[a][b]=-closure.reynoldsStress[a][b]+(a==b?2*q.k/3-q.pi:0);
        auto variant=m_robust;const double r=Distance(q.v.row.coordinate),kap=m_parameters.model.Kappa();
        if(m_variantName=="Original")variant=RANSRobustVariant{};
        else
        {
            variant.gamma=m_regularizationGamma+m_regularizationScale*kap*m_robustUstar*r;
            const double shear=m_robustUstar/(kap*r);
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)
                variant.referenceStrain[a][b]=shear*(m_referenceFlow[a]*m_referenceNormal[b]+m_referenceNormal[a]*m_referenceFlow[b]);
        }
        std::array<OuterAction,6> partial{};RANSTensor zero{};
        const bool strain=closure.strainInvariant>=closure.rotationInvariant,cls=m_parameters.model.Model()==RANSModel::CLSHighRe;
        outerBranches.push_back(strain?'S':'R');
        const long double lengthRatio=static_cast<long double>(q.k)*std::sqrt(static_cast<long double>(q.k))/
            (static_cast<long double>(m_parameters.model.Parameters().yapLength)*q.epsilon*r);
        yapBranches.push_back(cls && lengthRatio>1?'A':'I');
        partial[0]=OuterDerivative(q.k,q.epsilon,q.G,r,1,0,zero,m_parameters.model.Parameters(),cls,strain,variant);
        partial[1]=OuterDerivative(q.k,q.epsilon,q.G,r,0,1,zero,m_parameters.model.Parameters(),cls,strain,variant);
        for(int component=0;component<2;++component)for(int axis=0;axis<2;++axis)
        {
            auto h=zero;h[PhysicalAxis(component)][PhysicalAxis(axis)]=1;
            partial[2+2*component+axis]=OuterDerivative(q.k,q.epsilon,q.G,r,0,0,h,m_parameters.model.Parameters(),cls,strain,variant);
        }
        localKernels+=partial.size();
        auto evaluate=[&](const Variation &d)->Loads
        {
            RANSTensor stress{};double nt=partial[0].viscosity*d.k+partial[1].viscosity*d.epsilon;
            double sk=partial[0].kSource*d.k+partial[1].kSource*d.epsilon;
            double se=partial[0].epsilonSource*d.k+partial[1].epsilonSource*d.epsilon;
            for(int a=0;a<3;++a)for(int b=0;b<3;++b)stress[a][b]=partial[0].stress[a][b]*d.k+partial[1].stress[a][b]*d.epsilon;
            for(int component=0;component<2;++component)for(int axis=0;axis<2;++axis)
            {
                const auto &action=partial[2+2*component+axis];const double factor=d.G[PhysicalAxis(component)][PhysicalAxis(axis)];
                nt+=action.viscosity*factor;sk+=action.kSource*factor;se+=action.epsilonSource*factor;
                for(int a=0;a<3;++a)for(int b=0;b<3;++b)stress[a][b]+=action.stress[a][b]*factor;
            }
            Loads load;
            for(int component=0;component<2;++component)
            {
                const int a=PhysicalAxis(component);RANSVector flux{},change{};
                for(int b=0;b<3;++b)
                {flux[b]=q.stress[a][b]-q.u[a]*q.u[b];change[b]=stress[a][b]-d.u[a]*q.u[b]-q.u[a]*d.u[b];}
                load.Add(offset[component],q.v,0,flux,0,change,q.weight,0,0);
            }
            RANSVector kflux{},dkflux{},eflux{},deflux{};
            for(int a=0;a<3;++a)
            {
                kflux[a]=q.dk*q.gk[a]-q.u[a]*q.k;
                dkflux[a]=q.dk*d.gk[a]+nt/m_parameters.model.Parameters().sigmaK*q.gk[a]-d.u[a]*q.k-q.u[a]*d.k;
                eflux[a]=q.de*q.ge[a]-q.u[a]*q.epsilon;
                deflux[a]=q.de*d.ge[a]+nt/m_parameters.model.Parameters().sigmaEpsilon*q.ge[a]-d.u[a]*q.epsilon-q.u[a]*d.epsilon;
            }
            load.Add(offset[K()],q.v,-q.sk,kflux,-sk,dkflux,q.weight,0,0);
            load.Add(offset[Epsilon()],q.e,-q.se,eflux,-se,deflux,q.weight,0,0);
            load.Add(offset[OuterPressure()],q.p,q.G[0][0]+q.G[2][2],{},d.G[0][0]+d.G[2][2],{},q.weight,0,0);
            return load;
        };
        for(int component=0;component<2;++component)for(std::size_t j=0;j<q.v.row.ids.size();++j)
        {
            Variation d;const int a=PhysicalAxis(component);d.u[a]=q.v.value[j];for(int b=0;b<3;++b)d.G[a][b]=q.v.gradient[b][j];
            scatter(evaluate(d),offset[component]+q.v.row.ids[j]);
        }
        for(std::size_t j=0;j<q.v.row.ids.size();++j)
        {
            Variation d;d.k=q.v.value[j];for(int b=0;b<3;++b)d.gk[b]=q.v.gradient[b][j];
            scatter(evaluate(d),offset[K()]+q.v.row.ids[j]);
        }
        for(std::size_t j=0;j<q.e.row.ids.size();++j)
        {
            Variation d;d.epsilon=q.e.value[j];for(int b=0;b<3;++b)d.ge[b]=q.e.gradient[b][j];
            scatter(evaluate(d),offset[Epsilon()]+q.e.row.ids[j]);
        }
    }

    // External lifting has a fixed derivative. The Upstream matching trace
    // stays an independent unknown in all physical volume equations.
    const auto external=ExternalConstraints(time,data);
    for(auto [column,value]:external)
    {
        (void)value;
        if(column>=offset[Epsilon()] && column<offset[InnerPressure()] && m_traceIds.count(column-offset[Epsilon()]))continue;
        if(column<DifferentialSize())for(auto &row:jac)row[column]=0;
    }
    if(imposeConstraints)
    {
        for(auto [row,value]:constraints)
        {
            (void)value;std::fill(jac[row].begin(),jac[row].end(),0.L);
            if(row<DifferentialSize())jac[row][row]=1;
        }
        const auto face=[&](int axis,int side,const RANSWallGlobalField::Element &element)
            {return axis==2 && !side && std::abs(element.lower[2]-m_parameters.geometry.height)<1e-11;};
        std::set<Coordinate> samples;
        fields[Epsilon()]->Boundary([&](const Coordinate &p){samples.insert(p);return 0.;},face);
        struct TraceRows {Sparse baseline,current;};std::map<Coordinate,TraceRows> target;
        for(const auto &p:samples)
        {
            TraceRows rows;
            auto build=[&](bool baseline)
            {
                const auto &ray=rayLinear(p,baseline);const auto pre=ray.differential->At(m_parameters.geometry.height);
                Sparse result;Append(result,ray.sampled,ray.normalization.lambda*pre.EK);
                Append(result,ray.chi,ray.normalization.lambda*pre.Echi);Append(result,ray.lambda,pre.E);return result;
            };
            rows.baseline=build(true);rows.current=build(false);target.emplace(p,std::move(rows));++traceSamples;
        }
        for(const auto &sample:samples)
        {
            // The actual native constrained segment transform is a fixed
            // linear map. Applying it to each sample unit vector recovers
            // that map, including endpoints and shared corner priorities.
            const auto projection=fields[Epsilon()]->Boundary([&](const Coordinate &p){return p==sample?1.:0.;},face);
            for(auto [id,weight]:projection)
            {
                const auto &direction=m_baselineTraceIds.count(id)?target.at(sample).baseline:target.at(sample).current;
                for(auto [column,value]:direction)if(column<DifferentialSize())jac[offset[Epsilon()]+id][column]-=weight*value;
            }
        }
        // Projected targets do not depend on externally prescribed unknowns.
        for(auto [column,value]:external)
        {
            (void)value;
            if(column>=offset[Epsilon()] && column<offset[InnerPressure()] && m_traceIds.count(column-offset[Epsilon()]))continue;
            if(column<DifferentialSize())for(int id:m_traceIds)jac[offset[Epsilon()]+id][column]=0;
            if(column<DifferentialSize() && constraints.count(column))jac[column][column]=1;
        }
    }
    WallMatrix result(Size(),Vector(Size()));
    for(int row=0;row<Size();++row)for(int column=0;column<DifferentialSize();++column)
    {
        result[row][column]=static_cast<double>(jac[row][column]);
        Require(std::isfinite(result[row][column]),"nonfinite R3 nonpressure action coefficient");
    }
    if(metadata)
    {
        std::ostringstream text;text<<std::setprecision(17)
            <<"{\"policy\":\"R3Actions:v1\",\"nonpressure_columns\":"<<DifferentialSize()
            <<",\"volume_points\":"<<volumePoints<<",\"moving_points\":"<<movingPoints
            <<",\"local_directional_kernels\":"<<localKernels<<",\"profile_differentials\":"<<profileCount
            <<",\"capped_normalizations\":"<<normalizations<<",\"trace_samples\":"<<traceSamples
            <<",\"rays\":"<<cache.size()<<",\"residual_calls\":0,\"domain_calls\":0"
            <<",\"preparatory_constraints_calls\":1,\"preparatory_constraints_seconds\":"<<constraintSeconds
            <<",\"seconds\":"<<std::chrono::duration<double>(Clock::now()-started).count()
            <<",\"outer_branch_pattern\":\""<<outerBranches<<"\",\"outer_Yap_branch_pattern\":\""<<yapBranches<<"\""
            <<",\"invariant_tie\":\"selected strain; zero norm subgradient at origin\""
            <<",\"Yap_knot\":\"selected inactive\",\"cap_knot\":\"selected uncapped\""
            <<",\"moving_quadrature\":\"implicit discrete cutoff, moving panels/weights/field and test basis locations\""
            <<",\"profile_panel_ties\":\"stable ordered native break branch\"}";
        metadata->push_back(text.str());
        std::ostringstream profiles;profiles<<std::setprecision(17)<<"{\"component\":\"R3_profile_differentials\",\"rays\":[";
        bool comma=false;
        for(const auto &[key,value]:cache)
        {
            if(comma)profiles<<',';
            comma=true;
            profiles<<"{\"chart\":["<<key.first[0]<<','<<key.first[1]<<",0],\"baseline\":"<<(key.second?"true":"false")
                <<",\"sampled_K\":"<<value->profile->K()<<",\"chi\":"<<value->profile->Reference().weight
                <<",\"cutoff\":"<<value->profile->Reference().cutoff<<",\"cutoff_K\":"<<value->differential->CutoffK()<<'}';
        }
        profiles<<"]}";metadata->push_back(profiles.str());
    }
    return result;
}

RANSWallEvolutionPatch::WallMatrix RANSWallEvolutionPatch::NonpressureJacobianR3(
    const Vector &state,double time,const RANSWallEvolutionData &data,std::vector<std::string> *metadata) const
{return AssembleNonpressureR3(state,time,data,true,-1,metadata);}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::NonpressureDirectionalActionR3(
    const Vector &state,const Vector &direction,double time,const RANSWallEvolutionData &data,
    bool imposeConstraints,int region) const
{
    Require(direction.size()==state.size(),"R3 nonpressure direction has wrong shape");
    for(double value:direction)Require(std::isfinite(value),"nonfinite R3 nonpressure direction");
    for(int i=DifferentialSize();i<Size();++i)Require(direction[i]==0,"R3 nonpressure direction includes pressure/gauge entries");
    const auto matrix=AssembleNonpressureR3(state,time,data,imposeConstraints,region,nullptr);
    Vector result(Size());
    for(int i=0;i<Size();++i)
    {
        long double value=0;for(int j=0;j<DifferentialSize();++j)value+=static_cast<long double>(matrix[i][j])*direction[j];
        result[i]=static_cast<double>(value);
        Require(std::isfinite(result[i]),"nonfinite R3 nonpressure directional action");
    }
    return result;
}
} // namespace Nektar
