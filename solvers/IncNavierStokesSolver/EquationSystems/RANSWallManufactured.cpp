// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallManufactured.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallManufacturedJet.h>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <map>
#include <sstream>
#include <tuple>
namespace Nektar
{
using namespace SolverUtils;
using RANSWallAlgebra::Require;
namespace
{
using namespace CT6;
struct Geometry
{
    Jet a,J,gamma;
    Vec tangent{},normal{};
    Geometry(const Jet &q,const Jet &r,double amplitude,double length)
    {
        const Jet h=amplitude*Exp(-q*q/(length*length));
        const Jet hp=(-2/(length*length))*q*h;
        const Jet hpp=(4*q*q/(length*length*length*length)-2/(length*length))*h;
        a=Sqrt(1+hp*hp);gamma=-hpp/(a*a*a);J=1+r*gamma;
        tangent={1/a,0,hp/a};normal={-hp/a,0,1/a};
    }
    Vec Grad(const Jet &v) const
    {Vec g{};for(int i=0;i<3;++i)g[i]=tangent[i]*v.D(0)/(a*J)+normal[i]*v.D(1);return g;}
    Mat Gradient(const Vec &v) const{Mat g{};for(int i=0;i<3;++i)g[i]=Grad(v[i]);return g;}
    Jet Div(const Vec &v) const{Jet s;for(int i=0;i<3;++i)s=s+Grad(v[i])[i];return s;}
};
struct Closure {Mat R{},Q{};Jet nut,P,ep,ed,rad;};
// Independent original laws in the planar Cartesian sector. The analytic
// largest eigenvalue of the xz strain block replaces native Jacobi iteration.
Closure Outer(RANSModel model,const Jet &k,const Jet &e,const Mat &G)
{
    Mat S{},O{};Jet tr,ss,oo;const Jet tau=k/e;
    for(int i=0;i<3;++i)tr=tr+G[i][i];
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {S[i][j]=(G[i][j]+G[j][i])*.5-(i==j?tr/3:Jet{});O[i][j]=(G[i][j]-G[j][i])*.5;ss=ss+S[i][j]*S[i][j];oo=oo+O[i][j]*O[i][j];}
    const Jet eta=tau*Sqrt(2*ss),xi=tau*Sqrt(2*oo);Jet cmu,cq;Closure law;
    if(model==RANSModel::A93)cmu=(2./3)/(1.25+eta+.9*xi);
    else
    {
        const Jet lambda=(S[0][0]+S[2][2])*.5+Sqrt((S[0][0]-S[2][2])*(S[0][0]-S[2][2])*.25+S[0][2]*S[0][2]);
        cmu=1/((model==RANSModel::A95?6.5:4.)+3*lambda*tau*Sqrt((ss+oo)/ss));
        if(model==RANSModel::A95){law.rad=1-9*cmu*cmu*tau*tau*ss;cq=Sqrt(law.rad)/(1+6*tau*tau*Sqrt(ss*oo));}
    }
    law.nut=cmu*k*k/e;const auto s2=Product(S,S),o2=Product(O,O),os=Product(O,S),so=Product(S,O);
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)
    {
        if(model==RANSModel::A93)law.Q[i][j]=k*tau*tau/(1000+eta*eta*eta)*(3*(s2[i][j]-(i==j?ss/3:Jet{}))-19*(o2[i][j]+(i==j?oo/3:Jet{}))+15*(os[i][j]-so[i][j]));
        else if(model==RANSModel::A95)law.Q[i][j]=2*cq*k*tau*tau*(os[i][j]-so[i][j]);
        law.R[i][j]=(i==j?2*k/3:Jet{})-2*law.nut*S[i][j]+law.Q[i][j];law.P=law.P-law.R[i][j]*G[i][j];
    }
    if(model==RANSModel::E94E95){const Jet c1=(eta/(eta+5)).Value()>.43?eta/(eta+5):Jet(.43);law.ep=c1*Sqrt(2*ss)*e;law.ed=1.9*e*e/k;}
    else{law.ep=1.44*e*law.P/k;law.ed=1.92*e*e/k;}
    return law;
}
Vec Cartesian(const Jet &ut,const Jet &w,const Geometry &g)
{Vec u{};for(int i=0;i<3;++i)u[i]=ut*g.tangent[i]+w*g.normal[i];return u;}
Jet Normal(const Mat &a,const Vec &left,const Vec &right)
{Jet r;for(int i=0;i<3;++i)for(int j=0;j<3;++j)r=r+left[i]*a[i][j]*right[j];return r;}
}
struct RANSWallManufactured::Impl
{
    const RANSWallEvolutionPatch &patch;
    RANSNeutralReference model;
    double H,nu,hr,K0,amplitude,length,kap,qm;
    mutable std::map<std::tuple<double,double,double,bool>,RANSWallManufacturedSample> cache;
    struct Ray {Jet c,l,jump,cutoff;};
    mutable std::map<std::pair<double,double>,Ray> rays;
    explicit Impl(const RANSWallEvolutionPatch &p):patch(p),model(p.Parameters().model),H(p.Parameters().geometry.height),nu(p.Parameters().geometry.molecularNu),hr(p.Parameters().geometry.roughness),K0(model.Q()*.16),amplitude(p.Mapping().kind==RANSWallNormalMap::Kind::Flat?0:p.Mapping().amplitude),length(p.Mapping().lengthX),kap(model.Kappa()),qm(model.Q())
    {
        Require(p.CurvedManufactured() && p.OpenBoundary() && model.NewModel() && model.KappaPolicy()==RANSKappaPolicy::Natural,"CT6 continuous provider requires its named original Natural manufactured capability");
        Require(nu==1.5e-5 && hr==.05 && H==.5 && length==1. && (amplitude==0. || amplitude==.005),"CT6 manufactured physical inputs differ from bounded contract");
    }
    Jet K(const Jet &q,const Jet &t) const{return K0*(1+.02*q+.01*q*q)*(1+.05*t+.1*t*t);}
    Jet A(const Jet &k) const{return Pow(k,1.5)/(kap*std::pow(qm,1.5));}
    Jet B(const Jet &t) const{return .2*(1+.2*t+.5*t*t);}
    Jet Nut(const Jet &q,const Jet &r,const Jet &t) const
    {
        const Jet k=K(q,t),root=Sqrt(k),x=hr*root/(70*nu),m1=.5-.4*Pow(x,.7),m2=1-.79*Pow(x,-.28);
        const Jet m=m1.Value()>m2.Value()?m1:m2,b=10.7*(1-Exp(m*Log(x))),zv=nu*b/root;
        Require(zv.Value()<0,"CT6 independent rough primitive requires negative viscosity intercept");
        return (kap/std::sqrt(qm))*root*(r-zv);
    }
    Vec Base(const Jet &q,const Jet &r,const Jet &t,const Geometry &g) const
    {return Cartesian(B(t)*(1+.02*q)*r,-B(t)*.02*r*r/(2*g.a*g.J),g);}
    Mat Stress(const Closure &law,const Jet &k,const Jet &pi) const
    {Mat T{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)T[i][j]=-law.R[i][j]+(i==j?2*k/3-pi:Jet{});return T;}
    Ray AtH(double qv,double tv) const
    {
        const auto key=std::make_pair(qv,tv);auto found=rays.find(key);if(found!=rays.end())return found->second;
        const auto q=Jet::Variable(qv,0),r=Jet::Variable(H,1),t=Jet::Variable(tv,2);Geometry g(q,r,amplitude,length);
        const auto u=Base(q,r,t,g);const Mat G=g.Gradient(u);const Jet k=K(q,t),e=A(k)/H;const auto law=Outer(model.Model(),k,e,G);
        Mat F{};for(int i=0;i<3;++i)for(int j=0;j<3;++j)F[i][j]=-law.R[i][j]+(i==j?2*k/3:Jet{});
        const Jet D=nu+Nut(q,r,t);Ray v;
        // Interface coefficients are functions of q,t only. Their r Taylor
        // coefficients MUST be removed before extending along a ray.
        const auto freeze=[](Jet a){const auto &p=Jet::Table().p;for(int i=0;i<35;++i)if(p[i][1])a.c[i]=0;return a;};
        const Jet w=-B(t)*.02*r*r/(2*g.a*g.J),ut=B(t)*(1+.02*q)*r;
        v.c=freeze(Normal(F,g.tangent,g.normal)/D-(B(t)*(1+.02*q)+w.D(0)/(g.a*g.J)-g.gamma*ut/g.J));
        v.l=freeze(law.nut/D*(.03*K0/H));
        // w_r is common, so the inner normal viscous stress is unaffected by c.
        v.jump=freeze(Normal(F,g.normal,g.normal)-2*D*w.D(1));
        const Jet aa=(kap/std::sqrt(qm))*Sqrt(k),d=nu+Nut(q,Jet(0),t),V=d+aa*H,lg=Log(V/d),t0=k/qm,Ae=A(k);
        const Jet production=t0*t0/aa*(lg+nu/V-nu/d)+g.gamma*t0*t0/(aa*aa)*(aa*H-(d+nu)*lg+d*nu*(1/d-1/V));
        Jet logc(std::log(Cutoff(k.Value(),g.gamma.Value())));
        for(int i=0;i<5;++i){const Jet cut=Exp(logc);logc=logc-(Ae*(1+std::log(H)-logc+g.gamma*(H-cut/2))-production)/(-Ae*(1+g.gamma*cut/2));}
        v.cutoff=freeze(Exp(logc));
        if(rays.size()>4096)rays.clear();
        rays.emplace(key,v);return v;
    }
    double Cutoff(double k,double gamma) const
    {
        const double root=std::sqrt(k),x=hr*root/(70*nu),m=std::max(.5-.4*std::pow(x,.7),1-.79*std::pow(x,-.28));
        const double b=10.7*(1-std::pow(x,m)),a=kap*root/std::sqrt(qm),d=nu-a*nu*b/root,V=d+a*H,t0=k/qm,Ae=std::pow(k/qm,1.5)/kap;
        const double logarithm=std::log(V/d);
        const double p0=t0*t0/a*(logarithm+nu/V-nu/d);
        const double p1=t0*t0/(a*a)*(a*H-(d+nu)*logarithm+d*nu*(1/d-1/V));
        const double production=p0+gamma*p1;
        Require(production>=Ae*(1+gamma*H/2),"CT6 independent metric budget inadmissible");
        double lo=std::log(H)-64,hi=std::log(H);
        for(int i=0;i<120;++i){const double mid=(lo+hi)/2,c=std::exp(mid);if(Ae*(1+std::log(H/c)+gamma*(H-c/2))>production)lo=mid;else hi=mid;}
        return std::exp((lo+hi)/2);
    }
    RANSWallManufacturedSample Evaluate(const RANSWallCoordinate &p,double time,bool outer) const
    {
        const auto key=std::make_tuple(p[0],p[2],time,outer);auto old=cache.find(key);if(old!=cache.end())return old->second;
        const auto q=Jet::Variable(p[0],0),r=Jet::Variable(p[2],1),t=Jet::Variable(time,2);Geometry g(q,r,amplitude,length);const auto ray=AtH(p[0],time);
        const Jet phi=r*r*(r-H)*(r-H)/(2*H*H),phiR=(2*r*(r-H)*(r-H)+2*r*r*(r-H))/(2*H*H);
        const Jet ut=B(t)*(1+.02*q)*r+(outer?Jet{}:ray.c*phiR);
        const Jet w=-(B(t)*.02*r*r/2+(outer?Jet{}:ray.c.D(0)*phi))/(g.a*g.J);
        const auto u=Cartesian(ut,w,g);const auto G=g.Gradient(u);
        const Jet k=K(q,t)+(outer?(.03*K0/H)*(r-H):ray.l*(r*r-H*H)/(2*H));
        const Jet eps=A(K(q,t))/(outer || r.Value()>ray.cutoff.Value()?r:ray.cutoff);
        const Jet pi=.03*(q+.2*r)+.02*t+.01*t*t+(outer?ray.jump:Jet{});
        Closure law;Mat T{};Jet D,de;
        if(outer){law=Outer(model.Model(),k,eps,G);D=law.nut;de=D/model.Parameters().sigmaEpsilon;T=Stress(law,k,pi);}
        else
        {
            law.nut=Nut(q,r,t);D=nu+law.nut;
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            {law.R[i][j]=(i==j?2*k/3:Jet{})-law.nut*(G[i][j]+G[j][i]);law.P=law.P-law.R[i][j]*G[i][j];T[i][j]=D*(G[i][j]+G[j][i])-(i==j?pi:Jet{});}
        }
        Vec body{};for(int i=0;i<3;++i){Vec flux;for(int j=0;j<3;++j)flux[j]=u[i]*u[j]-T[i][j];body[i]=u[i].D(2)+g.Div(flux);}
        const auto kg=g.Grad(k),eg=g.Grad(eps);Vec kf{},ef{};
        for(int i=0;i<3;++i){kf[i]=u[i]*k-D*kg[i];ef[i]=u[i]*eps-de*eg[i];}
        RANSWallManufacturedSample v;const double sampled=K(q,t).Value();v.K=sampled;v.cutoff=ray.cutoff.Value();
        const double sink=eps.Value();
        v.kSource=(k.D(2)+g.Div(kf)-law.P).Value()+(outer?eps.Value():sink);
        v.epsilonSource=outer?(eps.D(2)+g.Div(ef)-law.ep+law.ed).Value():0;
        v.kTime=k.D(2).Value();v.epsilonTime=eps.D(2).Value();
        v.value={u[0].Value(),u[2].Value(),k.Value(),outer?eps.Value():sink,pi.Value()};
        const std::array<Jet,5> fields{u[0],u[2],k,eps,pi};
        for(int b=0;b<5;++b){const auto grad=g.Grad(fields[b]);for(int j=0;j<3;++j)v.gradient[b][j]=grad[j].Value();}
        for(int i=0;i<3;++i)
        {
            v.velocityTime[i]=u[i].D(2).Value();v.body[i]=body[i].Value();v.normal[i]=g.normal[i].Value();v.tangent[i]=g.tangent[i].Value();
            for(int j=0;j<3;++j){v.velocityGradient[i][j]=G[i][j].Value();v.reynolds[i][j]=law.R[i][j].Value();v.quadratic[i][j]=law.Q[i][j].Value();v.stress[i][j]=T[i][j].Value();}
        }
        v.kDiffusion=D.Value();v.epsilonDiffusion=de.Value();v.production=law.P.Value();v.epsilonProduction=law.ep.Value();v.epsilonDestruction=law.ed.Value();v.radicand=law.rad.Value();v.J=g.J.Value();v.a=g.a.Value();v.gamma=g.gamma.Value();
        if(cache.size()>30000)cache.clear();
        cache.emplace(key,v);return v;
    }
    RANSVector Rough(double q,double time) const
    {
        const auto v=Evaluate({q,0,H},time,true);const double ut=v.value[0]*v.tangent[0]+v.value[1]*v.tangent[2];
        double derivative=0,pressure=0,body=0;
        for(int a=0;a<3;++a)
        {
            double d=0;for(int i=0;i<3;++i)d+=v.tangent[i]*v.velocityGradient[i][a];
            derivative+=d*v.tangent[a];pressure+=v.tangent[a]*(v.gradient[4][a]-2*v.gradient[2][a]/3);body+=v.tangent[a]*v.body[a];
        }
        const double w=v.value[0]*v.normal[0]+v.value[1]*v.normal[2];derivative-=v.gamma/v.J*w;
        const double force=2*ut*derivative+pressure-body;RANSVector f{};for(int i=0;i<3;++i)f[i]=v.tangent[i]*force;return f;
    }
};
RANSWallManufactured::RANSWallManufactured(const RANSWallEvolutionPatch &p):impl(std::make_shared<Impl>(p)){}
RANSWallManufacturedSample RANSWallManufactured::At(const RANSWallCoordinate &p,double t,bool outer) const{return impl->Evaluate(p,t,outer);}
RANSWallEvolutionData RANSWallManufactured::Data() const
{
    RANSWallEvolutionData d(2);auto p=impl;d.openBoundaryIdentity=p->patch.OpenBoundaryIdentity();
    for(int a=0;a<2;++a)
    {
        const int axis=a?2:0;
        d.innerBody[a]=[p,axis](const RANSWallCoordinate &q,double t){const auto v=p->Evaluate(q,t,false);return v.body[axis]-(q[2]<=p->hr?p->Rough(q[0],t)[axis]:0);};
        d.outerBody[a]=[p,axis](const RANSWallCoordinate &q,double t){return p->Evaluate(q,t,true).body[axis];};
        d.boundary[a]=[p,a](const RANSWallCoordinate &q,double t){return p->Evaluate(q,t,q[2]>p->H).value[a];};
        for(int region=0;region<2;++region)d.outletTraction[region][a]=[p,axis,region](const RANSWallCoordinate &q,double t){const auto v=p->Evaluate(q,t,region);const auto area=p->patch.Mapping().Cofactor(q,0,1);const double norm=std::hypot(area[0],area[2]);double traction=0;for(int j=0;j<3;++j)traction+=v.stress[axis][j]*area[j]/norm;return traction;};
    }
    d.boundary[2]=[p](const RANSWallCoordinate &q,double t){return p->Evaluate(q,t,q[2]>p->H).value[2];};
    d.boundary[3]=[p](const RANSWallCoordinate &q,double t){return p->Evaluate(q,t,true).value[3];};
    d.innerKSource=[p](const RANSWallCoordinate &q,double t){return p->Evaluate(q,t,false).kSource;};
    d.outerKSource=[p](const RANSWallCoordinate &q,double t){return p->Evaluate(q,t,true).kSource;};
    d.outerEpsilonSource=[p](const RANSWallCoordinate &q,double t){return p->Evaluate(q,t,true).epsilonSource;};
    for(int region=0;region<2;++region)d.outletKDiffusion[region]=[p,region](const RANSWallCoordinate &q,double t){const auto v=p->Evaluate(q,t,region);const auto area=p->patch.Mapping().Cofactor(q,0,1);return v.kDiffusion*(v.gradient[2][0]*area[0]+v.gradient[2][2]*area[2])/std::hypot(area[0],area[2]);};
    d.outletEpsilonDiffusion=[p](const RANSWallCoordinate &q,double t){const auto v=p->Evaluate(q,t,true);const auto area=p->patch.Mapping().Cofactor(q,0,1);return v.epsilonDiffusion*(v.gradient[3][0]*area[0]+v.gradient[3][2]*area[2])/std::hypot(area[0],area[2]);};
    return d;
}
std::vector<RANSWallEvolutionFunction> RANSWallManufactured::InitialFunctions() const
{
    auto p=impl;std::vector<RANSWallEvolutionFunction> f;
    for(int b=0;b<p->patch.Gauge();++b)f.push_back([p,b](const RANSWallCoordinate &q,double t){const bool outer=b==p->patch.OuterPressure() || b==p->patch.Epsilon() || (b<=p->patch.K() && q[2]>p->H);const auto v=p->Evaluate(q,t,outer);return b<p->patch.InnerPressure()?v.value[b]:v.value[4];});
    return f;
}
bool RANSUsesAcceptedManufactured(const LibUtilities::SessionReaderSharedPtr &s)
{std::string v;s->LoadSolverInfo("RANSWallManufactured",v,"Off");Require(v=="Off" || v=="CT6ContinuousV1","unknown Accepted manufactured provider");return v=="CT6ContinuousV1";}
std::string RANSAcceptedManufacturedIdentity(const LibUtilities::SessionReaderSharedPtr &s)
{
    Require(RANSUsesAcceptedManufactured(s),"manufactured boundary requires CT6ContinuousV1");
    for(const std::string f:{"InitialConditionsInner","InitialConditionsOuter","BoundaryValuesInner","BoundaryValuesOuter","BodyForceInner","BodyForceOuter","ScalarSourceInner","ScalarSourceOuter","OutletTotalTractionInner","OutletTotalTractionOuter","OutletDiffusiveFluxInner","OutletDiffusiveFluxOuter"})Require(!s->DefinesFunction(f),"CT6 provider forbids competing physical production functions");
    std::ostringstream o;o<<std::setprecision(17)<<"AcceptedManufacturedOpenV1:CT6ContinuousV1:Cartesian-stress:chart-functions:Pi-absolute:no-mean-gauge:no-state-restoration";
    const std::set<std::string> work{"NUMSTEPS","TIMESTEP","TIME","IO_CHECKSTEPS","IO_INFOSTEPS","IO_TIMER_LEVEL","RANSWALLNORMALORDER","RANSWALLTANGENTIALORDER","RANSWALLMASSORDER","RANSWALLVOLUMEORDER","RANSWALLINTERFACEORDER"};
    for(const auto &[k,v]:s->GetParameters())if(!work.count(k))o<<'|'<<k<<'='<<v;
    return o.str();
}
} // namespace Nektar
