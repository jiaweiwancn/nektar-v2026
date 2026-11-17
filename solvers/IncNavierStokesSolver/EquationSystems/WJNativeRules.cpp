// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/WJNativeRules.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>
namespace Nektar::WJNative
{
namespace
{
void Require(bool v,const char*s){if(!v)throw std::runtime_error(s);}
std::pair<Vec,Vec>Gauss(int n)
{
    Vec x(n),w(n);double pi=std::acos(-1.);
    for(int i=0;i<n;++i){double z=std::cos(pi*(i+.75)/(n+.5)),d=0;
        for(int iteration=0;iteration<48;++iteration){double a=1,b=z;for(int k=2;k<=n;++k){double c=((2*k-1)*z*b-(k-1)*a)/k;a=b;b=c;}d=n*(z*b-a)/(z*z-1);double delta=b/d;z-=delta;if(std::abs(delta)<2e-16)break;}x[i]=z;w[i]=2/((1-z*z)*d*d);}
    return{x,w};
}
struct Factors
{
    Vec nodes,barycentric;
    Mat value,derivative;
    Factors(const LocalRegions::ExpansionSharedPtr&e,int direction)
    {
        auto basis=e->GetBasis(direction);auto z=basis->GetZ();int nq=z.size(),modes=e->GetBasisNumModes(direction);
        nodes=Vec(nq);barycentric=Vec::Ones(nq);value=Mat(nq,modes);derivative=Mat(nq,modes);
        const auto b=basis->GetBdata(),db=basis->GetDbdata();
        for(int i=0;i<nq;++i){nodes[i]=z[i];for(int j=0;j<modes;++j){value(i,j)=b[j*nq+i];derivative(i,j)=db[j*nq+i];}}
        for(int i=0;i<nq;++i)for(int j=0;j<nq;++j)if(i!=j)barycentric[i]/=nodes[i]-nodes[j];
    }
    std::pair<Vec,Vec>At(double x)const
    {
        for(int i=0;i<nodes.size();++i)if(std::abs(x-nodes[i])<4e-15)return{value.row(i).transpose(),derivative.row(i).transpose()};
        Vec l(nodes.size());for(int i=0;i<nodes.size();++i)l[i]=barycentric[i]/(x-nodes[i]);l/=l.sum();return{value.transpose()*l,derivative.transpose()*l};
    }
};
}
RuleFactory::RuleFactory(const Flow&flow,int nx,int nz):m_flow(flow),m_nx(nx),m_nz(nz)
{
    Require(flow.Mapping().curve==0&&nx>=4&&nz>=4,"split native rules initially support flat quads and Gauss order>=4");
    for(const auto&e:flow.Elements())
    {
        auto x=e.velocity->GetBasis(0)->GetZ(),z=e.velocity->GetBasis(1)->GetZ();
        for(int j=0;j<int(z.size());++j)for(int i=0;i<int(x.size());++i){int q=i+x.size()*j;double px=(e.xl+e.xr)/2+(e.xr-e.xl)*x[i]/2,pz=flow.Mapping().d0+(e.rl+e.rr)/2+(e.rr-e.rl)*z[j]/2;Require(std::abs(e.x[q]-px)<1e-10&&std::abs(e.z[q]-pz)<1e-10,"native rule geometry is not an affine flat rectangle");}
        m_static.push_back(std::make_shared<NativeRule>(Volume(e,{})));
    }
}
NativeRule RuleFactory::Evaluate(const Element&e,const Vec&xi,const Vec&eta,const Vec&weights)const
{
    Require(xi.size()==eta.size()&&xi.size()==weights.size(),"native rule coordinate/weight counts differ");
    NativeRule r;r.geometry=e.geometry;r.u=e.u;r.p=e.p;r.us=e.us;r.ps=e.ps;r.w=weights;r.x=Vec(xi.size());r.z=Vec(xi.size());
    r.B=Mat(xi.size(),16);r.Dx=r.B;r.Dz=r.B;r.P=Mat(xi.size(),9);r.Px=r.P;r.Pz=r.P;
    Factors ux(e.velocity,0),uz(e.velocity,1),px(e.pressure,0),pz(e.pressure,1);
    double sx=2/(e.xr-e.xl),sz=2/(e.rr-e.rl);
    for(int q=0;q<xi.size();++q)
    {
        auto x=ux.At(xi[q]),z=uz.At(eta[q]),p=px.At(xi[q]),v=pz.At(eta[q]);
        r.x[q]=(e.xl+e.xr)/2+(e.xr-e.xl)*xi[q]/2;r.z[q]=m_flow.Mapping().d0+(e.rl+e.rr)/2+(e.rr-e.rl)*eta[q]/2;
        for(int j=0;j<16;++j){int a=j%4,b=j/4;r.B(q,j)=x.first[a]*z.first[b]*e.us[j];r.Dx(q,j)=x.second[a]*z.first[b]*sx*e.us[j];r.Dz(q,j)=x.first[a]*z.second[b]*sz*e.us[j];}
        for(int j=0;j<9;++j){int a=j%3,b=j/3;r.P(q,j)=p.first[a]*v.first[b]*e.ps[j];r.Px(q,j)=p.second[a]*v.first[b]*sx*e.ps[j];r.Pz(q,j)=p.first[a]*v.second[b]*sz*e.ps[j];}
    }
    return r;
}
NativeRule RuleFactory::Volume(const Element&e,const std::function<double(const Element&,double)>&cutoff,
    const std::vector<double>&xBreaks)const
{
    auto x=Gauss(m_nx),z=Gauss(m_nz);std::vector<double>xs,ys,ws;const auto geo=m_flow.Mapping();
    const double xRoundoff=64*std::numeric_limits<double>::epsilon()*std::max({std::abs(e.xl),std::abs(e.xr),e.xr-e.xl});
    std::vector<double>xPanels{e.xl,e.xr};
    for(double point:xBreaks){
        Require(std::isfinite(point)&&point>=e.xl-xRoundoff&&point<=e.xr+xRoundoff,"cutoff x partition left its original owner");
        if(point>e.xl+xRoundoff&&point<e.xr-xRoundoff)xPanels.push_back(point);
    }
    std::sort(xPanels.begin(),xPanels.end());
    if(xPanels.size()>2){
        std::vector<double>distinct;
        for(double point:xPanels)if(distinct.empty()||point-distinct.back()>xRoundoff)distinct.push_back(point);
        xPanels=std::move(distinct);
    }
    for(int panel=0;panel+1<int(xPanels.size());++panel)for(int i=0;i<m_nx;++i)
    {
        const double left=xPanels[panel],right=xPanels[panel+1],physicalX=(left+right)/2+(right-left)*x.first[i]/2;
        // Preserve the original unsplit rule exactly when there is no x event.
        const double xi=xPanels.size()==2?x.first[i]:2*(physicalX-e.xl)/(e.xr-e.xl)-1;
        std::vector<double>cuts{e.rl,e.rr};
        // Native physical coordinates at an existing material edge can differ
        // from its supplied model radius by a few arithmetic ulps. Do not
        // create a spurious panel across that same edge (e.g. [H,H+2e-17]).
        // No state/cutoff is clamped; only numerically coincident geometric
        // integration boundaries are identified.
        const double roundoff=64*std::numeric_limits<double>::epsilon()*
            std::max({geo.H,std::abs(e.rl),std::abs(e.rr)});
        for(double a:{.4*geo.H,.75*geo.H,geo.H,geo.h})if(a>e.rl+roundoff&&a<e.rr-roundoff)cuts.push_back(a);
        if(cutoff){double c=cutoff(e,physicalX);Require(std::isfinite(c)&&c>0&&c<.4*geo.H,"current cutoff left certified below-collar branch");if(c>e.rl+roundoff&&c<e.rr-roundoff)cuts.push_back(c);}
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
        for(int s=0;s+1<int(cuts.size());++s)for(int j=0;j<m_nz;++j){double physicalR=(cuts[s]+cuts[s+1])/2+(cuts[s+1]-cuts[s])*z.first[j]/2;xs.push_back(xi);ys.push_back(2*(physicalR-e.rl)/(e.rr-e.rl)-1);ws.push_back(x.second[i]*z.second[j]*(right-left)*(cuts[s+1]-cuts[s])/4);}
    }
    Vec xi(xs.size()),eta(ys.size()),weight(ws.size());for(int i=0;i<int(xs.size());++i){xi[i]=xs[i];eta[i]=ys[i];weight[i]=ws[i];}return Evaluate(e,xi,eta,weight);
}
Rules RuleFactory::Cutoff(const std::function<double(const Element&,double)>&cutoff,const CutoffPartition*partition)const
{
    Require(bool(cutoff),"missing current cutoff function");Rules result=m_static;auto q=Gauss(m_nx);
    if(partition)Require(bool(partition->xBreaks)&&std::isfinite(partition->radiusRange[0])&&std::isfinite(partition->radiusRange[1])&&
        partition->radiusRange[0]>0&&partition->radiusRange[1]>=partition->radiusRange[0]&&partition->radiusRange[1]<.4*m_flow.Mapping().H,
        "invalid certified cutoff partition range");
    for(int i=0;i<int(m_flow.Elements().size());++i)
    {
        const auto&e=m_flow.Elements()[i];if(e.rl>=.4*m_flow.Mapping().H||e.rr<=0)continue;bool intersects=false;
        if(partition){
            // The whole certified range protects narrow intersections that
            // are absent at every original x Gauss node.
            if(e.rr<partition->radiusRange[0]||e.rl>partition->radiusRange[1])continue;
            result[i]=std::make_shared<NativeRule>(Volume(e,cutoff,partition->xBreaks(e)));continue;
        }
        for(int j=0;j<m_nx;++j){double x=(e.xl+e.xr)/2+(e.xr-e.xl)*q.first[j]/2,c=cutoff(e,x);Require(std::isfinite(c)&&c>0&&c<.4*m_flow.Mapping().H,"cutoff outside certified below-collar branch");intersects|=c>e.rl&&c<e.rr;}
        if(intersects)result[i]=std::make_shared<NativeRule>(Volume(e,cutoff));
    }
    return result;
}
NativeRule RuleFactory::Face(const Element&e,int side,int order,
    const std::function<double(const Element&,double)>&cutoff)const
{
    Require(side>=0&&side<4&&order>=4,"invalid native face rule");auto q=Gauss(order);const auto geo=m_flow.Mapping();
    std::vector<double>cuts{e.rl,e.rr};
    const double roundoff=64*std::numeric_limits<double>::epsilon()*std::max({geo.H,std::abs(e.rl),std::abs(e.rr)});
    if(side<2){
        for(double a:{.4*geo.H,.75*geo.H,geo.H,geo.h})if(a>e.rl+roundoff&&a<e.rr-roundoff)cuts.push_back(a);
        if(cutoff){double x=side==0?e.xl:e.xr,c=cutoff(e,x);Require(std::isfinite(c)&&c>0&&c<.4*geo.H,"face cutoff left certified below-collar branch");if(c>e.rl+roundoff&&c<e.rr-roundoff)cuts.push_back(c);}
        std::sort(cuts.begin(),cuts.end());cuts.erase(std::unique(cuts.begin(),cuts.end()),cuts.end());
    }else if(cutoff){double radius=side==2?e.rl:e.rr;Require(radius<=roundoff||radius>=.4*geo.H-roundoff,"moving-cutoff intersection on an interior-height horizontal face is not supported");}
    const int panels=side<2?int(cuts.size())-1:1;Vec xi(order*panels),eta(xi.size()),weight(xi.size());
    for(int panel=0;panel<panels;++panel)for(int i=0;i<order;++i){int j=panel*order+i;
        if(side<2){double r=(cuts[panel]+cuts[panel+1])/2+(cuts[panel+1]-cuts[panel])*q.first[i]/2;xi[j]=side==0?-1.:1.;eta[j]=2*(r-e.rl)/(e.rr-e.rl)-1;weight[j]=q.second[i]*(cuts[panel+1]-cuts[panel])/2;}
        else{xi[j]=q.first[i];eta[j]=side==2?-1.:1.;weight[j]=q.second[i]*(e.xr-e.xl)/2;}
    }
    NativeRule result=Evaluate(e,xi,eta,weight);result.normal=side==0?Vector{-1,0,0}:side==1?Vector{1,0,0}:side==2?Vector{0,0,-1}:Vector{0,0,1};return result;
}
NativeRule RuleFactory::Points(const Element&e,const Vec&x,const Vec&z)const
{
    Require(x.size()==z.size(),"native point coordinate counts differ");
    Vec xi=2*(x.array()-(e.xl+e.xr)/2)/(e.xr-e.xl);
    Vec eta=2*(z.array()-m_flow.Mapping().d0-(e.rl+e.rr)/2)/(e.rr-e.rl);
    Require((xi.array().abs()<=1+1e-11).all()&&(eta.array().abs()<=1+1e-11).all(),"native point lies outside its specified owner");
    return Evaluate(e,xi,eta,Vec::Zero(x.size()));
}
NativeRule RuleFactory::HorizontalTrace(const Element&e,const Vec&x,int side)const
{
    Require(side==2||side==3,"native horizontal trace requires lower/upper side");
    // Use the same endpoint-relative parameter as the direct edge polynomial.
    // The selected original owner and edge remain explicit; no state, cutoff
    // or arbitrary physical sampling point is snapped to this edge.
    Vec xi=2*(x.array()-e.xl)/(e.xr-e.xl)-1;
    Require((xi.array().abs()<=1+1e-11).all(),"native horizontal trace lies outside its specified owner");
    return Evaluate(e,xi,Vec::Constant(x.size(),side==2?-1.:1.),Vec::Zero(x.size()));
}
}
