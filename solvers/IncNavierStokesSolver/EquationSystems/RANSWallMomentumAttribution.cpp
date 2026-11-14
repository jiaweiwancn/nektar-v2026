// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <set>
#include <numeric>

namespace Nektar
{
using namespace RANSWallAlgebra;
using namespace SolverUtils;
namespace J=RANSWallSolveJson;
namespace
{
std::string AttributionRow(const RANSWallGlobalPoint &q)
{
    return J::Object({{"ids",J::Numbers(q.ids)},{"phi",J::Numbers(q.phi)},
        {"derivative",J::Rows(q.derivative)},{"geometry",J::Number(q.geometry)},
        {"element",J::Number(q.element)},{"coordinate",J::Numbers(q.coordinate)}});
}
}
std::string RANSWallEvolutionPatch::MomentumAttributionPacket(const Vector &state,double time,
    const RANSWallEvolutionData &data,bool exportPoints) const
{
    Require(SerialSparse() && !Mapped() && m_dimension==2 && DisplacedAccepted() && RichardsNorrisTop() &&
        m_parameters.model.Model()==RANSModel::A95 && m_eta==0.,"CT10 attribution requires the unanchored displaced A95 flat/RN sparse route");
    Require(!exportPoints || Size()<=12000,"CT10 full native point export is limited to12000 slots; aggregate sparse loads remain available");
    const int nm=offset[K()];const auto s=ExternalLift(state,time,data);Refresh(s,time,&data);
    double us,level;m_sparseSession->LoadParameter("RANSWallFrictionVelocity",us);m_sparseSession->LoadParameter("RANSOpenPressureLevel",level,0.);
    const auto &model=m_parameters.model;const auto &g=m_parameters.geometry;
    const double H=g.height,d=AerodynamicDisplacement(),K0=model.Q()*us*us,A=us*us*us/model.Kappa();
    const auto ref=Profile(K0,{lower[0],0.,0.});
    RANSTensor referenceG{};referenceG[0][2]=us/(model.Kappa()*(H-d));
    const auto referenceLaw=OuterClosure(K0,A/(H-d),referenceG,{lower[0],0.,H});
    const double piOuter=level-referenceLaw.quadraticStress[2][2];
    using Terms=std::map<std::string,Vector>;
    using Local=std::map<std::string,std::map<int,long double>>;
    Terms actual,analytic;std::map<std::string,Local> locals[2];
    std::vector<std::string> points,samplers,interfaces,facePoints;
    std::set<int> regionIds[2];
    for(int region=0;region<2;++region)for(const auto &e:fields[0]->Elements())
        if((e.lower[2]>=H)==bool(region))
        {
            auto mid=e.lower;for(int axis:{0,2})mid[axis]=(e.lower[axis]+e.upper[axis])/2;
            const auto q=GeometryPoint(0,e.geometry,mid);regionIds[region].insert(q.ids.begin(),q.ids.end());
        }
    std::vector<std::string> shared;
    for(int id:regionIds[0])if(regionIds[1].count(id))shared.push_back(J::Number(id));
    auto add=[&](int version,const std::string &key,const std::string &element,int block,
        const RANSWallGlobalPoint &q,double value,const RANSVector &flux,double weight)
    {
        auto &terms=version?analytic:actual;auto [at,inserted]=terms.try_emplace(key,Vector(nm));(void)inserted;
        for(int j=0;j<q.ids.size();++j)
        {
            long double v=q.phi[j]*value;
            for(int axis=0;axis<3;++axis)v+=(long double)q.derivative[axis][j]*flux[axis];
            const double increment=weight*static_cast<double>(v);const int row=offset[block]+q.ids[j];
            at->second[row]+=increment;locals[version][element][key][row]+=increment;
        }
    };
    auto atReference=[&](const RANSWallCoordinate &p,bool outer,RANSVector &u,RANSTensor &G,
        double &k,double &epsilon,double &pi)
    {
        const double z=p[2];u[0]=data.boundary[0](p,time);G[0][2]=outer?us/(model.Kappa()*(z-d)):us*us/ref->Original().TotalViscosity(z);
        k=outer?K0:ref->ReferenceK(z);epsilon=outer?A/(z-d):ref->Dissipation(z);pi=outer?piOuter:level;
    };
    auto stressParts=[&](bool outer,double k,double epsilon,double pi,const RANSTensor &G,
        double nt,const RANSWallCoordinate &p,RANSTensor &molecular,RANSTensor &linear,RANSTensor &nonlinear)
    {
        if(!outer)for(int i=0;i<3;++i)for(int j=0;j<3;++j)
        {molecular[i][j]=g.molecularNu*(G[i][j]+G[j][i]);linear[i][j]=nt*(G[i][j]+G[j][i]);}
        else
        {
            const auto law=OuterClosure(k,epsilon,G,p);
            for(int i=0;i<3;++i)for(int j=0;j<3;++j)
            {nonlinear[i][j]=-law.quadraticStress[i][j];linear[i][j]=-law.reynoldsStress[i][j]+(i==j?2*k/3:0)+law.quadraticStress[i][j];}
        }
        (void)pi;
    };
    auto volume=[&](const VolumePoint &q,bool outer,const EquilibriumAWFProfiles *profile,
        const RANSVector &fp,const RANSVector &fi,int sampler)
    {
        const auto &el=fields[0]->Elements()[q.v.element];const std::string element=std::to_string(q.v.geometry);
        const std::string region=outer?"outer/":"inner/";
        std::vector<std::string> fieldsAt;
        for(int version=0;version<2;++version)
        {
            RANSVector u{};RANSTensor G{},mol{},lin{},nl{};double k,ep,pi;
            if(version)atReference(q.v.coordinate,outer,u,G,k,ep,pi);
            else
            {
                for(int a=0;a<2;++a){u[PhysicalAxis(a)]=q.v.Value(s,offset[a]);G[PhysicalAxis(a)]=q.v.Gradient(s,offset[a]);}
                k=q.v.Value(s,offset[K()]);ep=outer?q.e.Value(s,offset[Epsilon()]):profile->Dissipation(q.v.coordinate[2]);
                pi=q.p.Value(s,offset[InnerPressure()+outer]);
            }
            const double nt=outer?0.:(version?ref.get():profile)->Original().EddyViscosity(q.v.coordinate[2]);
            stressParts(outer,k,ep,pi,G,nt,q.v.coordinate,mol,lin,nl);
            for(int a=0;a<2;++a)
            {
                const int aa=PhysicalAxis(a);RANSVector pressure{},adv{};pressure[aa]=-pi;
                for(int b=0;b<3;++b)adv[b]=-u[aa]*u[b];
                add(version,region+"volume/augmented_pressure",element,a,q.v,0,pressure,q.weight);
                add(version,region+"volume/advection",element,a,q.v,0,adv,q.weight);
                for(const auto &entry:std::vector<std::pair<std::string,const RANSTensor*>>{{"molecular",&mol},{"linear",&lin},{"nonlinear",&nl}})
                {
                    RANSVector normal{},shear{};
                    for(int b=0;b<3;++b)(b==aa?normal[b]:shear[b])=(*entry.second)[aa][b];
                    add(version,region+"volume/"+entry.first+"_normal",element,a,q.v,0,normal,q.weight);
                    add(version,region+"volume/"+entry.first+"_shear",element,a,q.v,0,shear,q.weight);
                }
                const bool active=!outer && q.v.coordinate[2]<=g.roughness;
                add(version,region+"volume/rough_pressure",element,a,q.v,active && !version?-fp[aa]:0,{},q.weight);
                add(version,region+"volume/rough_inertia",element,a,q.v,active && !version?-fi[aa]:0,{},q.weight);
            }
            if(exportPoints)fieldsAt.push_back(J::Object({{"u",J::Numbers(u)},{"G",J::Rows(G)},{"k",J::Number(k)},
                {"epsilon",J::Number(ep)},{"Pi",J::Number(pi)},{"nt_aux",J::Number(nt)},
                {"molecular",J::Rows(mol)},{"linear",J::Rows(lin)},{"nonlinear",J::Rows(nl)}}));
        }
        if(exportPoints)points.push_back(J::Object({{"region",J::Number(outer?1:0)},{"weight",J::Number(q.weight)},
            {"velocity_k",AttributionRow(q.v)},{"pressure",AttributionRow(q.p)},
            {"epsilon_row",outer?AttributionRow(q.e):"null"},{"lower",J::Numbers(el.lower)},{"upper",J::Numbers(el.upper)},
            {"sampler",J::Number(sampler)},{"actual_analytic",J::List(fieldsAt)}}));
    };
    int sampler=0;
    for(const auto &ray:m_rays)
    {
        auto p=ray.coordinate;p[2]=H;const auto pr=Point(OuterPressure(),p,1),kr=Point(K(),p,1),vr=Point(0,p,1);
        const auto gp=pr.Gradient(s,offset[OuterPressure()]),gk=kr.Gradient(s,offset[K()]);
        RANSVector fp{},fi{};fp[0]=gp[0]-2*gk[0]/3;const auto f=CartesianForcing(s,ray.coordinate,time,data,false);
        for(int b=0;b<3;++b)fi[b]=f[b]-fp[b];
        if(exportPoints)samplers.push_back(J::Object({{"coordinate",J::Numbers(p)},{"outer_velocity_k",AttributionRow(vr)},
            {"outer_pressure",AttributionRow(pr)},{"K",J::Number(ray.profile->K())},{"pressure_force",J::Numbers(fp)},{"inertial_force",J::Numbers(fi)}}));
        for(const auto &q:ray.points)volume(q,false,ray.profile.get(),fp,fi,sampler);
        ++sampler;
    }
    for(const auto &q:m_outer)volume(q,true,nullptr,{}, {},-1);
    const auto body=OwnedBodyLoad(time,data);Require(Maximum(body)==0.,"CT10 physical attribution forbids undeclared body forcing");
    for(int version=0;version<2;++version)(version?analytic:actual).emplace("prescribed_body",Vector(nm));
    for(const auto &name:{std::string("outlet"),std::string("top")})
        for(const auto &face:Surface(0,name=="outlet"?0:2,name=="outlet"?upper[0]:upper[2],1))
        {
            const auto &q=face.test;const bool outer=fields[0]->Elements()[q.element].lower[2]>=H;
            const std::string region=outer?"outer/":"inner/",element=std::to_string(q.geometry);const double area=std::hypot(face.area[0],face.area[2]);
            for(int version=0;version<2;++version)
            {
                RANSVector u{};RANSTensor G{};double k,ep,pi;
                if(version)atReference(q.coordinate,outer,u,G,k,ep,pi);
                else for(int a=0;a<2;++a)u[PhysicalAxis(a)]=q.Value(s,offset[a]);
                double un=0;for(int b=0;b<3;++b)un+=u[b]*face.area[b];
                for(int a=0;a<2;++a)
                {
                    const double traction=name=="outlet"?data.outletTraction[outer][a](q.coordinate,time):(a==0?data.topTraction[0](q.coordinate,time):0.);
                    // w at the roof is essential. There is no independent normal traction datum.
                    add(version,region+name+"/natural_traction",element,a,q,-area*traction,{},face.weight);
                    add(version,region+name+"/natural_advection",element,a,q,u[PhysicalAxis(a)]*un,{},face.weight);
                }
            }
            if(exportPoints)facePoints.push_back(J::Object({{"kind",J::Quote(name)},{"region",J::Number(outer)},
                {"row",AttributionRow(q)},{"weight",J::Number(face.weight)},{"area",J::Numbers(face.area)}}));
        }
    for(int side:{-1,1})for(const auto &face:Surface(0,2,H,side,64))
    {
        const auto &q=face.test;const bool outer=fields[0]->Elements()[q.element].lower[2]>=H;
        for(int version=0;version<2;++version)
        {
            RANSVector u{},traction{},flux{};RANSTensor G{},mol{},lin{},nl{};double k,ep,pi;
            if(version)atReference(q.coordinate,outer,u,G,k,ep,pi);
            else
            {
                for(int a=0;a<2;++a){u[PhysicalAxis(a)]=q.Value(s,offset[a]);G[PhysicalAxis(a)]=q.Gradient(s,offset[a]);}
                k=q.Value(s,offset[K()]);ep=outer?GeometryPoint(Epsilon(),q.geometry,q.coordinate).Value(s,offset[Epsilon()]):0.;
                pi=GeometryPoint(InnerPressure()+outer,q.geometry,q.coordinate).Value(s,offset[InnerPressure()+outer]);
            }
            const double nt=outer?0.:(version?ref:Profile(SampledK(s,q.coordinate),q.coordinate))->Original().EddyViscosity(H);
            stressParts(outer,k,ep,pi,G,nt,q.coordinate,mol,lin,nl);double un=0;for(int b=0;b<3;++b)un+=u[b]*face.area[b];
            for(int a:{0,2})for(int b=0;b<3;++b)traction[a]+=(mol[a][b]+lin[a][b]+nl[a][b]-(a==b?pi:0))*face.area[b];
            for(int a:{0,2})flux[a]=u[a]*un-traction[a];
            interfaces.push_back(J::Object({{"region",J::Number(outer)},{"analytic",J::Boolean(version)},
                {"side",J::Number(side)},{"row",AttributionRow(q)},{"weight",J::Number(face.weight)},
                {"area",J::Numbers(face.area)},{"u",J::Numbers(u)},{"traction",J::Numbers(traction)},{"outward_total_flux",J::Numbers(flux)}}));
        }
    }
    const auto native=Residual(s,s,0.,time,data,false,false);Vector total(nm),totalA(nm),error(nm);
    std::vector<std::pair<std::string,std::string>> actualJson,analyticJson;std::vector<std::string> localJson;
    for(const auto &[key,v]:actual){actualJson.push_back({key,J::Numbers(v)});for(int i=0;i<nm;++i)total[i]+=v[i];}
    for(const auto &[key,v]:analytic){analyticJson.push_back({key,J::Numbers(v)});for(int i=0;i<nm;++i)totalA[i]+=v[i];}
    for(int i=0;i<nm;++i)error[i]=total[i]-native[i];
    Require(Maximum(error)<=2e-10,"CT10 signed momentum components do not reconstruct original residual");
    for(int version=0;version<2;++version)for(const auto &[geometry,terms]:locals[version])
    {
        std::vector<std::pair<std::string,std::string>> entries;
        for(const auto &[key,load]:terms)
        {std::vector<int> ids;Vector values;for(auto [id,value]:load){ids.push_back(id);values.push_back(value);}entries.push_back({key,J::Object({{"ids",J::Numbers(ids)},{"values",J::Numbers(values)}})});}
        localJson.push_back(J::Object({{"geometry",J::Quote(geometry)},{"analytic",J::Boolean(version)},{"terms",J::Object(entries)}}));
    }
    auto constraints=Constraints(s,time,data);std::vector<int> free,fixed;Vector constrained(total),reactions(nm);
    for(int i=0;i<nm;++i)if(constraints.count(i)){fixed.push_back(i);reactions[i]=native[i];constrained[i]=s[i]-constraints.at(i);}else free.push_back(i);
    std::vector<std::string> elements;for(const auto &e:fields[0]->Elements())elements.push_back(J::Object({{"geometry",J::Number(e.geometry)},{"lower",J::Numbers(e.lower)},{"upper",J::Numbers(e.upper)}}));
    return J::Object({{"schema",J::Quote("CT10SignedMomentumV1")},{"identity",J::Quote(PreservationIdentity())},
        {"read_only","true"},{"momentum_size",J::Number(nm)},{"offsets",J::Numbers(offset)},
        {"actual",J::Object(actualJson)},{"analytic",J::Object(analyticJson)},{"total_actual",J::Numbers(total)},
        {"total_analytic",J::Numbers(totalA)},{"native_raw",J::Numbers(Vector(native.begin(),native.begin()+nm))},
        {"constrained",J::Numbers(constrained)},{"prescribed_reactions",J::Numbers(reactions)},
        {"free_rows",J::Numbers(free)},{"fixed_rows",J::Numbers(fixed)},{"shared_velocity_trace_ids",J::List(shared)},
        {"elements",J::List(elements)},{"local_sparse",J::List(localJson)},{"volume_points",J::List(points)},
        {"samplers",J::List(samplers)},{"natural_faces",J::List(facePoints)},{"interface_fluxes",J::List(interfaces)},
        {"maximum_reconstruction_error",J::Number(Maximum(error))},
        {"scope",J::Quote("same native generating owners,coordinates,weights and test rows; sparse element support only; interface fluxes are diagnostic and not added twice; no state/history/boundary change or governing solve")}});
}
}
