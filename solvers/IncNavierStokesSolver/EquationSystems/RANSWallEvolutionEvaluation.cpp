// SPDX-License-Identifier: MIT
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <cstring>
#include <cmath>
#include <list>
#include <tuple>
#include <unordered_map>

namespace Nektar
{
using namespace RANSWallAlgebra;
using namespace SolverUtils;
namespace J=RANSWallSolveJson;
namespace
{
using Bits=std::uint64_t;
using Key=std::vector<Bits>;
using RowKey=std::array<Bits,7>;
struct RowHash {std::size_t operator()(const RowKey &key) const noexcept
    {std::size_t h=1469598103934665603ULL;for(Bits x:key){h^=x;h*=1099511628211ULL;}return h;}};
Bits Bit(double x){Bits b;static_assert(sizeof b==sizeof x);std::memcpy(&b,&x,sizeof b);return b;}
bool Exact(const Vector &a,const Vector &b)
{return a.size()==b.size() && (a.empty() || std::memcmp(a.data(),b.data(),a.size()*sizeof(double))==0);}
std::string Fixed(const RANSWallEvolutionData &data)
{
    std::string s;
    for(const auto &[id,value]:data.fixedBoundary)
    {s.append(reinterpret_cast<const char*>(&id),sizeof id);s.append(reinterpret_cast<const char*>(&value),sizeof value);}
    return s;
}
Key DomainKey(bool only,bool coverage,const std::vector<RANSWallCoordinate> &anchors)
{
    Key k{Bits(only),Bits(coverage),Bits(anchors.size())};
    for(const auto &p:anchors)for(double x:p)k.push_back(Bit(x));
    return k;
}
}

struct RANSWallEvolutionPatch::EvaluationCache
{
    static constexpr std::size_t maximumStates=4,maximumRows=524288,maximumProfiles=4096,maximumQueries=32;
    const RANSWallEvolutionData *data=nullptr;
    std::uint64_t generation=0;Bits time=0;std::string fixed;bool ready=false;
    struct Domain
    {
        std::shared_ptr<const ScalarDomainResult> result;
        std::shared_ptr<const std::set<RANSWallCoordinate>> locations,functionals;
    };
    struct ResidualEntry {Key key;Vector history;std::shared_ptr<const Vector> value;};
    struct State
    {
        Vector state;
        std::map<Key,Domain> domains;
        std::map<Key,std::shared_ptr<const std::map<int,double>>> constraints;
        std::vector<ResidualEntry> residuals;
        std::map<std::tuple<const EquilibriumAWFProfiles*,Bits,Bits,Bits,bool>,
            std::shared_ptr<const RANSCappedDestructionValues>> normalization;
    };
    std::list<std::shared_ptr<State>> states;
    std::unordered_map<RowKey,std::shared_ptr<const RANSWallGlobalPoint>,RowHash> rows;
    std::map<const EquilibriumAWFProfiles*,std::shared_ptr<EquilibriumAWFProfiles>> profiles;
    std::map<std::pair<const EquilibriumAWFProfiles*,Bits>,std::shared_ptr<EquilibriumAWFProfiles>> reweighted;
    const void *coverageToken=nullptr;
    State &Find(const Vector &state,EvaluationWork &work)
    {
        for(auto it=states.begin();it!=states.end();++it)if(Exact((*it)->state,state))
        {auto result=*it;states.splice(states.begin(),states,it);return *result;}
        auto result=std::make_shared<State>();result->state=state;states.push_front(result);
        if(states.size()>maximumStates){states.pop_back();coverageToken=nullptr;++work.stateEvictions;}
        work.peakStates=std::max<std::uint64_t>(work.peakStates,states.size());return *result;
    }
    void Clear()
    {states.clear();rows.clear();profiles.clear();reweighted.clear();coverageToken=nullptr;}
};

RANSWallEvolutionPatch::EvaluationScope::EvaluationScope(const RANSWallEvolutionPatch &p,
    const RANSWallEvolutionData &data,double time):patch(p)
{
    if(!patch.m_evaluationReuse)return;
    if(patch.m_evaluation && patch.m_evaluation->data==&data)
    {patch.PrepareEvaluation(data,time);return;}
    previous=patch.m_evaluation;patch.m_evaluation=std::make_shared<EvaluationCache>();
    patch.m_evaluation->data=&data;owns=true;++patch.m_evaluationWork.scopes;
    patch.PrepareEvaluation(data,time);
}
RANSWallEvolutionPatch::EvaluationScope::~EvaluationScope()
{
    if(owns){patch.m_evaluation=previous;patch.InvalidateCoverageToken();}
}
void RANSWallEvolutionPatch::PrepareEvaluation(const RANSWallEvolutionData &data,double time) const
{
    if(!m_evaluation)return;
    Require(m_evaluation->data==&data && std::isfinite(time),"evaluation scope/input identity differs");
    const auto fixed=Fixed(data);auto &c=*m_evaluation;
    if(!c.ready || c.time!=Bit(time) || c.generation!=data.generation || c.fixed!=fixed)
    {
        if(c.ready)++m_evaluationWork.inputInvalidations;
        c.Clear();c.time=Bit(time);c.generation=data.generation;c.fixed=fixed;c.ready=true;
    }
}
void RANSWallEvolutionPatch::InvalidateCoverageToken() const
{if(m_evaluation)m_evaluation->coverageToken=nullptr;}

int RANSWallEvolutionPatch::GeometryElement(int block,int geometry) const
{
    RANSWallTiming::Scope timing("geometry.owner_catalog",true,true);
    if(m_evaluationReuse)
    {
        if(m_geometryCatalog.empty())
        {
            m_geometryCatalog.resize(fields.size());
            for(int b=0;b<fields.size();++b)for(int e=0;e<fields[b]->Elements().size();++e)
                Require(m_geometryCatalog[b].emplace(fields[b]->Elements()[e].geometry,e).second,
                    "duplicate native geometry in evaluation catalog");
        }
        const auto it=m_geometryCatalog.at(block).find(geometry);
        Require(it!=m_geometryCatalog.at(block).end(),"observation geometry absent from requested native field");return it->second;
    }
    for(int e=0;e<fields[block]->Elements().size();++e)if(fields[block]->Elements()[e].geometry==geometry)return e;
    throw std::runtime_error("observation geometry absent from requested native field");
}
std::shared_ptr<const RANSWallGlobalPoint> RANSWallEvolutionPatch::EvaluationRow(int block,int element,
    const RANSWallCoordinate &coordinate) const
{
    RANSWallTiming::Scope timing("geometry.evaluation_row",true,true);
    const auto &el=fields.at(block)->Elements().at(element);
    const int side=coordinate[2]==m_parameters.geometry.height?
        (el.upper[2]<=m_parameters.geometry.height?-1:1):0;
    // A patch owns an immutable field/mesh/map generation. Both geometry and
    // native element, actual radial side and all coordinate bits enter the key.
    const RowKey key{Bits(block),Bits(el.geometry),Bits(element),Bits(side+1),Bit(coordinate[0]),Bit(coordinate[1]),Bit(coordinate[2])};
    if(m_evaluation)
    {
        const auto it=m_evaluation->rows.find(key);
        if(it!=m_evaluation->rows.end()){++m_evaluationWork.rowHits;return it->second;}
    }
    ++m_evaluationWork.rowConstructs;
    auto row=std::make_shared<const RANSWallGlobalPoint>(Transform(fields[block]->ElementPoint(element,coordinate)));
    if(m_evaluation)
    {
        if(m_evaluation->rows.size()>=EvaluationCache::maximumRows)
        {m_evaluation->rows.clear();++m_evaluationWork.rowEvictions;}
        m_evaluation->rows.emplace(key,row);
        m_evaluationWork.peakRows=std::max<std::uint64_t>(m_evaluationWork.peakRows,m_evaluation->rows.size());
    }
    return row;
}

std::shared_ptr<const RANSWallEvolutionPatch::ScalarDomainResult> RANSWallEvolutionPatch::CachedDomain(
    const Vector &s,double time,const RANSWallEvolutionData &data,bool only,
    const std::vector<RANSWallCoordinate> &anchors) const
{
    if(!m_evaluation)return nullptr;
    PrepareEvaluation(data,time);auto &c=*m_evaluation;auto &entry=c.Find(s,m_evaluationWork);
    const auto key=DomainKey(only,m_assertKCoverage,anchors);
    auto it=entry.domains.find(key);
    if(it==entry.domains.end() && !NewModels() && AcceptedAWF() && !m_transitionEnabled && !m_cappedDestruction &&
        s.size()==std::size_t(Size()) && std::all_of(s.begin(),s.end(),[](double x){return std::isfinite(x);}))
    {
        // The accepted-AWF scalar inventory and dependent trace are functions
        // only of the complete k/epsilon blocks, immutable geometry/data and
        // time. It has no chi/sensor or capped-normalization dependence on U/p.
        // Full state contexts remain distinct; only proven-identical scalar
        // evaluations are shared, and both returned full states are rebound.
        const int begin=offset[K()],count=offset[InnerPressure()]-begin;
        for(const auto &other:c.states)if(other.get()!=&entry && other->state.size()==s.size() &&
            std::memcmp(s.data()+begin,other->state.data()+begin,count*sizeof(double))==0)
        {
            const auto old=other->domains.find(key);if(old==other->domains.end())continue;
            auto shared=old->second;auto result=std::make_shared<ScalarDomainResult>(*shared.result);
            if(c.coverageToken!=shared.locations.get())
            {m_coveredKLocations=*shared.locations;m_coveredKFunctionals=*shared.functionals;c.coverageToken=shared.locations.get();}
            result->effectiveState=ExternalLift(s,time,data);
            if(!only)
            {result->traceState=result->effectiveState;for(const auto &[row,value]:Constraints(s,time,data))result->traceState[row]=value;}
            shared.result=result;entry.domains.emplace(key,std::move(shared));++m_evaluationWork.scalarFieldReuses;
            break;
        }
        it=entry.domains.find(key);
    }
    if(it==entry.domains.end())return nullptr;
    if(c.coverageToken!=it->second.locations.get())
    {m_coveredKLocations=*it->second.locations;m_coveredKFunctionals=*it->second.functionals;c.coverageToken=it->second.locations.get();}
    return it->second.result;
}
void RANSWallEvolutionPatch::StoreDomain(const Vector &s,double time,const RANSWallEvolutionData &data,
    bool only,const std::vector<RANSWallCoordinate> &anchors,const ScalarDomainResult &result) const
{
    if(!m_evaluation || !result.admissible)return;
    PrepareEvaluation(data,time);auto &c=*m_evaluation;auto &entry=c.Find(s,m_evaluationWork);
    if(entry.domains.size()>=EvaluationCache::maximumQueries){entry.domains.clear();c.coverageToken=nullptr;}
    EvaluationCache::Domain value{std::make_shared<const ScalarDomainResult>(result),
        std::make_shared<const std::set<RANSWallCoordinate>>(m_coveredKLocations),
        std::make_shared<const std::set<RANSWallCoordinate>>(m_coveredKFunctionals)};
    c.coverageToken=value.locations.get();entry.domains[DomainKey(only,m_assertKCoverage,anchors)]=std::move(value);
}
RANSWallEvolutionPatch::ScalarDomainResult RANSWallEvolutionPatch::ScalarDomain(
    const Vector &s,double time,const RANSWallEvolutionData &data,bool retain,bool only,
    const std::vector<RANSWallCoordinate> &anchors) const
{
    RANSWallTiming::Scope timing("domain.scalar");EvaluationScope scope(*this,data,time);++m_evaluationWork.domainCalls;
    if(auto cached=CachedDomain(s,time,data,only,anchors))
    {++m_evaluationWork.domainHits;auto result=*cached;if(!retain)result.observations.clear();
    return result;}
    ++m_evaluationWork.domainComputes;
    auto result=EvaluateScalarDomain(s,time,data,retain || bool(m_evaluation),only,anchors);
    StoreDomain(s,time,data,only,anchors,result);
    if(!retain)result.observations.clear();
    return result;
}

std::shared_ptr<EquilibriumAWFProfiles> RANSWallEvolutionPatch::EvaluationProfile(
    const std::shared_ptr<EquilibriumAWFProfiles> &base,double weight) const
{
    auto construct=[&](){return weight==0 || weight==1?base:std::make_shared<EquilibriumAWFProfiles>(base->Reweighted(weight));};
    if(!m_evaluation)return construct();
    auto &c=*m_evaluation;
    if(c.profiles.size()+2>EvaluationCache::maximumProfiles)
    {
        for(auto &s:c.states)s->normalization.clear();
        c.reweighted.clear();c.profiles.clear();
    }
    c.profiles[base.get()]=base;
    std::shared_ptr<EquilibriumAWFProfiles> result;
    if(weight==0 || weight==1)result=base;
    else
    {
        auto &item=c.reweighted[{base.get(),Bit(weight)}];
        if(item){result=item;++m_evaluationWork.profileReuses;}else{result=construct();item=result;}
    }
    c.profiles[result.get()]=result;m_evaluationWork.peakProfiles=std::max<std::uint64_t>(m_evaluationWork.peakProfiles,c.profiles.size());return result;
}
std::shared_ptr<const RANSCappedDestructionValues> RANSWallEvolutionPatch::CachedNormalization(
    const Vector &s,const RANSWallCoordinate &p,const EquilibriumAWFProfiles &profile) const
{
    ++m_evaluationWork.normalizationCalls;
    if(!m_evaluation || Mapped() || !m_evaluation->profiles.count(&profile))return nullptr;
    auto &c=*m_evaluation;if(c.generation!=c.data->generation){PrepareEvaluation(*c.data,[] (Bits b){double t;std::memcpy(&t,&b,sizeof t);return t;}(c.time));return nullptr;}
    auto &entry=c.Find(s,m_evaluationWork);
    const auto key=std::make_tuple(&profile,Bit(p[0]),Bit(p[1]),Bit(m_mapping.At(p).scale),m_assertKCoverage);
    const auto it=entry.normalization.find(key);
    if(it==entry.normalization.end())return nullptr;
    // Coverage is a live acceptance assertion, even when numerical ray work
    // is reused. The owned immutable profile identifies the complete branch,
    // quadrature, coefficients and geometry, not merely sampled K.
    if(m_assertKCoverage)for(const auto &q:profile.Quadrature())
    {auto point=p;point[2]=q.r/m_mapping.At(p).scale;Require(m_coveredKLocations.count(point),"unprotected constitutive k evaluation: cached normalization ray");}
    ++m_evaluationWork.normalizationHits;return it->second;
}
void RANSWallEvolutionPatch::StoreNormalization(const Vector &s,const RANSWallCoordinate &p,
    const EquilibriumAWFProfiles &profile,const RANSCappedDestructionValues &value) const
{
    if(!m_evaluation || Mapped() || !m_evaluation->profiles.count(&profile))return;
    auto &entry=m_evaluation->Find(s,m_evaluationWork);
    entry.normalization[{&profile,Bit(p[0]),Bit(p[1]),Bit(m_mapping.At(p).scale),m_assertKCoverage}]=
        std::make_shared<const RANSCappedDestructionValues>(value);
}

std::shared_ptr<const Vector> RANSWallEvolutionPatch::CachedResidual(const Vector &s,const Vector &h,
    double alpha,double time,const RANSWallEvolutionData &data,bool anchored,bool constraints,int region) const
{
    ++m_evaluationWork.residualCalls;
    if(!m_evaluation || m_assertKCoverage || m_solveDiagnostics)return nullptr;
    PrepareEvaluation(data,time);const Key key{Bit(alpha),Bits(anchored),Bits(constraints),Bits(region+1)};
    auto &entry=m_evaluation->Find(s,m_evaluationWork);
    for(const auto &item:entry.residuals)if(item.key==key && Exact(item.history,h))
    {++m_evaluationWork.residualHits;return item.value;}
    return nullptr;
}
void RANSWallEvolutionPatch::StoreResidual(const Vector &s,const Vector &h,double alpha,double time,
    const RANSWallEvolutionData &data,bool anchored,bool constraints,int region,const Vector &r) const
{
    if(!m_evaluation || m_assertKCoverage || m_solveDiagnostics)return;
    PrepareEvaluation(data,time);auto &entry=m_evaluation->Find(s,m_evaluationWork);
    if(entry.residuals.size()>=EvaluationCache::maximumQueries)entry.residuals.clear();
    entry.residuals.push_back({{Bit(alpha),Bits(anchored),Bits(constraints),Bits(region+1)},h,std::make_shared<const Vector>(r)});
}
std::shared_ptr<const std::map<int,double>> RANSWallEvolutionPatch::CachedConstraints(const Vector &s,
    double time,const RANSWallEvolutionData &data) const
{
    ++m_evaluationWork.constraintCalls;
    if(!m_evaluation || m_assertKCoverage)return nullptr;
    PrepareEvaluation(data,time);auto &entry=m_evaluation->Find(s,m_evaluationWork);
    const auto it=entry.constraints.find({});if(it==entry.constraints.end())return nullptr;
    ++m_evaluationWork.constraintHits;return it->second;
}
void RANSWallEvolutionPatch::StoreConstraints(const Vector &s,double time,const RANSWallEvolutionData &data,
    const std::map<int,double> &value) const
{
    if(!m_evaluation || m_assertKCoverage)return;
    PrepareEvaluation(data,time);m_evaluation->Find(s,m_evaluationWork).constraints[{}]=std::make_shared<const std::map<int,double>>(value);
}

std::string RANSWallEvolutionPatch::EvaluationReuseJSON() const
{
    const auto &w=m_evaluationWork;std::size_t fixedRows=0,surfaceRows=0;
    for(const auto &[key,rows]:m_reconstructionCatalog){(void)key;fixedRows+=rows.size();}
    for(const auto &[key,rows]:m_surfaceCatalog){(void)key;surfaceRows+=rows.size();}
    return J::Object({{"mode",J::Quote(m_evaluationReuse?"ExactStepV1":"Off")},
        {"scopes",J::Number(w.scopes)},{"domain_calls",J::Number(w.domainCalls)},
        {"domain_computes",J::Number(w.domainComputes)},{"domain_hits",J::Number(w.domainHits)},
        {"prefix_hits",J::Number(w.prefixHits)},{"row_constructs",J::Number(w.rowConstructs)},{"row_hits",J::Number(w.rowHits)},
        {"row_evictions",J::Number(w.rowEvictions)},{"reconstruction_builds",J::Number(w.reconstructionBuilds)},
        {"reconstruction_hits",J::Number(w.reconstructionHits)},{"surface_builds",J::Number(w.surfaceBuilds)},{"surface_hits",J::Number(w.surfaceHits)},
        {"normalization_calls",J::Number(w.normalizationCalls)},{"normalization_computes",J::Number(w.normalizationComputes)},
        {"normalization_hits",J::Number(w.normalizationHits)},{"residual_calls",J::Number(w.residualCalls)},
        {"residual_hits",J::Number(w.residualHits)},{"constraint_calls",J::Number(w.constraintCalls)},
        {"constraint_hits",J::Number(w.constraintHits)},{"state_evictions",J::Number(w.stateEvictions)},
        {"input_invalidations",J::Number(w.inputInvalidations)},{"profile_reuses",J::Number(w.profileReuses)},
        {"scalar_field_reuses",J::Number(w.scalarFieldReuses)},
        {"fixed_reconstruction_rows",J::Number(fixedRows)},{"fixed_surface_rows",J::Number(surfaceRows)},
        {"row_counter_scope",J::Quote("scoped point/derivative requests; fixed reconstruction and surface catalogs counted separately")},
        {"peak_states",J::Number(w.peakStates)},{"peak_rows",J::Number(w.peakRows)},{"peak_profiles",J::Number(w.peakProfiles)},
        {"maximum_states",J::Number(EvaluationCache::maximumStates)},{"maximum_rows",J::Number(EvaluationCache::maximumRows)},
        {"maximum_profiles",J::Number(EvaluationCache::maximumProfiles)},
        {"lifetime",J::Quote("immutable fixed catalogs per patch; exact-state, moving-row and owned-profile caches per evaluation scope")}});
}
}
