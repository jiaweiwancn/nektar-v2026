///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <filesystem>
#include <numeric>

namespace Nektar
{
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;

RANSWallEvolutionPatch::SparseFlowResult RANSWallEvolutionPatch::SparseFlowBlocks(
    const PicardBlock &flow,int nv) const
{
    Require(SerialSparse() && flow.sparse,"invalid sparse flow blocks");
    const int n=flow.rows.size();SparseFlowResult r;
    r.a=flow.csr.Block(0,nv,0,nv);r.g=flow.csr.Block(0,nv,nv,n);
    r.d=flow.csr.Block(nv,n,0,nv);r.c=flow.csr.Block(nv,n,nv,n);
    r.m=SparseMatrix(nv,nv);std::vector<int> index(Size(),-1);
    for(int i=0;i<nv;++i)index[flow.rows[i]]=i;
    for(int i=0;i<nv;++i)
    {
        const int b=flow.rows[i]/size[0],local=flow.rows[i]-offset[b];
        m_sparseMass[b].Each(local,[&](int j,double v){const int k=index[offset[b]+j];if(k>=0)r.m.Entry(i,k)=v;});
    }
    r.m.Finalize();return r;
}

RANSWallEvolutionPatch::SparseFlowResult RANSWallEvolutionPatch::SparseFlowCorrection(
    const PicardBlock &flow,int nv,double alpha,const std::function<void(const std::string&)> &event) const
{
    RANSWallTiming::Scope timing("sparse.flow_linear_operation");
    Require(alpha>0,"invalid sparse flow alpha");
    const int n=flow.rows.size(),np=n-nv;auto r=SparseFlowBlocks(flow,nv);
    if(event)event("momentum_attempt");
    std::vector<int> groups;
    if(PredictorOrdering()==SparsePredictorOrdering::VectorCOLAMDV1)
    {groups.resize(nv);for(int i=0;i<nv;++i)groups[i]=flow.rows[i]%size[0];}
    SparseFactor momentum(r.a,false,PredictorOrdering(),groups),mass(r.m,true);
    r.predictorFill=momentum.FactorNonZeros();r.massFill=mass.FactorNonZeros();
    Vector rhs(nv);for(int i=0;i<nv;++i)rhs[i]=-flow.residual[i];
    r.predictorRhsMaximum=Maximum(rhs);
    r.predictor=momentum.Apply(rhs);r.predictorDefect=r.a.Apply(r.predictor);
    if(event)event("momentum_complete");
    for(int i=0;i<nv;++i)r.predictorDefect[i]-=rhs[i];
    auto pressureRhs=r.d.Apply(r.predictor);
    for(int i=0;i<np;++i)pressureRhs[i]+=flow.residual[nv+i];
    r.pressureRhsMaximum=Maximum(pressureRhs);
    // Local sparse product with DIAGONAL mass appears ONLY in the pressure
    // preconditioner. The governing action below always uses consistent M.
    SparseMatrix approximate(np,np);
    for(int i=0;i<np;++i)
    {
        r.d.Each(i,[&](int j,double v)
        {
            const double diagonal=r.m.Get(j,j);Require(diagonal>0,"nonpositive mass diagonal in pressure preconditioner");
            r.g.Each(j,[&](int k,double w){approximate.Entry(i,k)+=v*w/(alpha*diagonal);});
        });
        r.c.Each(i,[&](int j,double v){approximate.Entry(i,j)-=v;});
    }
    approximate.Finalize();SparseFactor prec(approximate);r.preconditionerFill=prec.FactorNonZeros();
    const auto action=[&](const Vector &p)
    {
        ++r.massApplications;if(event)event("mass");auto v=mass.Apply(r.g.Apply(p));for(auto &x:v)x/=alpha;
        auto out=r.d.Apply(v);const auto c=r.c.Apply(p);for(int i=0;i<np;++i)out[i]-=c[i];return out;
    };
    if(event)event("pressure_attempt");
    const auto solved=SparsePressureGMRES(m_sparseSession,m_comm,pressureRhs,action,[&](const Vector &x){return prec.Apply(x);});
    if(event)event("pressure_complete");
    r.dz=solved.solution;r.iterations=solved.iterations;r.pressureDefect=action(r.dz);
    for(int i=0;i<np;++i)r.pressureDefect[i]-=pressureRhs[i];
    r.Gz=r.g.Apply(r.dz);++r.massApplications;if(event)event("mass");r.bGz=mass.Apply(r.Gz);for(auto &x:r.bGz)x/=alpha;
    r.du=r.predictor;for(int i=0;i<nv;++i)r.du[i]-=r.bGz[i];
    auto update=r.du;update.insert(update.end(),r.dz.begin(),r.dz.end());r.defect=flow.Apply(update);
    for(int i=0;i<n;++i)r.defect[i]+=flow.residual[i];
    const auto ab=r.a.Apply(r.bGz);r.identity.resize(nv);r.constraints.resize(np);
    for(int i=0;i<nv;++i)r.identity[i]=r.defect[i]-(r.Gz[i]-ab[i]);
    for(int i=0;i<np;++i)r.constraints[i]=r.defect[nv+i];
    return r;
}

std::string RANSWallEvolutionPatch::PreFactorFrozenPacket(const Vector &raw,const Vector &history,
    double alpha,double time,const RANSWallEvolutionData &data,const std::string &prefix) const
{
    Require(SerialSparse() && !Mapped() && PimpleSingleSweep() && alpha>0 &&
        raw.size()==Size() && history.size()==Size() && m_exportReferenceOnly &&
        SparseFactorizationCount()==0,"invalid pre-factor frozen request");
    const auto state=ExternalLift(raw,time,data);
    auto external=ExternalConstraints(time,data);
    for(int id:m_traceIds)external.erase(offset[Epsilon()]+id);
    const auto flow=AssemblePicardBlock(0,state,history,alpha,time,data,external);
    const int nv=std::lower_bound(flow.rows.begin(),flow.rows.end(),DifferentialSize())-flow.rows.begin();
    const auto r=SparseFlowBlocks(flow,nv);
    const auto base=std::filesystem::absolute(prefix);std::filesystem::create_directories(base.parent_path());
    const auto save=[&](const SparseMatrix &matrix,const std::string &name)
    {
        const auto file=base.string()+"."+name+".csr";matrix.Write(file);
        return J::Object({{"file",J::Quote(file)},{"sha256",J::Quote(RANSWallR3::HashFile(file))},
            {"rows",J::Number(matrix.Rows())},{"columns",J::Number(matrix.Columns())},
            {"nnz",J::Number(matrix.NonZeros())},{"encoding",J::Quote("RANSCSR1-native-little-endian")}});
    };
    Vector lift=state;for(int i=0;i<Size();++i)lift[i]-=raw[i];
    std::vector<std::string> maps;
    for(int b=0;b<Gauge();++b)
    {
        const auto f=fields[b]->Field();const auto map=f->GetLocalToGlobalMap();
        std::vector<int> ids;Vector signs;std::vector<std::string> owners;
        for(int e=0;e<f->GetExpSize();++e)
        {
            const int start=f->GetCoeff_Offset(e),count=f->GetExp(e)->GetNcoeffs();
            for(int j=0;j<count;++j){ids.push_back(map->GetLocalToGlobalMap(start+j));signs.push_back(map->GetLocalToGlobalSign(start+j));}
            owners.push_back(J::Object({{"geometry",J::Number(fields[b]->Elements()[e].geometry)},
                {"local_start",J::Number(start)},{"coefficients",J::Number(count)}}));
        }
        maps.push_back(J::Object({{"block",J::Number(b)},{"offset",J::Number(offset[b])},
            {"local_to_global",J::Numbers(ids)},{"signs",J::Numbers(signs)},{"owners",J::List(owners)}}));
    }
    Require(SparseFactorizationCount()==0,"frozen assembly unexpectedly factorized a global system");
    std::vector<std::string> pressureMass;
    for(int b:{InnerPressure(),OuterPressure()})pressureMass.push_back(J::Object({
        {"block",J::Number(b)},{"mass",save(m_sparseMass[b],b==InnerPressure()?"pressure-inner-mass":"pressure-outer-mass")},
        {"constant",J::Numbers(m_constants[b])},{"physical_volume",J::Number(m_volumes[b])}}));
    return J::Object({{"schema",J::Quote("PreFactorFrozenCSRExportV2")},{"read_only","true"},
        {"whole_command_sparse_factorizations",J::Number(SparseFactorizationCount())},
        {"reference_mass_dual_norms",J::Quote("not evaluated: diagnostic capacity only; required separately for physical admission")},
        {"factorizations","0"},{"updates","0"},{"time",J::Number(time)},{"alpha",J::Number(alpha)},
        {"identity",J::Quote(PreservationIdentity())},{"descriptor",J::Quote(m_parameters.model.Signature())},
        {"stored_slots",J::Number(Size())},{"block_sizes",J::Numbers(size)},{"offsets",J::Numbers(offset)},
        {"ownership",J::Quote("GeneratedVolumeOwnerV1; canonical signed native maps; traces keep declared sides")},
        {"input",J::Numbers(raw)},{"effective_input",J::Numbers(state)},{"lifting",J::Numbers(lift)},
        {"history",J::Numbers(history)},{"history_sha256",J::Quote(RANSWallR3::Hash(history))},
        {"input_sha256",J::Quote(RANSWallR3::Hash(raw))},{"rows",J::Numbers(flow.rows)},
        {"velocity_size",J::Number(nv)},{"external_constraints",J::Constraints(external)},
        {"fixed_reference",J::Numbers(m_reference)},{"fixed_boundary",J::Constraints(data.fixedBoundary)},
        {"pressure_physical_mass",J::List(pressureMass)},
        {"frozen_flow_residual",J::Numbers(flow.residual)},{"signed_field_maps",J::List(maps)},
        {"frozen_flow_matrix",save(flow.csr,"system")},{"free_velocity_mass",save(r.m,"mass")},
        {"A",save(r.a,"A")},{"G",save(r.g,"G")},{"D",save(r.d,"D")},{"C",save(r.c,"C")},
        {"scope",J::Quote("actual native assembly/restriction and actual supplied history, before any factor or update; no solve, no history substitution")}});
}

std::string RANSWallEvolutionPatch::SparseBackendPacket(const Vector &state,double time,
    const RANSWallEvolutionData &data,double alpha) const
{
    Require(SerialSparse(),"sparse backend packet requires its explicit policy");
    auto external=ExternalConstraints(time,data);for(int id:m_traceIds)external.erase(offset[Epsilon()]+id);
    const auto flow=AssemblePicardBlock(0,state,state,alpha,time,data,external);
    const int nv=std::lower_bound(flow.rows.begin(),flow.rows.end(),DifferentialSize())-flow.rows.begin();
    const auto r=SparseFlowCorrection(flow,nv,alpha);
    Vector pressureDirection(Size());for(int j=0;j<r.dz.size();++j)pressureDirection[flow.rows[nv+j]]=r.dz[j];
    const auto nativeG=AffinePressureAction(pressureDirection,external);Vector actionDefect(nv);
    for(int i=0;i<nv;++i)actionDefect[i]=nativeG[flow.rows[i]]-r.Gz[i];
    Require(Maximum(actionDefect)<1e-10,"sparse capacity omitted complete native pressure action");
    auto massDefect=r.m.Apply(r.bGz);for(int i=0;i<nv;++i)massDefect[i]=alpha*massDefect[i]-r.Gz[i];
    std::vector<std::string> masses,scalars;
    for(int b=0;b<Gauge();++b)masses.push_back(J::Object({{"rows",J::Number(size[b])},{"nnz",J::Number(m_sparseMass[b].NonZeros())}}));
    for(int b:{Epsilon(),K()})
    {
        const auto block=AssemblePicardBlock(b,state,state,alpha,time,data,external);SparseFactor factor(block.csr);
        auto rhs=block.residual;for(auto &x:rhs)x=-x;const auto delta=factor.Apply(rhs);
        auto error=block.Apply(delta);for(int i=0;i<error.size();++i)error[i]-=rhs[i];
        scalars.push_back(J::Object({{"block",J::Number(b)},{"rows",J::Number(block.rows.size())},{"nnz",J::Number(block.csr.NonZeros())},
            {"factor_nnz",J::Number(factor.FactorNonZeros())},{"rhs_maximum",J::Number(Maximum(rhs))},
            {"maximum_true_residual",J::Number(Maximum(error))},{"diagnostic_increment",J::Numbers(delta)}}));
    }
    // Bounded aggregate reference errors on the actual native mass rows.
    // No per-volume dense packets or changed physical/reference coefficients.
    double us,level;m_sparseSession->LoadParameter("RANSWallFrictionVelocity",us);
    m_sparseSession->LoadParameter("RANSOpenPressureLevel",level,0.);
    const auto &model=m_parameters.model;const double H=m_parameters.geometry.height;
    const double K0=model.Q()*us*us,A=us*us*us/model.Kappa(),d=AerodynamicDisplacement();
    const auto referenceProfile=Profile(K0,{lower[0],0.,0.});
    SolverUtils::RANSTensor referenceG{};referenceG[0][2]=us/(model.Kappa()*(H-d));
    const auto referenceLaw=OuterClosure(K0,A/(H-d),referenceG,{lower[0],0.,H});
    std::array<long double,5> square{},gradientSquare{},volume{};
    Vector largest(5),scales{us,us,K0,A/(H-d),us*us};
    for(const auto &q:Quadrature(0))
    {
        const bool outer=q.coordinate[2]>=H;const double z=q.coordinate[2];
        const double du=outer?us/(model.Kappa()*(z-d)):us*us/referenceProfile->Original().TotalViscosity(z);
        const double dk=outer?0.:referenceProfile->ReferenceKGradient(z);
        for(int b=0;b<5;++b)
        {
            if(b==3 && !outer)continue;
            const int block=b==4?InnerPressure()+outer:b;
            const auto row=block==0 || block==1 || block==K()?q:GeometryPoint(block,q.geometry,q.coordinate);
            const double expected=b==4?level-(outer?referenceLaw.quadraticStress[2][2]:0.):data.boundary[b](q.coordinate,time);
            const RANSWallCoordinate gradient{0.,0.,b==0?du:b==K()?dk:b==3?-A/((z-d)*(z-d)):0.};
            const double delta=row.Value(state,offset[block])-expected;const auto actualGradient=row.Gradient(state,offset[block]);
            square[b]+=q.weight*delta*delta;volume[b]+=q.weight;largest[b]=std::max(largest[b],std::abs(delta));
            for(int axis=0;axis<3;++axis)gradientSquare[b]+=q.weight*(actualGradient[axis]-gradient[axis])*(actualGradient[axis]-gradient[axis]);
        }
    }
    std::vector<std::string> errors;for(int b=0;b<5;++b)errors.push_back(J::Object({{"block",J::Number(b)},
        {"weighted_RMS",J::Number(std::sqrt(square[b]/volume[b]))},{"true_L2",J::Number(std::sqrt(square[b]))},
        {"gradient_weighted_RMS",J::Number(std::sqrt(gradientSquare[b]/volume[b]))},{"maximum",J::Number(largest[b])},
        {"scale",J::Number(scales[b])},{"physical_volume",J::Number(volume[b])}}));
    return J::Object({{"schema",J::Quote("AcceptedSerialSparsePacketV1")},{"read_only","true"},{"alpha",J::Number(alpha)},
        {"identity",J::Quote(PreservationIdentity())},{"stored_slots",J::Number(Size())},{"masses",J::List(masses)},
        {"flow_rows",J::Number(flow.rows.size())},{"flow_nnz",J::Number(flow.csr.NonZeros())},
        {"velocity_rows",J::Number(nv)},{"pressure_rows",J::Number(r.dz.size())},{"pressure_C_nnz",J::Number(r.c.NonZeros())},
        {"predictor_factor_nnz",J::Number(r.predictorFill)},{"consistent_mass_factor_nnz",J::Number(r.massFill)},
        {"pressure_preconditioner_factor_nnz",J::Number(r.preconditionerFill)},{"pressure_iterations",J::Number(r.iterations)},
        {"mass_applications",J::Number(r.massApplications)},{"predictor_residual",J::Number(Maximum(r.predictorDefect))},
        {"predictor_rhs_maximum",J::Number(r.predictorRhsMaximum)},
        {"pressure_rhs_maximum",J::Number(r.pressureRhsMaximum)},
        {"full_native_G_action_defect",J::Number(Maximum(actionDefect))},
        {"consistent_mass_correction_rhs_maximum",J::Number(Maximum(r.Gz))},
        {"consistent_mass_correction_true_residual",J::Number(Maximum(massDefect))},
        {"true_pressure_residual",J::Number(Maximum(r.pressureDefect))},{"true_continuity",J::Number(Maximum(r.constraints))},
        {"splitting_identity",J::Number(Maximum(r.identity))},{"predictor",J::Numbers(r.predictor)},
        {"pressure_increment",J::Numbers(r.dz)},{"consistent_mass_correction",J::Numbers(r.bGz)},{"scalars",J::List(scalars)},
        {"continuous_reference_errors",J::List(errors)},
        {"scope",J::Quote("frozen initial-state diagnostic linear operations only; no evolution or extra governing sweep; consistent mass in action, diagonal mass only in named preconditioner; no dense Schur or M-inverse-G columns")}});
}
}
