///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <numeric>

namespace Nektar
{
using namespace RANSWallAlgebra;
namespace J=RANSWallSolveJson;

RANSWallEvolutionPatch::SparseFlowResult RANSWallEvolutionPatch::SparseFlowCorrection(
    const PicardBlock &flow,int nv,double alpha,const std::function<void(const std::string&)> &event) const
{
    RANSWallTiming::Scope timing("sparse.flow_linear_operation");
    Require(SerialSparse() && flow.sparse && alpha>0,"invalid sparse flow operation");
    const int n=flow.rows.size(),np=n-nv;SparseFlowResult r;
    r.a=flow.csr.Block(0,nv,0,nv);r.g=flow.csr.Block(0,nv,nv,n);
    r.d=flow.csr.Block(nv,n,0,nv);r.c=flow.csr.Block(nv,n,nv,n);
    r.m=SparseMatrix(nv,nv);std::vector<int> index(Size(),-1);
    for(int i=0;i<nv;++i)index[flow.rows[i]]=i;
    for(int i=0;i<nv;++i)
    {
        const int b=flow.rows[i]/size[0],local=flow.rows[i]-offset[b];
        m_sparseMass[b].Each(local,[&](int j,double v){const int k=index[offset[b]+j];if(k>=0)r.m.Entry(i,k)=v;});
    }
    if(event)event("momentum_attempt");
    r.m.Finalize();SparseFactor momentum(r.a),mass(r.m,true);
    r.predictorFill=momentum.FactorNonZeros();r.massFill=mass.FactorNonZeros();
    Vector rhs(nv);for(int i=0;i<nv;++i)rhs[i]=-flow.residual[i];
    r.predictor=momentum.Apply(rhs);r.predictorDefect=r.a.Apply(r.predictor);
    if(event)event("momentum_complete");
    for(int i=0;i<nv;++i)r.predictorDefect[i]-=rhs[i];
    auto pressureRhs=r.d.Apply(r.predictor);
    for(int i=0;i<np;++i)pressureRhs[i]+=flow.residual[nv+i];
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

std::string RANSWallEvolutionPatch::SparseBackendPacket(const Vector &state,double time,
    const RANSWallEvolutionData &data,double alpha) const
{
    Require(SerialSparse(),"sparse backend packet requires its explicit policy");
    auto external=ExternalConstraints(time,data);for(int id:m_traceIds)external.erase(offset[Epsilon()]+id);
    const auto flow=AssemblePicardBlock(0,state,state,alpha,time,data,external);
    const int nv=std::lower_bound(flow.rows.begin(),flow.rows.end(),DifferentialSize())-flow.rows.begin();
    const auto r=SparseFlowCorrection(flow,nv,alpha);
    std::vector<std::string> masses,scalars;
    for(int b=0;b<Gauge();++b)masses.push_back(J::Object({{"rows",J::Number(size[b])},{"nnz",J::Number(m_sparseMass[b].NonZeros())}}));
    for(int b:{Epsilon(),K()})
    {
        const auto block=AssemblePicardBlock(b,state,state,alpha,time,data,external);SparseFactor factor(block.csr);
        auto rhs=block.residual;for(auto &x:rhs)x=-x;const auto delta=factor.Apply(rhs);
        auto error=block.Apply(delta);for(int i=0;i<error.size();++i)error[i]-=rhs[i];
        scalars.push_back(J::Object({{"block",J::Number(b)},{"rows",J::Number(block.rows.size())},{"nnz",J::Number(block.csr.NonZeros())},
            {"factor_nnz",J::Number(factor.FactorNonZeros())},{"maximum_true_residual",J::Number(Maximum(error))},{"diagnostic_increment",J::Numbers(delta)}}));
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
        {"true_pressure_residual",J::Number(Maximum(r.pressureDefect))},{"true_continuity",J::Number(Maximum(r.constraints))},
        {"splitting_identity",J::Number(Maximum(r.identity))},{"predictor",J::Numbers(r.predictor)},
        {"pressure_increment",J::Numbers(r.dz)},{"consistent_mass_correction",J::Numbers(r.bGz)},{"scalars",J::List(scalars)},
        {"continuous_reference_errors",J::List(errors)},
        {"scope",J::Quote("frozen initial-state diagnostic linear operations only; no evolution or extra governing sweep; consistent mass in action, diagonal mass only in named preconditioner; no dense Schur or M-inverse-G columns")}});
}
}
