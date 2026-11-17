///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Opt-in H26 prescribed-data flow kernel. No transported turbulence state.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_WJ_NATIVE_FLOW_H
#define NEKTAR_WJ_NATIVE_FLOW_H
#include <MultiRegions/ContField.h>
#include <SolverUtils/Turbulence/WJEARSM.h>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include <functional>
#include <memory>
#include <vector>

namespace Nektar::WJNative
{
using Vec=Eigen::VectorXd;
using Mat=Eigen::MatrixXd;
using Sparse=Eigen::SparseMatrix<double>;
using Tensor=SolverUtils::RANSTensor;
using Vector=SolverUtils::RANSVector;
using Function=std::function<double(double,double)>;
struct Geometry
{
    double d0=0,H=.15,h=.10854503542535461,curve=0;
    double Ground(double x) const {return d0+curve*x*x*(1-x)*(1-x);}
    double Slope(double x) const {return curve*(2*x-6*x*x+4*x*x*x);}
};
struct WallStress
{
    Tensor R{},V{},remainder{};
    double mu=0,chi=0,production=0;
};
WallStress Adapt(double k,double epsilon,const Tensor &G,const Vector &normal,
                 double nu,double nuA,double r,double H);
struct SolveInfo
{
    int iterations=0,rebuilds=0;
    double momentum=0,continuity=0,gauge=0,correction=0;
    double factorSeconds=0,solveSeconds=0,residualSeconds=0;
};
struct Element
{
    int geometry=0;
    double xl=0,xr=0,rl=0,rr=0;
    std::vector<int> u,p;
    std::vector<double> us,ps;
    Mat B,Dx,Dz,P,Px,Pz;
    Vec x,z,w;
    LocalRegions::ExpansionSharedPtr velocity,pressure;
};
/** Native continuous mixed Galerkin blocks. Native quadrature/physical
 * derivatives and signed assembly maps define every operator. Only sparse
 * global matrices are used; full dense saddle solves belong to tiny oracles.
 */
class Flow
{
public:
    Flow(const std::vector<MultiRegions::ContFieldSharedPtr>&fields,Geometry geometry);
    int VelocitySize()const{return m_nu;}
    int PressureSize()const{return m_np;}
    int DirichletSize()const{return m_nd;}
    int FreeVelocitySize()const{return m_nf;}
    int Size()const{return 2*m_nf+m_np+(m_gauge?1:0);}
    bool Gauge()const{return m_gauge;}
    const std::vector<Element>&Elements()const{return m_elements;}
    const Sparse&Mass()const{return m_mass;}
    const Sparse&PressureMass()const{return m_qmass;}
    const Sparse&Divergence()const{return m_div;}
    const Sparse&Gradient()const{return m_grad;}
    const Sparse&AWF()const{return m_awf;}
    const Vec&PressureMean()const{return m_mean;}
    Geometry Mapping()const{return m_geometry;}
    Sparse Strain(const Function &mu,bool transpose=true)const;
    Sparse StrainSamples(const std::vector<Vec>&mu,bool transpose=true)const;
    Sparse System(double alpha,const Sparse &K,double awf)const;
    Vec Load(const std::array<Function,2>&f)const;
    Vec StressLoad(const std::function<Tensor(double,double)> &T)const;
    Vec StressSamples(const std::vector<std::vector<Tensor>>&T)const;
    Vec AWFLoad(const std::function<Vector(double,double)>&sampleAtH)const;
    Vec KnownKLoad(const Vec &prescribedK)const;
    Vec MassAction(const Vec &full)const;
    Vec Convection(const Vec &full)const;
    double PreparePreconditioner(const Sparse &a);
    Vec TractionLoad(int boundary,const std::function<Vector(double,double,const Vector&)>&t)const;
    Vec Project(const Function &f,int component)const;
    Vec Lift(const Vec &full)const;
    Vec Restrict(const Vec &full)const;
    Vec Expand(const Vec &free,const Vec &lift)const;
    Vec LiftedRhs(double alpha,const Sparse &K,double awf,const Vec &load,const Vec &lift)const;
    Vec ApplyFull(double alpha,const Sparse &K,double awf,const Vec &full)const;
    SolveInfo Solve(const Sparse &current,const Vec &rhs,Vec &x,double relativeChange=0,double tolerance=1e-11);
    SolveInfo Residual(const Sparse &current,const Vec &rhs,const Vec &x)const;
    std::array<double,3> Error(const Vec &full,const std::array<Function,3>&exact)const;
    double StrongDivergence(const Vec &full)const;
    double VelocityDual(const Vec &load)const;
    double PressureNorm(const Vec &pressure)const;
    Vec Physical(const Vec &full,int component)const;
private:
    void BuildInvariant();
    void BuildAWF();
    Vec RowScale()const;
    Vec Precondition(const Vec&rhs)const;
    std::vector<MultiRegions::ContFieldSharedPtr>m_fields;
    std::vector<Element>m_elements;
    Geometry m_geometry;
    int m_nu=0,m_np=0,m_nd=0,m_nf=0;
    bool m_gauge=true;
    double m_volume=0,m_length=1,m_forceScale=1;
    Sparse m_mass,m_qmass,m_div,m_grad,m_awf,m_freeMass;
    Vec m_mean,m_scale,m_pressureOne;
    int m_pressurePivot=0;
    Eigen::SimplicialLDLT<Sparse>m_velocityRiesz,m_pressureRiesz;
    std::unique_ptr<Eigen::SparseLU<Sparse>>m_factor;
    int m_initialIterations=0;
    int m_lastIterations=0;
};
}
#endif
