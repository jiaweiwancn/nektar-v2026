// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_MIXED_FLOW_H
#define NEKTAR_WJ_NATIVE_MIXED_FLOW_H
#include <IncNavierStokesSolver/EquationSystems/WJNativeFlow.h>

namespace Nektar::WJNative
{
/** H27 opt-in component-boundary mixed solver. Full native actions use the
 * u signed basis as a canonical vector basis. The actual v essential mask
 * and its signed native permutation are retained separately. Historical
 * H26 Flow and its constructor/default behavior are unchanged.
 */
class MixedFlow
{
public:
    MixedFlow(const std::vector<MultiRegions::ContFieldSharedPtr>&fields,
              Geometry geometry);
    int VelocitySize()const{return m_core.VelocitySize();}
    int PressureSize()const{return m_core.PressureSize();}
    int PressureOffset()const{return int(m_freeFull.size());}
    int Size()const{return PressureOffset()+PressureSize()+(m_gauge?1:0);}
    bool Gauge()const{return m_gauge;}
    double ConstantPressureAction()const{return m_constantPressureAction;}
    const Flow&Core()const{return m_core;}
    const std::array<std::vector<int>,2>&FreeVelocityIds()const{return m_free;}
    const std::vector<int>&NativeVToCanonical()const{return m_vToU;}
    const std::vector<double>&NativeVSigns()const{return m_vSign;}
    const std::vector<Element>&Elements()const{return m_core.Elements();}
    const Vec&PressureOne()const{return m_core.PressureOne();}
    const Vec&PressureMean()const{return m_core.PressureMean();}
    Sparse System(double alpha,const Sparse&K,double awf)const;
    Vec Project(const Function&f,int component)const;
    Vec NativeVelocity(const Vec&full,int component)const;
    Vec Lift(const Vec&full)const;
    Vec Restrict(const Vec&full)const;
    Vec Expand(const Vec&free,const Vec&lift)const;
    Vec LiftedRhs(double alpha,const Sparse&K,double awf,
                  const Vec&load,const Vec&lift)const;
    SolveInfo Solve(const Sparse&current,const Vec&rhs,Vec&x,
                    double coefficientChange=0,double tolerance=1e-11);
    SolveInfo Residual(const Sparse&current,const Vec&rhs,const Vec&x)const;
    /** Fixed, homogeneous norms of the raw free momentum/divergence/mean rows.
     * No time-mass or state-dependent denominator enters this action. */
    SolveInfo RawResidualNorm(const Vec&residual)const;
    const Vec&ResidualRowScale()const{return m_scale;}
    double PrepareLinearPreconditioner(const Sparse&current,bool amdOrdering=false)
    {return Prepare(current,amdOrdering);}
    Vec ApplyLinearPreconditioner(const Vec&rhs)const{return Precondition(rhs);}
    void ResetPreconditioner()
    {m_factor.reset();m_amdFactor.reset();m_initialIterations=m_lastIterations=0;}
    double VelocityDual(const Vec&load)const;
    void SetTransportedKSpace(const MultiRegions::ContFieldSharedPtr&field);
    Vec ProjectTransportedK(const Function&f)const;
    Vec KnownScalarKLoad(const Vec&nativeGlobalK)const;
    const Sparse&KnownScalarAction()const{return m_scalarAction;}
    const MultiRegions::ContFieldSharedPtr&TransportedKField()const{return m_kField;}
private:
    Vec ProjectNative(const Function&f,const MultiRegions::ContFieldSharedPtr&field)const;
    double Prepare(const Sparse&current,bool amdOrdering=false);
    Vec Precondition(const Vec&rhs)const;
    std::vector<MultiRegions::ContFieldSharedPtr>m_fields;
    Flow m_core;
    std::array<std::vector<int>,2>m_free;
    std::vector<int>m_freeFull,m_fullToFree,m_vToU;
    std::vector<double>m_vSign;
    bool m_gauge=false,m_sharedRiesz=false;
    double m_volume=0,m_constantPressureAction=0;
    int m_pressurePivot=0,m_initialIterations=0,m_lastIterations=0;
    Vec m_scale;
    std::array<std::unique_ptr<Eigen::SimplicialLDLT<Sparse>>,2>m_velocityRiesz;
    Eigen::SimplicialLDLT<Sparse>m_pressureRiesz;
    std::unique_ptr<Eigen::SparseLU<Sparse>>m_factor;
    std::unique_ptr<Eigen::SparseLU<Sparse,Eigen::AMDOrdering<int>>>m_amdFactor;
    MultiRegions::ContFieldSharedPtr m_kField;
    Sparse m_scalarAction;
};
}
#endif
