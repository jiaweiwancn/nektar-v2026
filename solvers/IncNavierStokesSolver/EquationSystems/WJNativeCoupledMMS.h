// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_COUPLED_MMS_H
#define NEKTAR_WJ_NATIVE_COUPLED_MMS_H

#include <IncNavierStokesSolver/EquationSystems/WJNativeFlow.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeReference.h>

namespace Nektar::WJNative
{
enum class MMSVariant { Historical, PolynomialScalarsV1 };

/** Continuous derivatives in physical x and flat wall distance r=z-d0.
 * dr is also the physical z derivative. These are analytic field derivatives,
 * not derivatives of projected native coefficients.
 */
struct MMSJet
{
    double value=0, dx=0, dr=0, dt=0, dxx=0, dxr=0, drr=0;
};

struct MMSFields
{
    MMSJet u, w, Pi, k;
    Tensor gradient{}, gradientX{}, gradientR{};
};

struct MMSForcing
{
    Vector momentum{}, awf{}, body{};
    double k=0, epsilon=0, production=0, destruction=0;
    double diffusionK=0, diffusionE=0;
    bool outer=false;
};

/** VERIFICATION ONLY: evolving continuous coupled MMS for the amended v2 law.
 * Domain x in [0,1], r in [0,.5], flat geometry, physical nu in both layers.
 * This helper never supplies application initial/inlet data or model balance
 * forcing. Its manufactured momentum load is DISTINCT from the physical body
 * force. The historical default bodyAmplitude=0 remains identically zero;
 * the optional analytic nonzero-body fixture exercises bulk and sampled-H
 * contributions separately, without changing the manufactured fields.
 * Constitutive admissibility remains a separate required native screen.
 */
class CoupledMMS
{
public:
    static constexpr const char *Identity = "WJH27ContinuousCoupledMMSV1";
    static constexpr const char *BodyIdentity = "WJH27PrescribedPhysicalBodyMMSV1";
    CoupledMMS(const Reference &reference, Geometry geometry,
               double topRadius=.5, double bodyAmplitude=0,
               MMSVariant variant=MMSVariant::Historical);

    Geometry Mapping() const { return m_geometry; }
    double TopRadius() const { return m_top; }
    double BodyAmplitude() const { return m_bodyAmplitude; }
    MMSVariant Variant() const { return m_variant; }
    const char *VariantName() const;
    MMSFields Fields(double x, double z, double time) const;
    /** Tangential/time jet of K(x,t)=k(x,H,t); r derivatives are zero. */
    MMSJet InterfaceK(double x, double time) const;
    /** Transported epsilon exists only outside H. Includes its outer Hessian.
     * The dependent inner Ew is supplied by Reference and has a cutoff kink;
     * it is not represented as a fictitious transported smooth scalar here.
     */
    MMSJet OuterEpsilon(double x, double z, double time) const;
    MMSForcing Forcing(double x, double z, double time) const;

    /** Unsuppressed H-sampled AWF amplitude; Forcing applies 1_(r<=h).
     * The return vector has components x,y,z; y and z are zero.
     */
    Vector AWFAtH(double x, double time) const;
    /** Prescribed physical body, independent of current/projected state:
     * a*(2e-4*(1+x+.5*r)*cos(.3*t), 0, 3e-4*(1+r)*sin(.2*t)).
     * Forcing.momentum is the verification load L(U,Pi)-Body-AWF;
     * it must be combined with bulk Body and the model's sampled body term.
     */
    Vector Body(double x, double z, double time) const;
    Vector BodyAtH(double x, double time) const;
    /** Exact outward diffusive flux -D grad(q).normal, never a gradient BC.
     * Use k on natural right/top/ground and epsilon on natural right/top;
     * inlet and outer epsilon's H trace remain essential conditions.
     */
    double KOutwardFlux(double x, double z, double time,
                        const Vector &normal) const;
    double EpsilonOutwardFlux(double x, double z, double time,
                              const Vector &normal) const;

private:
    double Radius(double x, double z, double time) const;
    const Reference &m_reference;
    Geometry m_geometry;
    double m_top, m_bodyAmplitude;
    MMSVariant m_variant;
};
}
#endif
