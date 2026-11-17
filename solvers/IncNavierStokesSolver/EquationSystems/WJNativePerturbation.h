// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_PERTURBATION_H
#define NEKTAR_WJ_NATIVE_PERTURBATION_H

#include <IncNavierStokesSolver/EquationSystems/WJNativeFlow.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeReference.h>
#include <array>
#include <string>
#include <vector>

namespace Nektar::WJNative
{
enum class PerturbationShape { PrimaryQuartile=0, SecondaryLowerSkew=1 };

/** Select once on the admitted base mesh, serialize, and reuse unchanged on
 * refinements. Coordinates are x and wall distance r, not physical z.
 */
struct PerturbationLayout
{
    std::array<double,5> xKnots{}, rKnots{};
    double H=.15, epsilonPeakRadius=0;
    PerturbationShape shape=PerturbationShape::PrimaryQuartile;
};

struct PerturbationSigns
{
    int velocity=1, k=1, epsilon=1;
};

/** External immutable Linux evidence for the GLOBAL peak of UnitVelocity.
 * This container is not itself a proof. The qualification runner must bind
 * the supplied SHA256 to the actual record, polynomial peak bounds and mesh.
 */
struct PerturbationNormalization
{
    PerturbationLayout layout;
    double unitPeakLower=0, unitPeakEstimate=0, unitPeakUpper=0;
    std::string evidenceSha256;
};

struct PerturbationScalarJet
{
    double value=0, dx=0, dr=0;
};

/** Prospective finite physical perturbation shapes; no engine/state mutation.
 * The cubic B-spline streamfunction gives exact native-Q3-representable,
 * compact, divergence-free velocity on a supported flat tensor mesh. Scalar k
 * is also piecewise Q3. The independent epsilon TARGET includes Ain/r and
 * requires separate qualification after native Q3 projection.
 *
 * Native edge vectors must describe common tensor-product layers across the
 * support. The runner must verify that condition from actual element topology;
 * a union of unrelated local edge coordinates is not sufficient evidence.
 */
class Perturbation
{
public:
    static constexpr const char *Identity="WJH27FinitePerturbationB3V1";
    static constexpr double VelocityPeak=1e-3;
    static constexpr double KRelativePeak=.01;
    static constexpr double EpsilonRelativePeak=.01;
    static constexpr double PeakBracketRelativeWidth=1e-8;

    static PerturbationLayout SelectLayout(
        const std::vector<double> &nativeXEdges,
        const std::vector<double> &nativeREdges, double H,
        PerturbationShape shape=PerturbationShape::PrimaryQuartile);

    Perturbation(const Reference &reference, const PerturbationLayout &layout,
                 const std::vector<double> &nativeXEdges,
                 const std::vector<double> &nativeREdges,
                 PerturbationSigns signs={});

    const PerturbationLayout &Layout() const { return m_layout; }
    PerturbationSigns Signs() const { return m_signs; }
    int StreamwiseIntervals() const { return m_xIntervals; }
    int NormalIntervals() const { return m_rIntervals; }
    int OuterNormalIntervals() const { return m_outerIntervals; }
    double SupportArea() const;
    double OuterSupportArea() const;
    double DomainArea() const { return m_domainArea; }
    double OuterDomainArea() const { return m_outerDomainArea; }
    double KIn() const { return m_KIn; }
    double AIn() const { return m_AIn; }
    double Displacement() const { return m_d0; }

    /** Load exactly once during initial-data construction, before histories or
     * pilot evolution. No guessed constant or subsequent amplitude tuning.
     */
    void SetNormalization(const PerturbationNormalization &normalization);
    bool HasNormalization() const { return m_hasNormalization; }
    const PerturbationNormalization &Normalization() const;
    double StreamfunctionCoefficient() const;

    // Public evaluation coordinates are physical x,z. Unit methods use
    // unsigned psi=B3_x*B3_r and require no peak normalization.
    double UnitStreamfunction(double x, double z) const;
    Vector UnitVelocity(double x, double z) const;
    Tensor UnitGradient(double x, double z) const;
    Vector Velocity(double x, double z) const;
    Tensor VelocityGradient(double x, double z) const;
    // Unit scalar shapes have peak one. The independent epsilon unit shape
    // is the reference-scaled bump, before multiplication by Ain/r.
    PerturbationScalarJet KShape(double x, double z) const;
    PerturbationScalarJet IndependentEpsilonShape(double x, double z) const;
    PerturbationScalarJet KIncrement(double x, double z) const;
    PerturbationScalarJet IndependentEpsilonIncrement(double x, double z) const;

private:
    double Radius(double x, double z) const;
    PerturbationLayout m_layout;
    PerturbationSigns m_signs;
    PerturbationNormalization m_normalization;
    bool m_hasNormalization=false;
    double m_KIn=0, m_AIn=0, m_d0=0, m_domainArea=0, m_outerDomainArea=0;
    int m_xIntervals=0, m_rIntervals=0, m_outerIntervals=0;
};
}
#endif
