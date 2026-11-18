// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_REFERENCE_H
#define NEKTAR_WJ_NATIVE_REFERENCE_H

#include <array>
#include <string>
#include <vector>

namespace Nektar::WJNative
{
struct ReferenceParameters
{
    double KIn, H, roughness, nu, Cd, sigmaK, sigmaE;
    double qe, Ne, kappa, cLength, d0, topRadius, tau0, uTau, z0;
    double KMinimum, KMaximum;
};

struct ReferenceValue
{
    double value, dK, dr;
};

enum class InnerKTransport
{
    ScalarV2,
    AuxiliaryHistoricalFixtureV1
};

enum class ReferencePolicy
{
    Historical,
    ExactInnerFiniteNuV1
};

/** Immutable, qualified H27-v2 counterfactual tables and prescribed inlet.
 * All heights are MODEL radius r=physical_z-d0. K is the current sampled
 * interface value, restricted to the qualified interval. These tables never
 * impose a zero actual-flow production/destruction budget.
 */
class Reference
{
public:
    static constexpr const char *QualifiedTableSha256 =
        "b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d";
    static constexpr const char *QualifiedExactInnerSha256 =
        "3eeb32e6c1f53d6045d9d03707344c65f0f1ec3964dc69116a87e04e104325e1";
    static ReferencePolicy ParsePolicy(const std::string &name);
    Reference(const std::string &json, ReferencePolicy policy);
    explicit Reference(const std::string &json,
        const std::string &expectedSha256 = QualifiedTableSha256,
        InnerKTransport transport = InnerKTransport::ScalarV2);
    const ReferenceParameters &Parameters() const { return m_p; }
    const ReferenceParameters &GetParameters() const { return m_p; }
    const std::string &Identity() const { return m_identity; }
    const std::string &BoundaryIdentity() const { return m_boundaryIdentity; }
    ReferencePolicy Policy() const { return m_policy; }
    const char *PolicyName() const;
    InnerKTransport Transport() const { return m_transport; }
    const char *TransportName() const;

    double Cutoff(double K) const;
    double CutoffDerivative(double K) const;
    std::array<double,2> LogCutoffDerivativeBounds() const { return m_logCutoffSlopeBounds; }
    std::array<double,2> CutoffRadiusRange() const;
    // Geometry-only inversion of the unchanged certified decreasing table.
    // This selects x integration panels, never a new cutoff calibration.
    double KAtCutoff(double radius) const;
    double Amplitude(double K) const;
    double AmplitudeDerivative(double K) const;
    double AmplitudeSecondDerivative(double K) const;
    ReferenceValue Auxiliary(double K, double r) const;
    ReferenceValue DependentEpsilon(double K, double r) const;
    double Aux(double K, double r) const { return Auxiliary(K,r).value; }
    double AuxDerivativeK(double K, double r) const { return Auxiliary(K,r).dK; }
    double AuxDerivativeR(double K, double r) const { return Auxiliary(K,r).dr; }
    double Ew(double K, double r) const { return DependentEpsilon(K,r).value; }
    double EwDerivativeK(double K, double r) const { return DependentEpsilon(K,r).dK; }
    double EwDerivativeR(double K, double r) const { return DependentEpsilon(K,r).dr; }
    double KDiffusivity(double k, double epsilon, double auxiliary, bool outer) const;

    // Prescribed inner profile and fixed outer logarithmic U and constant k.
    // InletE belongs to the transported outer field. The inner destruction is
    // always evaluated with Ew(current K,r).
    double InletU(double r) const;
    double InletK(double r) const;
    double InletE(double r) const;
    double InletUGradient(double r) const;
    double InletKGradient(double r) const;
    double InletEGradient(double r) const;
    double PrescribedInletCutoff() const { return m_inletCutoff; }
    double InletModifiedPressure(double r) const;
    double TopTangentialTraction() const;
    std::array<double,2> OutletTraction(double r) const; // components x,z
    double TopEpsilonOutwardFlux() const; // -D_epsilon * d epsilon_ref/dr

private:
    struct Polynomial
    {
        double lower, upper;
        std::vector<double> coefficients;
        double Evaluate(double argument) const;
        Polynomial Differentiate() const;
    };
    struct InletPanel
    {
        Polynomial u, k, pressure, du, dk;
    };
    Reference(const std::string &json, const std::string &expectedSha256,
              InnerKTransport transport, ReferencePolicy policy);
    const InletPanel &InnerPanel(double r) const;
    void CheckK(double K) const;
    void CheckRadius(double r) const;
    std::array<double,2> HistoricalLogCutoff(double K) const;
    void BindFixtureTransport();
    double InletQ(double r) const;
    ReferenceParameters m_p{};
    std::string m_identity, m_boundaryIdentity;
    ReferencePolicy m_policy;
    InnerKTransport m_transport = InnerKTransport::ScalarV2;
    Polynomial m_logCutoff, m_logCutoffDerivative;
    std::array<double,2> m_logCutoffSlopeBounds{};
    std::array<Polynomial,3> m_outlet;
    std::vector<InletPanel> m_exactInnerInlet;
    double m_referenceShear = 0.;
    double m_inletA=0, m_inletSlope=0, m_inletD=0, m_inletCutoff=0;
};
}
#endif
