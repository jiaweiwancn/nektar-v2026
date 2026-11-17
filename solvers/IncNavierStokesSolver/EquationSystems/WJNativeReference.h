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

/** Immutable, qualified H27-v2 counterfactual tables and historical inlet.
 * All heights are MODEL radius r=physical_z-d0. K is the current sampled
 * interface value, restricted to the qualified interval. These tables never
 * impose a zero actual-flow production/destruction budget.
 */
class Reference
{
public:
    static constexpr const char *QualifiedTableSha256 =
        "b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d";
    explicit Reference(const std::string &json,
        const std::string &expectedSha256 = QualifiedTableSha256);
    const ReferenceParameters &Parameters() const { return m_p; }
    const ReferenceParameters &GetParameters() const { return m_p; }
    const std::string &Identity() const { return m_identity; }

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

    // Fixed historical inlet: old finite-nu momentum and nonconstant k below
    // H, high-Re-derived logarithmic U and constant k above H. InletE is only
    // the transported OUTER epsilon; inner destruction is Ew(current K,r).
    double InletU(double r) const;
    double InletK(double r) const;
    double InletE(double r) const;
    double InletUGradient(double r) const;
    double InletKGradient(double r) const;
    double InletEGradient(double r) const;
    double HistoricalInletCutoff() const { return m_inletCutoff; }
    std::array<double,2> OutletTraction(double r) const; // components x,z
    double TopEpsilonOutwardFlux() const; // -D_epsilon * d epsilon_ref/dr

private:
    struct Polynomial
    {
        double lower, upper;
        std::vector<double> coefficients;
        double Evaluate(double argument) const;
    };
    void CheckK(double K) const;
    void CheckRadius(double r) const;
    double InletQ(double r) const;
    ReferenceParameters m_p{};
    std::string m_identity;
    Polynomial m_logCutoff, m_logCutoffDerivative;
    std::array<double,2> m_logCutoffSlopeBounds{};
    std::array<Polynomial,3> m_outlet;
    double m_inletA=0, m_inletSlope=0, m_inletD=0, m_inletCutoff=0;
};
}
#endif
