///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_REDUCED_MODEL_H
#define NEKTAR_RANS_WALL_REDUCED_MODEL_H
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Nektar::RANSWallReduced
{
using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;
inline double Dot(const Vector &a, const Vector &b)
{
    if (a.size() != b.size())
        throw std::logic_error("reduced model dot product shape mismatch");
    long double value = 0.;
    for (std::size_t i = 0; i < a.size(); ++i)
        value += static_cast<long double>(a[i]) * b[i];
    return static_cast<double>(value);
}
inline double Length(const Vector &a)
{
    return std::sqrt(std::max(0., Dot(a, a)));
}
// Relative progress remains meaningful near the original qualification
// boundary. An absolute O(1) scale can reject a large fractional restoration
// when the remaining normalized violation is small.
inline bool QualificationDecrease(double before, double after,
                                  double fraction, double rate)
{
    if (!std::isfinite(before) || !std::isfinite(after) || before < 0. || after < 0. ||
        !(fraction > 0. && fraction <= 1.) || !(rate > 0. && rate < 1.))
        throw std::logic_error("invalid qualification progress comparison");
    return after < before * (1. - rate * fraction);
}
struct QuadraticResult
{
    bool solved = false;
    Vector solution, multipliers;
    std::vector<int> active;
    int steps = 0, redundantRows = 0;
    double primal = 0., stationarity = 0., complementarity = 0.;
    std::string failure;
};

// Broyden update in the actual field metric. This is a finite-step response
// model, not a claim of an exact Jacobian. The caller must still evaluate its
// nonlinear trial. Leave changes below the declared signal scale untouched.
inline bool UpdateSecantActions(Matrix &actions, const Matrix &columns,
                                const Vector &coefficients,
                                const Vector &observedChange, double minimumNorm)
{
    const int n = columns.size();
    if (!n || actions.size() != columns.size() || coefficients.size() != columns.size() ||
        !(minimumNorm > 0.) || !std::isfinite(minimumNorm))
        throw std::logic_error("invalid reduced response secant shape or signal scale");
    Vector field(columns[0].size()), defect = observedChange, covector(n);
    for (int i = 0; i < n; ++i)
    {
        if (columns[i].size() != field.size() || actions[i].size() != defect.size())
            throw std::logic_error("inconsistent reduced response secant rows");
        for (int j = 0; j < field.size(); ++j) field[j] += coefficients[i]*columns[i][j];
        for (int j = 0; j < defect.size(); ++j) defect[j] -= coefficients[i]*actions[i][j];
    }
    const double denominator = Dot(field,field);
    if (!std::isfinite(denominator))
        throw std::runtime_error("nonfinite reduced response secant signal");
    if (denominator <= minimumNorm*minimumNorm) return false;
    for (int i = 0; i < n; ++i) covector[i] = Dot(columns[i],field)/denominator;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < defect.size(); ++j)
        {
            actions[i][j] += defect[j]*covector[i];
            if (!std::isfinite(actions[i][j]))
                throw std::runtime_error("nonfinite reduced response secant update");
        }
    return true;
}

// Small strictly convex subspace problem, starting at a feasible point.
// The factor type supplies the native dense factorization in production;
// this kernel has no field assembly or governing-equation solve.
template <typename Factor, typename Cut>
QuadraticResult SolveQuadratic(const Matrix &hessian, const Vector &gradient,
                              const std::vector<Cut> &cuts, Vector initial,
                              int maximumSteps = 120)
{
    const int n = gradient.size();
    if (!n || initial.size() != gradient.size() ||
        hessian.size() != gradient.size() || maximumSteps <= 0)
        throw std::logic_error("invalid reduced quadratic model shape");
    for (const auto &row : hessian)
        if (row.size() != gradient.size())
            throw std::logic_error("nonsquare reduced quadratic model");
    for (const auto &cut : cuts)
        if (cut.a.size() != gradient.size())
            throw std::logic_error("invalid reduced quadratic cut shape");
    QuadraticResult result;
    result.solution = std::move(initial);
    auto &solution = result.solution;
    auto &active = result.active;
    auto &multipliers = result.multipliers;
    auto finish = [&](const std::string &failure) {
        result.failure = failure;
        Vector g = gradient;
        for (int i = 0; i < n; ++i)
            g[i] += Dot(hessian[i], solution);
        for (int i = 0; i < active.size(); ++i)
        {
            const double lambda = i < multipliers.size() ? multipliers[i] : 0.;
            for (int j = 0; j < n; ++j)
                g[j] += lambda * cuts[active[i]].a[j];
            result.complementarity = std::max(result.complementarity,
                std::abs(lambda * (Dot(cuts[active[i]].a, solution) - cuts[active[i]].b)));
        }
        result.stationarity = Length(g);
        for (const auto &cut : cuts)
            result.primal = std::max(result.primal, Dot(cut.a, solution) - cut.b);
        result.solved = failure.empty() && result.primal <= 1e-8;
        if (failure.empty() && !result.solved)
            result.failure = "quadratic_primal_defect";
        return result;
    };
    for (const auto &cut : cuts)
        if (Dot(cut.a, solution) > cut.b + 1e-8)
            return finish("quadratic_infeasible_start");
    try
    {
        const Factor objective(hessian);
        Vector initialGradient(n);
        for (int i = 0; i < n; ++i)
            initialGradient[i] = -gradient[i] - Dot(hessian[i], solution);
        const double stepScale = std::max(Length(solution),
                                         Length(objective.Apply(initialGradient)));
        while (result.steps < maximumSteps)
        {
            ++result.steps;
            Vector g(n);
            for (int i = 0; i < n; ++i)
                g[i] = -gradient[i] - Dot(hessian[i], solution);
            // A_active^T = Q R. Work in its Euclidean null space rather
            // than forming A H^-1 A^T, which squares active-row conditioning.
            std::vector<Vector> face;
            Matrix upper(active.size(), Vector(active.size()));
            for (int i = 0; i < active.size(); ++i)
            {
                auto row = cuts[active[i]].a;
                for (int pass = 0; pass < 2; ++pass)
                    for (int j = 0; j < face.size(); ++j)
                    {
                        const double c = Dot(row, face[j]);
                        upper[j][i] += c;
                        for (int k = 0; k < n; ++k) row[k] -= c * face[j][k];
                    }
                const double length = Length(row);
                if (length < 1e-10 * Length(cuts[active[i]].a))
                    return finish("quadratic_active_rank_loss");
                upper[i][i] = length;
                for (double &x : row) x /= length;
                face.push_back(std::move(row));
            }
            Vector direction(n);
            if (active.empty())
                direction = objective.Apply(g);
            else
            {
                std::vector<Vector> null;
                for (int axis = 0; axis < n; ++axis)
                {
                    Vector v(n); v[axis] = 1.;
                    for (int pass = 0; pass < 2; ++pass)
                    {
                        for (const auto &q : face)
                        {
                            const double c = Dot(v, q);
                            for (int j = 0; j < n; ++j) v[j] -= c * q[j];
                        }
                        for (const auto &q : null)
                        {
                            const double c = Dot(v, q);
                            for (int j = 0; j < n; ++j) v[j] -= c * q[j];
                        }
                    }
                    const double length = Length(v);
                    if (length > 1e-8)
                    {
                        for (double &x : v) x /= length;
                        null.push_back(std::move(v));
                    }
                }
                if (null.size() + active.size() != gradient.size())
                    return finish("quadratic_null_rank_loss");
                if (!null.empty())
                {
                    Matrix reduced(null.size(), Vector(null.size()));
                    Vector rhs(null.size());
                    for (int i = 0; i < null.size(); ++i)
                    {
                        rhs[i] = Dot(null[i], g);
                        Vector hv(n);
                        for (int k = 0; k < n; ++k) hv[k] = Dot(hessian[k], null[i]);
                        for (int j = 0; j < null.size(); ++j)
                            reduced[i][j] = Dot(null[j], hv);
                    }
                    const auto z = Factor(reduced).Apply(rhs);
                    for (int i = 0; i < null.size(); ++i)
                        for (int j = 0; j < n; ++j) direction[j] += z[i] * null[i][j];
                }
                Vector stationarity(n), rhs(active.size());
                for (int i = 0; i < n; ++i)
                    stationarity[i] = g[i] - Dot(hessian[i], direction);
                for (int i = 0; i < active.size(); ++i) rhs[i] = Dot(face[i], stationarity);
                multipliers = Factor(upper).Apply(rhs);
            }
            if (Length(direction) <= 1e-9 * std::max(stepScale, Length(solution)))
            {
                if (active.empty() ||
                    *std::min_element(multipliers.begin(), multipliers.end()) >=
                        -1e-9 * Length(multipliers))
                    return finish("");
                active.erase(active.begin() + (std::min_element(multipliers.begin(),
                    multipliers.end()) - multipliers.begin()));
                continue;
            }
            double fraction = 1.;
            int entering = -1;
            for (int i = 0; i < cuts.size(); ++i)
            {
                if (std::find(active.begin(), active.end(), i) != active.end())
                    continue;
                const double advance = Dot(cuts[i].a, direction);
                if (advance <= 64. * std::numeric_limits<double>::epsilon() *
                                   Length(cuts[i].a) * Length(direction))
                    continue;
                const double bound = (cuts[i].b - Dot(cuts[i].a, solution)) / advance;
                if (bound < fraction)
                {
                    auto remainder = cuts[i].a;
                    for (int pass = 0; pass < 2; ++pass)
                        for (const auto &q : face)
                        {
                            const double c = Dot(remainder, q);
                            for (int j = 0; j < n; ++j) remainder[j] -= c * q[j];
                        }
                    if (Length(remainder) < 1e-10 * Length(cuts[i].a))
                    {
                        // A dependent row cannot block exact motion on this
                        // face. Ignore roundoff-sized advances only while the
                        // row is feasible; finish() still scans every row.
                        if (Dot(cuts[i].a, solution) - cuts[i].b > 1e-8)
                            return finish("quadratic_degenerate_primal_defect");
                        ++result.redundantRows;
                        continue;
                    }
                    fraction = std::max(0., bound);
                    entering = i;
                }
            }
            for (int i = 0; i < n; ++i)
                solution[i] += fraction * direction[i];
            if (entering >= 0) active.push_back(entering);
        }
        return finish("quadratic_iteration_limit");
    }
    catch (const std::runtime_error &error)
    {
        return finish(std::string("quadratic_factor_failure: ") + error.what());
    }
}

// Bound the final combined field correction, including the effects of cuts.
// Model predictions must be recomputed after this operation.
inline Vector BoundCorrection(Vector &coefficients, const Matrix &columns,
                              double radius)
{
    if (columns.empty() || columns.size() > coefficients.size() || !(radius > 0.))
        throw std::logic_error("invalid reduced correction radius or shape");
    Vector step(columns[0].size());
    for (int i = 0; i < columns.size(); ++i)
    {
        if (columns[i].size() != step.size())
            throw std::logic_error("reduced correction column shape mismatch");
        for (int j = 0; j < step.size(); ++j)
            step[j] += coefficients[i] * columns[i][j];
    }
    const double length = Length(step);
    if (length > radius)
    {
        const double fraction = radius / length;
        for (int i = 0; i < columns.size(); ++i) coefficients[i] *= fraction;
        for (double &x : step) x *= fraction;
    }
    return step;
}

struct RadiusResult
{
    QuadraticResult quadratic;
    Matrix hessian;
    double multiplier = 0., rawNorm = 0., norm = 0., complementarity = 0.;
    int solves = 0, steps = 0;
    std::string refinementFailure;
};

// Solve the linear-cut model with ||V z|| <= radius, rather than scaling an
// arbitrarily long minimizer afterwards. A nonnegative Lagrange multiplier
// adds lambda V^T V only to field coefficients (never the restoration slack).
// The conservative inside bracket gives an inexact radius solve; actual field
// bounds and nonlinear acceptance remain the caller's responsibility.
template <typename Factor, typename Cut>
RadiusResult SolveRadiusQuadratic(const Matrix &hessian, const Vector &gradient,
                                 const std::vector<Cut> &cuts, const Vector &initial,
                                 const Matrix &columns, double radius,
                                 int maximumSteps = 500)
{
    if (columns.empty() || columns.size() > gradient.size() ||
        !(radius > 0.) || !std::isfinite(radius))
        throw std::logic_error("invalid radius-constrained model shape or radius");
    RadiusResult result;
    struct Trial
    {
        QuadraticResult quadratic;
        Matrix hessian;
        double multiplier = 0., norm = 0.;
    };
    Matrix gram;
    auto solve = [&](double multiplier) {
        Trial t;
        t.multiplier = multiplier;
        t.hessian = hessian;
        if (multiplier > 0.)
            for (int i = 0; i < columns.size(); ++i)
                for (int j = 0; j < columns.size(); ++j)
                    t.hessian[i][j] += multiplier * gram[i][j];
        t.quadratic = SolveQuadratic<Factor>(t.hessian, gradient, cuts, initial, maximumSteps);
        ++result.solves;
        result.steps += t.quadratic.steps;
        if (t.quadratic.solved)
        {
            auto coefficients = t.quadratic.solution;
            const auto field = BoundCorrection(coefficients, columns,
                                               std::numeric_limits<double>::max());
            if (std::any_of(field.begin(), field.end(), [](double x) { return !std::isfinite(x); }))
            {
                t.quadratic.solved = false;
                t.quadratic.failure = "nonfinite_radius_model_correction";
            }
            else
            {
                t.norm = Length(field);
                if (!std::isfinite(t.norm))
                {
                    t.norm = 0.;
                    t.quadratic.solved = false;
                    t.quadratic.failure = "nonfinite_radius_model_norm";
                }
            }
        }
        return t;
    };
    auto finish = [&](Trial t) {
        result.multiplier = t.multiplier;
        result.norm = t.norm;
        result.complementarity = t.multiplier * std::abs(t.norm*t.norm-radius*radius);
        result.hessian = std::move(t.hessian);
        result.quadratic = std::move(t.quadratic);
        return result;
    };
    auto current = solve(0.);
    result.rawNorm = current.norm;
    if (!current.quadratic.solved || current.norm <= radius)
        return finish(std::move(current));
    gram.assign(columns.size(), Vector(columns.size()));
    double scale = 0.;
    for (int i = 0; i < columns.size(); ++i)
    {
        for (int j = 0; j < columns.size(); ++j)
            gram[i][j] = Dot(columns[i], columns[j]);
        scale = std::max(scale, std::abs(hessian[i][i]) / std::max(1e-30, gram[i][i]));
    }
    double lower = 0., upper = std::max(1e-3, 1e-3 * scale);
    bool bracketed = false;
    for (int bracket = 0; bracket < 14; ++bracket)
    {
        current = solve(upper);
        if (!current.quadratic.solved) return finish(std::move(current));
        if (current.norm <= radius) { bracketed = true; break; }
        lower = upper;
        upper *= 10.;
    }
    if (!bracketed)
    {
        current.quadratic.solved = false;
        current.quadratic.failure = "radius_multiplier_bracket_exhausted";
        return finish(std::move(current));
    }
    for (int refine = 0; refine < 16 && current.norm < .999 * radius; ++refine)
    {
        auto trial = solve(.5 * (lower + upper));
        if (!trial.quadratic.solved)
        {
            result.refinementFailure = trial.quadratic.failure;
            break; // Keep the converged inside-radius bracket, not this failed query.
        }
        if (trial.norm > radius) lower = trial.multiplier;
        else
        {
            upper = trial.multiplier;
            current = std::move(trial);
        }
    }
    return finish(std::move(current));
}
} // namespace Nektar::RANSWallReduced
#endif
