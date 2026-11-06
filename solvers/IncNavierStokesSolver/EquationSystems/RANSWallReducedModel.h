///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_REDUCED_MODEL_H
#define NEKTAR_RANS_WALL_REDUCED_MODEL_H
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace Nektar::RANSWallReduced
{
using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;
// Native 96-response models already spend as much time on repeated null-space
// reconstruction as on their separately converged field responses. Retain
// the established path through64 coordinates; use checked primal faces above.
inline constexpr int DirectFaceThreshold = 64;
struct WitnessKey
{
    int observation, family, region;
};
// Equivalent active faces can list the same physical observations in different
// orders. Choose representatives in original inventory order, so that this
// optimizer detail cannot decide which field responses enter the next model.
inline std::vector<int> SelectWitnesses(std::vector<WitnessKey> keys)
{
    std::sort(keys.begin(),keys.end(),[](const auto &a,const auto &b) {
        return a.observation<b.observation;
    });
    std::set<std::pair<int,int>> regions;
    std::vector<int> selected;
    for (int i=0;i<keys.size();++i)
    {
        const auto &key=keys[i];
        if (key.observation<0 || key.family<0 ||
            (i && keys[i-1].observation==key.observation &&
             (keys[i-1].family!=key.family || keys[i-1].region!=key.region)))
            throw std::logic_error("invalid canonical witness metadata");
        if (regions.emplace(key.family,key.region).second)
            selected.push_back(key.observation);
    }
    return selected;
}
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
    int augmentedFaceSolves = 0, augmentedFaceFallbacks = 0;
    double primal = 0., stationarity = 0., complementarity = 0.;
    std::string failure;
};

struct RaySelection
{
    double fraction = 0., merit = 0.;
    int evaluations = 0, invalid = 0;
    std::vector<std::pair<double,double>> samples;
};

// A bounded surrogate selector only. The caller must equilibrate and assess
// every actual trial against the original residual and qualification inventory.
template <typename Evaluate>
RaySelection SelectSourceRay(double baseMerit, Evaluate evaluate)
{
    if (!(baseMerit > 0.) || !std::isfinite(baseMerit))
        throw std::logic_error("invalid source-ray base merit");
    RaySelection result;
    result.merit = baseMerit;
    auto sample = [&](double f) {
        double value = 0.;
        ++result.evaluations;
        if (!evaluate(f,value)) { ++result.invalid; return; }
        if (!(value >= 0.) || !std::isfinite(value))
            throw std::logic_error("invalid source-ray surrogate value");
        result.samples.emplace_back(f,value);
        if (value < result.merit) { result.merit=value; result.fraction=f; }
    };
    for (int j=0;j<=10;++j) sample(std::ldexp(1.,-j));
    auto best=std::find_if(result.samples.begin(),result.samples.end(),
        [&](const auto &p) { return p.first==result.fraction; });
    if (best!=result.samples.begin() && best!=result.samples.end() && best+1!=result.samples.end())
    {
        const auto a=*(best+1), b=*best, c=*(best-1);
        const double slope1=(b.second*b.second-a.second*a.second)/(b.first-a.first);
        const double slope2=(c.second*c.second-b.second*b.second)/(c.first-b.first);
        const double quadratic=(slope2-slope1)/(c.first-a.first);
        const double linear=slope1-quadratic*(a.first+b.first);
        if (quadratic>0.)
        {
            const double f=-linear/(2.*quadratic);
            if (f>a.first && f<c.first && std::isfinite(f)) sample(f);
        }
    }
    if (result.merit >= baseMerit*(1.-1e-6)) result.fraction=0.;
    return result;
}

struct ModelDerivative
{
    bool available = false;
    Vector coarse, fine;
    double step = 0.;
};

// Surrogate differentiation only. Stay in the callback's evaluation domain,
// and return both scales so the caller can reject an unresolved derivative.
template <typename Evaluate>
ModelDerivative DifferentiateModel(const Vector &at, const Vector &base,
                                   int column, double step, Evaluate evaluate)
{
    if (column < 0 || column >= at.size() || base.empty() ||
        !(step > 0.) || !std::isfinite(step))
        throw std::logic_error("invalid nonlinear model derivative request");
    for (int shrink=0;shrink<4;++shrink)
        for (double sign : {1.,-1.})
        {
            const double h=sign*std::ldexp(step,-shrink);
            auto a=at,b=at;
            a[column]+=h; b[column]+=.5*h;
            Vector coarse,fine;
            if (!evaluate(a,coarse) || !evaluate(b,fine)) continue;
            if (coarse.size()!=base.size() || fine.size()!=base.size())
                throw std::logic_error("nonlinear model derivative output shape");
            for (int j=0;j<base.size();++j)
            {
                coarse[j]=(coarse[j]-base[j])/h;
                fine[j]=(fine[j]-base[j])/(.5*h);
                if (!std::isfinite(coarse[j]) || !std::isfinite(fine[j]))
                    throw std::logic_error("nonfinite nonlinear model derivative");
            }
            return {true,std::move(coarse),std::move(fine),h};
        }
    return {};
}

struct ModelLinearization
{
    struct Stencil
    {
        int column, halvings;
        double initialError, finalError, coarseStep;
    };
    bool available = false;
    Vector at, value;
    Matrix derivatives;
    int availableDerivatives = 0, anchorReductions = 0, attempts = 0;
    int refinedDerivatives = 0, derivativeRefinements = 0;
    double error = 0., initialError = 0., minimumStep = 0.;
    std::vector<Stencil> stencils;
};

// A constrained corner can forbid both signs of a coordinate probe even when
// the model is smooth on its domain. Retry its linearization toward the base
// of the ray. This changes only the internal model's anchor; the incumbent
// proposal and actual nonlinear acceptance tests remain the caller's authority.
template <typename Evaluate, typename RelativeError>
ModelLinearization LinearizeModelOnRay(const Vector &at, const Vector &value,
                                      double step, int maximumReductions,
                                      double maximumError, Evaluate evaluate,
                                      RelativeError relativeError)
{
    if (at.empty() || value.empty() || !(step > 0.) || !std::isfinite(step) ||
        maximumReductions < 0 || maximumReductions > 3 ||
        !(maximumError > 0.) || !std::isfinite(maximumError))
        throw std::logic_error("invalid source-model ray linearization request");
    ModelLinearization result;
    for (int anchor = 0; anchor <= maximumReductions; ++anchor)
    {
        ModelLinearization current;
        current.at = at;
        current.anchorReductions = anchor;
        current.attempts = anchor + 1;
        for (double &entry : current.at) entry *= std::ldexp(1., -anchor);
        if (anchor == 0) current.value = value;
        else if (!evaluate(current.at, current.value))
        {
            result.attempts = current.attempts;
            continue;
        }
        if (current.value.size() != value.size())
            throw std::logic_error("source-model anchor changed residual shape");
        for (int i = 0; i < at.size(); ++i)
        {
            auto d = DifferentiateModel(current.at, current.value, i, step, evaluate);
            if (!d.available) break;
            double error = relativeError(d);
            if (!(error >= 0.) || !std::isfinite(error))
                throw std::logic_error("invalid source-model derivative error");
            current.initialError = std::max(current.initialError, error);
            const double beforeError = error;
            int halvings = 0;
            bool refined = false;
            // Reuse the previous half-step sample as the new coarse action.
            // Each retry needs only one additional, guarded quarter sample.
            // A topology/roundoff error that does not resolve remains a
            // failure; it does not relax the caller's accuracy threshold.
            for (int retry = 0; retry < 3 && error > maximumError; ++retry)
            {
                auto atQuarter = current.at;
                atQuarter[i] += .25*d.step;
                Vector value;
                if (!evaluate(atQuarter, value)) break;
                if (value.size() != current.value.size())
                    throw std::logic_error("source-model quarter sample changed shape");
                for (int j = 0; j < value.size(); ++j)
                {
                    value[j] = (value[j]-current.value[j])/(.25*d.step);
                    if (!std::isfinite(value[j]))
                        throw std::logic_error("nonfinite refined source-model derivative");
                }
                d.coarse = std::move(d.fine);
                d.fine = std::move(value);
                d.step *= .5;
                error = relativeError(d);
                if (!(error >= 0.) || !std::isfinite(error))
                    throw std::logic_error("invalid refined source-model derivative error");
                ++current.derivativeRefinements;
                ++halvings;
                refined = true;
            }
            current.refinedDerivatives += refined;
            if (refined) current.stencils.push_back({i,halvings,beforeError,error,std::abs(d.step)});
            current.minimumStep = current.minimumStep > 0.
                ? std::min(current.minimumStep,std::abs(d.step)) : std::abs(d.step);
            current.error = std::max(current.error, error);
            current.derivatives.push_back(std::move(d.fine));
        }
        current.availableDerivatives = current.derivatives.size();
        current.available = current.availableDerivatives == at.size() &&
                            current.error <= maximumError;
        result = std::move(current);
        if (result.available) return result;
    }
    return result;
}

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
                              int maximumSteps = 120,
                              const std::vector<int> &initialActive = {})
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
    // Cut geometry is constant throughout this model solve. The complete
    // qualification inventory can contain more than100,000 rows.
    Vector cutLengths;
    cutLengths.reserve(cuts.size());
    for (const auto &cut : cuts) cutLengths.push_back(Length(cut.a));
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
    for (int id:initialActive)
    {
        if (id<0 || id>=cuts.size() ||
            std::find(active.begin(),active.end(),id)!=active.end())
            throw std::logic_error("invalid reduced quadratic warm face index");
        if (std::abs(Dot(cuts[id].a,solution)-cuts[id].b)>1e-8)
            return finish("quadratic_inconsistent_warm_face");
        active.push_back(id);
    }
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
                if (length < 1e-10 * cutLengths[active[i]])
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
                bool augmented = false;
                if (n > DirectFaceThreshold)
                {
                    // The orthonormal active face also defines a small saddle
                    // problem in response coordinates. This avoids rebuilding
                    // the complete null basis at every large-model QP step.
                    // The governing fields are not part of this matrix.
                    try
                    {
                        const int size = n + face.size();
                        double scale = 1.;
                        for (int i = 0; i < n; ++i)
                            scale = std::max(scale, std::sqrt(std::abs(hessian[i][i])));
                        Matrix saddle(size, Vector(size));
                        Vector rhs(size);
                        for (int i = 0; i < n; ++i)
                        {
                            std::copy(hessian[i].begin(), hessian[i].end(), saddle[i].begin());
                            rhs[i] = g[i];
                            for (int j = 0; j < face.size(); ++j)
                                saddle[i][n+j] = saddle[n+j][i] = scale * face[j][i];
                        }
                        const Factor factor(saddle);
                        auto answer = factor.Apply(rhs);
                        for (int refine = 0; refine < 2; ++refine)
                        {
                            auto residual = rhs;
                            for (int i = 0; i < size; ++i)
                                residual[i] -= Dot(saddle[i], answer);
                            const auto correction = factor.Apply(residual);
                            for (int i = 0; i < size; ++i) answer[i] += correction[i];
                        }
                        std::copy(answer.begin(), answer.begin()+n, direction.begin());
                        for (int pass = 0; pass < 2; ++pass)
                            for (const auto &q : face)
                            {
                                const double c = Dot(direction, q);
                                for (int i = 0; i < n; ++i) direction[i] -= c*q[i];
                            }
                        Vector defect(n), hd(n);
                        for (int i = 0; i < n; ++i)
                        {
                            hd[i] = Dot(hessian[i], direction);
                            defect[i] = g[i] - hd[i];
                        }
                        for (int pass = 0; pass < 2; ++pass)
                            for (const auto &q : face)
                            {
                                const double c = Dot(defect, q);
                                for (int i = 0; i < n; ++i) defect[i] -= c*q[i];
                            }
                        if (!std::isfinite(Length(defect)) ||
                            Length(defect) > 1e-10 * std::max({1.,Length(g),Length(hd)}))
                            throw std::runtime_error("inaccurate augmented face solve");
                        ++result.augmentedFaceSolves;
                        augmented = true;
                    }
                    catch (const std::runtime_error &)
                    {
                        ++result.augmentedFaceFallbacks;
                        std::fill(direction.begin(), direction.end(), 0.);
                    }
                }
                if (!augmented)
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
                }
                Vector stationarity(n), rhs(active.size());
                for (int i = 0; i < n; ++i)
                    stationarity[i] = g[i] - Dot(hessian[i], direction);
                for (int i = 0; i < active.size(); ++i) rhs[i] = Dot(face[i], stationarity);
                multipliers = Factor(upper).Apply(rhs);
            }
            const double directionLength = Length(direction);
            if (directionLength <= 1e-9 * std::max(stepScale, Length(solution)))
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
                                   cutLengths[i] * directionLength)
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
                    if (Length(remainder) < 1e-10 * cutLengths[i])
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
    int augmentedFaceSolves = 0, augmentedFaceFallbacks = 0;
    int ballRedundantCuts = 0;
    int warmFaceSolves = 0, warmFaceFallbacks = 0;
    bool unconstrainedAvailable = false;
    int failedPositiveQueries = 0;
    std::string unconstrainedFailure, refinementFailure;
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
    QuadraticResult previousQuadratic;
    struct Trial
    {
        QuadraticResult quadratic;
        Matrix hessian;
        double multiplier = 0., norm = 0.;
    };
    Matrix gram;
    auto buildGram = [&]() {
        if (!gram.empty()) return;
        gram.assign(columns.size(), Vector(columns.size()));
        for (int i = 0; i < columns.size(); ++i)
            for (int j = 0; j < columns.size(); ++j)
                gram[i][j] = Dot(columns[i], columns[j]);
    };
    // A conservative field-Gram bound certifies cuts that cannot bind
    // anywhere in this displacement ball. Keep slack-variable models intact.
    // Recheck every original cut at the final inside-ball solution below.
    std::vector<Cut> filtered;
    std::vector<int> originalIds;
    const std::vector<Cut> *modelCuts = &cuts;
    if (columns.size() > DirectFaceThreshold &&
        columns.size() == gradient.size() && cuts.size() > 4096)
    {
        buildGram();
        double lower = std::numeric_limits<double>::max(), upper = 0.;
        for (int i = 0; i < gram.size(); ++i)
        {
            long double off = 0.;
            for (int j = 0; j < gram.size(); ++j)
                if (i != j) off += std::abs(gram[i][j]);
            lower = std::min(lower, gram[i][i] - static_cast<double>(off));
            upper = std::max(upper, gram[i][i] + static_cast<double>(off));
        }
        if (std::isfinite(upper) && lower > .5*upper)
        {
            const double safeLower = lower - 64.*(columns.size()+columns[0].size()+1)*
                std::numeric_limits<double>::epsilon()*upper;
            if (safeLower > 0.)
            {
                const double coefficientRadius = radius/std::sqrt(safeLower);
                for (int i = 0; i < cuts.size(); ++i)
                {
                    if (cuts[i].a.size() != gradient.size())
                        throw std::logic_error("invalid radius model cut shape");
                    const double extent = coefficientRadius*Length(cuts[i].a);
                    if (cuts[i].b > extent + 1e-10*std::max({1.,std::abs(cuts[i].b),extent}))
                        ++result.ballRedundantCuts;
                    else
                    {
                        filtered.push_back(cuts[i]);
                        originalIds.push_back(i);
                    }
                }
                if (result.ballRedundantCuts) modelCuts = &filtered;
            }
        }
    }
    auto solve = [&](double multiplier) {
        Trial t;
        t.multiplier = multiplier;
        t.hessian = hessian;
        if (multiplier > 0.)
            for (int i = 0; i < columns.size(); ++i)
                for (int j = 0; j < columns.size(); ++j)
                    t.hessian[i][j] += multiplier * gram[i][j];
        auto record=[&](const QuadraticResult &q) {
            ++result.solves;
            result.steps+=q.steps;
            result.augmentedFaceSolves+=q.augmentedFaceSolves;
            result.augmentedFaceFallbacks+=q.augmentedFaceFallbacks;
        };
        // Radius queries change the objective multiplier, not the cut geometry.
        // Reuse a previously feasible point and its independent binding face
        // in large models. Wrong-sign multipliers remove obsolete active cuts
        // through the same active-set solve. Any numerical failure retries the
        // original feasible initial point with no assumed active face.
        if (columns.size()>DirectFaceThreshold && previousQuadratic.solved)
        {
            ++result.warmFaceSolves;
            t.quadratic=SolveQuadratic<Factor>(t.hessian,gradient,*modelCuts,
                previousQuadratic.solution,maximumSteps,previousQuadratic.active);
            record(t.quadratic);
            if (!t.quadratic.solved) ++result.warmFaceFallbacks;
        }
        if (!t.quadratic.solved)
        {
            t.quadratic=SolveQuadratic<Factor>(t.hessian,gradient,*modelCuts,initial,maximumSteps);
            record(t.quadratic);
        }
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
        if (t.quadratic.solved) previousQuadratic=t.quadratic;
        return t;
    };
    auto finish = [&](Trial t) {
        if (modelCuts != &cuts)
        {
            for (int &id : t.quadratic.active) id = originalIds.at(id);
            if (t.quadratic.solved && t.norm <= radius*(1.+1e-12))
            {
                for (const auto &cut : cuts)
                    t.quadratic.primal = std::max(t.quadratic.primal,
                        Dot(cut.a,t.quadratic.solution)-cut.b);
                if (t.quadratic.primal > 1e-8)
                {
                    t.quadratic.solved = false;
                    t.quadratic.failure = "radius_original_cut_defect";
                }
            }
        }
        result.multiplier = t.multiplier;
        result.norm = t.norm;
        result.complementarity = t.multiplier * std::abs(t.norm*t.norm-radius*radius);
        result.hessian = std::move(t.hessian);
        result.quadratic = std::move(t.quadratic);
        return result;
    };
    auto current = solve(0.);
    result.unconstrainedAvailable = current.quadratic.solved;
    result.unconstrainedFailure = current.quadratic.failure;
    result.rawNorm = current.norm;
    if ((current.quadratic.solved && current.norm <= radius) ||
        current.quadratic.failure == "quadratic_infeasible_start")
        return finish(std::move(current));
    // A failed radius-free query need not prevent a well-conditioned solve
    // inside the displacement ball. The native model can spend its budget
    // resolving a minimizer far outside that ball. Try positive multipliers
    // from the same feasible point; an unsolved query never supplies a norm
    // or a lower bracket, and no partial iterate becomes a solved model.
    buildGram();
    double scale = 0.;
    for (int i = 0; i < columns.size(); ++i)
    {
        scale = std::max(scale, std::abs(hessian[i][i]) / std::max(1e-30, gram[i][i]));
    }
    double lower = 0., upper = std::max(1e-3, 1e-3 * scale);
    bool bracketed = false;
    for (int bracket = 0; bracket < 14; ++bracket)
    {
        current = solve(upper);
        if (!current.quadratic.solved)
        {
            ++result.failedPositiveQueries;
            upper *= 10.;
            continue;
        }
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
