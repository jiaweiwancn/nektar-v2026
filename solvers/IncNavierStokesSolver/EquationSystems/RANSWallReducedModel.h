///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_WALL_REDUCED_MODEL_H
#define NEKTAR_RANS_WALL_REDUCED_MODEL_H
#include <algorithm>
#include <cmath>
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
struct QuadraticResult
{
    bool solved = false;
    Vector solution, multipliers;
    std::vector<int> active;
    int steps = 0;
    double primal = 0., stationarity = 0., complementarity = 0.;
    std::string failure;
};

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
        std::vector<Vector> inverse;
        for (const auto &cut : cuts)
            inverse.push_back(objective.Apply(cut.a));
        while (result.steps < maximumSteps)
        {
            ++result.steps;
            Vector g(n);
            for (int i = 0; i < n; ++i)
                g[i] = -gradient[i] - Dot(hessian[i], solution);
            auto direction = objective.Apply(g);
            if (!active.empty())
            {
                Matrix normal(active.size(), Vector(active.size()));
                Vector rhs(active.size());
                for (int i = 0; i < active.size(); ++i)
                {
                    rhs[i] = Dot(cuts[active[i]].a, direction);
                    for (int j = 0; j < active.size(); ++j)
                        normal[i][j] = Dot(cuts[active[i]].a, inverse[active[j]]);
                }
                multipliers = Factor(normal).Apply(rhs);
                for (int i = 0; i < active.size(); ++i)
                    for (int j = 0; j < n; ++j)
                        direction[j] -= multipliers[i] * inverse[active[i]][j];
            }
            if (Length(direction) < 1e-9 * std::max(1., Length(solution)))
            {
                if (active.empty() ||
                    *std::min_element(multipliers.begin(), multipliers.end()) >= -1e-8)
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
                if (advance <= 1e-12)
                    continue;
                const double bound = (cuts[i].b - Dot(cuts[i].a, solution)) / advance;
                if (bound < fraction)
                {
                    fraction = std::max(0., bound);
                    entering = i;
                }
            }
            for (int i = 0; i < n; ++i)
                solution[i] += fraction * direction[i];
            if (entering >= 0)
            {
                auto row = cuts[entering].a;
                std::vector<Vector> q;
                for (int id : active)
                {
                    auto a = cuts[id].a;
                    for (int pass = 0; pass < 2; ++pass)
                        for (const auto &v : q)
                        {
                            const double c = Dot(a, v);
                            for (int j = 0; j < n; ++j)
                                a[j] -= c * v[j];
                        }
                    if (Length(a) > 1e-8)
                    {
                        const double length = Length(a);
                        for (double &x : a) x /= length;
                        q.push_back(std::move(a));
                    }
                }
                for (int pass = 0; pass < 2; ++pass)
                    for (const auto &v : q)
                    {
                        const double c = Dot(row, v);
                        for (int j = 0; j < n; ++j) row[j] -= c * v[j];
                    }
                if (Length(row) < 1e-8)
                    return finish("quadratic_dependent_active_rows");
                active.push_back(entering);
            }
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
} // namespace Nektar::RANSWallReduced
#endif
