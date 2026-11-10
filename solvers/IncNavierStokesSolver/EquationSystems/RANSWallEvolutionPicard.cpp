///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSWallTiming.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallDurableProgress.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallSolveDiagnostics.h>
#include <algorithm>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <boost/property_tree/json_parser.hpp>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <numeric>

namespace Nektar
{
using namespace SolverUtils;
using RANSWallAlgebra::Maximum;
using RANSWallAlgebra::Multiply;
using RANSWallAlgebra::Require;
namespace J = RANSWallSolveJson;

namespace
{
template <class Factor>
double DualNorm(const Factor &mass, const RANSWallAlgebra::Vector &load)
{
    const auto dual = mass.Apply(load);
    return std::sqrt(std::max(
        0., std::inner_product(load.begin(), load.end(), dual.begin(), 0.)));
}
} // namespace

RANSWallEvolutionPatch::WallMatrix RANSWallEvolutionPatch::PicardKMass(
    const std::vector<int> &rows) const
{
    const int first = offset[K()];
    WallMatrix mass(rows.size(), Vector(rows.size()));
    for (int i = 0; i < rows.size(); ++i)
        for (int j = 0; j < rows.size(); ++j)
            mass[i][j] = m_mass[K()][rows[i] - first][rows[j] - first];
    return mass;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::PicardKResidual(
    const Vector &state, const Vector &history, double alpha, double time,
    const RANSWallEvolutionData &data, const std::vector<int> &rows) const
{
    const auto residual = Residual(state, history, alpha, time, data);
    Vector result;
    result.reserve(rows.size());
    for (int row : rows)
        result.push_back(residual[row]);
    return result;
}

RANSWallEvolutionPatch::PicardBlock RANSWallEvolutionPatch::AssemblePicardBlock(
    int block, const Vector &s, const Vector &history, double alpha,
    double time, const RANSWallEvolutionData &data,
    const std::map<int, double> &external) const
{
    RANSWallTiming::Scope timing(block==0?"flow.assembly":block==Epsilon()?"epsilon.assembly":"k.assembly");
    const bool flow = block == 0;
    PicardBlock system;
    std::vector<int> index(Size(), -1);
    for (int row = 0; row < Size(); ++row)
    {
        const bool selected =
            flow ? row < offset[K()] || row >= DifferentialSize()
                 : row >= offset[block] && row < offset[block + 1];
        if (selected && !external.count(row))
        {
            index[row] = system.rows.size();
            system.rows.push_back(row);
        }
    }
    const auto n = system.rows.size();
    system.matrix.assign(n, Vector(n));
    system.residual.assign(n, 0.);

    const double nu   = m_parameters.geometry.molecularNu;
    Refresh(s, time, &data);

    auto addScalar = [&](const RANSWallGlobalPoint &point,
                         const RANSVector &velocity, double weight,
                         double diffusion, double rate, double source) {
        const double value = point.Value(s, offset[block]);
        const auto grad    = point.Gradient(s, offset[block]);
        for (int i = 0; i < point.ids.size(); ++i)
        {
            const int row = index[offset[block] + point.ids[i]];
            if (row < 0)
                continue;
            const double test = point.phi[i];
            const double adv  = point.derivative[0][i] * velocity[0] +
                               point.derivative[2][i] * velocity[2];
            const double dx = point.derivative[0][i],
                         dz = point.derivative[2][i];
            system.residual[row] +=
                weight * (test * (rate * value - source) - adv * value +
                          diffusion * (dx * grad[0] + dz * grad[2]));
            for (int j = 0; j < point.ids.size(); ++j)
            {
                const int col = index[offset[block] + point.ids[j]];
                if (col >= 0)
                    system.matrix[row][col] +=
                        weight * ((test * rate - adv) * point.phi[j] +
                                  diffusion * (dx * point.derivative[0][j] +
                                               dz * point.derivative[2][j]));
            }
        }
    };
    auto evaluate = [&](const VolumePoint &q, const Ray *ray,
                        const RANSVector &force,
                        const RANSWallGlobalPoint *forcingPressure) {
        const bool inner = ray != nullptr;
        RANSVector u{};
        RANSTensor grad{};
        for (int a = 0; a < m_dimension; ++a)
        {
            const int axis = PhysicalAxis(a);
            u[axis]        = q.v.Value(s, offset[a]);
            grad[axis]     = q.v.Gradient(s, offset[a]);
        }
        const double k = q.v.Value(s, offset[K()]);
        double viscosity;
        RANSTensor stress{};
        const int pBlock = inner ? InnerPressure() : OuterPressure();
        if (inner)
        {
            const double nt = ray->profile->Original().EddyViscosity(
                Distance(q.v.coordinate));
            const auto law = InnerStress(k, nt, grad);
            if (!flow)
            {
                double rate = 0.;
                const auto source=WallSource(s, q.v.coordinate, law.production, *ray->profile,
                           &ray->destruction, &rate);
                addScalar(q.v, u, q.weight, nu + nt / SigmaK(), rate,
                          (CappedDestruction()?law.production:source.kSource) +
                              data.innerKSource(q.v.coordinate, time));
                return;
            }
            viscosity = nu + law.tangentialTangent;
            stress =
                RANSWallInnerMeanStress(law, nu, q.p.Value(s, offset[pBlock]));
        }
        else
        {
            const double epsilon = q.e.Value(s, offset[Epsilon()]);
            const auto law = OuterClosure(k, epsilon, grad, q.v.coordinate);
            if (!flow)
            {
                if (block == K())
                {
                    const double rate =
                        (epsilon - std::min(0., law.production)) / k;
                    addScalar(q.v, u, q.weight, law.kDiffusion, rate,
                              std::max(0., law.production) +
                                  data.outerKSource(q.v.coordinate, time));
                }
                else
                {
                    addScalar(
                        q.e, u, q.weight, law.epsilonDiffusion,
                        Epsilon2() * epsilon / k,
                        law.epsilonProduction + law.yap +
                            data.outerEpsilonSource(q.v.coordinate, time));
                }
                return;
            }
            viscosity = law.eddyViscosity + m_regularizationGamma +
                        m_regularizationScale * m_parameters.model.Kappa() *
                            m_robustUstar * Distance(q.v.coordinate);
            const double pi = q.p.Value(s, offset[pBlock]);
            for (int a = 0; a < 3; ++a)
                for (int b = 0; b < 3; ++b)
                    stress[a][b] = -law.reynoldsStress[a][b] +
                                   (a == b ? 2 * k / 3 - pi : 0.);
        }
        const double divergence = grad[0][0] + grad[2][2];
        for (int p = 0; p < q.p.ids.size(); ++p)
            system.residual[index[offset[pBlock] + q.p.ids[p]]] +=
                q.weight * q.p.phi[p] * divergence;

        // Oseen advection and a positive symmetric-strain part are implicit.
        // The load at s includes the complete stress, so its remainder
        // (including CLS nonlinear stresses and the reference-strain load)
        // is deferred. Both regions use the shared velocity test space.
        for (int a = 0; a < m_dimension; ++a)
        {
            const int aa = PhysicalAxis(a);
            for (int i = 0; i < q.v.ids.size(); ++i)
            {
                const int row               = index[offset[a] + q.v.ids[i]];
                const double dx             = q.v.derivative[0][i],
                             dz             = q.v.derivative[2][i];
                const double testDerivative = q.v.derivative[aa][i];
                if (row >= 0)
                {
                    system.residual[row] +=
                        q.weight * (dx * (stress[aa][0] - u[aa] * u[0]) +
                                    dz * (stress[aa][2] - u[aa] * u[2]) -
                                    q.v.phi[i] * force[aa]);
                    for (int b = 0; b < m_dimension; ++b)
                    {
                        const int bb = PhysicalAxis(b);
                        for (int j = 0; j < q.v.ids.size(); ++j)
                        {
                            const int col = index[offset[b] + q.v.ids[j]];
                            if (col < 0)
                                continue;
                            double value = viscosity * q.v.derivative[bb][i] *
                                           q.v.derivative[aa][j];
                            if (a == b)
                                value +=
                                    viscosity * (dx * q.v.derivative[0][j] +
                                                 dz * q.v.derivative[2][j]) -
                                    (dx * u[0] + dz * u[2]) * q.v.phi[j];
                            system.matrix[row][col] += q.weight * value;
                        }
                    }
                    if (forcingPressure && aa == 0)
                    {
                        // Flat tangential pressure gradient in the roughness
                        // load. The mechanical -2 grad(k)/3 and convection
                        // parts of CartesianForcing stay frozen in this solve.
                        const auto &fp = *forcingPressure;
                        for (int p = 0; p < fp.ids.size(); ++p)
                            system.matrix[row][index[offset[OuterPressure()] +
                                                     fp.ids[p]]] -=
                                q.weight * q.v.phi[i] * fp.derivative[0][p];
                    }
                }
                if (row < 0)
                    continue;
                for (int p = 0; p < q.p.ids.size(); ++p)
                {
                    const int pressure = index[offset[pBlock] + q.p.ids[p]];
                    const double entry = q.weight * q.p.phi[p] * testDerivative;
                    system.matrix[row][pressure] -= entry;
                    system.matrix[pressure][row] += entry;
                }
            }
        }
    };

    if (block != Epsilon())
        for (const auto &ray : m_rays)
        {
            RANSVector force{};
            if (flow)
                force = CartesianForcing(s, ray.coordinate, time, data, false);
            for (const auto &q : ray.points)
            {
                const bool rough = flow && Distance(q.v.coordinate) <=
                                               m_parameters.geometry.roughness;
                evaluate(q, &ray, rough ? force : RANSVector{},
                         rough ? &ray.outerP : nullptr);
            }
        }
    for (const auto &q : m_outer)
        evaluate(q, nullptr, {}, nullptr);

    const auto body = flow ? OwnedBodyLoad(time, data) : Vector{};
    for (int b = flow ? 0 : block; b < (flow ? K() : block + 1); ++b)
    {
        Vector change(size[b]);
        for (int j = 0; j < size[b]; ++j)
            change[j] = s[offset[b] + j] - history[offset[b] + j];
        const auto storage = Multiply(m_mass[b], change);
        for (int i = 0; i < size[b]; ++i)
        {
            const int row = index[offset[b] + i];
            if (row < 0)
                continue;
            system.residual[row] += alpha * storage[i];
            if (flow)
                system.residual[row] -= body[offset[b] + i];
            if (m_eta)
                system.residual[row] -= m_eta * m_defect[offset[b] + i];
            for (int j = 0; j < size[b]; ++j)
            {
                const int col = index[offset[b] + j];
                if (col >= 0)
                    system.matrix[row][col] +=
                        (alpha + m_picardPseudoAlpha) * m_mass[b][i][j];
            }
        }
    }
    if (flow)
    {
        const int gauge = index[offset[Gauge()]];
        for (int side = 0; side < 2; ++side)
            for (int i = 0; i < m_pressureMean[side].size(); ++i)
            {
                const int id              = offset[InnerPressure() + side] + i,
                          row             = index[id];
                const double w            = m_pressureMean[side][i];
                system.matrix[row][gauge] = system.matrix[gauge][row] = w;
                system.residual[row] += s[offset[Gauge()]] * w;
                system.residual[gauge] += s[id] * w;
            }
    }
    else if (block == Epsilon())
    {
        // Hold Phi fixed for this scalar solve. Free transport rows retain
        // the columns of the matching trace, including its update.
        const auto trace = DissipationTrace(s, time, data);
        for (const auto &[id, value] : trace[0])
        {
            const int row = index[offset[block] + id];
            Require(row >= 0,
                    "Picard matching trace overlaps external lifting");
            std::fill(system.matrix[row].begin(), system.matrix[row].end(), 0.);
            system.matrix[row][row] = 1.;
            system.residual[row] =
                s[offset[block] + id] - (value + trace[1].at(id));
        }
    }
    return system;
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::ProjectPicardGround(
    const Vector &state, const PicardBlock &block, const Factor &factor,
    const Vector &step, double relaxation, int &cuts,
    std::string *diagnostic) const
{
    // Project the proposed k increment in the frozen mass-dual residual
    // norm. The metric inverse is A^-1 M A^-T on the free k coefficients.
    // Ground cuts constrain the search only. They do not replace any k row.
    struct Cut
    {
        int geometry;
        double x, target;
        Vector row, inverse;
    };
    const int n = block.rows.size(), kOffset = offset[K()];
    std::vector<int> index(size[K()], -1);
    for (int i = 0; i < n; ++i)
        index[block.rows[i] - kOffset] = i;
    Vector proposed = step;
    for (double &entry : proposed)
        entry *= relaxation;
    Vector direction = proposed;
    std::vector<Cut> active;
    Vector multiplier;
    auto dot = [](const Vector &a, const Vector &b) {
        return std::inner_product(a.begin(), a.end(), b.begin(), 0.);
    };
    for (cuts = 0; cuts <= 20; ++cuts)
    {
        auto trial = state;
        for (int i = 0; i < n; ++i)
            trial[block.rows[i]] += direction[i];
        const auto extrema = WallTraceR4(trial);
        const auto worst   = std::min_element(
            extrema.begin(), extrema.end(), [](const auto &a, const auto &b) {
                return a.value - a.error < b.value - b.error;
            });
        if (worst->value >= worst->error)
        {
            if (diagnostic)
            {
                Vector change(n), dualRows(n);
                for (int i = 0; i < n; ++i)
                    change[i] = direction[i] - proposed[i];
                const auto dual = Factor(PicardKMass(block.rows))
                                      .Apply(Multiply(block.matrix, change));
                std::vector<std::string> constraints;
                double complementarity = 0., primal = 0.;
                for (int i = 0; i < active.size(); ++i)
                {
                    const auto &c    = active[i];
                    const double gap = dot(c.row, direction) - c.target;
                    primal           = std::max(primal, std::max(0., -gap));
                    complementarity  = std::max(complementarity,
                                                std::abs(gap * multiplier[i]));
                    for (int j = 0; j < n; ++j)
                        dualRows[j] -= multiplier[i] * c.row[j];
                    constraints.push_back(
                        J::Object({{"x", J::Number(c.x)},
                                   {"geometry", J::Number(c.geometry)},
                                   {"row", J::Numbers(c.row)},
                                   {"target", J::Number(c.target)},
                                   {"multiplier", J::Number(multiplier[i])},
                                   {"gap", J::Number(gap)}}));
                }
                for (int i = 0; i < n; ++i)
                    for (int j = 0; j < n; ++j)
                        dualRows[i] += block.matrix[j][i] * dual[j];
                *diagnostic =
                    J::Object({{"active", J::List(constraints)},
                               {"primal_defect", J::Number(primal)},
                               {"complementarity", J::Number(complementarity)},
                               {"stationarity", J::Number(Maximum(dualRows))},
                               {"ground_margin",
                                J::Number(worst->value - worst->error)}});
            }
            return direction;
        }
        const auto &point = worst->row;
        if (cuts == 20 ||
            std::any_of(active.begin(), active.end(), [&](const Cut &c) {
                return c.geometry == point.geometry &&
                       std::abs(c.x - point.coordinate[0]) < 1e-12;
            }))
            return {};

        Vector row(n);
        for (int j = 0; j < point.ids.size(); ++j)
        {
            const int col = index[point.ids[j]];
            if (col >= 0)
                row[col] = point.phi[j];
        }
        const auto adjoint = factor.Apply(row, true);
        Vector load(n);
        for (int i = 0; i < n; ++i)
            for (int j = 0; j < n; ++j)
                load[i] += m_mass[K()][block.rows[i] - kOffset]
                                 [block.rows[j] - kOffset] *
                           adjoint[j];
        auto inverse        = factor.Apply(load);
        const double square = dot(row, inverse);
        if (!(square > 0.))
            return {};
        const double scale = std::sqrt(square);
        for (int j = 0; j < n; ++j)
        {
            row[j] /= scale;
            inverse[j] /= scale;
        }
        const double base = point.Value(state, kOffset);
        const double reserve =
            std::max(8. * static_cast<double>(worst->error), .005 * base);
        active.push_back({point.geometry, point.coordinate[0],
                          (reserve - base) / scale, std::move(row),
                          std::move(inverse)});

        while (!active.empty())
        {
            const int count = active.size();
            WallMatrix gram(count, Vector(count));
            Vector rhs(count);
            for (int i = 0; i < count; ++i)
            {
                rhs[i] = active[i].target - dot(active[i].row, proposed);
                for (int j = 0; j < count; ++j)
                    gram[i][j] = dot(active[i].row, active[j].inverse);
            }
            try
            {
                multiplier = Factor(gram).Apply(rhs);
            }
            catch (const std::runtime_error &)
            {
                return {};
            }
            const auto leaving =
                std::min_element(multiplier.begin(), multiplier.end());
            if (*leaving < 0.)
            {
                active.erase(active.begin() + (leaving - multiplier.begin()));
                continue;
            }
            direction = proposed;
            for (int i = 0; i < count; ++i)
                for (int j = 0; j < n; ++j)
                    direction[j] += multiplier[i] * active[i].inverse[j];
            break;
        }
        if (active.empty())
            return {};
    }
    return {};
}

void RANSWallEvolutionPatch::DiagnosePicardK(
    const Vector &state, const Vector &history, double alpha, double time,
    const RANSWallEvolutionData &data,
    const std::map<int, double> &external) const
{
    Require(!std::ifstream(m_picardKDiagnostic).good(),
            "refusing to overwrite Picard k diagnostic");
    auto frozen =
        AssemblePicardBlock(K(), state, history, alpha, time, data, external);
    frozen.residual =
        PicardKResidual(state, history, alpha, time, data, frozen.rows);
    auto exact      = frozen;
    exact.matrix    = PicardKJacobian(state, alpha, time, data, frozen.rows);
    const auto mass = PicardKMass(frozen.rows);
    const Factor metric(mass);
    const double scale = m_parameters.model.Q() * std::pow(m_robustUstar, 3) /
                         m_parameters.geometry.height *
                         std::sqrt(m_volumes[K()]);
    auto norm  = [&](const Vector &r) { return DualNorm(metric, r) / scale; };
    auto trial = [&](const Vector &step, double fraction) {
        auto result = state;
        for (int i = 0; i < frozen.rows.size(); ++i)
            result[frozen.rows[i]] += fraction * step[i];
        return result;
    };
    auto extrema = [&](const Vector &q) {
        std::vector<std::string> observations;
        for (const auto &point : WallTraceR4(q))
            observations.push_back(
                J::Object({{"x", J::Number(point.row.coordinate[0])},
                           {"geometry", J::Number(point.row.geometry)},
                           {"value", J::Number(point.value)},
                           {"error", J::Number(point.error)}}));
        return J::List(observations);
    };
    const auto base = ScalarDomain(state, time, data);
    auto evaluate   = [&](const Vector &q, ScalarDomainResult &domain) {
        bool ground = false;
        WallTraceJSON(q, ground);
        if (!ground)
            return Vector{};
        domain = ScalarDomain(q, time, data);
        return domain.admissible
                     ? PicardKResidual(q, history, alpha, time, data, frozen.rows)
                     : Vector{};
    };
    auto dot = [](const Vector &a, const Vector &b) {
        return std::inner_product(a.begin(), a.end(), b.begin(), 0.);
    };
    const auto contacts = WallTraceR4(state);
    std::vector<Vector> active;
    Vector inward(frozen.rows.size()), sampled(inward.size()),
        normalization(inward.size());
    for (const auto &point : contacts)
    {
        if (point.value > 1e-7)
            continue;
        Vector row(inward.size());
        for (int i = 0; i < point.row.ids.size(); ++i)
        {
            const auto id = offset[K()] + point.row.ids[i];
            const auto col =
                std::lower_bound(frozen.rows.begin(), frozen.rows.end(), id);
            if (col != frozen.rows.end() && *col == id)
                row[col - frozen.rows.begin()] = point.row.phi[i];
        }
        if (Maximum(row) == 0.)
            continue;
        for (int i = 0; i < row.size(); ++i)
            inward[i] += row[i] / dot(row, row);
        active.push_back(std::move(row));
    }
    auto shape = [&](double z, Vector &direction) {
        const auto point = Point(K(), {1.57, 0., z});
        for (int j = 0; j < point.ids.size(); ++j)
        {
            const auto id = offset[K()] + point.ids[j];
            const auto col =
                std::lower_bound(frozen.rows.begin(), frozen.rows.end(), id);
            if (col != frozen.rows.end() && *col == id)
                direction[col - frozen.rows.begin()] = point.phi[j];
        }
    };
    if (active.empty())
        shape(0., inward);
    const auto sampleRow=SampleKRow({1.57,0.,0.});
    for(int j=0;j<sampleRow.ids.size();++j)
    {
        const auto id=offset[K()]+sampleRow.ids[j];
        const auto col=std::lower_bound(frozen.rows.begin(),frozen.rows.end(),id);
        if(col!=frozen.rows.end() && *col==id)sampled[col-frozen.rows.begin()]=sampleRow.phi[j];
    }
    shape(.04, normalization);
    std::vector<std::pair<std::string, Vector>> directions{
        {"ground_inward", inward},
        {"sampled_K", sampled},
        {"normalization_ray", normalization}};
    for (int i = 1; i <= 3; ++i)
    {
        Vector tangent(inward.size());
        for (int j = 0; j < tangent.size(); ++j)
            tangent[j] = std::sin((.31 + i * .67) * frozen.rows[j]);
        for (int pass = 0; pass < 3; ++pass)
            for (const auto &row : active)
            {
                const double component = dot(tangent, row) / dot(row, row);
                for (int j = 0; j < tangent.size(); ++j)
                    tangent[j] -= component * row[j];
            }
        for (int j = 0; j < tangent.size(); ++j)
            tangent[j] += .01 * inward[j];
        directions.emplace_back("near_tangent_" + std::to_string(i),
                                std::move(tangent));
    }
    std::vector<std::string> actions;
    for (auto &[name, direction] : directions)
    {
        const double amplitude = Maximum(direction);
        Require(amplitude > 0., "empty conditional k verification direction");
        for (auto &entry : direction)
            entry /= amplitude;
        const auto action = Multiply(exact.matrix, direction);
        auto difference   = Multiply(frozen.matrix, direction);
        for (int j = 0; j < difference.size(); ++j)
            difference[j] -= action[j];
        std::vector<std::string> samples;
        double best = 1e300;
        for (double h : {1e-4, 1e-5, 1e-6, 1e-7, 1e-8})
        {
            ScalarDomainResult plusDomain, minusDomain;
            const auto plus  = evaluate(trial(direction, h), plusDomain);
            const auto minus = evaluate(trial(direction, -h), minusDomain);
            if (plus.empty() && minus.empty())
                continue;
            const auto &lo     = minus.empty() ? frozen.residual : minus;
            const auto &hi     = plus.empty() ? frozen.residual : plus;
            const double width = !plus.empty() && !minus.empty() ? 2 * h : h;
            Vector error(direction.size());
            for (int j = 0; j < error.size(); ++j)
                error[j] = (hi[j] - lo[j]) / width - action[j];
            const double relative = norm(error) / std::max(1., norm(action));
            best                  = std::min(best, relative);
            const auto &observed  = plus.empty() ? minusDomain : plusDomain;
            samples.push_back(J::Object(
                {{"h", J::Number(h)},
                 {"scheme", J::Quote(width == 2 * h ? "central"
                                     : plus.empty() ? "backward"
                                                    : "forward")},
                 {"relative_error", J::Number(relative)},
                 {"b_change", J::Number(observed.minimumB - base.minimumB)},
                 {"sampled_minimum_change",
                  J::Number(observed.minimumSampledK - base.minimumSampledK)},
                 {"profile_branch_change",
                  J::Boolean(observed.profileBranches !=
                             base.profileBranches)}}));
        }
        actions.push_back(
            J::Object({{"name", J::Quote(name)},
                       {"direction", J::Numbers(direction)},
                       {"action_norm", J::Number(norm(action))},
                       {"frozen_difference", J::Number(norm(difference))},
                       {"best_relative_error", J::Number(best)},
                       {"samples", J::List(samples)}}));
    }
    auto projection = [&](const PicardBlock &block) {
        const Factor factor(block.matrix);
        auto rhs = block.residual;
        for (auto &entry : rhs)
            entry = -entry;
        const auto raw       = factor.Apply(rhs);
        int cuts             = 0;
        std::string kkt      = "null";
        const auto projected = ProjectPicardGround(
            state, block, factor, raw, m_picardRelaxation[1], cuts, &kkt);
        std::vector<std::string> samples;
        if (!projected.empty())
        {
            const auto change = Multiply(block.matrix, projected);
            for (double fraction : {1., .125, .015625})
            {
                ScalarDomainResult domain;
                const auto actual =
                    evaluate(trial(projected, fraction), domain);
                auto predicted = block.residual;
                for (int i = 0; i < predicted.size(); ++i)
                    predicted[i] += fraction * change[i];
                samples.push_back(J::Object(
                    {{"fraction", J::Number(fraction)},
                     {"predicted_k_norm", J::Number(norm(predicted))},
                     {"actual_k_norm",
                      actual.empty() ? "null" : J::Number(norm(actual))}}));
            }
        }
        return J::Object({{"unprojected", J::Numbers(raw)},
                          {"affine_ground_extrema", extrema(trial(raw, 1.))},
                          {"projected", J::Numbers(projected)},
                          {"cuts", J::Number(cuts)},
                          {"kkt", kkt},
                          {"trials", J::List(samples)}});
    };
    const auto qualified =
        CompleteRestrictedAudit(state, history, alpha, time, data, true);
    std::ofstream output(m_picardKDiagnostic);
    output << J::Object(
                  {{"alpha", J::Number(alpha)},
                   {"time", J::Number(time)},
                   {"guess_sha256", J::Quote(RANSWallR3::Hash(state))},
                   {"history_sha256", J::Quote(RANSWallR3::Hash(history))},
                   {"outer_range", J::Numbers(qualified.at("r2_outer_range"))},
                   {"outer_gradient_witness",
                    J::Numbers(qualified.at("r4_outer_gradient_witness"))},
                   {"k_norm", J::Number(norm(frozen.residual))},
                   {"free_rows", J::Numbers(frozen.rows)},
                   {"k_residual", J::Numbers(frozen.residual)},
                   {"mass", J::Rows(mass)},
                   {"frozen_matrix", J::Rows(frozen.matrix)},
                   {"k_jacobian", J::Rows(exact.matrix)},
                   {"ground_extrema", extrema(state)},
                   {"actions", J::List(actions)},
                   {"frozen_projection", projection(frozen)},
                   {"true_projection", projection(exact)}})
           << '\n';
    Require(output.good(), "cannot write Picard k diagnostic");
}

RANSWallEvolutionPatch::Vector RANSWallEvolutionPatch::SolveSegregatedPicard(
    Vector s, const Vector &history, double alpha, double time,
    const RANSWallEvolutionData &data, int &iterations, double tolerance,
    int maximumIterations) const
{
    Require(s.size() == Size() && history.size() == s.size() &&
                std::isfinite(alpha) && alpha > 0. && std::isfinite(time) &&
                std::isfinite(tolerance) && tolerance > 0. &&
                maximumIterations >= 0,
            "invalid Picard physical stage");
    iterations = 0;
    m_lastFeasibleState.clear();
    const auto incomingGuessHash = RANSWallR3::Hash(s);
    if (!m_picardGuessFile.empty())
    {
        boost::property_tree::ptree packet;
        boost::property_tree::read_json(m_picardGuessFile, packet);
        auto read = [&](const char *name) {
            Vector values;
            for (const auto &entry : packet.get_child(name))
                values.push_back(entry.second.get_value<double>());
            return values;
        };
        auto guess = read("state");
        Require(read("history") == history && guess.size() == s.size() &&
                    packet.get<double>("alpha") == alpha &&
                    packet.get<double>("time") == time,
                "saved Picard guess differs from native stage or history");
        Require(ExternalLift(guess, time, data) == guess,
                "saved Picard guess changes boundary data");
        s = std::move(guess);
    }
    s             = ExternalLift(s, time, data);
    auto external = ExternalConstraints(time, data);
    for (int id : m_traceIds)
        external.erase(offset[Epsilon()] + id);
    Require(external.empty() || external.rbegin()->first < DifferentialSize(),
            "Picard external lifting may constrain differential rows only");
    const auto start            = std::chrono::steady_clock::now();
    const bool nonK             = m_picardNonKTolerance > 0.;
    const double sweepTolerance = nonK ? m_picardNonKTolerance : tolerance;
    auto sweepNorms             = [&](Vector values) {
        if (nonK)
            values.at(K()) = 0.;
        return values;
    };
    int flowSolves = 0, kSolves = 0, epsilonSolves = 0;
    auto emit = [&](const std::string &event,
                    std::vector<std::pair<std::string, std::string>> values) {
        values.emplace_back("event", J::Quote(event));
        values.emplace_back("sweep", J::Number(iterations));
        values.emplace_back("alpha", J::Number(alpha));
        values.emplace_back("time", J::Number(time));
        values.emplace_back(
            "seconds", J::Number(std::chrono::duration<double>(
                                     std::chrono::steady_clock::now() - start)
                                     .count()));
        std::cout << "EP-AWF Picard " << J::Object(values) << std::endl;
    };
    auto admissible = [&](const Vector &state, ScalarDomainResult &domain,
                          std::string &ground) {
        domain     = ScalarDomain(state, time, data);
        bool valid = false;
        if (domain.admissible)
            ground = WallTraceJSON(domain.effectiveState, valid);
        return valid;
    };
    auto normsAt = [&](const Vector &state) {
        const auto residual = Residual(state, history, alpha, time, data);
        return Norms(residual, Constraints(state, time, data), m_robustUstar);
    };
    auto trialNorms = [&](const Vector &state, double fraction) {
        try
        {
            return normsAt(state);
        }
        catch (const std::exception &error)
        {
            emit("rejected_residual", {{"fraction", J::Number(fraction)},
                                       {"reason", J::Quote(error.what())}});
            return Vector{};
        }
    };
    auto blockCheck = [&](const PicardBlock &system, const Vector &state,
                          const char *name) {
        const auto native = Residual(state, history, alpha, time, data);
        Vector error(Size());
        for (int i = 0; i < system.rows.size(); ++i)
            error[system.rows[i]] = system.residual[i] - native[system.rows[i]];
        const auto constraints = Constraints(state, time, data);
        const auto norms       = Norms(error, constraints, m_robustUstar);
        emit("fixed_point", {{"block", J::Quote(name)},
                             {"scaled_difference", J::Numbers(norms)}});
        Require(
            Maximum(norms) <=
                1e-10 * std::max(1., Maximum(Norms(native, constraints,
                                                   m_robustUstar))),
            "frozen Picard block differs from the native same-state residual");
    };
    auto apply = [](Vector &state, const PicardBlock &block, const Vector &step,
                    double relaxation) {
        if(step.size() != block.rows.size())
            throw std::logic_error("Picard block increment shape mismatch");
        for (int i = 0; i < block.rows.size(); ++i)
            state.at(block.rows[i]) += relaxation * step[i];
    };
    auto scalarStep = [&](int b, Vector &state, double requested,
                          ScalarDomainResult &domain, std::string &ground) {
        const char *name = b == K() ? "k" : "epsilon_transport_trace";
        const auto block =
            AssemblePicardBlock(b, state, history, alpha, time, data, external);
        if (m_picardVerify && iterations == 0)
            blockCheck(block, state, name);
        if (b == K() && m_picardKNewton)
        {
            auto model = block;
            model.residual =
                PicardKResidual(state, history, alpha, time, data, model.rows);
            model.matrix =
                PicardKJacobian(state, alpha, time, data, model.rows);
            const Factor jacobian(model.matrix), mass(PicardKMass(model.rows));
            const double base = DualNorm(mass, model.residual);
            auto rhs          = model.residual;
            for (auto &value : rhs)
                value = -value;
            const auto newton   = jacobian.Apply(rhs);
            const double radius = m_picardKRadius * m_parameters.model.Q() *
                                  m_robustUstar * m_robustUstar;
            double fraction =
                std::min(requested, radius / std::max(radius, Maximum(newton)));
            ++kSolves;
            for (int retry = 0; retry < 9; ++retry, fraction *= .5)
            {
                int cuts        = 0;
                const auto step = ProjectPicardGround(state, model, jacobian,
                                                      newton, fraction, cuts);
                if (step.empty())
                    continue;
                auto predicted = Multiply(model.matrix, step);
                for (int i = 0; i < predicted.size(); ++i)
                    predicted[i] += model.residual[i];
                const double expected = DualNorm(mass, predicted);
                if (expected >= base)
                    continue;
                auto trial = state;
                apply(trial, model, step, 1.);
                if (!admissible(trial, domain, ground))
                    continue;
                const auto residual = PicardKResidual(trial, history, alpha,
                                                      time, data, model.rows);
                const double actual = DualNorm(mass, residual);
                const bool accept   = actual < base - 1e-4 * (base - expected);
                emit("conditional_k_trial",
                     {{"base", J::Number(base)},
                      {"predicted", J::Number(expected)},
                      {"actual", J::Number(actual)},
                      {"fraction", J::Number(fraction)},
                      {"cuts", J::Number(cuts)},
                      {"step_max", J::Number(Maximum(step))},
                      {"accepted", J::Boolean(accept)},
                      {"dimension", J::Number(step.size())}});
                if (accept)
                {
                    state = std::move(trial);
                    return fraction;
                }
            }
            emit("conditional_k_fallback", {{"base", J::Number(base)}});
        }
        auto rhs = block.residual;
        for (double &entry : rhs)
            entry = -entry;
        const Factor factor(block.matrix);
        const auto step = factor.Apply(rhs);
        auto defect     = Multiply(block.matrix, step);
        for (int i = 0; i < defect.size(); ++i)
            defect[i] -= rhs[i];
        Require(Maximum(defect) <= 1e-9 * std::max(1., Maximum(rhs)),
                "inaccurate Picard scalar block solve");
        ++(b == K() ? kSolves : epsilonSolves);
        emit("solve", {{"block", J::Quote(name)},
                       {"dimension", J::Number(step.size())},
                       {"frozen_scalar_defect", J::Number(Maximum(defect))}});
        for (int i = 0; i < 14; ++i, requested *= .5)
        {
            auto trial = state;
            if (b == K() && m_picardGroundProjection)
            {
                int cuts             = 0;
                const auto projected = ProjectPicardGround(
                    state, block, factor, step, requested, cuts);
                if (cuts)
                    emit("ground_projection",
                         {{"cuts", J::Number(cuts)},
                          {"relaxation", J::Number(requested)},
                          {"solved", J::Boolean(!projected.empty())}});
                if (projected.empty())
                    continue;
                apply(trial, block, projected, 1.);
            }
            else
                apply(trial, block, step, requested);
            if (admissible(trial, domain, ground))
            {
                state = std::move(trial);
                return requested;
            }
        }
        emit("rejected_scalar",
             {{"block", J::Quote(name)}, {"domain", ScalarDomainJSON(domain)}});
        return 0.;
    };
    bool outerViolationSeen = false;
    auto report             = [&](const char *event, const Vector &norms,
                      const ScalarDomainResult &domain,
                      const std::string &ground, const Vector &relaxation) {
        double traceDefect = 0.;
        for (int id : m_traceIds)
        {
            const int row = offset[Epsilon()] + id;
            traceDefect   = std::max(traceDefect,
                                                 std::abs(s[row] - domain.traceState[row]));
        }
        emit(event, {{"norms", J::Numbers(norms)},
                     {"state_sha256", J::Quote(RANSWallR3::Hash(s))},
                     {"relaxation", J::Numbers(relaxation)},
                     {"scalar_domain", ScalarDomainJSON(domain)},
                     {"ground", ground},
                     {"trace_coefficient_defect", J::Number(traceDefect)},
                     {"solve_counts",
                      J::Numbers(std::array<int, 4>{flowSolves, flowSolves,
                                                    kSolves, epsilonSolves})}});
        if ((m_picardKNewton || nonK) &&
            (!outerViolationSeen || iterations % 10 == 0))
        {
            const auto audit =
                CompleteRestrictedAudit(s, history, alpha, time, data, true);
            const auto &range = audit.at("r2_outer_range");
            outerViolationSeen |= range[0] > .001;
            emit(
                "qualification",
                {{"outer_range", J::Numbers(range)},
                 {"outer_gradient_witness",
                  J::Numbers(audit.at("r4_outer_gradient_witness"))},
                 {"qualified", J::Boolean(audit.at("r2_qualified")[0] == 1.)}});
        }
    };
    auto perturb = [&](Vector &guess, double amplitude) {
        const std::array<double, 4> scales{
            {m_robustUstar, m_robustUstar,
             m_parameters.model.Q() * m_robustUstar * m_robustUstar,
             std::pow(m_robustUstar, 3) /
                 (m_parameters.model.Kappa() * m_parameters.geometry.height)}};
        for (int b = 0; b < InnerPressure(); ++b)
            for (int i = 0; i < size[b]; ++i)
            {
                const int row = offset[b] + i;
                if (!external.count(row))
                    guess[row] += amplitude * scales[b] * std::sin(1. + row);
            }
    };

    ScalarDomainResult domain;
    std::string ground;
    if (m_picardGuessPerturbation)
        perturb(s, m_picardGuessPerturbation);
    Require(admissible(s, domain, ground),
            "Picard initial state fails scalar or ground domain: " +
                domain.failure);
    m_lastFeasibleState = s;
    auto norms          = normsAt(s);
    auto historyValues  = history;
    for (auto &value : historyValues)
        if (value == 0.)
            value = 0.;
    emit("stage_input",
         {{"unknowns", J::Number(Size())},
          {"tolerance", J::Number(tolerance)},
          {"history_sha256", J::Quote(RANSWallR3::Hash(history))},
          {"history_value_sha256", J::Quote(RANSWallR3::Hash(historyValues))},
          {"incoming_guess_sha256", J::Quote(incomingGuessHash)},
          {"guess_sha256", J::Quote(RANSWallR3::Hash(s))}});
    report("initial", norms, domain, ground, {});
    if (m_picardReduced)
        return SolveReducedPicard(std::move(s), history, alpha, time, data,
                                  iterations, tolerance, maximumIterations);
    if (!m_picardKDiagnostic.empty())
        DiagnosePicardK(s, history, alpha, time, data, external);
    if (m_picardVerify)
    {
        auto perturbed = s;
        perturb(perturbed, 1e-7);
        ScalarDomainResult perturbedDomain;
        std::string perturbedGround;
        Require(admissible(perturbed, perturbedDomain, perturbedGround),
                "Picard consistency perturbation leaves the scalar or ground "
                "domain");
        for (const int b : {0, K(), Epsilon()})
        {
            const auto baseName = b == 0     ? "initial_flow"
                                  : b == K() ? "initial_k"
                                             : "initial_epsilon";
            blockCheck(
                AssemblePicardBlock(b, s, history, alpha, time, data, external),
                s, baseName);
            const auto name = b == 0     ? "perturbed_flow"
                              : b == K() ? "perturbed_k"
                                         : "perturbed_epsilon";
            blockCheck(AssemblePicardBlock(b, perturbed, history, alpha, time,
                                           data, external),
                       perturbed, name);
        }
    }
    const Vector fixedK(s.begin() + offset[K()], s.begin() + offset[Epsilon()]);
    auto flowDirection = [&](const Vector &base, bool verify) {
        auto correction =
            PicardFlowCorrection(base, history, alpha, time, data, external);
        auto flow    = std::move(correction.block);
        const int nv = correction.velocitySize, np = correction.pressureSize;
        const auto &predictor = correction.predictor;
        auto step             = std::move(correction.step);
        ++flowSolves;
        emit("solve", {{"block", J::Quote("velocity_predictor")},
                       {"dimension", J::Number(nv)}});
        emit("solve",
             {{"block", J::Quote("pressure_correction")},
              {"dimension", J::Number(np)},
              {"schur_velocity_rhs", J::Number(np)},
              {"predictor_continuity_gauge_max",
               J::Number(correction.predictorDefect)},
              {"frozen_flow_defect", J::Number(Maximum(correction.defect))}});
        if (verify)
        {
            auto directRhs = flow.residual;
            for (double &entry : directRhs)
                entry = -entry;
            const auto direct = RANSWallAlgebra::Solve(flow.matrix, directRhs);
            Vector difference(step.size());
            double velocityCorrection = 0.;
            for (int i = 0; i < step.size(); ++i)
                difference[i] = step[i] - direct[i];
            for (int i = 0; i < nv; ++i)
                velocityCorrection = std::max(velocityCorrection,
                                              std::abs(step[i] - predictor[i]));
            auto predicted = base, corrected = base;
            Require(predictor.size() == nv,
                    "Picard velocity predictor shape mismatch");
            for (int i = 0; i < nv; ++i)
                predicted.at(flow.rows[i]) += predictor.at(i);
            apply(corrected, flow, step, 1.);
            auto pressureDirection = Vector(Size());
            for (int j = 0; j < np; ++j)
                pressureDirection[flow.rows[nv + j]] = step[nv + j];
            const auto nativePressure =
                AffinePressureAction(pressureDirection, external);
            Vector assembledPressure(nv);
            for (int i = 0; i < nv; ++i)
                for (int j = 0; j < np; ++j)
                    assembledPressure[i] +=
                        flow.matrix[i][nv + j] * step[nv + j];
            Vector pressureError(nv);
            for (int i = 0; i < nv; ++i)
                pressureError[i] =
                    assembledPressure[i] - nativePressure[flow.rows[i]];
            const auto predictedNorms = normsAt(predicted),
                       correctedNorms = normsAt(corrected);
            auto directState          = base;
            apply(directState, flow, direct, 1.);
            const auto correctedInterface =
                Audit(corrected, history, alpha, time, data, m_robustUstar)
                    .at("interface");
            const auto directInterface =
                Audit(directState, history, alpha, time, data, m_robustUstar)
                    .at("interface");
            Vector interfaceError(correctedInterface.size());
            for (int i = 0; i < interfaceError.size(); ++i)
                interfaceError[i] = correctedInterface[i] - directInterface[i];
            emit("frozen_saddle_comparison",
                 {{"coefficient_difference", J::Number(Maximum(difference))},
                  {"velocity_correction", J::Number(velocityCorrection)},
                  {"affine_pressure_difference",
                   J::Number(Maximum(pressureError))},
                  {"interface_corrected", J::Numbers(correctedInterface)},
                  {"interface_direct", J::Numbers(directInterface)},
                  {"original_predictor_continuity_gauge",
                   J::Numbers(Vector(predictedNorms.begin() + InnerPressure(),
                                     predictedNorms.begin() + Gauge() + 1))},
                  {"original_corrected_continuity_gauge",
                   J::Numbers(Vector(correctedNorms.begin() + InnerPressure(),
                                     correctedNorms.begin() + Gauge() + 1))}});
            Require(Maximum(difference) <= 1e-8 * std::max(1., Maximum(direct)),
                    "Picard pressure correction disagrees with the direct "
                    "frozen saddle system");
            Require(Maximum(pressureError) <=
                        1e-10 * std::max(1., Maximum(assembledPressure)),
                    "Picard affine pressure increment disagrees with the "
                    "native action");
            Require(Maximum(interfaceError) <=
                        1e-8 * std::max(1., Maximum(directInterface)),
                    "Picard interface observations disagree with the direct "
                    "frozen saddle state");
        }

        return std::make_pair(std::move(flow), std::move(step));
    };

    auto meritOf = [](const Vector &values) {
        return std::inner_product(values.begin(), values.end(), values.begin(),
                                  0.);
    };
    auto equilibrate = [&](Vector &state, Vector &values,
                           ScalarDomainResult &domain, std::string &ground) {
        auto h = [&](Vector n) {
            n.at(K()) = 0.;
            return n;
        };
        const auto before   = values;
        const double target = std::max(1e-11, 1e-3 * Maximum(h(values)));
        const Vector held(state.begin() + offset[K()],
                          state.begin() + offset[Epsilon()]);
        int updates = 0;
        while (updates < 12 && Maximum(h(values)) > target)
        {
            const auto [flow, step] = flowDirection(state, false);
            bool improved           = false;
            for (double fraction : {1., .5, .25, .125})
            {
                auto candidate = state;
                apply(candidate, flow, step, fraction);
                ScalarDomainResult candidateDomain;
                std::string candidateGround;
                if (!admissible(candidate, candidateDomain, candidateGround) ||
                    !scalarStep(Epsilon(), candidate, fraction, candidateDomain,
                                candidateGround))
                    continue;
                auto next = trialNorms(candidate, fraction);
                if (next.empty() || meritOf(h(next)) >= meritOf(h(values)))
                    continue;
                Require(std::equal(held.begin(), held.end(),
                                   candidate.begin() + offset[K()]),
                        "reduced Picard trial changed fixed k during "
                        "equilibration");
                state  = std::move(candidate);
                values = std::move(next);
                domain = std::move(candidateDomain);
                ground = std::move(candidateGround);
                ++updates;
                improved = true;
                break;
            }
            if (!improved)
                break;
        }
        emit("reduced_k_trial", {{"before_equilibration", J::Numbers(before)},
                                 {"after_equilibration", J::Numbers(values)},
                                 {"non_k_target", J::Number(target)},
                                 {"non_k_updates", J::Number(updates)}});
    };

    while (Maximum(sweepNorms(norms)) > sweepTolerance &&
           iterations < maximumIterations)
    {
        const auto [flow, step] =
            flowDirection(s, m_picardVerify && iterations == 0);
        const auto governingNorms = sweepNorms(norms);
        const double merit =
            std::inner_product(governingNorms.begin(), governingNorms.end(),
                               governingNorms.begin(), 0.);
        bool accepted   = false;
        double fraction = 1.;
        for (int retry = 0; retry < 12; ++retry, fraction *= .5)
        {
            auto trial = s;
            ScalarDomainResult trialDomain;
            std::string trialGround;
            const double flowRelax = m_picardRelaxation[0] * fraction;
            apply(trial, flow, step, flowRelax);
            if (!admissible(trial, trialDomain, trialGround))
            {
                emit("rejected_flow",
                     {{"fraction", J::Number(fraction)},
                      {"domain", ScalarDomainJSON(trialDomain)}});
                continue;
            }
            const auto flowNorms = trialNorms(trial, fraction);
            if (flowNorms.empty())
                continue;
            const double kRelax =
                nonK ? 0.
                     : scalarStep(K(), trial, m_picardRelaxation[1] * fraction,
                                  trialDomain, trialGround);
            if (!nonK && !kRelax)
                continue;
            // The k update may leave a trace mismatch. Assemble epsilon
            // with the new k and freshly evaluated mechanical wall forcing
            // and sensors. No other unknown is lifted in this scalar solve.
            const double epsRelax =
                scalarStep(Epsilon(), trial, m_picardRelaxation[2] * fraction,
                           trialDomain, trialGround);
            if (!epsRelax)
                continue;
            auto next = trialNorms(trial, fraction);
            if (next.empty())
                continue;
            const Vector relaxation{flowRelax, kRelax, epsRelax};
            emit("trial", {{"norms", J::Numbers(next)},
                           {"relaxation", J::Numbers(relaxation)},
                           {"post_flow_continuity_gauge",
                            J::Numbers(Vector{flowNorms[InnerPressure()],
                                              flowNorms[OuterPressure()],
                                              flowNorms[Gauge()]})}});
            if (m_picardEquilibrateK && retry == 0 && meritOf(next) >= merit)
                equilibrate(trial, next, trialDomain, trialGround);
            const auto governingNext = sweepNorms(next);
            if (Maximum(governingNext) > sweepTolerance &&
                std::inner_product(governingNext.begin(), governingNext.end(),
                                   governingNext.begin(), 0.) >= merit)
                continue;
            if (nonK)
                Require(std::equal(fixedK.begin(), fixedK.end(),
                                   trial.begin() + offset[K()]),
                        "non-k equilibration changed a fixed k coefficient");
            s                   = std::move(trial);
            m_lastFeasibleState = s;
            norms               = next;
            ++iterations;
            report("accepted_sweep", norms, trialDomain, trialGround,
                   relaxation);
            accepted = true;
            break;
        }
        Require(
            accepted,
            "Picard could not find an admissible residual-decreasing sweep");
    }
    if (nonK)
    {
        Require(Maximum(sweepNorms(norms)) <= sweepTolerance,
                "Picard diagnostic exhausted the non-k equilibration budget");
        std::vector<int> rows;
        for (int row = offset[K()]; row < offset[Epsilon()]; ++row)
            if (!external.count(row))
                rows.push_back(row);
        const auto audit =
            CompleteRestrictedAudit(s, history, alpha, time, data, true);
        emit("non_k_equilibrated",
             {{"norms", J::Numbers(norms)},
              {"non_k_norm", J::Number(Maximum(sweepNorms(norms)))},
              {"k_residual", J::Numbers(PicardKResidual(s, history, alpha, time,
                                                        data, rows))},
              {"outer_range", J::Numbers(audit.at("r2_outer_range"))},
              {"outer_gradient_witness",
               J::Numbers(audit.at("r4_outer_gradient_witness"))}});
        throw std::runtime_error("Picard diagnostic completed non-k "
                                 "equilibration without advancing history");
    }
    Require(Maximum(norms) <= tolerance,
            "Picard sweep limit before full stage convergence");
    emit("converged", {{"norms", J::Numbers(norms)}});
    return s;
}
} // namespace Nektar

namespace Nektar
{
std::string RANSWallEvolutionPatch::WJNativeChecks(const Vector &state,const RANSWallEvolutionData &data) const
{
    Require(WJ() && !Mapped() && m_dimension==2 && !CappedDestruction(),"WJ checks require the selected native pilot");
    const auto hash=RANSWallR3::Hash(state);const double alpha=1e6,time=0.;
    const auto constraints=Constraints(state,time,data),external=ExternalConstraints(time,data);
    const auto physical=Residual(state,state,alpha,time,data,false);
    double blockError=0,pressureColumnError=0,pressureFDError=0,traceColumn=0;
    std::vector<std::string> blockReports,pressureReports;
    for(int block:{0,K(),Epsilon()})
    {
        const auto frozen=AssemblePicardBlock(block,state,state,alpha,time,data,external);
        double error=0;
        for(int i=0;i<frozen.rows.size();++i)error=std::max(error,std::abs(frozen.residual[i]-physical[frozen.rows[i]]));
        blockError=std::max(blockError,error);
        blockReports.push_back(J::Object({{"block",J::Number(block)},{"same_state_max_defect",J::Number(error)},
            {"rows",J::Numbers(frozen.rows)},{"frozen_residual",J::Numbers(frozen.residual)}}));
        if(block==0)
        {
            for(int mode=0;mode<5;++mode)
            {
                Vector direction(Size());
                if(mode==2)direction[offset[Gauge()]]=1.;
                else if(mode==4)
                {
                    // Nonconstant outer pressure produces a nonzero sampled
                    // tangential roughness forcing with all sensors absent.
                    const auto values=fields[OuterPressure()]->Project([](const RANSWallCoordinate &p){return p[0]+.13*p[0]*p[2];});
                    std::copy(values.begin(),values.end(),direction.begin()+offset[OuterPressure()]);
                }
                else for(int side=0;side<2;++side)for(int i=0;i<size[InnerPressure()+side];++i)
                    if((mode==0 && side==0) || (mode==1 && side==1) || mode==3)
                        direction[offset[InnerPressure()+side]+i]=(mode==3 && side==1?-1:1)*m_constants[InnerPressure()+side][i];
                const auto action=AffinePressureAction(direction,constraints);
                Vector reduced;for(int row:frozen.rows)reduced.push_back(direction[row]);
                const auto column=RANSWallAlgebra::Multiply(frozen.matrix,reduced);
                double ce=0,fd=1e300;
                for(int i=0;i<frozen.rows.size();++i)ce=std::max(ce,std::abs(column[i]-action[frozen.rows[i]]));
                for(double h:{1e-3,1e-4,1e-5})
                {
                    auto plus=state,minus=state;
                    for(int i=0;i<Size();++i){plus[i]+=h*direction[i];minus[i]-=h*direction[i];}
                    const auto rp=Residual(plus,state,alpha,time,data,false),rm=Residual(minus,state,alpha,time,data,false);
                    double errorFD=0;for(int i=0;i<Size();++i)errorFD=std::max(errorFD,std::abs((rp[i]-rm[i])/(2*h)-action[i]));
                    fd=std::min(fd,errorFD);
                }
                pressureColumnError=std::max(pressureColumnError,ce);pressureFDError=std::max(pressureFDError,fd);
                pressureReports.push_back(J::Object({{"mode",J::Number(mode)},{"frozen_column_defect",J::Number(ce)},{"physical_difference_defect",J::Number(fd)}}));
            }
        }
        if(block==Epsilon())
        {
            for(int id:m_traceIds)
            {
                const auto found=std::find(frozen.rows.begin(),frozen.rows.end(),offset[Epsilon()]+id);
                Require(found!=frozen.rows.end(),"WJ epsilon trace row lost");
                const int col=found-frozen.rows.begin();
                for(int i=0;i<frozen.rows.size();++i)if(i!=col)traceColumn=std::max(traceColumn,std::abs(frozen.matrix[i][col]));
            }
        }
    }
    Require(blockError<1e-8 && pressureColumnError<1e-10 && pressureFDError<1e-7 && traceColumn>0,
        "WJ frozen same-state/pressure/trace check failed");
    // The resting-velocity fixture keeps the same transported scalars and
    // sampled K. Complete prescribed destruction must remain in the residual.
    auto resting=state;std::fill(resting.begin(),resting.begin()+offset[K()],0.);
    auto restingData=data;for(auto &[row,value]:restingData.fixedBoundary)if(row<offset[K()])value=0;
    const auto restingResidual=Residual(resting,resting,0,time,restingData,false);
    const auto restingFrozen=AssemblePicardBlock(K(),resting,resting,0,time,restingData,ExternalConstraints(time,restingData));
    double sinkBlockError=0,sinkPointError=0;
    for(int i=0;i<restingFrozen.rows.size();++i)sinkBlockError=std::max(sinkBlockError,std::abs(restingFrozen.residual[i]-restingResidual[restingFrozen.rows[i]]));
    Refresh(resting,time,&restingData);
    for(const auto &ray:m_rays)for(const auto &q:ray.points)
    {
        const auto source=WallSource(resting,q.v.coordinate,0,*ray.profile);
        const double expected=WJEP()?-ray.profile->ReferenceProduction(Distance(q.v.coordinate)):-ray.profile->Dissipation(Distance(q.v.coordinate));
        sinkPointError=std::max(sinkPointError,std::abs(source.kSource-expected));
    }
    Require(sinkBlockError<1e-9 && sinkPointError<1e-12,"uncapped WJ sink was dropped");
    const auto resolution=ResidualResolutionComponents(state,state,0,time,data,false);
    std::vector<std::pair<std::string,std::string>> terms;
    for(const auto &[name,values]:resolution.terms)terms.push_back({name,J::Numbers(values)});
    // Independent test quadrature: new higher-order Gauss rules evaluate the
    // unchanged represented fields, split at physical native/source knots.
    // These evaluations perform no projection, Lift or governing solve.
    std::vector<std::string> samples;
    auto append=[&](const VolumePoint &q,bool inner,const EquilibriumAWFProfiles *profile)
    {
        const auto &v=q.v;RANSTensor g{};for(int a=0;a<m_dimension;++a)g[PhysicalAxis(a)]=v.Gradient(state,offset[a]);
        Vector gradient;for(const auto &row:g)gradient.insert(gradient.end(),row.begin(),row.end());
        const double k=v.Value(state,offset[K()]),epsilon=inner?profile->Dissipation(v.coordinate[2]):q.e.Value(state,offset[Epsilon()]);
        const auto gk=v.Gradient(state,offset[K()]);RANSVector ge{};if(!inner)ge=q.e.Gradient(state,offset[Epsilon()]);
        std::vector<std::string> derivatives;for(int axis:{0,2})derivatives.push_back(J::Numbers(v.derivative[axis]));
        std::vector<std::string> epsilonDerivatives;if(!inner)for(int axis:{0,2})epsilonDerivatives.push_back(J::Numbers(q.e.derivative[axis]));
        samples.push_back(J::Object({{"inner",J::Boolean(inner)},{"geometry",J::Number(v.geometry)},{"element",J::Number(v.element)},
            {"coordinate",J::Numbers(v.coordinate)},{"weight",J::Number(q.weight)},{"k",J::Number(k)},{"epsilon",J::Number(epsilon)},
            {"G",J::Numbers(gradient)},{"grad_k",J::Numbers(gk)},{"grad_epsilon",J::Numbers(ge)},
            {"sampled_K",inner?J::Number(profile->K()):"null"},{"ids",J::Numbers(v.ids)},{"phi",J::Numbers(v.phi)},
            {"derivative_xz",J::List(derivatives)},{"epsilon_ids",inner?"[]":J::Numbers(q.e.ids)},
            {"epsilon_phi",inner?"[]":J::Numbers(q.e.phi)},{"epsilon_derivative_xz",J::List(epsilonDerivatives)}}));
    };
    const auto xrule=LibUtilities::PointsManager()[LibUtilities::PointsKey(10,LibUtilities::eGaussGaussLegendre)];
    for(const auto &element:innerVelocity->Elements())if(element.lower[2]==0)
        for(int j=0;j<10;++j)
        {
            RANSWallCoordinate anchor{element.lower[0]+(xrule->GetZ()[j]+1)*(element.upper[0]-element.lower[0])/2,0,0};
            const auto selected=StateProfile(state,anchor,time,data);auto options=selected->Options();options.quadratureOrder=32;
            EquilibriumAWFProfiles profile(selected->Geometry(),selected->K(),m_parameters.model,{},options,m_wallParameters);
            const double wx=xrule->GetW()[j]*(element.upper[0]-element.lower[0])/2;
            for(const auto &r:profile.Quadrature()){auto p=anchor;p[2]=r.r;append({Point(0,p),Point(InnerPressure(),p),{},wx*r.weight},true,&profile);}
        }
    for(const auto &v:outerVelocity->Quadrature(24))append({v,Point(OuterPressure(),v.coordinate),Point(Epsilon(),v.coordinate),v.weight},false,nullptr);
    // Extract actual native off-constraint states and test the common action
    // in a simultaneous k/epsilon/raw-G direction without changing the state.
    std::vector<std::string> constitutive;int extracted=0;
    for(const auto &q:m_outer)
    {
        if(extracted++>=4)break;
        RANSTensor g{},h{{{.2,.1,.3},{0,-.1,0},{-.2,0,.4}}};
        for(int a=0;a<m_dimension;++a)g[PhysicalAxis(a)]=q.v.Gradient(state,offset[a]);
        g[0][0]+=.07;g[2][0]-=.13;
        const double k=q.v.Value(state,offset[K()]),epsilon=q.e.Value(state,offset[Epsilon()]);
        const auto c=OuterClosure(k,epsilon,g,q.v.coordinate),d=ApplyWJClosureDirection(k,epsilon,g,.13,-.04,h,m_parameters.wj,WJGradientPolicy::DevStrainRawWorkV1);
        Vector G,H,R,dR;for(int a=0;a<3;++a)for(int b=0;b<3;++b){G.push_back(g[a][b]);H.push_back(h[a][b]);R.push_back(c.reynoldsStress[a][b]);dR.push_back(d.reynoldsStress[a][b]);}
        constitutive.push_back(J::Object({{"coordinate",J::Numbers(q.v.coordinate)},{"k",J::Number(k)},{"epsilon",J::Number(epsilon)},
            {"G",J::Numbers(G)},{"direction_G",J::Numbers(H)},{"R",J::Numbers(R)},{"direction_R",J::Numbers(dR)},
            {"nu_eff",J::Number(c.eddyViscosity)},{"nu_d",J::Number(c.scalarViscosity)},{"Praw",J::Number(c.production)},
            {"Pdev",J::Number(c.productionDev)},{"divergence",J::Number(c.divergence)},
            {"direction_sources_diffusion",J::Numbers(Vector{d.kSource,d.epsilonSource,d.kDiffusion,d.epsilonDiffusion})}}));
    }
    auto anchor=RANSWallCoordinate{(lower[0]+upper[0])/2,0,m_parameters.geometry.height};
    const double sample=SampleK(state,anchor),delta=sample*1e-5;
    auto plus=state,minus=state;
    for(int j=0;j<size[K()];++j){plus[offset[K()]+j]+=delta*m_constants[K()][j];minus[offset[K()]+j]-=delta*m_constants[K()][j];}
    const auto pp=StateProfile(plus,anchor,time,data),pm=StateProfile(minus,anchor,time,data);
    const double traceDirection=(pp->Dissipation(anchor[2])-pm->Dissipation(anchor[2]))/(2*delta);
    const double traceExpected=1.5*StateProfile(state,anchor,time,data)->Dissipation(anchor[2])/sample;
    Require(std::abs(traceDirection-traceExpected)<1e-8,"WJ endpoint-K dependent trace response failed");
    Require(hash==RANSWallR3::Hash(state),"WJ read-only check changed state");
    return J::Object({{"schema","1"},{"model",J::Quote("WJHighRe")},{"audit_revision",J::Quote(AuditRevision())},
        {"governing_solves","0"},{"state_unchanged","true"},{"state_sha256",J::Quote(hash)},
        {"frozen_blocks",J::List(blockReports)},{"pressure_modes",J::List(pressureReports)},
        {"max_same_state_defect",J::Number(blockError)},{"epsilon_retained_trace_column",J::Number(traceColumn)},
        {"zero_production_sink_point_error",J::Number(sinkPointError)},{"zero_production_sink_block_error",J::Number(sinkBlockError)},
        {"endpoint_K_cutoff_plus_minus",J::Numbers(Vector{pp->Reference().cutoff,pm->Reference().cutoff})},
        {"endpoint_K_trace_action_expected",J::Numbers(Vector{traceDirection,traceExpected})},
        {"terms",J::Object(terms)},{"independent_quadrature",J::List(samples)},
        {"extracted_constitutive",J::List(constitutive)}});
}
}
