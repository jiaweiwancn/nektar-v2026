///////////////////////////////////////////////////////////////////////////////
// File: RANSConstitutiveProbe.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <SolverUtils/Turbulence/RANSModels.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Nektar::SolverUtils;

namespace
{
// This executable only samples the production law. The independent indexed
// stress and its complex-step derivatives live in validate_cls_epawf_m1.py.
struct Case
{
    std::string name;
    RANSTensor gradient, direction;
    double k, epsilon, distance, nu;
    // 0: physical max; 1/2: formal strain/rotation limiting derivative.
    int branch;
    std::vector<double> steps;
};

void Tensor(std::ostream &out, const RANSTensor &value)
{
    out << '[';
    for (int i = 0; i < 3; ++i)
    {
        if (i) out << ',';
        out << '[';
        for (int j = 0; j < 3; ++j)
        {
            if (j) out << ',';
            if (!std::isfinite(value[i][j]))
                throw std::runtime_error("Nonfinite constitutive sample");
            out << value[i][j];
        }
        out << ']';
    }
    out << ']';
}

void State(std::ostream &out, const Case &c, double step)
{
    auto gradient = c.gradient;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            gradient[i][j] += step * c.direction[i][j];
    const auto v = c.branch == 0 ? EvaluateRANSClosure(RANSModel::CLSHighRe,
        c.k, c.epsilon, gradient, c.distance, {}, c.nu) :
        EvaluateRANSClosureLinearization(RANSModel::CLSHighRe, c.k, c.epsilon,
            gradient, c.distance, {}, c.nu, c.branch == 1);
    RANSTensor flux{};
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            flux[i][j] = c.nu * (gradient[i][j] + gradient[j][i]) -
                         v.reynoldsStress[i][j];
    out << "{\"gradient\":"; Tensor(out, gradient);
    out << ",\"stress\":"; Tensor(out, v.reynoldsStress);
    out << ",\"quadratic\":"; Tensor(out, v.quadraticStress);
    out << ",\"cubic\":"; Tensor(out, v.cubicStress);
    out << ",\"flux\":"; Tensor(out, flux);
    out << ",\"cMu\":" << v.cMu << ",\"nut\":" << v.eddyViscosity
        << ",\"strain\":" << v.strainInvariant
        << ",\"rotation\":" << v.rotationInvariant
        << ",\"production\":" << v.production
        << ",\"k_source\":" << v.kSource
        << ",\"epsilon_source\":" << v.epsilonSource
        << ",\"yap\":" << v.yap << '}';
}

void Write(std::ostream &out, const Case &c)
{
    out << "{\"name\":\"" << c.name << "\",\"branch\":\""
        << (c.branch == 0 ? "physical" : c.branch == 1 ? "strain" : "rotation")
        << "\",\"k\":" << c.k << ",\"epsilon\":" << c.epsilon
        << ",\"distance\":" << c.distance << ",\"nu\":" << c.nu
        << ",\"direction\":"; Tensor(out, c.direction);
    out << ",\"base\":"; State(out, c, 0);
    out << ",\"samples\":[";
    bool first = true;
    for (double step : c.steps)
    {
        if (!first) out << ',';
        first = false;
        out << "{\"step\":" << step << ",\"plus\":"; State(out, c, step);
        out << ",\"minus\":"; State(out, c, -step);
        out << '}';
    }
    out << "]}";
}
} // namespace

int main(int argc, char **argv)
{
    try
    {
        if (argc != 3 || std::string(argv[1]) != "--output")
            throw std::invalid_argument("Use --output PATH");
        const auto reference = RANSNeutralReference::CLS();
        const double q = reference.Q();
        RANSTensor base{}, transverseX{}, transverseZ{}, inplane{};
        base[0][2] = 1;
        transverseX[1][0] = 1;
        transverseZ[1][2] = 1;
        const RANSVector xi{.45, 0, std::sqrt(1 - .45 * .45)};
        const RANSVector a{-xi[2], 0, .45};
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j) inplane[i][j] = a[i] * xi[j];
        const std::vector<double> steps{1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
        std::vector<Case> cases;
        for (int branch : {0, 1, 2})
        {
            cases.push_back({"transverse_x", base, transverseX,
                             q, 1, 1, 0, branch, steps});
            cases.push_back({"transverse_z", base, transverseZ,
                             q, 1, 1, 0, branch, steps});
            if (branch != 0)
                cases.push_back({"reference_inplane", base, inplane,
                                 q, 1, 1, 0, branch, steps});
        }
        auto nearby = base;
        nearby[0][0] = 1e-4; nearby[2][2] = -1e-4;
        // The invariant gap is O(1e-8). These physical differences stay on
        // its strain branch; larger central differences would cross the tie.
        cases.push_back({"nearby_smooth_inplane", nearby, inplane,
                         q, 1, 1, 0, 0, {3e-9, 1e-9, 3e-10}});
        for (double offset : {1e-3, 1e-5})
            for (int sign : {-1, 1})
            {
                auto limit = base; limit[2][0] = sign * offset;
                cases.push_back({sign > 0 ? "strain_limit" : "rotation_limit",
                    limit, inplane, q, 1, 1, 0, 0,
                    {offset * .01, offset * .003, offset * .001}});
            }
        const RANSTensor general{{{.1, .4, .8}, {-.2, -.3, .5}, {.3, -.6, .2}}};
        cases.push_back({"general_3d", general, inplane, .7, .6, .3, 1.5e-5, 0, steps});
        auto twoDimensional = base;
        twoDimensional[0][0] = .2; twoDimensional[2][2] = -.2;
        twoDimensional[2][0] = .1;
        cases.push_back({"two_dimensional_three_stresses", twoDimensional,
                         inplane, .61, .73, .4, 0, 0, steps});
        cases.push_back({"transverse_molecular_control", base, transverseX,
                         q, 1, 1, .025, 0, steps});
        cases.push_back({"same_state_double_wall_distance", base, transverseX,
                         q, 1, 2, 0, 0, steps});
        std::ofstream out(argv[2]);
        if (!out) throw std::runtime_error("Cannot open constitutive report");
        out << std::setprecision(17) << "{\"schema\":1,\"q\":" << q
            << ",\"kappa\":" << reference.Kappa() << ",\"cMu\":" << reference.Cmu()
            << ",\"reference_profiles\":[";
        bool firstProfile = true;
        for (double z : {.1, 1., 2., 10., 100.})
        {
            const auto p = reference.AtHeight(z, .4, 2., 5.7);
            RANSTensor gradient{}; gradient[0][2] = p.velocityGradient;
            const Case profile{"profile", gradient, {}, p.k, p.epsilon, z, 0, 0, {}};
            if (!firstProfile) out << ',';
            firstProfile = false;
            out << "{\"z\":" << z << ",\"k\":" << p.k
                << ",\"epsilon\":" << p.epsilon << ",\"velocity\":" << p.velocity
                << ",\"epsilon_gradient\":" << p.epsilonGradient
                << ",\"epsilon_gradient_flux\":" << p.epsilonGradientFlux
                << ",\"state\":";
            State(out, profile, 0); out << '}';
        }
        out << "],\"cases\":[";
        for (std::size_t i = 0; i < cases.size(); ++i)
        {
            if (i) out << ',';
            Write(out, cases[i]);
        }
        out << "]}\n";
        out.close();
        if (!out) throw std::runtime_error("Cannot finish constitutive report");
        std::cout << "Sampled original CLS stress in " << cases.size()
                  << " constitutive cases; no stability verdict inferred.\n";
    }
    catch (const std::exception &error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
    return 0;
}
