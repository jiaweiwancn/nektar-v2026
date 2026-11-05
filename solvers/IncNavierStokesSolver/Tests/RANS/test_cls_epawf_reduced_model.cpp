///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
// Isolated known-answer checks of the small reduced quadratic model. These
// checks do not assemble native fields and do not certify M8 convergence.
#include <IncNavierStokesSolver/EquationSystems/RANSWallReducedModel.h>
#include <iostream>
#include <fstream>
#include <cstdint>
#include <iomanip>

namespace R = Nektar::RANSWallReduced;
namespace
{
struct Cut { R::Vector a; double b; };
class Factor
{
public:
    explicit Factor(const R::Matrix &matrix) : a(matrix) {}
    R::Vector Apply(R::Vector rhs) const
    {
        auto m = a;
        const int n = rhs.size();
        for (int j = 0; j < n; ++j)
        {
            int pivot = j;
            for (int i = j + 1; i < n; ++i)
                if (std::abs(m[i][j]) > std::abs(m[pivot][j])) pivot = i;
            if (std::abs(m[pivot][j]) < 1e-15)
                throw std::runtime_error("singular test model");
            std::swap(m[pivot], m[j]);
            std::swap(rhs[pivot], rhs[j]);
            for (int i = j + 1; i < n; ++i)
            {
                const double f = m[i][j] / m[j][j];
                for (int k = j; k < n; ++k) m[i][k] -= f * m[j][k];
                rhs[i] -= f * rhs[j];
            }
        }
        for (int j = n - 1; j >= 0; --j)
        {
            for (int k = j + 1; k < n; ++k) rhs[j] -= m[j][k] * rhs[k];
            rhs[j] /= m[j][j];
        }
        return rhs;
    }
private:
    R::Matrix a;
};
void Check(bool valid, const char *message)
{
    if (!valid) throw std::runtime_error(message);
}
void Close(double actual, double expected)
{
    Check(std::abs(actual - expected) < 1e-8, "wrong known-answer correction");
}
R::QuadraticResult Solve(const R::Matrix &h, const R::Vector &g,
                         const std::vector<Cut> &cuts, const R::Vector &start)
{
    const auto result = R::SolveQuadratic<Factor>(h, g, cuts, start);
    const double dualTolerance = std::max(1e-8, 1e-12 * R::Length(g));
    Check(result.solved && result.primal < 1e-9 &&
          result.stationarity < dualTolerance && result.complementarity < dualTolerance,
          "known-answer model failed its affine KKT checks");
    return result;
}
}

int main(int argc, char **argv)
{
    if (argc == 2)
    {
        std::ifstream in(argv[1], std::ios::binary);
        std::uint64_t shape[2]{};
        in.read(reinterpret_cast<char *>(shape), sizeof(shape));
        Check(in.good() && shape[0] > 0 && shape[0] <= 64 && shape[1] <= 1000000,
              "invalid native affine replay packet shape");
        const int n = shape[0];
        auto read = [&](R::Vector &v) {
            in.read(reinterpret_cast<char *>(v.data()), v.size() * sizeof(double));
        };
        R::Matrix h(n, R::Vector(n));
        R::Vector g(n), initial(n);
        for (auto &row : h) read(row);
        read(g); read(initial);
        std::vector<Cut> cuts;
        for (std::uint64_t i = 0; i < shape[1]; ++i)
        {
            Cut cut{R::Vector(n),0.}; read(cut.a);
            in.read(reinterpret_cast<char *>(&cut.b), sizeof(double));
            cuts.push_back(std::move(cut));
        }
        Check(in.good() && in.peek() == std::char_traits<char>::eof(), "truncated or extended model packet");
        const auto result = R::SolveQuadratic<Factor>(h, g, cuts, initial, 500);
        std::cout << std::setprecision(17) << "{\"solved\":" << (result.solved?"true":"false")
                  << ",\"failure\":\"" << result.failure << "\",\"steps\":" << result.steps
                  << ",\"primal\":" << result.primal << ",\"stationarity\":" << result.stationarity
                  << ",\"complementarity\":" << result.complementarity << ",\"redundant_rows\":"
                  << result.redundantRows << ",\"solution\":[";
        for (int i = 0; i < n; ++i) std::cout << (i?",":"") << result.solution[i];
        std::cout << "]}\n";
        return result.solved ? 0 : 1;
    }
    auto r = Solve({{2.,0.},{0.,4.}}, {-2.,-8.}, {}, {0.,0.});
    Close(r.solution[0], 1.); Close(r.solution[1], 2.);
    r = Solve({{2.,0.},{0.,4.}}, {-2.,-8.},
              {{{0.,1.}, .5}, {{0.,1.}, .5}}, {0.,0.});
    Close(r.solution[0], 1.); Close(r.solution[1], .5);
    Check(r.active.size() == 1, "duplicate cut entered the active set");
    // A full scan must discover the sum constraint and then ground contact.
    r = Solve({{1.,0.},{0.,1.}}, {-2.,-4.},
              {{{-1.,0.},0.}, {{0.,1.},10.}, {{1.,1.},1.}}, {0.,0.});
    Close(r.solution[0],0.); Close(r.solution[1],1.);
    Check(r.active.size() == 2, "missed a newly active observation");
    r = Solve({{1e10,0.},{0.,1e6}}, {-1e10,-1e6*(1.+1e-12)},
              {{{1.,0.},.5}, {{1.,1e-7},.5000001},
               {{2.,1e-7},1.0000001}}, {0.,0.});
    Close(r.solution[0],.5); Close(r.solution[1],1.);
    // Minimax restoration balances two opposing gradient witnesses.
    r = Solve({{.001,0.},{0.,1.}}, {0.,0.},
              {{{-1.,-1.},-2.}, {{1.,-1.},-.5}, {{0.,-1.},0.},
               {{-1.,0.},0.}}, {0.,2.});
    Close(r.solution[0],.75); Close(r.solution[1],1.25);
    // An independent lower-bound observation changes that restoration step.
    r = Solve({{.001,0.},{0.,1.}}, {0.,0.},
              {{{-1.,-1.},-2.}, {{1.,-1.},-.5}, {{0.,-1.},0.},
               {{-1.,0.},-.9}}, {.9,2.});
    Close(r.solution[0],.9); Close(r.solution[1],1.4);
    // H(y,k)=y-3k and F(k,y)=k+y-4 have a reduced slope of 4.
    // Feedback in B, rather than post-direction settling, selects dk=1.
    r = Solve({{16.}}, {-16.}, {}, {0.});
    Close(r.solution[0],1.);
    r = Solve({{16.}}, {-16e-15}, {}, {0.});
    Check(std::abs(r.solution[0]/1e-15 - 1.) < 1e-12,
          "an absolute inner step threshold erased a small Newton correction");
    R::Vector z{3.,4.,2.};
    const auto step = R::BoundCorrection(z, {{1.,0.},{1.,1.}}, .2);
    Close(R::Length(step),.2);
    Close(step[0],z[0]+z[1]); Close(step[1],z[1]); Close(z[2],2.);
    const auto infeasible = R::SolveQuadratic<Factor>(
        R::Matrix{{1.}}, R::Vector{0.}, std::vector<Cut>{{{1.},0.}}, R::Vector{1.});
    Check(!infeasible.solved && infeasible.failure == "quadratic_infeasible_start",
          "infeasible model start was silently accepted");
    const auto exhausted = R::SolveQuadratic<Factor>(
        R::Matrix{{1.}}, R::Vector{-1.}, std::vector<Cut>{}, R::Vector{0.}, 1);
    Check(!exhausted.solved && exhausted.failure == "quadratic_iteration_limit",
          "exhausted model solve was silently accepted");
    const auto singular = R::SolveQuadratic<Factor>(
        R::Matrix{{0.}}, R::Vector{0.}, std::vector<Cut>{}, R::Vector{0.});
    Check(!singular.solved && singular.failure.find("quadratic_factor_failure") == 0,
          "factorization failure was silently accepted");
    bool badShape = false;
    try { R::BoundCorrection(z, {{1.},{1.,2.}}, .2); }
    catch (const std::logic_error &) { badShape = true; }
    Check(badShape, "programming error was treated as a numerical retry");
    std::cout << "passed_reduced_model_checks (not M8 convergence)\n";
}
