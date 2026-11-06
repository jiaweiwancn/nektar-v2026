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
class RejectAugmentedFactor
{
public:
    explicit RejectAugmentedFactor(const R::Matrix &matrix) : a(matrix)
    {
        if (a.size() > 192) throw std::runtime_error("deliberate saddle factor failure");
    }
    R::Vector Apply(R::Vector rhs) const { return Factor(a).Apply(std::move(rhs)); }
private:
    R::Matrix a;
};
class RejectWeakObjectiveFactor
{
public:
    explicit RejectWeakObjectiveFactor(const R::Matrix &matrix) : a(matrix)
    {
        if (a[0][0] < .05) throw std::runtime_error("deliberate weak model rejection");
    }
    R::Vector Apply(R::Vector rhs) const { return Factor(a).Apply(std::move(rhs)); }
private:
    R::Matrix a;
};
class RejectEveryFactor
{
public:
    explicit RejectEveryFactor(const R::Matrix &)
    { throw std::runtime_error("deliberate numerical factor failure"); }
    R::Vector Apply(R::Vector) const { throw std::logic_error("unreachable factor"); }
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
    if (argc == 2 || (argc == 5 && std::string(argv[2]) == "--orthogonal-radius"))
    {
        std::ifstream in(argv[1], std::ios::binary);
        std::uint64_t shape[2]{};
        in.read(reinterpret_cast<char *>(shape), sizeof(shape));
        Check(in.good() && shape[0] > 0 && shape[0] <= 513 && shape[1] <= 1000000,
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
        R::QuadraticResult result;
        double replayRadius=0., replayMultiplier=0.;
        int radiusSteps=0,radiusSolves=0,warmFaces=0,warmFallbacks=0;
        if (argc == 5)
        {
            // Diagnostic only: the caller must have verified ordinary model
            // coordinates with orthonormal field columns. Remove the saved
            // radius multiplier, retaining the original merit/regularization
            // and every cut, then solve with the requested field radius.
            replayRadius=std::stod(argv[3]);
            const double oldMultiplier=std::stod(argv[4]);
            Check(n<=512 && replayRadius>0. && std::isfinite(replayRadius) &&
                  oldMultiplier>=0. && std::isfinite(oldMultiplier) &&
                  R::Length(initial)==0., "orthogonal-radius replay assumptions differ");
            R::Matrix columns(n,R::Vector(n));
            for (int i=0; i<n; ++i) { columns[i][i]=1.; h[i][i]-=oldMultiplier; }
            const auto bounded=R::SolveRadiusQuadratic<Factor>(
                h,g,cuts,initial,columns,replayRadius,500);
            replayMultiplier=bounded.multiplier;
            radiusSteps=bounded.steps;radiusSolves=bounded.solves;
            warmFaces=bounded.warmFaceSolves;warmFallbacks=bounded.warmFaceFallbacks;
            result=bounded.quadratic;
        }
        else result=R::SolveQuadratic<Factor>(h, g, cuts, initial, 500);
        std::cout << std::setprecision(17) << "{\"solved\":" << (result.solved?"true":"false")
                  << ",\"failure\":\"" << result.failure << "\",\"steps\":" << result.steps
                  << ",\"primal\":" << result.primal << ",\"stationarity\":" << result.stationarity
                  << ",\"complementarity\":" << result.complementarity << ",\"redundant_rows\":"
                  << result.redundantRows << ",\"diagnostic_radius\":" << replayRadius
                  << ",\"augmented_face_solves\":" << result.augmentedFaceSolves
                  << ",\"augmented_face_fallbacks\":" << result.augmentedFaceFallbacks
                  << ",\"radius_multiplier\":" << replayMultiplier
                  << ",\"radius_total_steps\":" << radiusSteps
                  << ",\"radius_model_solves\":" << radiusSolves
                  << ",\"warm_face_solves\":" << warmFaces
                  << ",\"warm_face_fallbacks\":" << warmFallbacks
                  << ",\"scope\":\"affine model replay only; not nonlinear stationarity or M8 convergence\",\"solution\":[";
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
    // Nonlinear destruction can make the full linear-model step ascend.
    auto ray=R::SelectSourceRay(1.,[](double f,double &value) {
        value=1.-f+4.*f*f; return true;
    });
    Close(ray.fraction,.125); Close(ray.merit,.9375);
    Check(ray.evaluations<=12,"unbounded surrogate sample budget");
    ray=R::SelectSourceRay(1.,[](double f,double &value) {
        if (f>.5) return false;
        value=std::abs(1.-2.*f); return true;
    });
    Close(ray.fraction,.5); Close(ray.merit,0.);
    Check(ray.invalid==1,"invalid surrogate state was not distinguished");
    ray=R::SelectSourceRay(1.,[](double f,double &value) { value=1.+f; return true; });
    Check(ray.fraction==0.,"nondecreasing surrogate replaced the control");
    auto derivative=R::DifferentiateModel(R::Vector{.125},R::Vector{.9375},0,1e-4,
        [](const R::Vector &z,R::Vector &value) {
            if (z[0]>.125) return false; // Only the inward one-sided branch is defined.
            value={1.-z[0]+4.*z[0]*z[0]}; return true;
        });
    Check(derivative.available && derivative.step<0.,"model derivative left its domain");
    Check(std::abs(derivative.coarse[0]+.0004)<1e-10 &&
          std::abs(derivative.fine[0]+.0002)<1e-10,"two-scale model derivative is incorrect");
    derivative=R::DifferentiateModel(R::Vector{0.},R::Vector{1.},0,1e-4,
        [](const R::Vector &,R::Vector &) { return false; });
    Check(!derivative.available,"undefined model derivative silently accepted");
    // Two active ground-like cuts forbid both signs of the first coordinate.
    // A fractional ray anchor restores a genuine model neighborhood without
    // modifying a physical state or evaluating outside the declared domain.
    int guarded=0;
    auto cornerModel=[&](const R::Vector &z,R::Vector &value) {
        const double a=1.-z[0]-z[1],b=1.+z[0]-3.*z[1];
        if (a<0. || b<0.) { ++guarded; return false; }
        value={a*a+2.*b*b+z[0]*z[0]+3.*z[1]*z[1]};
        return true;
    };
    auto derivativeError=[](const R::ModelDerivative &d) {
        return std::abs(d.coarse[0]-d.fine[0])/std::max(1.,std::abs(d.fine[0]));
    };
    R::Vector cornerValue;
    Check(cornerModel({.5,.5},cornerValue),"corner model was not defined");
    auto linearization=R::LinearizeModelOnRay(R::Vector{.5,.5},cornerValue,
        1e-4,0,.03,cornerModel,derivativeError);
    Check(!linearization.available && linearization.availableDerivatives==0 && guarded>0,
          "blocked corner derivative was silently used");
    linearization=R::LinearizeModelOnRay(R::Vector{.5,.5},cornerValue,
        1e-4,2,.03,cornerModel,derivativeError);
    Check(linearization.available && linearization.availableDerivatives==2 &&
          linearization.anchorReductions==1 && linearization.attempts==2,
          "interior source anchor did not recover the model derivatives");
    Close(linearization.at[0],.25); Close(linearization.at[1],.25);
    Check(std::abs(linearization.derivatives[0][0]-1.5)<.001 &&
          std::abs(linearization.derivatives[1][0]+5.5)<.002,
          "anchored derivatives differ from their independent quadratic values");
    auto interior=linearization;
    linearization=R::LinearizeModelOnRay(interior.at,interior.value,
        1e-4,2,.03,cornerModel,derivativeError);
    Check(linearization.available && linearization.anchorReductions==0 &&
          linearization.at==interior.at,"valid incumbent model was moved unnecessarily");
    linearization=R::LinearizeModelOnRay(R::Vector{.5,.5},cornerValue,
        1e-4,2,.03,[](const R::Vector &,R::Vector &) { return false; },derivativeError);
    Check(!linearization.available && linearization.attempts==3,
          "exhausted anchor attempts were silently accepted");
    // Coarse stencils can miss the declared accuracy at a smooth interior
    // point. The quadratic's independent derivative at zero is2; one reused
    // half stencil reduces its truncation error without moving the anchor.
    int precisionCalls=0;
    auto curvedModel=[&](const R::Vector &z,R::Vector &value) {
        ++precisionCalls;value={1.+2.*z[0]+1e4*z[0]*z[0]};return true;
    };
    linearization=R::LinearizeModelOnRay(R::Vector{0.},R::Vector{1.},
        2e-5,0,.03,curvedModel,derivativeError);
    Check(linearization.available && linearization.anchorReductions==0 &&
          linearization.refinedDerivatives==1 && linearization.derivativeRefinements==1 &&
          linearization.initialError>.03 && linearization.error<.03 && precisionCalls==3 &&
          linearization.minimumStep==1e-5 && linearization.stencils.size()==1,
          "an inaccurate coarse source stencil discarded a resolvable derivative");
    Check(std::abs(linearization.derivatives[0][0]-2.)<.06,
          "refined source action disagrees with the independent quadratic derivative");
    precisionCalls=0;
    linearization=R::LinearizeModelOnRay(R::Vector{0.},R::Vector{1.},
        2e-5,0,.03,[&](const R::Vector &z,R::Vector &value) {
            ++precisionCalls;value={z[0]==0.?1.:2.};return true;
        },derivativeError);
    Check(!linearization.available && linearization.derivativeRefinements==3 &&
          linearization.error>.03 && precisionCalls==5,
          "an unresolved source jump was silently accepted as a derivative");
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
    // A residual-reducing proposal must honor the same qualification funnel
    // that the nonlinear filter will enforce, even while outside the box.
    r = Solve({{16.,0.},{0.,1.}}, {-16.,0.},
              {{{1.,-1.},-2.}, {{0.,1.},2.02}, {{0.,-1.},0.}}, {0.,2.});
    Close(r.solution[0],.02); Close(r.solution[1],2.02);
    // A 50% restoration is useful even close to the original boundary and
    // when residual merit must rise; an absolute violation floor erases it.
    Check(R::QualificationDecrease(1e-5,5e-6,1.,1e-3), "small-violation restoration erased");
    Check(R::QualificationDecrease(1e-14,0.,1.,1e-3), "arithmetic-scale restoration erased");
    Check(!R::QualificationDecrease(1e-5,1e-5-1e-12,1.,1e-3), "insignificant restoration counted");
    Check(!R::QualificationDecrease(0.,0.,1.,1e-3), "unchanged qualification counted as progress");
    // A finite mixed penalty can stop outside a feasible qualification set.
    // At this exact soft-penalty fixed point, pure restoration still moves.
    const double w=1e-4, mu=1e-3, xstar=1./(1.+w), violation=1.-xstar;
    r = Solve({{w+mu,0.},{0.,1.}}, {w*xstar,0.},
              {{{-1.,-1.},-violation},{{0.,-1.},0.}}, {0.,violation});
    Close(r.solution[0],0.); Close(r.solution[1],violation);
    r = Solve({{mu,0.},{0.,1.}}, {0.,0.},
              {{{-1.,-1.},-violation},{{0.,-1.},0.}}, {0.,violation});
    Close(r.solution[0],violation/(1.+mu));
    Check(R::QualificationDecrease(violation,r.solution[1],1.,1e-3),
          "pure qualification restoration stalled at a mixed-penalty fixed point");
    // A tiny positive excess can otherwise remain outside after a soft
    // restoration. An inward model target crosses the original bound; the
    // physical qualification observation remains x <= 0, without an allowance.
    const double tinyExcess=1.8e-9, inwardTarget=1e-6;
    r = Solve({{mu,0.},{0.,1.}}, {0.,0.},
              {{{1.,-1.},-tinyExcess},{{0.,-1.},0.}}, {0.,tinyExcess});
    Check(tinyExcess+r.solution[0] > 0., "soft-restoration control unexpectedly interior");
    r = Solve({{mu,0.},{0.,1.}}, {0.,0.},
              {{{1.,-1.},-tinyExcess-inwardTarget},{{0.,-1.},0.}},
              {0.,tinyExcess+inwardTarget});
    Check(tinyExcess+r.solution[0] < -0.9*inwardTarget,
          "inward restoration model did not cross the unchanged original bound");
    // H(y,k)=y-3k and F(k,y)=k+y-4 have a reduced slope of 4.
    // Feedback in B, rather than post-direction settling, selects dk=1.
    r = Solve({{16.}}, {-16.}, {}, {0.});
    Close(r.solution[0],1.);
    // H(y,k)=y-3*k2, F=(k1+y-4,k2-1). The conditional
    // preconditioner is I. Its residual seed alone misses the reduced root;
    // -P^-1*(Jred-P)*v supplies an independent response direction and the
    // two-dimensional model recovers the known positive solution (1,1).
    const double root17=std::sqrt(17.);
    r = Solve({{50./17.}}, {-29./root17}, {}, {0.});
    const double restrictedK1=4.*r.solution[0]/root17,
                 restrictedK2=r.solution[0]/root17;
    Check(std::hypot(restrictedK1+3.*restrictedK2-4.,restrictedK2-1.) > .1,
          "single conditional seed unexpectedly solved the reduced problem");
    r = Solve({{50./17.,-7./root17},{-7./root17,1.}},
              {-29./root17,4.}, {}, {0.,0.});
    Close(4.*r.solution[0]/root17-r.solution[1],1.);
    Close(r.solution[0]/root17,1.);
    r = Solve({{16.}}, {-16e-15}, {}, {0.});
    Check(std::abs(r.solution[0]/1e-15 - 1.) < 1e-12,
          "an absolute inner step threshold erased a small Newton correction");
    R::Vector z{3.,4.,2.};
    const auto step = R::BoundCorrection(z, {{1.,0.},{1.,1.}}, .2);
    Close(R::Length(step),.2);
    Close(step[0],z[0]+z[1]); Close(step[1],z[1]); Close(z[2],2.);
    // Post-hoc radial scaling wastes a small radius on the weak direction of
    // an anisotropic model. The radius-constrained solve uses it on the useful
    // field direction, and leaves an independent restoration slack untouched.
    const R::Matrix ballH{{1.,0.,0.},{0.,100.,0.},{0.,0.,1.}},
                    ballV{{1.,0.},{0.,1.}};
    const R::Vector ballG{-10.,-100.,-7.};
    const auto ball = R::SolveRadiusQuadratic<Factor>(ballH, ballG,
        std::vector<Cut>{}, R::Vector(3), ballV, 1.);
    Check(ball.quadratic.solved && ball.multiplier > 0. &&
          ball.norm <= 1. && ball.norm >= .999 && ball.solves > 1,
          "radius-constrained model did not resolve its field bound");
    const auto &bz=ball.quadratic.solution;
    Close(bz[2],7.);
    const double ballObjective=.5*std::pow(bz[0]-10.,2)+50.*std::pow(bz[1]-1.,2);
    R::Vector radial{10.,1.,7.};
    R::BoundCorrection(radial,ballV,1.);
    const double radialObjective=.5*std::pow(radial[0]-10.,2)+50.*std::pow(radial[1]-1.,2);
    Check(ballObjective < 46.2 && radialObjective > 80.,
          "radius-aware direction lost its known advantage over radial scaling");
    const auto cutBall = R::SolveRadiusQuadratic<Factor>(ballH, ballG,
        std::vector<Cut>{{{1.,0.,0.},.3}}, R::Vector(3), ballV, 1.);
    Check(cutBall.quadratic.solved && cutBall.norm <= 1., "radius solve lost an active cut");
    Close(cutBall.quadratic.solution[0],.3);
    Check(std::abs(cutBall.quadratic.solution[1]-std::sqrt(.91)) < .002,
          "constrained circle model missed its known boundary solution");
    const auto interiorBall = R::SolveRadiusQuadratic<Factor>(
        R::Matrix{{1.}}, R::Vector{-.1}, std::vector<Cut>{}, R::Vector{0.},
        R::Matrix{{1.}}, 1.);
    Check(interiorBall.quadratic.solved && interiorBall.multiplier == 0. && interiorBall.solves == 1,
          "unnecessary radius penalty changed an interior model");
    Close(interiorBall.quadratic.solution[0],.1);
    // A radius larger than one remains a bound on the field displacement.
    // The identity model's exact constrained minimizer is (2.5,0).
    const auto largeBall = R::SolveRadiusQuadratic<Factor>(
        R::Matrix{{1.,0.},{0.,1.}}, R::Vector{-4.,0.}, std::vector<Cut>{},
        R::Vector(2), R::Matrix{{1.,0.},{0.,1.}}, 2.5);
    Check(largeBall.quadratic.solved && largeBall.multiplier > 0. &&
          largeBall.norm <= 2.5 && largeBall.norm >= 2.4975 &&
          std::abs(largeBall.quadratic.solution[0]-2.5) <= .0025,
          "large field radius lost its known constrained minimizer");
    Close(largeBall.quadratic.solution[1],0.);
    // A failed zero-multiplier objective is not a failed radius problem.
    // This weak scalar model has the independent constrained answer x=1.
    // The first two positive queries deliberately fail as well.
    const auto recoveredBall = R::SolveRadiusQuadratic<RejectWeakObjectiveFactor>(
        R::Matrix{{1e-12}}, R::Vector{-1.}, std::vector<Cut>{{{1.},2.}},
        R::Vector{0.}, R::Matrix{{1.}}, 1.);
    Check(recoveredBall.quadratic.solved && !recoveredBall.unconstrainedAvailable &&
          !recoveredBall.unconstrainedFailure.empty() &&
          recoveredBall.failedPositiveQueries==2 && recoveredBall.multiplier>0. &&
          recoveredBall.norm<=1. && recoveredBall.norm>=.999 &&
          recoveredBall.quadratic.primal<1e-9 && recoveredBall.quadratic.stationarity<1e-9,
          "failed unbounded query prevented a valid bounded model solve");
    Close(recoveredBall.quadratic.solution[0],1.);
    const auto exhaustedBall = R::SolveRadiusQuadratic<RejectEveryFactor>(
        R::Matrix{{1.}}, R::Vector{-1.}, std::vector<Cut>{},
        R::Vector{0.}, R::Matrix{{1.}}, 1.);
    Check(!exhaustedBall.quadratic.solved && !exhaustedBall.unconstrainedAvailable &&
          exhaustedBall.failedPositiveQueries==14 &&
          exhaustedBall.quadratic.failure=="radius_multiplier_bracket_exhausted",
          "exhausted positive queries silently accepted a partial model");
    const auto obliqueBall = R::SolveRadiusQuadratic<Factor>(
        R::Matrix{{1.,0.},{0.,1.}}, R::Vector{-3.,-4.}, std::vector<Cut>{},
        R::Vector(2), R::Matrix{{1.,0.},{1.,1.}}, .2);
    const auto &oz=obliqueBall.quadratic.solution;
    const double actualFieldNorm=std::hypot(oz[0]+oz[1],oz[1]);
    Check(obliqueBall.quadratic.solved && actualFieldNorm <= .2 && actualFieldNorm >= .1998,
          "oblique response coefficients were mistaken for the physical field norm");
    // H(y,k)=y-k1^2 has Y(k)=k1^2. At k=(0,0), the checked
    // qualification derivative of Y(k)-k2<=0 is (0,-1). The tangent
    // residual step (.1,0) leaves the true surface by .01. Correcting the
    // measured trial observation in the same response space restores it,
    // keeps ||delta||<.05, and still reduces F=(k1-.1,k2) from its base.
    const double curvedViolation=.01, margin=1e-6;
    const auto curved = R::SolveRadiusQuadratic<Factor>(
        R::Matrix{{1.,0.,0.},{0.,1.,0.},{0.,0.,1e8}}, R::Vector(3),
        std::vector<Cut>{{{0.,-1.,-1.},-curvedViolation-margin},{{0.,0.,-1.},0.}},
        R::Vector{0.,0.,curvedViolation+margin+1e-10},
        R::Matrix{{1.,0.},{0.,1.}}, .05);
    Check(curved.quadratic.solved && curved.norm <= .05,
          "nonlinear-qualification correction model failed");
    const auto &cz=curved.quadratic.solution;
    const double actualConstraint=std::pow(.1+cz[0],2)-cz[1],
                 actualResidual=std::hypot(cz[0],cz[1]);
    Check(actualConstraint < 0. && actualResidual < .02 &&
          std::hypot(.1+cz[0],cz[1]) < .11,
          "curved qualification was not restored with useful true residual progress");
    // A lower bound on the nonlinear H response Y(k)=k^2: the derivative
    // from base k=1.5 is too steep at trial k=.95. Three unchanged chord
    // corrections still violate Y>=1. A secant measured from two true H
    // responses restores the bound in two corrections, while reducing the
    // original F=k-.9 relative to its base value .6.
    double staleK=.95;
    for (int i=0; i<3; ++i) staleK += (1.-staleK*staleK+1e-6)/3.;
    Check(staleK*staleK<1., "stale chord unexpectedly restored the curved bound");
    const double secantStart=.95, firstChange=(1.-secantStart*secantStart+1e-6)/3.;
    const double afterFirst=secantStart+firstChange;
    R::Matrix response{{3.}};
    Check(R::UpdateSecantActions(response,{{1.}},{firstChange},
          {afterFirst*afterFirst-secantStart*secantStart},1e-6),
          "measured response secant was not used");
    Close(response[0][0],afterFirst+secantStart);
    const double restored=afterFirst+(1.-afterFirst*afterFirst+1e-6)/response[0][0];
    Check(restored*restored>=1. && std::abs(restored-.9)<.11 && std::abs(restored-1.5)<.55,
          "response secant lost actual qualification, residual progress or displacement bound");
    R::Matrix obliqueResponse{{2.,3.},{4.,5.}};
    Check(R::UpdateSecantActions(obliqueResponse,{{1.,0.},{1.,1.}}, {1.,-1.}, {7.,11.},1e-6),
          "oblique response secant failed");
    Close(obliqueResponse[0][0],2.); Close(obliqueResponse[0][1],3.);
    Close(obliqueResponse[0][0]-obliqueResponse[1][0],7.);
    Close(obliqueResponse[0][1]-obliqueResponse[1][1],11.);
    const auto beforeSmall=response;
    Check(!R::UpdateSecantActions(response,{{1.}},{1e-9},{1e-8},1e-6) && response==beforeSmall,
          "a below-signal secant changed the response model");
    // Identical active observations must supply the same next field witnesses
    // even when an equivalent optimization face reorders or repeats its rows.
    const std::vector<R::WitnessKey> witnesses{{5,0,3},{3,0,3},{2,2,3},{4,0,4},{2,2,3}};
    std::vector<int> witnessOrder{0,1,2,3,4};
    do
    {
        std::vector<R::WitnessKey> permuted;
        for (int i:witnessOrder) permuted.push_back(witnesses[i]);
        Check(R::SelectWitnesses(permuted)==std::vector<int>({2,3,4}),
              "equivalent active-row order changed next response witnesses");
    } while (std::next_permutation(witnessOrder.begin(),witnessOrder.end()));
    bool inconsistentWitness=false;
    try { R::SelectWitnesses({{2,0,1},{2,1,1}}); }
    catch (const std::logic_error &) { inconsistentWitness=true; }
    Check(inconsistentWitness,"conflicting observation metadata was treated as a valid witness");
    // Exercise both sides of the medium-model dispatch with closed answers.
    // The64-coordinate control keeps the original null-space path;96 uses
    // checked primal faces, warm radius queries and certified cut omission.
    for (int size : {64,96})
    {
        R::Matrix h(size,R::Vector(size)), v=h;
        R::Vector g(size), upper(size), lower(size);
        for (int i=0;i<size;++i) { h[i][i]=i+1.; v[i][i]=1.; }
        g[0]=-2.; g[1]=3.; upper[0]=1.; lower[1]=-1.;
        const std::vector<Cut> bounds{{upper,.25},{lower,.5}};
        const auto q=R::SolveQuadratic<Factor>(h,g,bounds,R::Vector(size));
        Check(q.solved && q.primal<1e-10 && q.stationarity<1e-8 &&
              q.complementarity<1e-8 &&
              (size==64 ? q.augmentedFaceSolves==0 : q.augmentedFaceSolves>0),
              "medium-model primal-face dispatch missed its checked solution");
        Close(q.solution[0],.25); Close(q.solution[1],-.5);
        g[0]=1.;
        const auto changed=R::SolveQuadratic<Factor>(h,g,bounds,R::Vector(size));
        const auto warm=R::SolveQuadratic<Factor>(h,g,bounds,q.solution,120,q.active);
        Check(changed.solved && warm.solved && warm.stationarity<1e-8,
              "medium warm face did not release an obsolete bound");
        for (int i=0;i<size;++i) Close(changed.solution[i],warm.solution[i]);
        Close(warm.solution[0],-1.); Close(warm.solution[1],-.5);
        for (int i=0;i<size;++i) h[i][i]=1.;
        h[1][1]=100.; g.assign(size,0.);
        g[0]=-271.3; g[1]=-1000.*std::sqrt(.91);
        R::Vector far(size), near(size); far[1]=1.; near[0]=1.;
        std::vector<Cut> cuts(5000,Cut{far,5.}); cuts.push_back({near,.3});
        const auto bounded=R::SolveRadiusQuadratic<Factor>(
            h,g,cuts,R::Vector(size),v,1.);
        Check(bounded.quadratic.solved && bounded.norm<=1. && bounded.norm>=.999 &&
              bounded.quadratic.primal<1e-10 &&
              bounded.quadratic.active==std::vector<int>{5000} &&
              (size==64 ? bounded.ballRedundantCuts==0 && bounded.warmFaceSolves==0
                        : bounded.ballRedundantCuts==5000 && bounded.warmFaceSolves>0),
              "medium radius solve lost original cuts, active indices or its field bound");
        Close(bounded.quadratic.solution[0],.3);
        Check(std::abs(bounded.quadratic.solution[1]-std::sqrt(.91))<.002,
              "medium radius model missed its known boundary answer");
    }
    // A closed-form large model exercises the response-coordinate saddle
    // path and its numerical-failure fallback with two independent bounds.
    R::Matrix largeH(192,R::Vector(192));
    R::Vector largeG(192), upperCut(192), lowerCut(192);
    for (int i=0;i<192;++i) largeH[i][i]=i+1.;
    largeG[0]=-2.; largeG[1]=3.; upperCut[0]=1.; lowerCut[1]=-1.;
    auto dependent=upperCut; dependent[2]=1e-14;
    const std::vector<Cut> largeCuts{{upperCut,.25},{lowerCut,.5},{dependent,.25}};
    const auto large=R::SolveQuadratic<Factor>(largeH,largeG,largeCuts,R::Vector(192));
    Check(large.solved && large.augmentedFaceSolves>0 && large.augmentedFaceFallbacks==0 &&
          large.primal<1e-10 && large.stationarity<1e-8 && large.complementarity<1e-8,
          "large orthonormal-face saddle model failed");
    Close(large.solution[0],.25); Close(large.solution[1],-.5);
    for (int i=2;i<192;++i) Close(large.solution[i],0.);
    auto changedG=largeG;changedG[0]=1.;
    const auto coldChanged=R::SolveQuadratic<Factor>(largeH,changedG,largeCuts,R::Vector(192));
    const auto warmChanged=R::SolveQuadratic<Factor>(
        largeH,changedG,largeCuts,large.solution,120,large.active);
    Check(coldChanged.solved && warmChanged.solved && warmChanged.primal<1e-10 &&
          warmChanged.stationarity<1e-8 && warmChanged.complementarity<1e-8,
          "warm face failed to remove an obsolete constraint");
    for (int i=0;i<192;++i) Close(warmChanged.solution[i],coldChanged.solution[i]);
    Close(warmChanged.solution[0],-1.);Close(warmChanged.solution[1],-.5);
    const auto staleFace=R::SolveQuadratic<Factor>(
        largeH,changedG,largeCuts,R::Vector(192),120,large.active);
    Check(!staleFace.solved && staleFace.failure=="quadratic_inconsistent_warm_face",
          "inconsistent warm face was silently accepted");
    const auto fallback=R::SolveQuadratic<RejectAugmentedFactor>(
        largeH,largeG,largeCuts,R::Vector(192));
    Check(fallback.solved && fallback.augmentedFaceFallbacks>0 && fallback.augmentedFaceSolves==0,
          "large face factor failure did not use the original null-space solve");
    for (int i=0;i<192;++i) Close(fallback.solution[i],large.solution[i]);
    R::Matrix largeBallH(192,R::Vector(192)), largeBallV=largeBallH;
    R::Vector largeBallG(192), far(192), near(192);
    for (int i=0;i<192;++i) largeBallH[i][i]=largeBallV[i][i]=1.;
    largeBallH[1][1]=100.;
    largeBallG[0]=-271.3; largeBallG[1]=-1000.*std::sqrt(.91);
    far[1]=1.; near[0]=1.;
    std::vector<Cut> manyCuts(5000,Cut{far,5.});
    manyCuts.push_back({near,.3});
    const auto pruned=R::SolveRadiusQuadratic<Factor>(
        largeBallH,largeBallG,manyCuts,R::Vector(192),largeBallV,1.);
    Check(pruned.quadratic.solved && pruned.ballRedundantCuts==5000 &&
          pruned.warmFaceSolves>0 && pruned.warmFaceFallbacks==0 &&
          pruned.norm<=1. && pruned.norm>=.999 &&
          pruned.quadratic.active==std::vector<int>{5000},
          "ball-inactive cuts or original active-row indices differ");
    Close(pruned.quadratic.solution[0],.3);
    Check(std::abs(pruned.quadratic.solution[1]-std::sqrt(.91))<.002,
          "pruned ball model missed its closed-form boundary solution");
    // A nearly cancelling field basis allows large coefficients inside the
    // physical ball. Its far coefficient cut must remain in the model.
    largeBallV[1][0]=1.; largeBallV[1][1]=1e-3;
    largeBallH[1][1]=1.; largeBallG.assign(192,0.);
    largeBallG[0]=10.; largeBallG[1]=-10.;
    const auto unpruned=R::SolveRadiusQuadratic<Factor>(
        largeBallH,largeBallG,manyCuts,R::Vector(192),largeBallV,1.);
    Check(unpruned.quadratic.solved && unpruned.ballRedundantCuts==0 && unpruned.norm<=1.,
          "ill-conditioned field geometry lost a necessary coefficient cut");
    Close(unpruned.quadratic.solution[1],5.);
    Check(std::abs(unpruned.quadratic.solution[0]+6.)<.002,
          "oblique field radius was confused with a coefficient radius");
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
