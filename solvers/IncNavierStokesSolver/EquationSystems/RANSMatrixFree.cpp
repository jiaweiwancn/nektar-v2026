///////////////////////////////////////////////////////////////////////////////
// File: RANSMatrixFree.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSMatrixFree.h>
#include <algorithm>
#include <cmath>

namespace Nektar
{
RANSMatrixFree::RANSMatrixFree(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    int dimension, NekDouble tolerance, RobinOperator robin)
    : m_fields(fields), m_dim(dimension), m_tolerance(tolerance), m_robin(std::move(robin))
{
    ASSERTL0(std::isfinite(tolerance) && tolerance > 0, "Invalid RANS Krylov tolerance");
}

RANSMatrixFree::~RANSMatrixFree()
{
    // The field owns its native systems. Destruction can follow MPI_Finalize;
    // do not invoke collective cache cleanup here.
}

int RANSMatrixFree::Solve(int f, const Array<OneD, NekDouble> &rhs,
    const Array<OneD, NekDouble> &diffusion, const Array<OneD, NekDouble> &mass,
    NekDouble alpha, Array<OneD, NekDouble> &output)
{
    auto field = std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[f]);
    const auto map = field->GetLocalToGlobalMap();
    ASSERTL0(map->GetGlobalSysSolnType() != MultiRegions::eIterativeFull,
        "RANSMatrixFree uses a native static-condensation preconditioner; do not select IterativeFull");
    const int nc = field->GetNcoeffs(), nq = field->GetTotPoints();
    const int ng = map->GetNumGlobalCoeffs(), nd = map->GetNumGlobalDirBndCoeffs();
    Preconditioner *preconditioner = nullptr;
    for (const auto &p : m_preconditioners)
        if (p->field == f && p->alpha == alpha) preconditioner = p.get();
    if (!preconditioner)
    {
        int count = 0;
        for (const auto &p : m_preconditioners) count += p->field == f;
        ASSERTL0(count < 2, "RANSMatrixFree requires constant-step IMEX1/2 (at most two inertia coefficients)");
        auto p = std::make_unique<Preconditioner>(); p->field = f; p->alpha = alpha;
        p->diffusion = Array<OneD, NekDouble>(nq); p->mass = Array<OneD, NekDouble>(nq);
        Vmath::Vcopy(nq, diffusion, 1, p->diffusion, 1); Vmath::Vcopy(nq, mass, 1, p->mass, 1);
        StdRegions::VarCoeffMap coefficients;
        coefficients[StdRegions::eVarCoeffD00] = p->diffusion;
        coefficients[StdRegions::eVarCoeffD11] = p->diffusion;
        if (m_dim == 3) coefficients[StdRegions::eVarCoeffD22] = p->diffusion;
        coefficients[StdRegions::eVarCoeffMass] = p->mass;
        StdRegions::ConstFactorMap factors; factors[StdRegions::eFactorLambda] = 1;
        p->key = std::make_unique<MultiRegions::GlobalLinSysKey>(StdRegions::eHelmholtz, map, factors, coefficients);
        preconditioner = p.get(); m_preconditioners.push_back(std::move(p));
    }
    Array<OneD, NekDouble> local(nc), applied(nc), physical(nq), b(ng), x(ng), ax(ng),
        r(ng), z(ng), p(ng), ap(ng), counts(ng), countLocal(nc);
    Array<OneD, NekDouble> massPhysical(nq), massCoefficients(nc);
    Fields flux(m_dim);
    for (auto &v : flux) v = Array<OneD, NekDouble>(nq);
    for (int i = 0; i < nc; ++i) countLocal[i] = map->GetLocalToGlobalSign(i);
    map->Assemble(countLocal, counts);
    Vmath::Smul(nc, -1.0, rhs, 1, local, 1);
    field->ImposeNeumannConditions(local); field->ImposeRobinConditions(local);
    map->Assemble(local, b);
    field->ImposeDirichletConditions(field->UpdateCoeffs());
    map->LocalToGlobal(field->GetCoeffs(), x, true);
    bool hasRobin = false;
    std::vector<std::pair<MultiRegions::ExpListSharedPtr, Array<OneD, NekDouble>>> saved;
    for (int boundary = 0; boundary < field->GetBndCondExpansions().size(); ++boundary)
    {
        const auto bc = field->GetBndConditions()[boundary];
        hasRobin |= bc->GetBoundaryConditionType() == SpatialDomains::eRobin;
        if (bc->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
        {
            auto exp = field->GetBndCondExpansions()[boundary];
            Array<OneD, NekDouble> copy(exp->GetNcoeffs());
            Vmath::Vcopy(copy.size(), exp->GetCoeffs(), 1, copy, 1);
            saved.emplace_back(exp, copy);
            Vmath::Zero(copy.size(), exp->UpdateCoeffs(), 1);
        }
    }
    auto restore = [&]() {
        for (const auto &item : saved)
            Vmath::Vcopy(item.second.size(), item.second, 1, item.first->UpdateCoeffs(), 1);
    };
    auto dot = [&](const Array<OneD, NekDouble> &a, const Array<OneD, NekDouble> &v) {
        double value = 0;
        for (int i = nd; i < ng; ++i)
            if (map->GetGlobalToUniversalMapUnique(i) > 0) value += a[i]*v[i];
        field->GetComm()->AllReduce(value, LibUtilities::ReduceSum);
        return value;
    };
    auto apply = [&](const Array<OneD, NekDouble> &v, Array<OneD, NekDouble> &av) {
        map->GlobalToLocal(v, local);
        field->BwdTrans(local, physical);
        // Assemble the current weak operator directly. Some native Helmholtz
        // collection backends accept diffusion maps but omit variable mass;
        // using the elementary CG operations covers both terms on every shape.
        for (int d = 0; d < m_dim; ++d)
        {
            field->PhysDeriv(MultiRegions::DirCartesianMap[d], physical, flux[d]);
            Vmath::Vmul(nq, diffusion, 1, flux[d], 1, flux[d], 1);
        }
        field->IProductWRTDerivBase(flux, applied);
        Vmath::Vmul(nq, mass, 1, physical, 1, massPhysical, 1);
        field->IProductWRTBase(massPhysical, massCoefficients);
        Vmath::Vadd(nc, applied, 1, massCoefficients, 1, applied, 1);
        if (hasRobin)
        {
            ASSERTL0(m_robin, "Missing RANS Robin matrix-free operator");
            Fields value(m_dim+2); value[f] = physical;
            m_robin(f, value, applied);
        }
        map->Assemble(applied, av);
    };
    auto precondition = [&]() {
        // The native iterative system smooths the magnitude of successive
        // right-hand sides. Normalize each application so a small Krylov
        // residual is not discarded relative to a much larger earlier RHS.
        const double scale = std::sqrt(dot(r,r));
        ASSERTL0(std::isfinite(scale) && scale > 0, "Invalid RANS preconditioner residual scale");
        map->GlobalToLocal(r, local);
        for (int i = 0; i < nc; ++i) local[i] /= scale*counts[map->GetLocalToGlobalMap(i)];
        Vmath::Zero(nc, applied, 1);
        field->GlobalSolve(*preconditioner->key, local, applied);
        map->LocalToGlobal(applied, z, true);
        Vmath::Smul(ng,scale,z,1,z,1);
        Vmath::Zero(nd, z, 1);
    };
    int iteration = 0;
    try
    {
        apply(x, ax);
        for (int i = 0; i < ng; ++i) r[i] = i < nd ? 0 : b[i] - ax[i];
        const double target = m_tolerance * std::max(1.0, std::sqrt(dot(b, b)));
        double rz = 0;
        bool restart = true;
        for (; iteration < 1000; ++iteration)
        {
            if (std::sqrt(dot(r, r)) <= target)
            {
                apply(x, ax);
                for (int i = nd; i < ng; ++i) r[i] = b[i] - ax[i];
                if (std::sqrt(dot(r, r)) <= 5*target) break;
                restart = true;
            }
            precondition();
            const double next = dot(r, z);
            ASSERTL0(std::isfinite(next) && next > 0, "RANS matrix-free preconditioner is not positive definite");
            const double beta = restart ? 0.0 : next/rz;
            for (int i = 0; i < ng; ++i) p[i] = i < nd ? 0 : (restart ? z[i] : z[i] + beta*p[i]);
            rz = next; restart = false;
            apply(p, ap);
            const double curvature = dot(p, ap);
            ASSERTL0(std::isfinite(curvature) && curvature > 0, "RANS matrix-free Helmholtz operator is not positive definite");
            const double step = rz/curvature;
            for (int i = nd; i < ng; ++i) { x[i] += step*p[i]; r[i] -= step*ap[i]; }
            if (iteration % 20 == 19)
            {
                apply(x, ax);
                for (int i = nd; i < ng; ++i) r[i] = b[i] - ax[i];
                restart = true;
            }
        }
        ASSERTL0(iteration < 1000, "RANS matrix-free CG failed to converge");
        map->GlobalToLocal(x, field->UpdateCoeffs());
        field->BwdTrans(field->GetCoeffs(), output);
    }
    catch (...) { restore(); throw; }
    restore();
    return iteration;
}
}
