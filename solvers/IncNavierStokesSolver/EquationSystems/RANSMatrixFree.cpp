///////////////////////////////////////////////////////////////////////////////
// File: RANSMatrixFree.cpp
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#include <IncNavierStokesSolver/EquationSystems/RANSMatrixFree.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace Nektar
{
RANSMatrixFree::RANSMatrixFree(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    int dimension, NekDouble tolerance, RobinOperator robin, NekDouble gradDiv)
    : m_fields(fields), m_dim(dimension), m_tolerance(tolerance),
      m_gradDiv(gradDiv), m_robin(std::move(robin))
{
    ASSERTL0(std::isfinite(tolerance) && tolerance > 0, "Invalid RANS Krylov tolerance");
    ASSERTL0(std::isfinite(gradDiv) && gradDiv >= 0, "Invalid RANS grad-div coefficient");
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
        // Scale by the assembled effective RHS, including Dirichlet lifting.
        // A dimensionful max(1,||b||) creates an accuracy floor on small 3D
        // domains and makes finite-difference Newton probes inconsistent.
        Array<OneD,NekDouble> lifting(ng,0.0),effective(ng,0.0);
        for (int i = 0; i < nd; ++i) lifting[i] = x[i];
        apply(lifting,ax);
        for (int i = nd; i < ng; ++i) effective[i] = b[i]-ax[i];
        const double norm = std::sqrt(dot(effective,effective));
        const double roundoff = 32*std::numeric_limits<double>::epsilon()*
            (std::sqrt(dot(b,b))+std::sqrt(dot(ax,ax)));
        const double target = std::max(m_tolerance*norm,roundoff);
        if (norm == 0) Vmath::Vcopy(ng,lifting,1,x,1);
        apply(x, ax);
        for (int i = 0; i < ng; ++i) r[i] = i < nd ? 0 : b[i] - ax[i];
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

int RANSMatrixFree::SolveVelocity(const Fields &rhs, const Fields &diffusion,
    NekDouble alpha, Fields &output)
{
    // Terrain walls couple Cartesian components through Pt and n*n^T. Solve
    // that SPD block implicitly: lagging its off-diagonals is coordinate
    // dependent and can destabilize a rigidly rotated version of a flat case.
    auto first = std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[0]);
    const auto map = first->GetLocalToGlobalMap();
    const int nc = first->GetNcoeffs(), nq = first->GetTotPoints();
    const int ng = map->GetNumGlobalCoeffs(), nd = map->GetNumGlobalDirBndCoeffs();
    std::vector<std::shared_ptr<MultiRegions::ContField>> fields;
    std::vector<Preconditioner *> preconditioners;
    for (int f = 0; f < m_dim; ++f)
    {
        auto field = std::dynamic_pointer_cast<MultiRegions::ContField>(m_fields[f]);
        const auto other = field->GetLocalToGlobalMap();
        ASSERTL0(other->GetNumGlobalCoeffs() == ng && other->GetNumGlobalDirBndCoeffs() == nd,
                 "Terrain velocity components need identical scalar assembly maps");
        for (int i = 0; i < nc; ++i)
            ASSERTL0(map->GetLocalToGlobalMap(i) == other->GetLocalToGlobalMap(i) &&
                     map->GetLocalToGlobalSign(i) == other->GetLocalToGlobalSign(i),
                     "Terrain velocity assembly numbering differs between components");
        fields.push_back(field);
        Preconditioner *preconditioner = nullptr;
        for (const auto &p : m_preconditioners)
            if (p->field == f && p->alpha == alpha) preconditioner = p.get();
        if (!preconditioner)
        {
            int count = 0;
            for (const auto &p : m_preconditioners) count += p->field == f;
            ASSERTL0(count < 2,"Terrain vector solve requires constant-step IMEX1/2");
            auto p = std::make_unique<Preconditioner>(); p->field = f; p->alpha = alpha;
            p->diffusion = Array<OneD,NekDouble>(nq); p->mass = Array<OneD,NekDouble>(nq,alpha);
            Vmath::Vcopy(nq,diffusion[f],1,p->diffusion,1);
            p->directional = Array<OneD,NekDouble>(nq);
            Vmath::Sadd(nq,m_gradDiv,p->diffusion,1,p->directional,1);
            StdRegions::VarCoeffMap coefficients;
            coefficients[StdRegions::eVarCoeffD00] = f == 0 ? p->directional : p->diffusion;
            coefficients[StdRegions::eVarCoeffD11] = f == 1 ? p->directional : p->diffusion;
            if (m_dim == 3) coefficients[StdRegions::eVarCoeffD22] = f == 2 ? p->directional : p->diffusion;
            coefficients[StdRegions::eVarCoeffMass] = p->mass;
            StdRegions::ConstFactorMap factors; factors[StdRegions::eFactorLambda] = 1;
            p->key = std::make_unique<MultiRegions::GlobalLinSysKey>(StdRegions::eHelmholtz,other,factors,coefficients);
            preconditioner = p.get(); m_preconditioners.push_back(std::move(p));
        }
        preconditioners.push_back(preconditioner);
    }
    auto arrays = [&](int count) {
        Fields a(m_dim);
        for (auto &v : a) v = Array<OneD,NekDouble>(count);
        return a;
    };
    auto b = arrays(ng), x = arrays(ng), ax = arrays(ng), r = arrays(ng), z = arrays(ng), p = arrays(ng), ap = arrays(ng);
    auto physical = arrays(nq);
    Array<OneD,NekDouble> divergence(nq),derivative(nq);
    Array<OneD,NekDouble> local(nc), applied(nc), coefficients(nc), counts(ng), countLocal(nc);
    for (int i = 0; i < nc; ++i) countLocal[i] = map->GetLocalToGlobalSign(i);
    map->Assemble(countLocal,counts);
    std::vector<std::pair<MultiRegions::ExpListSharedPtr,Array<OneD,NekDouble>>> saved;
    for (int f = 0; f < m_dim; ++f)
    {
        Vmath::Smul(nc,-1.0,rhs[f],1,local,1);
        fields[f]->ImposeNeumannConditions(local); fields[f]->ImposeRobinConditions(local);
        map->Assemble(local,b[f]);
        fields[f]->ImposeDirichletConditions(fields[f]->UpdateCoeffs());
        map->LocalToGlobal(fields[f]->GetCoeffs(),x[f],true);
        for (int boundary = 0; boundary < fields[f]->GetBndConditions().size(); ++boundary)
            if (fields[f]->GetBndConditions()[boundary]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
            {
                auto exp = fields[f]->GetBndCondExpansions()[boundary];
                Array<OneD,NekDouble> copy(exp->GetNcoeffs());
                Vmath::Vcopy(copy.size(),exp->GetCoeffs(),1,copy,1);
                saved.emplace_back(exp,copy); Vmath::Zero(copy.size(),exp->UpdateCoeffs(),1);
            }
    }
    auto restore = [&]() {
        for (const auto &v : saved) Vmath::Vcopy(v.second.size(),v.second,1,v.first->UpdateCoeffs(),1);
    };
    auto dot = [&](const Fields &a, const Fields &v) {
        double sum = 0;
        for (int f = 0; f < m_dim; ++f)
            for (int i = nd; i < ng; ++i)
                if (map->GetGlobalToUniversalMapUnique(i) > 0) sum += a[f][i]*v[f][i];
        first->GetComm()->AllReduce(sum,LibUtilities::ReduceSum);
        return sum;
    };
    auto apply = [&](const Fields &v, Fields &av) {
        for (int f = 0; f < m_dim; ++f)
        {
            map->GlobalToLocal(v[f],local); fields[f]->BwdTrans(local,physical[f]);
        }
        if (m_gradDiv > 0)
        {
            Vmath::Zero(nq,divergence,1);
            for (int f = 0; f < m_dim; ++f)
            {
                fields[f]->PhysDeriv(MultiRegions::DirCartesianMap[f],physical[f],derivative);
                Vmath::Vadd(nq,derivative,1,divergence,1,divergence,1);
            }
        }
        for (int f = 0; f < m_dim; ++f)
        {
            Fields flux(m_dim);
            for (int d = 0; d < m_dim; ++d)
            {
                flux[d] = Array<OneD,NekDouble>(nq);
                fields[f]->PhysDeriv(MultiRegions::DirCartesianMap[d],physical[f],flux[d]);
                Vmath::Vmul(nq,diffusion[f],1,flux[d],1,flux[d],1);
                if (m_gradDiv > 0 && d == f)
                    Vmath::Svtvp(nq,m_gradDiv,divergence,1,flux[d],1,flux[d],1);
            }
            fields[f]->IProductWRTDerivBase(flux,applied);
            fields[f]->IProductWRTBase(physical[f],coefficients);
            Vmath::Svtvp(nc,alpha,coefficients,1,applied,1,applied,1);
            m_robin(f,physical,applied);
            map->Assemble(applied,av[f]);
        }
    };
    auto precondition = [&]() {
        const double scale = std::sqrt(dot(r,r));
        ASSERTL0(std::isfinite(scale) && scale > 0,"Invalid terrain preconditioner residual");
        for (int f = 0; f < m_dim; ++f)
        {
            map->GlobalToLocal(r[f],local);
            for (int i = 0; i < nc; ++i) local[i] /= scale*counts[map->GetLocalToGlobalMap(i)];
            Vmath::Zero(nc,applied,1);
            fields[f]->GlobalSolve(*preconditioners[f]->key,local,applied);
            map->LocalToGlobal(applied,z[f],true);
            Vmath::Smul(ng,scale,z[f],1,z[f],1); Vmath::Zero(nd,z[f],1);
        }
    };
    int iteration = 0;
    try
    {
        auto lifting = arrays(ng),effective = arrays(ng);
        for (int f = 0; f < m_dim; ++f)
            for (int i = 0; i < ng; ++i) lifting[f][i] = i < nd ? x[f][i] : 0;
        apply(lifting,ax);
        for (int f = 0; f < m_dim; ++f)
            for (int i = 0; i < ng; ++i) effective[f][i] = i < nd ? 0 : b[f][i]-ax[f][i];
        const double norm = std::sqrt(dot(effective,effective));
        const double roundoff = 32*std::numeric_limits<double>::epsilon()*
            (std::sqrt(dot(b,b))+std::sqrt(dot(ax,ax)));
        const double target = std::max(m_tolerance*norm,roundoff);
        if (norm == 0)
            for (int f = 0; f < m_dim; ++f) Vmath::Vcopy(ng,lifting[f],1,x[f],1);
        apply(x,ax);
        for (int f = 0; f < m_dim; ++f)
            for (int i = 0; i < ng; ++i) r[f][i] = i < nd ? 0 : b[f][i]-ax[f][i];
        double rz = 0; bool restart = true;
        for (; iteration < 2000; ++iteration)
        {
            if (std::sqrt(dot(r,r)) <= target)
            {
                apply(x,ax);
                for (int f = 0; f < m_dim; ++f)
                    for (int i = nd; i < ng; ++i) r[f][i] = b[f][i]-ax[f][i];
                if (std::sqrt(dot(r,r)) <= 5*target) break;
                restart = true;
            }
            precondition();
            const double next = dot(r,z);
            ASSERTL0(std::isfinite(next) && next > 0,"Terrain vector preconditioner is not positive definite");
            const double beta = restart ? 0 : next/rz;
            for (int f = 0; f < m_dim; ++f)
                for (int i = 0; i < ng; ++i) p[f][i] = i < nd ? 0 : (restart ? z[f][i] : z[f][i]+beta*p[f][i]);
            rz = next; restart = false; apply(p,ap);
            const double curvature = dot(p,ap);
            ASSERTL0(std::isfinite(curvature) && curvature > 0,"Terrain vector operator is not positive definite");
            const double step = rz/curvature;
            for (int f = 0; f < m_dim; ++f)
                for (int i = nd; i < ng; ++i) { x[f][i] += step*p[f][i]; r[f][i] -= step*ap[f][i]; }
            if (iteration%20 == 19)
            {
                apply(x,ax);
                for (int f = 0; f < m_dim; ++f)
                    for (int i = nd; i < ng; ++i) r[f][i] = b[f][i]-ax[f][i];
                restart = true;
            }
        }
        ASSERTL0(iteration < 2000,"Terrain vector Krylov solve failed to converge");
        for (int f = 0; f < m_dim; ++f)
        {
            map->GlobalToLocal(x[f],fields[f]->UpdateCoeffs());
            fields[f]->BwdTrans(fields[f]->GetCoeffs(),output[f]);
        }
    }
    catch (...) { restore(); throw; }
    restore();
    return iteration;
}
}
