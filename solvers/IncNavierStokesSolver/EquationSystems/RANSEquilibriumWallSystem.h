///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_EQUILIBRIUM_WALL_SYSTEM_H
#define NEKTAR_RANS_EQUILIBRIUM_WALL_SYSTEM_H
#include <SolverUtils/EquationSystem.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallPatch.h>
namespace Nektar
{
/** Opt-in native flat wall-element stage solver. M3 solves a steady or a
 * supplied alpha*(state-history) stage. The native IMEX map/history/restart
 * qualification is M4 and is not inferred from these stage solves. */
class RANSEquilibriumWallSystem : public SolverUtils::EquationSystem
{
public:
    friend class MemoryManager<RANSEquilibriumWallSystem>;
    static SolverUtils::EquationSystemSharedPtr create(const LibUtilities::SessionReaderSharedPtr &session,
        const SpatialDomains::MeshGraphSharedPtr &graph)
    {
        auto p=MemoryManager<RANSEquilibriumWallSystem>::AllocateSharedPtr(session,graph);
        p->InitObject();return p;
    }
    static std::string className;
protected:
    RANSEquilibriumWallSystem(const LibUtilities::SessionReaderSharedPtr &session,
        const SpatialDomains::MeshGraphSharedPtr &graph):EquationSystem(session,graph){}
    void v_InitObject(bool declareFields=true) override;
    void v_DoInitialise(bool dumpInitialConditions=true) override;
    void v_DoSolve() override;
    void v_Output() override;
    void v_GenerateSummary(SolverUtils::SummaryList &summary) override;
    void v_EvaluateExactSolution(unsigned int field,Array<OneD,NekDouble> &output,NekDouble time) override;
private:
    RANSWallTimeFunction Function(const std::string &name,const std::string &variable,bool required=false) const;
    std::array<RANSWallTimeFunction,6> Functions(const std::string &name,bool required) const;
    void UpdateFields();
    std::shared_ptr<RANSWallFieldLayout> m_layout;
    std::shared_ptr<RANSWallPatch> m_patch;
    RANSWallPatchData m_data;
    RANSWallAlgebra::Vector m_state,m_history,m_residual,m_norms;
    double m_alpha=0,m_tolerance=2e-12;
    int m_iterations=0,m_maximumIterations=15;
    bool m_condensed=false;
};
}
#endif
