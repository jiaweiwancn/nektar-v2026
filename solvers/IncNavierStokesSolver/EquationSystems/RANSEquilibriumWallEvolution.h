///////////////////////////////////////////////////////////////////////////////
// SPDX-License-Identifier: MIT
// Copyright (c) 2026 Nektar++ contributors.
///////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_RANS_EQUILIBRIUM_WALL_EVOLUTION_H
#define NEKTAR_RANS_EQUILIBRIUM_WALL_EVOLUTION_H
#include <SolverUtils/EquationSystem.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallEvolutionPatch.h>
#include <IncNavierStokesSolver/EquationSystems/RANSWallTimeIntegration.h>
namespace Nektar
{
/** Native IMEX evolution of the revised flat multiregion wall element. */
class RANSEquilibriumWallEvolution : public SolverUtils::EquationSystem
{
public:
    friend class MemoryManager<RANSEquilibriumWallEvolution>;
    static SolverUtils::EquationSystemSharedPtr create(const LibUtilities::SessionReaderSharedPtr &session,
        const SpatialDomains::MeshGraphSharedPtr &graph)
    {auto p=MemoryManager<RANSEquilibriumWallEvolution>::AllocateSharedPtr(session,graph);p->InitObject();return p;}
    static std::string className;
protected:
    RANSEquilibriumWallEvolution(const LibUtilities::SessionReaderSharedPtr &session,
        const SpatialDomains::MeshGraphSharedPtr &graph):EquationSystem(session,graph){}
    void v_InitObject(bool declareFields=true) override;
    void v_DoInitialise(bool dump=true) override;
    void v_DoSolve() override;
    void v_Output() override;
    void v_GenerateSummary(SolverUtils::SummaryList &summary) override;
    void v_EvaluateExactSolution(unsigned int field,Array<OneD,NekDouble> &output,NekDouble time) override;
private:
    using Vector=RANSWallAlgebra::Vector;
    RANSWallEvolutionFunction Function(const std::string &name,const std::string &variable,bool required=false) const;
    std::vector<RANSWallEvolutionFunction> Functions(const std::string &name,bool required=true) const;
    Vector Implicit(const Vector &history,double time,double lambda);
    Vector Explicit(const Vector &state,double time);
    void UpdateFields();
    void WriteReport() const;
    void RecordAudit(const Vector &state,const Vector &history,double alpha,double time);
    void Checkpoint();
    void WriteMappedField(const std::string &name,const std::shared_ptr<RANSWallGlobalField> &field,const Vector &values) const;
    std::string PhysicsSignature() const;
    LibUtilities::SessionReaderSharedPtr m_serialSession;
    std::shared_ptr<RANSWallEvolutionPatch> m_patch;
    std::unique_ptr<RANSWallTimeIntegration> m_integrator;
    RANSWallEvolutionData m_data{2},m_implicitData{2};
    Vector m_state,m_initial,m_rawNorms,m_correctedNorms,m_lastHistory;
    std::vector<Vector> m_stageAudits,m_stepAudits,m_explicitAudits;
    double m_eta=0,m_ustar=.4,m_tolerance=2e-12,m_stageTime=0,m_stageAlpha=0;
    std::array<double,2> m_direction{1,0};
    int m_maximumIterations=20,m_order=2;
    bool m_condensed=false,m_neutral=true;
    std::string m_restart,m_checkpointBase,m_lastCheckpoint,m_signature;
    bool m_chartFunctions=false;
    Vector m_roughnessBudget,m_energyBudget,m_finalStageLoads;
    std::vector<Vector> m_regionalStageLoads;
    std::map<int,double> m_finalStageConstraints;
    bool m_diagnostics=false,m_historyReady=false;
    std::string m_failure;
    Vector m_rejectedState;
    std::vector<Vector> m_acceptedStates;
    std::vector<double> m_auditTimes;
    std::vector<RANSWallEvolutionPatch::Diagnostic> m_stateAudits;
    static std::string solverTypeLookupId;
};
}
#endif
