// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_COUPLED_H
#define NEKTAR_WJ_NATIVE_COUPLED_H
#include <SolverUtils/EquationSystem.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupling.h>
#include <fstream>

namespace Nektar
{
class WJNativeCoupled:public SolverUtils::EquationSystem
{
public:
    friend class MemoryManager<WJNativeCoupled>;
    static SolverUtils::EquationSystemSharedPtr create(const LibUtilities::SessionReaderSharedPtr&s,const SpatialDomains::MeshGraphSharedPtr&g)
    {auto eq=MemoryManager<WJNativeCoupled>::AllocateSharedPtr(s,g);eq->InitObject();return eq;}
    static std::string className;
protected:
    WJNativeCoupled(const LibUtilities::SessionReaderSharedPtr&s,const SpatialDomains::MeshGraphSharedPtr&g):EquationSystem(s,g){}
    void v_InitObject(bool declareFields=true)override;
    void v_DoInitialise(bool dumpInitialConditions=true)override;
    void v_DoSolve()override;
    void v_Output()override;
    void v_GenerateSummary(SolverUtils::SummaryList&summary)override;
    void v_EvaluateExactSolution(unsigned int field,Array<OneD,NekDouble>&output,NekDouble time)override;
private:
    std::unique_ptr<WJNative::Coupling>m_coupling;
    WJNative::CoupledHistory m_history;
    WJNative::CoupledStepInfo m_last;
    std::string m_prefix,m_restart,m_initialFields;
    int m_infoStride=20,m_recoveryStride=200,m_scienceStride=200,m_injectRetryStep=-1;
    double m_initialisationSeconds=0,m_constructionSeconds=0,m_completeAdvanceSeconds=0;
    bool m_coefficients=false;
    bool m_steadyMode=false;
    double m_steadyPreconditionerMass=0,m_steadySeed=0;
    WJNative::SteadyKrylovOptions m_steadyOptions;
    WJNative::SteadyKrylovResult m_steadyResult;
    std::vector<std::filesystem::path>m_recovery;
    std::vector<double>m_completeStepSeconds,m_recoverySeconds;
    std::vector<int>m_removedRecovery;
    std::ofstream m_stepsLog;
    static std::string solverTypeLookupId;
};
}
#endif
