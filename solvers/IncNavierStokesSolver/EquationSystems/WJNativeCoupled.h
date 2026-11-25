// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_COUPLED_H
#define NEKTAR_WJ_NATIVE_COUPLED_H
#include <SolverUtils/EquationSystem.h>
#include <IncNavierStokesSolver/EquationSystems/WJNativeCoupling.h>
#include <fstream>
#include <cstdint>

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
    WJNative::CoupledRecoveryTiming m_recoveryTiming;
    WJNative::CoupledStepInfo m_last;
    std::string m_prefix,m_restart,m_initialFields;
    std::string m_originalSteadyReceipt,m_originalSteadyReceiptSha256;
    double m_originalSteadyKPulse=0.;
    int m_originalSteadyKPulseCell=-1;
    std::optional<double> m_geometrySeedHeight;
    std::string m_perturbationCertificate,m_perturbationSha256;
    int m_infoStride=20,m_recoveryStride=200,m_scienceStride=200,m_injectRetryStep=-1;
    double m_initialisationSeconds=0,m_constructionSeconds=0,m_completeAdvanceSeconds=0;
    bool m_coefficients=false;
    bool m_steadyMode=false;
    bool m_diagnosticSteady=false;
    WJNative::DiagnosticSteadyOptions m_diagnosticSteadyOptions;
    std::optional<WJNative::CoupledState> m_diagnosticSeed;
    double m_steadyPreconditionerMass=0,m_steadySeed=0;
    WJNative::SteadyKrylovOptions m_steadyOptions;
    WJNative::SteadyKrylovResult m_steadyResult;
    std::vector<std::filesystem::path>m_recovery;
    std::vector<double>m_completeStepSeconds,m_recoverySeconds;
    std::vector<int>m_removedRecovery;
    std::ofstream m_stepsLog;
    enum ScientificReason:unsigned
    {
        Periodic=1,BeginLeft=2,BeginRight=4,BeginExact=8,
        EndLeft=16,EndRight=32,EndExact=64,Final=128
    };
    struct ScientificFrame
    {
        std::string path,sha256,origin;
        std::uint64_t step=0,bytes=0;
        double time=0,seconds=0;
        unsigned reasons=0;
    };
    bool m_terminalOutput=false,m_currentProcessAccepted=false,m_solveCompleted=false;
    double m_terminalBegin=0,m_terminalEnd=0,m_initialTime=0;
    double m_scientificFrameSeconds=0,m_scientificIndexSeconds=0;
    std::uint64_t m_initialStep=0;
    std::vector<ScientificFrame>m_scientificFrames;
    unsigned TerminalReasons(double left,double right,bool before)const;
    void WriteScientificFrame(unsigned reasons);
    void WriteScientificIndex();
    static std::string solverTypeLookupId;
};
}
#endif
