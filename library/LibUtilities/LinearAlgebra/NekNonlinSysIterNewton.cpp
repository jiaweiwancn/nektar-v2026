///////////////////////////////////////////////////////////////////////////////
//
// File: NekNonlinSysIterNewton.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: NekNonlinSysIterNewton definition
//
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/LinearAlgebra/NekNonlinSysIterNewton.h>

using namespace std;

namespace Nektar::LibUtilities
{
/**
 * @class  NekNonlinSysIterNewton
 *
 * Solves a nonlinear system using iterative methods.
 */
string NekNonlinSysIterNewton::className =
    LibUtilities::GetNekNonlinSysIterFactory().RegisterCreatorFunction(
        "Newton", NekNonlinSysIterNewton::create,
        "NekNonlinSysIterNewton solver.");

NekNonlinSysIterNewton::NekNonlinSysIterNewton(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const LibUtilities::CommSharedPtr &vRowComm, const int nscale,
    const NekSysKey &pKey)
    : NekNonlinSysIter(pSession, vRowComm, nscale, pKey)
{
}

void NekNonlinSysIterNewton::v_InitObject()
{
    NekNonlinSysIter::v_InitObject();
}

/**
 *
 **/
int NekNonlinSysIterNewton::v_SolveSystem(
    const int nGlobal, const Array<OneD, const NekDouble> &pInput,
    Array<OneD, NekDouble> &pOutput, [[maybe_unused]] const int nDir)
{
    ASSERTL0(nGlobal == m_SysDimen, "nGlobal != m_SysDimen");

    m_SourceVec     = pInput;
    m_Solution      = pOutput;
    m_NtotLinSysIts = 0;
    m_converged     = false;

    Vmath::Vcopy(nGlobal, pInput, 1, m_Solution, 1);

    NekDouble resnormOld = 0.0;
    int NttlNonlinIte    = 0;
    for (; NttlNonlinIte < m_NekNonlinSysMaxIterations; ++NttlNonlinIte)
    {
        m_operator.DoNekSysResEval(m_Solution, m_Residual, true);

        ConvergenceCheck(NttlNonlinIte, m_Residual);
        if (m_converged)
        {
            break;
        }

        NekDouble LinSysRelativeIteTol =
            CalcInexactNewtonForcing(NttlNonlinIte, resnormOld, m_SysResNorm);
        resnormOld = m_SysResNorm;
        m_linsol->SetRhsMagnitude(m_SysResNorm);
        m_linsol->SetNekLinSysTolerance(LinSysRelativeIteTol);
        int ntmpLinSysIts =
            m_linsol->SolveSystem(nGlobal, m_Residual, m_DeltSltn, 0);
        m_NtotLinSysIts += ntmpLinSysIts;

        Vmath::Vsub(nGlobal, m_Solution, 1, m_DeltSltn, 1, m_Solution, 1);
    }

    if ((!m_converged || m_verbose) && m_root && m_FlagWarnings)
    {
        int nwidthcolm = 11;

        WARNINGL0(m_converged,
                  "     # Nonlinear solver not converge in DoImplicitSolve");
        cout << right << scientific << setw(nwidthcolm)
             << setprecision(nwidthcolm - 6)
             << "     * Newton-Its converged (RES=" << sqrt(m_SysResNorm)
             << " Res/Res0= " << sqrt(m_SysResNorm / m_SysResNorm0)
             << " Res/DtRHS= " << sqrt(m_SysResNorm / m_rhs_magnitude)
             << " with " << setw(3) << NttlNonlinIte << " Non-Its)" << endl;
    }

    return NttlNonlinIte;
}

NekDouble NekNonlinSysIterNewton::CalcInexactNewtonForcing(
    const int &nIteration, const NekDouble &resnormOld,
    const NekDouble &resnorm)
{
    if (nIteration == 0 || !m_InexactNewtonForcing)
    {
        return m_NekLinSysTolerance;
    }
    else
    {
        static const NekDouble forcingGamma = 1.0;
        static const NekDouble forcingAlpha = 0.5 * (1.0 + sqrt(5.0));
        NekDouble tmpForc =
            forcingGamma * pow((resnorm / resnormOld), forcingAlpha);
        return max(min(m_NekLinSysTolerance, tmpForc), 1.0E-6);
    }
}

} // namespace Nektar::LibUtilities
