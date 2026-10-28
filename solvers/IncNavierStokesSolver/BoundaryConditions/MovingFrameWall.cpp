///////////////////////////////////////////////////////////////////////////////
//
// File: MovingFrameWall.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
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
// Description: Abstract base class for Extrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/MovingFrameWall.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{

std::string MovingFrameWall::className =
    GetIncBCFactory().RegisterCreatorFunction(
        "MovingFrameWall", MovingFrameWall::create,
        "Far field boundary condition of moving reference frame");

MovingFrameWall::MovingFrameWall(
    [[maybe_unused]] const LibUtilities::SessionReaderSharedPtr pSession,
    [[maybe_unused]] Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
    [[maybe_unused]] Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
    [[maybe_unused]] Array<OneD, MultiRegions::ExpListSharedPtr> exp,
    [[maybe_unused]] int nbnd, [[maybe_unused]] int spacedim,
    [[maybe_unused]] int bnddim)
    : StaticWall(pSession, pFields, cond, exp, nbnd, spacedim, bnddim)
{
    classname = "MovingFrameWall";
    m_hasVels = false;
    for (size_t i = 0; i < m_bnddim; ++i)
    {
        if (cond[i]->GetUserDefined() == classname)
        {
            m_hasVels   = true;
            m_BndExp[i] = exp[i];
        }
    }
    m_pressure    = cond.size() - 1;
    m_hasPressure = cond[m_pressure]->GetUserDefined() == classname;
    if (m_hasPressure)
    {
        m_BndExp[m_pressure] = exp[m_pressure];
    }
}

void MovingFrameWall::v_Initialise(
    const LibUtilities::SessionReaderSharedPtr &pSession)
{
    IncBaseCondition::v_Initialise(pSession);
    m_field->GetBndElmtExpansion(m_nbnd, m_bndElmtExps, false);
    if (m_hasPressure)
    {
        m_viscous =
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_intSteps);
        for (int n = 0; n < m_intSteps; ++n)
        {
            m_viscous[n] = Array<OneD, Array<OneD, NekDouble>>(m_bnddim);
            for (int i = 0; i < m_bnddim; ++i)
            {
                m_viscous[n][i] = Array<OneD, NekDouble>(m_npoints, 0.0);
            }
        }
    }
}

/// @brief v_Update set correct BCs (in wavespace)
/// @param fields
/// @param Adv is in wavespace for 3DH1D
/// @param params
void MovingFrameWall::v_Update(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &fields,
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &Adv,
    std::map<std::string, NekDouble> &params)
{
    int nptsPlane0 = 0;
    SetNumPointsOnPlane0(nptsPlane0);
    if (params.find("Omega_x") != params.end() ||
        params.find("Omega_y") != params.end() ||
        params.find("Omega_z") != params.end())
    {
        InitialiseCoords(params);
    }
    // pressure
    if (m_hasPressure && fields.size() > 0)
    {
        ++m_numCalls;

        Array<OneD, Array<OneD, NekDouble>> rhs(m_bnddim);
        for (int i = 0; i < m_bnddim; ++i)
        {
            rhs[i] = Array<OneD, NekDouble>(m_npoints, 0.);
        }
        // add viscous term
        AddVisPressureBCs(fields, rhs, params);
        // add DuDt term
        AddRigidBodyAcc(rhs, params, nptsPlane0);
        m_BndExp[m_pressure]->NormVectorIProductWRTBase(
            rhs, m_BndExp[m_pressure]->UpdateCoeffs());
    }
    // velocity
    if (m_hasVels && nptsPlane0)
    {
        Array<OneD, Array<OneD, NekDouble>> velocities(m_bnddim);
        for (size_t k = 0; k < m_bnddim; ++k)
        {
            if (m_BndExp.find(k) != m_BndExp.end())
            {
                velocities[k] = Array<OneD, NekDouble>(nptsPlane0, 0.0);
            }
        }
        RigidBodyVelocity(velocities, params, nptsPlane0);
        for (int k = 0; k < m_bnddim; ++k)
        {
            if (m_BndExp.find(k) != m_BndExp.end())
            {
                if (m_BndExp[k]->GetExpType() == MultiRegions::e2DH1D)
                {
                    m_BndExp[k]->GetPlane(0)->FwdTransBndConstrained(
                        velocities[k],
                        m_BndExp[k]->GetPlane(0)->UpdateCoeffs());
                }
                else
                {
                    m_BndExp[k]->FwdTransBndConstrained(
                        velocities[k], m_BndExp[k]->UpdateCoeffs());
                }
            }
        }
    }
}

} // namespace Nektar
