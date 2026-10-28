///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMovingReferenceFrame.cpp
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
// Description: Solving the absolute flow in a moving body frame,
// by adding (U0 + Omega X (x - x0)) . grad u - Omega X u
// as the body force.
// U0 is the translational velocity of the body frame.
// Omega is the angular velocity.
// x0 is the rotation pivot in the body frame.
// All vectors use the basis of the body frame.
// Translational motion is allowed for all dimensions.
// Rotation is not allowed for 1D, 2DH1D, 3DH2D.
// Rotation in z direction is allowed for 2D and 3DH1D.
// Rotation in 3 directions are allowed for 3D.
// TODO: add suport for 3D rotation using Quaternion
///////////////////////////////////////////////////////////////////////////////

#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/Vmath.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <SolverUtils/Forcing/ForcingMovingReferenceFrame.h>
#include <boost/format.hpp>

namespace Nektar::SolverUtils
{

std::string ForcingMovingReferenceFrame::classNameBody =
    GetForcingFactory().RegisterCreatorFunction(
        "MovingReferenceFrame", ForcingMovingReferenceFrame::create,
        "Moving Frame");

/**
 * @brief
 * @param pSession
 * @param pEquation
 */
ForcingMovingReferenceFrame::ForcingMovingReferenceFrame(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const std::weak_ptr<EquationSystem> &pEquation)
    : Forcing(pSession, pEquation)
{
}

ForcingMovingReferenceFrame::~ForcingMovingReferenceFrame(void)
{
    if (m_isRoot)
    {
        m_outputStream.close();
    }
}

/**
 * @brief Initialise the forcing module
 * @param pFields
 * @param pNumForcingFields
 * @param pForce
 */
void ForcingMovingReferenceFrame::v_InitObject(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    [[maybe_unused]] const unsigned int &pNumForcingFields,
    const TiXmlElement *pForce)
{
    m_session->MatchSolverInfo("Homogeneous", "1D", m_isH1d, false);
    m_session->MatchSolverInfo("Homogeneous", "2D", m_isH2d, false);
    bool singleMode, halfMode;
    m_session->MatchSolverInfo("ModeType", "SingleMode", singleMode, false);
    m_session->MatchSolverInfo("ModeType", "HalfMode", halfMode, false);
    if (singleMode || halfMode)
    {
        m_isH1d = false;
    }
    m_expdim    = m_isH2d ? 1 : pFields[0]->GetGraph()->GetMeshDimension();
    m_spacedim  = m_expdim + (m_isH1d ? 1 : 0) + (m_isH2d ? 2 : 0);
    m_isRoot    = pFields[0]->GetComm()->TreatAsRankZero();
    m_hasPlane0 = true;
    if (m_isH1d)
    {
        m_hasPlane0 = pFields[0]->GetZIDs()[0] == 0;
    }
    // initialize variables
    m_velxyz           = Array<OneD, NekDouble>(3, 0.0);
    m_omegaxyz         = Array<OneD, NekDouble>(3, 0.0);
    m_body.extForceXYZ = Array<OneD, NekDouble>(6, 0.0);
    m_hasVel           = Array<OneD, bool>(3, false);
    m_hasOmega         = Array<OneD, bool>(3, false);
    m_travelWave       = Array<OneD, NekDouble>(3, 0.);
    m_pivotPoint       = Array<OneD, NekDouble>(3, 0.0);
    m_index            = 0;
    m_currentTime      = -1.;
    LoadParameters(pForce);
    InitBodySolver(pForce);
    if (m_body.isCircular || m_expdim == 1)
    {
        for (int i = 0; i < 3; ++i)
        {
            m_hasOmega[i] = false;
        }
        m_hasRotation = false;
    }
    // account for the effect of rotation
    // Omega_X results in having v and w even if not defined by user
    // Omega_Y results in having u and w even if not defined by user
    // Omega_Z results in having u and v even if not defined by user
    for (int i = 0; i < 3; ++i)
    {
        int j = (i + 1) % 3;
        int k = (i + 2) % 3;
        if (m_hasOmega[i])
        {
            m_hasVel[j] = true;
            m_hasVel[k] = true;
        }
    }
    if (m_hasRotation)
    {
        int npoints = pFields[0]->GetNpoints();
        m_coords    = Array<OneD, Array<OneD, NekDouble>>(3);
        for (int j = 0; j < m_spacedim; ++j)
        {
            m_coords[j] = Array<OneD, NekDouble>(npoints);
        }
        pFields[0]->GetCoords(m_coords[0], m_coords[1], m_coords[2]);
        // move the origin to the pivot point
        for (int i = 0; i < m_spacedim; ++i)
        {
            Vmath::Sadd(npoints, -m_pivotPoint[i], m_coords[i], 1, m_coords[i],
                        1);
        }
    }
    // initialise pivot point for fluid interface
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer to the equation system is expired");
    auto FluidEq = std::dynamic_pointer_cast<FluidInterface>(equ);
    FluidEq->SetMovingFramePivot(m_pivotPoint);
    // initialize aeroforce filter
    if (m_body.hasFreeMotion)
    {
        InitialiseFilter(m_session, pFields, pForce);
    }
}

NekDouble ForcingMovingReferenceFrame::EvaluateExpression(
    std::string expression)
{
    NekDouble value = 0.;
    try
    {
        LibUtilities::Equation expession(m_session->GetInterpreter(),
                                         expression);
        value = expession.Evaluate();
    }
    catch (const std::runtime_error &)
    {
        NEKERROR(ErrorUtil::efatal, "Error evaluating expression" + expression);
    }
    return value;
}

void ForcingMovingReferenceFrame::LoadParameters(const TiXmlElement *pForce)
{
    const TiXmlElement *funcNameElmt;
    // body is a circular cylinder
    const TiXmlElement *mssgTagSpecial =
        pForce->FirstChildElement("CIRCULARCYLINDER");
    if (mssgTagSpecial)
    {
        m_body.isCircular = true;
    }
    else
    {
        m_body.isCircular = false;
    }
    // load frame velocity
    funcNameElmt = pForce->FirstChildElement("FRAMEVELOCITY");
    if (funcNameElmt)
    {
        std::string FuncName = funcNameElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(FuncName),
                 "Function '" + FuncName + "' is not defined in the session.");
        // linear velocity
        for (int i = 0; i < m_spacedim; ++i)
        {
            std::string var = m_session->GetVariable(i);
            if (m_session->DefinesFunction(FuncName, var))
            {
                m_frameVelFunction[i] = m_session->GetFunction(FuncName, var);
                m_hasVel[i]           = true;
            }
        }
        // linear displacement and acceleration
        std::vector<std::string> linearDispVar = {"X", "Y", "Z"};
        std::vector<std::string> linearAcceVar = {"A_x", "A_y", "A_z"};
        for (int i = 0; i < m_spacedim; ++i)
        {
            if (m_session->DefinesFunction(FuncName, linearDispVar[i]))
            {
                m_frameVelFunction[i + 6] =
                    m_session->GetFunction(FuncName, linearDispVar[i]);
            }
            if (m_session->DefinesFunction(FuncName, linearAcceVar[i]))
            {
                m_frameVelFunction[i + 12] =
                    m_session->GetFunction(FuncName, linearAcceVar[i]);
            }
        }
        // angular velocities
        m_hasRotation                       = false;
        std::vector<std::string> angularVar = {"Omega_x", "Omega_y", "Omega_z"};
        for (int i = 0; i < 3; ++i)
        {
            std::string var = angularVar[i];
            if (m_session->DefinesFunction(FuncName, var))
            {
                m_frameVelFunction[i + 3] =
                    m_session->GetFunction(FuncName, var);
                m_hasOmega[i] = true;
                m_hasRotation = true;
            }
        }
        // angular displacement and acceleration
        std::vector<std::string> angularDispVar  = {"Theta_x", "Theta_y",
                                                    "Theta_z"};
        std::vector<std::string> angularAccepVar = {"DOmega_x", "DOmega_y",
                                                    "DOmega_z"};
        for (int i = 0; i < 3; ++i)
        {
            if (m_session->DefinesFunction(FuncName, angularDispVar[i]))
            {
                m_frameVelFunction[i + 3 + 6] =
                    m_session->GetFunction(FuncName, angularDispVar[i]);
            }
            if (m_session->DefinesFunction(FuncName, angularAccepVar[i]))
            {
                m_frameVelFunction[i + 3 + 12] =
                    m_session->GetFunction(FuncName, angularAccepVar[i]);
            }
        }
        // TODO: add the support for general rotation
        for (int i = 0; i < 2; ++i)
        {
            ASSERTL0(!m_hasOmega[i], "Currently only Omega_z is supported");
        }
    }

    // load external force
    funcNameElmt = pForce->FirstChildElement("EXTERNALFORCE");
    if (funcNameElmt)
    {
        std::string FuncName = funcNameElmt->GetText();
        if (m_session->DefinesFunction(FuncName))
        {
            std::vector<std::string> forceVar = {"Fx", "Fy", "Fz"};
            for (int i = 0; i < m_spacedim; ++i)
            {
                std::string var = forceVar[i];
                if (m_session->DefinesFunction(FuncName, var))
                {
                    m_body.extForceFunction[i] =
                        m_session->GetFunction(FuncName, var);
                }
            }
            std::vector<std::string> momentVar = {"Mx", "My", "Mz"};
            for (int i = 0; i < 3; ++i)
            {
                std::string var = momentVar[i];
                if (m_session->DefinesFunction(FuncName, var))
                {
                    m_body.extForceFunction[i + 3] =
                        m_session->GetFunction(FuncName, var);
                }
            }
        }
    }

    // load pitching pivot
    const TiXmlElement *pivotElmt = pForce->FirstChildElement("PIVOTPOINT");
    if (pivotElmt)
    {
        std::vector<std::string> values;
        std::string mssgStr = pivotElmt->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        for (int j = 0; j < m_spacedim; ++j)
        {
            m_pivotPoint[j] = EvaluateExpression(values[j]);
        }
    }

    // load travelling wave speed
    const TiXmlElement *TWSElmt =
        pForce->FirstChildElement("TRAVELINGWAVESPEED");
    if (TWSElmt)
    {
        std::vector<std::string> values;
        std::string mssgStr = TWSElmt->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        for (int j = 0; j < m_spacedim; ++j)
        {
            NekDouble tmp = EvaluateExpression(values[j]);
            if (fabs(tmp) > NekConstants::kNekMachineEpsilon)
            {
                m_travelWave[j] = tmp;
                m_hasVel[j]     = true;
            }
        }
    }

    // OutputFile
    const TiXmlElement *mssgTag = pForce->FirstChildElement("OutputFile");
    std::string filename;
    if (mssgTag)
    {
        filename = mssgTag->GetText();
    }
    else
    {
        filename = m_session->GetSessionName();
    }
    if (!(filename.length() >= 4 &&
          filename.substr(filename.length() - 4) == ".mrf"))
    {
        filename += ".mrf";
    }
    if (m_isRoot)
    {
        m_outputStream.open(filename.c_str());
        if (m_spacedim == 2)
        {
            m_outputStream
                << "Variables = t, x, ux, ax, y, uy, ay, theta, omega, domega"
                << std::endl;
        }
        else if (m_spacedim == 3)
        {
            m_outputStream << "Variables = t, x, ux, ax, y, uy, ay, z, uz, az, "
                              "theta, omega, domega"
                           << std::endl;
        }
    }

    // output frequency
    m_outputFrequency = 1;
    mssgTag           = pForce->FirstChildElement("OutputFrequency");
    if (mssgTag)
    {
        std::vector<std::string> values;
        std::string mssgStr = mssgTag->GetText();
        m_outputFrequency   = round(EvaluateExpression(mssgStr));
    }
    ASSERTL0(m_outputFrequency > 0,
             "OutputFrequency should be greater than zero.");
}

void ForcingMovingReferenceFrame::InitBodySolver(const TiXmlElement *pForce)
{
    int NumDof = m_spacedim + 1;
    const TiXmlElement *mssgTag;
    std::string mssgStr;
    // allocate memory and initialise
    m_body.vel = Array<OneD, Array<OneD, NekDouble>>(3);
    for (size_t i = 0; i < 3; ++i)
    {
        m_body.vel[i] = Array<OneD, NekDouble>(NumDof, 0.);
    }
    SetInitialConditions();
    UpdateFluidInterface(m_body.vel, 1);
    // read free motion DoFs
    m_body.dirDoFs.clear();
    for (int i = 0; i < NumDof; ++i)
    {
        m_body.dirDoFs.insert(i);
    }
    mssgTag = pForce->FirstChildElement("MOTIONPRESCRIBED");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == NumDof,
                 "MOTIONPRESCRIBED vector should be of size " +
                     std::to_string(NumDof));
        for (int i = 0; i < NumDof; ++i)
        {
            if (EvaluateExpression(values[i]) == 0)
            {
                m_body.dirDoFs.erase(i);
                if (i < m_spacedim)
                {
                    m_hasVel[i] = true;
                }
                else if (i == m_spacedim)
                {
                    m_hasOmega[2] = true;
                    m_hasRotation = true;
                }
            }
        }
    }
    m_body.hasFreeMotion = m_body.dirDoFs.size() < NumDof;
    // read mass matrix
    m_body.M = Array<OneD, NekDouble>(NumDof * NumDof, 0.);
    mssgTag  = pForce->FirstChildElement("MASS");
    ASSERTL0(m_body.dirDoFs.size() == NumDof || mssgTag,
             "Mass matrix is required.");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == NumDof * NumDof,
                 "Mass matrix should be of size " + std::to_string(NumDof) +
                     "X" + std::to_string(NumDof));
        int count = 0;
        for (int i = 0; i < NumDof; ++i)
        {
            for (int j = 0; j < NumDof; ++j)
            {
                m_body.M[count] = EvaluateExpression(values[count]);
                ++count;
            }
        }
    }
    // read damping matrix
    m_body.C = Array<OneD, NekDouble>(NumDof * NumDof, 0.);
    mssgTag  = pForce->FirstChildElement("DAMPING");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == NumDof * NumDof,
                 "Damping matrix should be of size " + std::to_string(NumDof) +
                     "X" + std::to_string(NumDof));
        int count = 0;
        for (int i = 0; i < NumDof; ++i)
        {
            for (int j = 0; j < NumDof; ++j)
            {
                m_body.C[count] = EvaluateExpression(values[count]);
                ++count;
            }
        }
    }
    // read rigidity matrix
    m_body.K = Array<OneD, NekDouble>(NumDof * NumDof, 0.);
    mssgTag  = pForce->FirstChildElement("RIGIDITY");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == NumDof * NumDof,
                 "Rigidity matrix should be of size " + std::to_string(NumDof) +
                     "X" + std::to_string(NumDof));
        int count = 0;
        for (int i = 0; i < NumDof; ++i)
        {
            for (int j = 0; j < NumDof; ++j)
            {
                m_body.K[count] = EvaluateExpression(values[count]);
                ++count;
            }
        }
    }
    // read Newmark Beta paramters
    m_timestep      = m_session->GetParameter("TimeStep");
    NekDouble beta  = 0.25;
    NekDouble gamma = 0.51;
    if (m_session->DefinesParameter("NewmarkBeta"))
    {
        beta = m_session->GetParameter("NewmarkBeta");
    }
    if (m_session->DefinesParameter("NewmarkGamma"))
    {
        gamma = m_session->GetParameter("NewmarkGamma");
    }
    m_bodySolver.SetNewmarkBeta(beta, gamma, m_timestep, m_body.M, m_body.C,
                                m_body.K, m_body.dirDoFs);
}

void ForcingMovingReferenceFrame::UpdatePrescribed(
    const NekDouble &time, std::map<int, NekDouble> &Dirs)
{
    int NumDof = m_spacedim + 1;
    for (auto it : m_frameVelFunction)
    {
        if (it.first < 3)
        {
            Dirs[it.first] = it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first == 5)
        {
            Dirs[m_spacedim] = it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first < 9)
        {
            Dirs[NumDof + it.first - 6] = it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first == 11)
        {
            Dirs[NumDof + m_spacedim] = it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first < 15)
        {
            Dirs[(NumDof << 1) + it.first - 12] =
                it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first == 17)
        {
            Dirs[(NumDof << 1) + m_spacedim] =
                it.second->Evaluate(0., 0., 0., time);
        }
    }
    for (auto i : m_body.dirDoFs)
    {
        if (Dirs.find(i) == Dirs.end())
        {
            Dirs[i]                 = 0.;
            Dirs[i + NumDof]        = 0.;
            Dirs[i + (NumDof << 1)] = 0.;
        }
    }
}
/**
 * @brief Updates the forcing array with the current required forcing.
 * @param pFields
 * @param time
 */
void ForcingMovingReferenceFrame::UpdateFrameVelocity(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const NekDouble &time)
{
    if (m_currentTime >= time)
    {
        return;
    }
    m_currentTime = time;
    std::map<int, NekDouble> Dirs;
    UpdatePrescribed(time, Dirs);
    if (m_index == 0)
    {
        SetInitialConditions(Dirs);
    }
    else
    {
        // compute the velocites whoes functions are provided in inertial frame
        Array<OneD, NekDouble> forcebody(6, 0.); // fluid force
        if (m_body.hasFreeMotion)
        {
            for (auto it : m_body.extForceFunction)
            {
                m_body.extForceXYZ[it.first] =
                    it.second->Evaluate(0., 0., 0., time);
            }
            m_body.aeroforceFilter->GetForces(pFields, NullNekDouble1DArray,
                                              time);
            auto equ = m_equ.lock();
            ASSERTL0(equ, "Weak pointer to the equation system is expired");
            auto FluidEq = std::dynamic_pointer_cast<FluidInterface>(equ);
            FluidEq->GetAeroForce(forcebody);
        }
        SolveBodyMotion(m_body.vel, forcebody, Dirs);
    }
    if (m_isRoot && m_index % m_outputFrequency == 0)
    {
        m_outputStream << boost::format("%25.19e") % time << " ";
        for (size_t i = 0; i < m_body.vel[0].size(); ++i)
        {
            m_outputStream << boost::format("%25.19e") % m_body.vel[0][i] << " "
                           << boost::format("%25.19e") % m_body.vel[1][i] << " "
                           << boost::format("%25.19e") % m_body.vel[2][i]
                           << " ";
        }
        m_outputStream << std::endl;
    }
    // extract values and transform to body frame
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_velxyz[i] = m_body.vel[1][i];
    }
    if (m_hasRotation)
    {
        m_omegaxyz[2] = m_body.vel[1][m_spacedim];
        Array<OneD, NekDouble> angle(3, 0.);
        angle[2] = m_body.vel[0][m_spacedim];
        m_frame.SetAngle(angle);
        m_frame.IneritalToBody(3, m_velxyz, m_velxyz);
    }
    // set displacements at the current time step
    UpdateFluidInterface(m_body.vel, 0);
    ++m_index;
}

void ForcingMovingReferenceFrame::UpdateFluidInterface(
    Array<OneD, Array<OneD, NekDouble>> &bodyVel, const int step)
{
    // set displacements at the current or next time step
    Array<OneD, NekDouble> disp(6, 0.);
    Array<OneD, NekDouble> vel(12, 0.0);
    for (int i = 0; i < m_spacedim; ++i)
    {
        disp[i] = bodyVel[0][i];
    }
    if (!m_body.isCircular)
    {
        disp[5] = bodyVel[0][m_spacedim];
    }
    // to set the boundary condition of the next time step
    // update the frame velocities and accelerations
    for (int i = 0; i < m_spacedim; ++i)
    {
        vel[i]     = bodyVel[1][i];
        vel[i + 6] = bodyVel[2][i];
    }
    vel[5]  = bodyVel[1][m_spacedim];
    vel[11] = bodyVel[2][m_spacedim];
    Array<OneD, NekDouble> tmp;
    if (step && m_hasRotation)
    {
        m_frame.SetAngle(disp + 3);
        m_frame.IneritalToBody(3, vel, vel);
        m_frame.IneritalToBody(3, vel + 6, tmp = vel + 6);
    }
    auto equ = m_equ.lock();
    ASSERTL0(equ, "Weak pointer to the equation system is expired");
    auto FluidEq = std::dynamic_pointer_cast<FluidInterface>(equ);
    FluidEq->SetMovingFrameVelocities(vel, step);
    FluidEq->SetMovingFrameDisp(disp, step);
}

void ForcingMovingReferenceFrame::SolveBodyMotion(
    Array<OneD, Array<OneD, NekDouble>> &bodyVel,
    const Array<OneD, NekDouble> &forcebody, std::map<int, NekDouble> &Dirs)
{
    if (!m_body.hasFreeMotion)
    {
        m_bodySolver.SolvePrescribed(bodyVel, Dirs);
    }
    else if (!m_hasRotation || Dirs.find(m_spacedim) != Dirs.end())
    {
        m_bodySolver.SolvePrescribed(bodyVel, Dirs);
        Array<OneD, NekDouble> force(6, 0.), tmp;
        if (m_hasRotation && Dirs.find(m_spacedim) != Dirs.end())
        {
            Array<OneD, NekDouble> angle(3, 0.);
            angle[2] = bodyVel[0][m_spacedim];
            m_frame.SetAngle(angle);
            m_frame.BodyToInerital(3, forcebody, force);
            m_frame.BodyToInerital(3, forcebody + 3, tmp = force + 3);
        }
        else
        {
            Vmath::Vcopy(6, forcebody, 1, force, 1);
        }
        for (int i = 0; i < m_spacedim; ++i)
        {
            force[i] += m_body.extForceXYZ[i];
        }
        force[m_spacedim] = force[5] + m_body.extForceXYZ[5];
        m_bodySolver.SolveFree(bodyVel, force);
    }
    else
    {
        Array<OneD, Array<OneD, NekDouble>> tmpbodyVel(bodyVel.size());
        for (size_t i = 0; i < bodyVel.size(); ++i)
        {
            tmpbodyVel[i] = Array<OneD, NekDouble>(bodyVel[i].size());
        }
        Array<OneD, NekDouble> angle(3, 0.);
        angle[2] = bodyVel[0][m_spacedim];
        for (int iter = 0; iter < 2; ++iter)
        {
            // copy initial condition
            for (size_t i = 0; i < bodyVel.size(); ++i)
            {
                Vmath::Vcopy(bodyVel[i].size(), bodyVel[i], 1, tmpbodyVel[i],
                             1);
            }
            Array<OneD, NekDouble> force(6, 0.), tmp;
            m_frame.SetAngle(angle);
            m_frame.BodyToInerital(3, forcebody, force);
            m_frame.BodyToInerital(3, forcebody + 3, tmp = force + 3);
            for (int i = 0; i < m_spacedim; ++i)
            {
                force[i] += m_body.extForceXYZ[i];
            }
            force[m_spacedim] = force[5] + m_body.extForceXYZ[5];
            m_bodySolver.Solve(tmpbodyVel, force, Dirs);
            angle[2] = tmpbodyVel[0][m_spacedim];
        }
        // copy final results
        for (size_t i = 0; i < bodyVel.size(); ++i)
        {
            Vmath::Vcopy(bodyVel[i].size(), tmpbodyVel[i], 1, bodyVel[i], 1);
        }
    }
}

/**
 * @brief Adds the body force, -Omega X u.
 * @param fields
 * @param inarray
 * @param outarray
 * @param time
 */
void ForcingMovingReferenceFrame::v_Apply(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time)
{
    // If there is no rotation, body force is zero,
    // nothing needs to be done here.
    if (m_hasRotation)
    {
        int npoints = fields[0]->GetNpoints();
        addRotation(npoints, outarray, -1., inarray, outarray);
    }
    UpdateBoundaryConditions(time);
}

/**
 * @brief Set velocity boundary condition at the next time step
 */
void ForcingMovingReferenceFrame::UpdateBoundaryConditions(NekDouble time)
{
    time += m_timestep;
    // compute the velocities whose functions are provided in inertial frame
    std::map<int, NekDouble> Dirs;
    UpdatePrescribed(time, Dirs);
    Array<OneD, Array<OneD, NekDouble>> bodyVel(m_body.vel.size());
    for (size_t i = 0; i < m_body.vel.size(); ++i)
    {
        bodyVel[i] = Array<OneD, NekDouble>(m_body.vel[i].size());
        Vmath::Vcopy(m_body.vel[i].size(), m_body.vel[i], 1, bodyVel[i], 1);
    }
    m_bodySolver.SolvePrescribed(bodyVel, Dirs);
    // set displacements at the next time step
    UpdateFluidInterface(bodyVel, 1);
}

/**
 * @brief outarray = inarray0 + angVelScale Omega x inarray1
 */
void ForcingMovingReferenceFrame::addRotation(
    int nPnts, // number of points
    const Array<OneD, Array<OneD, NekDouble>> &inarray0, NekDouble angVelScale,
    const Array<OneD, Array<OneD, NekDouble>> &inarray1,
    Array<OneD, Array<OneD, NekDouble>> &outarray)
{
    ASSERTL0(&inarray1 != &outarray, "inarray1 and outarray "
                                     "should not be the same.");

    // TODO: In case of having support for all three components of Omega,
    // they should be transformed into the rotating frame first!

    // In case that the inarray0 and outarry are different, to avoid using
    // un-initialized array, copy the array first
    if (&inarray0 != &outarray)
    {
        ASSERTL0(inarray0.size() == outarray.size(),
                 "inarray0 and outarray must have same dimentions");
        for (int i = 0; i < inarray0.size(); ++i)
        {
            Vmath::Vcopy(nPnts, inarray0[i], 1, outarray[i], 1);
        }
    }

    if (m_spacedim >= 2 && m_hasOmega[2])
    {
        NekDouble cp = m_omegaxyz[2] * angVelScale;
        NekDouble cm = -1. * cp;

        Vmath::Svtvp(nPnts, cm, inarray1[1], 1, outarray[0], 1, outarray[0], 1);
        Vmath::Svtvp(nPnts, cp, inarray1[0], 1, outarray[1], 1, outarray[1], 1);
    }

    if (m_spacedim == 3 && m_hasOmega[0])
    {
        NekDouble cp = m_omegaxyz[0] * angVelScale;
        NekDouble cm = -1. * cp;

        Vmath::Svtvp(nPnts, cp, inarray1[1], 1, outarray[2], 1, outarray[2], 1);
        Vmath::Svtvp(nPnts, cm, inarray1[2], 1, outarray[1], 1, outarray[1], 1);
    }

    if (m_spacedim == 3 && m_hasOmega[1])
    {
        NekDouble cp = m_omegaxyz[1] * angVelScale;
        NekDouble cm = -1. * cp;

        Vmath::Svtvp(nPnts, cp, inarray1[2], 1, outarray[0], 1, outarray[0], 1);
        Vmath::Svtvp(nPnts, cm, inarray1[0], 1, outarray[2], 1, outarray[2], 1);
    }
}

/**
 * @brief Compute the moving frame velocity at given time
 */
void ForcingMovingReferenceFrame::v_PreApply(
    const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
    const Array<OneD, Array<OneD, NekDouble>> &inarray,
    Array<OneD, Array<OneD, NekDouble>> &outarray, const NekDouble &time)
{
    UpdateFrameVelocity(fields, time);
    int npoints = fields[0]->GetNpoints();
    if (m_isH2d && fields[0]->GetWaveSpace())
    {
        for (int i = 0; i < m_spacedim; ++i)
        {
            if (m_hasVel[i])
            {
                Array<OneD, NekDouble> tmpphys(npoints,
                                               -m_velxyz[i] - m_travelWave[i]);
                Array<OneD, NekDouble> tmpcoef(npoints);
                fields[0]->HomogeneousFwdTrans(npoints, tmpphys, tmpcoef);
                Vmath::Vadd(npoints, tmpcoef, 1, inarray[i], 1, outarray[i], 1);
            }
            else if (&inarray != &outarray)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
        }
    }
    else
    {
        int npoints0 = npoints;
        if (m_isH1d && fields[0]->GetWaveSpace())
        {
            npoints0 = m_hasPlane0 ? fields[0]->GetPlane(0)->GetNpoints() : 0;
        }
        for (int i = 0; i < m_spacedim; ++i)
        {
            if (m_hasVel[i])
            {
                Vmath::Sadd(npoints0, -m_velxyz[i] - m_travelWave[i],
                            inarray[i], 1, outarray[i], 1);
                if (&inarray != &outarray && npoints != npoints0)
                {
                    Array<OneD, NekDouble> tmp = outarray[i] + npoints0;
                    Vmath::Vcopy(npoints - npoints0, inarray[i] + npoints0, 1,
                                 tmp, 1);
                }
            }
            else if (&inarray != &outarray)
            {
                Vmath::Vcopy(npoints, inarray[i], 1, outarray[i], 1);
            }
        }
        if (m_hasRotation)
        {
            addRotation(npoints0, outarray, -1., m_coords, outarray);
        }
    }
}

void ForcingMovingReferenceFrame::SetInitialConditions()
{
    NekDouble time = 0.;
    std::map<std::string, std::string> fieldMetaDataMap;
    std::vector<std::string> strFrameData = {
        "X",   "Y",   "Z",   "Theta_x",  "Theta_y",  "Theta_z",
        "U",   "V",   "W",   "Omega_x",  "Omega_y",  "Omega_z",
        "A_x", "A_y", "A_z", "DOmega_x", "DOmega_y", "DOmega_z"};
    std::map<std::string, NekDouble> fileData;
    if (m_session->DefinesFunction("InitialConditions"))
    {
        for (int i = 0; i < m_session->GetVariables().size(); ++i)
        {
            if (m_session->GetFunctionType("InitialConditions",
                                           m_session->GetVariable(i)) ==
                LibUtilities::eFunctionTypeFile)
            {
                std::string filename = m_session->GetFunctionFilename(
                    "InitialConditions", m_session->GetVariable(i));
                fs::path pfilename(filename);
                // redefine path for parallel file which is in directory
                if (fs::is_directory(pfilename))
                {
                    fs::path metafile("Info.xml");
                    fs::path fullpath = pfilename / metafile;
                    filename          = LibUtilities::PortablePath(fullpath);
                }
                LibUtilities::FieldIOSharedPtr fld =
                    LibUtilities::FieldIO::CreateForFile(m_session, filename);
                fld->ImportFieldMetaData(filename, fieldMetaDataMap);

                // check to see if time is defined
                if (fieldMetaDataMap != LibUtilities::NullFieldMetaDataMap)
                {
                    if (fieldMetaDataMap.find("Time") != fieldMetaDataMap.end())
                    {
                        time = std::stod(fieldMetaDataMap["Time"]);
                    }
                    fileData.clear();
                    for (auto &var : strFrameData)
                    {
                        if (fieldMetaDataMap.find(var) !=
                            fieldMetaDataMap.end())
                        {
                            fileData[var] = std::stod(fieldMetaDataMap[var]);
                        }
                    }
                    if (fileData.size() == strFrameData.size())
                    {
                        break;
                    }
                }
            }
        }
    }
    if (m_session->DefinesCmdLineArgument("set-start-time"))
    {
        time = std::stod(
            m_session->GetCmdLineArgument<std::string>("set-start-time"));
    }
    if (fileData.size() == strFrameData.size())
    {
        int NumDofm1 = m_body.vel[0].size() - 1;
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_body.vel[0][i] = fileData[strFrameData[i]];
            m_body.vel[1][i] = fileData[strFrameData[i + 6]];
            m_body.vel[2][i] = fileData[strFrameData[i + 12]];
        }
        m_body.vel[0][NumDofm1] = fileData[strFrameData[5]];
        m_body.vel[1][NumDofm1] = fileData[strFrameData[11]];
        m_body.vel[2][NumDofm1] = fileData[strFrameData[17]];
    }
    std::map<int, NekDouble> Dirs;
    UpdatePrescribed(time, Dirs);
    SetInitialConditions(Dirs);
}

void ForcingMovingReferenceFrame::SetInitialConditions(
    std::map<int, NekDouble> &Dirs)
{
    for (auto it : Dirs)
    {
        int NumDof  = m_body.vel[0].size();
        int NumDof2 = NumDof << 1;
        if (it.first < m_body.vel[0].size())
        {
            m_body.vel[1][it.first] = it.second;
        }
        else if (it.first < NumDof2)
        {
            m_body.vel[0][it.first - NumDof] = it.second;
        }
        else
        {
            m_body.vel[2][it.first - NumDof2] = it.second;
        }
    }
}

void ForcingMovingReferenceFrame::InitialiseFilter(
    const LibUtilities::SessionReaderSharedPtr &pSession,
    const Array<OneD, MultiRegions::ExpListSharedPtr> &pFields,
    const TiXmlElement *pForce)
{
    std::map<std::string, std::string> vParams;
    vParams["OutputFile"]     = ".dummyMRFForceFile";
    const TiXmlElement *param = pForce->FirstChildElement("BOUNDARY");
    ASSERTL0(param, "Body surface should be assigned");

    vParams["Boundary"]           = param->GetText();
    const TiXmlElement *pivotElmt = pForce->FirstChildElement("PIVOTPOINT");
    if (pivotElmt)
    {
        std::string pstr = pivotElmt->GetText();
        std::replace(pstr.begin(), pstr.end(), ',', ' ');
        vParams["MomentPoint"] = pstr;
    }
    m_body.aeroforceFilter = MemoryManager<FilterAeroForces>::AllocateSharedPtr(
        pSession, m_equ.lock(), vParams);

    m_body.aeroforceFilter->Initialise(pFields, 0.0);
}

void Newmark_BetaSolver::SetNewmarkBeta(NekDouble beta, NekDouble gamma,
                                        NekDouble dt, Array<OneD, NekDouble> M,
                                        Array<OneD, NekDouble> C,
                                        Array<OneD, NekDouble> K,
                                        std::set<int> DirDoFs)
{
    m_coeffs    = Array<OneD, NekDouble>(5, 0.);
    m_coeffs[0] = 1. / (gamma * dt);
    m_coeffs[1] = 1. / gamma - 1.;
    m_coeffs[2] = beta * dt / gamma;
    m_coeffs[3] = dt * (1. - beta / gamma);
    m_coeffs[4] = (0.5 - beta / gamma) * dt * dt;

    m_rows = sqrt(M.size());
    m_index.resize(m_rows, -1);
    m_motionDofs = 0;
    for (int i = 0; i < m_rows; ++i)
    {
        if (DirDoFs.find(i) == DirDoFs.end())
        {
            m_index[m_motionDofs++] = i;
        }
    }
    for (int i = 0, count = m_motionDofs; i < m_rows; ++i)
    {
        if (DirDoFs.find(i) != DirDoFs.end())
        {
            m_index[count++] = i;
        }
    }
    if (m_motionDofs)
    {
        Array<OneD, NekDouble> temp;
        m_Matrix = Array<OneD, Array<OneD, NekDouble>>(m_motionDofs);
        m_M      = Array<OneD, Array<OneD, NekDouble>>(m_motionDofs);
        m_C      = Array<OneD, Array<OneD, NekDouble>>(m_motionDofs);
        m_K      = Array<OneD, Array<OneD, NekDouble>>(m_motionDofs);
        DNekMatSharedPtr inverseMatrix =
            MemoryManager<DNekMat>::AllocateSharedPtr(m_motionDofs,
                                                      m_motionDofs, 0.0, eFULL);
        for (int i = 0; i < m_motionDofs; ++i)
        {
            m_Matrix[i] = Array<OneD, NekDouble>(m_motionDofs, 0.);
            m_M[i]      = Array<OneD, NekDouble>(m_rows, 0.);
            m_C[i]      = Array<OneD, NekDouble>(m_rows, 0.);
            m_K[i]      = Array<OneD, NekDouble>(m_rows, 0.);
            int offset  = m_index[i] * m_rows;
            for (int j = 0; j < m_rows; ++j)
            {
                int ind   = offset + m_index[j];
                m_M[i][j] = M[ind];
                m_C[i][j] = C[ind];
                m_K[i][j] = K[ind];
                NekDouble value =
                    m_coeffs[0] * M[ind] + C[ind] + m_coeffs[2] * K[ind];
                if (j < m_motionDofs)
                {
                    inverseMatrix->SetValue(i, j, value);
                }
            }
        }
        inverseMatrix->Invert();
        for (int i = 0; i < m_motionDofs; ++i)
        {
            for (int j = 0; j < m_rows; ++j)
            {
                if (j < m_motionDofs)
                {
                    m_Matrix[i][j] = inverseMatrix->GetValue(i, j);
                }
            }
        }
    }
}

void Newmark_BetaSolver::SolvePrescribed(
    Array<OneD, Array<OneD, NekDouble>> u,
    std::map<int, NekDouble> motionPrescribed)
{
    for (int i = 0; i < m_rows; ++i)
    {
        if (motionPrescribed.find(i) != motionPrescribed.end())
        {
            int i0       = i + m_rows;
            int i2       = i0 + m_rows;
            NekDouble bm = 0., bk = 0.;
            if (motionPrescribed.find(i2) == motionPrescribed.end())
            {
                bm = m_coeffs[0] * u[1][i] + m_coeffs[1] * u[2][i];
            }
            if (motionPrescribed.find(i0) == motionPrescribed.end())
            {
                bk = u[0][i] + m_coeffs[3] * u[1][i] + m_coeffs[4] * u[2][i];
            }

            u[1][i] = motionPrescribed[i];
            if (motionPrescribed.find(i2) == motionPrescribed.end())
            {
                u[2][i] = m_coeffs[0] * u[1][i] - bm;
            }
            else
            {
                u[2][i] = motionPrescribed[i2];
            }
            if (motionPrescribed.find(i0) == motionPrescribed.end())
            {
                u[0][i] = m_coeffs[2] * u[1][i] + bk;
            }
            else
            {
                u[0][i] = motionPrescribed[i0];
            }
        }
    }
}

void Newmark_BetaSolver::Solve(Array<OneD, Array<OneD, NekDouble>> u,
                               Array<OneD, NekDouble> force,
                               std::map<int, NekDouble> motionPrescribed)
{
    SolvePrescribed(u, motionPrescribed);
    SolveFree(u, force);
}

void Newmark_BetaSolver::SolveFree(Array<OneD, Array<OneD, NekDouble>> u,
                                   Array<OneD, NekDouble> force)
{
    if (m_motionDofs)
    {
        Array<OneD, NekDouble> bm(m_motionDofs, 0.);
        Array<OneD, NekDouble> bk(m_motionDofs, 0.);
        for (int j = 0; j < m_motionDofs; ++j)
        {
            int j1 = m_index[j];
            bm[j]  = m_coeffs[0] * u[1][j1] + m_coeffs[1] * u[2][j1];
            bk[j]  = u[0][j1] + m_coeffs[3] * u[1][j1] + m_coeffs[4] * u[2][j1];
        }
        Array<OneD, NekDouble> rhs(m_motionDofs, 0.);
        for (int i = 0; i < m_motionDofs; ++i)
        {
            rhs[i] = force[m_index[i]];
            for (int j = 0; j < m_motionDofs; ++j)
            {
                rhs[i] += m_M[i][j] * bm[j] - m_K[i][j] * bk[j];
            }
            for (int j = m_motionDofs; j < m_rows; ++j)
            {
                int j1 = m_index[j];
                rhs[i] -= m_M[i][j] * u[2][j1] + m_C[i][j] * u[1][j1] +
                          m_K[i][j] * u[0][j1];
            }
        }
        for (int j = 0; j < m_motionDofs; ++j)
        {
            int j1   = m_index[j];
            u[1][j1] = Vmath::Dot(m_motionDofs, m_Matrix[j], 1, rhs, 1);
            u[0][j1] = m_coeffs[2] * u[1][j1] + bk[j];
            u[2][j1] = m_coeffs[0] * u[1][j1] - bm[j];
        }
    }
}

FrameTransform::FrameTransform()
{
    m_matrix = Array<OneD, NekDouble>(2, 0.);
}

void FrameTransform::SetAngle(const Array<OneD, NekDouble> theta)
{
    m_matrix[0] = cos(theta[2]);
    m_matrix[1] = sin(theta[2]);
}

void FrameTransform::BodyToInerital(const int dim,
                                    const Array<OneD, NekDouble> &body,
                                    Array<OneD, NekDouble> &inertial)
{
    NekDouble xi = body[0], eta = body[1];
    inertial[0] = xi * m_matrix[0] - eta * m_matrix[1];
    inertial[1] = eta * m_matrix[0] + xi * m_matrix[1];
    if (body != inertial && dim > 2)
    {
        for (int i = 2; i < dim; ++i)
        {
            inertial[i] = body[i];
        }
    }
}

void FrameTransform::IneritalToBody(const int dim,
                                    const Array<OneD, NekDouble> &inertial,
                                    Array<OneD, NekDouble> &body)
{
    NekDouble x = inertial[0], y = inertial[1];
    body[0] = x * m_matrix[0] + y * m_matrix[1];
    body[1] = y * m_matrix[0] - x * m_matrix[1];
    if (body != inertial && dim > 2)
    {
        for (int i = 2; i < dim; ++i)
        {
            body[i] = inertial[i];
        }
    }
}

} // namespace Nektar::SolverUtils
