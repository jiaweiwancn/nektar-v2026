///////////////////////////////////////////////////////////////////////////////
//
// File: FilterInterfaces.hpp
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
// Description: Interface class for solvers that support fluid physics
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_FILTERS_FILTERINTERFACES_HPP
#define NEKTAR_SOLVERUTILS_FILTERS_FILTERINTERFACES_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <SolverUtils/SolverUtilsDeclspec.h>

namespace Nektar::SolverUtils
{

class FluidInterface
{
public:
    virtual ~FluidInterface() = default;

    /// Extract array with velocity from physfield
    SOLVER_UTILS_EXPORT void GetVelocity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &velocity);

    SOLVER_UTILS_EXPORT bool HasConstantDensity();

    /// Extract array with density from physfield
    SOLVER_UTILS_EXPORT void GetDensity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &density);

    /// Extract array with pressure from physfield
    SOLVER_UTILS_EXPORT void GetPressure(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &pressure);

    // gave access and set to the moving frame velocity
    // for Moving reference frame formulation
    SOLVER_UTILS_EXPORT void SetMovingFrameVelocities(
        const Array<OneD, NekDouble> &vFrameVels, const int step);

    SOLVER_UTILS_EXPORT bool GetMovingFrameVelocities(
        Array<OneD, NekDouble> &vFrameVels, const int step);

    // gave access and set the displacement and angles between moving frame and
    // stationary one
    SOLVER_UTILS_EXPORT void SetMovingFrameDisp(
        const Array<OneD, NekDouble> &vFrameDisp, const int step);

    SOLVER_UTILS_EXPORT void SetMovingFramePivot(
        const Array<OneD, NekDouble> &vFramePivot);

    SOLVER_UTILS_EXPORT bool GetMovingFrameDisp(
        Array<OneD, NekDouble> &vFrameDisp, const int step);

    /// Set aerodynamic force and moment
    SOLVER_UTILS_EXPORT void SetAeroForce(Array<OneD, NekDouble> forces);

    /// Get aerodynamic force and moment
    SOLVER_UTILS_EXPORT void GetAeroForce(Array<OneD, NekDouble> forces);

protected:
    SOLVER_UTILS_EXPORT virtual void v_GetVelocity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, Array<OneD, NekDouble>> &velocity)      = 0;
    SOLVER_UTILS_EXPORT virtual bool v_HasConstantDensity() = 0;
    SOLVER_UTILS_EXPORT virtual void v_GetDensity(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &density) = 0;
    SOLVER_UTILS_EXPORT virtual void v_GetPressure(
        const Array<OneD, const Array<OneD, NekDouble>> &physfield,
        Array<OneD, NekDouble> &pressure) = 0;
    SOLVER_UTILS_EXPORT virtual void v_SetMovingFrameVelocities(
        [[maybe_unused]] const Array<OneD, NekDouble> &vFrameVels,
        [[maybe_unused]] const int step)
    {
    }
    SOLVER_UTILS_EXPORT virtual bool v_GetMovingFrameVelocities(
        [[maybe_unused]] Array<OneD, NekDouble> &vFrameVels,
        [[maybe_unused]] const int step)
    {
        return false;
    }
    SOLVER_UTILS_EXPORT virtual void v_SetMovingFrameDisp(
        [[maybe_unused]] const Array<OneD, NekDouble> &vFrameDisp,
        [[maybe_unused]] const int step)
    {
    }
    SOLVER_UTILS_EXPORT virtual void v_SetMovingFramePivot(
        [[maybe_unused]] const Array<OneD, NekDouble> &vFramePivot)
    {
    }
    SOLVER_UTILS_EXPORT virtual bool v_GetMovingFrameDisp(
        [[maybe_unused]] Array<OneD, NekDouble> &vFrameDisp,
        [[maybe_unused]] const int step)
    {
        return false;
    }

    SOLVER_UTILS_EXPORT virtual void v_SetAeroForce(
        [[maybe_unused]] Array<OneD, NekDouble> forces)
    {
    }

    SOLVER_UTILS_EXPORT virtual void v_GetAeroForce(
        [[maybe_unused]] Array<OneD, NekDouble> forces)
    {
    }
};

/**
 *
 */
inline void FluidInterface::GetVelocity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, Array<OneD, NekDouble>> &velocity)
{
    v_GetVelocity(physfield, velocity);
}

/**
 *
 */
inline bool FluidInterface::HasConstantDensity()
{
    return v_HasConstantDensity();
}

/**
 *
 */
inline void FluidInterface::GetDensity(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &density)
{
    v_GetDensity(physfield, density);
}

/**
 *
 */
inline void FluidInterface::GetPressure(
    const Array<OneD, const Array<OneD, NekDouble>> &physfield,
    Array<OneD, NekDouble> &pressure)
{
    v_GetPressure(physfield, pressure);
}

/**
 *
 */
inline void FluidInterface::SetMovingFrameVelocities(
    const Array<OneD, NekDouble> &vFrameVels, const int step)
{
    v_SetMovingFrameVelocities(vFrameVels, step);
}

/**
 *
 */
inline bool FluidInterface::GetMovingFrameVelocities(
    Array<OneD, NekDouble> &vFrameVels, const int step)
{
    return v_GetMovingFrameVelocities(vFrameVels, step);
}

/**
 *
 */
inline void FluidInterface::SetMovingFrameDisp(
    const Array<OneD, NekDouble> &vFrameDisp, const int step)
{
    v_SetMovingFrameDisp(vFrameDisp, step);
}

/**
 *
 */
inline void FluidInterface::SetMovingFramePivot(
    const Array<OneD, NekDouble> &vFramePivot)
{
    v_SetMovingFramePivot(vFramePivot);
}

inline bool FluidInterface::GetMovingFrameDisp(
    Array<OneD, NekDouble> &vFrameDisp, const int step)
{
    return v_GetMovingFrameDisp(vFrameDisp, step);
}

inline void FluidInterface::SetAeroForce(Array<OneD, NekDouble> forces)
{
    v_SetAeroForce(forces);
}

inline void FluidInterface::GetAeroForce(Array<OneD, NekDouble> forces)
{
    v_GetAeroForce(forces);
}

} // namespace Nektar::SolverUtils

#endif
