// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_ADMISSIBILITY_H
#define NEKTAR_WJ_NATIVE_ADMISSIBILITY_H
#include <IncNavierStokesSolver/EquationSystems/WJNativeWallV2.h>

namespace Nektar::WJNative
{
struct PlanarAdmissibility
{
    double covarianceOverK=0, principalOverMu=0, principalAngle=0, stationaryResidual=0;
    double principalRawOverMu=0, principalReserveOverMu=0;
    int stationaryCandidates=0;
    bool principalEvaluated=false;
};

/** Current flat x-z state, n_wall=e_z; inactive spanwise rows/columns must
 * vanish. The base gradient need not have pointwise zero divergence. The
 * acoustic perturbations alone are incompressible planar rank-one matrices.
 *
 * Full 3D covariance is checked. The principal result minimizes the PLANAR
 * angle continuum with frozen k,E,aux,r,H; it does not certify all 3D modes.
 * For chi>0 principalOverMu subtracts the reported arithmetic reserve from
 * the raw angular minimum. Negative scientific margins are returned normally.
 * principal=false skips the acoustic calculation and sets principalEvaluated
 * false; its zero principal fields must not be interpreted as a pass.
 */
PlanarAdmissibility CheckPlanarState(double k,double epsilon,const Tensor &gradient,
    double nu,double auxiliary,double r,double H,bool principal=true,
    const WallStress*currentValue=nullptr,
    ConstitutiveBackend backend=ConstitutiveBackend::Historical,WallEvaluationStats*stats=nullptr);
}
#endif
