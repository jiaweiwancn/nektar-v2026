// SPDX-License-Identifier: MIT
#ifndef NEKTAR_RANS_WALL_NUMERICAL_OPERATOR_H
#define NEKTAR_RANS_WALL_NUMERICAL_OPERATOR_H
namespace Nektar
{
// Numerical row ownership, independent of closure/wall physical descriptors.
// No legacy evolution history may be continued with this changed operator.
inline constexpr const char *RANSWallNumericalOperatorIdentity()
{ return "GeneratedVolumeOwnerV1"; }
}
#endif
