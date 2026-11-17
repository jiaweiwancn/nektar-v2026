// SPDX-License-Identifier: MIT
#ifndef NEKTAR_WJ_NATIVE_FAILURES_H
#define NEKTAR_WJ_NATIVE_FAILURES_H
#include <stdexcept>
#include <string>
namespace Nektar::WJNative
{
/** A recoverable factor/Krylov failure. Scientific input, positivity, root,
 * covariance and constitutive failures deliberately use other exceptions. */
class NativeLinearFailure:public std::runtime_error
{public:using std::runtime_error::runtime_error;};
/** A provisional current state violates a declared scientific margin. Time
 * advancement treats this as fatal without altering accepted histories. A
 * steady solver may explicitly classify it as an infeasible Newton trial. */
class NativeInadmissibleState:public std::runtime_error
{public:using std::runtime_error::runtime_error;};
inline void RequireState(bool condition,const std::string&message)
{if(!condition)throw NativeInadmissibleState(message);}
inline void RequireLinear(bool condition,const char*message)
{if(!condition)throw NativeLinearFailure(message);}
}
#endif
