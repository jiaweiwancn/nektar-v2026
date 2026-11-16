# H25R1 native prototype contracts

Frozen before measuring the prototypes. These are separate from the rejected
h350 representation and prescribed finite seed. Neither changes the WJ law,
root, physical input, numerical gates, IMEX1 sweep or ordinary restart identity.

## Exact assembly scatter

RANSWallAssemblyScatter=ElementExactV1, default Off, creates assembly-local
handles to the native sparse matrix entries for each element and velocity/scalar
block pair. Each original contribution retains its expression, visitation order,
zero-skip rule and binary64 +=. No matrix/factor/coefficient is reused. Handles
are discarded before exterior/constraint row erasure and never cross matrix,
free-row, geometry, owner, time, data or state generations. Extra zero map entries
are excluded by the original Each/Pattern/Nonzeros operations.

The native component compares coefficient bits, action, pattern, signed-zero/
cancellation inputs, dense/sparse sinks and row replacement after handle release.
Then matched Off two-step and On four-step advancing replays use the retained
H24 M16 start and W0 states/history/callbacks, with exact full endpoint native audit
parity. Timing trees retain the original disjoint taxonomy. No benefit is assumed.

## Sufficient-positivity prototype with exact fallback

The offline Bernstein heuristic covers all 700 k cells in the saved state, but
does not prove enclosure of native barycentric evaluation. The bounded native
prototype therefore uses the actual computed native basis values, avoiding
an unsupported analytic-to-native rounding claim. It is read-only: every exact
point evaluation and every production guard still runs, including full-audit
observations and original adverse first failures.

Let c be the native constant centre and delta_j the exact binary64 coefficient
differences actually used by RANSWallGlobalPoint::Value. Compute M=max|delta_j|
once per element for the current const effective lifted state. For the actual
computed double basis values bx_i and bz_j, let Sx and Sz be outward upper sums
of their absolute values. A cubic 2D row has 16 entries and one rounded product
per signed tensor basis value. With u=2^-53 and d=denorm_min(binary64), use

    P = up(up(Sx*Sz)*(1+2u) + 16d)
    D = up(M*P)
    A = up(abs(c)+D)
    gamma = up(35u/(1-35u))
    E = up(up(gamma*A)+64d)
    L = down(down(c-D)-E).

Every up/down operation applies nextafter to enclose its binary64 result.
The 16 long-double products and additions plus final double conversion need
at most 33 rounding operations; gamma35 uses the larger binary64 unit roundoff
when long double has at least 53 bits. The additive 64d covers underflow. The
rounded basis product sum is bounded by P. Therefore L bounds the actual native
constant-centred value from below, under IEEE binary64 and FE_TONEAREST.
Unsupported rounding, nonfinite/overflowed bounds or inconclusive L use exact
fallback. Native coordinate validation and the original clamped basis evaluation
remain live. This is a sufficient point positivity bound, not a cell polynomial
minimum, physical/principal certificate or replacement for the current-K budget.

Summaries are local to each current domain evaluation and bind its const effective
state address/size, field coefficient map, element/owner and offset. Fresh summaries
are built for each state/data/time generation; source/effective-state hashes, map
identity and time bits are recorded. Test A, changed k, changed epsilon, negative k
and return to A with fresh generations. The exact fallback is evaluated at every
prototype point and must satisfy L<=native_value. Certified lower bounds and
actually observed minima are reported separately.

Report measured native moving-point evaluation-cost-weighted coverage and bound
cost, with the observation inventory and zero skipped production guards. This
micro-cost coverage is distinct from the enclosing production guard time and
complete-step speedup. Integrating typed certificates into correction/minimum
diagnostics is not accepted by this read-only prototype. Full guards remain the
selected policy unless a later, separately qualified integration preserves all
failure, correction, minimum and coverage semantics.
