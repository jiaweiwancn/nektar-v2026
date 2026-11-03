# R2 pressure-action localization

The replacement pressure action passes the unchanged absolute, dimensionless
physical mass-dual/block/constraint threshold of `1e-7`. This conclusion is
supported by an independent exact affine assembly and independent analytic
sensor/profile/capped-source/projected-trace differentiation. It does not claim
that every finite-difference step passes that absolute threshold.

| Retained test | Independent exact-action error | Scope |
| --- | ---: | --- |
| Quiet native fixture, all four pressure/gauge directions | maximum `3.4191644960823145e-15` | Linear inner law; affine pressure/gauge action |
| Pressure-active native fixture, all four directions | maximum `9.787439875072821e-11` | SaturatedR2 inner law; interior pressure-driven transition, cap attenuation and nonzero projected epsilon-H action |
| Actual R1-B2 failing state and pressure/gauge direction | `3.4287813857863435e-15` | Archived coefficients used only for operator localization |

The quiet and active independent checks counted 260 and 360 assertions,
respectively. The active outer-pressure nonlinear source/trace action has
physical norm `1487.239431049782`, including a constraint norm
`48.548394581796074`; the action is therefore measurably active rather than an
affine surrogate. Its independent constraint-action error is
`6.120027527698541e-12`.

At the actual R1-B2 state, the independently reconstructed nonlinear pressure
source/trace action is zero. This was evaluated from the actual state and gates,
not assumed from a global chi setting. Relative to the independent exact
action, the retained cached and fresh R1 Jacobian actions have absolute errors
`2.4528642896580794e-7` and `1.5701508091699785e-6`. Both remain failures of
`1e-7`. These are errors against the exact action, distinct from the historical
`2.58e-6` fresh-Jacobian versus fine-difference discrepancy.

At that same state, the full-residual central-difference error worsens from
`9.741561763681852e-9` at dimensionless `h=2e-4`, through
`1.6276407242677775e-7` at `h=1e-5`, to `3.2441193813552027e-6` at `h=5e-7`.
The pressure contribution is affine on this actual inactive nonlinear branch,
and the independently assembled affine action agrees with the new native
action to roundoff. The retained ladder consequently localizes subtraction
error in the old whole-residual differencing; making that difference finer is
not evidence of a more accurate action. The new assembly avoids subtracting
the unchanged residual components while retaining the nonlinear chain where
it is active. This operator diagnosis does not establish that a changed
Jacobian policy cures the nonlinear stage failure.

In the pressure-active fixture, central differencing also has finite-step
nonlinear truncation error. At `h=5e-7`, the outer-pressure difference has
absolute error `2.9213449414673976e-6` and action-normalized relative error
`1.9642734589180027e-9`. The joint-pressure difference has absolute error
`6.326713562861098e-6`. These absolute errors fail `1e-7`; their small relative
errors are reported separately and do not replace the fixed absolute gate.
The exact analytic chain is the evidence that passes the gate.

Both fixtures are operator-verification cases, not physically qualified
trajectories or fixed-outer-box coverage claims. No native time trajectory or
nonlinear retry was run for this analysis.

## Producing and analysis identities

Native outputs were produced on Linux from clean source
`1ee3a8454c857c807bcb6911024186afb9d22721`. The probe executable SHA-256 is
`a818a70969de2cb7446d5b7f38b36ae7ee3b31e42b6798fee9673785a20714f2`.
The original validation attempt failed only at JSON serialization of an
explicit NumPy boolean; its reports and failed receipt were preserved.

The analysis-only replay used clean driver source
`10660ea1b2b1fd8c7f2884c0fe6aff4da2be84fd` and driver SHA-256
`f169abac5abb7af02adaa727c86e6ec4d91685fd6e951070fac955380c0231cd`.
It ran on CPU 48, one process, all six numerical thread settings equal to one.
The replay issued zero native commands, verified the original before/after
runtime pins, checked the probe image against those pins, and required the
exact native-output hashes recorded before the replay driver was written.

The canonical results remain at
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r2/pressure-reanalysis-002/validation.json`:
176,034 bytes, SHA-256
`27303ba92b3deb57b8f298cb2c074e9aa203ea710d7028dfe248dbf7c3b73021`.
That receipt retains every ladder result, imported-source hashes and the
historical runtime section hashes. The native reports remain in
`m8-r2/pressure-verification-001/results/`:

- `affine/probe.json`: 778,665 bytes, SHA-256
  `f14726aeca94ee69cf9616a5c09dc99fca87c989b72732dc31a5ff3e40f62760`.
- `pressure-active/probe.json`: 3,043,737 bytes, SHA-256
  `cf3a8844dff263b41d3ed39793c2300acd7d3b668610168fa4f26b8a450d33af`.

The R1-B2 nonlinear trace is unchanged, SHA-256
`2df975858f4d506b62030cdf4ea899fff5dd31740ea54289c166a10b156c4475`.
Its failing state was never used as a physical warm start.

The pressure suite was subsequently repeated at the actual producing solver
epoch `9f20b3053c303e8805f5183f668a8ec468d9656b`, with the same 260/360
independent fixture checks and the unchanged `1e-7` physical action gate.
That current-runtime result passes, with unchanged before/after fingerprints;
it is distinct from the historical analysis-only replay above. Its report is
`m8-r2/epoch-002/pressure-verification-001/results/validation.json`,
174,042 bytes, SHA256
`eb8e184f496c22bbec5713af6f82eeb4437011e20d6afdf02d9026e541d8c50a`.
The [inner evidence ledger](inner-verification-result.json) also binds this
separate pressure result. It does not certify nonpressure Jacobian columns.
