# Initial auxiliary J0 action diagnosis

The retained initial Jacobian fails the unchanged absolute, dimensionless
physical-norm threshold `1e-7` against a stable independent analytic reference:
the maximum error is `7.15996593112029e-7`. This localizes an initial numerical
action-accuracy failure beyond the original Jacobian-versus-finite-difference
comparison. It does not establish a cause of nonlinear convergence behavior,
a defect in the constitutive law, or a rigorous floating-point error theorem.
The compact [initial-action-result.json](initial-action-result.json) retains
the separate original and component gate/error fields and their source hashes;
large vectors and matrices remain on Linux.

The scope is the 2,165-unknown Linear-inner-law initial state, the first
auxiliary coefficient `alpha=3200`, and physical stage time `0.0025`. The state
binding uses identical source-bound preflight/stage geometry, expansions and
physical inputs, the producing initialization code, and the iteration-zero
native log. A completed nonlinear-history record was not used. No solved path,
target-stage status, physical trajectory or stability is qualified here.

The independent reference separates affine mechanical forcing, affine inner
scalar transport and exact time mass from the indexed outer constitutive and
scalar actions. The
roughness contribution to a momentum residual has the positive sign

\[
 \delta R_{\rm rough}[\phi]
 =\frac23\int_{r\le h}\phi\,\mathbf P\nabla\delta k_H\,dV.
\]

Its pressure counterpart is
`-∫(r≤h) phi P grad(delta pi_o)`, alongside the volume pressure term
`-∫ delta pi grad(phi)` and existing mean/gauge block. The initial exact
pressure/gauge action errors are at most `1.8470051735284574e-14`; the separate
weak-continuity/gauge column error is `1.1828409075072567e-11`. Native
nonpressure columns remain finite differences. These initial pressure checks
do not replace the separately retained pressure-active fixture.

The actual outer-k direction has zero overlap with the physical K-sampler
trace degrees of freedom and zero sampler value/tangential-gradient action.
The initial limiter plateau margin is `4.0`, and the state-gate squared
distance to activation is `0.0004` at the source/trace anchors. These checks
establish zero initial directional derivatives of chi, normalization lambda
and the projected epsilon-H target for this direction. No such assumption is
carried to later iterates.

The outer stress uses the existing indexed oracle with complex k, symbolic
cancellation of the isotropic `2k/3` mean-stress terms, and zero derivative of
the prescribed Gamma increment at fixed velocity gradient. Both formal
strain/rotation branches and complex steps `1e-20` and `1e-30` agree within
`8.393435159282324e-15`. At fixed velocity gradient and epsilon, a strain/rotation
tie does not create a k-direction ambiguity. The minimum initial Yap-switch
distance is `0.004339630861109223`.

| Comparison against the analytic reference | Absolute physical error |
| --- | ---: |
| Stored J0 action, streamwise momentum | `7.15996593112029e-7` |
| Stored J0, outer velocity load-support mask | `7.06580251680021e-7` |
| Stored J0, roughness velocity load-support mask | `1.157274987447134e-7` |
| Stored J0, k | `6.410811162823151e-8` |
| Stored J0, outer epsilon | `4.859374905274791e-8` |
| Recorded oracle FD, `h=2e-5` | `4.233260375824112e-8` |
| Recorded oracle FD, `h=2e-6` | `3.5412536912264136e-7` |
| Recorded oracle FD, `h=5e-7` | `1.4457706056678056e-6` |

The roughness and outer velocity masks contain disjoint coefficient load
rows. Their mass-dual norms are not orthogonal regional energies and are not
added. The decomposition separately verifies
`Jd - FD(h) = (Jd - A) - (FD(h) - A)`. It shows both a stored-J action mismatch
and increasing subtraction error in the finer whole-oracle differences. The
original failed gate remains false; small action-relative errors do not
replace the absolute threshold. No native intervention or retry was requested.

The producing source is `9f20b3053c303e8805f5183f668a8ec468d9656b`; solver
SHA-256 is `a90a196fe554d591fa39c6556166bf58ad43afa787b6b8e8d3d5cbf3355a5668`.
Retained J0 SHA-256 is
`469e86042f7c0d9e0291b0fa7f16011e3c7cfab53a714246778358a804128dee`,
and decoded initial h0 SHA-256 is
`f006292c518da35196f096cd36d46b062315aacc6700441cadd417d684fa69ca`.
The original action report remains unchanged, SHA-256
`52491ec3974757cf4b177b0502360171e5715448c68ab1ceb19b8bf564e38977`.

The component diagnosis used clean analysis source
`a268962a5dcc9adfe59fb412981b57242d671654`, driver SHA-256
`306c8bf37fce93beeb3217be6448ee9733bbc4f9786951b8d82f3b2ce01b9797`,
CPU 48 and all six numerical thread settings equal to one. It completed 63
checks and issued zero native calls. Authoritative Linux files remain under
`/media/jay/新加卷/nektar-v2026-runs/cls-epawf-revised/m8-r2/epoch-002/initial-k-diagnosis-001/`:

- `validation.json`: 161,857 bytes, SHA-256
  `2224a4e387a92d0a8a160e817bb9265534a3d7291384faf025d23ddea8b6ecf0`.
- `receipt.json`: 11,880 bytes, SHA-256
  `63c0506106cc85bdadb3f384f8478163bef516d45b780048d170ce13189b8fba`.

The next analysis helper adds a diagnostic-only signature and distance for
the original Yap switch `k^(3/2)/(2.5 epsilon r)=1` on the fixed actual outer
quadrature points. It changes no source, constant, field or admission rule.
M8 remains open.

## Binding to the actual original target

The retained first `alpha=400` target matrix independently confirms the same
`7.15996593112029e-7` momentum action error. The producing control flow and
the completed failed auxiliary record establish that the target starts from
the original h0: a failed auxiliary candidate never replaces the guess.
The log contains a unique target iteration-zero entry. This is a source and
record binding; no unpublished target vector or completed target path is
inferred.

The independent analytic reference transfers with the exact scalar/velocity
mass change, `A400=A3200-(7/8)*M3200`. The discrepancy in that mass-change
control is `1.3472529138953944e-11`, below its `1e-10` control threshold.
Pressure actions have zero target-minus-auxiliary difference; their analytic
triangle bounds remain at most `1.84701e-14`. The continuity-column error
remains `1.18284e-11`. These controls preserve the failed absolute outer-k
action gate and establish no causal explanation of the target's convergence.

The Linux check uses clean analysis source
`8d03c9a84e2dc1de118a327c1e17a2e93a8b4766`, CPU49 and all six numerical
thread settings equal to one: 72 checks, zero native calls and zero new
residual endpoints. Producing C++/runtime remain the same 9f20 epoch.
Target J0 SHA-256 is
`ae6e804b43d90be64ce77c3e44dd532970d35f7a37c449984cb28e6bfaf127e6`.
The authoritative files remain in `epoch-002/target-initial-jacobian-001`:

- `validation.json`: 158,084 bytes, SHA-256
  `29396c6623849f0257869e77624eafb6cae6637ab9c36f0709c1e90b98522d76`.
- `receipt.json`: 37,490 bytes, SHA-256
  `a33961cfc448cb8cd99059644a3ffbcf8d2550a5129bfd0e082936ddec236b4a`.
