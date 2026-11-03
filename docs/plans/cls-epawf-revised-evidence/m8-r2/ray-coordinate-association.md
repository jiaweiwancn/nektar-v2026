# Moving-ray coordinate association

The original full auxiliary-path audit remains **failed**. It finished after
9,683.398 seconds with five coordinate-identity exceptions: iteration bases
3, 5, 6 and 7, and the terminal state. No full-path pass is inferred from
the focused follow-up below, and that audit has not been rerun.

Its authoritative files remain under
`epoch-002/auxiliary-analysis-linear-002` on Linux:

- `validation.json`: 19,119,016 bytes, SHA256
  `9208cf90b8d6bd7f11f0f46d117c5e9b9d79cac640f8518fecf49855465d0aae`.
- `receipt.json`: SHA256
  `0e114d5fe6e9978bd5ed939226204ad9c94a0c96b3e3efe2678f639e98abf116`.

A single, bounded terminal-domain replay localizes that terminal exception.
One of 89,490 keys differs because these adjacent binary64 radii straddle
the decimal midpoint `0.1915526386605`:

| Quantity | Native | Independent |
| --- | --- | --- |
| Exact binary64 representation | `0x1.884cbff413279p-3` | `0x1.884cbff413278p-3` |
| Radius printed in decimal | `.191552638660500013...` | `.191552638660499985...` |
| Key rounded to 12 places | `.191552638661` | `.191552638660` |

Their separation is one ULP, `2.7755575615628914e-17`. The single unmatched
pair is unique within its exact six structural descriptors. All sorted rows
retain the same structural identities; their maximum raw coordinate
difference is `4.440892098500626e-16`, and their maximum weight difference
is `2.7755575615628914e-17`. The existing coordinate and weight accuracy
thresholds are both `2e-11`. Native and independent scalar minima, lambda,
b and cutoffs agree. No state or evaluation coordinate is replaced.

This is an independent-analyzer lookup defect, not an observed defect of
the native trial locations or a cause of nonlinear stagnation. The diagnosis
used clean source `42863271b58e34aef7b97d7240de6fdccecfb10c`, CPU49 and
six numerical thread settings equal to one: one domain reconstruction,
66.378 seconds, zero weak-residual, Jacobian, endpoint or native calls.
Its original failed identity gate remains false. Exact Linux bindings are
in `epoch-002/ray-identity-diagnosis-001`:

- `validation.json`: SHA256
  `1a590ae86fe3b33c0ad65b19f40d7055cfd4aa36d7dbe02e2d110ef055f137c8`.
- `receipt.json`: SHA256
  `e1577b4dbac6f4a59611b941af1162591b02a30937c8238f9b6d39d9a818f13d`.

Correspondence v2 keeps the original sorting, exact row count and exact six
structural fields. It compares the actual coordinates and weights at the
unchanged `2e-11` thresholds. If rounded keys differ, **every paired raw
coordinate additionally must agree within `1e-13`**. It records the old
rounded-key equality separately. Eight focused Linux fixtures pass,
including wrong labels, counts, coordinate separations and weights that
must reject. Their receipt is
`epoch-002/ray-correspondence-tests-001/receipt.json`, SHA256
`2850769ab07229518af06a1109a7f4ea71ca3a29e6abab12dea6ad312a80b0e0`.

Only the independent association code changes. The producing C++/binary,
constitutive equations, physical box, scalar evaluations, nonlinear residual
and all numerical acceptance tolerances remain unchanged. The initial
outer-k action failures and the original full-path failure are retained.
