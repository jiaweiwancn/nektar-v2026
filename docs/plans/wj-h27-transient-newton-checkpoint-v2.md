# H27 unchanged-v2 transient Newton checkpoint

The registered solver now has an opt-in `WJImplicitSolver=CoupledJFNKV1`
advancement path. It solves the same complete current BE/BDF2-EX2 equations
as the default PicardV1 control. The approved v2 tensor, scalar law, cutoff,
reference, physical geometry/inlet/boundaries and constants remain active.
The separately proposed v4 scalar/reference amendment is not implemented or
approved. This is a numerical checkpoint, not positive H27 closure.

## Source, runtime and reproduction

Tested numerical identity:
`1c0091bff8ead2247d11602203de3a742a49e1cad3ca8ef2b349312b8f2ca5d9`.
The immutable Linux source-v12 capsule is source-v11 plus 123 frozen overlay
files from Mac parent `45e21cac27fb763297cea56b1cafeaf7e630ff16` and the
uncommitted numerical changes. Current C++/headers match the tested capsule.
This containing commit binds those tested contents; the build preceded it.

The full normal-main executable was built on authoritative Linux in 203.561 s,
with four compiler workers. Its ELF SHA256 is
`2cafcb83456168871d85eddfc745affaed11588ce9708bf805fe1b2e313d8ff7`.
The evidence directory
[`wj-h27-linux/native-time-newton-v2d`](wj-h27-linux/native-time-newton-v2d/inventory.json)
contains exact source/runtime/shared-library identities, compile/link commands,
every new receipt and failure log, final XML/opt inputs, complete small MMS
states/checkpoints/retained outputs, independent reports and compact large
physical results. Complete raw fields remain under
`/home/jay/build/wj-h27-coupling-r1-20261117`, bound by each summary's path/hash.
The evidence archive is 4,264,571 bytes, SHA256
`404a899da15f47c6105605c2c04e770ef5a48679187de92d956e2ca95283b8e9`.

All numerical jobs used one MPI rank, with OMP/OpenBLAS/MKL/BLIS/VECLIB/NumExpr
threads one and the <=72 actual shared-host watchdog. No installed library or
unrelated checkout was changed. Mac work was source/Git/hash/AST/review only.
Build, native, oracle and restart commands are preserved verbatim in the
operation requests. A representative native command is:

```
/home/jay/build/wj-h27-coupling-r1-20261117/runtime-v12/IncNavierStokesSolver \
  --use-opt-file /home/jay/build/wj-h27-coupling-r1-20261117/science/coupled-fixtures-v12/physical-new.xml.opt \
  /home/jay/build/wj-h27-coupling-r1-20261117/science/coupled-fixtures-v12/physical-new.xml
```

That is the recorded completed command. Replay uses a fresh owned prefix;
the create-only outputs intentionally reject overwriting existing evidence.

## Algorithm and failure semantics

Newton uses exact physical BDF coefficients and complete accepted histories,
with conservative explicit EX advection fixed throughout each implicit solve.
The current full stress/scalar/trace residual is never replaced by its inverse
approximation. The latter uses the full frozen-scalar velocity tangent, mass,
positive scalar diffusion/reaction, transported-Q3-k action and dependent Gamma
derivative. It contains no current transport Jacobian for explicit advection.

Uncommitted trials retain typed state guards. Every physical step recomputes
full planar principal admissibility and the complete current residual before a
single history commit. Only typed linear failure permits the one same-dt cold
retry. One cooperative step deadline includes setup, both possible solves and
final checks. Failed accepted histories and disjoint work are preserved;
interrupted engine counts are explicitly marked incomplete. Inclusive Newton
trace observations are labelled as uncommitted, with FD probes excluded.

Three native readers received a test-only exact `-0` parsing repair outside the
immutable source capsule; no equation or threshold changed. The new history
checker is likewise a separately hashed validation overlay. The first physical
oracle copy lacked the retained Fig.3 CSV and failed before profile assessment.
Its failed report/40.853 s receipt remain intact. A new validation namespace
restores the exact original fixture; the repeated 51.888 s check meets the
accepted-step equation/history/quadrature subgates. Its top-level NOT_MET
correctly retains the separate physical-reference/full-H27 gates.

## Independent equations and algorithm comparisons

Both PolynomialScalarsV1/body=1 (four steps, dt=.02) and Historical MMS (eight
steps, dt=.0025) pass independent current equations at refined integration
orders. All paired Picard/Newton full states, original-owner fields, dependent
traces/cutoffs and complete explicit histories pass the fixed <=1e-8 gate.
The physical pair uses three actual unforced startup steps at dt=.005, the
same retained 700 cells, tolerances 1e-12 linear/3e-12 coupled, and no Anderson
acceleration in either algorithm. Pressure levels are retained.

For the physical endpoint, accepted coefficient MAX differences in the fixed
scales are U 1.834e-12, Pi 5.291e-9, k 1.470e-12 and epsilon 5.179e-13.
Independent current momentum/continuity/k/epsilon residuals are
9.366e-15 / 6.524e-16 / 3.405e-13 / 1.193e-14. Refined quadrature and
accepted-history boundary checks meet their declared subgates.

| Paired calculation | Picard whole operation s | Newton whole operation s |
| --- | ---: | ---: |
| Polynomial MMS, 4 steps | 2.782 | 2.770 |
| Historical MMS, 8 steps | 4.173 | 5.539 |
| Physical startup, 3 steps | 143.723 | 39.234 |

| Physical step | Picard complete step s | Newton complete step s |
| --- | ---: | ---: |
| BE, t=.005 | 38.654 | 11.534 |
| BDF2/EX2, t=.010 | 55.669 | 10.977 |
| BDF2/EX2, t=.015 | 43.540 | 10.988 |

The measured three-step means are 45.9545/11.1666 s, a 4.11535 ratio. These
are initial advancing prefixes, not three warmed blocks, a selected accurate
dt or a complete production forecast. The final Newton step has four accepted
Newton updates and seven coupled Krylov iterations. Its 10.9844 s native work
contains 2.2357 s constitutive, 1.7773 s sampling/assembly, 4.2660 s rebuild,
1.8515 s guards and 0.5112 s residual work; remaining disjoint categories are
reported exactly. Advancement attribution does not close whole-trajectory cost.

The Newton MMS foundation reproduces current/previous/older fields and explicit
histories exactly after cold continuation, injected retry and continuation from
retained recovery step 3. Four scheduled writes actually remove recovery 1/2,
retain 3/4, and retain scientific outputs 1..4. One-Newton-iteration and tiny
deadline failures preserve the complete step-zero accepted state. A Picard
restart of a Newton checkpoint rejects its scientific identity. These controls
pass independently; physical cold restart, terminal output and full H27-D remain
open.

## Physical reference and perturbation foundations

The additional unchanged-source-v11 steady reference doubles x resolution from
2800 to 5600 cells while keeping the 100 radial intervals. It converges in
708.726 s, five Newton updates and 115 Krylov iterations, native norm 2.400e-13;
peak RSS is 6,621,138,944 bytes. Independent equations meet their thresholds.
Inner-k retained-inlet RMS/MAX departures are .132643/.423656 against .05/.10.
The preceding 700/2800 values are .132331/.459519 and .132795/.440730.
The 2800/5600 inner-k differences .008029/.130882 still exceed .001/.001,
so neither a resolved reference nor a conclusive continuum barrier is claimed.

The independently checked source-v11 perturbation increments pass both fixed
shapes/signs, global peak/support/RMS/half-area and native projection checks.
The primary shape also passes on 2800 cells using its original certificate,
knots and normalization. This qualifies increment construction; no seeded,
admissible coupled state or evolved perturbation pilot is claimed.

## Closure and continuation

The completed OPEN account has 8715.746188972902 aggregate operation seconds,
290 entries and no unresolved operation at publication. It includes the prior
closed 703.6786409680499 s and 120 s bootstrap/control provision, all failed
experiments and expected negative controls. The same 43200 s ceiling and 3600 s
protected reserve remain; no allowance resets and the account is not closed.

All overall gates remain unclosed: model/reference and physical accuracy fail;
full coupled-numerical closure must bind the final selected model/recipe;
physical pilots, terminal retention/monitoring/storage, warmed complete cost and
the executable Task3 manifest remain outstanding. `TASK3_READY=NO`.
T=505.02466521278404, terminal [447.8483832973273,505.02466521278404] and
per-trajectory 1.25*Ctotal<=86400 s remain unchanged. No full Task3 run occurred.

The [audited v4 proposal](wj-h27-historical-inner-transport-proposal-v4.md)
is available for the separate scalar/reference model decision. It supplies
exact equations and a preliminary independent column control, including the
mandatory nonzero common-flux MMS after a diffusion jump. Approval is pending;
active v2 physics has not been changed to v4.
