# H27 CP3 short startup and complete-process cost screen

**H27_COMPLETE_COST=TARGET_NOT_MET.** The bounded source19 comparison and
startup control matrix finished. The measured startup settings give no
accuracy-qualified way to meet the complete per-trajectory budget. CP2
classification C also leaves the developed initial state unqualified.
`TASK3_READY=NO`.

The Linux experiments use the source19 executable, table and retained-inlet
700-cell 14×50 flat mesh, with Q3 velocity, k and outer epsilon, Q2 pressure
and coupled Gauss 8×12 per governing panel. The node rule for native field
construction remains GLL8×8. The coupled model is the approved finite-ν
v2 law with actual stress production, moving cutoff and complete current
scalar diffusion. CP1's coupled numerical foundation passes on this source.
[Numerical qualification](wj-h27-continuation-cp1.md)
[Reference decision](wj-h27-continuation-cp2.md)

The short startup uses BE and BDF2/EX2, the original applied boundary laws,
current-data inverse reuse and genuine accepted history. No developed or
finite-perturbation case is admitted while the physical reference gate is
open. The reference preparation costs already consumed in H27 are separate
from recurring launch work. [Frozen short-screen protocol](wj-h27-cp3-short-screen-admission-v1.md)

The two earlier conventional data sets keep their original scopes. H26
used GLL8×8, prescribed analytic forcing and matched boundary data within
its kernel comparison. H25R1's S0 measured the conventional no-RANS path on
the same cell count under its own integration and output configuration.
The actual source19 paired cost run uses stock `VCSWeakPressure`, continuous
Q2 pressure and Q3 velocity, constant molecular viscosity and tensor
Gauss-Legendre 8×12 per cell. The WJ and stock discretizations share the
physical mesh, historical analytic inlet speed, time step .005, 12-step
window, same ELF, one rank and six single-thread limits. The stock pressure
Poisson and component Helmholtz scheme uses high-order pressure conditions
and outlet pressure Dirichlet data with zero natural velocity-component
loads. The WJ mixed pressure and velocity solve retains the prescribed
total traction, the wall load, and both turbulence quantities. The velocity
projection operators and output contents differ. Both write each third
step and at the end. Whole process time includes each solver's actual output
work. [H26 comparison](wj-h26-closure.md) [H25R1 comparison](wj-h25r1-closure.md)

| 12-step run at .005 | Complete process wall s | Warm seconds per step | Timer scope |
| --- | ---: | --- | --- |
| Stock 1 | 1.48 | .010383, .010632, .010490 | Child process elapsed, warm internal advances |
| Stock 2 | 1.47 | .010382, .013116, .011122 | Child process elapsed, warm internal advances |
| Stock 3 | 1.49 | .010373, .010815, .010469 | Child process elapsed, warm internal advances |
| Coupled | 79.818 | 6.288725, 5.810171, 6.269996 | Runner process elapsed, warm complete steps |

All stock repeats complete the 12 advances with one initial checkpoint,
four scheduled checkpoints and a final field file, totaling 913,091 bytes
per block. Stock uses no SVV. Its internal timer excludes checkpoint
work and some post-step work. The source19 coupled endpoint reproduces
all four accepted current, previous and explicit-history SHA256 roles
from the previously independently audited .06 startup. Each advance
uses four Newton bases and the actual coupled Krylov counts are
9, 7, 8, 7, 7, 7, 7, 7, 7, 7, 7, 7. Seven inverse refreshes produce
21 factor preparations. There are zero rejected trials and zero retries.

The coupled native engine advances total 74.188 s. Its disjoint shares
are 29.9% guards, 25.0% assembly, 23.6% constitutive and wall work, 10.1%
inverse rebuilds, 8.7% residual work, 2.2% linear solves and .4%
engine orchestration, with small other components. Complete-process costs
include 1.102 s construction, 2.420 s initialization, .527 s scientific
frames, .101 s recovery and .040 s final export. Scientific frame and
recovery costs are already included in complete advancing wall time.
Exclusive measured categories account for 98.17% of the total after
leaving the residual runner, journal and engine remainder unattributed.
The recorded model/equation/linear work categories account for 55.30% of
the whole process, with guard and inverse rebuild time separate. Those
shares apply to the paired startup. The short .08 case has only two warm
observations and measured whole-process attribution of 91.72%.

## Startup accuracy and affordability

The three retained-inlet starts at dt=.02/.04/.08 accept 12/6/3 genuine
steps respectively and all reach t=.24. The independent equations pass
on both sides of all compared endpoints. The common-time physical field
comparisons retain the original element owners and fixed scales.

| Common-time comparison | Inner-k RMS / MAX in Kin | Inner-U RMS / MAX | Frozen comparison misses |
| --- | ---: | ---: | ---: |
| .04 vs .02 | .00042800 / .00113463 | 7.354e-7 / 2.232e-6 | 62 |
| .08 vs .04 | .00163905 / .00431980 | 2.937e-6 / 6.495e-6 | 116 |
| k and U limits | .001 / .001 | 1e-4 / 1e-4 | 0 |

Increasing the inverse reuse age from 8 to 24 at fixed dt=.04,
threshold .02 and tolerance 3e-12 saves one inverse refresh. Its complete
process is 46.102 s versus 47.298 s, and its largest observed warm step is
7.365 s versus 7.387 s. Every common-time physical numerical row passes.
The result is a short optimization measurement with nearly unchanged
worst warm cost. The independent tolerance control keeps inverse age 8,
changes only the nonlinear tolerance to 1e-10 and reduces four Newton
updates to three on steps 5 and 6. It leaves the largest measured warm
step essentially unchanged. The independent tolerance check passes every
comparison row, with inner-k RMS/MAX 4.58e-12/1.48e-11. The independent
looser endpoint k equation residual is 2.085e-11 against 1e-8.

| Short recipe | Largest observed warm s/step | s per physical time unit | 1.25 × advancement-only seconds to T |
| --- | ---: | ---: | ---: |
| .005, paired | 6.28873 | 1257.745 | 793,991 |
| .02 | 6.37074 | 318.537 | 201,092 |
| .04 | 7.38700 | 184.675 | 116,585 |
| .08 | 7.42046 | 92.756 | 58,557 |

The projection multiplies the largest observed warm mean by the step
count to reach or exceed T=505.02466521278404. It is conditional on rates
from at most .24 physical time. The .08 screen fits the advancement-only
proxy and fails the unchanged accuracy gate. The .04 screen misses both
inner-k maximum accuracy and the projected day. The .02 screen exceeds
the day before fixed costs. No timestep is qualified by this screen.
Spatial and polynomial accuracy and finite-perturbation timestep agreement
remain open. Production output, retry and cold restart, terminal assessment
and monitoring costs remain to be measured on a qualified recipe.

The 86,400 s limit with 1.25 provision allows 69,120 s complete wall time,
or under 136.865 wall seconds per physical unit before any fixed costs.
Across four trajectories, the accepted-advancement-only sum at the .04
startup rate would be 373,073 s, and 187,382 s at the .08 rate, before
provision. These are conditional sums using the same observed startup rate
for illustration. Separate developed-control and perturbation rates and
qualified campaign costs are still unknown. Four required pilot endpoints
at or beyond 26.13350407707765 would consume 33,306 s at the .02 rate,
19,324 s at the .04 rate and 9,706 s at the .08 rate for advancement alone.
The .02 proxy exceeds the remaining development account. At .04 it leaves
about 503 s for the rest of the required matrix. The faster .08 proxy
does not pass the existing k and stress budgets. The required long matrix
has no reference and affordability admission.

## Operational foundations and remaining gates

Source19 already exercises full-history cold starts and actual scheduled
recovery rotation in short controls. Both mathematical perturbation
certificates bind the unchanged 700-cell geometry and perturbation source
files. Source11 evidence passes the primary, finer-primary and secondary
native increment constructions. A source19 coupled perturbation requires
a qualified developed state, current dependent trace and accepted
admissibility at the original amplitude. [Perturbation binding](wj-h27-linux/continuation-cp3-source19/perturbation-foundation-reuse.json)

The source15 accumulator, terminal checker and native scientific frame
writer match their source19 bytes. Coupled actions, reference, moving rules
and the independent physical oracle have changed. This is a source reuse
audit for the matching pieces. Selected-recipe cold persisted statistics,
monitoring escalation, terminal reduction accuracy and storage are open.
The production terminal window remains [447.8483832973273,
505.02466521278404]. [Statistics source audit](wj-h27-linux/continuation-cp3-source19/persisted-statistics-source-reuse.json)

| Global gate | Current result |
| --- | --- |
| H27_MODEL_AND_REFERENCE | NOT_MET. CP2 classification C. |
| H27_COUPLED_NUMERICS | PASS on source19. |
| H27_PHYSICAL_PILOTS_AND_ACCURACY | NOT_MET. Four principal pilots and a jointly qualified h/p/q/dt recipe remain open. |
| H27_RESTART_MONITORING_STORAGE | NOT_MET on a selected recipe. |
| H27_COMPLETE_COST | TARGET_NOT_MET. Conditional rates and actual complete short blocks are measured. |
| TASK3_READY | NO. Launch commands and production forecasts are unset. |

## Provenance and resources

The first summary operation treated the native residual vector as a mapping
and stopped without publishing a report. It is preserved and charged.
The corrected report reads the vector and the actual coupled Krylov counters
from the exact archived step logs. All native runs and independent accuracy
results are unchanged. The whole CP3 phase charged 679.074828 s including
both reporting operations and publication. The OPEN aggregate account was
19,772.376549 of 43,200 s at settlement, with 19,827.623451 s available
before the 3,600 s evidence reserve and no reservations. CP3 used no build.
The peak was 12 host workers and 638,025,728 bytes process-tree RSS. Science
and source/build/runtime footprints were 3,487,350,585 and 3,007,184,492
bytes. Free space was 68,651,671,552 bytes.

The authoritative Linux archive `evidence/continuation-cp3-source19.tar.gz`
contains 87 checked members plus the inventory, 302,493 bytes, SHA256
`0fb39e82bffc0951569a97c37445eb047806d3ee8978a2d016df1b5f4dfbbc6d`.
The archive includes exact session/options files, launch and step logs and
runner requests. Raw coefficient fields, recovery and oracle outputs stay
in the live capsule. Reproduction uses `runtime-v19/IncNavierStokesSolver
--use-opt-file <session>.xml.opt <session>.xml` through `tools/run.py` with
new operation labels and output prefixes. The corrected analysis helper
can read the original raw data without starting the solver.
[Cost assessment](wj-h27-linux/continuation-cp3-source19/assessment-reviewed.json)
[Accuracy comparison summary](wj-h27-linux/continuation-cp3-source19/comparisons.json)
[Evidence inventory](wj-h27-linux/continuation-cp3-source19/inventory.json)
[Settled receipt and account](wj-h27-linux/continuation-cp3-source19/post-publication.json)
[Launch hold](wj-h27-task3-launch-hold-v2.json)
