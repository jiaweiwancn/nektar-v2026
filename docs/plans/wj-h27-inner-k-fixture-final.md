# WJ–AWF inner-k fixture closure

| Decision | Result |
| --- | --- |
| MECHANISM_DIAGNOSIS | PASS |
| COLUMN_REMEDY | SUPPORTED |
| INNER_K_2D_FIXTURE | FAIL |
| PRODUCTION_MODEL | v2_UNCHANGED |
| TASK3_READY | NO |

Auxiliary inner-k diffusion and the historical cutoff reduced the actual
2D inner-domain departure from .13279 to .08419 in Kin RMS on the matched
mesh. Refinement gave .08402. Both exceed .05. The largest refined owner
sample was .26277 Kin, above .10. Independently assessed steady equations
and current guards passed. Retain the v2 production law and keep row D
as a diagnostic.

## Signed mechanism and column result

At the ground the historical profile has $k'_A=0$ and $k''_A=-111.552029$.
The v2 diffusivity is 19.7362 times the historical auxiliary diffusivity.
The resulting diffusion divergence is −.252986841 and production minus
destruction is +.0142089175, giving a v2 tendency of −.238777923 in code units.
This is the interior tendency on the inherited profile. The imposed
inlet trace stays fixed.

| Fixed choice | Ground tendency | Kin-column k RMS / MAX | Column $Q(H)/A$ |
| --- | ---: | ---: | ---: |
| A scalar Cd, v2 cutoff | −.238777923 | .099384 / .326425 | ≈0 |
| B auxiliary, same v2 cutoff | +.001390506 | .187174 / .422194 | +.00031460 |
| C scalar Cd, historical cutoff | −.214167815 | .116530 / .486285 | +.04451675 |
| D auxiliary, historical cutoff | 0 | .018939 / .021467 | +.04487878 |

At the ground the fixed-v2-cutoff diffusion change contributes +.24016843.
The direct destruction change at fixed auxiliary diffusion contributes
−.00139051. The scalar diffusivity's response to the cutoff contributes
the four-choice interaction −.02600061. Row D recovers the lower-layer
historical balance to 7.49e-16 in the fixed action scale. The remaining
full-WJ production increment integrates to +.01731825 A on the same fixed
profile. The F1 tests retained physical constants, full current stress,
the outer scalar equations and outer-only epsilon. The finite column
root and guard inventory spans the admitted [.9,1.1] Kin interval.
[F0 signed diagnosis](wj-h27-inner-k-fixture-f0.md)
[F1 fixed-choice columns](wj-h27-inner-k-fixture-f1.md)

The diagnostic uses $D_{k,w}=\nu+\nu_A(K,r)$,
$c_A(K)=H\exp(1-J(K,H))$ and $E_w=A(K)/\max(r,c_A(K))$.
Here $a=\kappa\sqrt{K/q_e}$, $d=\nu-a z_v$ and
$J(K,r)=\log(1+ar/d)+\nu/(d+ar)-\nu/d$ on the original roughness
branch. This matches the depth-integrated auxiliary parallel production
at each K. The F1 full-stress column calculates $Q(H)$ and $U(H)$ as
outputs, giving $U(H)=.43092191$ for row D. The actual 2D interface solves
for common diffusive flux with distinct owner gradients, and the outer
epsilon trace remains $A(K)/H$.

## Native fixture results

The isolated WJH27InnerKFixture steady entry point uses the same coupled
operators as the ordinary application. F2 reproduced the source19
four field and advection history hashes under the default policy.
The changed signed actions passed within 7e-15 and the two H owners
matched a nonzero diffusive flux with unequal diffusivities. The
candidate's separate reference and scientific signature reject the
opposite policy's accepted histories. The source19 selected fields
entered the matched case as a fresh steady guess. [F2 parity and action checks](wj-h27-inner-k-fixture-f2.md)

Both 2D cases retained the source19 inlet, top, ground, outlet traction
and scalar boundary loads. The native solves advanced zero physical time.
The directional grid added 20 x faces near the early ground layer, the
inlet/H corner and the downstream width transitions. It retained all
original radial faces and limited adjacent x widths to a ratio of two.
All 29 frozen stations were assessed with their original owner rules.
[F3 admission and norms](wj-h27-inner-k-2d-admission.md)

| Physical comparison to the historical inlet | Inner-k RMS / Gauss MAX in Kin | Inner-U RMS / MAX | Outer epsilon RMS / MAX |
| --- | ---: | ---: | ---: |
| Source19 v2, 28×100 | .132795 / .441427 | .005338 / .011541 | .009219 / .071707 |
| Row D, 28×100 | .084193 / .302640 | .003584 / .007574 | .002949 / .020718 |
| Row D, 48×100 | .084020 / .262304 | .003580 / .006913 | .002950 / .020593 |
| Required limits | .05 / .10 | .01 / .02 | .10 / .20 |

The extra owner samples give .304268 at x=7.5 on the matched row D grid
and .262771 at the ground and outlet on the directional grid. The
directional case fails the inner-k envelope at 27 stations. Its
independent momentum, continuity, k and epsilon residuals are 4.88e-14,
3.02e-15, 4.78e-13 and 1.91e-13. The native minimum k and epsilon are
.00880770 and 1.74415e-5. The minimum covariance/k is .280839 and the
supported principal response divided by mu is .00428821. The independent
cubic-root error is 1.01e-15. U, W, outer epsilon and ground shear pass
their physical envelope checks.

The directional comparison gives inner-k RMS/MAX differences of
.007045/.124208 in Kin. Both miss .001, and the original velocity-gradient
and traction comparison limits also fail. The higher-quadrature q12×16
case and the tighter 3e-12 nonlinear and 3e-13 linear case each return the
same coefficient hash as the q8×12 directional case, with independent
current checks passed. Exact polynomial integration reproduces .007045.
Gauss20 and Gauss40 with both owner endpoints give the same difference
peak of .124208 at x=7.5, r=0. From x=0 to 15, band RMS remains .090–.104
on the refined mesh. The oscillations shrink with refinement and the
downstream departure persists. The sampled extrema remain spatially
unresolved. [Station RMS and MAX figure](wj-h27-linux/inner-k-f3/station-inner-k.png)

## Remaining transport imbalance

In the refined field the current parallel auxiliary budget integrates to
zero within 4.21e-15 in units of $25A(K_{in})$. Actual production minus
destruction is +.0106750 in the same units. The full-WJ excess over the
auxiliary stress contributes +.0590641, the current shear shift from the
auxiliary parallel state contributes −.0483922, and isotropic work
contributes +.00000315. Below .4H the shear shift alone contributes
+.0480846. Orders 20 and 32 agree within 1.18e-16 for this split.
Restoring the historical fixed-profile balance leaves these current
stress and transport changes in the coupled equations.

The refined inner layer exports +.0036576 by net streamwise advection
and +.0075596 by upward diffusion at H in the same units. H advection
is −.0005333. The outer H owner gives +.0075466 for upward diffusion.
With the column sign convention $Q=D_k k_r$, the actual inner H mean
is negative, whereas the full-stress reference column gives positive
$Q(H)/A=.04487878$ at Kin. The steady k weak-reaction balance is
1.62e-12 in fixed budget units, and the balance using all physical
owner gradients is 1.48e-5. The H traction jump has RMS .0002223 and
sampled MAX .015488 in its fixed shear scale. Both H owners and their
one-sided gradients remain visible in the evidence.
[2D assessment and source split](wj-h27-linux/inner-k-f3/assessment-summary.json)

## Tested identities and reproducibility

The reviewed branch started at 672baf2b803d3282618d147b40c16c8d9c7c0a77.
The separately verified F0, F1 and F2 checkpoints are 59c237bee23654f497cee4117abe60e96c519860,
faa39d77d537ff7f9af948397b97f0807aea5aa0, and ff10b980cfafad8c1f14ad4631f1a75ead970419.
The retained production source19 ID is
56c35390fdf092180913f60c3ae348ca10cdcdf4e6b889c3afd589eddb511acf and its ELF
hash is 37307b1662d694124326568fd1137050a56955928c48f6406989a171596fc52b.
The v21 diagnostic source ID is
f6000e775564bbdb665d8cce292511eeb39324ea68509309f75216de145ac90b and its ELF
hash is 148e33dc937beb6998509de194387f5714f958df867c91ec1baa33bc171d3706.
The historical fixture reference is
cdbabf20a178cf8d1ea89e9160517b8e2276050f1efd5dc0072e8df6a61cb2ae and retains
the boundary table b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d.

The Linux evidence records complete command arrays and receipts. The
Release target was WJH27InnerKFixture with four build workers. Native
physical runs used WJH27InnerKFixture steady-fixture CASE.xml NEW_DIR
from runtime-v21 in the existing H27 Linux root. Science ran serially
with one rank and all six numerical-library thread settings at one.
The matched run took 466.450 seconds at 2.874 GB peak tree RSS. The
directional run took 627.442 seconds at 5.506 GB peak tree RSS including
startup and field export. The observed worker peak for the package was 15.

The 2D package phase charged 1,927.847 aggregate seconds. Across F0–F3,
2,510.914 of the permitted 6,000 seconds were charged, including failed
checks and publication. The existing account remains open at
22,283.291 of 43,200 seconds with 17,316.709 non-protected seconds and
3,600 evidence-protected seconds. Settled source/build/runtime and
science storage were 3.260 GB and 3.775 GB, with 68.079 GB free at
settlement. The Linux primary worktree and the deployed source19 ELF
were verified in their original state. [Settled account](wj-h27-linux/inner-k-f3/settled.json)

The Linux archive retains the original failed input and signed-source
checks, their corrections, and the failed first matched-geometry check.
That comparison now checks original inventories before the independent
oracle snaps domain endpoints. Exact raw fields, owner observations,
commands and audits are bound by the per-phase inventory hashes.
The F3 archive SHA-256 is
b862f3229da94966c059a7187b24a285b72f386b5947566839b0cc350765a764.
[F0 evidence](wj-h27-linux/inner-k-f0/inventory.json)
[F1 evidence](wj-h27-linux/inner-k-f1/inventory.json)
[F2 evidence](wj-h27-linux/inner-k-f2/inventory.json)
[F3 evidence](wj-h27-linux/inner-k-f3/inventory.json)
