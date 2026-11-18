# Active-v2 inner inlet fixture closure

**BLOCKED at inlet construction.** The tracked increasing-shear branch
reaches its supported principal limit well before matching the fixed
outer velocity. Even its shear-cap column misses the target. Keep the
active v2 closure, calibrated outer inlet, external loads and deployed
historical-inlet preset. No matched inlet was admitted.

| Gate | Status |
| --- | --- |
| INLET_CONSTRUCTION | BLOCKED |
| NATIVE_INLET_BINDING | NOT_RUN |
| COUPLED_CURRENT_EQUATIONS | NOT_RUN |
| INNER_K_PHYSICAL_DEVELOPMENT | NOT_RUN |
| OTHER_PHYSICAL_ENVELOPES | NOT_RUN |
| INNER_K_NUMERICAL_RESOLUTION | NOT_RUN |
| FULL_COUPLED_SPATIAL_QUALIFICATION | NOT_RUN |
| TASK3_READY | NO |

The reviewed branch starts at ec2c2a1496730ed6530782721e5d6ad2ac40d2b3.
The prospective contract was pushed first at
5a16fafc8696c7fb28810659be6d808162c954d3. The revised inner profile
permission supersedes the earlier whole-height freeze. This test retains
the existing table SHA-256
b1bac2202d5991865a19c193528182954211f11f954fb5ca9b3f23998262250d,
ScalarCdV2, finite viscosity, the current v2 cutoff .02954486569124273
at Kin, auxiliary viscosity, roughness, the amended tensor and all
outer and exterior data. [Prospective contract](wj-h27-active-v2-inner-inlet-contract.md)

## Constant-shear compatibility

The boundary value problem varies $T_{\rm in}$ independently of
$a(K_{in})$ and solves $U'=G$, $k'=Q/D_k$ and
$Q'=E-(T_{\rm in}-\nu G)G$. It requires U(0)=Q(0)=0,
$k(H)=K_{in}$ and $U(H)=.43724424118668737$. The physical
increasing shear branch limits $q=kG/E$ to $72/13$ where the
collar is active. At H, $k=K_{in}$, $E=A(K_{in})/H$ and
the full-WJ shear function is at most .3. The exact total-shear
ceiling there is

$$T_H^{\rm cap}=.3K_{in}+\nu\,\frac{E(H)}{K_{in}}\frac{72}{13}
=.0027362151252511993.$$

Joint DOP853 shooting with analytic endpoint sensitivities started
from the canonical $U(H)=.42150691440$ and advanced to .43527707534.
The target correction hit the ceiling. The fixed-shear solver eliminated
the ground-flux equation by rederiving the lower-layer production integral
for arbitrary T and checking it with split Gauss64. It searched inside
the interval for Q(H) implied by $0\leq P\leq T^2/(4\nu)$, retaining
invalid shots. In 31 shear samples it found 24 bracketed positive full
columns. The nearest shear-cap sample was 1e-12 fractionally below the
ceiling. Its U(H) was .435826822102018, short by .001417419084669.
Independent two-ended RK45 and companion-root stress evaluation
reproduced U(H) within 1.16e-13 and compared scaled profiles and
derivatives within 8.62e-12. The flux and shear errors were below
1.21e-13 and 1.51e-15. [Shear ceiling and matched integration](wj-h27-linux/active-inlet/branch-assessment.json)

The supported planar check at H brackets a downward crossing of the
existing 1e-8 principal-response threshold at $q=3.69980211$ and
$T=.002678099344$. Two fixed-shear columns straddling it gave
minimum principal/mu +1.54e-6 and −1.52e-6 over 516 sampled radii.
The upper shear-cap column had minimum principal/mu −.104313.
The positive near-crossing witness reaches only .428649580 at H,
short of the immutable target by .008594661. Its independently
checked derivative-integrated polynomial has maximum scaled error
5.57e-9. Its distinct diagnostic schema records the unmatched value.
The lower-shear follow checked whether the velocity could return
to the target in the other direction. It found a turn near
$T/\tau_0=.92$ at $U(H)\simeq.41304$, then negative principal response
at .915. The .91 shot lost positivity below the collar. From .905
down to .84 the connected shots reached their admissibility limit
without a ground-flux bracket. The broader scan supplied coverage
down to .05. The branch, positivity and directional searches have
finite coverage and leave global uniqueness unclaimed.
[Principal-response checks](wj-h27-linux/active-inlet/principal-summary.json)
[Lower-shear continuation](wj-h27-linux/active-inlet/lower-branch-assessment.json)

| H-minus witness or retained H-plus datum | Near principal limit | Near shear cap | Retained outer |
| --- | ---: | ---: | ---: |
| U | .428649580 | .435826822 | .437244241 |
| k | .00880901515 | .00880901515 | .00880901515 |
| E | .00231552405 | .00231552405 | .00231552405 |
| U radial derivative | .972518986 | 1.455831408 | .890243773 |
| k radial derivative | −.013792982 | −.025905600 | 0 |
| Finite-nu total shear | .002678098221 | .002736215125 | .002658181905 |
| Q(H) | −4.118451093e-5 | −7.735161887e-5 | 0 |

Here $Q=D_k k_r$ and physical upward diffusion is $-Q$.
The positive principal witness exceeds the original outer total shear
by 1.99163e-5 and has $\int(P-E)dr=4.11845\times10^{-5}$.
Including Q(H) closes the scalar budget to 1.11e-13 in units of A.
Both H-side epsilon derivatives are −.015436827.
These limiting column values are diagnostics for compatibility.

## Coupled evidence and handoff

| 2D case and comparison reference | Inner-k RMS / sampled MAX in Kin |
| --- | ---: |
| Old active v2 against historical reference, 28×100 | .13279450 / .44142666 |
| New-inlet active v2 against the requested frozen reference | NOT_RUN |
| New-inlet active v2 against historical reference | NOT_RUN |

The separate row-D diagnostic remains at .084020/.262771 under the
complete historical inlet with an alternative transport identity.
The retained F2 native default parity verified all four source19
field and advection hashes. No native inlet provider, accepted-history
identity, projected native lift or 2D field was added at this blocked
checkpoint. The witness interpolation result measures only its
polynomial artifact. [Source19 CP2](wj-h27-continuation-cp2.md)
[Retained F2 parity](wj-h27-inner-k-fixture-f2.md)
[Earlier row-D fixture](wj-h27-inner-k-fixture-final.md)

## Resources and evidence

The Linux evidence binds the source, commands, fixed table and raw
shooting attempts. The initial near-cap forward RK stage crossed the
branch boundary. The independent two-ended solution enforces the exact
H value directly and avoids that issue. The first witness interpolation
missed the derivative gate. Integrating the sampled ODE derivatives and
checking both endpoint values produced the accepted diagnostic fit.
All rejected trials and operation exits remain archived.

| Disjoint operation group | Charged Linux seconds |
| --- | ---: |
| Joint shooting, fixed-shear scan and lower-branch follow | 112.577 |
| Independent edge, principal and witness checks including their column controls | 40.852 |
| Publication | 1.923 |
| Native build and 2D solves | 0 |
| New package total | 155.352 |

The original shared account is open at 22,438.642919 of 43,200
seconds. The package has 5,844.647801 of its 6,000 seconds left.
The global account retains 17,161.357081 non-protected and 3,600
protected seconds, with zero live reservations. The maximum tree
RSS was 66,043,904 bytes and the observed host worker peak was 5.
Each operation used one process and the six single-thread settings.
The Linux primary checkout, deployed source19 ELF and table were
rechecked at settlement. The Linux archive of 42 verified members
has SHA-256
2c1b3f6eac218a7dfb8b6b91b03ad16cdeef26d0773d63f019d8529d622c3586.
[Evidence inventory](wj-h27-linux/active-inlet/inventory.json)
[Settled account](wj-h27-linux/active-inlet/settled.json)

## Next task

Review the exact local stationarity requirement before another coupled
run. A separately specified value-matched developing inner inlet could
retain active v2 and the original outer calibration while permitting the
physical streamwise balance at the inlet. Check its supported-principal
margin and both H owner fluxes before freezing its comparison profile.
