# M8-R2 operational policy addendum, version 2

This additive record preserves the immutable candidate-contract-v2 and
operational-policy-v1 bytes. It documents the actual narrower implementation
guard and the boundary between residual-domain checks and extra constitutive
diagnostics before the next producing experiment epoch. It changes no
physical coefficient, residual, certificate bound, field value or tolerance.

The `R2LocalV1` implementation guard additionally requires **sampled
K/K0 in [.99,1.01]**, evaluated over its complete physical sampler inventory.
The mathematical R1 outer box still concerns the local outer k/K0 and actual
epsilon/epsilon0, not that sampled inner K. The mathematical inner law
permits positive sampled K, with [.5,1.5] used in its bounded kernel controls.
The narrower sampled-K condition is therefore an explicit operational
restriction, not a revision of the fixed R1 outer box or its interval proof.
Report the sampled-K range, strict scalar sign, and this extra guard margin
separately. The guard's finite covariance comparison uses the existing
`-1e-10*K0` tolerance; report the actual independent minimum and its strict
sign rather than reading an exact nonnegative value out of that tolerance.

Static review of the supported flat 2D Upstream+Capped, serial full system
confirms that ScalarDomain checks transported k at source, mass and
prescribed-load locations, and independent outer epsilon wherever these
locations are outer. It also checks interface/transfer locations;
native GLL and declared dense reconstruction; all physical matching and
sensor transfers; and the exact GLL evaluation nodes used by the dependent
trace projection. It checks sampled K at all associated x anchors. It then
rebuilds the actual candidate-dependent current volume rays and both baseline
and current dependent-trace normalization rays, checks every represented k
before normalization, recomputes b/lambda, and checks the projected target
trace as well as the independently retained actual EPS field. Positivity at
these finite points does not prove positivity throughout polynomial elements.

The actual coupled residual retains the Upstream epsilon-H coefficients as
independent unknowns. `trace_target_state` is a separate checking vector;
it is not substituted into the coupled stage. Constraints first apply the
same affine external lift as the residual, then compute the dependent trace.
The matched target keeps h0, target time and prescribed data unchanged at
every alpha; auxiliary solutions supply guesses only. The native integrator's
Advance() and accepted history do not run in this comparison mode.

The broader constitutive Audit additionally calls DependentEpsilon at inner
native/dense and other diagnostic reconstruction anchors. That diagnostic
path may create complete normalization rays beyond the current volume and
projected-trace rays required by the residual. Those extra rays must not be
mistaken for already checked residual evaluation sites. The independent
preflight reconstructs them as a declared diagnostic superset. If their
normalization fails after a coupled root converges, preserve the convergence
state, original residual and scalar-domain evidence, record the diagnostic
exception, and classify the root as **unqualified**. Do not replace fields,
expand the box, or discard the root by allowing a reporting exception to erase
its stage record. Root owns and verifies the corresponding exception capture
in the producing equation-system epoch; this document alone does not assert
that the implementation has passed that check.

The preflight supervisor prepares only outer6 with SaturatedR2 and Linear,
each at base/high source quadrature, and executes zero steps. It records
actual initial-state identity across laws, full scalar/ray/trace/covariance
and applicability margins, unanchored loads and source/energy/traction,
plus the actual transition sensor/chi/profile changes on common physical
anchors. Outer6 changes the normal support of the epsilon-H trace relative
to the fixed r=2.5 sensor; its effect is a measured representation change,
not a changed sensor definition or a presumed identical R1 history.

The supervisor records solver initialization elapsed time and any native
separate ScalarDomain/residual timing exported by the producing epoch.
Absent such native timing it reports that measurement unavailable, never
attributes the whole startup time to one domain evaluation. The first actual
auxiliary level within the bounded stage family supplies fresh-J cost; no
extra Jacobian probe is allocated. The earlier volume/dense scaling estimates
exclude repeated complete-domain screening. Current-iterate action checks are
recorded by the native solver and independently enforced in analysis against
the unchanged threshold; they are not a native hard-stop gate. The focused
kernel/operator/pressure verification precedes the stage family.
No preflight is launched until root records passing operator/action gates
and pins the next clean source/build/runtime. CPU49, one worker, six numerical
thread variables equal to one, shared validation lock and before/after full
runtime/source/input hashes apply. Historical builds and results remain intact.
