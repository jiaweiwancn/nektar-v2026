# H2h execution plan

13 November 2026, Asia/Shanghai. Reviewed parent f30c309d6944dc7a117f14d6d7e2b929bbcf0d54.

The standard-versus-CLS implementation/validation task remains open. H2g's two
calibrations, physical equations, mesh, boundary roles, references and runtimes
remain immutable. H2h freezes the exact archived standard-B velocity/pressure
and changes only its permitted outer scalar coefficients in the existing 8+8
conditioned scalar space. The origin is B; all envelopes refer to standard cold.

1. Verify indexed state/runtime/data identities and record policy before solving.
2. Export actual native B gradients, B/cold scalars and saved-basis responses
   at the complete native inventory and GL20. Verify zero/nonzero affine
   reconstruction, protected rows, profiles and the unchanged velocity.
3. Run exact zero-margin linear phase I with covariance slack alone and all
   scalar bounds hard. Verify primal/dual bounds, residuals and gap. If needed,
   use bounded supporting-cut generation for the convex positive margin, then
   an optional secondary coefficient-distance objective after feasibility.
4. Select and preserve one terminal candidate. Perform one fresh GL24 audit;
   GL24 cannot be used for search or tuned against afterward. A passing native,
   GL20/GL24 and exact-state/runtime receipt is required for H2h issuance.
5. Only a qualified standard seed enters 8/16 steps, coarse 4+4 restart and
   curved rollback. A first rejection is terminal. CLS remains archived/blocked.

Maximum optimizer work is 24 LP rounds, 80 complete scans, 300 seconds; at most
three native candidate checks. No basis enrichment, envelope widening, fit or
closure changes, joint velocity search, or long continuation is authorized here.

After H2h the comparison remains open for common flat spatial qualification,
bounded developed-flow continuation, temporal/spatial sensitivity and held-out
experimental validation. EP-AWF and unrelated models remain deferred.
