# H27 CP2: localize the approved v2 reference departure

CP1 is remotely verified at `971c50f848411994d4afceb643356821483b5cc6`.
The source19 coupled numerical foundation passes. The developed reference
remains unresolved. Its independent equations pass on 28×200 cells while its
inner-k RMS/MAX depart from the retained inlet by .13279468209091205/
.44168324293315253 in Kin units, exceeding .05/.10. The existing h comparisons
do not establish a selected-source error bound. [CP1](wj-h27-continuation-cp1.md)
[Physical checkpoint](wj-h27-cutoff-panels-checkpoint-v2.md)

Use the pre-existing physical source19 field and its exact independent
order48 report. The zero-time native reference holds the approved v2
equations, inlet, top, outlet, H and roughness. No solver advancement or
reference recalibration is needed to measure these balances.

For each original x column, independently integrate the inner-domain
production, dependent destruction, signed k departure and squared
departure. Retain radial bands bounded by the historical cutoff, .4H, h,
.75H and H. The independent volume rule also splits at the current cutoff
and its crossings of the original faces. Quantify each column's contribution
to the global squared departure and record the signed sampled peak and owner.

Use total current flux \(F=Uk-D_k\nabla k\), with \(U\) the vector velocity.
For each original cell, measure outward advection and diffusion separately
on every face. By column report the two x sections, ground and H fluxes,
the independent H+ flux, and the internal radial owner jumps. Keep the
actual sign of both fluxes at each face. Compute the defect between
outward physical flux and the actual integrated \(P-E_w\); a finite owner
gradient jump is an approximation diagnostic. Compare both Gauss order20
and order32 at fixed coefficients. Require accurate area and a repeatable
RMS, and bind the existing independent boundary and global weak budgets.
Use the measurements as localization and quadrature evidence on this mesh.

Evaluate the historically retained parallel inlet as a local diagnostic
under the active v2 scalar law. Use independently sampled gradients and
full collar production. Compare its historical and v2 cutoff, destruction,
diffusion and cumulative required \(Q(r)=\int_0^r(E_w-P)\,dr\) against
the gradient flux \(D_k k'_\mathrm{in}\). Keep the existing inlet profile
and actual H flux in the PDE. This test isolates source and normal diffusion
at fixed inlet data. Actual flow may also change through streamwise transport
and boundary reactions.

Admit the combined fixed-field order20/order32 analysis for up to 360 seconds
and evidence for up to 60 seconds on the live account. Keep one rank, all
six library threads one, 72 actual host workers, 8 GiB science RSS, 8 GiB
science/source-build storage and 15 GiB free. The account was OPEN at
17116.471991412807 seconds with 3600 protected. All outputs are create-only.
Use the results and the live balance to freeze a separate, complete bounded
spatial/quadrature/tolerance matrix before a native solve. The locality of
the physical departure alone does not establish which x cells need refinement.
Use previous directional numerical defects and owner gradient evidence too.
