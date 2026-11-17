# Approved-v2 streamwise reference refinement admission

The 700/2800-cell steady equations converge, but common-point inner-k
differences remain large at interior stations, including x=0 and 7.5. The
physical departures fail the fixed RMS/MAX envelope. The present evidence does
not separate streamwise Galerkin resolution from a limiting model adjustment.
Do not declare a resolved continuum obstruction or activate new physics.

Starting from the actual 2800-cell XML axes (28 x intervals, 100 radial
intervals), bisect every x interval only. The new 5600-cell calculation retains
all radial coordinates, Q3/Q2 spaces, native 8x12 split integration, v2 physics,
historical inlet, physical boundaries and reference table. Transfer the 2800-cell
field only as a fresh steady initial guess, rebuilding dependent Gamma and all
derived data. Source-v11 and FullGradientV1 remain pinned. Use linear 1e-12,
nonlinear 1e-10, Krylov dimension 80; physical time remains zero.

Assess the current equations independently, physical envelopes in every fixed
region/station, and all common-owner numerical fields against 2800 cells.
Locate signed ground/inner k variation and the original owner-gradient jumps;
do not average gradients, remove difficult stations or relax any envelope.
Any numerical stabilization is a separate prospective implementation and is
not authorized by this run's identity.

Measured 2800-cell native operation was 417.867 s with 2.716 GB peak tree RSS.
Forecast roughly 1000--1500 s and 5.5--7.5 GB for this x refinement; an actual
8 GiB watchdog remains binding. Admit generation <=30 s, native <=2100 s
(internal 1800 s), independent equations <=400 s, common-point comparison
<=450 s and localization <=120 s. The entire additional matrix is <=3100
aggregate operation seconds, including failures, in the continuing 43200 s
account. The last completed snapshot is 6947.866607507807 s; 3600 s remains
protected for evidence. No allowance is reset.

Use a single science rank, all six numerical-library thread counts one,
actual shared-host workers <=72, science/source/storage/free-space caps as
before. Keep outputs in the isolated healthy ext4 H27 namespace. No production
pilot or full Task3 horizon is admitted by this diagnostic.

A separate inactive mathematical factorial control is admitted for <=120 s
analysis (110 s internal), bringing this entire admission to <=3220 s. At Kin
only, compare the active scalar/auxiliary inner diffusivities and the selected
full-WJ zero-endpoint-flux/historical cutoffs. When holding the historical
cutoff, solve the counterfactual Qref(H) in the bounded [0,.4] A interval for
Qref(0)=0. Keep the complete finite-nu v2 shear/collar and k(H)=Kin. Measure
column controls and departures from the unchanged inlet. Record invalid
trials. This is a prospective mechanism study only: no runtime table,
new-model identity, executable equation or actual-flow boundary is changed.
Approval remains required for any further active scientific amendment.
