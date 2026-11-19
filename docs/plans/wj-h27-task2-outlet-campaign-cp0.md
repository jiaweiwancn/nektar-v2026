# Steady flat outlet B qualification campaign

The campaign resumes from `d28a8f9c46d600e70cee18172f6452faab19f08d`. The fetch found the same tip and no later commits. The registered Linux application, reference table and every audited mesh and typed state match the preceding hold. The prior work established physical and current-equation acceptance on Bf, with 12 spatial misses against Bn2 and no misses in the experimental region. The terminal x control passed the spatial comparison and exceeded the independent momentum limit. [Prior outcome](wj-h27-task2-outlet-qualification-decision.md) · [Pinned identities and gates](wj-h27-linux/task2-outlet-campaign/cp0/outlet-campaign-plan.json)

The native model remains WJNativeCoupledV2 with ExactInnerFiniteNuV1, the historical outer inlet, ScalarCdV2, finite molecular viscosity, current-K cutoff, sampled wall force, full production and the implicit positive viscosity with its complete stress remainder. H=.15 code represents 6.0 mm above the modeled wall. The base uses Q3/Q2, native 6×8 quadrature and linear/coupled tolerances of 1e−12/1e−10. Each prepared session uses a typed zero-time seed, fresh lifts and fresh histories. FullTractionV1 is the default. ZeroTransverseV1 is the explicit flat candidate. [Model and acceptance contract](wj-h27-task2-outlet-qualification-contract.md)

## Resources and opening experiment

The locked amendment raises the cumulative ceiling from 64,800 s to 108,000 s. The old settlement charged 59,341.329 s. The amendment preserves every charge and the 3,600 s protected reserve. Two invocations with the exact authorization leave the same ledger bytes after the first apply. The watchdog keeps one science rank, six single-thread library settings, 80 total shared workers, four build jobs, 16 GiB build RSS, both 8 GiB disk limits and 15 GiB free space. A named physical solve with forecast peak above 10 GiB uses the smallest whole-GiB limit covering its forecast, up to 20 GiB, with live host and cgroup headroom. Each staged forecast includes all planned segments, audit, comparison and evidence, while reserving 300 s for final publication. [Amendment](wj-h27-linux/task2-outlet-campaign/cp0/outlet-campaign-amendment.json) · [Opening settlement](wj-h27-linux/task2-outlet-campaign/cp0/accounting.json)

Bfl reproduces the archived 111×94 design exactly, adding r=.0621875 and .0665625 to Bf. It preserves every streamwise edge and the companion interval above r=.06875. The first comparison isolates the previously untouched lower blend interval. Its prospectively tested prediction is contraction of the shear and stress peaks at r=.06 and the axial derivative at r=.06875. It also tracks the near-inlet radial derivative within the held lower owner. The signed shear diagnostic at χ=0 uses an exact centered product difference in Ur, Wx and ν+ν_A. A separate coefficient-based polynomial evaluator verifies the signed derivatives on both sides of the faces and at interior points. Spatial acceptance uses the original raw gradients and includes endpoints in MAX. [First design](wj-h27-linux/task2-outlet-campaign/cp0/outlet-campaign-bfl-design.json)

The full initial forecast is 3,036 charged seconds and 12.21 GiB, admitted against a named 13 GiB cap. The starts use native time limits of 257, 1,059 and 380 s, with a one-update first segment, continuation from a typed iterate and a bounded retry. The cost uses measured Bf and Bn2 receipts and reserves 552 MB of additional run storage. [Exact axes, resource anchors and staged reservations](wj-h27-linux/task2-outlet-campaign/cp0/outlet-qual-bfl-forecast.json)

| Work | Planning grid from Bfl | Planning charge | Forecast peak RSS |
| --- | --- | ---: | ---: |
| Full blend and both radial flanks, one split per cell | 111×132 | 4,264 s | 17.15 GiB |
| All experimental and terminal owners, with inlet if needed | 139×94 | 3,802 s | 15.30 GiB |
| Two successive H-only splits on the same original support | 111×96 then 111×100 | 3,101 s and 3,230 s | 12.47 and 12.99 GiB |
| Separate higher quadrature and tighter tolerance controls | 111×94 for each | 4,294 s and 3,036 s | 13.17 and 12.21 GiB |
| x=17.5 matched solve on the Bfl radial grid | 118×94 | 3,228 s plus build and parity | 12.98 GiB |

The table scales the opening cell count and measured normalized RSS. It reserves extra native storage and a quadrature cost factor of √2 for 8×12. The exact designs, costs and resource admissions are regenerated from completed measurements. The candidate must directly pass each applicable witness. If a newly changed band exposes a failing maximum, the next forecast binds the exact affected owner and the preceding failed comparison. A complete stage that fits the measured ledger advances under a fresh label.

## Gate matrix

| Gate | Acceptance | Initial state |
| --- | --- | --- |
| Full spatial and experimental window | All 690 retained signed rows and fixed x∈[−2.5,8] at the original RMS and MAX limits. Finer witnesses directly bind the proposed candidate. | Open |
| Radial coverage | Blend [.06,.1125], flanks [.05,.06] and [.1125,.11541061001735], previously exposed lower bands and both H-adjacent cells. | Open |
| Streamwise coverage | All 25 columns touching the 23 in-window stations and the two terminal columns. Include the inlet column if the local derivative remains sensitive. | Terminal comparison measured. Independent momentum and station coverage open. |
| Current equations and guards | Momentum and continuity ≤1e−9 and k/ε ≤1e−8. Weak accounting, boundary work, reaction, positivity, root, covariance, principal and current-K checks. | Bf passes |
| Physical and inlet projection | Original field limits, separate wall shear versus tau0 and Tref(Kin), and original U/k/gradient projection limits. | Bf passes |
| H fluxes | Shared mass and k-advection jumps ≤1e−10 and successive one-sided refinement per the frozen rule below. | Shared terms pass. One-sided convergence open. |
| Native quadrature | Exact selected mesh at 8×12 and base tolerances. Inner-k RMS and MAX ≤1e−4 and inherited limits elsewhere. | Pending spatial selection |
| Algebraic | Exact selected mesh at 6×8 and 3e−13 linear / 3e−12 coupled. Same field gates as the quadrature control. | Pending spatial selection |
| Outlet location | Native and independent x=15 parity, matched x=17.5 and frozen fixed-window physical limits with numerical uncertainty. | Pending CP2 |

At H, use the same upward normal for both owners. Compare k diffusion scaled by the reference A and the two total-traction components scaled by tau0. For each side and component, integrate on the common x intervals and retain endpoints in MAX. Keep the x grid fixed and halve both original H-adjacent widths twice while retaining their support. Require both RMS and MAX between the candidate and its finest H witness to be ≤.001. Require the second successive difference to be at most 80% of the first for each norm, with an absolute resolved floor of 1e−5 when both differences are below that floor. Require the fine one-sided flux jumps to be ≤.001 and report separately how the jumps change at both levels. Shared mass and k advection retain the original 1e−10 gate. A mesh change outside the frozen support opens a separate convergence sequence.

## Resume

Use the existing ext4 experiment root on `jay@100.88.70.60`. The operation receipts bind the independent evaluator and admission tools. Inspect the installed driver log and settled receipts before continuing. Keep prior operations immutable and use the existing typed-packet recovery path. The opening plan uses the deployed `outlet-control-v13/task2r_run.py` and `outlet-python-v28/h27/task2_outlet_campaign.py`. The same amendment command can reconcile a resumed session without new authority or duplicate credit. [Commands and archive identity](wj-h27-linux/task2-outlet-campaign/cp0/operations.json)

```bash
W=/media/jay/cb7dc460-f22f-4533-8329-245d10620f47/nektar-v2026-workspace
R="$W/runs/wj-h27-task2r-20261118"
source "$W/env.sh"
S="$R/tools/outlet-control-v13/task2r_run.py"
D="$R/tools/outlet-python-v28/h27"
F="$R/evidence"
export PYTHONPATH="$R/tools/outlet-control-v13:$D:$R/tools"
python3 -B "$S" "$R" balance
python3 -B "$S" "$R" amend-outlet-campaign "$F/outlet-campaign-authorization.txt"
```

**Campaign OPEN. TASK2_COMPLETE=NO. TASK3_READY=NO.**
