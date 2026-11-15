# CT12 read-only reconstruction method

For each of the six captured CT11 corrections, restrict the recorded sparse
frozen matrix and consistent velocity mass with the actual native free IDs.
Check its four blocks against the separately captured G, D and C. Require C=0
and exclusion of the fixed-zero unused gauge slot. The physical pressure tests
remain free; the outlet total traction sets the absolute level.

With the captured residual convention the predictor is
`du*=−A⁻¹ r_m`, and the native pressure increment solves

```
S_alpha dpi = D du* + r_c
S_alpha p   = D [M_free⁻¹ (G p)] / alpha − C p.
```

Sparse direct mass/predictor actions use extended-precision CSR accumulation and
up to three original-system residual refinements. Diagonal mass is used only in
the pressure preconditioner. Restarted GMRES is bounded to at most 960 Arnoldi
iterations; its stopping flag is recorded independently from recomputed true
residuals. A work-limit flag is never called relative convergence. Tight repeats,
physical image differences, signed sums and native-increment reconstruction
provide the additional physical error tests. No dense Schur matrix or array of
all mass-inverse pressure columns is formed.

Map momentum loads with `P_alpha(r)=−S_alpha⁻¹ D A⁻¹ r`. Separate the directly
captured continuity image `S_alpha⁻¹ r_c` and any backward-storage/lifting
remainder from the spatial residual. On the same native integration rows, use

```
P_alpha(current spatial load)
  = P_alpha(sum continuous components)
  + sum_j P_alpha(initial represented_j − continuous_j)
  + sum_j P_alpha(current_j − initial represented_j).
```

This is a diagnostic linearity identity. No reference contribution is removed
from evolution. Signed native component vectors, local sparse supports,
reactions, natural top/outlet terms and owner/side metadata remain in the heavy
packets. The reorganization avoids separately solving and subtracting large
cancelling equilibrium pressure images. Component pressure differences and the
complete signed sum are checked in physical norms on the actual target's fixed
pressure observation space. Whole/inner/outer means and variations, and the
one-sided interface pressure jump, retain the absolute level.

Sum the captured increments from the unchanged initial absolute pressure and
compare with every corresponding accepted state. Independently check that
pressure is unchanged across corrected-flow, post-flow, pre-k and accepted
scalar stages. Report imposed dependent-epsilon trace and refreshed post-k lag
separately. A representative exact frozen coupled correction uses `D A⁻¹G−C`
only as a read-only comparison. It never replaces the governing mass correction.

For the initial representation's radial localization, stream each band's sparse
native supports into one global momentum vector and map it through the fixed
consistent-mass small-step image. Signed alignments, not unsigned stress norms,
identify reinforcing and cancelling supports. Interface flux observations have
opposite outward orientations; they are diagnostics, not extra shared-row loads.

All load/block reconstruction limits remain 2e−10, and the pressure-image,
native-increment, tighter-repeat and signed-sum uncertainty limit is 1e−8 on the
u_tau² scale. Held-state governing-quadrature controls preserve coefficients and
the pressure observation operator. Any failed capacity control is a blocked
quadrature result, not a changed physical trajectory or permission to continue.
