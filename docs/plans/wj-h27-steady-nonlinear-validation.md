# Optional complete nonlinear-state validation

`SteadyKrylovCallbacks::ValidateNonlinearState` has signature
`void(const SteadyVector&)`. When present, the solver calls it after `ValidTrial`
and before `Residual` for the initial state, each re-evaluated accepted Newton
base, and each nonlinear line-search candidate. A failed full guard therefore
cannot reach `TrueNorm`, `Observe`, best-state storage, or Armijo acceptance.
The derivative-probe base is also fully validated.

Finite-difference samples continue to call `ValidTrial` and `Residual` and retain
their existing central/one-sided/reduction and typed domain-rejection policy.
They do not call the additional complete nonlinear guard. Their recorded
`plusAdmissible`, `minusAdmissible`, and `halfAdmissible` flags certify domain
evaluation only; they do not claim complete principal admissibility or accepted
physical-state status. FD points never enter best-state storage.

A `SteadyInadmissibleTrial` from the complete guard rejects a nonlinear candidate
and records its reason, attempted/passed flags, and separate rejection count.
The same exception is fatal at an initial/re-evaluated base. Every unrelated
exception propagates unchanged. The generic solver has no validation cache.
The application may reuse complete validation only for exactly the same full
unknown vector; domain-only FD evaluations must not overwrite that certificate.
Final application-level physical validation remains the caller's responsibility.

A best candidate returned after bounded failure may still be Armijo-unaccepted,
as before. With the callback enabled it must have passed complete validation;
`bestFromAcceptedIterate` continues to distinguish these two provenances. A
failed initial validation produces no evaluated/best-state evidence. Absence of
the optional callback preserves the prior numerical algorithm and required
callback behavior.

The standalone probe retains all prior controls and adds rejected-best-state,
guarded half-step, fully validated but unaccepted best-candidate, no-op callback
equivalence, and three domain-only FD controls. Fatal controls cover initial,
derivative-probe base, re-evaluated accepted base, and unrelated candidate errors.
The independent checker uses exact scalar cubic stencils and nonlinear equations
to verify their evidence. This is a generic Newton evaluation optimization; it
does not change equations, scientific admissibility thresholds, or H27 model
qualification. Linux build/execution is required; Mac checks are static only.
