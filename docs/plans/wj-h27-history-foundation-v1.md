# Coupled history, retry and retention foundation

While the physical reference is refined, exercise unchanged source-v11 history
and failure handling on the tiny, genuinely evolving PolynomialScalarsV1 MMS
with separately active body forcing. This is not physical restart qualification
or H27-D closure.

Use nx=nr=2, dt=.02, body amplitude 1, linear tolerance 3e-13 and coupled
tolerance 3e-12. Compare a seven-step uninterrupted process with a fresh one-step
prefix followed by a COLD six-step process. Also compare a seven-step run with
the supported typed linear-failure injection at step two. It must preserve the
accepted first-step history and recover at the same dt exactly once.

The uninterrupted run writes recovery and scientific files at steps 2,4,6.
Observe actual deletion of recovery2 and retained recovery4/6, then continue
each retained checkpoint in a separate cold process to step7. Retain all three
scientific files and compare their actual states to the checkpoint payloads.
Info stride one records every accepted step's residual and admissibility.

An independent checker parses the binary payload, source/signature/checksum,
all four U/k/epsilon state/history sets and their raw hashes. Compare current
and previous fields plus both explicit histories at <=1e-8 using fixed U=1,
Pi=tau0, k=Kin and epsilon=Ain/.5 scales; report bit identity separately.
Incompatible-dt and deliberately corrupted owned checkpoint copies must reject
with their declared reasons. Independent current coupled equations on the
continuous/body, cold and retry endpoints remain separate checks.

Admit <=600 aggregate operation seconds for fixtures, eight small native
processes, independent equations and checker work, within the same continuing
H27 allowance. One rank, all numerical-library threads one, shared actual
workers<=72 and existing RSS/storage safeguards apply. Production physical
restart, terminal statistics, denser monitoring, storage and complete cost
remain open even if this foundation passes.
