# H20 compact evidence: configuration stop

Authoritative Linux root: `/home/jay/build/wj-h20-three-level-20261115`.
Start with [component gates](publication-v2/component-gates.json),
[first failure](publication-v2/first-failure.json),
[planned/executed counts](publication-v2/planned-versus-executed.csv),
[time pairs](publication-v2/physical-time-pairs.csv),
[old/new lineage](publication-v2/initial-and-anchor-lineage.json), and
[closed resources](closure-v1/verification.json).

One M governor launch rejected ACK timeout 180 > native limit 120 before restore
or advance. Zero governing attempts, zero new CP0/build and one protected strict
read of the old H18 M step4/T=.01 anchor occurred. Four reserved/unknown slots stay
charged; later branches were not launched. Corrected timeout/preflight passed
passive tests only. All H20 temporal results are NOT_RUN.

[Index-v1](authoritative-index-v1.json) binds publication-v2 and separate closure.
[Local verification](local-verification-v1.json) records exact archive/member
hashes. Publication-v1 is a retained failed prefix, not an authoritative result.
Its original prefix manifest and publication-code/source metadata are in
[transport metadata](transport-metadata). Copied evidence and derived compact
views are labelled; the optional raw strict gzip remains remote-only. No active
operation log was sealed as completed evidence.

The publication accounting copy is an OPEN snapshot. Final accounting is CLOSED
at 714.056254596/7200 s. Earlier ledgers and heavy scientific files are unchanged.
No scientific contraction/performance figure was produced from missing results.
