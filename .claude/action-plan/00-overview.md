---
plan: swift-audit-fixes
status: pending
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# Action Plan — Swift Audit Fixes (spmMathTools)

**Goal:** close every finding in
[`../review-for-fixes/2026-07-11-swift-audit.md`](../review-for-fixes/2026-07-11-swift-audit.md)
(F1–F13, O1–O5). The audit is authoritative; chunks back-reference its
finding IDs. Headline: the spectral surface (F1) is silently wrong today —
chunks 01→02→03 are the critical path. **Decision on record: Linux is a
required build target** (same decision as spmFoundationTools).

## Conventions every chunk inherits

1. **Gate:** `make build test` green at the end of every chunk (plus
   `make format`).
2. **TDD:** failing tests first for behavioral findings; regression tests
   named `test<FindingID>_...`.
3. **New tests use swift-testing** (`import Testing`) — audit O3 policy.
4. **Stay in scope:** only the chunk's file list; adjacent issues →
   completion notes.
5. **Reconcile `.claude/specs/known-issues.md`:** if a chunk fixes
   something that file records, update/remove the entry in the same chunk.
6. **Frontmatter:** tracking triad + `status: pending`; update as you work.
7. Faithful-numeric caveat: OTG and root-solver behavior is pinned by the
   existing test suites (OTGTruthTable/Comprehensive/Continuity) — a chunk
   that changes any of their outcomes is wrong until proven otherwise.

## Deferred (recorded, not chunked)

- **O5** TODO triage beyond what F1/chunk 02 unblocks — file follow-ups.

## Dependency graph

```
01 complex-arithmetic ─► 02 fft-butterfly ─► 03 spectral-tests
04 otg-debug-prints ──┐
05 crash-path-hardening ─┼─► 06 otg-fixme-characterization
07 packaging-and-linux   │        (06 after 04/05: clean solver first)
08 dead-scaffolding      │
09 qmul-dedup            │
10 dsp-tests-1  (after 02 — envelope/phase lean on working spectra)
11 dsp-tests-2  (independent of 02)
12 truth-table-hardening (independent)
13 float-double-collapse (LAST — needs 01–12 green as its harness)
14 otg-benchmark (optional, measure-first)
```

## Chunk index

| # | Chunk | Findings | Depends on |
|---|---|---|---|
| 01 | [complex-arithmetic-generic](01-complex-arithmetic-generic.md) | F6 | — |
| 02 | [fft-butterfly-restore](02-fft-butterfly-restore.md) | F1, F4(fft) | 01 |
| 03 | [spectral-tests](03-spectral-tests.md) | O1(spectral) | 02 |
| 04 | [otg-debug-prints](04-otg-debug-prints.md) | F2 | — |
| 05 | [crash-path-hardening](05-crash-path-hardening.md) | F3, F4, F5 | — |
| 06 | [otg-fixme-characterization](06-otg-fixme-characterization.md) | F7 | 04, 05 |
| 07 | [packaging-and-linux](07-packaging-and-linux.md) | F8, F10, F11 | — |
| 08 | [dead-scaffolding](08-dead-scaffolding.md) | F9 | — |
| 09 | [qmul-dedup](09-qmul-dedup.md) | F12 | — |
| 10 | [dsp-tests-1](10-dsp-tests-1.md) | O1(filter/envelope/phase) | 02 |
| 11 | [dsp-tests-2](11-dsp-tests-2.md) | O1(resample/window/trigger/zerox/generators) | — |
| 12 | [truth-table-hardening](12-truth-table-hardening.md) | O2 | — |
| 13 | [float-double-collapse](13-float-double-collapse.md) | F13 | 01–12 |
| 14 | [otg-benchmark](14-otg-benchmark.md) | O4 | 04 |
