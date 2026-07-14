---
chunk: 13-float-double-collapse
status: pending
depends_on: [01, 02, 03, 04, 05, 06, 07, 08, 09, 10, 11, 12]
audit: ../review-for-fixes/2026-07-11-swift-audit.md §F13
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 13 — Collapse Float/Double operator duplication (mechanical pass)

**Deliverable:** scalar-specialized operator extensions (`where T ==
Double` + `where T == Float` twins) become single generic
`where T: BinaryFloatingPoint & SIMDScalar` implementations. Runs LAST —
chunks 01–12 (especially the new test suites) are its safety harness.

## Files

- Edit: the twin-extension files under
  `FoundationMathTypes/Spatial/Operators/`,
  `FoundationMathTypes/Waveform1D/Operators/` (Double/Float pairs only —
  `where T == Int` extensions stay), `Complex/` is already done (chunk 01).
- No test files change (that is the point).

## Design constraints

1. **Mechanical rule:** a pair collapses ONLY if the two bodies are
   textually identical modulo the scalar type. Diff each pair first; any
   pair that differs beyond the scalar is skipped and listed in completion
   notes (it hides a real divergence).
2. Where a body needs scalar-specific literals (`1e-10` vs `1e-6`
   tolerances), keep the specialization for that member only — do not
   parameterize tolerances speculatively.
3. Work file-by-file, gate after each file (`make test`), so a miscollapse
   is bisectable.
4. Preserve `@inlinable`; expect no performance change (generic
   specialization handles it) — if the OTG performance suite regresses
   > 10%, revert that file and note it.

## TDD steps

1. Record full-suite pass state. 2. Collapse file-by-file with a gate run
   between files. 3. `make build test` green, results identical.

## Acceptance criteria

- [ ] `grep -rn "where T == Double" spm/Sources/spmMathTools/FoundationMathTypes/{Spatial,Waveform1D}/Operators` → only members justified in completion notes
- [ ] Zero test-file diffs; identical test outcomes
- [ ] `make build test` passes

## Out of scope

Int-specialized extensions; DSP `Extensions/` internals; API additions.
