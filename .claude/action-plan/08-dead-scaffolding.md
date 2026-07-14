---
chunk: 08-dead-scaffolding
status: pending
depends_on: []
audit: ../review-for-fixes/2026-07-11-swift-audit.md §F9
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 08 — Delete dead scaffolding

**Deliverable:** no stub/empty files shipping in the target.

## Files (all deletions)

- `spm/Sources/spmMathTools/FoundationMathTypes/Waveform1D/Operators/Waveform1D+Custom.swift`
  (18 TODO placeholders, zero API)
- `spm/Sources/spmMathTools/FoundationMathTypes/SpatialWaveforms/Operators/WaveformPosition+Arithmetic.swift`
- `…/SpatialWaveforms/Operators/WaveformQuaternion+Arithmetic.swift`
- `…/SpatialWaveforms/Operators/WaveformSpatialPose+Arithmetic.swift`
  (empty extension bodies)
- `spm/Sources/spmMathTools/FoundationMathTypes/PrecisionTime/Extensions/PrecisionTimestamp+Extensions.swift`
- `…/PrecisionTime/Extensions/PrecisionTimeInterval+Extensions.swift`
  (zero-byte)
- Edit: `.claude/specs/known-issues.md` — record that spatial-waveform
  arithmetic and precision-time extensions are intentionally
  not-yet-implemented (the stubs were the only marker of that intent).

## Design constraints

Verify each file defines zero usable API before deleting
(`grep -c "func\|var\|let" <file>` — comments only). If any acquired real
content since the audit, leave it and note it. Delete empty parent
directories left behind.

## TDD steps

1. Delete. 2. `make build test` green (nothing referenced the stubs).

## Acceptance criteria

- [ ] All six files gone; `find spm/Sources -name "*.swift" -size -10c` → no hits
- [ ] `grep -rn "TODO: Add custom" spm/Sources` → no hits
- [ ] known-issues.md carries the intent note
- [ ] `make build test` passes

## Out of scope

Implementing any of the stubbed features; test-side stubs.
