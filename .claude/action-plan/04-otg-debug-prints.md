---
chunk: 04-otg-debug-prints
status: pending
depends_on: []
audit: ../review-for-fixes/2026-07-11-swift-audit.md §F2
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 04 — Remove the 60 OTG debug prints

**Deliverable:** a silent OTG hot path.

## Files

- Edit: `spm/Sources/spmMathTools/FoundationMathTypes/OTG/CalculatorTarget.swift`,
  `OTG/position/PositionThirdOrderStep2.swift`, and every other file
  `grep -rn "print(" spm/Sources/spmMathTools/FoundationMathTypes/OTG` hits.

## Design constraints

1. Default action: **delete** each `print("[DEBUG] …")` line.
2. Exception: a print that carries genuine diagnostic value at a failure
   branch (e.g. the synchronize-failed dumps at `CalculatorTarget.swift:
   217-218, 766-769`) may instead move behind a single internal hook:

```swift
// OTG/OTGDiagnostics.swift
internal enum OTGDiagnostics {
    /// Caller-injectable; nil by default — zero cost in production.
    nonisolated(unsafe) static var handler: (@Sendable (String) -> Void)?
    @inline(__always) static func emit(_ message: @autoclosure () -> String) {
        handler?(message())
    }
}
```

   Keep the hook TINY; if `nonisolated(unsafe)` displeases the compiler
   settings, make handler a `let` configured via an init-time parameter on
   `OTG` instead — decide by whichever compiles cleanly first.
3. No behavior change: the full existing OTG test suite is the harness.

## TDD steps

1. Run the OTG suites, record pass state. 2. Remove/convert prints.
3. Suites unchanged; `make build test` green.

## Acceptance criteria

- [ ] `grep -rn "print(" spm/Sources/spmMathTools/FoundationMathTypes/OTG --include="*.swift"` → no hits
- [ ] OTG test results identical before/after
- [ ] `make build test` passes

## Out of scope

Fixing the branches the prints were annotating (chunks 05, 06); logging
frameworks (no swift-log dependency exists here — do not add one).
