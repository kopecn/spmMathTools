---
chunk: 05-crash-path-hardening
status: pending
depends_on: []
audit: ../review-for-fixes/2026-07-11-swift-audit.md §F3, §F4, §F5
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 05 — Crash-path hardening

**Deliverable:** no trap reachable from public API on runtime data.

## Files

- Edit: `spm/Sources/spmMathTools/FoundationMathTypes/OTG/CalculatorTarget.swift`
  (:208-210), `OTG/Trajectory.swift` (:55, :67, :122),
  `OTG/OutputParameter.swift` (:100 doc),
  `Waveform1D/Extensions/Waveform1D-Envelope.swift` (:285)
- Edit/Create tests beside the existing OTG + envelope tests.

## Design constraints

1. **F3** `blocks[divRem].a!.profile!` → `guard let` chain; on absence the
   calculate path returns `Result.errorSynchronizationCalculation` exactly
   as the surrounding error branches do (copy the adjacent early-return
   pattern in the same function). Same for `.b!`.
2. **F4** `Trajectory` `preconditionFailure`s → the not-computed state
   becomes representable: `stateToIntegrateFrom` (and the :122 sibling)
   return an optional or a documented safe zero-state consistent with how
   `atTime` treats an empty trajectory — read `atTime` first and match its
   convention; whichever you pick, the doc comment states it. NO trap.
3. **F4** `OutputParameter.swift:100`: fix the doc comment — it claims
   "Throws: preconditionFailure"; document the actual behavior after
   change (2).
4. **F5** `distances.min()!` → `guard let minValue = distances.min(),
   let index = distances.firstIndex(of: minValue) else { …empty-input
   behavior consistent with the function's existing empty-waveform
   handling… }`.
5. Existing OTG suites must be untouched by these changes (they never hit
   the trap paths — that is why the traps survived).

## TDD steps

1. Failing tests: (a) a synchronize state with an empty block interval
   returns the error Result (construct via the public calculate path with
   a crafted degenerate input if reachable; if not reachable publicly,
   test the internal function directly and say so); (b) `Trajectory`
   sampled before computation returns the documented safe value — no
   crash; (c) envelope on an empty/1-sample waveform returns the
   documented empty result — no crash.
2. Implement. 3. `make build test` green.

## Acceptance criteria

- [ ] `grep -rn "preconditionFailure\|fatalError" spm/Sources/spmMathTools/FoundationMathTypes/OTG` → no hits
- [ ] `grep -n '\.a!\|\.b!\|min()!' <the three edited source files>` → no hits
- [ ] All three regression tests pass; existing OTG suites unchanged
- [ ] `make build test` passes

## Out of scope

The FIXME branches (chunk 06); debug prints (04); API redesign of
Trajectory.
