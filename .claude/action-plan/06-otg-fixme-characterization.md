---
chunk: 06-otg-fixme-characterization
status: pending
depends_on: [04, 05]
audit: ../review-for-fixes/2026-07-11-swift-audit.md §F7
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 06 — Characterize the ThirdOrderStep2 FIXMEs

**Deliverable:** the two bare `// MARK: - FIXME` markers
(`OTG/position/PositionThirdOrderStep2.swift:1355, 1362`) each become
either (a) a passing test + explanatory comment, or (b) a precisely
documented known defect with a failing (skipped, linked) test.

## Files

- Edit: `spm/Sources/spmMathTools/FoundationMathTypes/OTG/position/PositionThirdOrderStep2.swift`
  (comments; code ONLY if a bug is proven)
- Create: `spm/Tests/spmMathToolsTests/OTGTests/OTGStep2FixmeTests.swift`
- Edit: `.claude/specs/known-issues.md` (add/close entries)

## Design constraints

1. First, archaeology: `git log -L1340,1380:…PositionThirdOrderStep2.swift`
   (adjust path) to find what the FIXMEs marked when introduced.
2. Construct inputs that exercise those exact branches (instrument
   temporarily with a counter — not print — to prove the branch fires;
   remove instrumentation after).
3. If output for those inputs violates an OTG invariant (limits exceeded,
   target missed, duration wrong vs Step1 optimal): that is a real bug —
   fix it ONLY if the fix is local to the branch; otherwise document
   precisely (inputs, expected, actual) in known-issues.md and leave a
   linked, `.disabled(...)`-annotated swift-testing case.
4. If output is correct: replace the FIXME with a comment stating what was
   verified and by which test.

## Acceptance criteria

- [ ] `grep -n "FIXME" spm/Sources/spmMathTools/FoundationMathTypes/OTG/position/PositionThirdOrderStep2.swift` → no hits
- [ ] Each former marker maps to a named test (passing, or disabled+documented)
- [ ] Truth-table/comprehensive/continuity suites unchanged
- [ ] `make build test` passes

## Out of scope

Rewriting solver branches beyond a proven-local fix; performance.
