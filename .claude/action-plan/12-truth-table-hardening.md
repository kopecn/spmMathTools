---
chunk: 12-truth-table-hardening
status: pending
depends_on: []
audit: ../review-for-fixes/2026-07-11-swift-audit.md §O2
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 12 — Truth-table resource hardening + numeric-golden export

**Deliverable:** missing truth-table resources FAIL instead of skipping,
and the 32-case numeric golden becomes a committed JSON other ports can
consume (py-MathTools' OTG oracle chunk expects exactly this).

## Files

- Edit: `spm/Tests/spmMathToolsTests/OTGTests/OTGComprehensiveTests.swift`
  (:406, :447)
- Edit: `spm/Tests/spmMathToolsTests/OTGTests/OTGTruthTableTests.swift`
- Create: `spm/Tests/spmMathToolsTests/OTGTests/truthTables/otg_numeric_truth.json`

## Design constraints

1. The two `XCTSkip("… Run generateRandomTrajectories.py first.")` become
   `XCTFail`/`Issue.record` — the JSONs ARE committed under
   `OTGTests/truthTables/`; a missing resource means bundling broke, which
   must be loud. Verify the resources load via `Bundle.module` (add
   `resources:` to the test target in `Package.swift` if they're currently
   loaded by path — check first).
2. Export the hardcoded 32-case array (inputs + `expectedDuration` +
   `expectedTimeIntervals`) to `otg_numeric_truth.json` with a `_meta` key
   naming the source (file + array symbol). Then make
   `OTGTruthTableTests` LOAD from the JSON and assert identical results —
   the JSON becomes the single source, the Swift array is deleted
   (Define Once).
3. JSON schema: `{"_meta": {...}, "cases": [{"input": {…InputParameter
   codable shape…}, "expectedDuration": …, "expectedTimeIntervals":
   [7 floats]}]}` — reuse the existing Codable keys so py-MathTools'
   `InputParameter.from_dict` reads it unmodified.

## TDD steps

1. Write a test asserting the JSON loads and has 32 cases; export; rewire
   the truth-table test to the JSON; delete the array.
2. `make build test` green — truth-table results identical.

## Acceptance criteria

- [ ] `grep -n "XCTSkip" spm/Tests/spmMathToolsTests/OTGTests/OTGComprehensiveTests.swift` → no hits
- [ ] `otg_numeric_truth.json` committed, 32 cases, `_meta` present
- [ ] The hardcoded Swift case array is gone; truth-table assertions unchanged in count and outcome
- [ ] `make build test` passes

## Out of scope

Adding new truth cases; regenerating the random corpora; py-MathTools.
