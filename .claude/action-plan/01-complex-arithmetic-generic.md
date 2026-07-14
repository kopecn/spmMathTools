---
chunk: 01-complex-arithmetic-generic
status: pending
depends_on: []
audit: ../review-for-fixes/2026-07-11-swift-audit.md §F6
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 01 — Generic Complex arithmetic

**Deliverable:** `Complex<T>` arithmetic for every
`T: BinaryFloatingPoint & SIMDScalar` — the prerequisite for the FFT fix.

## Files

- Edit: `spm/Sources/spmMathTools/FoundationMathTypes/Complex/Complex+Arithmetic.swift`
- Create: `spm/Tests/spmMathToolsTests/ComplexArithmeticTests.swift`

## Design constraints

1. Re-gate the existing Float-only extension (`:7`,
   `extension Complex where T == Float`) to
   `extension Complex where T: BinaryFloatingPoint & SIMDScalar` — the
   SIMD-backed bodies (`SIMD2` adds/multiplies) are already generic-safe.
   Keep `@inlinable`.
2. Surface unchanged: `+ - * /` (complex⊕complex, complex⊕scalar both
   orders where they exist today), prefix `-`, compound assignments. Do
   not add new operators.
3. Behavior for `Complex<Float>` must be bit-identical (existing operator
   semantics are the contract).

## TDD steps

1. Failing tests: `Complex<Double>` multiply — `(1+2i)(3+4i) == -5+10i`
   exactly; division round-trip `(a*b)/b == a` (rtol 1e-12); conjugate
   product `z * z.conjugate == |z|²` real; the same cases for Float pinned
   equal to current behavior (write these against Float FIRST, before the
   re-gate, so they document the incumbent).
2. Re-gate the extension. 3. `make build test` green.

## Acceptance criteria

- [ ] `grep -n "where T == Float" spm/Sources/spmMathTools/FoundationMathTypes/Complex/Complex+Arithmetic.swift` → no hits
- [ ] Double and Float arithmetic tests pass
- [ ] `make build test` passes

## Out of scope

The FFT (chunk 02); `Complex+Analysis.swift` (already both scalars);
adding Complex features Swift's base type lacks.
