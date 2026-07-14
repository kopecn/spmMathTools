---
chunk: 09-qmul-dedup
status: pending
depends_on: []
audit: ../review-for-fixes/2026-07-11-swift-audit.md §F12
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 09 — Deduplicate the quaternion kernels

**Deliverable:** `_qmul`/`_qrot` defined once.

## Files

- Create: `spm/Sources/spmMathTools/FoundationMathTypes/Spatial/Operators/QuaternionKernels.swift`
- Edit: `Spatial/Operators/Quaternion+InlinableQuaternion.swift`,
  `Spatial/Operators/SpatialPose+InlinableQuaternion.swift` (delete the
  copies; forward call sites)

## Design constraints

1. New home: internal `enum QuaternionKernels` with `@inlinable
   @inline(__always) static func qmul(_:_:)` and `qrot(_:_:)` for
   `SIMD4<Double>`/`SIMD4<Float>` (or one generic
   `where T: BinaryFloatingPoint & SIMDScalar` if it compiles to the same
   code — prefer generic, fall back to the two overloads if SIMD generics
   fight back; note which).
2. Bodies are moved verbatim from `Quaternion+InlinableQuaternion.swift:13,41`
   (the `SpatialPose` copies at `:12,:30` must be textually identical —
   diff them first; if they differ, STOP: that's a correctness finding).
3. Call sites in the operator extensions switch to
   `QuaternionKernels.qmul(...)`; the old private statics are deleted.

## TDD steps

1. The existing Quaternion/SpatialPose operator suites (~56 + 8 tests) are
   the harness — record pass state before, identical after.
2. `make build test` green.

## Acceptance criteria

- [ ] `grep -rn "static func _qmul\|static func _qrot" spm/Sources` → no hits
- [ ] Exactly one definition site (`QuaternionKernels.swift`)
- [ ] Operator test results identical before/after
- [ ] `make build test` passes

## Out of scope

The broader Float/Double collapse (chunk 13); changing operator semantics.
