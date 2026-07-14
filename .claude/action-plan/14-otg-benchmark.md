---
chunk: 14-otg-benchmark
status: pending
depends_on: [04]
audit: ../review-for-fixes/2026-07-11-swift-audit.md §O4
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 14 — OTG per-cycle benchmark (measure first, optional)

**Deliverable:** the measurement that decides whether O4's allocation
concern is real. NO optimization in this chunk (decision framework:
performance work is measurement-driven).

## Files

- Edit/Create: `spm/Tests/spmMathToolsTests/OTGTests/OTGPerformanceTests.swift`
  (extend the existing 4 tests)
- Create: `.claude/review-for-fixes/otg-benchmark-results.md`

## Design constraints

1. Benchmark a representative control loop: 3-DOF position interface,
   1 kHz control cycle, 10k `update` calls — report wall time per cycle
   (p50/p99) via `ContinuousClock`, on `make build` (release) binaries,
   not debug.
2. Approximate allocation pressure: run the same loop under
   `-Xswiftc -sanitize=address` timing delta, or simply count
   `Trajectory.atTime` array creations via a temporary counter — state the
   method used.
3. The report ends with a go/no-go on the two candidate optimizations the
   audit named (inout sampling buffers into `OutputParameter`; struct-ify
   step solvers) — with a projected benefit each, or "not worth it".

## Acceptance criteria

- [ ] Benchmark runs reproducibly (two runs within 20% of each other)
- [ ] Results doc committed with p50/p99 per cycle and the go/no-go
- [ ] Zero src changes (`git diff --stat` touches only Tests/ and .claude/)
- [ ] `make build test` passes

## Out of scope

Implementing any optimization — that is a follow-up scope decision made on
this chunk's numbers.
