---
chunk: 03-spectral-tests
status: pending
depends_on: [02]
audit: ../review-for-fixes/2026-07-11-swift-audit.md §O1 (spectral)
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 03 — Spectral surface tests (PSD / spectrogram / mel / features)

**Deliverable:** behavioral pins for everything downstream of the restored
FFT. Test-only chunk — `src` edits only where a test exposes a real defect
(then fix minimally and note it).

## Files

- Create: `spm/Tests/spmMathToolsTests/Waveform1DSpectralTests.swift`

## Design constraints (the analytic assertions)

1. `powerSpectralDensity` of seeded white noise: flat within ±6 dB of the
   mean across the passband; total power ≈ signal variance (rtol 0.1).
2. `powerSpectralDensity` of a 50 Hz sine: single dominant peak at 50 Hz.
3. `spectrogram` of a linear chirp: per-column argmax frequency is
   monotonically non-decreasing.
4. `melSpectrogram`: filterbank output non-negative; energy concentrated
   in ascending mel bins for the chirp.
5. `extractSpectralFeatures` of the 50 Hz sine: spectral centroid ≈ 50 Hz
   (rtol 5e-2); flatness near 0 for the sine, near 1 for white noise.
6. Every test seeds its noise deterministically; sizes are powers of two
   (the non-pow2 path is chunk 02's test).

## TDD steps

1. Write all tests; run. Passing = pins. Failing = a real downstream
   defect (likely: scaling conventions) — fix minimally in the owning
   file, cite which test forced it.
2. `make build test` green.

## Acceptance criteria

- [ ] All six assertion families present and passing
- [ ] Any `src` diff is justified by a named failing test in the completion notes
- [ ] `make build test` passes

## Out of scope

Filtering/envelope/phase (chunk 10); API changes to the spectral surface.
