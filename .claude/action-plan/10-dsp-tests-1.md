---
chunk: 10-dsp-tests-1
status: pending
depends_on: [02]
audit: ../review-for-fixes/2026-07-11-swift-audit.md §O1 (filtering / envelope / phase)
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 10 — DSP tests, part 1: filtering, envelope, phase

**Deliverable:** analytic behavioral pins for three untested families.
Test-only; `src` edits only when a test exposes a defect (fix minimally,
cite the test).

## Files

- Create: `spm/Tests/spmMathToolsTests/Waveform1DFilteringTests.swift`,
  `Waveform1DEnvelopeTests.swift`, `Waveform1DPhaseTests.swift`

## Design constraints (assertions per family)

- **Filtering:** low-pass on a 5 Hz + 200 Hz mix (fs 2 kHz) attenuates
  200 Hz ≥ 20 dB while 5 Hz survives (rtol 0.1, interior); high-pass
  mirror; moving average of a constant is identity; Savitzky–Golay on a
  noiseless cubic reproduces it (atol 1e-6); Whittaker–Henderson λ→0 ≈
  identity.
- **Envelope:** `amplitudeEnvelope` of `A·sin` ≈ A interior (rtol 0.1);
  damped sinusoid envelope tracks `A·e^{−λt}` (rtol 0.15 interior); upper
  ≥ lower everywhere; each `WaveformInstantaneousMethod` returns full
  length.
- **Phase:** unwrapped instantaneous phase of a chirp is monotone;
  instantaneous frequency of an f-Hz sine ≈ f interior (rtol 5e-2); phase
  difference of sin vs cos ≈ π/2 interior; PLV(self) == 1, PLV(independent
  seeded noise) < 0.5.
- Tolerances are stated loose deliberately (these pin *behavior classes*,
  not implementations); tighten only with evidence.

## TDD steps

1. Write all tests; run. 2. Fix any genuine defects minimally.
3. `make build test` green.

## Acceptance criteria

- [ ] Every listed assertion exists and passes
- [ ] Any src diff justified by a named failing test
- [ ] `make build test` passes

## Out of scope

Chunk 11's families; spectral (03); API changes.
