---
chunk: 11-dsp-tests-2
status: pending
depends_on: []
audit: ../review-for-fixes/2026-07-11-swift-audit.md §O1 (resampling / windowing / triggers / zero-crossings / generators)
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 11 — DSP tests, part 2: resampling, windowing, triggers, zero-crossings, generators

**Deliverable:** analytic behavioral pins for the remaining untested
families. Test-only; minimal src fixes only when a test exposes a defect.

## Files

- Create: `spm/Tests/spmMathToolsTests/Waveform1DResamplingTests.swift`,
  `Waveform1DWindowingTests.swift`, `Waveform1DTriggerTests.swift`,
  `Waveform1DZeroCrossingTests.swift`, `Waveform1DGeneratorTests.swift`

## Design constraints (assertions per family)

- **Generators** (cheapest, do first — other families consume them):
  `sine` matches `sin(2πft)` on the time axis (atol 1e-12); `constant`,
  `impulse`, `linearRamp`, `heaviside` shape checks; `whiteNoise(seed:)`
  reproducible for equal seeds, different for different seeds; `counter`
  is 0..<n.
- **Resampling:** `decimated(by: 4)` → count n/4 and `dt` ×4;
  `interpolated` ×2 then `decimated` ×2 round-trips a smooth signal
  interior (rtol 1e-2); polyphase variants preserve a pure tone's
  frequency (argmax bin via chunk 02's fft).
- **Windowing:** rectangular `windowed` is identity; Hann endpoints ≈ 0,
  midpoint ≈ 1; `windowCoherentGain(.hann)` ≈ 0.5 (rtol 2e-2); every
  `WaveformWindowType` case generates without error.
- **Triggers:** square wave — rising-edge count == cycles; falling ==
  cycles; window enter/exit pair on a sine crossing a band; no-hit →
  empty.
- **Zero-crossings:** f-Hz sine over integer periods → count == 2·periods
  (±1); positive+negative == both; empty/constant waveform → empty.

## TDD steps

1. Generators first, then the other four files. 2. Minimal fixes with
   citations. 3. `make build test` green.

## Acceptance criteria

- [ ] Every listed assertion exists and passes
- [ ] Any src diff justified by a named failing test
- [ ] `make build test` passes

## Out of scope

TimeAlignment (blocked on `Waveform1D-TimeAlignment.swift:93` TODO — file a
known-issues.md entry instead); chunk 10's families.
