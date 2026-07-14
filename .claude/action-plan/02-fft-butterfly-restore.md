---
chunk: 02-fft-butterfly-restore
status: pending
depends_on: [01]
audit: ../review-for-fixes/2026-07-11-swift-audit.md §F1, §F4(fft)
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# 02 — Restore the FFT butterfly

**Deliverable:** a working Cooley–Tukey FFT behind `fft()`, PSD, and
spectrogram — plus the `fatalError` removed from the helper.

## Files

- Edit: `spm/Sources/spmMathTools/FoundationMathTypes/Waveform1D/Extensions/Waveform1D-Spectrogram+FFTHelpers.swift`
- Edit (only if the helper's signature change forces it):
  `Waveform1D-FFT.swift`, `Waveform1D-Spectrogram.swift`
- Create: `spm/Tests/spmMathToolsTests/Waveform1DFFTTests.swift`

## Recipe

The commented-out butterfly (`FFTHelpers.swift:21-37`) is structurally
right; restore it with the twiddle recurrence hoisted (avoids `cos`/`sin`
per inner iteration):

```swift
var length = 2
while length <= n {
    let half = length / 2
    let angleStep = -2.0 * Double.pi / Double(length)
    let wStep = Complex<T>(real: T(cos(angleStep)), imaginary: T(sin(angleStep)))
    for start in stride(from: 0, to: n, by: length) {
        var w = Complex<T>(real: 1, imaginary: 0)
        for j in 0..<half {
            let u = x[start + j]
            let v = w * x[start + j + half]      // needs chunk 01
            x[start + j] = u + v
            x[start + j + half] = u - v
            w = w * wStep
        }
    }
    length *= 2
}
```

`fatalError("FFT input size must be a power of 2")` (`:11`) becomes a
`guard … else { return nil }` / thrown error consistent with what the
public callers already do on bad input (`fft()` returns an Optional —
propagate nil). Delete the FIXME comment block and the now-obsolete entry
in `.claude/specs/known-issues.md` if one exists.

## TDD steps

1. Failing tests FIRST (they fail loudly against the broken helper):
   (a) `fft()` of a 10 Hz sine (fs 1 kHz, n 1024) → magnitude argmax at
   the 10 Hz bin and > 100× the median bin; (b) Parseval:
   `Σ|X|²/n == Σx²` (rtol 1e-9); (c) impulse → flat magnitude spectrum
   (max/min < 1+1e-9); (d) linearity: `fft(a+b) == fft(a)+fft(b)`
   elementwise (atol 1e-9); (e) non-power-of-2 input → nil/throws, no
   crash.
2. Restore the butterfly per the recipe. 3. `make build test` green.

## Acceptance criteria

- [ ] All five FFT tests pass; the sine-peak test FAILED before the fix
      (record the before-run in completion notes — proves the oracle bites)
- [ ] `grep -n "fatalError" spm/Sources/spmMathTools/FoundationMathTypes/Waveform1D/` → no hits
- [ ] `grep -n "FIXME" …/Waveform1D-Spectrogram+FFTHelpers.swift` → no hits
- [ ] `make build test` passes

## Out of scope

PSD/spectrogram-level tests (chunk 03); performance beyond the twiddle
recurrence; the phase-analysis TODO fallbacks.
