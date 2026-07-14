---
type: audit
name: swift-audit-math-tools
purpose: Repo-wide Swift audit — findings ranked by severity with minimal fixes
last_updated: 2026-07-11
semver: 0.0.1
author: Nicholas Bergantz
---

# Swift Audit — spmMathTools

Scope: the `spmMathTools` target (~17,400 lines under
`FoundationMathTypes/` + `Extensions/`). Reviewed against the Swift specs
(`swift.md` + concurrency / error-handling / performance / cross-platform).
Each finding: severity → evidence → minimal fix. Paired audit:
`spmFoundationTools/.claude/review-for-fixes/2026-07-11-swift-audit.md`.

## Fix — BLOCKER

### F1. The FFT is non-functional — every spectral API returns garbage
`Waveform1D/Extensions/Waveform1D-Spectrogram+FFTHelpers.swift:17`: the
Cooley–Tukey **butterfly passes are commented out**
(`FIXME: … broken upstream`); `fft_cooleyTukey` performs only the
bit-reversal permutation. That helper is the engine for the public surface:

- `Waveform1D-FFT.swift:22` — `fft()`
- `Waveform1D-FFT.swift:119, 191` — `powerSpectralDensity` / single-segment PSD
- `Waveform1D-Spectrogram.swift:254` — `spectrogram` (and thus mel + features)

So `fft/PSD/spectrogram/melSpectrogram/spectralFeatures` all return
bit-reversed time samples, not a Fourier transform — and there are **zero
tests** on the spectral surface to catch it. Likely root cause: the
butterfly needs `Complex<T>` multiplication, but complex `*` exists only
for `T == Float` (see F6), so the generic code couldn't compile and was
stubbed out.
**Fix (in order):** (1) implement generic complex multiply (F6); (2)
restore the butterfly (the commented code is nearly correct — hoist the
twiddle recurrence out of the inner loop); (3) add analytic tests before
trusting anything: 10 Hz sine → single peak bin, Parseval energy check,
impulse → flat spectrum, PSD of white noise ≈ flat. Until then, consider
making the spectral methods return `nil`/throw so callers cannot consume
wrong numbers silently.

## Fix — correctness & error handling

### F2. 60 `print("[DEBUG] …")` statements in the OTG hot path — MAJOR
`OTG/CalculatorTarget.swift:217-218, 452, 766-769`,
`OTG/position/PositionThirdOrderStep2.swift:480-518`, and ~50 more. A
real-time trajectory library writing to stdout every failed sync branch is
both an error-handling violation (print ≠ logging) and a per-cycle
performance hit. **Fix:** delete them, or route through `swift-log` at
`.trace` behind a logger the caller injects (foundation already ships
`Logging`).

### F3. Force unwraps on runtime solver state — MAJOR
`OTG/CalculatorTarget.swift:208-210`:
`profiles[divRem] = blocks[divRem].a!.profile!` (and `.b!.profile!`). If the
block interval is absent at that branch, the whole process traps
(swift.md: no force unwraps outside tests/proven boundaries).
**Fix:** `guard let` and return `Result.errorSynchronizationCalculation` —
the enum exists precisely for this.

### F4. `fatalError`/`preconditionFailure` on recoverable library states — MAJOR
- `Waveform1D-Spectrogram+FFTHelpers.swift:11` — `fatalError("FFT input
  size must be a power of 2")` on runtime data reachable from public APIs.
- `OTG/Trajectory.swift:55, 67, 122` — `preconditionFailure("… Trajectory
  was not computed before use")` — misuse, but a host app dies for it.
- `OTG/OutputParameter.swift:100` — doc comment claims
  "`Throws: preconditionFailure`" — documentation of a crash as a throw.
**Fix:** FFT helper returns `nil`/throws (callers already return
optionals); Trajectory returns a typed error or optional state; fix the
OutputParameter doc to say what actually happens.

### F5. `min()!` traps on the empty-array path — MINOR
`Waveform1D/Extensions/Waveform1D-Envelope.swift:285`:
`distances.firstIndex(of: distances.min()!) ?? 0` — the `?? 0` guards the
wrong call; `min()!` itself traps when `distances` is empty.
**Fix:** `guard let minValue = distances.min(), let idx = …` early-return.

### F6. Complex arithmetic exists only for `T == Float` — MAJOR
`Complex/Complex+Arithmetic.swift:7` (`extension Complex where T == Float`).
`Complex<Double>` — the type every `DoubleWaveform1D` spectral path needs —
has no `* / + -`. This is the probable root cause of F1.
**Fix:** implement once generically
(`extension Complex where T: BinaryFloatingPoint & SIMDScalar`), delete the
Float-only copy.

### F7. Unresolved `FIXME` markers inside the core solver — MAJOR (investigate)
`OTG/position/PositionThirdOrderStep2.swift:1355, 1362` — bare
`// MARK: - FIXME` in root-selection branches of the largest solver file.
Unknown known-issues in the numeric core. **Fix:** characterize each with a
failing test (the OTGTruthTable/FailureFix suites are the harness) or
delete the marker with a comment saying why the branch is correct.

## Fix — packaging & structure

### F8. Branch pin `branch: "dev"` on spmFoundationTools — MAJOR
`Package.swift:19`. Builds are non-reproducible: any push to foundation's
`dev` changes what this package builds against (repo BKM: deployables pin
exact refs). **Fix:** tag foundation releases and pin
`from:`/`exact:` a version; keep branch pins only in local development via
a path-based override.

### F9. Dead scaffolding ships in the target — MINOR
- `Waveform1D/Operators/Waveform1D+Custom.swift` — 18 TODO placeholders,
  zero API.
- `SpatialWaveforms/Operators/Waveform{Position,Quaternion,SpatialPose}
  +Arithmetic.swift` — empty extension bodies.
- `PrecisionTime/Extensions/PrecisionTime{stamp,Interval}+Extensions.swift`
  — zero-byte files.
**Fix:** delete them all; git remembers the scaffolds when the features
arrive (dead-code baseline).

### F10. `depermaid` declared, used by no target — MINOR
`Package.swift:20` — same as the foundation repo. **Fix:** remove.

### F11. Platform mismatch with the dependency — MINOR
This package: `.macOS(.v14)` only; foundation supports iOS 16 / tvOS 16 /
watchOS 9. The math layer silently forbids every non-macOS consumer of the
foundation types. **Fix:** match foundation's platform list (nothing here
is macOS-specific), or record the macOS-only decision in `.claude/CLAUDE.md`.

## Fix — duplication (Define Once)

### F12. `_qmul` / `_qrot` implemented four times — MINOR
`Spatial/Operators/Quaternion+InlinableQuaternion.swift:13,41` and
`Spatial/Operators/SpatialPose+InlinableQuaternion.swift:12,30` carry
identical SIMD quaternion kernels (×2 files ×2 scalar types).
**Fix:** one internal `enum QuaternionKernels` (or free `@inlinable`
functions) both extensions call.

### F13. Float/Double operator file duplication — MINOR (structural)
Nearly every operator file is duplicated per scalar
(`where T == Double` + `where T == Float`). With F5's `sincos` helper on
the foundation side, most collapse to single
`where T: BinaryFloatingPoint & SIMDScalar` implementations. Large diff,
mechanical win — schedule as its own pass with the existing operator tests
as the harness.

## Optimize / test debt

### O1. The untested majority of Waveform1D — MAJOR (test debt)
Tested: calc, correlation, peaks (~105 tests). Untested: **FFT, PSD,
spectrogram (see F1), filtering, envelope, phase, resampling, windowing,
triggers, zero-crossings, generators** — public API with no behavioral
pins. **Fix:** analytic-signal tests per family (the py-MathTools port
plan's chunk criteria in
`py-MathTools/.claude/action-plan/18–29` enumerate exactly which analytic
assertions to use — reuse them here).

### O2. Truth-table JSONs skip silently — MINOR
`OTGTests/OTGComprehensiveTests.swift:406,447` — `XCTSkip("… Run
generateRandomTrajectories.py first.")` when the JSONs are missing. The
files ARE committed under `OTGTests/truthTables/`, so the skip masks a
resource-bundling failure if the path breaks. Also the JSONs record inputs
only; the numeric golden lives solely in the 32-case hardcoded array in
`OTGTruthTableTests.swift`. **Fix:** make missing-resource a failure, not a
skip; consider exporting the 32-case array to a committed JSON so other
ports (py-MathTools) can consume it directly.

### O3. Mixed XCTest (12) / swift-testing (4) — MINOR
Chamber Match: converge on swift-testing for new tests.

### O4. OTG per-cycle allocations — MINOR (measure first)
`Trajectory.atTime` returns three fresh `[Double]` per call inside the
control loop; step-solver classes are reference types churned per
calculation. Do not refactor on speculation (decision framework:
performance work is measurement-driven) — add a benchmark first
(`OTGPerformanceTests` exists; extend it with an allocations counter),
then consider `inout` buffers writing straight into `OutputParameter`.

### O5. TODO backlog worth triaging — MINOR
`Waveform1D-TimeAlignment.swift:93` (trimToCommonTimeBase blocked on
PrecisionTimestamp support), `Waveform1D-PhaseAnalysis.swift:211`
("Requires FFT extension — fallback to hilbert" — F1 again),
`Waveform1D-Spectrogram.swift:369` (spectral features want PrecisionTime
types). Each is a small design decision; file them or do them.

## Suggested order

F1+F6 (the spectral surface is silently wrong today) → F2/F3/F4 (library
crash/noise discipline) → F8 (reproducible builds) → O1 (test the fixed
spectral code first) → F9/F10/F12 → F13/O4 as measured, mechanical passes.
