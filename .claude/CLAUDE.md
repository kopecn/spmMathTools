# spmMathTools

## Project Overview

A **greenfield** Swift math and signal-processing library intended as a clean-room implementation to be **backported into Python** once stabilized. The library targets macOS 14+ with Swift 6.1 (strict concurrency) and ships as a single SPM library product.

> **Status: Early / Alpha** — The README warns: _"Dragons be here. Do not use this package yet. I have only validated 5% of the AI slop it spit out."_ Treat all numerics as unverified until covered by truth-table or property-based tests.

## Architecture

```
spmMathTools (library)
├── Extensions/            # Date, simd_float4x4, simd_double4x4 (Denavit-Hartenberg)
└── FoundationMathTypes/
    ├── Complex/           # SIMD-accelerated complex arithmetic (Float only)
    ├── Errors/            # MathErrors enum
    ├── Functional/        # Polynomial (zeroth–decic), root solvers, numeric utils
    ├── OTG/               # Optimal Trajectory Generation (Ruckig port, ~3 300 LOC)
    │   ├── enums/         # ControlInterface, Result, Synchronization, etc.
    │   ├── extensions/    # InputParameter+codable
    │   ├── position/      # 1st/2nd/3rd-order position step solvers
    │   ├── velocity/      # 2nd/3rd-order velocity step solvers
    │   └── support/
    ├── PrecisionTime/     # Extensions on FoundationTypes.PrecisionTimestamp/Interval
    ├── Spatial/           # SpatialPose, Position, Quaternion operators & constructors
    ├── SpatialWaveforms/  # Operators on WaveformPosition, WaveformQuaternion, WaveformSpatialPose
    ├── Waveform1D/        # 1-D signal type: operators, generators, calc, 40+ support types
    ├── WaveformPosition/  # Position-based waveform mutators
    ├── WaveformQuaternion/# Quaternion-based waveform mutators
    └── WaveformSpatialPose/ # SpatialPose-based waveform mutators
```

## Key Dependencies

| Package | Import name | Role |
|---------|------------|------|
| `spmFoundationTools` (dev branch) | `FoundationTypes` | Core types: `Waveform1D`, `Position`, `Quaternion`, `SpatialPose`, `PrecisionTimestamp`, `PrecisionTimeInterval`, `Complex` |
| `kvSIMD.swift` | `kvSIMD` | Cross-platform SIMD math (`atan2`, trig on vectors) |
| `depermaid` | _(plugin only)_ | Dependency graph visualization |

**Important:** The foundational _types_ live in `spmFoundationTools`. This repo provides **math operations, operators, and algorithms** on those types.

## Build & Test

```bash
# Build
swift build

# Test (includes 262 144 OTG truth-table trajectories)
swift test

# Run a single test class
swift test --filter spmMathToolsTests.OTGTruthTableTests
```

Tests live under `spm/Tests/spmMathToolsTests/` and mirror the source layout.

## Subspecs

Consult these when the task requires deeper context:

| Spec | When to read |
|------|-------------|
| [`specs/conventions.md`](specs/conventions.md) | Writing new files, adding types, naming questions, OTG or Polynomial module work |
| [`specs/known-issues.md`](specs/known-issues.md) | Debugging, adding features near PrecisionTime, Spatial, or Waveform types |
| [`specs/solid-assessment.md`](specs/solid-assessment.md) | Refactoring, architecture decisions, OTG restructuring |
| [`specs/python-backport.md`](specs/python-backport.md) | Any Python backport work |
