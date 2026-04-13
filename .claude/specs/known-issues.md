# spmMathTools — Known Issues & TODOs

1. **SpatialPose+Arithmetic** should use primitives and avoid creating intermediary types (per README)
2. **Lexicographical Comparison** for Waveform types not yet implemented (per README)
3. **PrecisionTime resolution** — Extensions in `PrecisionTime/` are currently empty stubs. `PrecisionTimestamp` arithmetic (`.seconds`, `.attoseconds`, `.sign`, `attosecondsPerSecond`) may not be fully surfaced from `spmFoundationTools` yet. The `Waveform1D+subset.swift` file works around this using `.secondsAsDouble` and `addingTimeInterval(add:)`.
4. **Spatial operator duplication** — Double and Float overloads are copy-pasted; could collapse via `where T: FloatingPoint` generic constraints.
