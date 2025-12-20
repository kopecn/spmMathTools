import Foundation

/// Time alignment methods
public enum WaveformAlignmentMethod<U: BinaryFloatingPoint & Sendable> {
    /// Use cross-correlation to find optimal alignment
    case crossCorrelation
    /// Align based on absolute timestamps (requires t0)
    case timeStamp
    /// Manual time offset in seconds
    case manualOffset(U)
    /// Align based on peak positions
    case peakAlignment
    /// Align based on energy envelopes
    case energyAlignment
}
