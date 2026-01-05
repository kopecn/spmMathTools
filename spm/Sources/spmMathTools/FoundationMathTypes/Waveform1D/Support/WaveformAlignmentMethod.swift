import Foundation
import FoundationTypes

/// Time alignment methods
public enum WaveformAlignmentMethod {
    /// Use cross-correlation to find optimal alignment
    case crossCorrelation
    /// Align based on absolute timestamps (requires t0)
    case timeStamp
    /// Manual time offset in seconds
    case manualOffset(PrecisionTimeInterval)
    /// Align based on peak positions
    case peakAlignment
    /// Align based on energy envelopes
    case energyAlignment
}
