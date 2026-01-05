import Foundation
import FoundationTypes

/// Represents a zero crossing event

public struct WaveformZeroCrossing<T: Numeric> {
    /// Interpolated sample index where crossing occurs
    public let sampleIndex: Int

    /// Absolute time of crossing (if t0 is available)
    public let time: PrecisionTimestamp?

    /// Time offset from waveform start in seconds
    public let timeOffset: PrecisionTimeInterval

    /// Type of crossing (rising or falling)
    public let type: WaveformZeroCrossingType

    /// Magnitude of the crossing (difference between adjacent samples)
    public let magnitude: T
}
