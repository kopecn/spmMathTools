import Foundation

/// Represents a zero crossing event

public struct WaveformZeroCrossing<T: Numeric, U: BinaryFloatingPoint> {
    /// Interpolated sample index where crossing occurs
    public let sampleIndex: U

    /// Absolute time of crossing (if t0 is available)
    public let time: Date?

    /// Time offset from waveform start in seconds
    public let timeOffset: U

    /// Type of crossing (rising or falling)
    public let type: WaveformZeroCrossingType

    /// Magnitude of the crossing (difference between adjacent samples)
    public let magnitude: T
}
