import Foundation

/// Represents time lag information between two waveforms
public struct WaveformTimeLag<T: BinaryFloatingPoint> {
    /// Lag in number of samples
    public let lagSamples: Int

    /// Lag in time (seconds)
    public let lagTime: T

    /// Cross-correlation value at the lag
    public let correlation: T

    /// Confidence in the lag estimate (0.0 to 1.0)
    public let confidence: T
}
