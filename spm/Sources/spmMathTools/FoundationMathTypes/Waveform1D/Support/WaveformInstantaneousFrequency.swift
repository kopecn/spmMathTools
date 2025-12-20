import Foundation

/// Instantaneous frequency analysis result
public struct WaveformInstantaneousFrequency<T: BinaryFloatingPoint> {
    /// Time frames
    public let timeFrames: [T]

    /// Base frequency bins
    public let frequencies: [T]

    /// Instantaneous frequency data [time][frequency]
    public let instantaneousFrequencies: [[T]]
}
