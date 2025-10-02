import Foundation

/// Instantaneous frequency analysis result
public struct WaveformInstantaneousFrequency<T> {
    /// Time frames
    public let timeFrames: [TimeInterval]

    /// Base frequency bins
    public let frequencies: [T]

    /// Instantaneous frequency data [time][frequency]
    public let instantaneousFrequencies: [[T]]
}
