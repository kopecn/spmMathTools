import Foundation

/// Mel-scale spectrogram data structure
public struct WaveformMelSpectrogram<T> {
    /// Time frames
    public let timeFrames: [TimeInterval]

    /// Mel-frequency bins
    public let melFrequencies: [T]

    /// Mel spectrogram data [time][mel_frequency]
    public let data: [[T]]

    /// Analysis parameters
    public let windowSize: Int
    public let hopSize: Int
    public let numMelBins: Int
    public let samplingRate: Double
}
