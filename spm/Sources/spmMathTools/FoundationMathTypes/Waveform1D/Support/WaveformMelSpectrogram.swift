import Foundation

/// Mel-scale spectrogram data structure
public struct WaveformMelSpectrogram<U: BinaryFloatingPoint & Sendable> {
    /// Time frames
    public let timeFrames: [U]

    /// Mel-frequency bins
    public let melFrequencies: [U]

    /// Mel spectrogram data [time][mel_frequency]
    public let data: [[U]]

    /// Analysis parameters
    public let windowSize: Int
    public let hopSize: Int
    public let numMelBins: Int
    public let samplingRate: U
}
