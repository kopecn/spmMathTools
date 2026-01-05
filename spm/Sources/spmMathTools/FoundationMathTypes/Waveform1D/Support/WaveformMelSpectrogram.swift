import Foundation
import FoundationTypes

/// Mel-scale spectrogram data structure
public struct WaveformMelSpectrogram {
    /// Time frames
    public let timeFrames: [PrecisionTimestamp]

    /// Mel-frequency bins
    public let melFrequencies: [Double]

    /// Mel spectrogram data [time][mel_frequency]
    public let data: [[Double]]

    /// Analysis parameters
    public let windowSize: Int
    public let hopSize: Int
    public let numMelBins: Int
    public let samplingRate: PrecisionTimeInterval
}
