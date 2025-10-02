/// Frequency range specification
public struct WaveformFrequencyRange {
    public let minFreq: Double
    public let maxFreq: Double

    public init(minFreq: Double, maxFreq: Double) {
        self.minFreq = minFreq
        self.maxFreq = maxFreq
    }
}
