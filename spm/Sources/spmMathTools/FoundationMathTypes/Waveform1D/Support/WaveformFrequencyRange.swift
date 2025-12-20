/// Frequency range specification
public struct WaveformFrequencyRange<T: BinaryFloatingPoint> {
    public let minFreq: T
    public let maxFreq: T

    public init(minFreq: T, maxFreq: T) {
        self.minFreq = minFreq
        self.maxFreq = maxFreq
    }
}
