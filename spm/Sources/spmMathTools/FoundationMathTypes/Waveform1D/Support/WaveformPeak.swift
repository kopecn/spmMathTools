import Foundation

public struct WaveformPeak<T> {
    /// Index in the values array where the peak occurs
    public let index: Int

    /// Value at the peak
    public let value: T

    /// Absolute time of the peak (if t0 is available)
    public let time: Date?

    /// Time offset from waveform start in seconds
    public let timeOffset: TimeInterval
}


