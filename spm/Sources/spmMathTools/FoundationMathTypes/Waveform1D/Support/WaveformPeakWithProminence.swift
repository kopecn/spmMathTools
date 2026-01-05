import Foundation

public struct WaveformPeakWithProminence<T: Numeric & Sendable> {
    let peak: WaveformPeak<T>
    let prominence: T
}
