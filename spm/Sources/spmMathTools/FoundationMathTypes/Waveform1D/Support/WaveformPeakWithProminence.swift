import Foundation

public struct WaveformPeakWithProminence<T: Numeric & Sendable, U: BinaryFloatingPoint & Sendable> {
    let peak: WaveformPeak<T,U>
    let prominence: T
}
