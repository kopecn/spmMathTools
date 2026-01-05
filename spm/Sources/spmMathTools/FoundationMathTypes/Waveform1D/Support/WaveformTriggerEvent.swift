import Foundation
import FoundationTypes

/// Trigger event information
public struct WaveformTriggerEvent<T: Numeric & Sendable> {
    /// Sample index where trigger occurred
    public let index: Int

    /// Value at trigger point
    public let value: T

    /// Absolute time of trigger (if t0 available)
    public let time: Date?

    /// Time offset from waveform start
    public let timeOffset: PrecisionTimeInterval

    /// Type of trigger that fired
    public let type: WaveformTriggerType<T>
}
