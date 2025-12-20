import Foundation

/// Trigger configuration
public struct WaveformTrigger<T: Numeric & Sendable, U: BinaryFloatingPoint & Sendable> {
    /// Type of trigger
    public let type: WaveformTriggerType<T>

    /// Minimum time interval between triggers
    public let minimumInterval: U?

    public init(type: WaveformTriggerType<T>, minimumInterval: U? = nil) {
        self.type = type
        self.minimumInterval = minimumInterval
    }
}
