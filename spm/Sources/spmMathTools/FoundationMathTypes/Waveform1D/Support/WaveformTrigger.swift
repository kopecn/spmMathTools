import Foundation

/// Trigger configuration
public struct WaveformTrigger<T> {
    /// Type of trigger
    public let type: WaveformTriggerType<T>

    /// Minimum time interval between triggers
    public let minimumInterval: TimeInterval?

    public init(type: WaveformTriggerType<T>, minimumInterval: TimeInterval? = nil) {
        self.type = type
        self.minimumInterval = minimumInterval
    }
}
