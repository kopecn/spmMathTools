import Foundation
import FoundationTypes

/// Trigger configuration
public struct WaveformTrigger<T: Numeric & Sendable> {
    /// Type of trigger
    public let type: WaveformTriggerType<T>

    /// Minimum time interval between triggers
    public let minimumInterval: PrecisionTimeInterval?

    public init(type: WaveformTriggerType<T>, minimumInterval: PrecisionTimeInterval? = nil) {
        self.type = type
        self.minimumInterval = minimumInterval
    }
}
