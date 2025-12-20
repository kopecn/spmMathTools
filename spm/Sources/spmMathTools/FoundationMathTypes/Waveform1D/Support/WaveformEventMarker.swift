import Foundation

/// Event marker for visualization
public struct WaveformEventMarker<U: BinaryFloatingPoint & Sendable> {
    /// Sample index of the event
    public let index: Int

    /// Absolute time of event (if available)
    public let time: Date?

    /// Time offset from waveform start
    public let timeOffset: U

    /// Annotation type
    public let annotation: WaveformEventAnnotation

    /// Additional metadata
    public let metadata: [String: String]
}
