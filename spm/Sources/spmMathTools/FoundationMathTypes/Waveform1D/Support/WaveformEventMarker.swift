import Foundation

/// Event marker for visualization
public struct WaveformEventMarker {
    /// Sample index of the event
    public let index: Int

    /// Absolute time of event (if available)
    public let time: Date?

    /// Time offset from waveform start
    public let timeOffset: TimeInterval

    /// Annotation type
    public let annotation: WaveformEventAnnotation

    /// Additional metadata
    public let metadata: [String: String]
}
