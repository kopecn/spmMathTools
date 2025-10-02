import Foundation

/// Waveform with event markers
public struct WaveformWithEvents<T: Numeric> {
    /// The original waveform
    public let waveform: Waveform1D<T>

    /// Array of event markers
    public let events: [WaveformEventMarker]

    /// Get events within a time range
    /// - Parameters:
    ///   - start: Start time
    ///   - end: End time
    /// - Returns: Filtered events
    public func events(from start: TimeInterval, to end: TimeInterval) -> [WaveformEventMarker] {
        return events.filter { event in
            event.timeOffset >= start && event.timeOffset <= end
        }
    }

    /// Get events of a specific annotation type
    /// - Parameter annotation: Annotation type to filter by
    /// - Returns: Filtered events
    public func events(ofType annotation: WaveformEventAnnotation) -> [WaveformEventMarker] {
        return events.filter { event in
            switch (event.annotation, annotation) {
            case (.marker, .marker), (.peak, .peak), (.valley, .valley):
                return true
            case (.custom(let a), .custom(let b)):
                return a == b
            default:
                return false
            }
        }
    }
}
