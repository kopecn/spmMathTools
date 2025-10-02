import Foundation

/// Spectrogram data structure
public struct WaveformSpectrogram<T> {
    /// Time frames (center time of each window)
    public let timeFrames: [TimeInterval]

    /// Frequency bins
    public let frequencies: [T]

    /// Spectrogram data [time][frequency]
    public let data: [[T]]

    /// Analysis parameters
    public let windowSize: Int
    public let hopSize: Int
    public let windowType: WaveformWindowType
    public let scaling: WaveformSpectrogramScaling
    public let samplingRate: Double
    public let originalDuration: TimeInterval

    /// Get spectrogram value at specific time and frequency indices
    /// - Parameters:
    ///   - timeIndex: Time frame index
    ///   - frequencyIndex: Frequency bin index
    /// - Returns: Spectrogram value or nil if indices are out of bounds
    public func value(at timeIndex: Int, frequencyIndex: Int) -> T? {
        guard timeIndex >= 0 && timeIndex < data.count && frequencyIndex >= 0 && frequencyIndex < frequencies.count
        else {
            return nil
        }
        return data[timeIndex][frequencyIndex]
    }

    /// Get time slice at specific frequency
    /// - Parameter frequencyIndex: Frequency bin index
    /// - Returns: Time series at the specified frequency
    public func timeSlice(at frequencyIndex: Int) -> [T]? {
        guard frequencyIndex >= 0 && frequencyIndex < frequencies.count else { return nil }
        return data.map { $0[frequencyIndex] }
    }

    /// Get frequency slice at specific time
    /// - Parameter timeIndex: Time frame index
    /// - Returns: Frequency spectrum at the specified time
    public func frequencySlice(at timeIndex: Int) -> [T]? {
        guard timeIndex >= 0 && timeIndex < data.count else { return nil }
        return data[timeIndex]
    }
}
