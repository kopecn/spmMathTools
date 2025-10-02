/// Envelope detection methods
public enum WaveformEnvelopeMethod {
    /// Hilbert transform-based envelope
    case hilbert
    /// Peak detection with sliding window
    case peakDetection(windowSize: Int)
    /// RMS envelope with sliding window
    case rms(windowSize: Int)
    /// Simple absolute value envelope
    case absoluteValue
}
