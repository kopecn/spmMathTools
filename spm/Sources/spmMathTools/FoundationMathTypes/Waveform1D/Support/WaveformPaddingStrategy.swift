/// Padding strategies for windowing
public enum WaveformPaddingStrategy {
    /// No padding, keep incomplete windows as-is
    case none
    /// Pad with zeros
    case zeros
    /// Pad with last value
    case lastValue
    /// Mirror padding
    case mirror
}
