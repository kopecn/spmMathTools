
/// Correlation computation modes
public enum WaveformCorrelationMode {
    /// Full correlation: output size = len(x) + len(y) - 1
    case full
    /// Valid correlation: output only where signals fully overlap
    case valid
    /// Same correlation: output same size as first input
    case same
}