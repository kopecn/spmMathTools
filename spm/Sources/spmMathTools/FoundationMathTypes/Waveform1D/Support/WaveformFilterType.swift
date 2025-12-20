/// Digital filter types
public enum WaveformFilterType<U: BinaryFloatingPoint & Sendable> {
    /// Low-pass filter with cutoff frequency
    case lowPass(cutoffFrequency: U)
    /// High-pass filter with cutoff frequency
    case highPass(cutoffFrequency: U)
    /// Band-pass filter with low and high frequencies
    case bandPass(lowFrequency: U, highFrequency: U)
    /// Band-stop (notch) filter with low and high frequencies
    case bandStop(lowFrequency: U, highFrequency: U)
    /// Simple moving average filter
    case movingAverage(windowSize: Int)
    /// Exponential smoothing filter
    case exponential(alpha: U)
}
