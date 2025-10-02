    /// Digital filter types
    public enum WaveformFilterType {
        /// Low-pass filter with cutoff frequency
        case lowPass(cutoffFrequency: Double)
        /// High-pass filter with cutoff frequency
        case highPass(cutoffFrequency: Double)
        /// Band-pass filter with low and high frequencies
        case bandPass(lowFrequency: Double, highFrequency: Double)
        /// Band-stop (notch) filter with low and high frequencies
        case bandStop(lowFrequency: Double, highFrequency: Double)
        /// Simple moving average filter
        case movingAverage(windowSize: Int)
        /// Exponential smoothing filter
        case exponential(alpha: Double)
    }