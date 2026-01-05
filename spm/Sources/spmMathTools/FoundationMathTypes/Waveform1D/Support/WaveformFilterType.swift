import FoundationTypes

/// Digital filter types
public enum WaveformFilterType {
    /// Low-pass filter with cutoff frequency (represented by the interval timestep e.g. 1/frequency)
    case lowPass(cutoffFrequency: PrecisionTimeInterval)
    /// High-pass filter with cutoff frequency (represented by the interval timestep e.g. 1/frequency)
    case highPass(cutoffFrequency: PrecisionTimeInterval)
    /// Band-pass filter with low and high frequencies (represented by the interval timestep e.g. 1/frequency)
    case bandPass(lowFrequency: PrecisionTimeInterval, highFrequency: PrecisionTimeInterval)
    /// Band-stop (notch) filter with low and high frequencies (represented by the interval timestep e.g. 1/frequency)
    case bandStop(lowFrequency: PrecisionTimeInterval, highFrequency: PrecisionTimeInterval)
    /// Simple moving average filter
    case movingAverage(windowSize: Int)
    /// Exponential smoothing filter (represented by the interval timestep e.g. 1/frequency)
    case exponential(alpha: PrecisionTimeInterval)
}
