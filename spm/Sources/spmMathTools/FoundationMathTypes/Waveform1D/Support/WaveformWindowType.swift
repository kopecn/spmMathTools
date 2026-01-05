import FoundationTypes

/// Available window function types
public enum WaveformWindowType {
    /// Rectangular window (no windowing)
    case rectangular
    /// Hanning window (raised cosine)
    case hanning
    /// Hamming window
    case hamming
    /// Blackman window
    case blackman
    /// Blackman-Harris window
    case blackmanHarris
    /// Kaiser window with beta parameter
    case kaiser(beta: PrecisionTimeInterval)
    /// Tukey window with taper ratio
    case tukey(taperRatio: PrecisionTimeInterval)
    /// Bartlett (triangular) window
    case bartlett
    /// Welch (parabolic) window
    case welch
}
