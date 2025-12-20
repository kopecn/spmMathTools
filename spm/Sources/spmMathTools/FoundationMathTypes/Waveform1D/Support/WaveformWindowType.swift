/// Available window function types
public enum WaveformWindowType<U: BinaryFloatingPoint & Sendable> {
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
    case kaiser(beta: U)
    /// Tukey window with taper ratio
    case tukey(taperRatio: U)
    /// Bartlett (triangular) window
    case bartlett
    /// Welch (parabolic) window
    case welch
}
