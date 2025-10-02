/// Spectrogram scaling types
public enum WaveformSpectrogramScaling {
    /// Linear magnitude
    case magnitude
    /// Power (magnitude squared)
    case power
    /// Power spectral density
    case powerSpectralDensity
    /// Decibel (20*log10(magnitude))
    case decibel
    /// Complex representation (magnitude and phase)
    case complex
}
