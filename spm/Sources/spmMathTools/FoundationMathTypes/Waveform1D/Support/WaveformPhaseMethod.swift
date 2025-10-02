/// Phase calculation methods
public enum WaveformPhaseMethod {
    /// Hilbert transform method
    case hilbert
    /// FFT-based method
    case fft
    /// Derivative-based method (for sinusoidal signals)
    case derivative
}
