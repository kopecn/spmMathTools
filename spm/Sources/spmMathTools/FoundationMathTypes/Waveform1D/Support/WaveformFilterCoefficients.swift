/// Filter coefficients structure
public struct WaveformFilterCoefficients<T: Numeric & Sendable> {
    let b: [T]  // Numerator coefficients (feedforward)
    let a: [T]  // Denominator coefficients (feedback)
}
