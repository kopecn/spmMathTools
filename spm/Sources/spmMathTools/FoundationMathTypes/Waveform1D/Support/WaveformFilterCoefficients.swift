/// Filter coefficients structure
public struct WaveformFilterCoefficients<T> {
    let b: [T]  // Numerator coefficients (feedforward)
    let a: [T]  // Denominator coefficients (feedback)
}
