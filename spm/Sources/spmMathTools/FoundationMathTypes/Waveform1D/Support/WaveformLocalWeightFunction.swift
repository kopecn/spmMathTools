/// Local weight functions for polynomial regression
public enum WaveformLocalWeightFunction<T: BinaryFloatingPoint> {
    /// Uniform weights (unweighted)
    case uniform
    /// Tricube weight function
    case tricube
    /// Gaussian weights with specified sigma
    case gaussian(sigma: T)
    /// Epanechnikov weight function
    case epanechnikov
}
