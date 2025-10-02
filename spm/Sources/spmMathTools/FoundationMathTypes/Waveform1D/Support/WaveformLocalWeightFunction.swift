/// Local weight functions for polynomial regression
public enum WaveformLocalWeightFunction {
    /// Uniform weights (unweighted)
    case uniform
    /// Tricube weight function
    case tricube
    /// Gaussian weights with specified sigma
    case gaussian(sigma: Double)
    /// Epanechnikov weight function
    case epanechnikov
}
