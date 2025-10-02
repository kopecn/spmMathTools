/// Interpolation methods for resampling operations
public enum WaveformInterpolationMethod {
    /// Linear interpolation between points
    case linear
    /// Cubic interpolation (Hermite)
    case cubic
    /// Nearest neighbor (no interpolation)
    case nearestNeighbor
}
