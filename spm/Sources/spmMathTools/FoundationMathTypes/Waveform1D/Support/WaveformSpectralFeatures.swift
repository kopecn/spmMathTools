import Foundation

/// Spectral features extracted from spectrogram
public struct WaveformSpectralFeatures<T> {
    /// Spectral centroid for each time frame
    public let spectralCentroids: [T]

    /// Spectral rolloff for each time frame
    public let spectralRolloffs: [T]

    /// Spectral flux for each time frame
    public let spectralFluxes: [T]

    /// Corresponding time frames
    public let timeFrames: [TimeInterval]
}
