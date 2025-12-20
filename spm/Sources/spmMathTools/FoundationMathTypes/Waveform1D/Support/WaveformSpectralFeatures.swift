import Foundation

/// Spectral features extracted from spectrogram
public struct WaveformSpectralFeatures<U: BinaryFloatingPoint & Sendable> {
    /// Spectral centroid for each time frame
    public let spectralCentroids: [U]

    /// Spectral rolloff for each time frame
    public let spectralRolloffs: [U]

    /// Spectral flux for each time frame
    public let spectralFluxes: [U]

    /// Corresponding time frames
    public let timeFrames: [U]
}
