import Foundation
import FoundationTypes

/// Spectral features extracted from spectrogram
public struct WaveformSpectralFeatures {
    /// Spectral centroid for each time frame
    public let spectralCentroids: [PrecisionTimestamp]

    /// Spectral rolloff for each time frame
    public let spectralRolloffs: [PrecisionTimeInterval]

    /// Spectral flux for each time frame
    public let spectralFluxes: [Double]

    /// Corresponding time frames
    public let timeFrames: [PrecisionTimestamp]
}
