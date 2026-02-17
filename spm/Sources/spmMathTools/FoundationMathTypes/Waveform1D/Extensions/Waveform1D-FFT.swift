import Foundation
import FoundationTypes

// MARK: - FFT Mathematical Operations
extension Waveform1D where T: BinaryFloatingPoint & SIMDScalar {

    /// Compute the Fast Fourier Transform using Cooley-Tukey algorithm
    /// - Returns: Tuple containing frequency array, magnitudes, and phases
    public func fft() -> (frequencies: [T], magnitudes: [T], phases: [T])? {
        guard values.count > 1 else { return nil }

        // Convert to complex numbers
        var complex = values.map { Complex<T>(real: $0, imaginary: T.zero) }

        // Pad to next power of 2 for efficiency
        let fftSize = nextPowerOfTwo(complex.count)
        while complex.count < fftSize {
            complex.append(Complex<T>(real: T.zero, imaginary: T.zero))
        }

        // Perform FFT
        fft_cooleyTukey(&complex)

        // Calculate frequencies (only positive frequencies)
        let samplingRate = 1.0 / dt.secondsAsDouble
        let nyquistBin = fftSize / 2
        let frequencies = (0..<nyquistBin).map { i in
            T(Double(i) * samplingRate / Double(fftSize))
        }

        // Calculate magnitudes and phases
        let magnitudes = (0..<nyquistBin).map { i in
            complexMagnitude(complex[i])
        }

        let phases = (0..<nyquistBin).map { i in
            complexPhase(complex[i])
        }

        return (frequencies, magnitudes, phases)
    }

    // MARK: - Private FFT Implementation

    private func nextPowerOfTwo(_ n: Int) -> Int {
        guard n > 1 else { return 1 }
        return 1 << Int(ceil(log2(Double(n))))
    }

    private func bitReversePermutation(_ x: inout [Complex<T>]) {
        let n = x.count
        var j = 0

        for i in 1..<n {
            var bit = n >> 1
            while j & bit != 0 {
                j ^= bit
                bit >>= 1
            }
            j ^= bit

            if i < j {
                x.swapAt(i, j)
            }
        }
    }
}

// MARK: - Power Spectral Density Operations
extension Waveform1D where T: BinaryFloatingPoint & SIMDScalar {

    /// Compute the Power Spectral Density using Welch's method
    /// - Parameters:
    ///   - windowType: Window function to apply (default: Hanning)
    ///   - windowSize: Size of each segment (default: signal length / 8)
    ///   - overlap: Overlap between segments as fraction (0.0 to 1.0, default: 0.5)
    ///   - scaling: PSD scaling type (default: density)
    /// - Returns: Tuple containing frequency array and PSD values
    public func powerSpectralDensity(
        windowType: WaveformWindowType = .hanning,
        windowSize: Int? = nil,
        overlap: T = 0.5,
        scaling: WaveformPSDScaling = .density
    ) -> (frequencies: [T], psd: [T])? {

        guard !values.isEmpty else { return nil }

        let actualWindowSize = windowSize ?? max(256, values.count / 8)
        let clampedOverlap = max(T(0.0), min(T(0.99), overlap))
        let hopSize = Int(T(actualWindowSize) * (T(1.0) - clampedOverlap))

        guard actualWindowSize <= values.count && hopSize > 0 else {
            // Fall back to single segment
            return singleSegmentPSD(windowType: windowType, scaling: scaling)
        }

        // Generate window using Windowing extension
        let window = Self.generateWindow(type: windowType, length: actualWindowSize)
        let windowPower = window.reduce(T(0.0)) { $0 + $1 * $1 }

        var psdAccumulator: [T] = []
        var segmentCount = 0

        // Process overlapping segments
        var startIndex = 0
        while startIndex + actualWindowSize <= values.count {
            let segment = Array(values[startIndex..<(startIndex + actualWindowSize)])

            // Apply window
            let windowedSegment = zip(segment, window).map { $0 * $1 }

            // Convert to complex and perform FFT
            var complex = windowedSegment.map { Complex<T>(real: $0, imaginary: T.zero) }
            let fftSize = nextPowerOfTwo(complex.count)
            while complex.count < fftSize {
                complex.append(Complex<T>(real: T.zero, imaginary: T.zero))
            }

            fft_cooleyTukey(&complex)

            // Calculate power spectrum for this segment
            let nyquistBin = fftSize / 2
            let segmentPSD = (0..<nyquistBin).map { i -> T in
                let magnitude = complexMagnitude(complex[i])
                return magnitude * magnitude
            }

            // Accumulate PSD
            if psdAccumulator.isEmpty {
                psdAccumulator = segmentPSD
            } else {
                for i in 0..<min(psdAccumulator.count, segmentPSD.count) {
                    psdAccumulator[i] += segmentPSD[i]
                }
            }

            segmentCount += 1
            startIndex += hopSize
        }

        guard segmentCount > 0 && !psdAccumulator.isEmpty else { return nil }

        // Average across segments
        let averagedPSD = psdAccumulator.map { $0 / T(segmentCount) }

        // Apply scaling
        let samplingRate = T(1.0 / dt.secondsAsDouble)
        let scaledPSD = applyWaveformPSDScaling(
            averagedPSD,
            samplingRate: samplingRate,
            windowPower: windowPower,
            windowSize: actualWindowSize,
            scaling: scaling
        )

        // Generate frequency array
        let samplingRateDouble = 1.0 / dt.secondsAsDouble
        let frequencies = (0..<scaledPSD.count).map { i in
            T(Double(i) * samplingRateDouble / Double(actualWindowSize))
        }

        return (frequencies, scaledPSD)
    }

    /// Compute PSD for a single segment (no windowing/averaging)
    /// - Parameters:
    ///   - windowType: Window function to apply
    ///   - scaling: PSD scaling type
    /// - Returns: Tuple containing frequency array and PSD values
    public func singleSegmentPSD(
        windowType: WaveformWindowType = .hanning,
        scaling: WaveformPSDScaling = .density
    ) -> (frequencies: [T], psd: [T])? {

        guard !values.isEmpty else { return nil }

        // Generate window using Windowing extension
        let window = Self.generateWindow(type: windowType, length: values.count)
        let windowPower = window.reduce(T(0.0)) { $0 + $1 * $1 }

        // Apply window
        let windowedValues = zip(values, window).map { $0 * $1 }

        // Convert to complex and perform FFT
        var complex = windowedValues.map { Complex<T>(real: $0, imaginary: T.zero) }
        let fftSize = nextPowerOfTwo(complex.count)
        while complex.count < fftSize {
            complex.append(Complex<T>(real: T.zero, imaginary: T.zero))
        }

        fft_cooleyTukey(&complex)

        // Calculate power spectrum
        let nyquistBin = fftSize / 2
        let powerSpectrum = (0..<nyquistBin).map { i -> T in
            let magnitude = complexMagnitude(complex[i])
            return magnitude * magnitude
        }

        // Apply scaling
        let samplingRate = T(1.0 / dt.secondsAsDouble)
        let scaledPSD = applyWaveformPSDScaling(
            powerSpectrum,
            samplingRate: samplingRate,
            windowPower: windowPower,
            windowSize: values.count,
            scaling: scaling
        )

        // Generate frequency array
        let samplingRateDouble = 1.0 / dt.secondsAsDouble
        let frequencies = (0..<scaledPSD.count).map { i in
            T(Double(i) * samplingRateDouble / Double(fftSize))
        }

        return (frequencies, scaledPSD)
    }

    /// Apply appropriate scaling to PSD values
    private func applyWaveformPSDScaling(
        _ powerSpectrum: [T],
        samplingRate: T,
        windowPower: T,
        windowSize: Int,
        scaling: WaveformPSDScaling
    ) -> [T] {

        switch scaling {
        case .density:
            // Power spectral density: units^2 / Hz
            let scale = T(1.0) / (samplingRate * windowPower)
            return powerSpectrum.map { $0 * scale }

        case .spectrum:
            // Power spectrum: units^2
            let scale = T(1.0) / (windowPower * windowPower)
            return powerSpectrum.map { $0 * scale }
        }
    }
}
