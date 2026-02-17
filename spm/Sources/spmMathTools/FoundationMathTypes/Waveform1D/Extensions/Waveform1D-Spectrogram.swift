import Foundation
import FoundationTypes

// MARK: - Spectrogram and Time-Frequency Analysis
extension Waveform1D where T: BinaryFloatingPoint & SIMDScalar {

    /// Compute the spectrogram using Short-Time Fourier Transform (STFT)
    /// - Parameters:
    ///   - windowSize: Size of the analysis window (default: 256)
    ///   - hopSize: Number of samples to advance between windows (default: windowSize/4)
    ///   - windowType: Window function to apply (default: hanning)
    ///   - scaling: Spectrogram scaling type (default: magnitude)
    ///   - frequencyRange: Frequency range to include (default: full range)
    /// - Returns: Spectrogram data structure
    public func spectrogram(
        windowSize: Int = 256,
        hopSize: Int? = nil,
        windowType: WaveformWindowType = .hanning,
        scaling: WaveformSpectrogramScaling = .magnitude,
        frequencyRange: WaveformFrequencyRange<T>? = nil
    ) -> WaveformSpectrogram<T>? {

        guard !values.isEmpty && windowSize > 0 && windowSize <= values.count else { return nil }

        let actualHopSize = hopSize ?? (windowSize / 4)
        guard actualHopSize > 0 else { return nil }

        // Calculate number of time frames
        let numFrames = max(1, (values.count - windowSize) / actualHopSize + 1)

        // Generate window coefficients
        let window = Self.generateWindow(type: windowType, length: windowSize)

        // Calculate frequency bins
        let frequencyBins = calculateFrequencyBins(windowSize: windowSize)
        let (startBin, endBin) = getFrequencyBinRange(frequencyBins: frequencyBins, range: frequencyRange)
        let selectedFrequencies = Array(frequencyBins[startBin..<endBin]).map { Double($0) }

        // Initialize spectrogram matrix
        var spectrogramData: [[T]] = Array(
            repeating: Array(repeating: T.zero, count: selectedFrequencies.count),
            count: numFrames
        )

        // Calculate time frames as PrecisionTimestamps
        var timeFrames: [PrecisionTimestamp] = []

        // Process each time frame
        for frameIndex in 0..<numFrames {
            let startSample = frameIndex * actualHopSize

            // Calculate time for this frame (center of window)
            let centerSample = startSample + windowSize / 2
            if let t0 = t0 {
                timeFrames.append(t0 + dt * centerSample)
            } else {
                // Use epoch-based timestamp if no t0
                let frameInterval = dt * centerSample
                timeFrames.append(PrecisionTimestamp(interval: frameInterval))
            }

            let endSample = min(startSample + windowSize, values.count)

            // Extract and pad segment if necessary
            var segment = Array(values[startSample..<endSample])
            if segment.count < windowSize {
                let paddingNeeded = windowSize - segment.count
                segment.append(contentsOf: Array(repeating: T.zero, count: paddingNeeded))
            }

            // Apply window function
            let windowedSegment = zip(segment, window).map { $0 * $1 }

            // Perform FFT
            if let fftResult = performSTFT(windowedSegment) {
                // Extract the desired frequency range and apply scaling
                for (binIndex, freqIndex) in (startBin..<endBin).enumerated() {
                    let magnitude = fftResult.magnitudes[freqIndex]
                    let phase = fftResult.phases[freqIndex]

                    spectrogramData[frameIndex][binIndex] = applySpectrogramScaling(
                        magnitude: magnitude,
                        phase: phase,
                        scaling: scaling
                    )
                }
            }
        }

        let sampRate: Double = samplingFrequencyInHz()

        return WaveformSpectrogram(
            timeFrames: timeFrames,
            frequencies: selectedFrequencies,
            data: spectrogramData,
            windowSize: windowSize,
            hopSize: actualHopSize,
            windowType: windowType,
            scaling: scaling,
            samplingRate: sampRate,
            originalDuration: duration
        )
    }

    /// Compute a mel-scale spectrogram for audio analysis
    /// - Parameters:
    ///   - windowSize: Size of the analysis window
    ///   - hopSize: Hop size between windows
    ///   - numMelBins: Number of mel-frequency bins (default: 80)
    ///   - minFreq: Minimum frequency for mel scale (default: 0 Hz)
    ///   - maxFreq: Maximum frequency for mel scale (default: Nyquist frequency)
    /// - Returns: Mel-scale spectrogram
    public func melSpectrogram(
        windowSize: Int = 512,
        hopSize: Int? = nil,
        numMelBins: Int = 80,
        minFreq: Double = 0.0,
        maxFreq: Double? = nil
    ) -> WaveformMelSpectrogram? {

        let nyquist: Double = nyquistFrequencyInHz()
        let actualMaxFreq = maxFreq ?? nyquist

        // First compute regular spectrogram
        guard
            let spec = self.spectrogram(
                windowSize: windowSize,
                hopSize: hopSize,
                scaling: .powerSpectralDensity
            )
        else { return nil }

        // Create mel filter bank
        let melFilters = createMelFilterBank(
            numBins: numMelBins,
            frequencyBins: spec.frequencies,
            minFreq: minFreq,
            maxFreq: actualMaxFreq
        )

        // Apply mel filters to spectrogram
        var melData: [[Double]] = []

        for timeFrame in spec.data {
            var melFrame: [Double] = []

            for melFilter in melFilters {
                var melValue = 0.0
                for (freqIndex, filterValue) in melFilter.enumerated() {
                    if freqIndex < timeFrame.count {
                        melValue += Double(timeFrame[freqIndex]) * filterValue
                    }
                }
                melFrame.append(melValue)
            }

            melData.append(melFrame)
        }

        // Create mel frequency axis
        let melFrequencies = (0..<numMelBins).map { i -> Double in
            let melMin = hzToMel(minFreq)
            let melMax = hzToMel(actualMaxFreq)
            let mel = melMin + (melMax - melMin) * Double(i) / Double(numMelBins - 1)
            return melToHz(mel)
        }

        let sampRate: Double = samplingFrequencyInHz()

        return WaveformMelSpectrogram(
            timeFrames: spec.timeFrames,
            melFrequencies: melFrequencies,
            data: melData,
            windowSize: windowSize,
            hopSize: hopSize ?? (windowSize / 4),
            numMelBins: numMelBins,
            samplingRate: PrecisionTimeInterval(seconds: sampRate)
        )
    }

    /// Compute instantaneous frequency using spectrogram phase derivatives
    /// - Parameters:
    ///   - windowSize: Analysis window size
    ///   - hopSize: Hop size between windows
    /// - Returns: Instantaneous frequency matrix
    public func spectrogramInstantaneousFrequency(
        windowSize: Int = 256,
        hopSize: Int? = nil
    ) -> WaveformInstantaneousFrequency<T>? {

        // Compute spectrograms with phase information
        guard
            let spec = self.spectrogram(
                windowSize: windowSize,
                hopSize: hopSize,
                scaling: .complex
            )
        else { return nil }

        let actualHopSize = hopSize ?? (windowSize / 4)
        let numFrames = spec.timeFrames.count
        let numFreqs = spec.frequencies.count

        var instFreqData: [[T]] = Array(
            repeating: Array(repeating: T.zero, count: numFreqs),
            count: max(0, numFrames - 1)
        )

        // Calculate instantaneous frequency from phase derivatives
        let dtSeconds = T(dt.secondsAsDouble)
        for frameIndex in 1..<numFrames {
            for freqIndex in 0..<numFreqs {
                let currentPhase = getPhaseFromComplex(spec.data[frameIndex][freqIndex])
                let previousPhase = getPhaseFromComplex(spec.data[frameIndex - 1][freqIndex])

                // Unwrap phase difference
                var phaseDiff = currentPhase - previousPhase
                while phaseDiff > T.pi { phaseDiff -= T(2.0 * Double.pi) }
                while phaseDiff < -T.pi { phaseDiff += T(2.0 * Double.pi) }

                // Convert to instantaneous frequency
                let deltaT = T(actualHopSize) * dtSeconds
                let instFreq = T(spec.frequencies[freqIndex]) + phaseDiff / (T(2.0 * Double.pi) * deltaT)

                instFreqData[frameIndex - 1][freqIndex] = instFreq
            }
        }

        // Convert timeFrames from PrecisionTimestamp to T (seconds offset)
        let timeFrameValues: [T] = spec.timeFrames.dropLast().map { ts in
            T(ts.interval.secondsAsDouble)
        }

        return WaveformInstantaneousFrequency(
            timeFrames: timeFrameValues,
            frequencies: spec.frequencies.map { T($0) },
            instantaneousFrequencies: instFreqData
        )
    }

    // MARK: - Private Helper Methods

    private func performSTFT(_ windowedSegment: [T]) -> (magnitudes: [T], phases: [T])? {
        // Create complex array for FFT
        var complex = windowedSegment.map { Complex<T>(real: $0, imaginary: T.zero) }

        // Pad to next power of 2
        let fftSize = spectrogramNextPowerOfTwo(complex.count)
        while complex.count < fftSize {
            complex.append(Complex<T>(real: T.zero, imaginary: T.zero))
        }

        // Perform FFT
        fft_cooleyTukey(&complex)

        // Extract magnitudes and phases (only positive frequencies)
        let nyquistBin = fftSize / 2
        let magnitudes = (0..<nyquistBin).map { complexMagnitude(complex[$0]) }
        let phases = (0..<nyquistBin).map { complexPhase(complex[$0]) }

        return (magnitudes, phases)
    }

    private func calculateFrequencyBins(windowSize: Int) -> [T] {
        let fftSize = spectrogramNextPowerOfTwo(windowSize)
        let nyquistBin = fftSize / 2
        let sampRate = 1.0 / dt.secondsAsDouble

        return (0..<nyquistBin).map { i in
            T(Double(i) * sampRate / Double(fftSize))
        }
    }

    private func getFrequencyBinRange(frequencyBins: [T], range: WaveformFrequencyRange<T>?) -> (start: Int, end: Int) {
        guard let range = range else {
            return (0, frequencyBins.count)
        }

        let startBin = frequencyBins.firstIndex { $0 >= range.minFreq } ?? 0
        let endBin = frequencyBins.lastIndex { $0 <= range.maxFreq }?.advanced(by: 1) ?? frequencyBins.count

        return (startBin, min(endBin, frequencyBins.count))
    }

    private func applySpectrogramScaling(magnitude: T, phase: T, scaling: WaveformSpectrogramScaling) -> T {
        switch scaling {
        case .magnitude:
            return magnitude
        case .power:
            return magnitude * magnitude
        case .powerSpectralDensity:
            return magnitude * magnitude
        case .decibel:
            return magnitude > T.zero ? T(20.0 * Foundation.log10(Double(magnitude))) : T(-Double.infinity)
        case .complex:
            // For complex scaling, we'd typically store both magnitude and phase
            // This is a simplified representation
            return magnitude
        }
    }

    private func getPhaseFromComplex(_ value: T) -> T {
        // This is a simplified implementation
        // In practice, you'd extract phase from stored complex representation
        return T.zero
    }

    private func createMelFilterBank(
        numBins: Int,
        frequencyBins: [Double],
        minFreq: Double,
        maxFreq: Double
    ) -> [[Double]] {
        let melMin = hzToMel(minFreq)
        let melMax = hzToMel(maxFreq)

        // Create mel-spaced center frequencies
        let melCenters = (0..<(numBins + 2)).map { i in
            melMin + (melMax - melMin) * Double(i) / Double(numBins + 1)
        }

        let hzCenters = melCenters.map { melToHz($0) }

        var filters: [[Double]] = []

        for i in 1..<(numBins + 1) {
            var filter = Array(repeating: 0.0, count: frequencyBins.count)

            let leftHz = hzCenters[i - 1]
            let centerHz = hzCenters[i]
            let rightHz = hzCenters[i + 1]

            for (freqIndex, freq) in frequencyBins.enumerated() {
                if freq >= leftHz && freq <= centerHz {
                    let denom = centerHz - leftHz
                    filter[freqIndex] = denom > 0 ? (freq - leftHz) / denom : 0.0
                } else if freq > centerHz && freq <= rightHz {
                    let denom = rightHz - centerHz
                    filter[freqIndex] = denom > 0 ? (rightHz - freq) / denom : 0.0
                }
            }

            filters.append(filter)
        }

        return filters
    }

    private func hzToMel(_ hz: Double) -> Double {
        return 2595.0 * Foundation.log10(1.0 + hz / 700.0)
    }

    private func melToHz(_ mel: Double) -> Double {
        return 700.0 * (pow(10.0, mel / 2595.0) - 1.0)
    }

    private func spectrogramNextPowerOfTwo(_ n: Int) -> Int {
        guard n > 1 else { return 1 }
        return 1 << Int(ceil(log2(Double(n))))
    }
}

// MARK: - Spectrogram Analysis and Utilities
extension Waveform1D where T: BinaryFloatingPoint & SIMDScalar {

    /// Extract spectral features from spectrogram
    /// - Parameter spectrogram: Input spectrogram
    /// - Returns: Spectral features including centroid, rolloff, etc.
    // TODO: WaveformSpectralFeatures uses PrecisionTimestamp for centroids and PrecisionTimeInterval
    // for rolloffs, which is semantically incorrect for spectral values. This is a design issue
    // in the support type. For now we convert to match the type definition.
    public static func extractSpectralFeatures(from spectrogram: WaveformSpectrogram<T>) -> WaveformSpectralFeatures {
        var spectralCentroids: [PrecisionTimestamp] = []
        var spectralRolloffs: [PrecisionTimeInterval] = []
        var spectralFluxes: [Double] = []

        for (timeIndex, timeFrame) in spectrogram.data.enumerated() {
            // Spectral centroid
            let totalEnergy = timeFrame.reduce(T.zero, +)
            if totalEnergy > T.zero {
                let weightedSum = zip(timeFrame, spectrogram.frequencies).reduce(0.0) { $0 + Double($1.0) * $1.1 }
                let centroidHz = weightedSum / Double(totalEnergy)
                spectralCentroids.append(PrecisionTimestamp(interval: PrecisionTimeInterval(seconds: centroidHz)))
            } else {
                spectralCentroids.append(PrecisionTimestamp(interval: .zero))
            }

            // Spectral rolloff (95% energy point)
            let rolloffThreshold = totalEnergy * T(0.95)
            var cumulativeEnergy = T.zero
            var rolloffFreq = spectrogram.frequencies.last ?? 0.0

            for (freqIndex, energy) in timeFrame.enumerated() {
                cumulativeEnergy += energy
                if cumulativeEnergy >= rolloffThreshold {
                    rolloffFreq = spectrogram.frequencies[freqIndex]
                    break
                }
            }
            spectralRolloffs.append(PrecisionTimeInterval(seconds: rolloffFreq))

            // Spectral flux (change from previous frame)
            if timeIndex > 0 {
                let previousFrame = spectrogram.data[timeIndex - 1]
                let flux = zip(timeFrame, previousFrame).reduce(0.0) { sum, pair in
                    let diff = Double(pair.0) - Double(pair.1)
                    return sum + max(0.0, diff)
                }
                spectralFluxes.append(flux)
            } else {
                spectralFluxes.append(0.0)
            }
        }

        return WaveformSpectralFeatures(
            spectralCentroids: spectralCentroids,
            spectralRolloffs: spectralRolloffs,
            spectralFluxes: spectralFluxes,
            timeFrames: spectrogram.timeFrames
        )
    }
}
