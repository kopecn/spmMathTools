import Foundation
import FoundationTypes

// MARK: - Spectrogram and Time-Frequency Analysis
extension Waveform1D where T: BinaryFloatingPoint {

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
        frequencyRange: WaveformFrequencyRange? = nil
    ) -> WaveformSpectrogram<T>? {

        guard !values.isEmpty && windowSize > 0 && windowSize <= values.count else { return nil }

        let actualHopSize = hopSize ?? (windowSize / 4)
        guard actualHopSize > 0 else { return nil }

        // Calculate number of time frames
        let numFrames = max(1, (values.count - windowSize) / actualHopSize + 1)

        // Generate window coefficients
        let window = generateWindow(type: windowType, length: windowSize)

        // Calculate frequency bins
        let frequencyBins = calculateFrequencyBins(windowSize: windowSize)
        let (startBin, endBin) = getFrequencyBinRange(frequencyBins: frequencyBins, range: frequencyRange)
        let selectedFrequencies = Array(frequencyBins[startBin..<endBin])

        // Initialize spectrogram matrix
        var spectrogramData: [[T]] = Array(
            repeating: Array(repeating: T.zero, count: selectedFrequencies.count),
            count: numFrames
        )

        // Calculate time frames
        var timeFrames: [TimeInterval] = []

        // Process each time frame
        for frameIndex in 0..<numFrames {
            let startSample = frameIndex * actualHopSize
            let endSample = min(startSample + windowSize, values.count)

            // Calculate time for this frame (center of window)
            let frameTime = TimeInterval(startSample + windowSize / 2) * dt
            timeFrames.append(frameTime)

            // Extract and pad segment if necessary
            var segment = Array(values[startSample..<endSample])
            if segment.count < windowSize {
                let paddingNeeded = windowSize - segment.count
                segment.append(contentsOf: Array(repeating: T.zero, count: paddingNeeded))
            }

            // Apply window function
            let windowedSegment = zip(segment, window).map { T($0.0) * T($0.1) }

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

        return WaveformSpectrogram(
            timeFrames: timeFrames,
            frequencies: selectedFrequencies,
            data: spectrogramData,
            windowSize: windowSize,
            hopSize: actualHopSize,
            windowType: windowType,
            scaling: scaling,
            samplingRate: samplingFrequency,
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
    ) -> WaveformMelSpectrogram<T>? {

        let actualMaxFreq = maxFreq ?? nyquistFrequency

        // First compute regular spectrogram
        guard
            let spectrogram = self.spectrogram(
                windowSize: windowSize,
                hopSize: hopSize,
                scaling: .powerSpectralDensity
            )
        else { return nil }

        // Create mel filter bank
        let melFilters = createMelFilterBank(
            numBins: numMelBins,
            frequencyBins: spectrogram.frequencies.map { Double($0) },
            minFreq: minFreq,
            maxFreq: actualMaxFreq
        )

        // Apply mel filters to spectrogram
        var melData: [[T]] = []

        for timeFrame in spectrogram.data {
            var melFrame: [T] = []

            for melFilter in melFilters {
                var melValue = T.zero
                for (freqIndex, filterValue) in melFilter.enumerated() {
                    if freqIndex < timeFrame.count {
                        melValue += timeFrame[freqIndex] * T(filterValue)
                    }
                }
                melFrame.append(melValue)
            }

            melData.append(melFrame)
        }

        // Create mel frequency axis
        let melFrequencies = (0..<numMelBins).map { i in
            let melMin = hzToMel(minFreq)
            let melMax = hzToMel(actualMaxFreq)
            let mel = melMin + (melMax - melMin) * Double(i) / Double(numMelBins - 1)
            return T(melToHz(mel))
        }

        return WaveformMelSpectrogram(
            timeFrames: spectrogram.timeFrames,
            melFrequencies: melFrequencies,
            data: melData,
            windowSize: windowSize,
            hopSize: hopSize ?? (windowSize / 4),
            numMelBins: numMelBins,
            samplingRate: samplingFrequency
        )
    }

    /// Compute instantaneous frequency using spectrogram phase derivatives
    /// - Parameters:
    ///   - windowSize: Analysis window size
    ///   - hopSize: Hop size between windows
    /// - Returns: Instantaneous frequency matrix
    public func instantaneousFrequency(
        windowSize: Int = 256,
        hopSize: Int? = nil
    ) -> WaveformInstantaneousFrequency<T>? {

        // Compute spectrograms with phase information
        guard
            let spectrogram = self.spectrogram(
                windowSize: windowSize,
                hopSize: hopSize,
                scaling: .complex
            )
        else { return nil }

        let actualHopSize = hopSize ?? (windowSize / 4)
        let numFrames = spectrogram.timeFrames.count
        let numFreqs = spectrogram.frequencies.count

        var instFreqData: [[T]] = Array(
            repeating: Array(repeating: T.zero, count: numFreqs),
            count: max(0, numFrames - 1)
        )

        // Calculate instantaneous frequency from phase derivatives
        for frameIndex in 1..<numFrames {
            for freqIndex in 0..<numFreqs {
                let currentPhase = getPhaseFromComplex(spectrogram.data[frameIndex][freqIndex])
                let previousPhase = getPhaseFromComplex(spectrogram.data[frameIndex - 1][freqIndex])

                // Unwrap phase difference
                var phaseDiff = currentPhase - previousPhase
                while phaseDiff > T.pi { phaseDiff -= T(2.0 * Double.pi) }
                while phaseDiff < -T.pi { phaseDiff += T(2.0 * Double.pi) }

                // Convert to instantaneous frequency
                let deltaT = T(actualHopSize) * T(dt)
                let instFreq = spectrogram.frequencies[freqIndex] + phaseDiff / (T(2.0 * Double.pi) * deltaT)

                instFreqData[frameIndex - 1][freqIndex] = instFreq
            }
        }

        return WaveformInstantaneousFrequency(
            timeFrames: Array(spectrogram.timeFrames.dropLast()),
            frequencies: spectrogram.frequencies,
            instantaneousFrequencies: instFreqData
        )
    }

    // MARK: - Private Helper Methods

    private func performSTFT(_ windowedSegment: [T]) -> (magnitudes: [T], phases: [T])? {
        // Convert to doubles for FFT processing
        let doubleSegment = windowedSegment.map { Double($0) }

        // Create complex array for FFT
        var complex = doubleSegment.map { Complex(real: $0, imaginary: 0.0) }

        // Pad to next power of 2
        let fftSize = nextPowerOfTwo(complex.count)
        while complex.count < fftSize {
            complex.append(Complex(real: 0.0, imaginary: 0.0))
        }

        // Perform FFT (reuse from FFT extension)
        fft_cooleyTukey(&complex)

        // Extract magnitudes and phases (only positive frequencies)
        let nyquistBin = fftSize / 2
        let magnitudes = (0..<nyquistBin).map { T(complex[$0].magnitude) }
        let phases = (0..<nyquistBin).map { T(complex[$0].phase) }

        return (magnitudes, phases)
    }

    private func calculateFrequencyBins(windowSize: Int) -> [T] {
        let fftSize = nextPowerOfTwo(windowSize)
        let nyquistBin = fftSize / 2

        return (0..<nyquistBin).map { i in
            T(Double(i) * samplingFrequency / Double(fftSize))
        }
    }

    private func getFrequencyBinRange(frequencyBins: [T], range: WaveformFrequencyRange?) -> (start: Int, end: Int) {
        guard let range = range else {
            return (0, frequencyBins.count)
        }

        let startBin = frequencyBins.firstIndex { $0 >= T(range.minFreq) } ?? 0
        let endBin = frequencyBins.lastIndex { $0 <= T(range.maxFreq) }?.advanced(by: 1) ?? frequencyBins.count

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
                    filter[freqIndex] = (freq - leftHz) / (centerHz - leftHz)
                } else if freq > centerHz && freq <= rightHz {
                    filter[freqIndex] = (rightHz - freq) / (rightHz - centerHz)
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

    // Reuse from FFT extension
    private func nextPowerOfTwo(_ n: Int) -> Int {
        guard n > 1 else { return 1 }
        return 1 << Int(ceil(log2(Double(n))))
    }

    // Reuse from existing windowing functionality
    private func generateWindow(type: WaveformWindowType, length: Int) -> [Double] {
        // This would call the existing generateWindow function from the windowing extension
        return Waveform1D<Double>.generateWindow(type: type, length: length)
    }
}

// MARK: - Spectrogram Analysis and Utilities
extension Waveform1D where T: BinaryFloatingPoint {

    /// Extract spectral features from spectrogram
    /// - Parameter spectrogram: Input spectrogram
    /// - Returns: Spectral features including centroid, rolloff, etc.
    public static func extractSpectralFeatures(from spectrogram: WaveformSpectrogram<T>) -> WaveformSpectralFeatures<T>
    {
        var spectralCentroids: [T] = []
        var spectralRolloffs: [T] = []
        var spectralFluxes: [T] = []

        for (timeIndex, timeFrame) in spectrogram.data.enumerated() {
            // Spectral centroid
            let totalEnergy = timeFrame.reduce(T.zero, +)
            if totalEnergy > T.zero {
                let weightedSum = zip(timeFrame, spectrogram.frequencies).reduce(T.zero) { $0 + $1.0 * $1.1 }
                spectralCentroids.append(weightedSum / totalEnergy)
            } else {
                spectralCentroids.append(T.zero)
            }

            // Spectral rolloff (95% energy point)
            let rolloffThreshold = totalEnergy * T(0.95)
            var cumulativeEnergy = T.zero
            var rolloffFreq = spectrogram.frequencies.last ?? T.zero

            for (freqIndex, energy) in timeFrame.enumerated() {
                cumulativeEnergy += energy
                if cumulativeEnergy >= rolloffThreshold {
                    rolloffFreq = spectrogram.frequencies[freqIndex]
                    break
                }
            }
            spectralRolloffs.append(rolloffFreq)

            // Spectral flux (change from previous frame)
            if timeIndex > 0 {
                let previousFrame = spectrogram.data[timeIndex - 1]
                let flux = zip(timeFrame, previousFrame).reduce(T.zero) { sum, pair in
                    let diff = pair.0 - pair.1
                    return sum + max(T.zero, diff)
                }
                spectralFluxes.append(flux)
            } else {
                spectralFluxes.append(T.zero)
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
