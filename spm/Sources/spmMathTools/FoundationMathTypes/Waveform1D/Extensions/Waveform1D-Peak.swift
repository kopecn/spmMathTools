import Foundation

// MARK: - Peak Detection
extension Waveform1D where T: BinaryFloatingPoint & Comparable {

    /// Detect peaks in the waveform with threshold and prominence filtering
    /// - Parameters:
    ///   - threshold: Minimum value to consider as a peak (default: no threshold)
    ///   - prominence: Minimum prominence required (default: no prominence filtering)
    ///   - minDistance: Minimum distance between peaks in samples (default: 1)
    ///   - edgePeaks: Whether to consider edge values as potential peaks (default: false)
    /// - Returns: Array of detected peaks
    public func detectPeaks(
        threshold: T? = nil,
        prominence: T? = nil,
        minDistance: Int = 1,
        edgePeaks: Bool = false
    ) -> [WaveformPeak<T>] {
        guard values.count >= 3 else { return [] }

        var candidatePeaks: [WaveformPeak<T>] = []
        let startIndex = edgePeaks ? 0 : 1
        let endIndex = edgePeaks ? values.count : values.count - 1

        // Find local maxima
        for i in startIndex..<endIndex {
            let isLocalMaximum: Bool

            if i == 0 {
                // First element - only check right neighbor
                isLocalMaximum = edgePeaks && i + 1 < values.count && values[i] > values[i + 1]
            } else if i == values.count - 1 {
                // Last element - only check left neighbor
                isLocalMaximum = edgePeaks && values[i] > values[i - 1]
            } else {
                // Interior element - check both neighbors with plateau handling
                let leftCondition = values[i] < values[i - 1]   // < (strict, mirroring peak's >)
                let rightCondition = values[i] <= values[i + 1] // <= (inclusive, mirroring peak's >=)
                
                // Must be at least as low as both neighbors, and strictly lower than at least one
                isLocalMaximum = leftCondition && rightCondition && 
                               (values[i] < values[i - 1] || values[i] < values[i + 1])
            }

            if isLocalMaximum {
                // Apply threshold filter
                if let threshold = threshold, values[i] < threshold {
                    continue
                }

                let peak = WaveformPeak<T>(
                    index: i,
                    value: values[i],
                    time: t0?.addingTimeInterval(TimeInterval(i) * dt),
                    timeOffset: TimeInterval(i) * dt
                )
                candidatePeaks.append(peak)
            }
        }

        // Apply minimum distance filtering
        var filteredPeaks = applyMinimumDistance(candidatePeaks, minDistance: minDistance)

        // Apply prominence filtering
        if let prominence = prominence {
            filteredPeaks = applyProminenceFilter(filteredPeaks, minimumProminence: prominence)
        }

        return filteredPeaks
    }

    /// Detect valleys (negative peaks) in the waveform
    /// - Parameters:
    ///   - threshold: Maximum value to consider as a valley (default: no threshold)
    ///   - prominence: Minimum prominence required (default: no prominence filtering)
    ///   - minDistance: Minimum distance between valleys in samples (default: 1)
    ///   - edgeValleys: Whether to consider edge values as potential valleys (default: false)
    /// - Returns: Array of detected valleys
    public func detectValleys(
        threshold: T? = nil,
        prominence: T? = nil,
        minDistance: Int = 1,
        edgeValleys: Bool = false
    ) -> [WaveformPeak<T>] {
        guard values.count >= 3 else { return [] }

        var candidateValleys: [WaveformPeak<T>] = []
        let startIndex = edgeValleys ? 0 : 1
        let endIndex = edgeValleys ? values.count : values.count - 1

        // Find local minima
        for i in startIndex..<endIndex {
            let isLocalMinimum: Bool

            if i == 0 {
                // First element - only check right neighbor
                isLocalMinimum = edgeValleys && i + 1 < values.count && values[i] < values[i + 1]
            } else if i == values.count - 1 {
                // Last element - only check left neighbor
                isLocalMinimum = edgeValleys && values[i] < values[i - 1]
            } else {
                // Interior element - check both neighbors with plateau handling
                let leftCondition = values[i] <= values[i - 1]  // <= instead of <
                let rightCondition = values[i] < values[i + 1]  // < instead of <

                // Must be at least as low as both neighbors, and strictly lower than at least one
                isLocalMinimum =
                    leftCondition && rightCondition && (values[i] < values[i - 1] || values[i] < values[i + 1])
            }

            if isLocalMinimum {
                // Apply threshold filter
                if let threshold = threshold, values[i] > threshold {
                    continue
                }

                let valley = WaveformPeak<T>(
                    index: i,
                    value: values[i],
                    time: t0?.addingTimeInterval(TimeInterval(i) * dt),
                    timeOffset: TimeInterval(i) * dt
                )
                candidateValleys.append(valley)
            }
        }

        // Apply minimum distance filtering
        var filteredValleys = applyMinimumDistance(candidateValleys, minDistance: minDistance)

        // Apply prominence filtering (for valleys, prominence is depth below surrounding peaks)
        if let prominence = prominence {
            filteredValleys = applyValleyProminenceFilter(filteredValleys, minimumProminence: prominence)
        }

        return filteredValleys
    }

    /// Find the most prominent peaks in the waveform
    /// - Parameters:
    ///   - count: Maximum number of peaks to return
    ///   - minDistance: Minimum distance between peaks in samples
    /// - Returns: Array of peaks sorted by prominence (highest first)
    public func findMostProminentPeaks(count: Int, minDistance: Int = 1) -> [WaveformPeak<T>] {
        let allPeaks = detectPeaks(minDistance: minDistance)
        let peaksWithProminence = allPeaks.map { peak in
            let prominence = calculateProminence(at: peak.index)
            return WaveformPeakWithProminence(peak: peak, prominence: prominence)
        }

        let sortedPeaks = peaksWithProminence.sorted { $0.prominence > $1.prominence }
        return Array(sortedPeaks.prefix(count).map { $0.peak })
    }

    // MARK: - Private Helper Methods

    private func applyMinimumDistance(_ peaks: [WaveformPeak<T>], minDistance: Int) -> [WaveformPeak<T>] {
        guard minDistance > 1 else { return peaks }

        var filteredPeaks: [WaveformPeak<T>] = []

        for peak in peaks.sorted(by: { $0.value > $1.value }) {  // Sort by height, highest first
            let tooClose = filteredPeaks.contains { existingPeak in
                abs(existingPeak.index - peak.index) < minDistance
            }

            if !tooClose {
                filteredPeaks.append(peak)
            }
        }

        return filteredPeaks.sorted { $0.index < $1.index }  // Sort back by index
    }

    private func applyProminenceFilter(_ peaks: [WaveformPeak<T>], minimumProminence: T) -> [WaveformPeak<T>] {
        return peaks.filter { peak in
            let prominence = calculateProminence(at: peak.index)
            return prominence >= minimumProminence
        }
    }

    private func applyValleyProminenceFilter(_ valleys: [WaveformPeak<T>], minimumProminence: T) -> [WaveformPeak<T>] {
        return valleys.filter { valley in
            let prominence = calculateValleyProminence(at: valley.index)
            return prominence >= minimumProminence
        }
    }

    private func calculateProminence(at peakIndex: Int) -> T {
        guard peakIndex >= 0 && peakIndex < values.count else { return T.zero }

        let peakValue = values[peakIndex]

        // Find the lowest point between this peak and higher peaks on both sides
        let leftBase = findLeftBase(from: peakIndex, peakValue: peakValue)
        let rightBase = findRightBase(from: peakIndex, peakValue: peakValue)

        let baseLevel = max(leftBase, rightBase)
        return peakValue - baseLevel
    }

    private func calculateValleyProminence(at valleyIndex: Int) -> T {
        guard valleyIndex >= 0 && valleyIndex < values.count else { return T.zero }

        let valleyValue = values[valleyIndex]

        // Find the highest point between this valley and lower valleys on both sides
        let leftBase = findLeftValleyBase(from: valleyIndex, valleyValue: valleyValue)
        let rightBase = findRightValleyBase(from: valleyIndex, valleyValue: valleyValue)

        let baseLevel = min(leftBase, rightBase)
        return baseLevel - valleyValue
    }

    private func findLeftBase(from peakIndex: Int, peakValue: T) -> T {
        var minValue = peakValue

        for i in stride(from: peakIndex - 1, through: 0, by: -1) {
            if values[i] < minValue {
                minValue = values[i]
            }
            if values[i] > peakValue {
                break  // Found a higher peak, stop here
            }
        }

        return minValue
    }

    private func findRightBase(from peakIndex: Int, peakValue: T) -> T {
        var minValue = peakValue

        for i in (peakIndex + 1)..<values.count {
            if values[i] < minValue {
                minValue = values[i]
            }
            if values[i] > peakValue {
                break  // Found a higher peak, stop here
            }
        }

        return minValue
    }

    private func findLeftValleyBase(from valleyIndex: Int, valleyValue: T) -> T {
        var maxValue = valleyValue

        for i in stride(from: valleyIndex - 1, through: 0, by: -1) {
            if values[i] > maxValue {
                maxValue = values[i]
            }
            if values[i] < valleyValue {
                break  // Found a lower valley, stop here
            }
        }

        return maxValue
    }

    private func findRightValleyBase(from valleyIndex: Int, valleyValue: T) -> T {
        var maxValue = valleyValue

        for i in (valleyIndex + 1)..<values.count {
            if values[i] > maxValue {
                maxValue = values[i]
            }
            if values[i] < valleyValue {
                break  // Found a lower valley, stop here
            }
        }

        return maxValue
    }
}

// MARK: - Integer Support
extension Waveform1D where T: BinaryInteger & Comparable {

    /// Detect peaks in integer waveforms
    public func detectPeaks(
        threshold: T? = nil,
        minDistance: Int = 1,
        edgePeaks: Bool = false
    ) -> [WaveformPeak<T>] {
        guard values.count >= 3 else { return [] }

        var candidatePeaks: [WaveformPeak<T>] = []
        let startIndex = edgePeaks ? 0 : 1
        let endIndex = edgePeaks ? values.count : values.count - 1

        // Find local maxima
        for i in startIndex..<endIndex {
            let isLocalMaximum: Bool

            if i == 0 {
                isLocalMaximum = i + 1 < values.count && values[i] > values[i + 1]
            } else if i == values.count - 1 {
                isLocalMaximum = values[i] > values[i - 1]
            } else {
                isLocalMaximum = values[i] > values[i - 1] && values[i] > values[i + 1]
            }

            if isLocalMaximum {
                if let threshold = threshold, values[i] < threshold {
                    continue
                }

                let peak = WaveformPeak(
                    index: i,
                    value: values[i],
                    time: t0?.addingTimeInterval(TimeInterval(i) * dt),
                    timeOffset: TimeInterval(i) * dt
                )
                candidatePeaks.append(peak)
            }
        }

        // Apply minimum distance filtering
        return applyMinimumDistanceInteger(candidatePeaks, minDistance: minDistance)
    }

    private func applyMinimumDistanceInteger(_ peaks: [WaveformPeak<T>], minDistance: Int) -> [WaveformPeak<T>] {
        guard minDistance > 1 else { return peaks }

        var filteredPeaks: [WaveformPeak<T>] = []

        for peak in peaks.sorted(by: { $0.value > $1.value }) {
            let tooClose = filteredPeaks.contains { existingPeak in
                abs(existingPeak.index - peak.index) < minDistance
            }

            if !tooClose {
                filteredPeaks.append(peak)
            }
        }

        return filteredPeaks.sorted { $0.index < $1.index }
    }
}
