import Foundation

// MARK: - Zero-Crossing Detection
extension Waveform1D where T: BinaryFloatingPoint & Comparable {

    /// Detect zero crossings in the waveform
    /// - Parameters:
    ///   - threshold: Minimum threshold to consider as zero (default: 1e-10)
    ///   - direction: Type of crossings to detect (default: all)
    /// - Returns: Array of zero crossing information
    public func zeroCrossings(
        threshold: T = T(1e-10),
        direction: WaveformZeroCrossingDirection = .all
    ) -> [WaveformZeroCrossing<T>] {
        guard values.count >= 2 else { return [] }

        var crossings: [WaveformZeroCrossing<T>] = []

        for i in 0..<(values.count - 1) {
            let currentValue = values[i]
            let nextValue = values[i + 1]

            // Check if there's a sign change (accounting for threshold)
            let currentSign = signWithThreshold(currentValue, threshold: threshold)
            let nextSign = signWithThreshold(nextValue, threshold: threshold)

            // Skip if either value is within threshold (considered zero)
            guard currentSign != .zero && nextSign != .zero else { continue }

            // Check for crossing
            if currentSign != nextSign {
                let crossingType: WaveformZeroCrossingType

                switch (currentSign, nextSign) {
                case (.negative, .positive):
                    crossingType = .rising
                case (.positive, .negative):
                    crossingType = .falling
                default:
                    continue  // Should not happen given our checks above
                }

                // Check if this crossing type is requested
                guard direction.includes(crossingType) else { continue }

                // Calculate interpolated crossing point
                let interpolatedIndex = interpolateZeroCrossing(
                    value1: currentValue,
                    value2: nextValue,
                    index1: i,
                    threshold: threshold
                )

                let crossingTime = t0?.addingTimeInterval(interpolatedIndex * dt)

                let crossing = WaveformZeroCrossing(
                    sampleIndex: interpolatedIndex,
                    time: crossingTime,
                    timeOffset: interpolatedIndex * dt,
                    type: crossingType,
                    magnitude: abs(nextValue - currentValue)
                )

                crossings.append(crossing)
            }
        }

        return crossings
    }

    /// Count zero crossings in the waveform
    /// - Parameters:
    ///   - threshold: Minimum threshold to consider as zero
    ///   - direction: Type of crossings to count
    /// - Returns: Number of zero crossings
    public func zeroCrossingCount(
        threshold: T = T(1e-10),
        direction: WaveformZeroCrossingDirection = .all
    ) -> Int {
        return zeroCrossings(threshold: threshold, direction: direction).count
    }

    /// Calculate zero crossing rate (crossings per unit time)
    /// - Parameters:
    ///   - threshold: Minimum threshold to consider as zero
    ///   - direction: Type of crossings to count
    /// - Returns: Zero crossing rate in crossings per second
    public func zeroCrossingRate(
        threshold: T = T(1e-10),
        direction: WaveformZeroCrossingDirection = .all
    ) -> Double {
        let crossingCount = zeroCrossingCount(threshold: threshold, direction: direction)
        let totalDuration = duration

        guard totalDuration > 0 else { return 0.0 }
        return Double(crossingCount) / totalDuration
    }

    /// Get segments between zero crossings
    /// - Parameters:
    ///   - threshold: Minimum threshold to consider as zero
    ///   - includePartial: Include partial segments at start/end
    /// - Returns: Array of waveform segments between crossings
    public func segmentsBetweenZeroCrossings(
        threshold: T = T(1e-10),
        includePartial: Bool = false
    ) -> [Waveform1D<T>] {
        let crossings = zeroCrossings(threshold: threshold, direction: .all)
        guard !crossings.isEmpty else {
            return includePartial ? [self] : []
        }

        var segments: [Waveform1D<T>] = []
        var segmentIndices: [Int] = []

        // Add start index if including partial segments
        if includePartial {
            segmentIndices.append(0)
        }

        // Add crossing indices
        segmentIndices.append(contentsOf: crossings.map { Int(round($0.sampleIndex)) })

        // Add end index if including partial segments
        if includePartial {
            segmentIndices.append(values.count - 1)
        }

        // Create segments
        for i in 0..<(segmentIndices.count - 1) {
            let startIdx = segmentIndices[i]
            let endIdx = segmentIndices[i + 1]

            guard startIdx < endIdx && endIdx < values.count else { continue }

            let segmentValues = Array(values[startIdx...endIdx])
            let segmentT0 = t0?.addingTimeInterval(TimeInterval(startIdx) * dt)

            let segment = Waveform1D(values: segmentValues, dt: dt, t0: segmentT0)
            segments.append(segment)
        }

        return segments
    }

    // MARK: - Private Helper Methods

    private func signWithThreshold(_ value: T, threshold: T) -> NumericSign {
        if abs(value) <= threshold {
            return .zero
        } else if value > threshold {
            return .positive
        } else {
            return .negative
        }
    }

    private func interpolateZeroCrossing(value1: T, value2: T, index1: Int, threshold: T) -> Double {
        // Linear interpolation to find exact crossing point
        let denominator = value2 - value1

        guard abs(denominator) > threshold else {
            return Double(index1) + 0.5  // Midpoint if values are too close
        }

        let fraction = -Double(value1) / Double(denominator)
        return Double(index1) + max(0.0, min(1.0, fraction))
    }
}

// MARK: - Integer Support
extension Waveform1D where T: BinaryInteger & Comparable {

    /// Detect zero crossings for integer waveforms
    /// - Parameter direction: Type of crossings to detect
    /// - Returns: Array of zero crossing information
    public func zeroCrossings(direction: WaveformZeroCrossingDirection = .all) -> [WaveformZeroCrossing<T>] {
        guard values.count >= 2 else { return [] }

        var crossings: [WaveformZeroCrossing<T>] = []

        for i in 0..<(values.count - 1) {
            let currentValue = values[i]
            let nextValue = values[i + 1]

            // For integers, zero is exactly zero
            let currentIsPositive = currentValue > 0
            let currentIsNegative = currentValue < 0
            let nextIsPositive = nextValue > 0
            let nextIsNegative = nextValue < 0

            // Skip if either value is zero
            guard currentValue != 0 && nextValue != 0 else { continue }

            // Check for crossing
            var crossingType: WaveformZeroCrossingType?

            if currentIsNegative && nextIsPositive {
                crossingType = .rising
            } else if currentIsPositive && nextIsNegative {
                crossingType = .falling
            }

            guard let type = crossingType, direction.includes(type) else { continue }

            // For integers, crossing occurs at midpoint
            let interpolatedIndex = Double(i) + 0.5
            let crossingTime = t0?.addingTimeInterval(interpolatedIndex * dt)

            let crossing = WaveformZeroCrossing(
                sampleIndex: interpolatedIndex,
                time: crossingTime,
                timeOffset: interpolatedIndex * dt,
                type: type,
                magnitude: T(abs(Int(nextValue) - Int(currentValue)))
            )

            crossings.append(crossing)
        }

        return crossings
    }

    /// Count zero crossings for integer waveforms
    public func zeroCrossingCount(direction: WaveformZeroCrossingDirection = .all) -> Int {
        return zeroCrossings(direction: direction).count
    }

    /// Calculate zero crossing rate for integer waveforms
    public func zeroCrossingRate(direction: WaveformZeroCrossingDirection = .all) -> Double {
        let crossingCount = zeroCrossingCount(direction: direction)
        let totalDuration = duration

        guard totalDuration > 0 else { return 0.0 }
        return Double(crossingCount) / totalDuration
    }
}
