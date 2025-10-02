import Foundation

// MARK: - Phase Analysis and Unwrapping
extension Waveform1D where T: BinaryFloatingPoint {

    /// Compute the instantaneous phase using Hilbert transform
    /// - Parameters:
    ///   - unwrap: Whether to unwrap the phase (default: true)
    ///   - method: Phase calculation method
    /// - Returns: New waveform containing instantaneous phase values
    public func instantaneousPhase(
        unwrap: Bool = true,
        method: WaveformPhaseMethod = .hilbert
    ) -> Waveform1D<T> {

        let phaseValues: [T]

        switch method {
        case .hilbert:
            phaseValues = hilbertPhase()
        case .fft:
            phaseValues = fftPhase()
        case .derivative:
            phaseValues = derivativePhase()
        }

        let finalPhase = unwrap ? unwrapPhase(phaseValues) : phaseValues

        return Waveform1D(values: finalPhase, dt: dt, t0: t0)
    }

    /// Unwrap phase values to remove 2π discontinuities
    /// - Parameters:
    ///   - phaseData: Array of phase values to unwrap
    ///   - threshold: Threshold for detecting phase jumps (default: π)
    /// - Returns: Unwrapped phase values
    public func unwrapPhase(_ phaseData: [T], threshold: T = T.pi) -> [T] {
        guard phaseData.count > 1 else { return phaseData }

        var unwrapped = phaseData
        var correction = T.zero

        for i in 1..<unwrapped.count {
            let difference = unwrapped[i] - unwrapped[i - 1]

            // Detect jumps greater than threshold
            if difference > threshold {
                correction -= T(2.0 * Double.pi)
            } else if difference < -threshold {
                correction += T(2.0 * Double.pi)
            }

            unwrapped[i] += correction
        }

        return unwrapped
    }

    /// Compute the instantaneous frequency from phase
    /// - Parameters:
    ///   - unwrappedPhase: Optional pre-computed unwrapped phase
    ///   - method: Phase calculation method if phase not provided
    /// - Returns: New waveform containing instantaneous frequency in Hz
    public func instantaneousFrequency(
        unwrappedPhase: [T]? = nil,
        method: WaveformPhaseMethod = .hilbert
    ) -> Waveform1D<T> {

        let phase = unwrappedPhase ?? instantaneousPhase(unwrap: true, method: method).values
        guard phase.count > 1 else {
            return Waveform1D(values: [], dt: dt, t0: t0)
        }

        // Compute frequency as derivative of phase divided by 2π
        var frequency: [T] = []
        let twoPi = T(2.0 * Double.pi)
        let dtValue = T(dt)

        // Forward difference for first point
        frequency.append((phase[1] - phase[0]) / (twoPi * dtValue))

        // Central difference for middle points
        for i in 1..<(phase.count - 1) {
            frequency.append((phase[i + 1] - phase[i - 1]) / (T(2.0) * twoPi * dtValue))
        }

        // Backward difference for last point
        if phase.count > 1 {
            frequency.append((phase[phase.count - 1] - phase[phase.count - 2]) / (twoPi * dtValue))
        }

        return Waveform1D(values: frequency, dt: dt, t0: t0)
    }

    /// Compute phase difference between two signals
    /// - Parameters:
    ///   - other: Other waveform to compare with
    ///   - method: Phase calculation method
    ///   - unwrap: Whether to unwrap the phase difference
    /// - Returns: Phase difference waveform
    public func phaseDifference(
        with other: Waveform1D<T>,
        method: WaveformPhaseMethod = .hilbert,
        unwrap: Bool = true
    ) -> Waveform1D<T>? {

        guard values.count == other.values.count else { return nil }

        let phase1 = instantaneousPhase(unwrap: false, method: method).values
        let phase2 = other.instantaneousPhase(unwrap: false, method: method).values

        let phaseDiff = zip(phase1, phase2).map { $0 - $1 }
        let finalPhaseDiff = unwrap ? unwrapPhase(phaseDiff) : phaseDiff

        return Waveform1D(values: finalPhaseDiff, dt: dt, t0: t0)
    }

    /// Compute group delay from phase response
    /// - Parameters:
    ///   - frequencies: Frequency array from FFT
    ///   - phases: Phase array from FFT
    /// - Returns: Group delay values in samples
    public func groupDelay(frequencies: [T], phases: [T]) -> [T] {
        guard frequencies.count == phases.count && frequencies.count > 1 else {
            return []
        }

        let unwrappedPhases = unwrapPhase(phases)
        var groupDelay: [T] = []

        // Compute negative derivative of phase with respect to frequency
        for i in 1..<(unwrappedPhases.count - 1) {
            let dPhase = unwrappedPhases[i + 1] - unwrappedPhases[i - 1]
            let dFreq = frequencies[i + 1] - frequencies[i - 1]

            if dFreq > T.zero {
                groupDelay.append(-dPhase / (T(2.0 * Double.pi) * dFreq))
            } else {
                groupDelay.append(T.zero)
            }
        }

        // Handle endpoints
        if !groupDelay.isEmpty {
            groupDelay.insert(groupDelay.first!, at: 0)
            groupDelay.append(groupDelay.last!)
        }

        return groupDelay
    }

    /// Compute phase coherence between two signals
    /// - Parameters:
    ///   - other: Other waveform to compare with
    ///   - windowSize: Window size for coherence calculation
    ///   - overlap: Overlap between windows (0.0 to 1.0)
    /// - Returns: Phase coherence values
    public func phaseCoherence(
        with other: Waveform1D<T>,
        windowSize: Int? = nil,
        overlap: Double = 0.5
    ) -> Waveform1D<T>? {

        guard values.count == other.values.count && values.count > 1 else { return nil }

        let actualWindowSize = windowSize ?? min(256, values.count / 4)
        let hopSize = Int(Double(actualWindowSize) * (1.0 - max(0.0, min(0.99, overlap))))

        var coherenceValues: [T] = []
        var startIndex = 0

        while startIndex + actualWindowSize <= values.count {
            let segment1 = Array(values[startIndex..<(startIndex + actualWindowSize)])
            let segment2 = Array(other.values[startIndex..<(startIndex + actualWindowSize)])

            let coherence = calculatePhaseCoherence(segment1, segment2)
            coherenceValues.append(coherence)

            startIndex += hopSize
        }

        // Create time axis for coherence values
        let coherenceDt = dt * Double(hopSize)
        let coherenceT0 = t0?.addingTimeInterval(TimeInterval(actualWindowSize / 2) * dt)

        return Waveform1D(values: coherenceValues, dt: coherenceDt, t0: coherenceT0)
    }

    // MARK: - Private Implementation Methods

    private func hilbertPhase() -> [T] {
        guard values.count > 2 else {
            return values.map { _ in T.zero }
        }

        var phases: [T] = []

        for i in 0..<values.count {
            let real = values[i]
            let imaginary = quadratureComponent(at: i)
            let phase = T(atan2(Double(imaginary), Double(real)))
            phases.append(phase)
        }

        return phases
    }

    private func fftPhase() -> [T] {
        guard let fftResult = fft() else {
            return Array(repeating: T.zero, count: values.count)
        }

        // Use phases from FFT, but we need to interpolate back to original length
        let phases = fftResult.phases

        // Simple linear interpolation to match original signal length
        var interpolatedPhases: [T] = []
        let scale = Double(phases.count - 1) / Double(values.count - 1)

        for i in 0..<values.count {
            let index = Double(i) * scale
            let lowerIndex = Int(index)
            let upperIndex = min(lowerIndex + 1, phases.count - 1)
            let fraction = T(index - Double(lowerIndex))

            if lowerIndex == upperIndex {
                interpolatedPhases.append(phases[lowerIndex])
            } else {
                let interpolated = phases[lowerIndex] + fraction * (phases[upperIndex] - phases[lowerIndex])
                interpolatedPhases.append(interpolated)
            }
        }

        return interpolatedPhases
    }

    private func derivativePhase() -> [T] {
        // Phase from derivative approximation (for sinusoidal signals)
        let deriv = derivative()

        return zip(values, deriv.values).map { value, derivative in
            if abs(value) > T(1e-10) {
                return T(atan2(Double(derivative), Double(value)))
            } else {
                return T.zero
            }
        }
    }

    private func quadratureComponent(at index: Int) -> T {
        // Simple quadrature filter approximation for Hilbert transform
        let windowSize = min(5, values.count / 4)
        guard windowSize > 1 else { return T.zero }

        var sum = T.zero
        var weightSum = T.zero

        for i in max(0, index - windowSize)...min(values.count - 1, index + windowSize) {
            if i != index {
                let distance = abs(i - index)
                let weight = T(1.0 / Double(distance))
                let sign = T((i > index) ? 1.0 : -1.0)
                sum += values[i] * weight * sign
                weightSum += weight
            }
        }

        return weightSum > T.zero ? sum / weightSum : T.zero
    }

    private func calculatePhaseCoherence(_ x: [T], _ y: [T]) -> T {
        guard x.count == y.count && !x.isEmpty else { return T.zero }

        // Calculate instantaneous phases
        let waveform1 = Waveform1D(values: x, dt: dt)
        let waveform2 = Waveform1D(values: y, dt: dt)

        let phase1 = waveform1.hilbertPhase()
        let phase2 = waveform2.hilbertPhase()

        // Calculate phase difference
        let phaseDiff = zip(phase1, phase2).map { $0 - $1 }

        // Calculate coherence as absolute value of mean complex exponential
        var sumReal = T.zero
        var sumImag = T.zero

        for diff in phaseDiff {
            sumReal += T(cos(Double(diff)))
            sumImag += T(sin(Double(diff)))
        }

        let meanReal = sumReal / T(phaseDiff.count)
        let meanImag = sumImag / T(phaseDiff.count)

        return (meanReal * meanReal + meanImag * meanImag).squareRoot()
    }
}

// MARK: - Phase-based Signal Analysis
extension Waveform1D where T: BinaryFloatingPoint {

    /// Detect phase-locked events
    /// - Parameters:
    ///   - referencePhase: Reference phase value (in radians)
    ///   - tolerance: Phase tolerance for detection (in radians)
    ///   - minInterval: Minimum time between events
    /// - Returns: Array of phase-locked events
    public func detectPhaseLockEvents(
        referencePhase: T,
        tolerance: T = T.pi / T(6.0),  // 30 degrees
        minInterval: TimeInterval? = nil
    ) -> [WaveformTriggerEvent<T>] {

        let phase = instantaneousPhase(unwrap: true).values
        var events: [WaveformTriggerEvent<T>] = []
        var lastEventIndex: Int?

        for (index, phaseValue) in phase.enumerated() {
            // Normalize phase difference to [-π, π]
            let phaseDiff = normalizePhase(phaseValue - referencePhase)

            if abs(phaseDiff) <= tolerance {
                // Check minimum interval
                if let lastIndex = lastEventIndex,
                    let minInterval = minInterval
                {
                    let timeSinceLastEvent = TimeInterval(index - lastIndex) * dt
                    if timeSinceLastEvent < minInterval {
                        continue
                    }
                }

                let event = WaveformTriggerEvent(
                    index: index,
                    value: values[index],
                    time: t0?.addingTimeInterval(TimeInterval(index) * dt),
                    timeOffset: TimeInterval(index) * dt,
                    type: .phase(referencePhase, tolerance: tolerance)
                )

                events.append(event)
                lastEventIndex = index
            }
        }

        return events
    }

    /// Calculate phase synchronization index between signals
    /// - Parameter other: Other waveform to compare with
    /// - Returns: Phase synchronization index (0 to 1)
    public func phaseSynchronizationIndex(with other: Waveform1D<T>) -> T? {
        guard let phaseDiff = phaseDifference(with: other, method: .hilbert, unwrap: false) else {
            return nil
        }

        // Calculate PLV (Phase Locking Value)
        var sumReal = T.zero
        var sumImag = T.zero

        for phase in phaseDiff.values {
            sumReal += T(cos(Double(phase)))
            sumImag += T(sin(Double(phase)))
        }

        let meanReal = sumReal / T(phaseDiff.values.count)
        let meanImag = sumImag / T(phaseDiff.values.count)

        return (meanReal * meanReal + meanImag * meanImag).squareRoot()
    }

    private func normalizePhase(_ phase: T) -> T {
        var normalized = phase
        let twoPi = T(2.0 * Double.pi)

        while normalized > T.pi {
            normalized -= twoPi
        }
        while normalized < -T.pi {
            normalized += twoPi
        }

        return normalized
    }
}
