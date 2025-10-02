import Foundation

// MARK: - Time Alignment and Synchronization
extension Waveform1D where T: BinaryFloatingPoint {

    /// Align this waveform to match another waveform's time characteristics
    /// - Parameters:
    ///   - reference: Reference waveform to align to
    ///   - method: Alignment method to use
    ///   - maxLag: Maximum lag to search for optimal alignment (in samples)
    /// - Returns: Aligned waveform
    public func aligned(
        to reference: Waveform1D<T>,
        method: WaveformAlignmentMethod = .crossCorrelation,
        maxLag: Int? = nil
    ) -> Waveform1D<T> {

        switch method {
        case .crossCorrelation:
            return alignByCrossCorrelation(to: reference, maxLag: maxLag)
        case .timeStamp:
            return alignByTimeStamp(to: reference)
        case .manualOffset(let offset):
            return alignByManualOffset(offset)
        case .peakAlignment:
            return alignByPeaks(to: reference)
        case .energyAlignment:
            return alignByEnergy(to: reference, maxLag: maxLag)
        }
    }

    /// Synchronize multiple waveforms to a common time base
    /// - Parameters:
    ///   - waveforms: Array of waveforms to synchronize
    ///   - reference: Reference waveform (if nil, uses first waveform)
    ///   - method: Synchronization method
    ///   - resampleToCommonRate: Whether to resample all to common sampling rate
    /// - Returns: Array of synchronized waveforms
    public static func synchronize(
        _ waveforms: [Waveform1D<T>],
        to reference: Waveform1D<T>? = nil,
        method: WaveformAlignmentMethod = .crossCorrelation,
        resampleToCommonRate: Bool = true
    ) -> [Waveform1D<T>] {

        guard !waveforms.isEmpty else { return [] }

        let referenceWaveform = reference ?? waveforms[0]
        var synchronized: [Waveform1D<T>] = []

        for waveform in waveforms {
            var alignedWaveform = waveform

            // Resample to common sampling rate if requested
            if resampleToCommonRate && abs(waveform.samplingFrequency - referenceWaveform.samplingFrequency) > 1e-10 {
                alignedWaveform = waveform.resampled(to: referenceWaveform.samplingFrequency)
            }

            // Align to reference
            alignedWaveform = alignedWaveform.aligned(to: referenceWaveform, method: method)

            synchronized.append(alignedWaveform)
        }

        return synchronized
    }

    /// Calculate time lag between two waveforms using cross-correlation
    /// - Parameters:
    ///   - other: Other waveform to compare with
    ///   - maxLag: Maximum lag to search (in samples)
    /// - Returns: Time lag information
    public func timeLag(to other: Waveform1D<T>, maxLag: Int? = nil) -> WaveformTimeLag<T>? {
        let searchRange = maxLag.map { -$0..<$0 }

        guard let correlation = findMaxCorrelation(with: other, searchRange: searchRange) else {
            return nil
        }

        return WaveformTimeLag(
            lagSamples: correlation.lagSamples,
            lagTime: correlation.lagTime,
            correlation: correlation.correlation,
            confidence: calculateLagConfidence(correlation.correlation)
        )
    }

    /// Trim waveforms to their common time overlap
    /// - Parameter waveforms: Array of waveforms to trim
    /// - Returns: Array of trimmed waveforms with common time base
    public static func trimToCommonTimeBase(_ waveforms: [Waveform1D<T>]) -> [Waveform1D<T>] {
        guard waveforms.count > 1 else { return waveforms }

        // Find common time range
        var latestStart: Date?
        var earliestEnd: Date?

        for waveform in waveforms {
            guard let start = waveform.t0, let end = waveform.endTime else { continue }

            if latestStart == nil || start > latestStart! {
                latestStart = start
            }

            if earliestEnd == nil || end < earliestEnd! {
                earliestEnd = end
            }
        }

        guard let commonStart = latestStart, let commonEnd = earliestEnd,
            commonStart < commonEnd
        else {
            return waveforms
        }

        // Trim each waveform to common time range
        return waveforms.compactMap { waveform in
            waveform.subset(from: commonStart, to: commonEnd)
        }
    }

    // MARK: - Private Alignment Methods

    private func alignByCrossCorrelation(to reference: Waveform1D<T>, maxLag: Int?) -> Waveform1D<T> {
        guard let lag = timeLag(to: reference, maxLag: maxLag) else { return self }

        return alignByManualOffset(lag.lagTime)
    }

    private func alignByTimeStamp(to reference: Waveform1D<T>) -> Waveform1D<T> {
        guard let referenceT0 = reference.t0, let selfT0 = t0 else { return self }

        let timeOffset = selfT0.timeIntervalSince(referenceT0)
        return alignByManualOffset(-timeOffset)
    }

    private func alignByManualOffset(_ offsetSeconds: TimeInterval) -> Waveform1D<T> {
        let newT0 = t0?.addingTimeInterval(offsetSeconds)
        return Waveform1D(values: values, dt: dt, t0: newT0)
    }

    private func alignByPeaks(to reference: Waveform1D<T>) -> Waveform1D<T> {
        // Find the most prominent peak in each waveform
        let selfPeaks = detectPeaks(minDistance: max(1, values.count / 20))
        let refPeaks = reference.detectPeaks(minDistance: max(1, reference.values.count / 20))

        guard let selfMainPeak = selfPeaks.max(by: { $0.value < $1.value }),
            let refMainPeak = refPeaks.max(by: { $0.value < $1.value })
        else {
            return self
        }

        // Calculate time offset needed to align peaks
        let selfPeakTime = selfMainPeak.timeOffset
        let refPeakTime = refMainPeak.timeOffset
        let timeOffset = refPeakTime - selfPeakTime

        return alignByManualOffset(timeOffset)
    }

    private func alignByEnergy(to reference: Waveform1D<T>, maxLag: Int?) -> Waveform1D<T> {
        // Create energy envelopes
        let selfEnergy = amplitudeEnvelope(method: .rms(windowSize: max(10, values.count / 100)))
        let refEnergy = reference.amplitudeEnvelope(method: .rms(windowSize: max(10, reference.values.count / 100)))

        // Use cross-correlation on energy envelopes
        guard let lag = selfEnergy.timeLag(to: refEnergy, maxLag: maxLag) else { return self }

        return alignByManualOffset(lag.lagTime)
    }

    private func calculateLagConfidence(_ correlation: T) -> Double {
        // Simple confidence measure based on correlation strength
        let absCorr = abs(Double(correlation))
        return min(1.0, max(0.0, (absCorr - 0.1) / 0.9))
    }
}

// MARK: - Time-based Windowing and Segmentation
extension Waveform1D where T: BinaryFloatingPoint {

    /// Create overlapping time windows from the waveform
    /// - Parameters:
    ///   - windowDuration: Duration of each window in seconds
    ///   - overlap: Overlap between windows as fraction (0.0 to 1.0)
    ///   - padding: Padding strategy for incomplete windows
    /// - Returns: Array of windowed waveforms
    public func timeWindows(
        duration windowDuration: TimeInterval,
        overlap: Double = 0.0,
        padding: WaveformPaddingStrategy = .none
    ) -> [Waveform1D<T>] {

        guard windowDuration > 0 && windowDuration <= duration else { return [] }

        let windowSamples = Int(windowDuration / dt)
        let stepSamples = Int(Double(windowSamples) * (1.0 - max(0.0, min(0.99, overlap))))

        var windows: [Waveform1D<T>] = []
        var startIndex = 0

        while startIndex < values.count {
            let endIndex = min(startIndex + windowSamples, values.count)
            var windowValues = Array(values[startIndex..<endIndex])

            // Handle incomplete windows
            if windowValues.count < windowSamples {
                switch padding {
                case .none:
                    break  // Keep incomplete window
                case .zeros:
                    let paddingCount = windowSamples - windowValues.count
                    windowValues.append(contentsOf: Array(repeating: T.zero, count: paddingCount))
                case .lastValue:
                    if let lastValue = windowValues.last {
                        let paddingCount = windowSamples - windowValues.count
                        windowValues.append(contentsOf: Array(repeating: lastValue, count: paddingCount))
                    }
                case .mirror:
                    while windowValues.count < windowSamples && windowValues.count > 1 {
                        let mirrorIndex = windowValues.count - 2
                        if mirrorIndex >= 0 {
                            windowValues.append(windowValues[mirrorIndex])
                        } else {
                            break
                        }
                    }
                }
            }

            // Create windowed waveform with adjusted t0
            let windowT0 = t0?.addingTimeInterval(TimeInterval(startIndex) * dt)
            let window = Waveform1D(values: windowValues, dt: dt, t0: windowT0)
            windows.append(window)

            startIndex += stepSamples

            // Break if we've covered the signal and don't want partial windows
            if padding == .none && startIndex >= values.count {
                break
            }
        }

        return windows
    }

    /// Create non-overlapping time segments of equal duration
    /// - Parameters:
    ///   - segmentDuration: Duration of each segment in seconds
    ///   - discardIncomplete: If true, discards the last segment if incomplete
    /// - Returns: Array of time-segmented waveforms
    public func timeSegments(
        duration segmentDuration: TimeInterval,
        discardIncomplete: Bool = false
    ) -> [Waveform1D<T>] {

        let padding: WaveformPaddingStrategy = discardIncomplete ? .none : .lastValue
        return timeWindows(duration: segmentDuration, overlap: 0.0, padding: padding)
    }
}

// MARK: - Integer Support
extension Waveform1D where T: BinaryInteger {

    /// Simple time alignment for integer waveforms using timestamps
    /// - Parameter reference: Reference waveform to align to
    /// - Returns: Aligned waveform
    public func aligned(to reference: Waveform1D<T>) -> Waveform1D<T> {
        guard let referenceT0 = reference.t0, let selfT0 = t0 else { return self }

        let timeOffset = selfT0.timeIntervalSince(referenceT0)
        let newT0 = t0?.addingTimeInterval(-timeOffset)

        return Waveform1D(values: values, dt: dt, t0: newT0)
    }

    /// Create time windows for integer waveforms
    /// - Parameters:
    ///   - windowDuration: Duration of each window in seconds
    ///   - overlap: Overlap between windows as fraction
    /// - Returns: Array of windowed waveforms
    public func timeWindows(
        duration windowDuration: TimeInterval,
        overlap: Double = 0.0
    ) -> [Waveform1D<T>] {

        guard windowDuration > 0 && windowDuration <= duration else { return [] }

        let windowSamples = Int(windowDuration / dt)
        let stepSamples = Int(Double(windowSamples) * (1.0 - max(0.0, min(0.99, overlap))))

        var windows: [Waveform1D<T>] = []
        var startIndex = 0

        while startIndex + windowSamples <= values.count {
            let windowValues = Array(values[startIndex..<(startIndex + windowSamples)])
            let windowT0 = t0?.addingTimeInterval(TimeInterval(startIndex) * dt)
            let window = Waveform1D(values: windowValues, dt: dt, t0: windowT0)
            windows.append(window)

            startIndex += stepSamples
        }

        return windows
    }
}
