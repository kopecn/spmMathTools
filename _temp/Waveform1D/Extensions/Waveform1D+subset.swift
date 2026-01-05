import Foundation
import FoundationTypes

extension Waveform1D {

    /// Extract a subset of the waveform using a date range
    /// - Parameters:
    ///   - startDate: Start date for the subset
    ///   - endDate: End date for the subset
    ///   - paddingBefore: Number of samples to pad before start with first value
    ///   - paddingAfter: Number of samples to pad after end with last value
    ///   - retainT0: If true, keeps original t0. If false, sets t0 to nil
    /// - Returns: New waveform subset or nil if date range is outside waveform bounds
    public func subset(
        from startDate: PrecisionTimestamp,
        to endDate: PrecisionTimestamp,
        paddingBefore: Int = 0,
        paddingAfter: Int = 0,
        retainT0: Bool = true
    ) -> Waveform1D<T>? {

        guard let t0 = self.t0 else { return nil }
        guard startDate <= endDate else { return nil }

        // Convert PrecisionTimeInterval to U (BinaryFloatingPoint)
        let startInterval = startDate - t0
        let startTime = {
            let seconds = U(startInterval.seconds) + U(startInterval.attoseconds) / U(PrecisionTimeInterval.attosecondsPerSecond)
            switch startInterval.sign {
            case .positive: return seconds
            case .negative: return -seconds
            case .zero: return 0
            }
        }()

        let endInterval = endDate - t0
        let endTime = {
            let seconds = U(endInterval.seconds) + U(endInterval.attoseconds) / U(PrecisionTimeInterval.attosecondsPerSecond)
            switch endInterval.sign {
            case .positive: return seconds
            case .negative: return -seconds
            case .zero: return 0
            }
        }()

        return subset(
            fromTime: startTime,
            toTime: endTime,
            paddingBefore: paddingBefore,
            paddingAfter: paddingAfter,
            retainT0: retainT0,
            WaveformTimeReference: .waveformStart
        )
    }

    /// Extract a subset of the waveform using a time range
    /// - Parameters:
    ///   - startTime: Start time for the subset
    ///   - endTime: End time for the subset
    ///   - paddingBefore: Number of samples to pad before start with first value
    ///   - paddingAfter: Number of samples to pad after end with last value
    ///   - retainT0: If true, keeps original t0. If false, sets t0 to nil
    ///   - WaveformTimeReference: Whether time is relative to waveform start or Unix epoch
    /// - Returns: New waveform subset or nil if time range is outside waveform bounds
    public func subset(
        fromTime startTime: T,
        toTime endTime: T,
        paddingBefore: Int = 0,
        paddingAfter: Int = 0,
        retainT0: Bool = true,
        WaveformTimeReference: WaveformTimeReference = .waveformStart
    ) -> Waveform1D<T>? {

        guard !values.isEmpty else { return nil }
        guard startTime <= endTime else { return nil }

        let adjustedStartTime: T
        let adjustedEndTime: T

        switch WaveformTimeReference {
        case .waveformStart:
            adjustedStartTime = startTime
            adjustedEndTime = endTime
        case .epoch:
            guard let t0 = self.t0 else { return nil }
            // Convert t0's interval from epoch to U
            let t0Seconds = T(t0.interval.seconds) + T(t0.interval.attoseconds) / T(PrecisionTimeInterval.attosecondsPerSecond)
            let t0Interval: T
            switch t0.interval.sign {
            case .positive: t0Interval = t0Seconds
            case .negative: t0Interval = -t0Seconds
            case .zero: t0Interval = 0
            }
            adjustedStartTime = startTime - t0Interval
            adjustedEndTime = endTime - t0Interval
        }

        let waveformDuration = U(values.count - 1) * dt

        // Check if range is completely outside waveform bounds
        if adjustedEndTime < 0 || adjustedStartTime > waveformDuration {
            return nil
        }

        // Calculate nearest sample indices
        let startIndex = Int(round(adjustedStartTime / dt))
        let endIndex = Int(round(adjustedEndTime / dt))

        // Clamp indices to valid range
        let clampedStartIndex = max(0, min(startIndex, values.count - 1))
        let clampedEndIndex = max(0, min(endIndex, values.count - 1))

        // Ensure we have a valid range
        guard clampedStartIndex <= clampedEndIndex else { return nil }

        // Extract subset values
        var subsetValues = Array(values[clampedStartIndex...clampedEndIndex])

        // Add padding if requested
        if paddingBefore > 0, let firstValue = subsetValues.first {
            let paddingValues = Array(repeating: firstValue, count: paddingBefore)
            subsetValues = paddingValues + subsetValues
        }

        if paddingAfter > 0, let lastValue = subsetValues.last {
            let paddingValues = Array(repeating: lastValue, count: paddingAfter)
            subsetValues = subsetValues + paddingValues
        }

        // Calculate new t0 if retaining
        let newT0: PrecisionTimestamp?
        if retainT0, let originalT0 = self.t0 {
            let timeOffset = T(clampedStartIndex - paddingBefore) * dt
            newT0 = originalT0 + timeOffset
        } else {
            newT0 = nil
        }

        return Waveform1D(values: subsetValues, dt: dt, t0: newT0)
    }
}
