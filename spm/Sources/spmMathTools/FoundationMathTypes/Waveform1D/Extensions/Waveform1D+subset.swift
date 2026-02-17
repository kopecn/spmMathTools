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

        let startTime: Double = (startDate - t0).secondsAsDouble
        let endTime: Double = (endDate - t0).secondsAsDouble

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
    ///   - startTime: Start time for the subset (in seconds)
    ///   - endTime: End time for the subset (in seconds)
    ///   - paddingBefore: Number of samples to pad before start with first value
    ///   - paddingAfter: Number of samples to pad after end with last value
    ///   - retainT0: If true, keeps original t0. If false, sets t0 to nil
    ///   - WaveformTimeReference: Whether time is relative to waveform start or Unix epoch
    /// - Returns: New waveform subset or nil if time range is outside waveform bounds
    public func subset(
        fromTime startTime: Double,
        toTime endTime: Double,
        paddingBefore: Int = 0,
        paddingAfter: Int = 0,
        retainT0: Bool = true,
        WaveformTimeReference: WaveformTimeReference = .waveformStart
    ) -> Waveform1D<T>? {

        guard !values.isEmpty else { return nil }
        guard startTime <= endTime else { return nil }

        var adjustedStartTime: Double = startTime
        var adjustedEndTime: Double = endTime

        switch WaveformTimeReference {
        case .waveformStart:
            break
        case .epoch:
            guard let t0 = self.t0 else { return nil }
            let t0Interval: Double = t0.interval.secondsAsDouble
            adjustedStartTime = startTime - t0Interval
            adjustedEndTime = endTime - t0Interval
        }

        let dtSeconds = dt.secondsAsDouble
        let waveformDuration = Double(values.count - 1) * dtSeconds

        if adjustedEndTime < 0 || adjustedStartTime > waveformDuration {
            return nil
        }

        let startIndex = Int(round(adjustedStartTime / dtSeconds))
        let endIndex = Int(round(adjustedEndTime / dtSeconds))

        let clampedStartIndex = max(0, min(startIndex, values.count - 1))
        let clampedEndIndex = max(0, min(endIndex, values.count - 1))

        guard clampedStartIndex <= clampedEndIndex else { return nil }

        var subsetValues = Array(values[clampedStartIndex...clampedEndIndex])

        if paddingBefore > 0, let firstValue = subsetValues.first {
            let paddingValues = Array(repeating: firstValue, count: paddingBefore)
            subsetValues = paddingValues + subsetValues
        }

        if paddingAfter > 0, let lastValue = subsetValues.last {
            let paddingValues = Array(repeating: lastValue, count: paddingAfter)
            subsetValues = subsetValues + paddingValues
        }

        let newT0: PrecisionTimestamp?
        if retainT0, let originalT0 = self.t0 {
            let timeOffset = Double(clampedStartIndex - paddingBefore) * dtSeconds
            newT0 = originalT0.addingTimeInterval(add: timeOffset)
        } else {
            newT0 = nil
        }

        return Waveform1D(values: subsetValues, dt: dt, t0: newT0)
    }
}
