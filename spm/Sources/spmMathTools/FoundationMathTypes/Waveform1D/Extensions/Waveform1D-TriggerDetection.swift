import Foundation
import FoundationTypes

// MARK: - Trigger Detection and Event Marking
extension Waveform1D where T: BinaryFloatingPoint & Comparable {

    /// Detect trigger events in the waveform
    /// - Parameters:
    ///   - trigger: Trigger configuration
    ///   - hysteresis: Optional hysteresis to prevent false triggers
    /// - Returns: Array of detected trigger events
    public func detectTriggers(
        trigger: WaveformTrigger<T>,
        hysteresis: T? = nil
    ) -> [WaveformTriggerEvent<T>] {
        guard !values.isEmpty else { return [] }

        var events: [WaveformTriggerEvent<T>] = []
        var triggerState = false
        var lastTriggerIndex: Int?

        let hysteresisValue = hysteresis ?? T.zero

        for (index, value) in values.enumerated() {
            let triggerCondition = evaluateTrigger(value: value, trigger: trigger)

            if !triggerState && triggerCondition {
                // Trigger activated
                if let lastIndex = lastTriggerIndex,
                    let minInterval = trigger.minimumInterval
                {
                    let timeSinceLastTrigger = TimeInterval(index - lastIndex) * dt
                    if timeSinceLastTrigger < minInterval {
                        continue  // Skip due to minimum interval
                    }
                }

                let event = WaveformTriggerEvent(
                    index: index,
                    value: value,
                    time: t0?.addingTimeInterval(TimeInterval(index) * dt),
                    timeOffset: TimeInterval(index) * dt,
                    type: trigger.type
                )

                events.append(event)
                triggerState = true
                lastTriggerIndex = index

            } else if triggerState && hysteresisValue > T.zero {
                // Check for hysteresis reset
                let resetCondition = evaluateHysteresisReset(
                    value: value,
                    trigger: trigger,
                    hysteresis: hysteresisValue
                )
                if resetCondition {
                    triggerState = false
                }
            } else if triggerState && hysteresisValue == T.zero {
                // Reset immediately if no hysteresis
                if !triggerCondition {
                    triggerState = false
                }
            }
        }

        return events
    }

    /// Detect edge triggers (rising/falling edges)
    /// - Parameters:
    ///   - edgeType: Type of edge to detect
    ///   - threshold: Threshold value for edge detection
    ///   - minInterval: Minimum time between triggers
    /// - Returns: Array of edge trigger events
    public func detectEdgeTriggers(
        edgeType: WaveformEdgeType,
        threshold: T,
        minInterval: TimeInterval? = nil
    ) -> [WaveformTriggerEvent<T>] {

        let trigger = WaveformTrigger<T>(
            type: .edge(edgeType, threshold: threshold),
            minimumInterval: minInterval
        )

        return detectTriggers(trigger: trigger)
    }

    /// Detect level triggers (above/below threshold)
    /// - Parameters:
    ///   - levelType: Type of level trigger
    ///   - threshold: Threshold value
    ///   - minInterval: Minimum time between triggers
    ///   - hysteresis: Hysteresis value to prevent chattering
    /// - Returns: Array of level trigger events
    public func detectLevelTriggers(
        levelType: WaveformLevelType,
        threshold: T,
        minInterval: TimeInterval? = nil,
        hysteresis: T? = nil
    ) -> [WaveformTriggerEvent<T>] {

        let trigger = WaveformTrigger<T>(
            type: .level(levelType, threshold: threshold),
            minimumInterval: minInterval
        )

        return detectTriggers(trigger: trigger, hysteresis: hysteresis)
    }

    /// Detect window triggers (value enters/exits a range)
    /// - Parameters:
    ///   - windowType: Type of window trigger
    ///   - lowerBound: Lower threshold
    ///   - upperBound: Upper threshold
    ///   - minInterval: Minimum time between triggers
    /// - Returns: Array of window trigger events
    public func detectWindowTriggers(
        windowType: WaveformWindowTriggerType,
        lowerBound: T,
        upperBound: T,
        minInterval: TimeInterval? = nil
    ) -> [WaveformTriggerEvent<T>] {

        let trigger = WaveformTrigger<T>(
            type: .window(windowType, lower: lowerBound, upper: upperBound),
            minimumInterval: minInterval
        )

        return detectTriggers(trigger: trigger)
    }

    /// Detect pattern-based triggers using template matching
    /// - Parameters:
    ///   - pattern: Template pattern to match
    ///   - threshold: Correlation threshold for pattern matching
    ///   - minInterval: Minimum time between triggers
    /// - Returns: Array of pattern trigger events
    public func detectPatternTriggers(
        pattern: [T],
        threshold: T,
        minInterval: TimeInterval? = nil
    ) -> [WaveformTriggerEvent<T>] {

        guard pattern.count > 0 && pattern.count <= values.count else { return [] }

        var events: [WaveformTriggerEvent<T>] = []
        var lastTriggerIndex: Int?

        let halfPatternLength = pattern.count / 2

        for i in 0...(values.count - pattern.count) {
            let segment = Array(values[i..<(i + pattern.count)])
            let correlation = calculateNormalizedCorrelation(segment, pattern)

            if correlation >= threshold {
                // Check minimum interval
                if let lastIndex = lastTriggerIndex,
                    let minInterval = minInterval
                {
                    let timeSinceLastTrigger = TimeInterval(i - lastIndex) * dt
                    if timeSinceLastTrigger < minInterval {
                        continue
                    }
                }

                let triggerIndex = i + halfPatternLength
                let event = WaveformTriggerEvent(
                    index: triggerIndex,
                    value: values[triggerIndex],
                    time: t0?.addingTimeInterval(TimeInterval(triggerIndex) * dt),
                    timeOffset: TimeInterval(triggerIndex) * dt,
                    type: .pattern(pattern, threshold: threshold)
                )

                events.append(event)
                lastTriggerIndex = i
            }
        }

        return events
    }

    /// Mark events on the waveform with annotations
    /// - Parameters:
    ///   - events: Array of trigger events to mark
    ///   - annotation: Annotation type for the events
    /// - Returns: Waveform with event markers
    public func withEventMarkers(
        events: [WaveformTriggerEvent<T>],
        annotation: WaveformEventAnnotation = .marker
    ) -> WaveformWithEvents<T> {

        let eventMarkers = events.map { event in
            WaveformEventMarker(
                index: event.index,
                time: event.time,
                timeOffset: event.timeOffset,
                annotation: annotation,
                metadata: ["triggerType": String(describing: event.type)]
            )
        }

        return WaveformWithEvents<T>(
            waveform: self,
            events: eventMarkers
        )
    }

    // MARK: - Private Helper Methods

    private func evaluateTrigger(value: T, trigger: WaveformTrigger<T>) -> Bool {
        switch trigger.type {
        case .edge(let edgeType, let threshold):
            return evaluateEdgeTrigger(value: value, edgeType: edgeType, threshold: threshold)

        case .level(let levelType, let threshold):
            return evaluateLevelTrigger(value: value, levelType: levelType, threshold: threshold)

        case .window(let windowType, let lower, let upper):
            return evaluateWindowTrigger(value: value, windowType: windowType, lower: lower, upper: upper)

        case .pattern(_, _):
            return false  // Pattern triggers are handled separately
        }
    }

    private func evaluateEdgeTrigger(value: T, edgeType: WaveformEdgeType, threshold: T) -> Bool {
        // Edge triggers require previous value for comparison
        // This is a simplified implementation - full implementation would track state
        switch edgeType {
        case .rising:
            return value > threshold
        case .falling:
            return value < threshold
        case .both:
            return true  // Simplified - would need state tracking
        }
    }

    private func evaluateLevelTrigger(value: T, levelType: WaveformLevelType, threshold: T) -> Bool {
        switch levelType {
        case .above:
            return value > threshold
        case .below:
            return value < threshold
        case .equal:
            return abs(value - threshold) < T(1e-10)  // Small epsilon for floating point comparison
        }
    }

    private func evaluateWindowTrigger(value: T, windowType: WaveformWindowTriggerType, lower: T, upper: T) -> Bool {
        let inWindow = value >= lower && value <= upper

        switch windowType {
        case .enter:
            return inWindow
        case .exit:
            return !inWindow
        }
    }

    private func evaluateHysteresisReset(value: T, trigger: WaveformTrigger<T>, hysteresis: T) -> Bool {
        switch trigger.type {
        case .level(.above, let threshold):
            return value < (threshold - hysteresis)
        case .level(.below, let threshold):
            return value > (threshold + hysteresis)
        default:
            return false
        }
    }

    private func calculateNormalizedCorrelation(_ x: [T], _ y: [T]) -> T {
        guard x.count == y.count && !x.isEmpty else { return T.zero }

        let xMean = x.reduce(T.zero, +) / T(x.count)
        let yMean = y.reduce(T.zero, +) / T(y.count)

        var numerator = T.zero
        var xSumSquares = T.zero
        var ySumSquares = T.zero

        for i in 0..<x.count {
            let xDev = x[i] - xMean
            let yDev = y[i] - yMean

            numerator += xDev * yDev
            xSumSquares += xDev * xDev
            ySumSquares += yDev * yDev
        }

        let denominator = (xSumSquares * ySumSquares).squareRoot()

        guard denominator > T.zero else { return T.zero }

        return numerator / denominator
    }
}

// MARK: - Integer Support
extension Waveform1D where T: BinaryInteger & Comparable {

    /// Detect simple threshold triggers for integer waveforms
    /// - Parameters:
    ///   - threshold: Threshold value
    ///   - direction: Trigger direction (above/below)
    ///   - minInterval: Minimum time between triggers
    /// - Returns: Array of trigger events
    public func detectThresholdTriggers(
        threshold: T,
        direction: WaveformLevelType,
        minInterval: TimeInterval? = nil
    ) -> [WaveformTriggerEvent<T>] {

        var events: [WaveformTriggerEvent<T>] = []
        var lastTriggerIndex: Int?

        for (index, value) in values.enumerated() {
            let triggered: Bool

            switch direction {
            case .above:
                triggered = value > threshold
            case .below:
                triggered = value < threshold
            case .equal:
                triggered = value == threshold
            }

            if triggered {
                // Check minimum interval
                if let lastIndex = lastTriggerIndex,
                    let minInterval = minInterval
                {
                    let timeSinceLastTrigger = TimeInterval(index - lastIndex) * dt
                    if timeSinceLastTrigger < minInterval {
                        continue
                    }
                }

                let event = WaveformTriggerEvent(
                    index: index,
                    value: value,
                    time: t0?.addingTimeInterval(TimeInterval(index) * dt),
                    timeOffset: TimeInterval(index) * dt,
                    type: .level(direction, threshold: threshold)
                )

                events.append(event)
                lastTriggerIndex = index
            }
        }

        return events
    }
}
