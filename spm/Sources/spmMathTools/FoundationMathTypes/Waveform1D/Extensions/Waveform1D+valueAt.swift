import Foundation
import FoundationTypes

extension Waveform1D where T: BinaryFloatingPoint {

    /// Access value at a given Date with linear interpolation
    /// - Parameters:
    ///   - date: The date to sample at
    ///   - clamp: If true, returns edge values when date is outside range. If false, returns nil.
    ///   - interpolationWindow: Maximum time difference (seconds) for interpolation. Beyond this, nearest sample is used.
    /// - Returns: Interpolated value or nil if date is invalid/outside range and clamp is false
    public func value(at date: Date, clamp: Bool = false, interpolationWindow: TimeInterval = 0.5) -> T? {
        guard let t0 = self.t0 else { return nil }

        let timeOffset = date.timeIntervalSince(t0)
        return value(
            atTime: timeOffset,
            clamp: clamp,
            interpolationWindow: interpolationWindow,
            WaveformTimeReference: .waveformStart
        )
    }

    /// Access value at a given time with linear interpolation
    /// - Parameters:
    ///   - time: Time in seconds
    ///   - clamp: If true, returns edge values when time is outside range. If false, returns nil.
    ///   - interpolationWindow: Maximum time difference (seconds) for interpolation. Beyond this, nearest sample is used.
    ///   - WaveformTimeReference: Whether time is relative to waveform start or Unix epoch
    /// - Returns: Interpolated value or nil if time is invalid/outside range and clamp is false
    public func value(
        atTime time: TimeInterval,
        clamp: Bool = false,
        interpolationWindow: TimeInterval = 0.5,
        WaveformTimeReference: WaveformTimeReference = .waveformStart
    ) -> T? {

        guard !values.isEmpty else { return nil }

        let adjustedTime: TimeInterval
        switch WaveformTimeReference {
        case .waveformStart:
            adjustedTime = time
        case .epoch:
            guard let t0 = self.t0 else { return nil }
            adjustedTime = time - t0.timeIntervalSince1970
        }

        let sampleIndex = adjustedTime / dt
        let waveformDuration = TimeInterval(values.count - 1) * dt

        // Check bounds
        if sampleIndex < 0 {
            return clamp ? values.first : nil
        }
        if adjustedTime > waveformDuration {
            return clamp ? values.last : nil
        }

        // Exact sample match
        let floorIndex = Int(sampleIndex)
        if floorIndex >= values.count - 1 {
            return values.last
        }

        let fractionalPart = sampleIndex - Double(floorIndex)

        // If very close to a sample point or outside interpolation window, return nearest
        if fractionalPart < 1e-10 || abs(fractionalPart * dt) > interpolationWindow {
            return fractionalPart < 0.5 ? values[floorIndex] : values[floorIndex + 1]
        }

        // Linear interpolation
        let v1 = values[floorIndex]
        let v2 = values[floorIndex + 1]
        let weight = T(fractionalPart)

        return v1 + weight * (v2 - v1)
    }
}

// Extension for integer types (no interpolation, nearest neighbor)
extension Waveform1D where T: BinaryInteger {

    /// Access value at a given Date (nearest neighbor for integer types)
    public func value(at date: Date, clamp: Bool = false) -> T? {
        guard let t0 = self.t0 else { return nil }

        let timeOffset = date.timeIntervalSince(t0)
        return value(atTime: timeOffset, clamp: clamp, WaveformTimeReference: .waveformStart)
    }

    /// Access value at a given time (nearest neighbor for integer types)
    public func value(
        atTime time: TimeInterval,
        clamp: Bool = false,
        WaveformTimeReference: WaveformTimeReference = .waveformStart
    ) -> T? {

        guard !values.isEmpty else { return nil }

        let adjustedTime: TimeInterval
        switch WaveformTimeReference {
        case .waveformStart:
            adjustedTime = time
        case .epoch:
            guard let t0 = self.t0 else { return nil }
            adjustedTime = time - t0.timeIntervalSince1970
        }

        let sampleIndex = adjustedTime / dt
        let waveformDuration = TimeInterval(values.count - 1) * dt

        // Check bounds
        if sampleIndex < 0 {
            return clamp ? values.first : nil
        }
        if adjustedTime > waveformDuration {
            return clamp ? values.last : nil
        }

        // Nearest neighbor
        let nearestIndex = Int(round(sampleIndex))
        let clampedIndex = min(max(nearestIndex, 0), values.count - 1)

        return values[clampedIndex]
    }
}
