import Foundation
import FoundationTypes

// MARK: - Generic Numeric Extension
extension Waveform1D {

    /// Helper function to convert PrecisionTimeInterval to BinaryFloatingPoint type U
    private func intervalToU<U: BinaryFloatingPoint>(_ interval: PrecisionTimeInterval) -> U {
        let seconds = U(interval.seconds)
        let fractionalSeconds = U(interval.attoseconds) / U(PrecisionTimeInterval.attosecondsPerSecond)
        switch interval.sign {
        case .positive: return seconds + fractionalSeconds
        case .negative: return -(seconds + fractionalSeconds)
        case .zero: return 0
        }
    }

    /// Access value at a given Date
    /// - Parameters:
    ///   - date: The date to sample at
    ///   - clamp: If true, returns edge values when date is outside range. If false, returns nil.
    ///   - interpolationWindow: Maximum time difference (seconds) for interpolation. Beyond this, nearest sample is used. Only applies to BinaryFloatingPoint types.
    /// - Returns: Value at the date or nil if date is invalid/outside range and clamp is false
    public func value<U: BinaryFloatingPoint>(
        at date: PrecisionTimestamp,
        clamp: Bool = false,
        interpolationWindow: U = 0.5
    ) -> T? {
        guard let t0 = self.t0 else { return nil }

        // Convert PrecisionTimeInterval to U
        let interval = date - t0
        let timeOffset: U = intervalToU(interval)
        return value(
            atTime: timeOffset,
            clamp: clamp,
            interpolationWindow: interpolationWindow,
            WaveformTimeReference: .waveformStart
        )
    }

    /// Access value at a given time
    /// - Parameters:
    ///   - time: Time in seconds
    ///   - clamp: If true, returns edge values when time is outside range. If false, returns nil.
    ///   - interpolationWindow: Maximum time difference (seconds) for interpolation. Beyond this, nearest sample is used. Only applies to BinaryFloatingPoint types.
    ///   - WaveformTimeReference: Whether time is relative to waveform start or Unix epoch
    /// - Returns: Value at the time or nil if time is invalid/outside range and clamp is false
    public func value<U: BinaryFloatingPoint>(
        atTime time: U,
        clamp: Bool = false,
        interpolationWindow: U = 0.5,
        WaveformTimeReference: WaveformTimeReference = .waveformStart
    ) -> T? {

        guard !values.isEmpty else { return nil }

        let adjustedTime: U
        switch WaveformTimeReference {
        case .waveformStart:
            adjustedTime = time
        case .epoch:
            guard let t0 = self.t0 else { return nil }
            // Convert t0's interval from epoch to U
            let t0Interval: U = intervalToU(t0.interval)
            adjustedTime = time - t0Interval
        }

        let dtAsU: U = intervalToU(dt)
        let sampleIndex = adjustedTime / dtAsU
        let waveformDuration = U(values.count - 1) * dtAsU

        // Check bounds
        if sampleIndex < 0 {
            return clamp ? values.first : nil
        }
        if adjustedTime > waveformDuration {
            return clamp ? values.last : nil
        }

        // Get indices
        let floorIndex = Int(sampleIndex)
        if floorIndex >= values.count - 1 {
            return values.last
        }

        let fractionalPart = sampleIndex - U(floorIndex)

        // Determine if we should use interpolation or nearest neighbor
        // For BinaryFloatingPoint types, we can do interpolation
        // For other types (like integers), we use nearest neighbor
        if T.self is any BinaryFloatingPoint.Type {
            // If very close to a sample point or outside interpolation window, return nearest
            if fractionalPart < 1e-10 || abs(fractionalPart * dtAsU) > interpolationWindow {
                return fractionalPart < 0.5 ? values[floorIndex] : values[floorIndex + 1]
            }

            // Linear interpolation for BinaryFloatingPoint types
            let v1 = values[floorIndex]
            let v2 = values[floorIndex + 1]

            // Convert fractionalPart to T safely
            if let weight = convertToNumeric(fractionalPart, targetType: T.self) {
                return v1 + weight * (v2 - v1)
            } else {
                // Fallback to nearest neighbor if conversion fails
                return fractionalPart < 0.5 ? values[floorIndex] : values[floorIndex + 1]
            }
        } else {
            // Nearest neighbor for non-BinaryFloatingPoint types (e.g., integers)
            let nearestIndex = Int(round(Double(sampleIndex)))
            let clampedIndex = min(max(nearestIndex, 0), values.count - 1)
            return values[clampedIndex]
        }
    }

    /// Helper function to convert a BinaryFloatingPoint value to a Numeric type
    /// - Parameters:
    ///   - value: The floating point value to convert
    ///   - targetType: The target numeric type
    /// - Returns: Converted value or nil if conversion is not possible
    private func convertToNumeric<F: BinaryFloatingPoint, N: Numeric>(_ value: F, targetType: N.Type) -> N? {
        // Handle conversion based on target type
        switch N.self {
        case is Double.Type:
            return Double(value) as? N
        case is Float.Type:
            return Float(value) as? N
        case is CGFloat.Type:
            return CGFloat(value) as? N
        case is any BinaryInteger.Type:
            return N(exactly: Int(value))
        default:
            // For other numeric types, try integer conversion
            return N(exactly: Int(value))
        }
    }
}
