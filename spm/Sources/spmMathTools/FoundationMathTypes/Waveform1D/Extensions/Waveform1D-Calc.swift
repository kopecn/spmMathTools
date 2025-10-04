import Foundation
import FoundationTypes

// Note: kvSIMD can be imported when available for performance optimization

// MARK: - Calc Mathematical Operations
extension Waveform1D where T: BinaryFloatingPoint {

    /// Compute the numerical integral of the waveform using trapezoidal rule
    /// - Parameter initialValue: Starting value for integration (default: 0)
    /// - Returns: New waveform representing the integral
    public func integrate(initialValue: T = T.zero) -> Waveform1D<T> {
        guard !values.isEmpty else { return Waveform1D(values: [], dt: dt, t0: t0) }

        var integral: [T] = [initialValue]
        let dtValue = T(dt)

        for i in 1..<values.count {
            let trapezoidArea = (values[i - 1] + values[i]) * dtValue / T(2.0)
            integral.append(integral.last! + trapezoidArea)
        }

        return Waveform1D(values: integral, dt: dt, t0: t0)
    }

    /// Compute the numerical derivative of the waveform using central difference
    /// - Returns: New waveform representing the derivative
    public func derivative() -> Waveform1D<T> {
        guard values.count >= 2 else { return Waveform1D(values: [], dt: dt, t0: t0) }

        var derivative: [T] = []
        let dtValue = T(dt)

        // Forward difference for first point
        derivative.append((values[1] - values[0]) / dtValue)

        // Central difference for middle points
        for i in 1..<(values.count - 1) {
            derivative.append((values[i + 1] - values[i - 1]) / (T(2.0) * dtValue))
        }

        // Backward difference for last point
        if values.count > 1 {
            derivative.append((values[values.count - 1] - values[values.count - 2]) / dtValue)
        }

        return Waveform1D(values: derivative, dt: dt, t0: t0)
    }
}
