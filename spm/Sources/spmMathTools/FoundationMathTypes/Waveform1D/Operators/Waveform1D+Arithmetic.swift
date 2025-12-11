import Foundation
import FoundationTypes

// MARK: - Arithmetic Operators for Double Waveforms

extension Waveform1D where T == Double {

    // MARK: Addition

    public static func + (lhs: Waveform1D<Double>, rhs: Waveform1D<Double>) -> Waveform1D<Double> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] + rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func + (lhs: Waveform1D<Double>, rhs: Double) -> Waveform1D<Double> {
        let result = lhs.values.map { $0 + rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func + (lhs: Double, rhs: Waveform1D<Double>) -> Waveform1D<Double> {
        return rhs + lhs
    }

    // MARK: Subtraction

    public static func - (lhs: Waveform1D<Double>, rhs: Waveform1D<Double>) -> Waveform1D<Double> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] - rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func - (lhs: Waveform1D<Double>, rhs: Double) -> Waveform1D<Double> {
        let result = lhs.values.map { $0 - rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func - (lhs: Double, rhs: Waveform1D<Double>) -> Waveform1D<Double> {
        let result = rhs.values.map { lhs - $0 }
        return Waveform1D(values: result, dt: rhs.dt, t0: rhs.t0)
    }

    // MARK: Multiplication

    public static func * (lhs: Waveform1D<Double>, rhs: Waveform1D<Double>) -> Waveform1D<Double> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] * rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func * (lhs: Waveform1D<Double>, rhs: Double) -> Waveform1D<Double> {
        let result = lhs.values.map { $0 * rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func * (lhs: Double, rhs: Waveform1D<Double>) -> Waveform1D<Double> {
        return rhs * lhs
    }

    // MARK: Division

    public static func / (lhs: Waveform1D<Double>, rhs: Waveform1D<Double>) -> Waveform1D<Double> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] / rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func / (lhs: Waveform1D<Double>, rhs: Double) -> Waveform1D<Double> {
        let result = lhs.values.map { $0 / rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func / (lhs: Double, rhs: Waveform1D<Double>) -> Waveform1D<Double> {
        let result = rhs.values.map { lhs / $0 }
        return Waveform1D(values: result, dt: rhs.dt, t0: rhs.t0)
    }
}

// MARK: - Arithmetic Operators for Float Waveforms

extension Waveform1D where T == Float {

    // MARK: Addition

    public static func + (lhs: Waveform1D<Float>, rhs: Waveform1D<Float>) -> Waveform1D<Float> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] + rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func + (lhs: Waveform1D<Float>, rhs: Float) -> Waveform1D<Float> {
        let result = lhs.values.map { $0 + rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func + (lhs: Float, rhs: Waveform1D<Float>) -> Waveform1D<Float> {
        return rhs + lhs
    }

    // MARK: Subtraction

    public static func - (lhs: Waveform1D<Float>, rhs: Waveform1D<Float>) -> Waveform1D<Float> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] - rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func - (lhs: Waveform1D<Float>, rhs: Float) -> Waveform1D<Float> {
        let result = lhs.values.map { $0 - rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func - (lhs: Float, rhs: Waveform1D<Float>) -> Waveform1D<Float> {
        let result = rhs.values.map { lhs - $0 }
        return Waveform1D(values: result, dt: rhs.dt, t0: rhs.t0)
    }

    // MARK: Multiplication

    public static func * (lhs: Waveform1D<Float>, rhs: Waveform1D<Float>) -> Waveform1D<Float> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] * rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func * (lhs: Waveform1D<Float>, rhs: Float) -> Waveform1D<Float> {
        let result = lhs.values.map { $0 * rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func * (lhs: Float, rhs: Waveform1D<Float>) -> Waveform1D<Float> {
        return rhs * lhs
    }

    // MARK: Division

    public static func / (lhs: Waveform1D<Float>, rhs: Waveform1D<Float>) -> Waveform1D<Float> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] / rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func / (lhs: Waveform1D<Float>, rhs: Float) -> Waveform1D<Float> {
        let result = lhs.values.map { $0 / rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func / (lhs: Float, rhs: Waveform1D<Float>) -> Waveform1D<Float> {
        let result = rhs.values.map { lhs / $0 }
        return Waveform1D(values: result, dt: rhs.dt, t0: rhs.t0)
    }
}

// MARK: - Arithmetic Operators for Int Waveforms

extension Waveform1D where T == Int {

    // MARK: Addition

    public static func + (lhs: Waveform1D<Int>, rhs: Waveform1D<Int>) -> Waveform1D<Int> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] + rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func + (lhs: Waveform1D<Int>, rhs: Int) -> Waveform1D<Int> {
        let result = lhs.values.map { $0 + rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func + (lhs: Int, rhs: Waveform1D<Int>) -> Waveform1D<Int> {
        return rhs + lhs
    }

    // MARK: Subtraction

    public static func - (lhs: Waveform1D<Int>, rhs: Waveform1D<Int>) -> Waveform1D<Int> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] - rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func - (lhs: Waveform1D<Int>, rhs: Int) -> Waveform1D<Int> {
        let result = lhs.values.map { $0 - rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func - (lhs: Int, rhs: Waveform1D<Int>) -> Waveform1D<Int> {
        let result = rhs.values.map { lhs - $0 }
        return Waveform1D(values: result, dt: rhs.dt, t0: rhs.t0)
    }

    // MARK: Multiplication

    public static func * (lhs: Waveform1D<Int>, rhs: Waveform1D<Int>) -> Waveform1D<Int> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] * rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func * (lhs: Waveform1D<Int>, rhs: Int) -> Waveform1D<Int> {
        let result = lhs.values.map { $0 * rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func * (lhs: Int, rhs: Waveform1D<Int>) -> Waveform1D<Int> {
        return rhs * lhs
    }

    // MARK: Division

    public static func / (lhs: Waveform1D<Int>, rhs: Waveform1D<Int>) -> Waveform1D<Int> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] / rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func / (lhs: Waveform1D<Int>, rhs: Int) -> Waveform1D<Int> {
        let result = lhs.values.map { $0 / rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func / (lhs: Int, rhs: Waveform1D<Int>) -> Waveform1D<Int> {
        let result = rhs.values.map { lhs / $0 }
        return Waveform1D(values: result, dt: rhs.dt, t0: rhs.t0)
    }

    // MARK: Remainder

    public static func % (lhs: Waveform1D<Int>, rhs: Waveform1D<Int>) -> Waveform1D<Int> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] % rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func % (lhs: Waveform1D<Int>, rhs: Int) -> Waveform1D<Int> {
        let result = lhs.values.map { $0 % rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func % (lhs: Int, rhs: Waveform1D<Int>) -> Waveform1D<Int> {
        let result = rhs.values.map { lhs % $0 }
        return Waveform1D(values: result, dt: rhs.dt, t0: rhs.t0)
    }
}
