import Foundation
import FoundationTypes

// MARK: - Bitwise Operators for Int Waveforms

extension Waveform1D where T == Int {

    // MARK: Bitwise AND

    public static func & (lhs: Waveform1D<T>, rhs: Waveform1D<T>) -> Waveform1D<T> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] & rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func & (lhs: Waveform1D<T>, rhs: Int) -> Waveform1D<T> {
        let result = lhs.values.map { $0 & rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func & (lhs: Int, rhs: Waveform1D<T>) -> Waveform1D<T> {
        return rhs & lhs
    }

    public static func &= (lhs: inout Waveform1D<T>, rhs: Waveform1D<T>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] & rhs.values[$0] }
    }

    public static func &= (lhs: inout Waveform1D<T>, rhs: Int) {
        for i in lhs.values.indices {
            lhs.values[i] &= rhs
        }
    }

    // MARK: Bitwise OR

    public static func | (lhs: Waveform1D<T>, rhs: Waveform1D<T>) -> Waveform1D<T> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] | rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func | (lhs: Waveform1D<T>, rhs: Int) -> Waveform1D<T> {
        let result = lhs.values.map { $0 | rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func | (lhs: Int, rhs: Waveform1D<T>) -> Waveform1D<T> {
        return rhs | lhs
    }

    public static func |= (lhs: inout Waveform1D<T>, rhs: Waveform1D<T>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] | rhs.values[$0] }
    }

    public static func |= (lhs: inout Waveform1D<T>, rhs: Int) {
        for i in lhs.values.indices {
            lhs.values[i] |= rhs
        }
    }

    // MARK: Bitwise XOR

    public static func ^ (lhs: Waveform1D<T>, rhs: Waveform1D<T>) -> Waveform1D<T> {
        let count = min(lhs.values.count, rhs.values.count)
        let result = (0..<count).map { lhs.values[$0] ^ rhs.values[$0] }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func ^ (lhs: Waveform1D<T>, rhs: Int) -> Waveform1D<T> {
        let result = lhs.values.map { $0 ^ rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func ^ (lhs: Int, rhs: Waveform1D<T>) -> Waveform1D<T> {
        return rhs ^ lhs
    }

    public static func ^= (lhs: inout Waveform1D<T>, rhs: Waveform1D<T>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] ^ rhs.values[$0] }
    }

    public static func ^= (lhs: inout Waveform1D<T>, rhs: Int) {
        for i in lhs.values.indices {
            lhs.values[i] ^= rhs
        }
    }

    // MARK: Bitwise Shift

    public static func << (lhs: Waveform1D<T>, rhs: Int) -> Waveform1D<T> {
        let result = lhs.values.map { $0 << rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func <<= (lhs: inout Waveform1D<T>, rhs: Int) {
        for i in lhs.values.indices {
            lhs.values[i] <<= rhs
        }
    }

    public static func >> (lhs: Waveform1D<T>, rhs: Int) -> Waveform1D<T> {
        let result = lhs.values.map { $0 >> rhs }
        return Waveform1D(values: result, dt: lhs.dt, t0: lhs.t0)
    }

    public static func >>= (lhs: inout Waveform1D<T>, rhs: Int) {
        for i in lhs.values.indices {
            lhs.values[i] >>= rhs
        }
    }
}
