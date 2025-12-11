import Foundation
import FoundationTypes

// MARK: - Compound Assignment Operators for Double Waveforms

extension Waveform1D where T == Double {

    public static func += (lhs: inout Waveform1D<Double>, rhs: Waveform1D<Double>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] + rhs.values[$0] }
    }

    public static func += (lhs: inout Waveform1D<Double>, rhs: Double) {
        for i in lhs.values.indices {
            lhs.values[i] += rhs
        }
    }

    public static func -= (lhs: inout Waveform1D<Double>, rhs: Waveform1D<Double>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] - rhs.values[$0] }
    }

    public static func -= (lhs: inout Waveform1D<Double>, rhs: Double) {
        for i in lhs.values.indices {
            lhs.values[i] -= rhs
        }
    }

    public static func *= (lhs: inout Waveform1D<Double>, rhs: Waveform1D<Double>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] * rhs.values[$0] }
    }

    public static func *= (lhs: inout Waveform1D<Double>, rhs: Double) {
        for i in lhs.values.indices {
            lhs.values[i] *= rhs
        }
    }

    public static func /= (lhs: inout Waveform1D<Double>, rhs: Waveform1D<Double>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] / rhs.values[$0] }
    }

    public static func /= (lhs: inout Waveform1D<Double>, rhs: Double) {
        for i in lhs.values.indices {
            lhs.values[i] /= rhs
        }
    }
}

// MARK: - Compound Assignment Operators for Float Waveforms

extension Waveform1D where T == Float {

    public static func += (lhs: inout Waveform1D<Float>, rhs: Waveform1D<Float>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] + rhs.values[$0] }
    }

    public static func += (lhs: inout Waveform1D<Float>, rhs: Float) {
        for i in lhs.values.indices {
            lhs.values[i] += rhs
        }
    }

    public static func -= (lhs: inout Waveform1D<Float>, rhs: Waveform1D<Float>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] - rhs.values[$0] }
    }

    public static func -= (lhs: inout Waveform1D<Float>, rhs: Float) {
        for i in lhs.values.indices {
            lhs.values[i] -= rhs
        }
    }

    public static func *= (lhs: inout Waveform1D<Float>, rhs: Waveform1D<Float>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] * rhs.values[$0] }
    }

    public static func *= (lhs: inout Waveform1D<Float>, rhs: Float) {
        for i in lhs.values.indices {
            lhs.values[i] *= rhs
        }
    }

    public static func /= (lhs: inout Waveform1D<Float>, rhs: Waveform1D<Float>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] / rhs.values[$0] }
    }

    public static func /= (lhs: inout Waveform1D<Float>, rhs: Float) {
        for i in lhs.values.indices {
            lhs.values[i] /= rhs
        }
    }
}

// MARK: - Compound Assignment Operators for Int Waveforms

extension Waveform1D where T == Int {

    public static func += (lhs: inout Waveform1D<Int>, rhs: Waveform1D<Int>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] + rhs.values[$0] }
    }

    public static func += (lhs: inout Waveform1D<Int>, rhs: Int) {
        for i in lhs.values.indices {
            lhs.values[i] += rhs
        }
    }

    public static func -= (lhs: inout Waveform1D<Int>, rhs: Waveform1D<Int>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] - rhs.values[$0] }
    }

    public static func -= (lhs: inout Waveform1D<Int>, rhs: Int) {
        for i in lhs.values.indices {
            lhs.values[i] -= rhs
        }
    }

    public static func *= (lhs: inout Waveform1D<Int>, rhs: Waveform1D<Int>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] * rhs.values[$0] }
    }

    public static func *= (lhs: inout Waveform1D<Int>, rhs: Int) {
        for i in lhs.values.indices {
            lhs.values[i] *= rhs
        }
    }

    public static func /= (lhs: inout Waveform1D<Int>, rhs: Waveform1D<Int>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] / rhs.values[$0] }
    }

    public static func /= (lhs: inout Waveform1D<Int>, rhs: Int) {
        for i in lhs.values.indices {
            lhs.values[i] /= rhs
        }
    }

    public static func %= (lhs: inout Waveform1D<Int>, rhs: Waveform1D<Int>) {
        let count = min(lhs.values.count, rhs.values.count)
        lhs.values = (0..<count).map { lhs.values[$0] % rhs.values[$0] }
    }

    public static func %= (lhs: inout Waveform1D<Int>, rhs: Int) {
        for i in lhs.values.indices {
            lhs.values[i] %= rhs
        }
    }
}
