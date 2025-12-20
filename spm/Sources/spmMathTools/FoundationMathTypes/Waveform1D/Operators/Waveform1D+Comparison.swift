import Foundation
import FoundationTypes

// MARK: - Comparison Methods for Double Waveforms

extension Waveform1D where T == Double {

    /// Element-wise equality
    public func elementsEqual(to other: Waveform1D<T,U>) -> [Bool] {
        let count = min(values.count, other.values.count)
        return (0..<count).map { values[$0] == other.values[$0] }
    }

    public func elementsEqual(to value: Double) -> [Bool] {
        return values.map { $0 == value }
    }

    /// Element-wise less than
    public func elementsLessThan(_ other: Waveform1D<T,U>) -> [Bool] {
        let count = min(values.count, other.values.count)
        return (0..<count).map { values[$0] < other.values[$0] }
    }

    public func elementsLessThan(_ value: Double) -> [Bool] {
        return values.map { $0 < value }
    }

    /// Element-wise greater than
    public func elementsGreaterThan(_ other: Waveform1D<T,U>) -> [Bool] {
        let count = min(values.count, other.values.count)
        return (0..<count).map { values[$0] > other.values[$0] }
    }

    public func elementsGreaterThan(_ value: Double) -> [Bool] {
        return values.map { $0 > value }
    }

    /// Element-wise approximate equality
    public func elementsApproximatelyEqual(to other: Waveform1D<T,U>, tolerance: Double = 1e-10) -> [Bool] {
        let count = min(values.count, other.values.count)
        return (0..<count).map { abs(values[$0] - other.values[$0]) <= tolerance }
    }

    public func elementsApproximatelyEqual(to value: Double, tolerance: Double = 1e-10) -> [Bool] {
        return values.map { abs($0 - value) <= tolerance }
    }
}

// MARK: - Comparison Methods for Float Waveforms

extension Waveform1D where T == Float {

    /// Element-wise equality
    public func elementsEqual(to other: Waveform1D<T,U>) -> [Bool] {
        let count = min(values.count, other.values.count)
        return (0..<count).map { values[$0] == other.values[$0] }
    }

    public func elementsEqual(to value: Float) -> [Bool] {
        return values.map { $0 == value }
    }

    /// Element-wise less than
    public func elementsLessThan(_ other: Waveform1D<T,U>) -> [Bool] {
        let count = min(values.count, other.values.count)
        return (0..<count).map { values[$0] < other.values[$0] }
    }

    public func elementsLessThan(_ value: Float) -> [Bool] {
        return values.map { $0 < value }
    }

    /// Element-wise greater than
    public func elementsGreaterThan(_ other: Waveform1D<T,U>) -> [Bool] {
        let count = min(values.count, other.values.count)
        return (0..<count).map { values[$0] > other.values[$0] }
    }

    public func elementsGreaterThan(_ value: Float) -> [Bool] {
        return values.map { $0 > value }
    }

    /// Element-wise approximate equality
    public func elementsApproximatelyEqual(to other: Waveform1D<T,U>, tolerance: Float = 1e-6) -> [Bool] {
        let count = min(values.count, other.values.count)
        return (0..<count).map { abs(values[$0] - other.values[$0]) <= tolerance }
    }

    public func elementsApproximatelyEqual(to value: Float, tolerance: Float = 1e-6) -> [Bool] {
        return values.map { abs($0 - value) <= tolerance }
    }
}

// MARK: - Comparison Methods for Int Waveforms

extension Waveform1D where T == Int, U == Double {

    /// Element-wise equality
    public func elementsEqual(to other: Waveform1D<T,U>) -> [Bool] {
        let count = min(values.count, other.values.count)
        return (0..<count).map { values[$0] == other.values[$0] }
    }

    public func elementsEqual(to value: Int) -> [Bool] {
        return values.map { $0 == value }
    }

    /// Element-wise less than
    public func elementsLessThan(_ other: Waveform1D<T,U>) -> [Bool] {
        let count = min(values.count, other.values.count)
        return (0..<count).map { values[$0] < other.values[$0] }
    }

    public func elementsLessThan(_ value: Int) -> [Bool] {
        return values.map { $0 < value }
    }

    /// Element-wise greater than
    public func elementsGreaterThan(_ other: Waveform1D<T,U>) -> [Bool] {
        let count = min(values.count, other.values.count)
        return (0..<count).map { values[$0] > other.values[$0] }
    }

    public func elementsGreaterThan(_ value: Int) -> [Bool] {
        return values.map { $0 > value }
    }
}
