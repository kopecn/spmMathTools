import Foundation
import FoundationTypes

// MARK: - Unary Operators for Double Waveforms

extension Waveform1D where T == Double {

    /// Unary plus (returns copy)
    public static prefix func + (waveform: Waveform1D<T>) -> Waveform1D<T> {
        return waveform
    }

    /// Negate all elements
    public static prefix func - (waveform: Waveform1D<T>) -> Waveform1D<T> {
        let negated = waveform.values.map { -$0 }
        return Waveform1D<T>(values: negated, dt: waveform.dt, t0: waveform.t0)
    }
}

// MARK: - Unary Operators for Float Waveforms

extension Waveform1D where T == Float {

    /// Unary plus (returns copy)
    public static prefix func + (waveform: Waveform1D<T>) -> Waveform1D<T> {
        return waveform
    }

    /// Negate all elements
    public static prefix func - (waveform: Waveform1D<T>) -> Waveform1D<T> {
        let negated = waveform.values.map { -$0 }
        return Waveform1D<T>(values: negated, dt: waveform.dt, t0: waveform.t0)
    }
}

// MARK: - Unary Operators for Int Waveforms

extension Waveform1D where T == Int {

    /// Unary plus (returns copy)
    public static prefix func + (waveform: Waveform1D<Int>) -> Waveform1D<Int> {
        return waveform
    }

    /// Negate all elements
    public static prefix func - (waveform: Waveform1D<Int>) -> Waveform1D<Int> {
        let negated = waveform.values.map { -$0 }
        return Waveform1D(values: negated, dt: waveform.dt, t0: waveform.t0)
    }

    /// Bitwise NOT
    public static prefix func ~ (waveform: Waveform1D<Int>) -> Waveform1D<Int> {
        let complemented = waveform.values.map { ~$0 }
        return Waveform1D(values: complemented, dt: waveform.dt, t0: waveform.t0)
    }
}
