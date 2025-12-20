import Foundation
import FoundationTypes

// MARK: - Unary Operators for Double Waveforms

extension Waveform1D where T == Double {

    /// Unary plus (returns copy)
    public static prefix func + (waveform: Waveform1D<T,U>) -> Waveform1D<T,U> {
        return waveform
    }

    /// Negate all elements
    public static prefix func - (waveform: Waveform1D<T,U>) -> Waveform1D<T,U> {
        let negated = waveform.values.map { -$0 }
        return Waveform1D<T,U>(values: negated, dt: waveform.dt, t0: waveform.t0)
    }
}

// MARK: - Unary Operators for Float Waveforms

extension Waveform1D where T == Float {

    /// Unary plus (returns copy)
    public static prefix func + (waveform: Waveform1D<T,U>) -> Waveform1D<T,U> {
        return waveform
    }

    /// Negate all elements
    public static prefix func - (waveform: Waveform1D<T,U>) -> Waveform1D<T,U> {
        let negated = waveform.values.map { -$0 }
        return Waveform1D<T,U>(values: negated, dt: waveform.dt, t0: waveform.t0)
    }
}

// MARK: - Unary Operators for Int Waveforms

extension Waveform1D where T == Int, U == Double {

    /// Unary plus (returns copy)
    public static prefix func + (waveform: Waveform1D<Int,U>) -> Waveform1D<Int,U> {
        return waveform
    }

    /// Negate all elements
    public static prefix func - (waveform: Waveform1D<Int,U>) -> Waveform1D<Int,U> {
        let negated = waveform.values.map { -$0 }
        return Waveform1D(values: negated, dt: waveform.dt, t0: waveform.t0)
    }

    /// Bitwise NOT
    public static prefix func ~ (waveform: Waveform1D<Int,U>) -> Waveform1D<Int,U> {
        let complemented = waveform.values.map { ~$0 }
        return Waveform1D(values: complemented, dt: waveform.dt, t0: waveform.t0)
    }
}

// MARK: - Unary Operators for Int Waveforms

extension Waveform1D where T == Int, U == Float {

    /// Unary plus (returns copy)
    public static prefix func + (waveform: Waveform1D<Int,U>) -> Waveform1D<Int,U> {
        return waveform
    }

    /// Negate all elements
    public static prefix func - (waveform: Waveform1D<Int,U>) -> Waveform1D<Int,U> {
        let negated = waveform.values.map { -$0 }
        return Waveform1D(values: negated, dt: waveform.dt, t0: waveform.t0)
    }

    /// Bitwise NOT
    public static prefix func ~ (waveform: Waveform1D<Int,U>) -> Waveform1D<Int,U> {
        let complemented = waveform.values.map { ~$0 }
        return Waveform1D(values: complemented, dt: waveform.dt, t0: waveform.t0)
    }
}