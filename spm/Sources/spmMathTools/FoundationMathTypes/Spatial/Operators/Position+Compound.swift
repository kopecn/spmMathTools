import Foundation
import FoundationTypes
import simd

// MARK: - Compound Operators for Double Position

extension Position where T == Double {

    /// Add another position to this position (in-place)
    @inlinable
    public static func += (lhs: inout Position<Double>, rhs: Position<Double>) {
        lhs.vector += rhs.vector
    }

    /// Add a scalar to all components (in-place)
    @inlinable
    public static func += (lhs: inout Position<Double>, rhs: Double) {
        lhs.vector += SIMD3<Double>(repeating: rhs)
    }

    /// Subtract another position from this position (in-place)
    @inlinable
    public static func -= (lhs: inout Position<Double>, rhs: Position<Double>) {
        lhs.vector -= rhs.vector
    }

    /// Subtract a scalar from all components (in-place)
    @inlinable
    public static func -= (lhs: inout Position<Double>, rhs: Double) {
        lhs.vector -= SIMD3<Double>(repeating: rhs)
    }

    /// Scale this position by a scalar (in-place)
    @inlinable
    public static func *= (lhs: inout Position<Double>, rhs: Double) {
        lhs.vector *= rhs
    }

    /// Divide this position by a scalar (in-place)
    @inlinable
    public static func /= (lhs: inout Position<Double>, rhs: Double) {
        lhs.vector /= rhs
    }
}

// MARK: - Compound Operators for Float Position

extension Position where T == Float {

    /// Add another position to this position (in-place)
    @inlinable
    public static func += (lhs: inout Position<Float>, rhs: Position<Float>) {
        lhs.vector += rhs.vector
    }

    /// Add a scalar to all components (in-place)
    @inlinable
    public static func += (lhs: inout Position<Float>, rhs: Float) {
        lhs.vector += SIMD3<Float>(repeating: rhs)
    }

    /// Subtract another position from this position (in-place)
    @inlinable
    public static func -= (lhs: inout Position<Float>, rhs: Position<Float>) {
        lhs.vector -= rhs.vector
    }

    /// Subtract a scalar from all components (in-place)
    @inlinable
    public static func -= (lhs: inout Position<Float>, rhs: Float) {
        lhs.vector -= SIMD3<Float>(repeating: rhs)
    }

    /// Scale this position by a scalar (in-place)
    @inlinable
    public static func *= (lhs: inout Position<Float>, rhs: Float) {
        lhs.vector *= rhs
    }

    /// Divide this position by a scalar (in-place)
    @inlinable
    public static func /= (lhs: inout Position<Float>, rhs: Float) {
        lhs.vector /= rhs
    }
}
