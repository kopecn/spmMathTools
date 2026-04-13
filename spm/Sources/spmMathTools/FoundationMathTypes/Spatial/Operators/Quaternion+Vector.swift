import Foundation
import FoundationTypes
import simd

// MARK: - Quaternion Operators for Double Quaternion

extension Quaternion where T == Double {

    /// Compute the conjugate of this quaternion
    @inlinable
    public var conjugate: Quaternion<Double> {
        Quaternion(vector: vector * SIMD4<Double>(-1, -1, -1, 1))
    }

    /// Compute the inverse of this quaternion
    @inlinable
    public var inverse: Quaternion<Double> {
        let lenSq = simd_length_squared(vector)
        if lenSq > Double.ulpOfOne {
            return Quaternion(vector: vector * SIMD4<Double>(-1, -1, -1, 1) / lenSq)
        }
        return .identity
    }

    /// Compute dot product of two quaternions
    @inlinable
    public func dot(_ other: Quaternion<Double>) -> Double {
        simd_dot(vector, other.vector)
    }

    /// Compute dot product using operator
    @inlinable
    public static func • (lhs: Quaternion<Double>, rhs: Quaternion<Double>) -> Double {
        simd_dot(lhs.vector, rhs.vector)
    }

    /// Rotate a position by this quaternion
    @inlinable
    public func rotate(_ position: Position<Double>) -> Position<Double> {
        Position(vector: Self._qrot(vector, position.vector))
    }

    /// Rotate a position using operator
    @inlinable
    public static func * (quaternion: Quaternion<Double>, position: Position<Double>) -> Position<Double> {
        quaternion.rotate(position)
    }
}

// MARK: - Quaternion Operators for Float Quaternion

extension Quaternion where T == Float {

    /// Compute the conjugate of this quaternion
    @inlinable
    public var conjugate: Quaternion<Float> {
        Quaternion(vector: vector * SIMD4<Float>(-1, -1, -1, 1))
    }

    /// Compute the inverse of this quaternion
    @inlinable
    public var inverse: Quaternion<Float> {
        let lenSq = simd_length_squared(vector)
        if lenSq > Float.ulpOfOne {
            return Quaternion(vector: vector * SIMD4<Float>(-1, -1, -1, 1) / lenSq)
        }
        return .identity
    }

    /// Compute dot product of two quaternions
    @inlinable
    public func dot(_ other: Quaternion<Float>) -> Float {
        simd_dot(vector, other.vector)
    }

    /// Compute dot product using operator
    @inlinable
    public static func • (lhs: Quaternion<Float>, rhs: Quaternion<Float>) -> Float {
        simd_dot(lhs.vector, rhs.vector)
    }

    /// Rotate a position by this quaternion
    @inlinable
    public func rotate(_ position: Position<Float>) -> Position<Float> {
        Position(vector: Self._qrot(vector, position.vector))
    }

    /// Rotate a position using operator
    @inlinable
    public static func * (quaternion: Quaternion<Float>, position: Position<Float>) -> Position<Float> {
        quaternion.rotate(position)
    }
}
