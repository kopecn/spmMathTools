import Foundation
import FoundationTypes
import simd

// MARK: - Compound Operators for Double Quaternion

extension Quaternion where T == Double {

    /// Multiply this quaternion by another (in-place composition)
    public static func *= (lhs: inout Quaternion<Double>, rhs: Quaternion<Double>) {
        lhs = lhs * rhs
    }

    /// Multiply this quaternion by a scalar (in-place)
    public static func *= (lhs: inout Quaternion<Double>, rhs: Double) {
        lhs.vector *= rhs
    }

    /// Divide this quaternion by a scalar (in-place)
    public static func /= (lhs: inout Quaternion<Double>, rhs: Double) {
        lhs.vector /= rhs
    }

    /// Add another quaternion to this quaternion (in-place)
    public static func += (lhs: inout Quaternion<Double>, rhs: Quaternion<Double>) {
        lhs.vector += rhs.vector
    }

    /// Add scalar to all components (in-place)
    public static func += (lhs: inout Quaternion<Double>, rhs: Double) {
        lhs.vector += SIMD4<Double>(repeating: rhs)
    }

    /// Subtract another quaternion from this quaternion (in-place)
    public static func -= (lhs: inout Quaternion<Double>, rhs: Quaternion<Double>) {
        lhs.vector -= rhs.vector
    }

    /// Subtract scalar from all components (in-place)
    public static func -= (lhs: inout Quaternion<Double>, rhs: Double) {
        lhs.vector -= SIMD4<Double>(repeating: rhs)
    }
}

// MARK: - Compound Operators for Float Quaternion

extension Quaternion where T == Float {

    /// Multiply this quaternion by another (in-place composition)
    public static func *= (lhs: inout Quaternion<Float>, rhs: Quaternion<Float>) {
        lhs = lhs * rhs
    }

    /// Multiply this quaternion by a scalar (in-place)
    public static func *= (lhs: inout Quaternion<Float>, rhs: Float) {
        lhs.vector *= rhs
    }

    /// Divide this quaternion by a scalar (in-place)
    public static func /= (lhs: inout Quaternion<Float>, rhs: Float) {
        lhs.vector /= rhs
    }

    /// Add another quaternion to this quaternion (in-place)
    public static func += (lhs: inout Quaternion<Float>, rhs: Quaternion<Float>) {
        lhs.vector += rhs.vector
    }

    /// Add scalar to all components (in-place)
    public static func += (lhs: inout Quaternion<Float>, rhs: Float) {
        lhs.vector += SIMD4<Float>(repeating: rhs)
    }

    /// Subtract another quaternion from this quaternion (in-place)
    public static func -= (lhs: inout Quaternion<Float>, rhs: Quaternion<Float>) {
        lhs.vector -= rhs.vector
    }

    /// Subtract scalar from all components (in-place)
    public static func -= (lhs: inout Quaternion<Float>, rhs: Float) {
        lhs.vector -= SIMD4<Float>(repeating: rhs)
    }
}
