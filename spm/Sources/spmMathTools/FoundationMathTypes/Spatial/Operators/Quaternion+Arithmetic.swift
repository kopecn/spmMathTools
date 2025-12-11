import Foundation
import FoundationTypes
import simd

// MARK: - Arithmetic Operators for Double Quaternion

extension Quaternion where T == Double {

    // MARK: Multiplication (Quaternion Composition)

    /// Multiply two quaternions (quaternion composition/rotation chaining)
    /// - Note: Quaternion multiplication is NOT commutative (q1 * q2 ≠ q2 * q1)
    @inlinable
    public static func * (lhs: Quaternion<Double>, rhs: Quaternion<Double>) -> Quaternion<Double> {
        Quaternion(vector: _qmul(lhs.vector, rhs.vector))
    }

    /// Multiply quaternion by scalar (scaling)
    @inlinable
    public static func * (lhs: Quaternion<Double>, rhs: Double) -> Quaternion<Double> {
        Quaternion(vector: lhs.vector * rhs)
    }

    /// Multiply scalar by quaternion (scaling)
    @inlinable
    public static func * (lhs: Double, rhs: Quaternion<Double>) -> Quaternion<Double> {
        Quaternion(vector: lhs * rhs.vector)
    }

    // MARK: Division

    /// Divide quaternion by scalar
    @inlinable
    public static func / (lhs: Quaternion<Double>, rhs: Double) -> Quaternion<Double> {
        Quaternion(vector: lhs.vector / rhs)
    }

    // MARK: Addition

    /// Add two quaternions (component-wise, useful for interpolation)
    @inlinable
    public static func + (lhs: Quaternion<Double>, rhs: Quaternion<Double>) -> Quaternion<Double> {
        Quaternion(vector: lhs.vector + rhs.vector)
    }

    /// Add scalar to all components
    @inlinable
    public static func + (lhs: Quaternion<Double>, rhs: Double) -> Quaternion<Double> {
        Quaternion(vector: lhs.vector + SIMD4<Double>(repeating: rhs))
    }

    /// Add quaternion to scalar
    @inlinable
    public static func + (lhs: Double, rhs: Quaternion<Double>) -> Quaternion<Double> {
        rhs + lhs
    }

    // MARK: Subtraction

    /// Subtract two quaternions (component-wise, useful for interpolation)
    @inlinable
    public static func - (lhs: Quaternion<Double>, rhs: Quaternion<Double>) -> Quaternion<Double> {
        Quaternion(vector: lhs.vector - rhs.vector)
    }

    /// Subtract scalar from all components
    @inlinable
    public static func - (lhs: Quaternion<Double>, rhs: Double) -> Quaternion<Double> {
        Quaternion(vector: lhs.vector - SIMD4<Double>(repeating: rhs))
    }

    /// Subtract quaternion from scalar
    @inlinable
    public static func - (lhs: Double, rhs: Quaternion<Double>) -> Quaternion<Double> {
        Quaternion(vector: SIMD4<Double>(repeating: lhs) - rhs.vector)
    }
}

// MARK: - Arithmetic Operators for Float Quaternion

extension Quaternion where T == Float {

    // MARK: Multiplication (Quaternion Composition)

    /// Multiply two quaternions (quaternion composition/rotation chaining)
    /// - Note: Quaternion multiplication is NOT commutative (q1 * q2 ≠ q2 * q1)
    @inlinable
    public static func * (lhs: Quaternion<Float>, rhs: Quaternion<Float>) -> Quaternion<Float> {
        Quaternion(vector: _qmul(lhs.vector, rhs.vector))
    }

    /// Multiply quaternion by scalar (scaling)
    @inlinable
    public static func * (lhs: Quaternion<Float>, rhs: Float) -> Quaternion<Float> {
        Quaternion(vector: lhs.vector * rhs)
    }

    /// Multiply scalar by quaternion (scaling)
    @inlinable
    public static func * (lhs: Float, rhs: Quaternion<Float>) -> Quaternion<Float> {
        Quaternion(vector: lhs * rhs.vector)
    }

    // MARK: Division

    /// Divide quaternion by scalar
    @inlinable
    public static func / (lhs: Quaternion<Float>, rhs: Float) -> Quaternion<Float> {
        Quaternion(vector: lhs.vector / rhs)
    }

    // MARK: Addition

    /// Add two quaternions (component-wise, useful for interpolation)
    @inlinable
    public static func + (lhs: Quaternion<Float>, rhs: Quaternion<Float>) -> Quaternion<Float> {
        Quaternion(vector: lhs.vector + rhs.vector)
    }

    /// Add scalar to all components
    @inlinable
    public static func + (lhs: Quaternion<Float>, rhs: Float) -> Quaternion<Float> {
        Quaternion(vector: lhs.vector + SIMD4<Float>(repeating: rhs))
    }

    /// Add quaternion to scalar
    @inlinable
    public static func + (lhs: Float, rhs: Quaternion<Float>) -> Quaternion<Float> {
        rhs + lhs
    }

    // MARK: Subtraction

    /// Subtract two quaternions (component-wise, useful for interpolation)
    @inlinable
    public static func - (lhs: Quaternion<Float>, rhs: Quaternion<Float>) -> Quaternion<Float> {
        Quaternion(vector: lhs.vector - rhs.vector)
    }

    /// Subtract scalar from all components
    @inlinable
    public static func - (lhs: Quaternion<Float>, rhs: Float) -> Quaternion<Float> {
        Quaternion(vector: lhs.vector - SIMD4<Float>(repeating: rhs))
    }

    /// Subtract quaternion from scalar
    @inlinable
    public static func - (lhs: Float, rhs: Quaternion<Float>) -> Quaternion<Float> {
        Quaternion(vector: SIMD4<Float>(repeating: lhs) - rhs.vector)
    }
}
