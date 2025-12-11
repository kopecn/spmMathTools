import Foundation
import FoundationTypes
import simd

// MARK: - Arithmetic Operators for Double Quaternion

extension Quaternion where T == Double {

    // MARK: Multiplication (Quaternion Composition)

    /// Multiply two quaternions (quaternion composition/rotation chaining)
    /// - Note: Quaternion multiplication is NOT commutative (q1 * q2 ≠ q2 * q1)
    public static func * (lhs: Quaternion<Double>, rhs: Quaternion<Double>) -> Quaternion<Double> {
        let w = lhs.w * rhs.w - lhs.x * rhs.x - lhs.y * rhs.y - lhs.z * rhs.z
        let x = lhs.w * rhs.x + lhs.x * rhs.w + lhs.y * rhs.z - lhs.z * rhs.y
        let y = lhs.w * rhs.y - lhs.x * rhs.z + lhs.y * rhs.w + lhs.z * rhs.x
        let z = lhs.w * rhs.z + lhs.x * rhs.y - lhs.y * rhs.x + lhs.z * rhs.w
        return Quaternion(x: x, y: y, z: z, w: w)
    }

    /// Multiply quaternion by scalar (scaling)
    public static func * (lhs: Quaternion<Double>, rhs: Double) -> Quaternion<Double> {
        return Quaternion(vector: lhs.vector * rhs)
    }

    /// Multiply scalar by quaternion (scaling)
    public static func * (lhs: Double, rhs: Quaternion<Double>) -> Quaternion<Double> {
        return Quaternion(vector: lhs * rhs.vector)
    }

    // MARK: Division

    /// Divide quaternion by scalar
    public static func / (lhs: Quaternion<Double>, rhs: Double) -> Quaternion<Double> {
        return Quaternion(vector: lhs.vector / rhs)
    }

    // MARK: Addition

    /// Add two quaternions (component-wise, useful for interpolation)
    public static func + (lhs: Quaternion<Double>, rhs: Quaternion<Double>) -> Quaternion<Double> {
        return Quaternion(vector: lhs.vector + rhs.vector)
    }

    /// Add scalar to all components
    public static func + (lhs: Quaternion<Double>, rhs: Double) -> Quaternion<Double> {
        return Quaternion(vector: lhs.vector + SIMD4<Double>(repeating: rhs))
    }

    /// Add quaternion to scalar
    public static func + (lhs: Double, rhs: Quaternion<Double>) -> Quaternion<Double> {
        return rhs + lhs
    }

    // MARK: Subtraction

    /// Subtract two quaternions (component-wise, useful for interpolation)
    public static func - (lhs: Quaternion<Double>, rhs: Quaternion<Double>) -> Quaternion<Double> {
        return Quaternion(vector: lhs.vector - rhs.vector)
    }

    /// Subtract scalar from all components
    public static func - (lhs: Quaternion<Double>, rhs: Double) -> Quaternion<Double> {
        return Quaternion(vector: lhs.vector - SIMD4<Double>(repeating: rhs))
    }

    /// Subtract quaternion from scalar
    public static func - (lhs: Double, rhs: Quaternion<Double>) -> Quaternion<Double> {
        return Quaternion(vector: SIMD4<Double>(repeating: lhs) - rhs.vector)
    }
}

// MARK: - Arithmetic Operators for Float Quaternion

extension Quaternion where T == Float {

    // MARK: Multiplication (Quaternion Composition)

    /// Multiply two quaternions (quaternion composition/rotation chaining)
    /// - Note: Quaternion multiplication is NOT commutative (q1 * q2 ≠ q2 * q1)
    public static func * (lhs: Quaternion<Float>, rhs: Quaternion<Float>) -> Quaternion<Float> {
        let w = lhs.w * rhs.w - lhs.x * rhs.x - lhs.y * rhs.y - lhs.z * rhs.z
        let x = lhs.w * rhs.x + lhs.x * rhs.w + lhs.y * rhs.z - lhs.z * rhs.y
        let y = lhs.w * rhs.y - lhs.x * rhs.z + lhs.y * rhs.w + lhs.z * rhs.x
        let z = lhs.w * rhs.z + lhs.x * rhs.y - lhs.y * rhs.x + lhs.z * rhs.w
        return Quaternion(x: x, y: y, z: z, w: w)
    }

    /// Multiply quaternion by scalar (scaling)
    public static func * (lhs: Quaternion<Float>, rhs: Float) -> Quaternion<Float> {
        return Quaternion(vector: lhs.vector * rhs)
    }

    /// Multiply scalar by quaternion (scaling)
    public static func * (lhs: Float, rhs: Quaternion<Float>) -> Quaternion<Float> {
        return Quaternion(vector: lhs * rhs.vector)
    }

    // MARK: Division

    /// Divide quaternion by scalar
    public static func / (lhs: Quaternion<Float>, rhs: Float) -> Quaternion<Float> {
        return Quaternion(vector: lhs.vector / rhs)
    }

    // MARK: Addition

    /// Add two quaternions (component-wise, useful for interpolation)
    public static func + (lhs: Quaternion<Float>, rhs: Quaternion<Float>) -> Quaternion<Float> {
        return Quaternion(vector: lhs.vector + rhs.vector)
    }

    /// Add scalar to all components
    public static func + (lhs: Quaternion<Float>, rhs: Float) -> Quaternion<Float> {
        return Quaternion(vector: lhs.vector + SIMD4<Float>(repeating: rhs))
    }

    /// Add quaternion to scalar
    public static func + (lhs: Float, rhs: Quaternion<Float>) -> Quaternion<Float> {
        return rhs + lhs
    }

    // MARK: Subtraction

    /// Subtract two quaternions (component-wise, useful for interpolation)
    public static func - (lhs: Quaternion<Float>, rhs: Quaternion<Float>) -> Quaternion<Float> {
        return Quaternion(vector: lhs.vector - rhs.vector)
    }

    /// Subtract scalar from all components
    public static func - (lhs: Quaternion<Float>, rhs: Float) -> Quaternion<Float> {
        return Quaternion(vector: lhs.vector - SIMD4<Float>(repeating: rhs))
    }

    /// Subtract quaternion from scalar
    public static func - (lhs: Float, rhs: Quaternion<Float>) -> Quaternion<Float> {
        return Quaternion(vector: SIMD4<Float>(repeating: lhs) - rhs.vector)
    }
}
