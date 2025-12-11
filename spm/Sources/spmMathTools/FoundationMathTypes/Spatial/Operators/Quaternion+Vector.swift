import Foundation
import FoundationTypes
import simd

// MARK: - Quaternion Operators for Double Quaternion

extension Quaternion where T == Double {

    /// Compute the conjugate of this quaternion
    public var conjugate: Quaternion<Double> {
        return Quaternion(x: -x, y: -y, z: -z, w: w)
    }

    /// Compute the inverse of this quaternion
    public var inverse: Quaternion<Double> {
        let magSq = simd_length_squared(vector)
        if magSq > Double.ulpOfOne {
            let conj = conjugate
            return Quaternion(vector: conj.vector / magSq)
        }
        return .identity
    }

    /// Compute dot product of two quaternions
    public func dot(_ other: Quaternion<Double>) -> Double {
        return simd_dot(self.vector, other.vector)
    }

    /// Compute dot product using operator
    public static func • (lhs: Quaternion<Double>, rhs: Quaternion<Double>) -> Double {
        return simd_dot(lhs.vector, rhs.vector)
    }

    /// Rotate a position by this quaternion
    public func rotate(_ position: Position<Double>) -> Position<Double> {
        let qv = SIMD3<Double>(x, y, z)
        let uv = simd_cross(qv, position.vector)
        let uuv = simd_cross(qv, uv)
        return Position(vector: position.vector + ((uv * w) + uuv) * 2)
    }

    /// Rotate a position using operator
    public static func * (quaternion: Quaternion<Double>, position: Position<Double>) -> Position<Double> {
        return quaternion.rotate(position)
    }
}

// MARK: - Quaternion Operators for Float Quaternion

extension Quaternion where T == Float {
    
    /// Compute the conjugate of this quaternion
    public var conjugate: Quaternion<Float> {
        return Quaternion(x: -x, y: -y, z: -z, w: w)
    }

    /// Compute the inverse of this quaternion
    public var inverse: Quaternion<Float> {
        let magSq = simd_length_squared(vector)
        if magSq > Float.ulpOfOne {
            let conj = conjugate
            return Quaternion(vector: conj.vector / magSq)
        }
        return .identity
    }

    /// Compute dot product of two quaternions
    public func dot(_ other: Quaternion<Float>) -> Float {
        return simd_dot(self.vector, other.vector)
    }

    /// Compute dot product using operator
    public static func • (lhs: Quaternion<Float>, rhs: Quaternion<Float>) -> Float {
        return simd_dot(lhs.vector, rhs.vector)
    }

    /// Rotate a position by this quaternion
    public func rotate(_ position: Position<Float>) -> Position<Float> {
        let qv = SIMD3<Float>(x, y, z)
        let uv = simd_cross(qv, position.vector)
        let uuv = simd_cross(qv, uv)
        return Position(vector: position.vector + ((uv * w) + uuv) * 2)
    }

    /// Rotate a position using operator
    public static func * (quaternion: Quaternion<Float>, position: Position<Float>) -> Position<Float> {
        return quaternion.rotate(position)
    }
}
