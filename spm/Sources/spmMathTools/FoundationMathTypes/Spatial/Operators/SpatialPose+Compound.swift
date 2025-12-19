import Foundation
import FoundationTypes
import simd

// MARK: - Arithmetic Operators for Double SpatialPose

extension SpatialPose where T == Double {

    /// Compose this pose with another (in-place)
    @inlinable
    public static func *= (lhs: inout SpatialPose<Double>, rhs: SpatialPose<Double>) {
        lhs = lhs * rhs
    }

    /// Offset the pose position by a vector (in-place)
    @inlinable
    public static func += (lhs: inout SpatialPose<Double>, rhs: Position<Double>) {
        lhs._pos += rhs
    }

    /// Offset the pose position by a negative vector (in-place)
    @inlinable
    public static func -= (lhs: inout SpatialPose<Double>, rhs: Position<Double>) {
        lhs._pos -= rhs
    }
}

// MARK: - Compound Operators for Float SpatialPose

extension SpatialPose where T == Float {

    /// Compose this pose with another (in-place)
    @inlinable
    public static func *= (lhs: inout SpatialPose<Float>, rhs: SpatialPose<Float>) {
        lhs = lhs * rhs
    }

    /// Offset the pose position by a vector (in-place)
    @inlinable
    public static func += (lhs: inout SpatialPose<Float>, rhs: Position<Float>) {
        lhs._pos += rhs
    }

    /// Offset the pose position by a negative vector (in-place)
    @inlinable
    public static func -= (lhs: inout SpatialPose<Float>, rhs: Position<Float>) {
        lhs._pos -= rhs
    }
}
