import Foundation
import FoundationTypes
import simd

// MARK: - Arithmetic Operators for Double SpatialPose

extension SpatialPose where T == Double {

    /// Compose two poses (apply lhs then rhs)
    /// - Note: Pose composition is NOT commutative (p1 * p2 ≠ p2 * p1)
    @inlinable
    public static func * (lhs: SpatialPose<Double>, rhs: SpatialPose<Double>) -> SpatialPose<Double> {
        SpatialPose(
            position: lhs._pos + _qrot(lhs._rot, rhs._pos),
            rotation: _qmul(lhs._rot, rhs._rot)
        )
    }

    /// Transform a position using operator
    @inlinable
    public static func * (pose: SpatialPose<Double>, position: Position<Double>) -> Position<Double> {
        pose.transform(position)
    }

    // MARK: Position Offset

    /// Offset the pose position by a vector (no rotation change)
    @inlinable
    public static func + (lhs: SpatialPose<Double>, rhs: Position<Double>) -> SpatialPose<Double> {
        SpatialPose(position: lhs._pos + rhs.vector, rotation: lhs._rot)
    }

    /// Offset the pose position by a negative vector
    @inlinable
    public static func - (lhs: SpatialPose<Double>, rhs: Position<Double>) -> SpatialPose<Double> {
        SpatialPose(position: lhs._pos - rhs.vector, rotation: lhs._rot)
    }
}

// MARK: - Arithmetic Operators for Float SpatialPose

extension SpatialPose where T == Float {
    
    /// Compose two poses (apply lhs then rhs)
    /// - Note: Pose composition is NOT commutative (p1 * p2 ≠ p2 * p1)
    @inlinable
    public static func * (lhs: SpatialPose<Float>, rhs: SpatialPose<Float>) -> SpatialPose<Float> {
        SpatialPose(
            position: lhs._pos + _qrot(lhs._rot, rhs._pos),
            rotation: _qmul(lhs._rot, rhs._rot)
        )
    }


    /// Transform a position using operator
    @inlinable
    public static func * (pose: SpatialPose<Float>, position: Position<Float>) -> Position<Float> {
        pose.transform(position)
    }

    // MARK: Position Offset

    /// Offset the pose position by a vector (no rotation change)
    @inlinable
    public static func + (lhs: SpatialPose<Float>, rhs: Position<Float>) -> SpatialPose<Float> {
        SpatialPose(position: lhs._pos + rhs.vector, rotation: lhs._rot)
    }


    /// Offset the pose position by a negative vector
    @inlinable
    public static func - (lhs: SpatialPose<Float>, rhs: Position<Float>) -> SpatialPose<Float> {
        SpatialPose(position: lhs._pos - rhs.vector, rotation: lhs._rot)
    }
}
