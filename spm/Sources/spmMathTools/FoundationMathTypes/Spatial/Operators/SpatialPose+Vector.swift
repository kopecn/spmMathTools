import Foundation
import FoundationTypes
import simd

// MARK: - Vector Operators for Double SpatialPose

extension SpatialPose where T == Double {

    // MARK: Position Transformation

    /// Transform a position by this pose (rotate then translate)
    @inlinable
    public func transform(_ position: Position<Double>) -> Position<Double> {
        Position<Double>(vector: _pos + Self._qrot(_rot, position.vector))
    }

    // MARK: Pose Operations

    /// Compute the inverse of this pose
    @inlinable
    public var inverse: SpatialPose<Double> {
        SpatialPose(
            position: Self._qrot(SIMD4<Double>(-_rot.x, -_rot.y, -_rot.z, _rot.w), -_pos),
            rotation: SIMD4<Double>(-_rot.x, -_rot.y, -_rot.z, _rot.w)
        )
    }

    /// Compute the relative pose from this pose to another
    /// - Returns: The pose that transforms from self to other
    @inlinable
    public func relativePose(to other: SpatialPose<Double>) -> SpatialPose<Double> {
        self.inverse * other
    }
}

// MARK: - Vector Operators for Float SpatialPose

extension SpatialPose where T == Float {

    // MARK: Position Transformation

    /// Transform a position by this pose (rotate then translate)
    @inlinable
    public func transform(_ position: Position<Float>) -> Position<Float> {
        Position<Float>(vector: _pos + Self._qrot(_rot, position.vector))
    }

    // MARK: Pose Operations

    /// Compute the inverse of this pose
    @inlinable
    public var inverse: SpatialPose<Float> {
        SpatialPose(
            position: Self._qrot(SIMD4<Float>(-_rot.x, -_rot.y, -_rot.z, _rot.w), -_pos),
            rotation: SIMD4<Float>(-_rot.x, -_rot.y, -_rot.z, _rot.w)
        )
    }

    /// Compute the relative pose from this pose to another
    /// - Returns: The pose that transforms from self to other
    @inlinable
    public func relativePose(to other: SpatialPose<Float>) -> SpatialPose<Float> {
        self.inverse * other
    }
}
