import Foundation
import FoundationTypes
import simd

// MARK: - Vector Operators for Double SpatialPose

extension SpatialPose where T == Double {

    // MARK: Position Transformation

    /// Transform a position by this pose (rotate then translate)
    @inlinable
    public func transform(_ position: Position<Double>) -> Position<Double> {
        _pos + _rot.rotate(position)
    }

    // MARK: Pose Operations

    /// Compute the inverse of this pose
    @inlinable
    public var inverse: SpatialPose<Double> {
        let _rotInv = _rot.inverse
        return SpatialPose(
            position: _rotInv.rotate(-_pos),
            rotation: _rotInv
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
        _pos + _rot.rotate(position)
    }

    // MARK: Pose Operations

    /// Compute the inverse of this pose
    @inlinable
    public var inverse: SpatialPose<Float> {
        let _rotInv = _rot.inverse
        return SpatialPose(
            position: _rotInv.rotate(-_pos),
            rotation: _rotInv
        )
    }

    /// Compute the relative pose from this pose to another
    /// - Returns: The pose that transforms from self to other
    @inlinable
    public func relativePose(to other: SpatialPose<Float>) -> SpatialPose<Float> {
        self.inverse * other
    }
}
