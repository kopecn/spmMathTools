import Foundation
import FoundationTypes
import simd

// MARK: - Arithmetic Operators for Double SpatialPose

extension SpatialPose where T == Double {

    // MARK: Pose Composition

    /// Compose two poses (apply lhs then rhs)
    /// - Note: Pose composition is NOT commutative (p1 * p2 ≠ p2 * p1)
    public static func * (lhs: SpatialPose<Double>, rhs: SpatialPose<Double>) -> SpatialPose<Double> {
        let newRotation = lhs.quaternion * rhs.quaternion
        let newPosition = lhs.quaternion.rotate(rhs.position) + lhs.position
        return SpatialPose(position: newPosition, rotation: newRotation)
    }

    /// Compose this pose with another (in-place)
    public static func *= (lhs: inout SpatialPose<Double>, rhs: SpatialPose<Double>) {
        lhs = lhs * rhs
    }

    // MARK: Position Transformation

    /// Transform a position by this pose (rotate then translate)
    public func transform(_ position: Position<Double>) -> Position<Double> {
        return quaternion.rotate(position) + self.position
    }

    /// Transform a position using operator
    public static func * (pose: SpatialPose<Double>, position: Position<Double>) -> Position<Double> {
        return pose.transform(position)
    }

    // MARK: Position Offset

    /// Offset the pose position by a vector (no rotation change)
    public static func + (lhs: SpatialPose<Double>, rhs: Position<Double>) -> SpatialPose<Double> {
        return SpatialPose(position: lhs.position + rhs, rotation: lhs.quaternion)
    }

    /// Offset the pose position by a vector (in-place)
    public static func += (lhs: inout SpatialPose<Double>, rhs: Position<Double>) {
        lhs = SpatialPose(position: lhs.position + rhs, rotation: lhs.quaternion)
    }

    /// Offset the pose position by a negative vector
    public static func - (lhs: SpatialPose<Double>, rhs: Position<Double>) -> SpatialPose<Double> {
        return SpatialPose(position: lhs.position - rhs, rotation: lhs.quaternion)
    }

    /// Offset the pose position by a negative vector (in-place)
    public static func -= (lhs: inout SpatialPose<Double>, rhs: Position<Double>) {
        lhs = SpatialPose(position: lhs.position - rhs, rotation: lhs.quaternion)
    }

    // MARK: Pose Operations

    /// Compute the inverse of this pose
    public var inverse: SpatialPose<Double> {
        let invRotation = quaternion.inverse
        let invPosition = invRotation.rotate(-position)
        return SpatialPose(position: invPosition, rotation: invRotation)
    }

    /// Compute the relative pose from this pose to another
    /// - Returns: The pose that transforms from self to other
    public func relativePose(to other: SpatialPose<Double>) -> SpatialPose<Double> {
        return self.inverse * other
    }

    /// Interpolate between two poses using SLERP for rotation and linear for position
    /// - Parameter other: Target pose
    /// - Parameter t: Interpolation parameter (0 = self, 1 = other)
    /// - Returns: Interpolated pose
    public func interpolate(to other: SpatialPose<Double>, t: Double) -> SpatialPose<Double> {
        // Linear interpolation for position
        let interpPosition = position + (other.position - position) * t

        // SLERP for rotation
        let dot = quaternion.dot(other.quaternion)
        let adjustedOther = dot < 0 ? -other.quaternion : other.quaternion
        let dotClamped = min(max(abs(dot), -1.0), 1.0)

        let theta = acos(dotClamped)
        if theta < Double.ulpOfOne {
            // Quaternions are too close, use linear interpolation
            let interpRotation = (quaternion + (adjustedOther - quaternion) * t)
            let magnitude = simd_length(interpRotation.vector)
            return SpatialPose(position: interpPosition, rotation: Quaternion(vector: interpRotation.vector / magnitude))
        } else {
            let sinTheta = sin(theta)
            let w1 = sin((1.0 - t) * theta) / sinTheta
            let w2 = sin(t * theta) / sinTheta
            let interpRotation = quaternion * w1 + adjustedOther * w2
            return SpatialPose(position: interpPosition, rotation: interpRotation)
        }
    }

    // MARK: Distance

    /// Compute the positional distance between two poses
    public func positionDistance(to other: SpatialPose<Double>) -> Double {
        return position.distance(to: other.position)
    }

    /// Compute the squared positional distance between two poses (more efficient)
    public func positionDistanceSquared(to other: SpatialPose<Double>) -> Double {
        return position.distanceSquared(to: other.position)
    }

    /// Compute the angular distance between two poses (in radians)
    public func angularDistance(to other: SpatialPose<Double>) -> Double {
        let relativeRotation = quaternion.inverse * other.quaternion
        return 2.0 * acos(min(max(abs(relativeRotation.w), -1.0), 1.0))
    }
}

// MARK: - Arithmetic Operators for Float SpatialPose

extension SpatialPose where T == Float {

    // MARK: Pose Composition

    /// Compose two poses (apply lhs then rhs)
    /// - Note: Pose composition is NOT commutative (p1 * p2 ≠ p2 * p1)
    public static func * (lhs: SpatialPose<Float>, rhs: SpatialPose<Float>) -> SpatialPose<Float> {
        let newRotation = lhs.quaternion * rhs.quaternion
        let newPosition = lhs.quaternion.rotate(rhs.position) + lhs.position
        return SpatialPose(position: newPosition, rotation: newRotation)
    }

    /// Compose this pose with another (in-place)
    public static func *= (lhs: inout SpatialPose<Float>, rhs: SpatialPose<Float>) {
        lhs = lhs * rhs
    }

    // MARK: Position Transformation

    /// Transform a position by this pose (rotate then translate)
    public func transform(_ position: Position<Float>) -> Position<Float> {
        return quaternion.rotate(position) + self.position
    }

    /// Transform a position using operator
    public static func * (pose: SpatialPose<Float>, position: Position<Float>) -> Position<Float> {
        return pose.transform(position)
    }

    // MARK: Position Offset

    /// Offset the pose position by a vector (no rotation change)
    public static func + (lhs: SpatialPose<Float>, rhs: Position<Float>) -> SpatialPose<Float> {
        return SpatialPose(position: lhs.position + rhs, rotation: lhs.quaternion)
    }

    /// Offset the pose position by a vector (in-place)
    public static func += (lhs: inout SpatialPose<Float>, rhs: Position<Float>) {
        lhs = SpatialPose(position: lhs.position + rhs, rotation: lhs.quaternion)
    }

    /// Offset the pose position by a negative vector
    public static func - (lhs: SpatialPose<Float>, rhs: Position<Float>) -> SpatialPose<Float> {
        return SpatialPose(position: lhs.position - rhs, rotation: lhs.quaternion)
    }

    /// Offset the pose position by a negative vector (in-place)
    public static func -= (lhs: inout SpatialPose<Float>, rhs: Position<Float>) {
        lhs = SpatialPose(position: lhs.position - rhs, rotation: lhs.quaternion)
    }

    // MARK: Pose Operations

    /// Compute the inverse of this pose
    public var inverse: SpatialPose<Float> {
        let invRotation = quaternion.inverse
        let invPosition = invRotation.rotate(-position)
        return SpatialPose(position: invPosition, rotation: invRotation)
    }

    /// Compute the relative pose from this pose to another
    /// - Returns: The pose that transforms from self to other
    public func relativePose(to other: SpatialPose<Float>) -> SpatialPose<Float> {
        return self.inverse * other
    }

    /// Interpolate between two poses using SLERP for rotation and linear for position
    /// - Parameter other: Target pose
    /// - Parameter t: Interpolation parameter (0 = self, 1 = other)
    /// - Returns: Interpolated pose
    public func interpolate(to other: SpatialPose<Float>, t: Float) -> SpatialPose<Float> {
        // Linear interpolation for position
        let interpPosition = position + (other.position - position) * t

        // SLERP for rotation
        let dot = quaternion.dot(other.quaternion)
        let adjustedOther = dot < 0 ? -other.quaternion : other.quaternion
        let dotClamped = min(max(abs(dot), -1.0), 1.0)

        let theta = acos(dotClamped)
        if theta < Float.ulpOfOne {
            // Quaternions are too close, use linear interpolation
            let interpRotation = (quaternion + (adjustedOther - quaternion) * t)
            let magnitude = simd_length(interpRotation.vector)
            return SpatialPose(position: interpPosition, rotation: Quaternion(vector: interpRotation.vector / magnitude))
        } else {
            let sinTheta = sin(theta)
            let w1 = sin((1.0 - t) * theta) / sinTheta
            let w2 = sin(t * theta) / sinTheta
            let interpRotation = quaternion * w1 + adjustedOther * w2
            return SpatialPose(position: interpPosition, rotation: interpRotation)
        }
    }

    // MARK: Distance

    /// Compute the positional distance between two poses
    public func positionDistance(to other: SpatialPose<Float>) -> Float {
        return position.distance(to: other.position)
    }

    /// Compute the squared positional distance between two poses (more efficient)
    public func positionDistanceSquared(to other: SpatialPose<Float>) -> Float {
        return position.distanceSquared(to: other.position)
    }

    /// Compute the angular distance between two poses (in radians)
    public func angularDistance(to other: SpatialPose<Float>) -> Float {
        let relativeRotation = quaternion.inverse * other.quaternion
        return 2.0 * acos(min(max(abs(relativeRotation.w), -1.0), 1.0))
    }
}
