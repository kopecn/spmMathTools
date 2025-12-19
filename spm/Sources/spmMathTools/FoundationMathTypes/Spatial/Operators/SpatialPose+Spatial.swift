import Foundation
import FoundationTypes
import simd

// MARK: - Spatial Operators for Double SpatialPose

extension SpatialPose where T == Double {

    /// Interpolate between two poses using SLERP for rotation and linear for position
    /// - Parameter other: Target pose
    /// - Parameter t: Interpolation parameter (0 = self, 1 = other)
    /// - Returns: Interpolated pose
    public func interpolate(to other: SpatialPose<Double>, t: Double) -> SpatialPose<Double> {
        SpatialPose(
            position: _pos + (other._pos - _pos) * t,
            rotation: {
                if (_rot • other._rot) < 0 {
                    if acos(min(max( _rot • -other._rot, -1.0), 1.0)) < Double.ulpOfOne {
                        return (_rot + (-other._rot - _rot).normalized * t)
                    } else {
                        return _rot
                            * (sin((1.0 - t) * acos(min(max( _rot • -other._rot , -1.0), 1.0)))
                                / sin(acos(min(max( _rot • -other._rot , -1.0), 1.0)))) + -other._rot
                            * (sin(t * acos(min(max( _rot • -other._rot, -1.0), 1.0)))
                                / sin(acos(min(max( _rot • -other._rot, -1.0), 1.0))))
                    }
                } else {
                    if acos(min(max(_rot • other._rot, -1.0), 1.0)) < Double.ulpOfOne {
                        return (_rot + (other._rot - _rot).normalized * t)
                    } else {
                        return _rot
                            * (sin((1.0 - t) * acos(min(max( _rot • other._rot, -1.0), 1.0)))
                                / sin(acos(min(max( _rot • other._rot, -1.0), 1.0)))) + other._rot
                            * (sin(t * acos(min(max( _rot • other._rot, -1.0), 1.0)))
                                / sin(acos(min(max( _rot • other._rot, -1.0), 1.0))))
                    }
                }
            }()
        )
    }

    // MARK: Distance

    /// Compute the positional distance between two poses
    @inlinable
    public func positionDistance(to other: SpatialPose<Double>) -> Double {
        _pos.distance(to: other._pos)
    }

    /// Compute the squared positional distance between two poses (more efficient)
    @inlinable
    public func positionDistanceSquared(to other: SpatialPose<Double>) -> Double {
        _pos.distanceSquared(to: other._pos)
    }

    /// Compute the angular distance between two poses (in radians)
    @inlinable
    public func angularDistance(to other: SpatialPose<Double>) -> Double {
        2.0 * acos(min(max(abs(Self._qmul(SIMD4<Double>(-_rot.x, -_rot.y, -_rot.z, _rot.w), other._rot.vector).w), -1.0), 1.0))
    }
}

// MARK: - Spatial Operators for Float SpatialPose

extension SpatialPose where T == Float {

    /// Interpolate between two poses using SLERP for rotation and linear for position
    /// - Parameter other: Target pose
    /// - Parameter t: Interpolation parameter (0 = self, 1 = other)
    /// - Returns: Interpolated pose
    public func interpolate(to other: SpatialPose<Float>, t: Float) -> SpatialPose<Float> {
        SpatialPose(
            position: _pos + (other._pos - _pos) * t,
            rotation: {
                if (_rot • other._rot) < 0 {
                    if acos(min(max( _rot • -other._rot, -1.0), 1.0)) < Float.ulpOfOne {
                        return (_rot + (-other._rot - _rot).normalized * t)
                    } else {
                        return _rot
                            * (sin((1.0 - t) * acos(min(max( _rot • -other._rot , -1.0), 1.0)))
                                / sin(acos(min(max( _rot • -other._rot , -1.0), 1.0)))) + -other._rot
                            * (sin(t * acos(min(max( _rot • -other._rot, -1.0), 1.0)))
                                / sin(acos(min(max( _rot • -other._rot, -1.0), 1.0))))
                    }
                } else {
                    if acos(min(max(_rot • other._rot, -1.0), 1.0)) < Float.ulpOfOne {
                        return (_rot + (other._rot - _rot).normalized * t)
                    } else {
                        return _rot
                            * (sin((1.0 - t) * acos(min(max( _rot • other._rot, -1.0), 1.0)))
                                / sin(acos(min(max( _rot • other._rot, -1.0), 1.0)))) + other._rot
                            * (sin(t * acos(min(max( _rot • other._rot, -1.0), 1.0)))
                                / sin(acos(min(max( _rot • other._rot, -1.0), 1.0))))
                    }
                }
            }()
        )
    }

    // MARK: Distance

    /// Compute the positional distance between two poses
    @inlinable
    public func positionDistance(to other: SpatialPose<Float>) -> Float {
        _pos.distance(to: other._pos)
    }

    /// Compute the squared positional distance between two poses (more efficient)
    @inlinable
    public func positionDistanceSquared(to other: SpatialPose<Float>) -> Float {
        _pos.distanceSquared(to: other._pos)
    }

    /// Compute the angular distance between two poses (in radians)
    @inlinable
    public func angularDistance(to other: SpatialPose<Float>) -> Float {
        2.0 * acos(min(max(abs(Self._qmul(SIMD4<Float>(-_rot.x, -_rot.y, -_rot.z, _rot.w), other._rot.vector).w), -1.0), 1.0))
    }
}
