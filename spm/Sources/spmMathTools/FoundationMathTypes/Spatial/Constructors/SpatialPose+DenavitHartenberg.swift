import Foundation
import simd
import FoundationTypes

// MARK: - Denavit-Hartenberg Initializers

extension SpatialPose where T == Double {
    /// Initialize a SpatialPose using standard Denavit-Hartenberg parameters.
    ///
    /// This is an optimized initializer that directly computes position and quaternion
    /// without creating intermediate transformation matrices.
    ///
    /// [See DH Parameters](https://en.wikipedia.org/wiki/Denavit–Hartenberg_parameters)
    ///
    /// - Parameters:
    ///   - a: Link length (distance along x axis)
    ///   - alpha: Link twist (angle in radians around x axis)
    ///   - d: Link offset (distance along z axis)
    ///   - theta: Joint angle (angle in radians around z axis)
    @inlinable
    public init(denavitHartenberg a: Double, alpha: Double, d: Double, theta: Double) {
        // Precompute trigonometric values
        let ct = cos(theta)
        let st = sin(theta)
        let ca = cos(alpha)
        let sa = sin(alpha)

        // Position from DH transformation: (a*cos(θ), a*sin(θ), d)
        let pos = SIMD3<Double>(a * ct, a * st, d)

        // Convert rotation matrix to quaternion
        // Matrix: [ ct, -st*ca,  st*sa ]
        //         [ st,  ct*ca, -ct*sa ]
        //         [  0,     sa,     ca ]
        let trace = ct + ct * ca + ca
        let rot: SIMD4<Double>

        if trace > 0 {
            let s = sqrt(trace + 1.0) * 2.0
            rot = SIMD4<Double>(
                (sa - (-ct * sa)) / s,
                (st * sa - 0.0) / s,
                (st - (-st * ca)) / s,
                0.25 * s
            )
        } else if ct > ct * ca && ct > ca {
            let s = sqrt(1.0 + ct - ct * ca - ca) * 2.0
            rot = SIMD4<Double>(
                0.25 * s,
                (-st * ca + st) / s,
                (st * sa + 0.0) / s,
                (sa - (-ct * sa)) / s
            )
        } else if ct * ca > ca {
            let s = sqrt(1.0 + ct * ca - ct - ca) * 2.0
            rot = SIMD4<Double>(
                (-st * ca + st) / s,
                0.25 * s,
                (-ct * sa + sa) / s,
                (st * sa - 0.0) / s
            )
        } else {
            let s = sqrt(1.0 + ca - ct - ct * ca) * 2.0
            rot = SIMD4<Double>(
                (st * sa + 0.0) / s,
                (-ct * sa + sa) / s,
                0.25 * s,
                (st - (-st * ca)) / s
            )
        }

        self.init(position: pos, rotation: rot)
    }

    /// Initialize a SpatialPose using Denavit-Hartenberg parameters with precomputed link twist.
    ///
    /// Use this optimized version when the link twist (alpha) is constant and its sine/cosine
    /// can be precomputed. This saves two trigonometric function calls per invocation.
    ///
    /// [See DH Parameters](https://en.wikipedia.org/wiki/Denavit–Hartenberg_parameters)
    ///
    /// - Parameters:
    ///   - a: Link length (distance along x axis)
    ///   - ca: Cosine of the link twist
    ///   - sa: Sine of the link twist
    ///   - d: Link offset (distance along z axis)
    ///   - theta: Joint angle (angle in radians around z axis)
    @inlinable
    public init(denavitHartenberg a: Double, ca: Double, sa: Double, d: Double, theta: Double) {
        // Compute trigonometric values for theta only
        let ct = cos(theta)
        let st = sin(theta)

        // Position from DH transformation
        let pos = SIMD3<Double>(a * ct, a * st, d)

        // Convert to quaternion
        // Matrix: [ ct, -st*ca,  st*sa ]
        //         [ st,  ct*ca, -ct*sa ]
        //         [  0,     sa,     ca ]
        let trace = ct + ct * ca + ca
        let rot: SIMD4<Double>

        if trace > 0 {
            let s = sqrt(trace + 1.0) * 2.0
            rot = SIMD4<Double>(
                (sa - (-ct * sa)) / s,
                (st * sa - 0.0) / s,
                (st - (-st * ca)) / s,
                0.25 * s
            )
        } else if ct > ct * ca && ct > ca {
            let s = sqrt(1.0 + ct - ct * ca - ca) * 2.0
            rot = SIMD4<Double>(
                0.25 * s,
                (-st * ca + st) / s,
                (st * sa + 0.0) / s,
                (sa - (-ct * sa)) / s
            )
        } else if ct * ca > ca {
            let s = sqrt(1.0 + ct * ca - ct - ca) * 2.0
            rot = SIMD4<Double>(
                (-st * ca + st) / s,
                0.25 * s,
                (-ct * sa + sa) / s,
                (st * sa - 0.0) / s
            )
        } else {
            let s = sqrt(1.0 + ca - ct - ct * ca) * 2.0
            rot = SIMD4<Double>(
                (st * sa + 0.0) / s,
                (-ct * sa + sa) / s,
                0.25 * s,
                (st - (-st * ca)) / s
            )
        }

        self.init(position: pos, rotation: rot)
    }
}

extension SpatialPose where T == Float {
    /// Initialize a SpatialPose using standard Denavit-Hartenberg parameters.
    ///
    /// This is an optimized initializer that directly computes position and quaternion
    /// without creating intermediate transformation matrices.
    ///
    /// [See DH Parameters](https://en.wikipedia.org/wiki/Denavit–Hartenberg_parameters)
    ///
    /// - Parameters:
    ///   - a: Link length (distance along x axis)
    ///   - alpha: Link twist (angle in radians around x axis)
    ///   - d: Link offset (distance along z axis)
    ///   - theta: Joint angle (angle in radians around z axis)
    @inlinable
    public init(denavitHartenberg a: Float, alpha: Float, d: Float, theta: Float) {
        // Precompute trigonometric values
        let ct = cos(theta)
        let st = sin(theta)
        let ca = cos(alpha)
        let sa = sin(alpha)

        // Position from DH transformation
        let pos = SIMD3<Float>(a * ct, a * st, d)

        // Convert rotation matrix to quaternion
        // Matrix: [ ct, -st*ca,  st*sa ]
        //         [ st,  ct*ca, -ct*sa ]
        //         [  0,     sa,     ca ]
        let trace = ct + ct * ca + ca
        let rot: SIMD4<Float>

        if trace > 0 {
            let s = sqrt(trace + 1.0) * 2.0
            rot = SIMD4<Float>(
                (sa - (-ct * sa)) / s,
                (st * sa - 0.0) / s,
                (st - (-st * ca)) / s,
                0.25 * s
            )
        } else if ct > ct * ca && ct > ca {
            let s = sqrt(1.0 + ct - ct * ca - ca) * 2.0
            rot = SIMD4<Float>(
                0.25 * s,
                (-st * ca + st) / s,
                (st * sa + 0.0) / s,
                (sa - (-ct * sa)) / s
            )
        } else if ct * ca > ca {
            let s = sqrt(1.0 + ct * ca - ct - ca) * 2.0
            rot = SIMD4<Float>(
                (-st * ca + st) / s,
                0.25 * s,
                (-ct * sa + sa) / s,
                (st * sa - 0.0) / s
            )
        } else {
            let s = sqrt(1.0 + ca - ct - ct * ca) * 2.0
            rot = SIMD4<Float>(
                (st * sa + 0.0) / s,
                (-ct * sa + sa) / s,
                0.25 * s,
                (st - (-st * ca)) / s
            )
        }

        self.init(position: pos, rotation: rot)
    }

    /// Initialize a SpatialPose using Denavit-Hartenberg parameters with precomputed link twist.
    ///
    /// Use this optimized version when the link twist (alpha) is constant and its sine/cosine
    /// can be precomputed. This saves two trigonometric function calls per invocation.
    ///
    /// [See DH Parameters](https://en.wikipedia.org/wiki/Denavit–Hartenberg_parameters)
    ///
    /// - Parameters:
    ///   - a: Link length (distance along x axis)
    ///   - ca: Cosine of the link twist
    ///   - sa: Sine of the link twist
    ///   - d: Link offset (distance along z axis)
    ///   - theta: Joint angle (angle in radians around z axis)
    @inlinable
    public init(denavitHartenberg a: Float, ca: Float, sa: Float, d: Float, theta: Float) {
        // Compute trigonometric values for theta only
        let ct = cos(theta)
        let st = sin(theta)

        // Position from DH transformation
        let pos = SIMD3<Float>(a * ct, a * st, d)

        // Convert to quaternion
        // Matrix: [ ct, -st*ca,  st*sa ]
        //         [ st,  ct*ca, -ct*sa ]
        //         [  0,     sa,     ca ]
        let trace = ct + ct * ca + ca
        let rot: SIMD4<Float>

        if trace > 0 {
            let s = sqrt(trace + 1.0) * 2.0
            rot = SIMD4<Float>(
                (sa - (-ct * sa)) / s,
                (st * sa - 0.0) / s,
                (st - (-st * ca)) / s,
                0.25 * s
            )
        } else if ct > ct * ca && ct > ca {
            let s = sqrt(1.0 + ct - ct * ca - ca) * 2.0
            rot = SIMD4<Float>(
                0.25 * s,
                (-st * ca + st) / s,
                (st * sa + 0.0) / s,
                (sa - (-ct * sa)) / s
            )
        } else if ct * ca > ca {
            let s = sqrt(1.0 + ct * ca - ct - ca) * 2.0
            rot = SIMD4<Float>(
                (-st * ca + st) / s,
                0.25 * s,
                (-ct * sa + sa) / s,
                (st * sa - 0.0) / s
            )
        } else {
            let s = sqrt(1.0 + ca - ct - ct * ca) * 2.0
            rot = SIMD4<Float>(
                (st * sa + 0.0) / s,
                (-ct * sa + sa) / s,
                0.25 * s,
                (st - (-st * ca)) / s
            )
        }

        self.init(position: pos, rotation: rot)
    }
}
