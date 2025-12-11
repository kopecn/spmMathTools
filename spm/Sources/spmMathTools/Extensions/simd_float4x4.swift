import Foundation
import simd

extension simd_float4x4 {
    /// Initialize a transformation matrix using standard Denavit-Hartenberg parameters.  [see DH](https://en.wikipedia.org/wiki/Denavit–Hartenberg_parameters)
    /// - Parameters:
    ///   - a: Link length (distance along x axis)
    ///   - alpha: Link twist (angle in radians around x axis)
    ///   - d: Link offset (distance along z axis)
    ///   - theta: Joint angle (angle in radians around z axis)
    public init(denavitHartenberg a: Float, alpha: Float, d: Float, theta: Float) {
        let ct: Float = cos(theta)
        let st: Float = sin(theta)
        let ca: Float = cos(alpha)
        let sa: Float = sin(alpha)

        self.init(
            columns: (
                SIMD4<Float>(ct, st, 0, 0),
                SIMD4<Float>(-st * ca, ct * ca, sa, 0),
                SIMD4<Float>(st * sa, -ct * sa, ca, 0),
                SIMD4<Float>(a * ct, a * st, d, 1)
            )
        )
    }

    /// Initialize a transformation matrix using standard Denavit-Hartenberg parameters.  [see DH](https://en.wikipedia.org/wiki/Denavit–Hartenberg_parameters)
    /// - Notes: Use this to save compute since cos and sin for the twist link does not need to be recomputed
    /// - Parameters:
    ///   - a: Link length (distance along x axis)
    ///   - sa: sine of the link twist (angle in radians around x axis)
    ///   - ca: cosine of the link twist (angle in radians around x axis)
    ///   - d: Link offset (distance along z axis)
    ///   - theta: Joint angle (angle in radians around z axis)
    public init(denavitHartenberg a: Float, ca: Float, sa: Float, d: Float, theta: Float) {
        let ct: Float = cos(theta)
        let st: Float = sin(theta)

        self.init(
            columns: (
                SIMD4<Float>(ct, st, 0, 0),
                SIMD4<Float>(-st * ca, ct * ca, sa, 0),
                SIMD4<Float>(st * sa, -ct * sa, ca, 0),
                SIMD4<Float>(a * ct, a * st, d, 1)
            )
        )
    }
}
