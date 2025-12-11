import Foundation
import FoundationTypes
import simd

// MARK: - InlinableQuaternion Operators for Double SpatialPose

extension SpatialPose where T == Double {

    // MARK: - Inline Quaternion Operations

    @inlinable
    @inline(__always)
    static func _qmul(_ q1: SIMD4<Double>, _ q2: SIMD4<Double>) -> SIMD4<Double> {
        SIMD4<Double>(
            q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
            q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
            q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w,
            q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z
        )
    }

    @inlinable
    @inline(__always)
    static func _qrot(_ q: SIMD4<Double>, _ v: SIMD3<Double>) -> SIMD3<Double> {
        v + 2.0
            * simd_cross(
                SIMD3<Double>(q.x, q.y, q.z),
                simd_cross(SIMD3<Double>(q.x, q.y, q.z), v) + q.w * v
            )
    }
}

// MARK: - InlinableQuaternion Operators for Float SpatialPose

extension SpatialPose where T == Float {

    // MARK: - Inline Quaternion Operations

    @inlinable
    @inline(__always)
    static func _qmul(_ q1: SIMD4<Float>, _ q2: SIMD4<Float>) -> SIMD4<Float> {
        SIMD4<Float>(
            q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y,
            q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x,
            q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w,
            q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z
        )
    }

    @inlinable
    @inline(__always)
    static func _qrot(_ q: SIMD4<Float>, _ v: SIMD3<Float>) -> SIMD3<Float> {
        v + 2.0
            * simd_cross(
                SIMD3<Float>(q.x, q.y, q.z),
                simd_cross(SIMD3<Float>(q.x, q.y, q.z), v) + q.w * v
            )
    }
}
