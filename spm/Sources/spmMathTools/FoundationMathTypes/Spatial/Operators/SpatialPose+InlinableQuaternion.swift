import Foundation
import FoundationTypes
import simd

// MARK: - InlinableQuaternion Operators for Double SpatialPose
// Delegates to Quaternion<Double>'s canonical implementations to avoid duplication.

extension SpatialPose where T == Double {

    @inlinable
    @inline(__always)
    static func _qmul(_ q1: SIMD4<Double>, _ q2: SIMD4<Double>) -> SIMD4<Double> {
        Quaternion<Double>._qmul(q1, q2)
    }

    @inlinable
    @inline(__always)
    static func _qrot(_ q: SIMD4<Double>, _ v: SIMD3<Double>) -> SIMD3<Double> {
        Quaternion<Double>._qrot(q, v)
    }
}

// MARK: - InlinableQuaternion Operators for Float SpatialPose
// Delegates to Quaternion<Float>'s canonical implementations to avoid duplication.

extension SpatialPose where T == Float {

    @inlinable
    @inline(__always)
    static func _qmul(_ q1: SIMD4<Float>, _ q2: SIMD4<Float>) -> SIMD4<Float> {
        Quaternion<Float>._qmul(q1, q2)
    }

    @inlinable
    @inline(__always)
    static func _qrot(_ q: SIMD4<Float>, _ v: SIMD3<Float>) -> SIMD3<Float> {
        Quaternion<Float>._qrot(q, v)
    }
}
