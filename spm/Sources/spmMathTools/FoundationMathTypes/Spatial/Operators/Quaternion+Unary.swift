import Foundation
import FoundationTypes
import simd

// MARK: - Unary Operators for Double Quaternion

extension Quaternion where T == Double {

    /// Negate a quaternion (flip all components)
    @inlinable
    public static prefix func - (quaternion: Quaternion<Double>) -> Quaternion<Double> {
        Quaternion(vector: -quaternion.vector)
    }

    /// Unary plus (returns copy)
    @inlinable
    public static prefix func + (quaternion: Quaternion<Double>) -> Quaternion<Double> {
        quaternion
    }

}

// MARK: - Unary Operators for Float Quaternion

extension Quaternion where T == Float {

    /// Negate a quaternion (flip all components)
    @inlinable
    public static prefix func - (quaternion: Quaternion<Float>) -> Quaternion<Float> {
        Quaternion(vector: -quaternion.vector)
    }

    /// Unary plus (returns copy)
    @inlinable
    public static prefix func + (quaternion: Quaternion<Float>) -> Quaternion<Float> {
        quaternion
    }
}
