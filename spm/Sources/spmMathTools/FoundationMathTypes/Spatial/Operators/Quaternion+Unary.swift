import Foundation
import FoundationTypes
import simd

// MARK: - Unary Operators for Double Quaternion

extension Quaternion where T == Double {

    /// Negate a quaternion (flip all components)
    public static prefix func - (quaternion: Quaternion<Double>) -> Quaternion<Double> {
        return Quaternion(vector: -quaternion.vector)
    }

    /// Unary plus (returns copy)
    public static prefix func + (quaternion: Quaternion<Double>) -> Quaternion<Double> {
        return quaternion
    }

}

// MARK: - Unary Operators for Float Quaternion

extension Quaternion where T == Float {

    /// Negate a quaternion (flip all components)
    public static prefix func - (quaternion: Quaternion<Float>) -> Quaternion<Float> {
        return Quaternion(vector: -quaternion.vector)
    }

    /// Unary plus (returns copy)
    public static prefix func + (quaternion: Quaternion<Float>) -> Quaternion<Float> {
        return quaternion
    }
}
