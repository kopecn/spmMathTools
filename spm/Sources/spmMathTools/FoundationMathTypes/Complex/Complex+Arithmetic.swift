import Foundation
import simd

import FoundationTypes

// MARK: - Arithmetic Operations (SIMD-accelerated)
extension Complex where T == Float {

    /// Adds two complex numbers using SIMD acceleration.
    @inlinable
    public static func + (lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
        return Complex(vector: lhs.storage + rhs.storage)
    }

    /// Subtracts two complex numbers using SIMD acceleration.
    @inlinable
    public static func - (lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
        return Complex(vector: lhs.storage - rhs.storage)
    }

    /// Multiplies two complex numbers: (a+bi)(c+di) = (ac-bd) + (ad+bc)i
    @inlinable
    public static func * (lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
        // Using SIMD operations where possible
        let ac = lhs.real * rhs.real
        let bd = lhs.imaginary * rhs.imaginary
        let ad = lhs.real * rhs.imaginary
        let bc = lhs.imaginary * rhs.real
        return Complex(real: ac - bd, imaginary: ad + bc)
    }

    /// Multiplies a complex number by a scalar.
    @inlinable
    public static func * (lhs: Complex<T>, rhs: T) -> Complex<T> {
        let result = lhs.storage * SIMD2(repeating: rhs)
        return Complex(real: result.x, imaginary: result.y)
    }

    /// Multiplies a scalar by a complex number.
    @inlinable
    public static func * (lhs: T, rhs: Complex<T>) -> Complex<T> {
        return rhs * lhs
    }

    /// Divides a complex number by a scalar.
    @inlinable
    public static func / (lhs: Complex<T>, rhs: T) -> Complex<T> {
        let result = lhs.storage / SIMD2(repeating: rhs)
        return Complex(real: result.x, imaginary: result.y)
    }

    /// Negates a complex number.
    @inlinable
    public static prefix func - (operand: Complex<T>) -> Complex<T> {
        let result = -operand.storage
        return Complex(real: result.x, imaginary: result.y)
    }

    // MARK: - Compound Assignment Operators

    /// Adds a complex number to this complex number.
    @inlinable
    public static func += (lhs: inout Complex<T>, rhs: Complex<T>) {
        lhs.storage += rhs.storage
    }

    /// Subtracts a complex number from this complex number.
    @inlinable
    public static func -= (lhs: inout Complex<T>, rhs: Complex<T>) {
        lhs.storage -= rhs.storage
    }

    /// Multiplies this complex number by another complex number.
    @inlinable
    public static func *= (lhs: inout Complex<T>, rhs: Complex<T>) {
        lhs = lhs * rhs
    }

    /// Multiplies this complex number by a scalar.
    @inlinable
    public static func *= (lhs: inout Complex<T>, rhs: T) {
        lhs.storage *= SIMD2(repeating: rhs)
    }

    /// Divides this complex number by a scalar.
    @inlinable
    public static func /= (lhs: inout Complex<T>, rhs: T) {
        lhs.storage /= SIMD2(repeating: rhs)
    }

    /// Divides two complex numbers: (a+bi)/(c+di) = [(ac+bd) + (bc-ad)i] / (c²+d²)
    @inlinable
    public static func / (lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
        let denominator = rhs.magnitudeSquared
        let ac = lhs.real * rhs.real
        let bd = lhs.imaginary * rhs.imaginary
        let bc = lhs.imaginary * rhs.real
        let ad = lhs.real * rhs.imaginary
        return Complex(real: (ac + bd) / denominator, imaginary: (bc - ad) / denominator)
    }

    /// Divides this complex number by another complex number.
    @inlinable
    public static func /= (lhs: inout Complex<T>, rhs: Complex<T>) {
        lhs = lhs / rhs
    }
}
