import Foundation
import simd

import FoundationTypes

// MARK: - Arithmetic Operations (SIMD-accelerated)
extension Complex where T == Float {

    /// Adds two complex numbers using SIMD acceleration.
    @inlinable
    public static func + (lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
        Complex(vector: lhs.storage + rhs.storage)
    }

    /// Subtracts two complex numbers using SIMD acceleration.
    @inlinable
    public static func - (lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
        Complex(vector: lhs.storage - rhs.storage)
    }

    /// Multiplies two complex numbers: (a+bi)(c+di) = (ac-bd) + (ad+bc)i
    /// Uses SIMD broadcast + shuffle: a*[c,d] + b*[-d,c]
    @inlinable
    public static func * (lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
        let a = SIMD2<T>(repeating: lhs.storage.x)
        let b = SIMD2<T>(repeating: lhs.storage.y)
        let rNeg = SIMD2<T>(-rhs.storage.y, rhs.storage.x)
        return Complex(vector: a * rhs.storage + b * rNeg)
    }

    /// Multiplies a complex number by a scalar.
    @inlinable
    public static func * (lhs: Complex<T>, rhs: T) -> Complex<T> {
        Complex(vector: lhs.storage * rhs)
    }

    /// Multiplies a scalar by a complex number.
    @inlinable
    public static func * (lhs: T, rhs: Complex<T>) -> Complex<T> {
        rhs * lhs
    }

    /// Divides a complex number by a scalar.
    @inlinable
    public static func / (lhs: Complex<T>, rhs: T) -> Complex<T> {
        Complex(vector: lhs.storage / rhs)
    }

    /// Negates a complex number.
    @inlinable
    public static prefix func - (operand: Complex<T>) -> Complex<T> {
        Complex(vector: -operand.storage)
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
        lhs.storage *= rhs
    }

    /// Divides this complex number by a scalar.
    @inlinable
    public static func /= (lhs: inout Complex<T>, rhs: T) {
        lhs.storage /= rhs
    }

    /// Divides two complex numbers: (a+bi)/(c+di) = [(ac+bd) + (bc-ad)i] / (c²+d²)
    /// Computed as lhs * conjugate(rhs) / |rhs|² using SIMD broadcast + shuffle.
    @inlinable
    public static func / (lhs: Complex<T>, rhs: Complex<T>) -> Complex<T> {
        let denom = rhs.magnitudeSquared
        let rhsConj = SIMD2<T>(rhs.storage.x, -rhs.storage.y)
        let a = SIMD2<T>(repeating: lhs.storage.x)
        let b = SIMD2<T>(repeating: lhs.storage.y)
        let rNeg = SIMD2<T>(-rhsConj.y, rhsConj.x)   // [d, c]
        return Complex(vector: (a * rhsConj + b * rNeg) / denom)
    }

    /// Divides this complex number by another complex number.
    @inlinable
    public static func /= (lhs: inout Complex<T>, rhs: Complex<T>) {
        lhs = lhs / rhs
    }
}
