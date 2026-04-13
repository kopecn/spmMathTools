import Foundation
import simd

import FoundationTypes

// MARK: - Double-specific properties
extension Complex where T == Double {

    /// The complex conjugate (real - imaginary*i).
    @inlinable
    public var conjugate: Complex<T> {
        Complex(vector: storage * SIMD2<Double>(1, -1))
    }

    /// The phase (argument) of the complex number in radians.
    @inlinable
    public var phase: T {
        simd.atan2(imaginary, real)
    }
}

// MARK: - Float-specific properties
extension Complex where T == Float {

    /// The complex conjugate (real - imaginary*i).
    @inlinable
    public var conjugate: Complex<T> {
        Complex(vector: storage * SIMD2<Float>(1, -1))
    }

    /// The phase (argument) of the complex number in radians.
    @inlinable
    public var phase: T {
        simd.atan2(imaginary, real)
    }
}
