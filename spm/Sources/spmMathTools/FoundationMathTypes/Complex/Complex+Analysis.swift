import Foundation
import simd

import FoundationTypes

extension Complex {

    /// The complex conjugate (real - imaginary*i).
    @inlinable
    public var conjugate: Complex<T> {
        return Complex(real: real, imaginary: -imaginary)
    }
}


// MARK: - Double-specific properties
extension Complex where T == Double {

    /// The phase (argument) of the complex number in radians.
    @inlinable
    public var phase: T {
        return simd.atan2(imaginary, real)
    }
}

// MARK: - Float-specific properties
extension Complex where T == Float {

    /// The phase (argument) of the complex number in radians.
    @inlinable
    public var phase: T {
        return simd.atan2(imaginary, real)
    }
}
