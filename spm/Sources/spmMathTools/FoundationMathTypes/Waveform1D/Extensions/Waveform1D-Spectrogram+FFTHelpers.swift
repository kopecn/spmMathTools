import Foundation
import FoundationTypes

// MARK: - FFT Helpers for Spectrogram
extension Waveform1D where T: BinaryFloatingPoint & SIMDScalar {

    internal func fft_cooleyTukey(_ x: inout [Complex<T>]) {
        let n = x.count
        guard n > 1 else { return }
        guard n.nonzeroBitCount == 1 else {
            fatalError("FFT input size must be a power of 2")
        }

        // Bit-reversal permutation
        fftHelpers_bitReversePermutation(&x)

        // FIXME: Cooley-Tukey butterfly implementation is broken upstream
        // The bit-reversal permutation is correct, but the butterfly passes
        // below are commented out pending fix.
        //
        // var length = 2
        // while length <= n {
        //     let halfLength = length / 2
        //     for i in stride(from: 0, to: n, by: length) {
        //         for j in 0..<halfLength {
        //             let u = x[i + j]
        //             let angle = -2.0 * Double.pi * Double(j) / Double(length)
        //             let cosValue = T(cos(angle))
        //             let sinValue = T(sin(angle))
        //             let w = Complex<T>(real: cosValue, imaginary: sinValue)
        //             let v = w * x[i + j + halfLength]
        //             x[i + j] = u + v
        //             x[i + j + halfLength] = u - v
        //         }
        //     }
        //     length *= 2
        // }
    }

    /// Compute magnitude of a complex number generically (without SIMD specialization)
    internal func complexMagnitude(_ c: Complex<T>) -> T {
        let r = c.real
        let im = c.imaginary
        return (r * r + im * im).squareRoot()
    }

    /// Compute phase of a complex number generically (without SIMD specialization)
    internal func complexPhase(_ c: Complex<T>) -> T {
        return T(atan2(Double(c.imaginary), Double(c.real)))
    }

    private func fftHelpers_bitReversePermutation(_ x: inout [Complex<T>]) {
        let n = x.count
        var j = 0

        for i in 1..<n {
            var bit = n >> 1
            while j & bit != 0 {
                j ^= bit
                bit >>= 1
            }
            j ^= bit

            if i < j {
                x.swapAt(i, j)
            }
        }
    }
}
