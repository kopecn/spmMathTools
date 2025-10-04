import Foundation
import FoundationTypes

// MARK: - FFT Helpers for Spectrogram (if not available from main FFT extension)
extension Waveform1D where T: BinaryFloatingPoint {

    // This file provides FFT helpers if they're not accessible from the main FFT extension
    // These should match the implementations in Waveform1D-FFT.swift

    internal func fft_cooleyTukey(_ x: inout [Complex]) {
        let n = x.count
        guard n > 1 else { return }
        guard n.nonzeroBitCount == 1 else {
            fatalError("FFT input size must be a power of 2")
        }

        // Bit-reversal permutation
        bitReversePermutation(&x)

        // Cooley-Tukey FFT
        var length = 2
        while length <= n {
            for i in stride(from: 0, to: n, by: length) {
                for j in 0..<(length / 2) {
                    let u = x[i + j]
                    let angle = -2.0 * Double.pi * Double(j) / Double(length)
                    let w = Complex(real: cos(angle), imaginary: sin(angle))
                    let v = w * x[i + j + length / 2]

                    x[i + j] = u + v
                    x[i + j + length / 2] = u - v
                }
            }
            length *= 2
        }
    }

    private func bitReversePermutation(_ x: inout [Complex]) {
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
