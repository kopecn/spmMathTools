import Foundation

// MARK: - FFT Mathematical Operations
extension Waveform1D where T: BinaryFloatingPoint {

    /// Compute the Fast Fourier Transform using Cooley-Tukey algorithm
    /// - Returns: Tuple containing frequency array, magnitudes, and phases
    public func fft() -> (frequencies: [T], magnitudes: [T], phases: [T])? {
        guard values.count > 1 else { return nil }

        // Convert to complex numbers
        var complex = values.map { Complex(real: Double($0), imaginary: 0.0) }

        // Pad to next power of 2 for efficiency
        let fftSize = nextPowerOfTwo(complex.count)
        while complex.count < fftSize {
            complex.append(Complex(real: 0.0, imaginary: 0.0))
        }

        // Perform FFT
        fft_cooleyTukey(&complex)

        // Calculate frequencies (only positive frequencies)
        let samplingRate = 1.0 / dt
        let nyquistBin = fftSize / 2
        let frequencies = (0..<nyquistBin).map { i in
            T(Double(i) * samplingRate / Double(fftSize))
        }

        // Calculate magnitudes and phases
        let magnitudes = (0..<nyquistBin).map { i in
            T(complex[i].magnitude)
        }

        let phases = (0..<nyquistBin).map { i in
            T(complex[i].phase)
        }

        return (frequencies, magnitudes, phases)
    }

    // MARK: - Private FFT Implementation

    private func nextPowerOfTwo(_ n: Int) -> Int {
        guard n > 1 else { return 1 }
        return 1 << Int(ceil(log2(Double(n))))
    }

    private func fft_cooleyTukey(_ x: inout [Complex]) {
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
            let step = n / length
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

// MARK: - Complex Number Helper
private struct Complex {
    let real: Double
    let imaginary: Double

    var magnitude: Double {
        return sqrt(real * real + imaginary * imaginary)
    }

    var phase: Double {
        return atan2(imaginary, real)
    }

    static func + (lhs: Complex, rhs: Complex) -> Complex {
        return Complex(real: lhs.real + rhs.real, imaginary: lhs.imaginary + rhs.imaginary)
    }

    static func - (lhs: Complex, rhs: Complex) -> Complex {
        return Complex(real: lhs.real - rhs.real, imaginary: lhs.imaginary - rhs.imaginary)
    }

    static func * (lhs: Complex, rhs: Complex) -> Complex {
        return Complex(
            real: lhs.real * rhs.real - lhs.imaginary * rhs.imaginary,
            imaginary: lhs.real * rhs.imaginary + lhs.imaginary * rhs.real
        )
    }
}

// MARK: - Statistical and Utility Operations
extension Waveform1D where T: BinaryFloatingPoint {

    /// Calculate the root mean square (RMS) value
    public var rms: T {
        guard !values.isEmpty else { return T.zero }
        let sumSquares = values.reduce(T.zero) { $0 + $1 * $1 }
        return (sumSquares / T(values.count)).squareRoot()
    }

    /// Calculate the mean value
    public var mean: T {
        guard !values.isEmpty else { return T.zero }
        return values.reduce(T.zero, +) / T(values.count)
    }

    /// Calculate the peak-to-peak amplitude
    public var peakToPeak: T? {
        guard !values.isEmpty else { return nil }
        let min = values.min()!
        let max = values.max()!
        return max - min
    }

    /// Get the total duration of the waveform
    public var duration: TimeInterval {
        guard values.count > 1 else { return 0 }
        return TimeInterval(values.count - 1) * dt
    }

    /// Get the sampling frequency (Hz)
    public var samplingFrequency: Double {
        return 1.0 / dt
    }

    /// Get the Nyquist frequency (Hz)
    public var nyquistFrequency: Double {
        return samplingFrequency / 2.0
    }
}
