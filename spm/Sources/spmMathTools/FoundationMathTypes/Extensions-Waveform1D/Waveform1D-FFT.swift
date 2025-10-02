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



// MARK: - Power Spectral Density Operations
extension Waveform1D where T: BinaryFloatingPoint {
    
    /// Compute the Power Spectral Density using Welch's method
    /// - Parameters:
    ///   - WaveformWindowType: Window function to apply (default: Hanning)
    ///   - windowSize: Size of each segment (default: signal length / 8)
    ///   - overlap: Overlap between segments as fraction (0.0 to 1.0, default: 0.5)
    ///   - scaling: PSD scaling type (default: density)
    /// - Returns: Tuple containing frequency array and PSD values
    public func powerSpectralDensity(WaveformWindowType: WaveformWindowType = .hanning,
                                   windowSize: Int? = nil,
                                   overlap: Double = 0.5,
                                   scaling: WaveformPSDScaling = .density) -> (frequencies: [T], psd: [T])? {
        
        guard !values.isEmpty else { return nil }
        
        let actualWindowSize = windowSize ?? max(256, values.count / 8)
        let clampedOverlap = max(0.0, min(0.99, overlap))
        let hopSize = Int(Double(actualWindowSize) * (1.0 - clampedOverlap))
        
        guard actualWindowSize <= values.count && hopSize > 0 else {
            // Fall back to single segment
            return singleSegmentPSD(WaveformWindowType: WaveformWindowType, scaling: scaling)
        }
        
        // Generate window
        let window = generateWindow(type: WaveformWindowType, length: actualWindowSize)
        let windowPower = window.reduce(0.0) { $0 + $1 * $1 }
        
        var psdAccumulator: [Double] = []
        var segmentCount = 0
        
        // Process overlapping segments
        var startIndex = 0
        while startIndex + actualWindowSize <= values.count {
            let segment = Array(values[startIndex..<(startIndex + actualWindowSize)])
            
            // Apply window
            let windowedSegment = zip(segment, window).map { Double($0.0) * $0.1 }
            
            // Convert to complex and perform FFT
            var complex = windowedSegment.map { Complex(real: $0, imaginary: 0.0) }
            let fftSize = nextPowerOfTwo(complex.count)
            while complex.count < fftSize {
                complex.append(Complex(real: 0.0, imaginary: 0.0))
            }
            
            fft_cooleyTukey(&complex)
            
            // Calculate power spectrum for this segment
            let nyquistBin = fftSize / 2
            let segmentPSD = (0..<nyquistBin).map { i -> Double in
                let magnitude = complex[i].magnitude
                return magnitude * magnitude
            }
            
            // Accumulate PSD
            if psdAccumulator.isEmpty {
                psdAccumulator = segmentPSD
            } else {
                for i in 0..<min(psdAccumulator.count, segmentPSD.count) {
                    psdAccumulator[i] += segmentPSD[i]
                }
            }
            
            segmentCount += 1
            startIndex += hopSize
        }
        
        guard segmentCount > 0 && !psdAccumulator.isEmpty else { return nil }
        
        // Average across segments
        let averagedPSD = psdAccumulator.map { $0 / Double(segmentCount) }
        
        // Apply scaling
        let samplingRate = 1.0 / dt
        let scaledPSD = applyWaveformPSDScaling(averagedPSD, 
                                       samplingRate: samplingRate, 
                                       windowPower: windowPower, 
                                       windowSize: actualWindowSize, 
                                       scaling: scaling)
        
        // Generate frequency array
        let frequencies = (0..<scaledPSD.count).map { i in
            T(Double(i) * samplingRate / Double(actualWindowSize))
        }
        
        return (frequencies, scaledPSD.map { T($0) })
    }
    
    /// Compute PSD for a single segment (no windowing/averaging)
    /// - Parameters:
    ///   - WaveformWindowType: Window function to apply
    ///   - scaling: PSD scaling type
    /// - Returns: Tuple containing frequency array and PSD values
    public func singleSegmentPSD(WaveformWindowType: WaveformWindowType = .hanning,
                                scaling: WaveformPSDScaling = .density) -> (frequencies: [T], psd: [T])? {
        
        guard !values.isEmpty else { return nil }
        
        // Generate window
        let window = generateWindow(type: WaveformWindowType, length: values.count)
        let windowPower = window.reduce(0.0) { $0 + $1 * $1 }
        
        // Apply window
        let windowedValues = zip(values, window).map { Double($0.0) * $0.1 }
        
        // Convert to complex and perform FFT
        var complex = windowedValues.map { Complex(real: $0, imaginary: 0.0) }
        let fftSize = nextPowerOfTwo(complex.count)
        while complex.count < fftSize {
            complex.append(Complex(real: 0.0, imaginary: 0.0))
        }
        
        fft_cooleyTukey(&complex)
        
        // Calculate power spectrum
        let nyquistBin = fftSize / 2
        let powerSpectrum = (0..<nyquistBin).map { i -> Double in
            let magnitude = complex[i].magnitude
            return magnitude * magnitude
        }
        
        // Apply scaling
        let samplingRate = 1.0 / dt
        let scaledPSD = applyWaveformPSDScaling(powerSpectrum,
                                       samplingRate: samplingRate,
                                       windowPower: windowPower,
                                       windowSize: values.count,
                                       scaling: scaling)
        
        // Generate frequency array
        let frequencies = (0..<scaledPSD.count).map { i in
            T(Double(i) * samplingRate / Double(fftSize))
        }
        
        return (frequencies, scaledPSD.map { T($0) })
    }
    
    /// Apply appropriate scaling to PSD values
    private func applyWaveformPSDScaling(_ powerSpectrum: [Double],
                                 samplingRate: Double,
                                 windowPower: Double,
                                 windowSize: Int,
                                 scaling: WaveformPSDScaling) -> [Double] {
        
        switch scaling {
        case .density:
            // Power spectral density: units^2 / Hz
            let scale = 1.0 / (samplingRate * windowPower)
            return powerSpectrum.map { $0 * scale }
            
        case .spectrum:
            // Power spectrum: units^2
            let scale = 1.0 / (windowPower * windowPower)
            return powerSpectrum.map { $0 * scale }
        }
    }
    
    /// Generate window coefficients (moved from main file if needed)
    private func generateWindow(type: WaveformWindowType, length: Int) -> [Double] {
        guard length > 0 else { return [] }
        guard length > 1 else { return [1.0] }
        
        let n = Double(length)
        
        switch type {
        case .rectangular:
            return Array(repeating: 1.0, count: length)
            
        case .hanning:
            return (0..<length).map { i in
                0.5 * (1.0 - cos(2.0 * .pi * Double(i) / (n - 1.0)))
            }
            
        case .hamming:
            return (0..<length).map { i in
                0.54 - 0.46 * cos(2.0 * .pi * Double(i) / (n - 1.0))
            }
            
        case .blackman:
            return (0..<length).map { i in
                let factor = 2.0 * .pi * Double(i) / (n - 1.0)
                return 0.42 - 0.5 * cos(factor) + 0.08 * cos(2.0 * factor)
            }
            
        case .blackmanHarris:
            return (0..<length).map { i in
                let factor = 2.0 * .pi * Double(i) / (n - 1.0)
                return 0.35875 - 0.48829 * cos(factor) + 0.14128 * cos(2.0 * factor) - 0.01168 * cos(3.0 * factor)
            }
            
        case .kaiser(beta: let beta):
            let alpha = (n - 1.0) / 2.0
            let i0Beta = modifiedBesselI0(beta)
            
            return (0..<length).map { i in
                let x = (Double(i) - alpha) / alpha
                let arg = beta * sqrt(max(0.0, 1.0 - x * x)) // Clamp to avoid sqrt of negative
                return modifiedBesselI0(arg) / i0Beta
            }
            
        case .tukey(taperRatio: let taperRatio):
            let clampedR = max(0.0, min(1.0, taperRatio))
            
            return (0..<length).map { i in
                let x = Double(i) / (n - 1.0)
                
                if x < clampedR / 2.0 {
                    // Taper up
                    return 0.5 * (1.0 + cos(.pi * (2.0 * x / clampedR - 1.0)))
                } else if x > 1.0 - clampedR / 2.0 {
                    // Taper down
                    return 0.5 * (1.0 + cos(.pi * (2.0 * (x - 1.0 + clampedR / 2.0) / clampedR)))
                } else {
                    // Flat top
                    return 1.0
                }
            }
            
        case .bartlett:
            return (0..<length).map { i in
                let x = Double(i) / (n - 1.0)
                return 1.0 - 2.0 * abs(x - 0.5)
            }
            
        case .welch:
            return (0..<length).map { i in
                let x = (Double(i) - (n - 1.0) / 2.0) / ((n - 1.0) / 2.0)
                return 1.0 - x * x
            }
        }
    }
}

/// Modified Bessel function of the first kind, order 0 (for Kaiser window)
private func modifiedBesselI0(_ x: Double) -> Double {
    let ax = abs(x)
    var ans: Double
    
    if ax < 3.75 {
        let y = x / 3.75
        let y2 = y * y
        ans = 1.0 + y2 * (3.5156229 + y2 * (3.0899424 + y2 * (1.2067492 + y2 * (0.2659732 + y2 * (0.360768e-1 + y2 * 0.45813e-2)))))
    } else {
        let y = 3.75 / ax
        ans = (exp(ax) / sqrt(ax)) * (0.39894228 + y * (0.1328592e-1 + y * (0.225319e-2 + y * (-0.157565e-2 + y * (0.916281e-2 + y * (-0.2057706e-1 + y * (0.2635537e-1 + y * (-0.1647633e-1 + y * 0.392377e-2))))))))
    }
    
    return ans
}
