import Foundation

// MARK: - Window Functions
extension Waveform1D where T: BinaryFloatingPoint {
    
    /// Apply a window function to the waveform
    /// - Parameter WaveformWindowType: The type of window to apply
    /// - Returns: New waveform with the window function applied
    public func windowed(with WaveformWindowType: WaveformWindowType) -> Waveform1D<T> {
        guard !values.isEmpty else { return self }
        
        let windowCoefficients = Self.generateWindow(type: WaveformWindowType, length: values.count)
        let windowedValues = zip(values, windowCoefficients).map { $0 * T($1) }
        
        return Waveform1D(values: windowedValues, dt: dt, t0: t0)
    }
    
    /// Generate window coefficients for a given length and type
    /// - Parameters:
    ///   - type: The window function type
    ///   - length: The number of samples in the window
    /// - Returns: Array of window coefficients
    public static func generateWindow(type: WaveformWindowType, length: Int) -> [Double] {
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
            
        case .kaiser(let beta):
            let alpha = (n - 1.0) / 2.0
            let i0Beta = modifiedBesselI0(beta)
            
            return (0..<length).map { i in
                let x = (Double(i) - alpha) / alpha
                let arg = beta * sqrt(1.0 - x * x)
                return modifiedBesselI0(arg) / i0Beta
            }
            
        case .tukey(let r):
            let clampedR = max(0.0, min(1.0, r))
            
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
    
    /// Calculate the coherent gain of a window function
    /// - Parameter WaveformWindowType: The window function type
    /// - Returns: The coherent gain factor
    public func windowCoherentGain(for WaveformWindowType: WaveformWindowType) -> Double {
        let window = Self.generateWindow(type: WaveformWindowType, length: values.count)
        return window.reduce(0.0, +) / Double(window.count)
    }
    
    /// Calculate the processing gain of a window function
    /// - Parameter WaveformWindowType: The window function type
    /// - Returns: The processing gain factor
    public func windowProcessingGain(for WaveformWindowType: WaveformWindowType) -> Double {
        let window = Self.generateWindow(type: WaveformWindowType, length: values.count)
        let sumSquares = window.reduce(0.0) { $0 + $1 * $1 }
        return sqrt(sumSquares / Double(window.count))
    }
    
    // MARK: - Private Helper Functions
    
    /// Modified Bessel function of the first kind, order 0
    private static func modifiedBesselI0(_ x: Double) -> Double {
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
}
