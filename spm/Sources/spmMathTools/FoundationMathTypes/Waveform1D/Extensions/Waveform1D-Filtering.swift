import Foundation

// MARK: - Digital Filtering
extension Waveform1D where T: BinaryFloatingPoint {
    
    /// Apply a digital filter to the waveform
    /// - Parameters:
    ///   - filterType: Type of filter to apply
    ///   - order: Filter order (higher order = steeper rolloff)
    /// - Returns: Filtered waveform
    public func filtered(with filterType: WaveformFilterType, order: Int = 4) -> Waveform1D<T> {
        guard !values.isEmpty && order > 0 else { return self }
        
        switch filterType {
        case .lowPass(let cutoffFrequency):
            return applyButterworthFilter(type: .lowPass, cutoffFrequency: cutoffFrequency, order: order)
        case .highPass(let cutoffFrequency):
            return applyButterworthFilter(type: .highPass, cutoffFrequency: cutoffFrequency, order: order)
        case .bandPass(let lowFrequency, let highFrequency):
            return applyBandPassFilter(lowFrequency: lowFrequency, highFrequency: highFrequency, order: order)
        case .bandStop(let lowFrequency, let highFrequency):
            return applyBandStopFilter(lowFrequency: lowFrequency, highFrequency: highFrequency, order: order)
        case .movingAverage(let windowSize):
            return applyMovingAverage(windowSize: windowSize)
        case .exponential(let alpha):
            return applyExponentialFilter(alpha: alpha)
        }
    }
    
    /// Apply a simple moving average filter
    /// - Parameter windowSize: Size of the moving average window
    /// - Returns: Filtered waveform
    public func movingAverageFilter(windowSize: Int) -> Waveform1D<T> {
        return filtered(with: .movingAverage(windowSize: windowSize))
    }
    
    /// Apply an exponential smoothing filter
    /// - Parameter alpha: Smoothing factor (0 < alpha <= 1, higher = less smoothing)
    /// - Returns: Filtered waveform
    public func exponentialFilter(alpha: Double) -> Waveform1D<T> {
        return filtered(with: .exponential(alpha: alpha))
    }
    
    /// Apply a Butterworth low-pass filter
    /// - Parameters:
    ///   - cutoffFrequency: Cutoff frequency in Hz
    ///   - order: Filter order
    /// - Returns: Filtered waveform
    public func lowPassFilter(cutoffFrequency: Double, order: Int = 4) -> Waveform1D<T> {
        return filtered(with: .lowPass(cutoffFrequency: cutoffFrequency), order: order)
    }
    
    /// Apply a Butterworth high-pass filter
    /// - Parameters:
    ///   - cutoffFrequency: Cutoff frequency in Hz
    ///   - order: Filter order
    /// - Returns: Filtered waveform
    public func highPassFilter(cutoffFrequency: Double, order: Int = 4) -> Waveform1D<T> {
        return filtered(with: .highPass(cutoffFrequency: cutoffFrequency), order: order)
    }
    
    /// Apply a band-pass filter
    /// - Parameters:
    ///   - lowFrequency: Lower cutoff frequency in Hz
    ///   - highFrequency: Upper cutoff frequency in Hz
    ///   - order: Filter order
    /// - Returns: Filtered waveform
    public func bandPassFilter(lowFrequency: Double, highFrequency: Double, order: Int = 4) -> Waveform1D<T> {
        return filtered(with: .bandPass(lowFrequency: lowFrequency, highFrequency: highFrequency), order: order)
    }
    
    // MARK: - Private Implementation Methods
    
    private func applyButterworthFilter(type: WaveformFilterButterworthType, cutoffFrequency: Double, order: Int) -> Waveform1D<T> {
        let nyquist = samplingFrequency / 2.0
        let normalizedCutoff = cutoffFrequency / nyquist
        
        // Clamp normalized frequency to valid range
        let clampedCutoff = max(0.001, min(0.999, normalizedCutoff))
        
        // Calculate filter coefficients
        let coefficients = calculateButterworthCoefficients(type: type, normalizedCutoff: clampedCutoff, order: order)
        
        // Apply IIR filter
        return applyIIRFilter(coefficients: coefficients)
    }
    
    private func applyBandPassFilter(lowFrequency: Double, highFrequency: Double, order: Int) -> Waveform1D<T> {
        // Apply high-pass first, then low-pass
        let highPassed = highPassFilter(cutoffFrequency: lowFrequency, order: order)
        return highPassed.lowPassFilter(cutoffFrequency: highFrequency, order: order)
    }
    
    private func applyBandStopFilter(lowFrequency: Double, highFrequency: Double, order: Int) -> Waveform1D<T> {
        // Create band-pass filter and subtract from original
        let bandPass = bandPassFilter(lowFrequency: lowFrequency, highFrequency: highFrequency, order: order)
        
        let filtered = zip(values, bandPass.values).map { original, bandPassValue in
            original - bandPassValue
        }
        
        return Waveform1D(values: filtered, dt: dt, t0: t0)
    }
    
    private func applyMovingAverage(windowSize: Int) -> Waveform1D<T> {
        guard windowSize > 0 && windowSize <= values.count else { return self }
        
        var filtered: [T] = []
        let halfWindow = windowSize / 2
        
        for i in 0..<values.count {
            let start = max(0, i - halfWindow)
            let end = min(values.count, i + halfWindow + 1)
            
            let windowSum = values[start..<end].reduce(T.zero, +)
            let average = windowSum / T(end - start)
            filtered.append(average)
        }
        
        return Waveform1D(values: filtered, dt: dt, t0: t0)
    }
    
    private func applyExponentialFilter(alpha: Double) -> Waveform1D<T> {
        guard !values.isEmpty && alpha > 0 && alpha <= 1 else { return self }
        
        var filtered: [T] = []
        let alphaT = T(alpha)
        let oneMinusAlpha = T(1.0 - alpha)
        
        // Initialize with first value
        var smoothed = values[0]
        filtered.append(smoothed)
        
        // Apply exponential smoothing
        for i in 1..<values.count {
            smoothed = alphaT * values[i] + oneMinusAlpha * smoothed
            filtered.append(smoothed)
        }
        
        return Waveform1D(values: filtered, dt: dt, t0: t0)
    }
    
    private func calculateButterworthCoefficients(type: WaveformFilterButterworthType, normalizedCutoff: Double, order: Int) -> WaveformFilterCoefficients<T> {
        // Simplified Butterworth coefficient calculation
        // For a more complete implementation, you'd use proper pole-zero placement
        
        let wc = tan(.pi * normalizedCutoff / 2.0)
        let wc2 = wc * wc
        
        switch type {
        case .lowPass:
            return calculateLowPassCoefficients(wc: wc, wc2: wc2, order: order)
        case .highPass:
            return calculateHighPassCoefficients(wc: wc, wc2: wc2, order: order)
        }
    }
    
    private func calculateLowPassCoefficients(wc: Double, wc2: Double, order: Int) -> WaveformFilterCoefficients<T> {
        // Simplified 2nd order Butterworth low-pass
        let k1 = sqrt(2.0) * wc
        let k2 = wc2
        let denominator = k2 + k1 + 1.0
        
        let b = [
            T(k2 / denominator),
            T(2.0 * k2 / denominator),
            T(k2 / denominator)
        ]
        
        let a = [
            T(1.0),
            T(2.0 * (k2 - 1.0) / denominator),
            T((k2 - k1 + 1.0) / denominator)
        ]
        
        return WaveformFilterCoefficients(b: b, a: a)
    }
    
    private func calculateHighPassCoefficients(wc: Double, wc2: Double, order: Int) -> WaveformFilterCoefficients<T> {
        // Simplified 2nd order Butterworth high-pass
        let k1 = sqrt(2.0) * wc
        let k2 = wc2
        let denominator = k2 + k1 + 1.0
        
        let b = [
            T(1.0 / denominator),
            T(-2.0 / denominator),
            T(1.0 / denominator)
        ]
        
        let a = [
            T(1.0),
            T(2.0 * (k2 - 1.0) / denominator),
            T((k2 - k1 + 1.0) / denominator)
        ]
        
        return WaveformFilterCoefficients(b: b, a: a)
    }
    
    private func applyIIRFilter(coefficients: WaveformFilterCoefficients<T>) -> Waveform1D<T> {
        guard !values.isEmpty else { return self }
        
        let b = coefficients.b
        let a = coefficients.a
        let orderB = b.count
        let orderA = a.count
        
        var filtered: [T] = []
        var inputHistory: [T] = Array(repeating: T.zero, count: orderB)
        var outputHistory: [T] = Array(repeating: T.zero, count: orderA)
        
        for input in values {
            // Shift input history
            for i in stride(from: orderB - 1, to: 0, by: -1) {
                inputHistory[i] = inputHistory[i - 1]
            }
            inputHistory[0] = input
            
            // Calculate output
            var output = T.zero
            
            // Feedforward (numerator)
            for i in 0..<orderB {
                output += b[i] * inputHistory[i]
            }
            
            // Feedback (denominator, excluding a[0])
            for i in 1..<orderA {
                output -= a[i] * outputHistory[i]
            }
            
            // Shift output history
            for i in stride(from: orderA - 1, to: 0, by: -1) {
                outputHistory[i] = outputHistory[i - 1]
            }
            outputHistory[0] = output
            
            filtered.append(output)
        }
        
        return Waveform1D(values: filtered, dt: dt, t0: t0)
    }
}


// MARK: - Frequency Response Analysis
extension Waveform1D where T: BinaryFloatingPoint {
    
    /// Calculate the frequency response of the signal
    /// - Parameters:
    ///   - frequencies: Array of frequencies to evaluate
    ///   - filterType: Filter type to analyze
    ///   - order: Filter order
    /// - Returns: Tuple containing magnitude and phase response
    public func frequencyResponse(at frequencies: [Double], 
                                 filterType: WaveformFilterType, 
                                 order: Int = 4) -> (magnitude: [T], phase: [T]) {
        
        let nyquist = samplingFrequency / 2.0
        var magnitudes: [T] = []
        var phases: [T] = []
        
        for freq in frequencies {
            let normalizedFreq = freq / nyquist
            let omega = .pi * normalizedFreq
            
            let response = calculateFilterResponse(omega: omega, filterType: filterType, order: order)
            magnitudes.append(T(response.magnitude))
            phases.append(T(response.phase))
        }
        
        return (magnitudes, phases)
    }
    
    private func calculateFilterResponse(omega: Double, filterType: WaveformFilterType, order: Int) -> (magnitude: Double, phase: Double) {
        // Simplified frequency response calculation
        // For demonstration purposes - a complete implementation would use proper transfer functions
        
        switch filterType {
        case .lowPass(let cutoffFreq):
            let wc = 2.0 * .pi * cutoffFreq / samplingFrequency
            let ratio = omega / wc
            let magnitude = 1.0 / sqrt(1.0 + pow(ratio, 2.0 * Double(order)))
            let phase = -Double(order) * atan(ratio)
            return (magnitude, phase)
            
        case .highPass(let cutoffFreq):
            let wc = 2.0 * .pi * cutoffFreq / samplingFrequency
            let ratio = wc / omega
            let magnitude = 1.0 / sqrt(1.0 + pow(ratio, 2.0 * Double(order)))
            let phase = Double(order) * atan(ratio)
            return (magnitude, phase)
            
        default:
            return (1.0, 0.0) // Placeholder for other filter types
        }
    }
}