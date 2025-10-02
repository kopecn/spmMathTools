import Foundation

// MARK: - Envelope Detection
extension Waveform1D where T: BinaryFloatingPoint & Comparable {
    
    /// Detect the amplitude envelope using the Hilbert transform approximation
    /// - Parameters:
    ///   - method: Envelope detection method (default: hilbert)
    ///   - smoothing: Optional smoothing window size for post-processing
    /// - Returns: New waveform representing the amplitude envelope
    public func amplitudeEnvelope(method: WaveformEnvelopeMethod = .hilbert, 
                                 smoothing: Int? = nil) -> Waveform1D<T> {
        guard !values.isEmpty else { return Waveform1D(values: [], dt: dt, t0: t0) }
        
        let envelope: [T]
        
        switch method {
        case .hilbert:
            envelope = hilbertEnvelope()
        case .peakDetection(let windowSize):
            envelope = peakDetectionEnvelope(windowSize: windowSize)
        case .rms(let windowSize):
            envelope = rmsEnvelope(windowSize: windowSize)
        case .absoluteValue:
            envelope = absoluteValueEnvelope()
        }
        
        // Apply smoothing if requested
        let finalEnvelope = smoothing.map { smoothEnvelope(envelope, windowSize: $0) } ?? envelope
        
        return Waveform1D(values: finalEnvelope, dt: dt, t0: t0)
    }
    
    /// Detect upper and lower envelopes using peak and valley detection
    /// - Parameters:
    ///   - windowSize: Window size for local extrema detection
    ///   - interpolationMethod: Method for interpolating between extrema
    /// - Returns: Tuple containing upper and lower envelope waveforms
    public func upperLowerEnvelopes(windowSize: Int = 5, 
                                   interpolationMethod: WaveformInterpolationMethod = .linear) 
    -> (upper: Waveform1D<T>, lower: Waveform1D<T>) {
        
        guard values.count > windowSize else {
            let constantUpper = Waveform1D(values: Array(repeating: maximum ?? T.zero, count: values.count), dt: dt, t0: t0)
            let constantLower = Waveform1D(values: Array(repeating: minimum ?? T.zero, count: values.count), dt: dt, t0: t0)
            return (constantUpper, constantLower)
        }
        
        // Find local maxima and minima
        let peaks = findLocalExtrema(type: .maxima, windowSize: windowSize)
        let valleys = findLocalExtrema(type: .minima, windowSize: windowSize)
        
        // Ensure we have boundary points
        var upperPoints = peaks
        var lowerPoints = valleys
        
        // Add boundary points if they're not already extrema
        if upperPoints.first?.index != 0 {
            upperPoints.insert((index: 0, value: values[0]), at: 0)
        }
        if upperPoints.last?.index != values.count - 1 {
            upperPoints.append((index: values.count - 1, value: values[values.count - 1]))
        }
        
        if lowerPoints.first?.index != 0 {
            lowerPoints.insert((index: 0, value: values[0]), at: 0)
        }
        if lowerPoints.last?.index != values.count - 1 {
            lowerPoints.append((index: values.count - 1, value: values[values.count - 1]))
        }
        
        // Interpolate envelopes
        let upperEnvelope = interpolateEnvelope(points: upperPoints, method: interpolationMethod)
        let lowerEnvelope = interpolateEnvelope(points: lowerPoints, method: interpolationMethod)
        
        return (
            Waveform1D(values: upperEnvelope, dt: dt, t0: t0),
            Waveform1D(values: lowerEnvelope, dt: dt, t0: t0)
        )
    }
    
    /// Calculate the instantaneous amplitude using multiple methods
    /// - Parameter method: Method for calculating instantaneous amplitude
    /// - Returns: New waveform representing instantaneous amplitude
    public func instantaneousAmplitude(method: WaveformInstantaneousMethod = .magnitude) -> Waveform1D<T> {
        switch method {
        case .magnitude:
            return amplitudeEnvelope(method: .hilbert)
        case .power:
            let envelope = amplitudeEnvelope(method: .hilbert)
            let powerValues = envelope.values.map { $0 * $0 }
            return Waveform1D(values: powerValues, dt: dt, t0: t0)
        case .logMagnitude:
            let envelope = amplitudeEnvelope(method: .hilbert)
            let logValues = envelope.values.map { value in
                value > T.zero ? T(log(Double(value))) : T(-Double.infinity)
            }
            return Waveform1D(values: logValues, dt: dt, t0: t0)
        }
    }
    
    // MARK: - Private Implementation Methods
    
    private func hilbertEnvelope() -> [T] {
        // Simplified Hilbert transform approximation using quadrature filter
        // For a more accurate implementation, use FFT-based Hilbert transform
        
        guard values.count > 2 else { return values.map { abs($0) } }
        
        var envelope: [T] = []
        
        // Use a simple quadrature approximation
        for i in 0..<values.count {
            let real = values[i]
            let imaginary = quadratureComponent(at: i)
            let magnitude = (real * real + imaginary * imaginary).squareRoot()
            envelope.append(magnitude)
        }
        
        return envelope
    }
    
    private func quadratureComponent(at index: Int) -> T {
        // Simple quadrature filter approximation
        let windowSize = min(5, values.count / 4)
        guard windowSize > 1 else { return T.zero }
        
        var sum = T.zero
        var count = 0
        
        for i in max(0, index - windowSize)...min(values.count - 1, index + windowSize) {
            if i != index {
                let weight = T(1.0 / Double(abs(i - index)))
                sum += values[i] * weight
                count += 1
            }
        }
        
        return count > 0 ? sum / T(count) : T.zero
    }
    
    private func peakDetectionEnvelope(windowSize: Int) -> [T] {
        var envelope: [T] = []
        let halfWindow = windowSize / 2
        
        for i in 0..<values.count {
            let start = max(0, i - halfWindow)
            let end = min(values.count, i + halfWindow + 1)
            
            let windowMax = values[start..<end].max() ?? values[i]
            envelope.append(abs(windowMax))
        }
        
        return envelope
    }
    
    private func rmsEnvelope(windowSize: Int) -> [T] {
        var envelope: [T] = []
        let halfWindow = windowSize / 2
        
        for i in 0..<values.count {
            let start = max(0, i - halfWindow)
            let end = min(values.count, i + halfWindow + 1)
            
            let windowValues = Array(values[start..<end])
            let rms = (windowValues.reduce(T.zero) { $0 + $1 * $1 } / T(windowValues.count)).squareRoot()
            envelope.append(rms)
        }
        
        return envelope
    }
    
    private func absoluteValueEnvelope() -> [T] {
        return values.map { abs($0) }
    }
    
    private func smoothEnvelope(_ envelope: [T], windowSize: Int) -> [T] {
        guard windowSize > 1 else { return envelope }
        
        var smoothed: [T] = []
        let halfWindow = windowSize / 2
        
        for i in 0..<envelope.count {
            let start = max(0, i - halfWindow)
            let end = min(envelope.count, i + halfWindow + 1)
            
            let windowSum = envelope[start..<end].reduce(T.zero, +)
            let average = windowSum / T(end - start)
            smoothed.append(average)
        }
        
        return smoothed
    }
    
    private func findLocalExtrema(type: WaveformExtremaType, windowSize: Int) -> [(index: Int, value: T)] {
        var extrema: [(index: Int, value: T)] = []
        let halfWindow = windowSize / 2
        
        for i in halfWindow..<(values.count - halfWindow) {
            let start = i - halfWindow
            let end = i + halfWindow + 1
            let window = Array(values[start..<end])
            let centerValue = values[i]
            
            let isExtremum: Bool
            switch type {
            case .maxima:
                isExtremum = window.allSatisfy { $0 <= centerValue }
            case .minima:
                isExtremum = window.allSatisfy { $0 >= centerValue }
            }
            
            if isExtremum {
                extrema.append((index: i, value: centerValue))
            }
        }
        
        return extrema
    }
    
    private func interpolateEnvelope(points: [(index: Int, value: T)], 
                                   method: WaveformInterpolationMethod) -> [T] {
        guard points.count >= 2 else {
            return Array(repeating: points.first?.value ?? T.zero, count: values.count)
        }
        
        var envelope: [T] = Array(repeating: T.zero, count: values.count)
        
        for i in 0..<values.count {
            switch method {
            case .linear:
                envelope[i] = linearInterpolate(at: i, points: points)
            case .cubic:
                envelope[i] = cubicInterpolate(at: i, points: points)
            case .nearestNeighbor:
                envelope[i] = nearestNeighborInterpolate(at: i, points: points)
            }
        }
        
        return envelope
    }
    
    private func linearInterpolate(at index: Int, points: [(index: Int, value: T)]) -> T {
        // Find surrounding points
        var leftPoint = points[0]
        var rightPoint = points[points.count - 1]
        
        for i in 0..<(points.count - 1) {
            if points[i].index <= index && points[i + 1].index >= index {
                leftPoint = points[i]
                rightPoint = points[i + 1]
                break
            }
        }
        
        guard leftPoint.index != rightPoint.index else { return leftPoint.value }
        
        let t = T(index - leftPoint.index) / T(rightPoint.index - leftPoint.index)
        return leftPoint.value + t * (rightPoint.value - leftPoint.value)
    }
    
    private func cubicInterpolate(at index: Int, points: [(index: Int, value: T)]) -> T {
        // Simplified cubic interpolation - fall back to linear for now
        return linearInterpolate(at: index, points: points)
    }
    
    private func nearestNeighborInterpolate(at index: Int, points: [(index: Int, value: T)]) -> T {
        let distances = points.map { abs($0.index - index) }
        let minIndex = distances.firstIndex(of: distances.min()!) ?? 0
        return points[minIndex].value
    }
}

    
