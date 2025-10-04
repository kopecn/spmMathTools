import Foundation
import FoundationTypes

// MARK: - Decimation and Interpolation
extension Waveform1D where T: BinaryFloatingPoint {

    /// Decimate the waveform by reducing the sampling rate
    /// - Parameters:
    ///   - factor: Decimation factor (must be positive integer)
    ///   - antiAliasFilter: Apply anti-aliasing filter before decimation
    ///   - filterOrder: Order of the anti-aliasing filter
    /// - Returns: Decimated waveform
    public func decimated(
        by factor: Int,
        antiAliasFilter: Bool = true,
        filterOrder: Int = 4
    ) -> Waveform1D<T> {
        guard factor > 0 && !values.isEmpty else { return self }
        guard factor > 1 else { return self }  // No decimation needed

        var processedWaveform = self

        // Apply anti-aliasing filter if requested
        if antiAliasFilter {
            let cutoffFrequency = samplingFrequency / (2.0 * Double(factor))
            processedWaveform = lowPassFilter(cutoffFrequency: cutoffFrequency, order: filterOrder)
        }

        // Perform decimation by keeping every factor-th sample
        let decimatedValues = stride(from: 0, to: processedWaveform.values.count, by: factor)
            .map { processedWaveform.values[$0] }

        let newDt = dt * Double(factor)

        return Waveform1D(values: decimatedValues, dt: newDt, t0: t0)
    }

    /// Interpolate the waveform to increase the sampling rate
    /// - Parameters:
    ///   - factor: Interpolation factor (must be positive integer)
    ///   - method: Interpolation method to use
    ///   - antiAliasFilter: Apply anti-aliasing filter after interpolation
    ///   - filterOrder: Order of the anti-aliasing filter
    /// - Returns: Interpolated waveform
    public func interpolated(
        by factor: Int,
        method: WaveformInterpolationMethod = .linear,
        antiAliasFilter: Bool = true,
        filterOrder: Int = 4
    ) -> Waveform1D<T> {
        guard factor > 0 && !values.isEmpty else { return self }
        guard factor > 1 else { return self }  // No interpolation needed

        let interpolatedValues = performInterpolation(factor: factor, method: method)
        let newDt = dt / Double(factor)

        var interpolatedWaveform = Waveform1D(values: interpolatedValues, dt: newDt, t0: t0)

        // Apply anti-aliasing filter if requested
        if antiAliasFilter {
            let cutoffFrequency = samplingFrequency / 2.0  // Original Nyquist frequency
            interpolatedWaveform = interpolatedWaveform.lowPassFilter(
                cutoffFrequency: cutoffFrequency,
                order: filterOrder
            )
        }

        return interpolatedWaveform
    }

    /// Resample the waveform to a new sampling rate using rational resampling
    /// - Parameters:
    ///   - newSamplingRate: Target sampling rate in Hz
    ///   - method: Interpolation method for fractional resampling
    ///   - antiAliasFilter: Apply anti-aliasing filters
    /// - Returns: Resampled waveform
    public func resampled(
        to newSamplingRate: Double,
        method: WaveformInterpolationMethod = .linear,
        antiAliasFilter: Bool = true
    ) -> Waveform1D<T> {
        guard newSamplingRate > 0 && !values.isEmpty else { return self }

        let currentRate = samplingFrequency
        guard abs(newSamplingRate - currentRate) > 1e-10 else { return self }  // Already at target rate

        if newSamplingRate < currentRate {
            // Downsampling - use decimation
            let ratio = currentRate / newSamplingRate
            if ratio == floor(ratio) {
                // Integer decimation
                return decimated(by: Int(ratio), antiAliasFilter: antiAliasFilter)
            } else {
                // Fractional decimation
                return fractionalResample(targetRate: newSamplingRate, method: method, antiAliasFilter: antiAliasFilter)
            }
        } else {
            // Upsampling - use interpolation
            let ratio = newSamplingRate / currentRate
            if ratio == floor(ratio) {
                // Integer interpolation
                return interpolated(by: Int(ratio), method: method, antiAliasFilter: antiAliasFilter)
            } else {
                // Fractional interpolation
                return fractionalResample(targetRate: newSamplingRate, method: method, antiAliasFilter: antiAliasFilter)
            }
        }
    }

    /// Resample to match another waveform's sampling characteristics
    /// - Parameters:
    ///   - other: Reference waveform to match
    ///   - method: Interpolation method
    /// - Returns: Resampled waveform with matching sampling rate
    public func resampledToMatch<U>(
        _ other: Waveform1D<U>,
        method: WaveformInterpolationMethod = .linear
    ) -> Waveform1D<T> {
        return resampled(to: other.samplingFrequency, method: method)
    }

    // MARK: - Private Implementation Methods

    private func performInterpolation(factor: Int, method: WaveformInterpolationMethod) -> [T] {
        guard values.count > 1 else { return values }

        var interpolatedValues: [T] = []

        for i in 0..<(values.count - 1) {
            interpolatedValues.append(values[i])

            // Insert interpolated values between samples
            for j in 1..<factor {
                let t = T(j) / T(factor)
                let interpolatedValue = interpolateValue(from: values[i], to: values[i + 1], t: t, method: method)
                interpolatedValues.append(interpolatedValue)
            }
        }

        // Add the last sample
        interpolatedValues.append(values.last!)

        return interpolatedValues
    }

    private func interpolateValue(from start: T, to end: T, t: T, method: WaveformInterpolationMethod) -> T {
        switch method {
        case .linear:
            return start + t * (end - start)
        case .cubic:
            // Simplified cubic interpolation (Hermite interpolation with zero derivatives)
            let t2 = t * t
            let t3 = t2 * t
            let h00 = T(2) * t3 - T(3) * t2 + T(1)
            _ = t3 - T(2) * t2 + t
            let h01 = -T(2) * t3 + T(3) * t2
            _ = t3 - t2

            // Assuming zero derivatives at endpoints for simplicity
            return h00 * start + h01 * end
        case .nearestNeighbor:
            return t < T(0.5) ? start : end
        }
    }

    private func fractionalResample(
        targetRate: Double,
        method: WaveformInterpolationMethod,
        antiAliasFilter: Bool
    ) -> Waveform1D<T> {
        let currentRate = samplingFrequency
        let ratio = targetRate / currentRate
        let newLength = Int(Double(values.count) * ratio)

        var resampledValues: [T] = []
        resampledValues.reserveCapacity(newLength)

        for i in 0..<newLength {
            let sourceIndex = Double(i) / ratio
            let floorIndex = Int(sourceIndex)
            let fractionalPart = sourceIndex - Double(floorIndex)

            if floorIndex >= values.count - 1 {
                resampledValues.append(values.last!)
            } else {
                let interpolatedValue = interpolateValue(
                    from: values[floorIndex],
                    to: values[floorIndex + 1],
                    t: T(fractionalPart),
                    method: method
                )
                resampledValues.append(interpolatedValue)
            }
        }

        let newDt = 1.0 / targetRate
        var resampledWaveform = Waveform1D(values: resampledValues, dt: newDt, t0: t0)

        // Apply anti-aliasing filter if requested
        if antiAliasFilter {
            let cutoffFrequency = min(currentRate, targetRate) / 2.0
            resampledWaveform = resampledWaveform.lowPassFilter(cutoffFrequency: cutoffFrequency)
        }

        return resampledWaveform
    }
}

// MARK: - Polyphase Filtering for Efficient Resampling
extension Waveform1D where T: BinaryFloatingPoint {

    /// Efficient polyphase decimation with built-in filtering
    /// - Parameters:
    ///   - factor: Decimation factor
    ///   - filterTaps: Number of filter taps (default: 64)
    /// - Returns: Decimated waveform with polyphase filtering
    public func polyphaseDecimated(by factor: Int, filterTaps: Int = 64) -> Waveform1D<T> {
        guard factor > 1 && !values.isEmpty else { return self }

        // Design a simple low-pass FIR filter
        let cutoffFrequency = 0.5 / Double(factor)  // Normalized cutoff frequency
        let filterCoeffs = designLowPassFIR(taps: filterTaps, cutoff: cutoffFrequency)

        // Apply polyphase decimation
        let decimatedValues = polyphaseFilter(coefficients: filterCoeffs, decimationFactor: factor)
        let newDt = dt * Double(factor)

        return Waveform1D(values: decimatedValues, dt: newDt, t0: t0)
    }

    /// Efficient polyphase interpolation with built-in filtering
    /// - Parameters:
    ///   - factor: Interpolation factor
    ///   - filterTaps: Number of filter taps (default: 64)
    /// - Returns: Interpolated waveform with polyphase filtering
    public func polyphaseInterpolated(by factor: Int, filterTaps: Int = 64) -> Waveform1D<T> {
        guard factor > 1 && !values.isEmpty else { return self }

        // Design a simple low-pass FIR filter
        let cutoffFrequency = 0.5 / Double(factor)  // Normalized cutoff frequency
        let filterCoeffs = designLowPassFIR(taps: filterTaps, cutoff: cutoffFrequency)

        // Zero-pad for interpolation
        var interpolatedValues: [T] = []
        for value in values {
            interpolatedValues.append(value * T(factor))  // Scale by interpolation factor
            for _ in 1..<factor {
                interpolatedValues.append(T.zero)
            }
        }

        // Apply low-pass filter
        let filteredValues = firFilter(interpolatedValues, coefficients: filterCoeffs)
        let newDt = dt / Double(factor)

        return Waveform1D(values: filteredValues, dt: newDt, t0: t0)
    }

    // MARK: - Private Filter Implementation

    private func designLowPassFIR(taps: Int, cutoff: Double) -> [T] {
        let halfTaps = taps / 2
        var coefficients: [T] = []

        for i in 0..<taps {
            let n = i - halfTaps
            if n == 0 {
                coefficients.append(T(2.0 * cutoff))
            } else {
                let sinc = sin(2.0 * .pi * cutoff * Double(n)) / (.pi * Double(n))
                // Apply Hamming window
                let window = 0.54 - 0.46 * cos(2.0 * .pi * Double(i) / Double(taps - 1))
                coefficients.append(T(sinc * window))
            }
        }

        return coefficients
    }

    private func polyphaseFilter(coefficients: [T], decimationFactor: Int) -> [T] {
        let outputLength = values.count / decimationFactor
        var output: [T] = []
        output.reserveCapacity(outputLength)

        for outputIndex in 0..<outputLength {
            let inputIndex = outputIndex * decimationFactor
            var sum = T.zero

            for coeffIndex in 0..<coefficients.count {
                let sampleIndex = inputIndex - coeffIndex
                if sampleIndex >= 0 && sampleIndex < values.count {
                    sum += values[sampleIndex] * coefficients[coeffIndex]
                }
            }

            output.append(sum)
        }

        return output
    }

    private func firFilter(_ input: [T], coefficients: [T]) -> [T] {
        var output: [T] = []
        output.reserveCapacity(input.count)

        for i in 0..<input.count {
            var sum = T.zero

            for j in 0..<coefficients.count {
                let inputIndex = i - j
                if inputIndex >= 0 && inputIndex < input.count {
                    sum += input[inputIndex] * coefficients[j]
                }
            }

            output.append(sum)
        }

        return output
    }
}

// MARK: - Integer Support
extension Waveform1D where T: BinaryInteger {

    /// Decimate integer waveform using nearest neighbor
    /// - Parameter factor: Decimation factor
    /// - Returns: Decimated waveform
    public func decimated(by factor: Int) -> Waveform1D<T> {
        guard factor > 1 && !values.isEmpty else { return self }

        let decimatedValues = stride(from: 0, to: values.count, by: factor)
            .map { values[$0] }

        let newDt = dt * Double(factor)

        return Waveform1D(values: decimatedValues, dt: newDt, t0: t0)
    }

    /// Interpolate integer waveform using nearest neighbor or linear interpolation
    /// - Parameters:
    ///   - factor: Interpolation factor
    ///   - method: Interpolation method (for integers, only nearestNeighbor makes sense)
    /// - Returns: Interpolated waveform
    public func interpolated(
        by factor: Int,
        method: WaveformInterpolationMethod = .nearestNeighbor
    ) -> Waveform1D<T> {
        guard factor > 1 && !values.isEmpty else { return self }

        var interpolatedValues: [T] = []

        for i in 0..<(values.count - 1) {
            interpolatedValues.append(values[i])

            // Insert interpolated values
            for _ in 1..<factor {
                switch method {
                case .nearestNeighbor:
                    interpolatedValues.append(values[i])
                case .linear, .cubic:
                    // For integers, approximate linear interpolation
                    let midpoint = (Int(values[i]) + Int(values[i + 1])) / 2
                    interpolatedValues.append(T(midpoint))
                }
            }
        }

        interpolatedValues.append(values.last!)
        let newDt = dt / Double(factor)

        return Waveform1D(values: interpolatedValues, dt: newDt, t0: t0)
    }
}
