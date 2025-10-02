import Foundation

// MARK: - Correlation Operations
extension Waveform1D where T: BinaryFloatingPoint {

    /// Compute the auto-correlation of the waveform
    /// - Parameters:
    ///   - maxLag: Maximum lag to compute (default: half the signal length)
    ///   - normalized: If true, normalizes the correlation (default: true)
    /// - Returns: New waveform containing auto-correlation values with lags as time indices
    public func autoCorrelation(maxLag: Int? = nil, normalized: Bool = true) -> Waveform1D<T> {
        guard !values.isEmpty else { return Waveform1D(values: [], dt: dt, t0: nil) }

        let actualMaxLag = maxLag ?? (values.count / 2)
        let clampedMaxLag = min(actualMaxLag, values.count - 1)

        var correlations: [T] = []
        let signalMean = normalized ? mean : T.zero

        // Calculate normalization factor for normalized correlation
        let normalizationFactor: T
        if normalized {
            let variance = values.reduce(T.zero) { sum, value in
                let deviation = value - signalMean
                return sum + deviation * deviation
            }
            normalizationFactor = variance > T.zero ? variance : T(1.0)
        } else {
            normalizationFactor = T(1.0)
        }

        // Compute correlation for each lag
        for lag in 0...clampedMaxLag {
            var correlation = T.zero
            let validSamples = values.count - lag

            for i in 0..<validSamples {
                let val1 = values[i] - signalMean
                let val2 = values[i + lag] - signalMean
                correlation += val1 * val2
            }

            if normalized {
                correlation = correlation / normalizationFactor
            }

            correlations.append(correlation)
        }

        return Waveform1D(values: correlations, dt: dt, t0: nil)
    }

    /// Compute the cross-correlation between this waveform and another
    /// - Parameters:
    ///   - other: The other waveform to correlate with
    ///   - mode: Correlation mode (full, valid, or same)
    ///   - normalized: If true, normalizes the correlation (default: true)
    /// - Returns: New waveform containing cross-correlation values
    public func crossCorrelation(
        with other: Waveform1D<T>,
        mode: WaveformCorrelationMode = .full,
        normalized: Bool = true
    ) -> Waveform1D<T>? {

        guard !values.isEmpty && !other.values.isEmpty else { return nil }

        let x = values
        let y = other.values
        let xMean = normalized ? mean : T.zero
        let yMean = normalized ? other.mean : T.zero

        // Calculate normalization factors
        let xNorm: T
        let yNorm: T

        if normalized {
            let xVariance = x.reduce(T.zero) { sum, value in
                let deviation = value - xMean
                return sum + deviation * deviation
            }
            let yVariance = y.reduce(T.zero) { sum, value in
                let deviation = value - yMean
                return sum + deviation * deviation
            }
            xNorm = xVariance > T.zero ? xVariance.squareRoot() : T(1.0)
            yNorm = yVariance > T.zero ? yVariance.squareRoot() : T(1.0)
        } else {
            xNorm = T(1.0)
            yNorm = T(1.0)
        }

        var correlations: [T] = []

        switch mode {
        case .full:
            // Full correlation: output size = len(x) + len(y) - 1
            let outputSize = x.count + y.count - 1
            let startLag = -(y.count - 1)

            for lag in 0..<outputSize {
                let actualLag = startLag + lag
                var correlation = T.zero

                for i in 0..<x.count {
                    let j = i + actualLag
                    if j >= 0 && j < y.count {
                        let xVal = (x[i] - xMean) / xNorm
                        let yVal = (y[j] - yMean) / yNorm
                        correlation += xVal * yVal
                    }
                }

                correlations.append(correlation)
            }

        case .valid:
            // Valid correlation: output only where signals fully overlap
            guard y.count <= x.count else { return nil }
            let outputSize = x.count - y.count + 1

            for lag in 0..<outputSize {
                var correlation = T.zero

                for i in 0..<y.count {
                    let xVal = (x[lag + i] - xMean) / xNorm
                    let yVal = (y[i] - yMean) / yNorm
                    correlation += xVal * yVal
                }

                correlations.append(correlation)
            }

        case .same:
            // Same correlation: output same size as first input
            let outputSize = x.count
            let startLag = -(y.count - 1) / 2

            for lag in 0..<outputSize {
                let actualLag = startLag + lag
                var correlation = T.zero

                for i in 0..<x.count {
                    let j = i + actualLag
                    if j >= 0 && j < y.count {
                        let xVal = (x[i] - xMean) / xNorm
                        let yVal = (y[j] - yMean) / yNorm
                        correlation += xVal * yVal
                    }
                }

                correlations.append(correlation)
            }
        }

        // Use the sampling rate of the first signal
        return Waveform1D(values: correlations, dt: dt, t0: nil)
    }

    /// Find the lag corresponding to maximum cross-correlation
    /// - Parameters:
    ///   - other: The other waveform to correlate with
    ///   - searchRange: Range of lags to search (default: full range)
    /// - Returns: Tuple containing (lag in samples, lag in time, correlation value)
    public func findMaxCorrelation(
        with other: Waveform1D<T>,
        searchRange: Range<Int>? = nil
    ) -> (lagSamples: Int, lagTime: TimeInterval, correlation: T)? {

        guard let crossCorr = crossCorrelation(with: other, mode: .full, normalized: true) else {
            return nil
        }

        let searchStart = searchRange?.lowerBound ?? 0
        let searchEnd = searchRange?.upperBound ?? crossCorr.values.count

        guard searchStart >= 0 && searchEnd <= crossCorr.values.count && searchStart < searchEnd else {
            return nil
        }

        var maxCorrelation = crossCorr.values[searchStart]
        var maxIndex = searchStart

        for i in searchStart..<searchEnd {
            if crossCorr.values[i] > maxCorrelation {
                maxCorrelation = crossCorr.values[i]
                maxIndex = i
            }
        }

        // Convert to actual lag (accounting for full correlation indexing)
        let actualLag = maxIndex - (other.values.count - 1)
        let lagTime = TimeInterval(actualLag) * dt

        return (actualLag, lagTime, maxCorrelation)
    }
}
