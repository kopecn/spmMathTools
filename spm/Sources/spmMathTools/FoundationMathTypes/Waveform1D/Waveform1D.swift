import Foundation

// Convenience typealiases for common use cases
public typealias DoubleWaveform1D = Waveform1D<Double>
public typealias FloatWaveform1D = Waveform1D<Float>
public typealias IntWaveform1D = Waveform1D<Int>

/// A one-dimensional waveform data structure representing time-series data.
///
/// `Waveform1D` stores uniformly sampled data points with their temporal characteristics,
/// making it suitable for signal processing, scientific measurements, and time-series analysis.
///
/// Example usage:
/// ```swift
/// let startTime = Date()
/// let samplingInterval = 0.001 // 1ms sampling
/// let samples: [Double] = [1.0, 2.0, 3.0, 4.0, 5.0]
/// var waveform = Waveform1D(values: samples, dt: samplingInterval, t0: startTime)
///
/// // Integer waveform
/// var intWaveform = Waveform1D(values: [1, 2, 3, 4, 5])
///
/// // Float waveform
/// var floatWaveform = Waveform1D<Float>(values: [1.0, 2.0, 3.0])
/// ```
public struct Waveform1D<T: Numeric> {
    /// The sampled data values of the waveform
    public var values: [T]

    /// The time interval between consecutive samples in seconds
    public var dt: TimeInterval

    /// The absolute start time of the first sample
    public var t0: Date?

    /// Initialize with all parameters
    public init(values: [T], dt: TimeInterval, t0: Date?) {
        self.values = values
        self.dt = dt
        self.t0 = t0
    }

    /// Initialize with values only, using default dt=1.0 and t0=nil
    public init(values: [T]) {
        self.init(values: values, dt: 1.0, t0: nil)
    }

    /// Initialize with values and dt, using default t0=nil
    public init(values: [T], dt: TimeInterval) {
        self.init(values: values, dt: dt, t0: nil)
    }

    // MARK: - Computed Properties (Available to all numeric types)

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

    /// Get the end time of the waveform
    public var endTime: Date? {
        guard let t0 = t0 else { return nil }
        return t0.addingTimeInterval(duration)
    }

    /// Get the number of samples
    public var sampleCount: Int {
        return values.count
    }
}

// MARK: - Computed Properties for Comparable Types
extension Waveform1D where T: Comparable {

    /// Calculate the peak-to-peak amplitude
    public var peakToPeak: T? {
        guard !values.isEmpty else { return nil }
        let min = values.min()!
        let max = values.max()!
        return max - min
    }

    /// Get the minimum value
    public var minimum: T? {
        return values.min()
    }

    /// Get the maximum value
    public var maximum: T? {
        return values.max()
    }
}

// MARK: - Computed Properties for Floating Point Types
extension Waveform1D where T: BinaryFloatingPoint {

    /// Calculate the mean value
    public var mean: T {
        guard !values.isEmpty else { return T.zero }
        return values.reduce(T.zero, +) / T(values.count)
    }

    /// Calculate the root mean square (RMS) value
    public var rms: T {
        guard !values.isEmpty else { return T.zero }
        let sumSquares = values.reduce(T.zero) { $0 + $1 * $1 }
        return (sumSquares / T(values.count)).squareRoot()
    }

    /// Calculate the standard deviation
    public var standardDeviation: T {
        guard values.count > 1 else { return T.zero }
        let meanValue = mean
        let variance =
            values.reduce(T.zero) { sum, value in
                let deviation = value - meanValue
                return sum + deviation * deviation
            } / T(values.count - 1)
        return variance.squareRoot()
    }

    /// Calculate the variance
    public var variance: T {
        guard values.count > 1 else { return T.zero }
        let meanValue = mean
        return values.reduce(T.zero) { sum, value in
            let deviation = value - meanValue
            return sum + deviation * deviation
        } / T(values.count - 1)
    }
}

// MARK: - Computed Properties for Integer Types
extension Waveform1D where T: BinaryInteger {

    /// Calculate the mean value (integer division)
    public var mean: T {
        guard !values.isEmpty else { return T.zero }
        return values.reduce(T.zero, +) / T(values.count)
    }

    /// Calculate the sum of all values
    public var sum: T {
        return values.reduce(T.zero, +)
    }
}

// MARK: - Computed Properties for Signed Integer Types
extension Waveform1D where T: SignedInteger {

    /// Calculate the absolute sum of all values
    public var absoluteSum: T {
        return values.reduce(T.zero) { $0 + abs($1) }
    }
}
