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
}
