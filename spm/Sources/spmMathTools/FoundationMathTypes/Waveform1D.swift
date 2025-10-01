import Foundation

/// A one-dimensional waveform data structure representing time-series data.  Based upon [Labview 1d Waveform](https://www.ni.com/docs/en-US/bundle/labview-nxg-data-types/page/waveforms.html)
///
/// `Waveform1D` stores uniformly sampled data points with their temporal characteristics,
/// making it suitable for signal processing, scientific measurements, and time-series analysis.
///
/// Example usage:
/// ```swift
/// let startTime = Date()
/// let samplingInterval = 0.001 // 1ms sampling
/// let samples = [1.0, 2.0, 3.0, 4.0, 5.0]
/// var waveform = Waveform1D(values: samples, dt: samplingInterval, t0: startTime)
/// 
/// // Simple initialization with defaults
/// var simpleWaveform = Waveform1D(values: [1.0, 2.0, 3.0])
/// simpleWaveform.dt = 0.5 // Modify sampling interval
/// ```
public struct Waveform1D {
    /// The sampled data values of the waveform
    public var values: [Double]
    
    /// The time interval between consecutive samples in seconds
    public var dt: TimeInterval
    
    /// The absolute start time of the first sample
    public var t0: Date?
    
    /// Initialize with all parameters
    public init(values: [Double], dt: TimeInterval, t0: Date?) {
        self.values = values
        self.dt = dt
        self.t0 = t0
    }
    
    /// Initialize with values only, using default dt=1.0 and t0=nil
    public init(values: [Double]) {
        self.init(values: values, dt: 1.0, t0: nil)
    }
    
    /// Initialize with values and dt, using default t0=nil
    public init(values: [Double], dt: TimeInterval) {
        self.init(values: values, dt: dt, t0: nil)
    }
}