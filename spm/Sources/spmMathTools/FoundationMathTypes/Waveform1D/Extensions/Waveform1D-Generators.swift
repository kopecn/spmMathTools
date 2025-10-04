import Foundation
import FoundationTypes

// MARK: - Waveform Generators
extension Waveform1D where T: BinaryFloatingPoint {

    // MARK: - Sinusoidal Waveforms

    /// Generate a sine wave
    /// - Parameters:
    ///   - frequency: Frequency in Hz
    ///   - amplitude: Peak amplitude (default: 1.0)
    ///   - phase: Phase offset in radians (default: 0.0)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Sine wave waveform
    public static func sine(
        frequency: Double,
        amplitude: T = T(1.0),
        phase: Double = 0.0,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)
        let omega = 2.0 * Double.pi * frequency

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            return amplitude * T(sin(omega * t + phase))
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    /// Generate a cosine wave
    /// - Parameters:
    ///   - frequency: Frequency in Hz
    ///   - amplitude: Peak amplitude (default: 1.0)
    ///   - phase: Phase offset in radians (default: 0.0)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Cosine wave waveform
    public static func cosine(
        frequency: Double,
        amplitude: T = T(1.0),
        phase: Double = 0.0,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)
        let omega = 2.0 * Double.pi * frequency

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            return amplitude * T(cos(omega * t + phase))
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    /// Generate a square wave
    /// - Parameters:
    ///   - frequency: Frequency in Hz
    ///   - amplitude: Peak amplitude (default: 1.0)
    ///   - dutyCycle: Duty cycle (0.0 to 1.0, default: 0.5)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Square wave waveform
    public static func square(
        frequency: Double,
        amplitude: T = T(1.0),
        dutyCycle: Double = 0.5,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)
        let period = 1.0 / frequency
        let clampedDutyCycle = max(0.0, min(1.0, dutyCycle))

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            let phaseInPeriod = (t.truncatingRemainder(dividingBy: period)) / period
            return phaseInPeriod < clampedDutyCycle ? amplitude : -amplitude
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    /// Generate a triangle wave
    /// - Parameters:
    ///   - frequency: Frequency in Hz
    ///   - amplitude: Peak amplitude (default: 1.0)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Triangle wave waveform
    public static func triangle(
        frequency: Double,
        amplitude: T = T(1.0),
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)
        let period = 1.0 / frequency

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            let phaseInPeriod = (t.truncatingRemainder(dividingBy: period)) / period
            let triangleValue =
                phaseInPeriod < 0.5
                ? 4.0 * phaseInPeriod - 1.0
                : 3.0 - 4.0 * phaseInPeriod
            return amplitude * T(triangleValue)
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    /// Generate a sawtooth wave
    /// - Parameters:
    ///   - frequency: Frequency in Hz
    ///   - amplitude: Peak amplitude (default: 1.0)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Sawtooth wave waveform
    public static func sawtooth(
        frequency: Double,
        amplitude: T = T(1.0),
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)
        let period = 1.0 / frequency

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            let phaseInPeriod = (t.truncatingRemainder(dividingBy: period)) / period
            let sawtoothValue = 2.0 * phaseInPeriod - 1.0
            return amplitude * T(sawtoothValue)
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    // MARK: - Chirp and Frequency-Modulated Signals

    /// Generate a linear chirp (frequency sweep)
    /// - Parameters:
    ///   - startFrequency: Starting frequency in Hz
    ///   - endFrequency: Ending frequency in Hz
    ///   - amplitude: Peak amplitude (default: 1.0)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Linear chirp waveform
    public static func chirp(
        startFrequency: Double,
        endFrequency: Double,
        amplitude: T = T(1.0),
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)
        let sweepRate = (endFrequency - startFrequency) / duration

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            let phase = 2.0 * Double.pi * (startFrequency * t + 0.5 * sweepRate * t * t)
            return amplitude * T(sin(phase))
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    // MARK: - Exponential Functions

    /// Generate an exponential decay
    /// - Parameters:
    ///   - amplitude: Initial amplitude (default: 1.0)
    ///   - timeConstant: Time constant (tau) in seconds
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Exponential decay waveform
    public static func exponentialDecay(
        amplitude: T = T(1.0),
        timeConstant: Double,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            return amplitude * T(exp(-t / timeConstant))
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    /// Generate an exponential growth
    /// - Parameters:
    ///   - amplitude: Initial amplitude (default: 1.0)
    ///   - timeConstant: Time constant (tau) in seconds
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Exponential growth waveform
    public static func exponentialGrowth(
        amplitude: T = T(1.0),
        timeConstant: Double,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            return amplitude * T(exp(t / timeConstant))
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    // MARK: - Polynomial Functions

    /// Generate a polynomial waveform
    /// - Parameters:
    ///   - coefficients: Polynomial coefficients [a₀, a₁, a₂, ...] for a₀ + a₁t + a₂t² + ...
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Polynomial waveform
    public static func polynomial(
        coefficients: [Double],
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        guard !coefficients.isEmpty else {
            return Waveform1D(values: [], dt: 1.0 / samplingRate, t0: t0)
        }

        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            var result = 0.0
            for (power, coeff) in coefficients.enumerated() {
                result += coeff * pow(t, Double(power))
            }
            return T(result)
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    /// Generate a linear ramp
    /// - Parameters:
    ///   - startValue: Starting value
    ///   - endValue: Ending value
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Linear ramp waveform
    public static func linearRamp(
        startValue: T,
        endValue: T,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)

        let values = (0..<sampleCount).map { i in
            let progress = T(Double(i) / Double(sampleCount - 1))
            return startValue + (endValue - startValue) * progress
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    // MARK: - Logarithmic Functions

    /// Generate a natural logarithm waveform
    /// - Parameters:
    ///   - amplitude: Amplitude scaling factor (default: 1.0)
    ///   - offset: Time offset to avoid log(0) (default: 0.01)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Natural logarithm waveform
    public static func logarithm(
        amplitude: T = T(1.0),
        offset: Double = 0.01,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt + offset
            return amplitude * T(log(t))
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    /// Generate a base-10 logarithm waveform
    /// - Parameters:
    ///   - amplitude: Amplitude scaling factor (default: 1.0)
    ///   - offset: Time offset to avoid log(0) (default: 0.01)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Base-10 logarithm waveform
    public static func log10(
        amplitude: T = T(1.0),
        offset: Double = 0.01,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt + offset
            return amplitude * T(_math.log10(t))
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    // MARK: - Square Root Functions

    /// Generate a square root waveform
    /// - Parameters:
    ///   - amplitude: Amplitude scaling factor (default: 1.0)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Square root waveform
    public static func squareRoot(
        amplitude: T = T(1.0),
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            return amplitude * T(sqrt(t))
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    // MARK: - Step Functions

    /// Generate a Heaviside step function
    /// - Parameters:
    ///   - amplitude: Step amplitude (default: 1.0)
    ///   - stepTime: Time at which step occurs (default: duration/2)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Heaviside step waveform
    public static func heaviside(
        amplitude: T = T(1.0),
        stepTime: Double? = nil,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)
        let actualStepTime = stepTime ?? (duration / 2.0)

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            return t >= actualStepTime ? amplitude : T.zero
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    /// Generate a ReLU (Rectified Linear Unit) function
    /// - Parameters:
    ///   - slope: Slope of the linear portion (default: 1.0)
    ///   - threshold: Threshold below which output is zero (default: 0.0)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: ReLU waveform
    public static func relu(
        slope: Double = 1.0,
        threshold: Double = 0.0,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            let x = t - threshold
            return x > 0 ? T(slope * x) : T.zero
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    /// Generate a sigmoid function
    /// - Parameters:
    ///   - amplitude: Maximum amplitude (default: 1.0)
    ///   - steepness: Steepness parameter (default: 1.0)
    ///   - center: Center point of sigmoid (default: duration/2)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Sigmoid waveform
    public static func sigmoid(
        amplitude: T = T(1.0),
        steepness: Double = 1.0,
        center: Double? = nil,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)
        let actualCenter = center ?? (duration / 2.0)

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            let x = steepness * (t - actualCenter)
            let sigmoidValue = 1.0 / (1.0 + exp(-x))
            return amplitude * T(sigmoidValue)
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    // MARK: - Noise and Random Signals

    /// Generate white noise
    /// - Parameters:
    ///   - amplitude: RMS amplitude (default: 1.0)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - seed: Random seed for reproducibility
    ///   - t0: Optional start time
    /// - Returns: White noise waveform
    public static func whiteNoise(
        amplitude: T = T(1.0),
        duration: TimeInterval,
        samplingRate: Double,
        seed: UInt64? = nil,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)

        let values = (0..<sampleCount).map { _ in
            let u1 = Double.random(in: 0...1)
            let u2 = Double.random(in: 0...1)
            // Box-Muller transform for Gaussian noise
            let gaussian = sqrt(-2.0 * log(u1)) * cos(2.0 * Double.pi * u2)
            return amplitude * T(gaussian)
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    /// Generate a constant (DC) signal
    /// - Parameters:
    ///   - value: Constant value
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Constant waveform
    public static func constant(
        value: T,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)
        let values = Array(repeating: value, count: sampleCount)

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    // MARK: - Impulse and Delta Functions

    /// Generate a unit impulse (delta function approximation)
    /// - Parameters:
    ///   - amplitude: Impulse amplitude (default: 1.0)
    ///   - impulseTime: Time at which impulse occurs (default: 0.0)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Unit impulse waveform
    public static func impulse(
        amplitude: T = T(1.0),
        impulseTime: Double = 0.0,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)
        let impulseIndex = Int(impulseTime * samplingRate)

        var values = Array(repeating: T.zero, count: sampleCount)
        if impulseIndex >= 0 && impulseIndex < sampleCount {
            values[impulseIndex] = amplitude
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    // MARK: - Composite Waveforms

    /// Generate a damped sinusoid (exponentially decaying sine wave)
    /// - Parameters:
    ///   - frequency: Frequency in Hz
    ///   - amplitude: Initial amplitude (default: 1.0)
    ///   - dampingConstant: Damping time constant in seconds
    ///   - phase: Phase offset in radians (default: 0.0)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Damped sinusoid waveform
    public static func dampedSinusoid(
        frequency: Double,
        amplitude: T = T(1.0),
        dampingConstant: Double,
        phase: Double = 0.0,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)
        let omega = 2.0 * Double.pi * frequency

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            let envelope = exp(-t / dampingConstant)
            return amplitude * T(envelope * sin(omega * t + phase))
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }
}

// MARK: - Integer Type Generators
extension Waveform1D where T: BinaryInteger {

    /// Generate a square wave with integer values
    /// - Parameters:
    ///   - frequency: Frequency in Hz
    ///   - highValue: High level value
    ///   - lowValue: Low level value
    ///   - dutyCycle: Duty cycle (0.0 to 1.0, default: 0.5)
    ///   - duration: Duration in seconds
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Integer square wave waveform
    public static func digitalSquare(
        frequency: Double,
        highValue: T,
        lowValue: T,
        dutyCycle: Double = 0.5,
        duration: TimeInterval,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate
        let sampleCount = Int(duration * samplingRate)
        let period = 1.0 / frequency
        let clampedDutyCycle = max(0.0, min(1.0, dutyCycle))

        let values = (0..<sampleCount).map { i in
            let t = Double(i) * dt
            let phaseInPeriod = (t.truncatingRemainder(dividingBy: period)) / period
            return phaseInPeriod < clampedDutyCycle ? highValue : lowValue
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }

    /// Generate a counter/ramp with integer values
    /// - Parameters:
    ///   - startValue: Starting count value
    ///   - increment: Increment per sample
    ///   - sampleCount: Number of samples
    ///   - samplingRate: Sampling rate in Hz
    ///   - t0: Optional start time
    /// - Returns: Integer counter waveform
    public static func counter(
        startValue: T,
        increment: T,
        sampleCount: Int,
        samplingRate: Double,
        t0: Date? = nil
    ) -> Waveform1D<T> {
        let dt = 1.0 / samplingRate

        var values: [T] = []
        var currentValue = startValue

        for _ in 0..<sampleCount {
            values.append(currentValue)
            currentValue += increment
        }

        return Waveform1D(values: values, dt: dt, t0: t0)
    }
}
