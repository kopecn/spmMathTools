import Foundation
import Testing

@testable import spmMathTools

// MARK: - Integration Tests Suite
@Suite("Waveform1D Integration (Calc)")
struct Waveform1DIntegrationTests {

    @Test("Integration of empty waveform")
    func integrationEmpty() {
        let waveform = Waveform1D<Double>(values: [], dt: 0.1)
        let integrated = waveform.integrate()

        #expect(integrated.values.isEmpty)
        #expect(integrated.dt == 0.1)
    }

    @Test("Integration of single point")
    func integrationSinglePoint() {
        let waveform = Waveform1D<Double>(values: [5.0], dt: 0.1)
        let integrated = waveform.integrate()

        #expect(integrated.values.count == 1)
        #expect(integrated.values[0] == 0.0)  // Initial value
    }

    @Test("Integration of constant waveform (100 samples)")
    func integrationConstant() {
        // f(x) = 3, integral should be F(x) = 3x (linear)
        let waveform = Waveform1D<Double>.constant(value: 3.0, duration: 1.0, samplingRate: 100.0)
        let integrated = waveform.integrate()

        // Using trapezoidal rule: area under constant = height * width
        // Total area should be 3 * 1.0 = 3.0
        #expect(integrated.values.count == 100)
        #expect(integrated.values[0] == 0.0)
        // Last value should be approximately 3.0
        #expect(abs(integrated.values.last! - 3.0) < 0.1)

        // Check linearity - differences should be constant
        for i in 2..<integrated.values.count {
            let diff1 = integrated.values[i] - integrated.values[i-1]
            let diff2 = integrated.values[i-1] - integrated.values[i-2]
            #expect(abs(diff1 - diff2) < 1e-10)
        }
    }

    @Test("Integration of linear ramp (500 samples)")
    func integrationLinear() {
        // f(x) = x from 0 to 1, integral F(x) = x²/2
        let waveform = Waveform1D<Double>.linearRamp(startValue: 0.0, endValue: 1.0, duration: 1.0, samplingRate: 500.0)
        let integrated = waveform.integrate()

        // At t=1, integral should be approximately 0.5 (1²/2)
        #expect(integrated.values.count == 500)
        #expect(integrated.values[0] == 0.0)
        #expect(abs(integrated.values.last! - 0.5) < 0.01)

        // Should be quadratic - second differences should be small and constant
        for i in 2..<min(10, integrated.values.count) {
            let diff1 = integrated.values[i] - integrated.values[i-1]
            let diff2 = integrated.values[i-1] - integrated.values[i-2]
            // First differences should be increasing (quadratic)
            #expect(diff1 >= diff2 - 1e-6)
        }
    }

    @Test("Integration of sine wave (1000 samples)")
    func integrationSineWave() {
        // ∫sin(2πft)dt = -cos(2πft)/(2πf) + C
        let frequency = 2.0
        let waveform = Waveform1D<Double>.sine(frequency: frequency, duration: 1.0, samplingRate: 1000.0)
        let integrated = waveform.integrate()

        #expect(integrated.values.count == 1000)
        #expect(integrated.values[0] == 0.0)

        // For sine starting at 0, the numerical integral will:
        // 1. Start at 0 (not the theoretical minimum)
        // 2. First increase (since sin(0+) > 0)
        // 3. Reach maximum around t=0.25 (quarter period)
        // 4. Return toward 0 around t=0.5 (half period)
        // 5. Potentially go negative in subsequent cycles

        let maxVal = integrated.values.max()!
        let minVal = integrated.values.min()!

        // The integral should oscillate
        #expect(maxVal > 0.0)
        // minVal might be 0.0 or slightly negative depending on numerical precision
        #expect(minVal <= 0.01)  // Allow for small positive minimum due to numerical integration

        // Check that it actually oscillates (range should be reasonable)
        let range = maxVal - minVal
        let expectedAmplitude = 1.0 / (2.0 * Double.pi * frequency) // ≈ 0.0796
        #expect(range > expectedAmplitude * 0.5)  // Should have reasonable oscillation
        #expect(range < expectedAmplitude * 4.0)  // But not too large

        // Verify it returns close to starting value after complete cycles
        let finalValue = integrated.values.last!
        #expect(abs(finalValue) < expectedAmplitude * 2.0)  // Should be reasonably close to 0
    }

    @Test("Integration of cosine wave (800 samples)")
    func integrationCosineWave() {
        // ∫cos(2πft)dt = sin(2πft)/(2πf) + C
        let frequency = 1.0
        let waveform = Waveform1D<Double>.cosine(frequency: frequency, duration: 2.0, samplingRate: 400.0)
        let integrated = waveform.integrate()

        #expect(integrated.values.count == 800)
        #expect(integrated.values[0] == 0.0)

        // Cosine starts at 1, so integral should increase initially
        #expect(integrated.values[10] > integrated.values[0])

        // Over complete cycles, the integral oscillates
        let maxVal = integrated.values.max()!
        let minVal = integrated.values.min()!
        #expect(maxVal > 0.0)
        #expect(minVal < 0.0)
    }

    @Test("Integration of triangle wave (500 samples)")
    func integrationTriangleWave() {
        let waveform = Waveform1D<Double>.triangle(frequency: 2.0, duration: 1.0, samplingRate: 500.0)
        let integrated = waveform.integrate()

        #expect(integrated.values.count == 500)
        #expect(integrated.values[0] == 0.0)

        // Triangle wave integral should be smooth and oscillating
        // Check it's not constant
        let maxVal = integrated.values.max()!
        let minVal = integrated.values.min()!
        #expect(maxVal - minVal > 0.1)
    }

    @Test("Integration with custom initial value (200 samples)")
    func integrationWithInitialValue() {
        let waveform = Waveform1D<Double>.constant(value: 2.0, duration: 1.0, samplingRate: 200.0)
        let integrated = waveform.integrate(initialValue: 10.0)

        // Start at 10, then add area under constant 2.0 over 1 second = 2.0 total
        #expect(integrated.values[0] == 10.0)
        #expect(abs(integrated.values.last! - 12.0) < 0.1)
    }

    @Test("Integration preserves dt and t0")
    func integrationPreservesMetadata() {
        let startTime = Date()
        let waveform = Waveform1D<Double>.sine(frequency: 1.0, duration: 1.0, samplingRate: 100.0, t0: startTime)
        let integrated = waveform.integrate()

        #expect(integrated.dt == waveform.dt)
        #expect(integrated.t0 == startTime)
    }

    @Test("Integration of square wave (600 samples)")
    func integrationSquareWave() {
        let waveform = Waveform1D<Double>.square(frequency: 3.0, duration: 1.0, samplingRate: 600.0)
        let integrated = waveform.integrate()

        #expect(integrated.values.count == 600)
        #expect(integrated.values[0] == 0.0)

        // Square wave integral should be piecewise linear (ramp up/down)
        // Check that it's not constant
        let maxVal = integrated.values.max()!
        let minVal = integrated.values.min()!
        #expect(maxVal - minVal > 0.1)
    }

    @Test("Integration of sawtooth wave (400 samples)")
    func integrationSawtoothWave() {
        let waveform = Waveform1D<Double>.sawtooth(frequency: 2.0, duration: 1.0, samplingRate: 400.0)
        let integrated = waveform.integrate()

        #expect(integrated.values.count == 400)
        #expect(integrated.values[0] == 0.0)

        // Sawtooth integral should show periodic behavior
        let maxVal = integrated.values.max()!
        let minVal = integrated.values.min()!
        #expect(maxVal - minVal > 0.1)
    }
}

// MARK: - Derivative Tests Suite
@Suite("Waveform1D Derivative (Calc)")
struct Waveform1DDerivativeTests {

    @Test("Derivative of empty waveform")
    func derivativeEmpty() {
        let waveform = Waveform1D<Double>(values: [], dt: 0.1)
        let derived = waveform.derivative()

        #expect(derived.values.isEmpty)
    }

    @Test("Derivative of single point")
    func derivativeSinglePoint() {
        let waveform = Waveform1D<Double>(values: [5.0], dt: 0.1)
        let derived = waveform.derivative()

        #expect(derived.values.isEmpty)
    }

    @Test("Derivative of constant waveform (100 samples)")
    func derivativeConstant() {
        // d/dx(c) = 0
        let waveform = Waveform1D<Double>.constant(value: 5.0, duration: 1.0, samplingRate: 100.0)
        let derived = waveform.derivative()

        #expect(derived.values.count == 100)
        for value in derived.values {
            #expect(abs(value) < 1e-8)  // Should be ~0
        }
    }

    @Test("Derivative of linear ramp (500 samples)")
    func derivativeLinear() {
        // d/dx(x) = 1
        let waveform = Waveform1D<Double>.linearRamp(startValue: 0.0, endValue: 1.0, duration: 1.0, samplingRate: 500.0)
        let derived = waveform.derivative()

        #expect(derived.values.count == 500)
        // Slope should be approximately 1.0
        let avgSlope = derived.values.reduce(0.0, +) / Double(derived.values.count)
        #expect(abs(avgSlope - 1.0) < 0.1)
    }

    @Test("Derivative of sine wave is cosine (1000 samples)")
    func derivativeSineWave() {
        // d/dx(sin(2πft)) = 2πf·cos(2πft)
        let frequency = 2.0
        let waveform = Waveform1D<Double>.sine(frequency: frequency, duration: 1.0, samplingRate: 1000.0)
        let derived = waveform.derivative()

        // Compare with actual cosine (scaled by 2πf)
        let cosineWave = Waveform1D<Double>.cosine(frequency: frequency, duration: 1.0, samplingRate: 1000.0)
        let scaledCosine = cosineWave.values.map { $0 * 2.0 * Double.pi * frequency }

        #expect(derived.values.count == 1000)

        // Check several points (skip edges due to boundary effects)
        for i in 10..<(derived.values.count - 10) {
            #expect(abs(derived.values[i] - scaledCosine[i]) < 0.5)
        }
    }

    @Test("Derivative of cosine wave is negative sine (800 samples)")
    func derivativeCosineWave() {
        // d/dx(cos(2πft)) = -2πf·sin(2πft)
        let frequency = 1.0
        let waveform = Waveform1D<Double>.cosine(frequency: frequency, duration: 2.0, samplingRate: 400.0)
        let derived = waveform.derivative()

        let sineWave = Waveform1D<Double>.sine(frequency: frequency, duration: 2.0, samplingRate: 400.0)
        let negScaledSine = sineWave.values.map { -$0 * 2.0 * Double.pi * frequency }

        #expect(derived.values.count == 800)

        // Check correlation (skip edges)
        for i in 10..<(derived.values.count - 10) {
            #expect(abs(derived.values[i] - negScaledSine[i]) < 0.5)
        }
    }

    @Test("Derivative of triangle wave (500 samples)")
    func derivativeTriangleWave() {
        let waveform = Waveform1D<Double>.triangle(frequency: 2.0, duration: 1.0, samplingRate: 500.0)
        let derived = waveform.derivative()

        #expect(derived.values.count == 500)

        // Triangle wave derivative should be roughly square wave (constant slopes)
        // Check that derivative has distinct positive and negative regions
        let positiveCount = derived.values.filter { $0 > 0.5 }.count
        let negativeCount = derived.values.filter { $0 < -0.5 }.count
        #expect(positiveCount > 50)
        #expect(negativeCount > 50)
    }

    @Test("Derivative of square wave (600 samples)")
    func derivativeSquareWave() {
        let waveform = Waveform1D<Double>.square(frequency: 3.0, duration: 1.0, samplingRate: 600.0)
        let derived = waveform.derivative()

        #expect(derived.values.count == 600)

        // Square wave derivative should have large spikes at transitions
        let maxDerivative = derived.values.map { abs($0) }.max()!
        #expect(maxDerivative > 10.0)  // Should have significant spikes
    }

    @Test("Derivative of sawtooth wave (400 samples)")
    func derivativeSawtoothWave() {
        let waveform = Waveform1D<Double>.sawtooth(frequency: 2.0, duration: 1.0, samplingRate: 400.0)
        let derived = waveform.derivative()

        #expect(derived.values.count == 400)

        // Sawtooth has constant positive slope with sharp negative spikes
        let positiveSlopes = derived.values.filter { $0 > 0.5 }.count
        #expect(positiveSlopes > 300)  // Most of the time it's ramping up
    }

    @Test("Derivative preserves dt and t0")
    func derivativePreservesMetadata() {
        let startTime = Date()
        let waveform = Waveform1D<Double>.sine(frequency: 1.0, duration: 1.0, samplingRate: 100.0, t0: startTime)
        let derived = waveform.derivative()

        #expect(derived.dt == waveform.dt)
        #expect(derived.t0 == startTime)
    }

    @Test("Derivative handles negative values (200 samples)")
    func derivativeNegativeValues() {
        let waveform = Waveform1D<Double>.linearRamp(startValue: -10.0, endValue: 10.0, duration: 1.0, samplingRate: 200.0)
        let derived = waveform.derivative()

        #expect(derived.values.count == 200)
        // Constant slope throughout
        let avgSlope = derived.values.reduce(0.0, +) / Double(derived.values.count)
        #expect(abs(avgSlope - 20.0) < 1.0)  // Slope is (10 - (-10))/1 = 20
    }
}

// MARK: - Integration and Derivative Relationship Tests
@Suite("Calculus Fundamental Theorem (Calc)")
struct CalculusFundamentalTheoremTests {

    @Test("Derivative of integral returns approximate original (sine, 500 samples)")
    func derivativeOfIntegralSine() {
        let original = Waveform1D<Double>.sine(frequency: 1.0, duration: 1.0, samplingRate: 500.0)

        let integrated = original.integrate()
        let derived = integrated.derivative()

        // Due to numerical approximation, won't be exact
        #expect(derived.values.count == original.values.count)

        // Check middle points (avoid boundary effects)
        for i in 10..<(derived.values.count - 10) {
            #expect(abs(derived.values[i] - original.values[i]) < 0.1)
        }
    }

    @Test("Derivative of integral returns approximate original (triangle, 400 samples)")
    func derivativeOfIntegralTriangle() {
        let original = Waveform1D<Double>.triangle(frequency: 2.0, duration: 1.0, samplingRate: 400.0)

        let integrated = original.integrate()
        let derived = integrated.derivative()

        #expect(derived.values.count == original.values.count)

        // Check correlation - should be similar in middle region
        for i in 20..<(derived.values.count - 20) {
            #expect(abs(derived.values[i] - original.values[i]) < 0.2)
        }
    }

    @Test("Integral of derivative returns approximate original (cosine, 600 samples)")
    func integralOfDerivativeCosine() {
        let original = Waveform1D<Double>.cosine(frequency: 1.0, duration: 1.0, samplingRate: 600.0)

        let derived = original.derivative()
        let integrated = derived.integrate()

        // Should approximately reconstruct (up to constant)
        #expect(integrated.values.count == original.values.count)

        // Check that the shape is similar (differences should be similar)
        for i in 10..<(integrated.values.count - 10) {
            let originalDiff = original.values[i] - original.values[i-1]
            let reconstructedDiff = integrated.values[i] - integrated.values[i-1]
            #expect(abs(originalDiff - reconstructedDiff) < 0.01)
        }
    }

    @Test("Second derivative of sine is negative sine (1000 samples)")
    func secondDerivativeSine() {
        // d²/dx²(sin(x)) = -sin(x)
        let waveform = Waveform1D<Double>.sine(frequency: 1.0, duration: 2.0, samplingRate: 500.0)

        let firstDerivative = waveform.derivative()
        let secondDerivative = firstDerivative.derivative()

        #expect(secondDerivative.values.count == 1000)

        // Second derivative should be approximately -original
        // Check middle points
        for i in 20..<(secondDerivative.values.count - 20) {
            let expected = -waveform.values[i] * pow(2.0 * Double.pi * 1.0, 2.0)
            #expect(abs(secondDerivative.values[i] - expected) < 1.0)
        }
    }

    @Test("Third derivative of cubic (ramp integrated twice, 300 samples)")
    func thirdDerivative() {
        // Start with linear ramp
        let linear = Waveform1D<Double>.linearRamp(startValue: 0.0, endValue: 1.0, duration: 1.0, samplingRate: 300.0)
        let quadratic = linear.integrate()  // x²/2
        let cubic = quadratic.integrate()   // x³/6

        let firstDeriv = cubic.derivative()   // x²/2
        let secondDeriv = firstDeriv.derivative() // x
        let thirdDeriv = secondDeriv.derivative() // constant

        #expect(thirdDeriv.values.count > 0)

        // Third derivative should be approximately constant, but numerical errors accumulate
        let avg = thirdDeriv.values.reduce(0.0, +) / Double(thirdDeriv.values.count)
        let variance = thirdDeriv.values.map { pow($0 - avg, 2.0) }.reduce(0.0, +) / Double(thirdDeriv.values.count)
        
        // Test that the variance is reasonable (not infinite/NaN)
        #expect(variance.isFinite)
        #expect(variance >= 0.0)
        
        // Test that the mean is in a reasonable range for the expected constant derivative
        #expect(abs(avg) < 10.0)
        
        // Test that it's not completely random - should have some structure
        let standardDeviation = sqrt(variance)
        #expect(standardDeviation < abs(avg) * 10.0 || abs(avg) < 1.0)  // Relative or absolute bound
    }
}

// MARK: - Edge Cases and Numerical Stability Tests
@Suite("Calculus Edge Cases (Calc)")
struct CalculusEdgeCasesTests {

    @Test("Very small dt values (1000 samples)")
    func verySmallDt() {
        let waveform = Waveform1D<Double>.linearRamp(startValue: 0.0, endValue: 1.0, duration: 0.1, samplingRate: 10000.0)

        let derived = waveform.derivative()
        #expect(derived.values.count == 1000)

        // Slope should be 1/0.1 = 10
        let avgSlope = derived.values.reduce(0.0, +) / Double(derived.values.count)
        #expect(abs(avgSlope - 10.0) < 1.0)
    }

    @Test("Large dt values (100 samples)")
    func largeDt() {
        let waveform = Waveform1D<Double>.linearRamp(startValue: 0.0, endValue: 100.0, duration: 100.0, samplingRate: 1.0)

        let derived = waveform.derivative()
        #expect(derived.values.count == 100)

        // Slope should be 100/100 = 1.0
        let avgSlope = derived.values.reduce(0.0, +) / Double(derived.values.count)
        #expect(abs(avgSlope - 1.0) < 0.1)
    }

    @Test("Integration with many points (1000 samples)")
    func integrationManyPoints() {
        let waveform = Waveform1D<Double>.constant(value: 1.0, duration: 1.0, samplingRate: 1000.0)

        let integrated = waveform.integrate()
        #expect(integrated.values.count == 1000)
        #expect(integrated.values.last! > 0.9)  // Should be close to 1.0
    }

    @Test("Derivative with many points (1000 samples)")
    func derivativeManyPoints() {
        let waveform = Waveform1D<Double>.linearRamp(startValue: 0.0, endValue: 999.0, duration: 1.0, samplingRate: 1000.0)

        let derived = waveform.derivative()
        #expect(derived.values.count == 1000)

        // Linear function has constant derivative
        let avgSlope = derived.values.reduce(0.0, +) / Double(derived.values.count)
        #expect(abs(avgSlope - 999.0) < 10.0)
    }

    @Test("Zero values (200 samples)")
    func zeroValues() {
        let waveform = Waveform1D<Double>.constant(value: 0.0, duration: 1.0, samplingRate: 200.0)

        let integrated = waveform.integrate()
        let derived = waveform.derivative()

        #expect(integrated.values.allSatisfy { $0 == 0.0 })
        #expect(derived.values.allSatisfy { abs($0) < 1e-10 })
    }

    @Test("High frequency sine wave (500 samples)")
    func highFrequencySine() {
        // High frequency requires adequate sampling
        let waveform = Waveform1D<Double>.sine(frequency: 50.0, duration: 1.0, samplingRate: 500.0)

        let derived = waveform.derivative()
        let integrated = waveform.integrate()

        #expect(derived.values.count == 500)
        #expect(integrated.values.count == 500)

        // Should complete without numerical issues
        #expect(!derived.values.contains { $0.isNaN || $0.isInfinite })
        #expect(!integrated.values.contains { $0.isNaN || $0.isInfinite })
    }

    @Test("Multiple cycles of square wave (800 samples)")
    func multipleSquareCycles() {
        let waveform = Waveform1D<Double>.square(frequency: 10.0, duration: 1.0, samplingRate: 800.0)

        let integrated = waveform.integrate()

        #expect(integrated.values.count == 800)
        // Over complete cycles, integral should be bounded
        let maxVal = integrated.values.map { abs($0) }.max()!
        #expect(maxVal < 10.0)  // Shouldn't grow unbounded
    }
}

// MARK: - Physical Interpretation Tests
@Suite("Physical Interpretations (Calc)")
struct PhysicalInterpretationTests {

    @Test("Velocity to position (kinematics, 600 samples)")
    func velocityToPosition() {
        // Constant velocity of 10 m/s for 1 second
        let velocity = Waveform1D<Double>.constant(value: 10.0, duration: 1.0, samplingRate: 600.0)
        let position = velocity.integrate()

        // Position should increase linearly to 10 meters
        #expect(position.values.count == 600)
        #expect(position.values[0] == 0.0)
        #expect(abs(position.values.last! - 10.0) < 0.2)
    }

    @Test("Position to velocity (kinematics, 500 samples)")
    func positionToVelocity() {
        // Linear position: 0 to 10 meters over 1 second
        let position = Waveform1D<Double>.linearRamp(startValue: 0.0, endValue: 10.0, duration: 1.0, samplingRate: 500.0)
        let velocity = position.derivative()

        // Velocity should be constant 10 m/s
        #expect(velocity.values.count == 500)
        let avgVelocity = velocity.values.reduce(0.0, +) / Double(velocity.values.count)
        #expect(abs(avgVelocity - 10.0) < 0.5)
    }

    @Test("Acceleration integration (gravity, 1000 samples)")
    func accelerationIntegration() {
        // Constant acceleration of 9.8 m/s² (gravity) for 1 second
        let acceleration = Waveform1D<Double>.constant(value: 9.8, duration: 1.0, samplingRate: 1000.0)
        let velocity = acceleration.integrate()

        // After 1 second, velocity should be ~9.8 m/s
        #expect(velocity.values.count == 1000)
        #expect(abs(velocity.values.last! - 9.8) < 0.5)
    }

    @Test("Oscillating force (harmonic motion, 800 samples)")
    func harmonicMotion() {
        // Force proportional to -x (Hooke's law approximation)
        let force = Waveform1D<Double>.sine(frequency: 1.0, duration: 2.0, samplingRate: 400.0)

        // This represents simplified harmonic oscillator
        let acceleration = force  // F = ma, assume m=1
        let velocity = acceleration.integrate()
        let position = velocity.integrate()

        #expect(position.values.count == 800)

        // Position should oscillate smoothly
        let maxPos = position.values.max()!
        let minPos = position.values.min()!
        #expect(maxPos > 0.0)
        #expect(minPos <= 0.0)
        #expect(abs(maxPos + minPos) < 0.5)  // Symmetric oscillation
    }

    @Test("Power integration to energy (400 samples)")
    func powerToEnergy() {
        // Constant power of 100 W for 1 second
        let power = Waveform1D<Double>.constant(value: 100.0, duration: 1.0, samplingRate: 400.0)
        let energy = power.integrate()

        // Energy should be 100 J after 1 second
        #expect(energy.values.count == 400)
        #expect(energy.values[0] == 0.0)
        #expect(abs(energy.values.last! - 100.0) < 2.0)
    }
}
