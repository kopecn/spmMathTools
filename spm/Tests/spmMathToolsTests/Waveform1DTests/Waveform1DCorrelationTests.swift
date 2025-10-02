import Foundation
import Testing

@testable import spmMathTools

// MARK: - Basic Auto-Correlation Suite
@Suite("Auto-Correlation Tests")
struct AutoCorrelationTests {

    @Test("Basic auto-correlation")
    func basicAutoCorrelation() {
        let values = [1.0, 2.0, 3.0, 2.0, 1.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let autoCorr = waveform.autoCorrelation()

        // Auto-correlation should have maximum at lag 0
        #expect(!autoCorr.values.isEmpty)
        #expect(autoCorr.values[0] > 0.0)  // Lag 0 should be positive
        #expect(autoCorr.dt == waveform.dt)
        #expect(autoCorr.t0 == nil)
    }

    @Test("Auto-correlation of sine wave")
    func autoCorrelationSineWave() {
        let waveform = Waveform1D<Double>.sine(frequency: 1.0, duration: 2.0, samplingRate: 100.0)
        let autoCorr = waveform.autoCorrelation(maxLag: 50)

        #expect(autoCorr.values.count == 51)  // 0 to 50 lags
        #expect(autoCorr.values[0] > 0.0)  // Maximum at lag 0

        // For sine wave, auto-correlation should be periodic
        // Find secondary peaks that indicate periodicity
        let maxVal = autoCorr.values[0]
        let secondaryPeaks = autoCorr.values.enumerated().filter { index, value in
            index > 10 && value > maxVal * 0.5
        }
        #expect(!secondaryPeaks.isEmpty)
    }

    @Test("Auto-correlation with custom max lag")
    func autoCorrelationCustomMaxLag() {
        let values = Array(0..<20).map { Double($0) }
        let waveform = Waveform1D(values: values, dt: 0.1)

        let autoCorr = waveform.autoCorrelation(maxLag: 5)

        #expect(autoCorr.values.count == 6)  // 0 to 5 lags
    }

    @Test("Auto-correlation normalized vs non-normalized")
    func autoCorrelationNormalization() {
        let values = [1.0, 4.0, 2.0, 3.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let normalizedCorr = waveform.autoCorrelation(normalized: true)
        let unnormalizedCorr = waveform.autoCorrelation(normalized: false)

        #expect(normalizedCorr.values.count == unnormalizedCorr.values.count)
        #expect(normalizedCorr.values[0] <= 1.0)  // Normalized should be ≤ 1
        #expect(unnormalizedCorr.values[0] > normalizedCorr.values[0])  // Unnormalized should be larger
    }

    @Test("Auto-correlation preserves metadata")
    func autoCorrelationMetadata() {
        let values = [1.0, 2.0, 3.0, 4.0]
        let dt = 0.25
        let waveform = Waveform1D(values: values, dt: dt)

        let autoCorr = waveform.autoCorrelation()

        #expect(autoCorr.dt == dt)
        #expect(autoCorr.t0 == nil)
    }
}

// MARK: - Cross-Correlation Suite
@Suite("Cross-Correlation Tests")
struct CrossCorrelationTests {

    @Test("Basic cross-correlation")
    func basicCrossCorrelation() {
        let values1 = [1.0, 2.0, 3.0]
        let values2 = [1.0, 1.0, 1.0]
        let waveform1 = Waveform1D(values: values1, dt: 0.1)
        let waveform2 = Waveform1D(values: values2, dt: 0.1)

        let crossCorr = waveform1.crossCorrelation(with: waveform2, mode: .full)

        #expect(crossCorr != nil)
        #expect(crossCorr!.values.count == 5)  // 3 + 3 - 1
        #expect(crossCorr!.dt == waveform1.dt)
    }

    @Test("Cross-correlation of identical signals")
    func crossCorrelationIdentical() {
        let values = [1.0, 2.0, 3.0, 2.0, 1.0]
        let waveform1 = Waveform1D(values: values, dt: 0.1)
        let waveform2 = Waveform1D(values: values, dt: 0.1)

        let crossCorr = waveform1.crossCorrelation(with: waveform2, mode: .full, normalized: true)

        #expect(crossCorr != nil)
        // Maximum correlation should be at the center (zero lag)
        let centerIndex = crossCorr!.values.count / 2
        let maxValue = crossCorr!.values.max()!
        #expect(abs(crossCorr!.values[centerIndex] - maxValue) < 1e-10)
        #expect(maxValue <= 1.0)  // Normalized should be ≤ 1
    }

    @Test("Cross-correlation with time shift")
    func crossCorrelationTimeShift() {
        let baseSignal = [0.0, 1.0, 2.0, 1.0, 0.0]
        let shiftedSignal = [1.0, 2.0, 1.0, 0.0, 0.0]  // Shifted left by 1
        
        let waveform1 = Waveform1D(values: baseSignal, dt: 0.1)
        let waveform2 = Waveform1D(values: shiftedSignal, dt: 0.1)

        let crossCorr = waveform1.crossCorrelation(with: waveform2, mode: .full, normalized: true)

        #expect(crossCorr != nil)
        
        // Find the maximum correlation
        let maxIndex = crossCorr!.values.enumerated().max { $0.element < $1.element }!.offset
        let expectedCenterIndex = baseSignal.count - 1
        
        // The shift should be detectable in the correlation
        #expect(maxIndex != expectedCenterIndex)
    }

    @Test("Cross-correlation with sine waves")
    func crossCorrelationSineWaves() {
        let waveform1 = Waveform1D<Double>.sine(frequency: 2.0, duration: 1.0, samplingRate: 100.0)
        let waveform2 = Waveform1D<Double>.sine(frequency: 2.0, phase: Double.pi / 2, duration: 1.0, samplingRate: 100.0)

        let crossCorr = waveform1.crossCorrelation(with: waveform2, mode: .same, normalized: true)

        #expect(crossCorr != nil)
        #expect(crossCorr!.values.count == waveform1.values.count)
        
        // Sine and cosine should have specific correlation properties
        let maxCorrelation = crossCorr!.values.map { abs($0) }.max()!
        #expect(maxCorrelation > 0.5)  // Should have reasonable correlation
    }
}

// MARK: - Correlation Modes Suite
@Suite("Correlation Modes")
struct CorrelationModesTests {

    @Test("Full correlation mode")
    func fullCorrelationMode() {
        let values1 = [1.0, 2.0, 3.0]
        let values2 = [1.0, 1.0]
        let waveform1 = Waveform1D(values: values1, dt: 0.1)
        let waveform2 = Waveform1D(values: values2, dt: 0.1)

        let crossCorr = waveform1.crossCorrelation(with: waveform2, mode: .full)

        #expect(crossCorr != nil)
        #expect(crossCorr!.values.count == 4)  // 3 + 2 - 1
    }

    @Test("Valid correlation mode")
    func validCorrelationMode() {
        let values1 = [1.0, 2.0, 3.0, 4.0, 5.0]
        let values2 = [1.0, 1.0, 1.0]
        let waveform1 = Waveform1D(values: values1, dt: 0.1)
        let waveform2 = Waveform1D(values: values2, dt: 0.1)

        let crossCorr = waveform1.crossCorrelation(with: waveform2, mode: .valid)

        #expect(crossCorr != nil)
        #expect(crossCorr!.values.count == 3)  // 5 - 3 + 1
    }

    @Test("Same correlation mode")
    func sameCorrelationMode() {
        let values1 = [1.0, 2.0, 3.0, 4.0]
        let values2 = [1.0, 1.0]
        let waveform1 = Waveform1D(values: values1, dt: 0.1)
        let waveform2 = Waveform1D(values: values2, dt: 0.1)

        let crossCorr = waveform1.crossCorrelation(with: waveform2, mode: .same)

        #expect(crossCorr != nil)
        #expect(crossCorr!.values.count == 4)  // Same as first input
    }

    @Test("Valid mode with incompatible lengths")
    func validModeIncompatibleLengths() {
        let values1 = [1.0, 2.0]
        let values2 = [1.0, 1.0, 1.0, 1.0]  // Longer than first signal
        let waveform1 = Waveform1D(values: values1, dt: 0.1)
        let waveform2 = Waveform1D(values: values2, dt: 0.1)

        let crossCorr = waveform1.crossCorrelation(with: waveform2, mode: .valid)

        #expect(crossCorr == nil)  // Should return nil for incompatible lengths
    }
}

// MARK: - Maximum Correlation Finding Suite
@Suite("Maximum Correlation Finding")
struct MaxCorrelationTests {

    @Test("Find max correlation with identical signals")
    func findMaxCorrelationIdentical() {
        let values = [1.0, 3.0, 2.0, 4.0, 1.0]
        let waveform1 = Waveform1D(values: values, dt: 0.1)
        let waveform2 = Waveform1D(values: values, dt: 0.1)

        let result = waveform1.findMaxCorrelation(with: waveform2)

        #expect(result != nil)
        #expect(result!.lagSamples == 0)  // No lag for identical signals
        #expect(result!.lagTime == 0.0)
        #expect(result!.correlation > 0.9)  // Should be close to 1.0
    }

    @Test("Find max correlation with shifted signal")
    func findMaxCorrelationShifted() {
        let original = [0.0, 0.0, 1.0, 2.0, 1.0, 0.0, 0.0]
        let shifted = [0.0, 1.0, 2.0, 1.0, 0.0, 0.0, 0.0]  // Shifted left by 1
        
        let waveform1 = Waveform1D(values: original, dt: 0.2)
        let waveform2 = Waveform1D(values: shifted, dt: 0.2)

        let result = waveform1.findMaxCorrelation(with: waveform2)

        #expect(result != nil)
        #expect(result!.lagSamples == -1)  // Shifted left by 1 sample
        #expect(result!.lagTime == -0.2)   // Shifted left by 0.2 seconds
        #expect(result!.correlation > 0.5)
    }

    @Test("Find max correlation with search range")
    func findMaxCorrelationSearchRange() {
        let values1 = [1.0, 2.0, 3.0, 4.0, 5.0]
        let values2 = [3.0, 4.0, 5.0]
        let waveform1 = Waveform1D(values: values1, dt: 0.1)
        let waveform2 = Waveform1D(values: values2, dt: 0.1)

        let result = waveform1.findMaxCorrelation(with: waveform2, searchRange: 1..<4)

        #expect(result != nil)
        // Verify the result is within the search range
        let crossCorr = waveform1.crossCorrelation(with: waveform2, mode: .full)!
        let actualLag = result!.lagSamples + (values2.count - 1)
        #expect(actualLag >= 1 && actualLag < 4)
    }

    @Test("Find max correlation with invalid search range")
    func findMaxCorrelationInvalidRange() {
        let values = [1.0, 2.0, 3.0]
        let waveform1 = Waveform1D(values: values, dt: 0.1)
        let waveform2 = Waveform1D(values: values, dt: 0.1)

        let result = waveform1.findMaxCorrelation(with: waveform2, searchRange: 10..<20)

        #expect(result == nil)  // Invalid range should return nil
    }
}

// MARK: - Edge Cases Suite
@Suite("Correlation Edge Cases")
struct CorrelationEdgeCasesTests {

    @Test("Auto-correlation with empty signal")
    func autoCorrelationEmpty() {
        let waveform = Waveform1D(values: [Double](), dt: 0.1)

        let autoCorr = waveform.autoCorrelation()

        #expect(autoCorr.values.isEmpty)
        #expect(autoCorr.dt == 0.1)
    }

    @Test("Cross-correlation with empty signals")
    func crossCorrelationEmpty() {
        let waveform1 = Waveform1D(values: [Double](), dt: 0.1)
        let waveform2 = Waveform1D(values: [1.0, 2.0], dt: 0.1)

        let crossCorr1 = waveform1.crossCorrelation(with: waveform2)
        let crossCorr2 = waveform2.crossCorrelation(with: waveform1)

        #expect(crossCorr1 == nil)
        #expect(crossCorr2 == nil)
    }

    @Test("Auto-correlation with single value")
    func autoCorrelationSingle() {
        let waveform = Waveform1D(values: [5.0], dt: 0.1)

        let autoCorr = waveform.autoCorrelation()

        #expect(autoCorr.values.count == 1)
        #expect(autoCorr.values[0] > 0.0)
    }

    @Test("Cross-correlation with single values")
    func crossCorrelationSingle() {
        let waveform1 = Waveform1D(values: [3.0], dt: 0.1)
        let waveform2 = Waveform1D(values: [2.0], dt: 0.1)

        let crossCorr = waveform1.crossCorrelation(with: waveform2, mode: .full)

        #expect(crossCorr != nil)
        #expect(crossCorr!.values.count == 1)
    }

    @Test("Auto-correlation with constant signal")
    func autoCorrelationConstant() {
        let waveform = Waveform1D(values: [5.0, 5.0, 5.0, 5.0], dt: 0.1)

        let autoCorr = waveform.autoCorrelation(normalized: true)

        #expect(!autoCorr.values.isEmpty)
        // For normalized constant signal, correlation should be 1.0 at all lags
        for value in autoCorr.values {
            #expect(abs(value - 1.0) < 1e-10)
        }
    }

    @Test("Cross-correlation with zero variance signal")
    func crossCorrelationZeroVariance() {
        let waveform1 = Waveform1D(values: [1.0, 2.0, 3.0], dt: 0.1)
        let waveform2 = Waveform1D(values: [5.0, 5.0, 5.0], dt: 0.1)

        let crossCorr = waveform1.crossCorrelation(with: waveform2, normalized: true)

        #expect(crossCorr != nil)
        #expect(!crossCorr!.values.isEmpty)
    }

    @Test("Auto-correlation with zero values")
    func autoCorrelationZeros() {
        let waveform = Waveform1D(values: [0.0, 0.0, 0.0, 0.0], dt: 0.1)

        let autoCorr = waveform.autoCorrelation()

        #expect(!autoCorr.values.isEmpty)
        for value in autoCorr.values {
            #expect(value == 0.0)
        }
    }
}

// MARK: - Mathematical Properties Suite
@Suite("Correlation Mathematical Properties")
struct CorrelationMathematicalPropertiesTests {

    @Test("Auto-correlation symmetry property")
    func autoCorrelationSymmetry() {
        let values = [1.0, 3.0, 2.0, 4.0, 1.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        // Create a symmetric auto-correlation by computing both directions
        let fullAutoCorr = waveform.crossCorrelation(with: waveform, mode: .full, normalized: true)!
        
        let centerIndex = fullAutoCorr.values.count / 2
        
        // Check symmetry around center
        for i in 1...min(2, centerIndex) {
            let leftValue = fullAutoCorr.values[centerIndex - i]
            let rightValue = fullAutoCorr.values[centerIndex + i]
            #expect(abs(leftValue - rightValue) < 1e-10)
        }
    }

    @Test("Cross-correlation maximum at zero lag for identical signals")
    func crossCorrelationMaxAtZeroLag() {
        let values = [1.0, 4.0, 2.0, 3.0, 1.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let crossCorr = waveform.crossCorrelation(with: waveform, mode: .full, normalized: true)!
        
        let centerIndex = crossCorr.values.count / 2
        let maxValue = crossCorr.values.max()!
        
        #expect(abs(crossCorr.values[centerIndex] - maxValue) < 1e-10)
    }

    @Test("Normalized correlation bounds")
    func normalizedCorrelationBounds() {
        let values1 = [1.0, 3.0, 2.0, 4.0, 1.0]
        let values2 = [2.0, 1.0, 3.0, 1.0, 2.0]
        let waveform1 = Waveform1D(values: values1, dt: 0.1)
        let waveform2 = Waveform1D(values: values2, dt: 0.1)

        let crossCorr = waveform1.crossCorrelation(with: waveform2, normalized: true)!

        // All values should be between -1 and 1
        for value in crossCorr.values {
            #expect(value >= -1.0)
            #expect(value <= 1.0)
        }
    }

    @Test("Auto-correlation at zero lag equals signal energy")
    func autoCorrelationZeroLagEnergy() {
        let values = [1.0, 2.0, 3.0, 2.0, 1.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let autoCorr = waveform.autoCorrelation(normalized: false)
        let signalEnergy = values.map { $0 * $0 }.reduce(0.0, +)

        #expect(abs(autoCorr.values[0] - signalEnergy) < 1e-10)
    }
}

// MARK: - Floating Point Types Suite
@Suite("Correlation Floating Point Types")
struct CorrelationFloatingPointTypesTests {

    @Test("Auto-correlation with Float")
    func autoCorrelationFloat() {
        let values: [Float] = [1.0, 2.0, 3.0, 2.0, 1.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let autoCorr = waveform.autoCorrelation()

        #expect(!autoCorr.values.isEmpty)
        #expect(autoCorr.values[0] > 0.0)
    }

    @Test("Cross-correlation with Float")
    func crossCorrelationFloat() {
        let values1: [Float] = [1.0, 2.0, 3.0]
        let values2: [Float] = [1.0, 1.0, 1.0]
        let waveform1 = Waveform1D(values: values1, dt: 0.1)
        let waveform2 = Waveform1D(values: values2, dt: 0.1)

        let crossCorr = waveform1.crossCorrelation(with: waveform2)

        #expect(crossCorr != nil)
        #expect(!crossCorr!.values.isEmpty)
    }

    @Test("Find max correlation with Float")
    func findMaxCorrelationFloat() {
        let values: [Float] = [1.0, 2.0, 3.0, 2.0, 1.0]
        let waveform1 = Waveform1D(values: values, dt: 0.1)
        let waveform2 = Waveform1D(values: values, dt: 0.1)

        let result = waveform1.findMaxCorrelation(with: waveform2)

        #expect(result != nil)
        #expect(result!.lagSamples == 0)
        #expect(result!.correlation > 0.9)
    }
}

// MARK: - Performance Tests Suite
@Suite("Correlation Performance")
struct CorrelationPerformanceTests {

    @Test("Auto-correlation performance with large dataset", .timeLimit(.minutes(1)))
    func autoCorrelationPerformance() {
        let values = (0..<5000).map { i in sin(Double(i) * 0.01) }
        let waveform = Waveform1D(values: values, dt: 0.01)

        let autoCorr = waveform.autoCorrelation(maxLag: 1000)

        #expect(autoCorr.values.count == 1001)
        #expect(autoCorr.values[0] > 0.0)
    }

    @Test("Cross-correlation performance with large datasets", .timeLimit(.minutes(1)))
    func crossCorrelationPerformance() {
        let values1 = (0..<2000).map { i in sin(Double(i) * 0.02) }
        let values2 = (0..<1000).map { i in cos(Double(i) * 0.02) }
        let waveform1 = Waveform1D(values: values1, dt: 0.01)
        let waveform2 = Waveform1D(values: values2, dt: 0.01)

        let crossCorr = waveform1.crossCorrelation(with: waveform2, mode: .valid)

        #expect(crossCorr != nil)
        #expect(crossCorr!.values.count == 1001)  // 2000 - 1000 + 1
    }

    @Test("Find max correlation performance", .timeLimit(.minutes(1)))
    func findMaxCorrelationPerformance() {
        let baseSignal = (0..<1000).map { i in sin(Double(i) * 0.05) + 0.1 * Double.random(in: -1...1) }
        let shiftedSignal = Array(baseSignal[50...]) + Array(repeating: 0.0, count: 50)
        
        let waveform1 = Waveform1D(values: baseSignal, dt: 0.01)
        let waveform2 = Waveform1D(values: shiftedSignal, dt: 0.01)

        let result = waveform1.findMaxCorrelation(with: waveform2)

        #expect(result != nil)
        #expect(abs(result!.lagSamples - (-50)) <= 5)  // Should detect the 50-sample shift within tolerance
    }
}

// MARK: - Real-World Applications Suite
@Suite("Correlation Real-World Applications")
struct CorrelationApplicationTests {

    @Test("Signal delay detection")
    func signalDelayDetection() {
        // Create a reference signal and a delayed version
        let reference = Waveform1D<Double>.triangle(frequency: 5.0, duration: 1.0, samplingRate: 100.0)
        
        // Create delayed signal by padding with zeros
        let delayInSamples = 10
        let delayedValues = Array(repeating: 0.0, count: delayInSamples) + reference.values
        let delayed = Waveform1D(values: delayedValues, dt: reference.dt)

        let result = reference.findMaxCorrelation(with: delayed)

        #expect(result != nil)
        #expect(abs(result!.lagSamples - delayInSamples) <= 2)  // Should detect the delay
        #expect(result!.correlation > 0.7)  // Should have good correlation
    }

    @Test("Periodic signal correlation")
    func periodicSignalCorrelation() {
        let waveform1 = Waveform1D<Double>.sine(frequency: 2.0, duration: 2.0, samplingRate: 50.0)
        let waveform2 = Waveform1D<Double>.sine(frequency: 2.0, duration: 2.0, samplingRate: 50.0)

        let autoCorr = waveform1.autoCorrelation(maxLag: 50)

        // For periodic signals, auto-correlation should show periodicity
        let period = Int(50.0 / 2.0)  // Expected period in samples
        
        // Check that correlation at period lag is high
        if period < autoCorr.values.count {
            #expect(autoCorr.values[period] > autoCorr.values[0] * 0.8)
        }
    }

    @Test("Noise correlation properties")
    func noiseCorrelationProperties() {
        // Generate white noise
        let noiseValues = (0..<200).map { _ in Double.random(in: -1...1) }
        let noiseWaveform = Waveform1D(values: noiseValues, dt: 0.01)

        let autoCorr = noiseWaveform.autoCorrelation(maxLag: 50, normalized: true)

        // For white noise, auto-correlation should decay quickly
        #expect(autoCorr.values[0] > 0.9)  // High at lag 0
        
        // Check that correlation decreases for higher lags
        if autoCorr.values.count > 10 {
            let avgHighLag = autoCorr.values[5...].reduce(0.0, +) / Double(autoCorr.values.count - 5)
            #expect(abs(avgHighLag) < 0.3)  // Should be much smaller than lag 0
        }
    }

    @Test("Mixed signal correlation")
    func mixedSignalCorrelation() {
        // Create a signal with both sine and noise components
        let pure = Waveform1D<Double>.sine(frequency: 3.0, duration: 1.0, samplingRate: 100.0)
        let noise = (0..<100).map { _ in 0.2 * Double.random(in: -1...1) }
        let mixed = pure.values.enumerated().map { index, value in value + noise[index] }
        let mixedWaveform = Waveform1D(values: mixed, dt: pure.dt)

        let crossCorr = pure.crossCorrelation(with: mixedWaveform, normalized: true)!

        // Should still have reasonable correlation despite noise
        let maxCorrelation = crossCorr.values.map { abs($0) }.max()!
        #expect(maxCorrelation > 0.7)
    }
}