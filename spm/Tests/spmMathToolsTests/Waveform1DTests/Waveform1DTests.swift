import Foundation
import Testing

@testable import spmMathTools

// MARK: - Initialization Tests Suite
@Suite("Waveform1D Initialization")
struct Waveform1DInitializationTests {

    @Test("Basic initialization with all parameters")
    func initializationWithAllParameters() {
        let startTime = Date()
        let values = [1.0, 2.0, 3.0, 4.0, 5.0]
        let dt = 0.001

        let waveform = Waveform1D(values: values, dt: dt, t0: startTime)

        #expect(waveform.values == values)
        #expect(waveform.dt == dt)
        #expect(waveform.t0 == startTime)
    }

    @Test("Initialization with values only")
    func initializationWithValuesOnly() {
        let values = [1, 2, 3, 4, 5]
        let waveform = Waveform1D(values: values)

        #expect(waveform.values == values)
        #expect(waveform.dt == 1.0)
        #expect(waveform.t0 == nil)
    }

    @Test("Initialization with values and dt")
    func initializationWithValuesAndDt() {
        let values = [1.0, 2.0, 3.0]
        let dt = 0.5
        let waveform = Waveform1D(values: values, dt: dt)

        #expect(waveform.values == values)
        #expect(waveform.dt == dt)
        #expect(waveform.t0 == nil)
    }

    @Test("Type aliases work correctly")
    func typeAliases() {
        let doubleWaveform = DoubleWaveform1D(values: [1.0, 2.0, 3.0])
        let floatWaveform = FloatWaveform1D(values: [1.0, 2.0, 3.0])
        let intWaveform = IntWaveform1D(values: [1, 2, 3])

        #expect(doubleWaveform.values.count == 3)
        #expect(floatWaveform.values.count == 3)
        #expect(intWaveform.values.count == 3)
    }
}

// MARK: - Computed Properties Suite
@Suite("Waveform1D Computed Properties")
struct Waveform1DComputedPropertiesTests {

    @Test("Duration calculation with multiple samples")
    func durationWithMultipleSamples() {
        let waveform = Waveform1D(values: [1.0, 2.0, 3.0, 4.0, 5.0], dt: 0.1)
        let expectedDuration = TimeInterval(4) * 0.1  // (5-1) * 0.1

        #expect(waveform.duration == expectedDuration)
    }

    @Test("Duration calculation with single sample")
    func durationWithSingleSample() {
        let waveform = Waveform1D(values: [1.0], dt: 0.1)

        #expect(waveform.duration == 0)
    }

    @Test("Duration calculation with empty array")
    func durationWithEmptyArray() {
        let waveform = Waveform1D(values: [Double](), dt: 0.1)

        #expect(waveform.duration == 0)
    }

    @Test("Sampling frequency calculation")
    func samplingFrequency() {
        let waveform = Waveform1D(values: [1.0, 2.0, 3.0], dt: 0.001)
        let expectedFrequency = 1.0 / 0.001

        #expect(waveform.samplingFrequency == expectedFrequency)
    }

    @Test("Nyquist frequency calculation")
    func nyquistFrequency() {
        let waveform = Waveform1D(values: [1.0, 2.0, 3.0], dt: 0.002)
        let expectedNyquist = (1.0 / 0.002) / 2.0

        #expect(waveform.nyquistFrequency == expectedNyquist)
    }

    @Test("End time calculation with t0")
    func endTimeWithT0() {
        let startTime = Date()
        let waveform = Waveform1D(values: [1.0, 2.0, 3.0, 4.0], dt: 0.1, t0: startTime)
        let expectedEndTime = startTime.addingTimeInterval(0.3)  // 3 * 0.1

        #expect(waveform.endTime == expectedEndTime)
    }

    @Test("End time calculation without t0")
    func endTimeWithoutT0() {
        let waveform = Waveform1D(values: [1.0, 2.0, 3.0], dt: 0.1)

        #expect(waveform.endTime == nil)
    }

    @Test("Sample count")
    func sampleCount() {
        let waveform = Waveform1D(values: [1, 2, 3, 4, 5, 6])

        #expect(waveform.sampleCount == 6)
    }
}

// MARK: - Comparable Extension Suite
@Suite("Waveform1D Comparable Operations")
struct Waveform1DComparableTests {

    @Test("Peak-to-peak calculation for comparable types")
    func peakToPeakComparable() {
        let waveform = Waveform1D(values: [1.0, 5.0, 2.0, 8.0, 3.0])
        let expectedPeakToPeak = 8.0 - 1.0

        #expect(waveform.peakToPeak == expectedPeakToPeak)
    }

    @Test("Peak-to-peak with empty array")
    func peakToPeakEmpty() {
        let waveform = Waveform1D(values: [Double]())

        #expect(waveform.peakToPeak == nil)
    }

    @Test("Minimum value calculation")
    func minimumValue() {
        let waveform = Waveform1D(values: [3, 1, 4, 1, 5])

        #expect(waveform.minimum == 1)
    }

    @Test("Maximum value calculation")
    func maximumValue() {
        let waveform = Waveform1D(values: [3, 1, 4, 1, 5])

        #expect(waveform.maximum == 5)
    }

    @Test("Min/Max with empty array")
    func minMaxEmpty() {
        let waveform = Waveform1D(values: [Int]())

        #expect(waveform.minimum == nil)
        #expect(waveform.maximum == nil)
    }
}

// MARK: - BinaryFloatingPoint Extension Suite
@Suite("Waveform1D Floating Point Operations")
struct Waveform1DFloatingPointTests {

    @Test("Mean calculation for floating point")
    func meanFloatingPoint() {
        let waveform = Waveform1D(values: [1.0, 2.0, 3.0, 4.0, 5.0])
        let expectedMean = 3.0

        #expect(waveform.mean == expectedMean)
    }

    @Test("Mean with empty array")
    func meanEmpty() {
        let waveform = Waveform1D(values: [Double]())

        #expect(waveform.mean == 0.0)
    }

    @Test("RMS calculation")
    func rmsCalculation() {
        let waveform = Waveform1D(values: [3.0, 4.0])  // 3-4-5 triangle
        let expectedRMS = sqrt((9.0 + 16.0) / 2.0)  // sqrt(25/2) = sqrt(12.5)

        #expect(abs(waveform.rms - expectedRMS) < 1e-10)
    }

    @Test("RMS with empty array")
    func rmsEmpty() {
        let waveform = Waveform1D(values: [Float]())

        #expect(waveform.rms == 0.0)
    }

    @Test("Standard deviation calculation")
    func standardDeviation() {
        let waveform = Waveform1D(values: [2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0])
        // Expected std dev ≈ 2.138 (using sample standard deviation)

        #expect(waveform.standardDeviation > 2.0)
        #expect(waveform.standardDeviation < 3.0)
    }

    @Test("Standard deviation with single value")
    func standardDeviationSingle() {
        let waveform = Waveform1D(values: [5.0])

        #expect(waveform.standardDeviation == 0.0)
    }

    @Test("Variance calculation")
    func variance() {
        let waveform = Waveform1D(values: [1.0, 2.0, 3.0])
        let expectedVariance = ((1.0 - 2.0) * (1.0 - 2.0) + (2.0 - 2.0) * (2.0 - 2.0) + (3.0 - 2.0) * (3.0 - 2.0)) / 2.0

        #expect(abs(waveform.variance - expectedVariance) < 1e-10)
    }

    @Test("Variance with single value")
    func varianceSingle() {
        let waveform = Waveform1D(values: [5.0])

        #expect(waveform.variance == 0.0)
    }
}

// MARK: - BinaryInteger Extension Suite
@Suite("Waveform1D Integer Operations")
struct Waveform1DIntegerTests {

    @Test("Mean calculation for integers")
    func meanInteger() {
        let waveform = Waveform1D(values: [1, 2, 3, 4, 5])
        let expectedMean = 15 / 5  // Integer division

        #expect(waveform.mean == expectedMean)
    }

    @Test("Sum calculation for integers")
    func sumInteger() {
        let waveform = Waveform1D(values: [1, 2, 3, 4, 5])
        let expectedSum = 15

        #expect(waveform.sum == expectedSum)
    }

    @Test("Sum with empty array")
    func sumEmpty() {
        let waveform = Waveform1D(values: [Int]())

        #expect(waveform.sum == 0)
    }
}

// MARK: - SignedInteger Extension Suite
@Suite("Waveform1D Signed Integer Operations")
struct Waveform1DSignedIntegerTests {

    @Test("Absolute sum calculation")
    func absoluteSum() {
        let waveform = Waveform1D(values: [-2, 3, -4, 5, -1])
        let expectedAbsoluteSum = 2 + 3 + 4 + 5 + 1

        #expect(waveform.absoluteSum == expectedAbsoluteSum)
    }

    @Test("Absolute sum with all positive values")
    func absoluteSumPositive() {
        let waveform = Waveform1D(values: [1, 2, 3, 4, 5])
        let expectedAbsoluteSum = 15

        #expect(waveform.absoluteSum == expectedAbsoluteSum)
    }

    @Test("Absolute sum with empty array")
    func absoluteSumEmpty() {
        let waveform = Waveform1D(values: [Int]())

        #expect(waveform.absoluteSum == 0)
    }
}

// MARK: - Edge Cases and Integration Suite
@Suite("Waveform1D Edge Cases")
struct Waveform1DEdgeCasesTests {

    @Test("Large dataset performance", .timeLimit(.minutes(1)))
    func largeDataset() {
        let largeArray = Array(repeating: 1.0, count: 100_000)
        let waveform = Waveform1D(values: largeArray, dt: 0.001)

        #expect(waveform.sampleCount == 100_000)
        #expect(waveform.mean == 1.0)
        #expect(waveform.duration == 99.999)  // (100000-1) * 0.001
    }

    @Test("Floating point precision")
    func floatingPointPrecision() {
        let waveform = Waveform1D(values: [0.1, 0.2, 0.3])
        let expectedMean = 0.2

        #expect(abs(waveform.mean - expectedMean) < 1e-15)
    }

    @Test("Negative values handling")
    func negativeValues() {
        let waveform = Waveform1D(values: [-5.0, -2.0, 3.0, 7.0])

        #expect(waveform.minimum == -5.0)
        #expect(waveform.maximum == 7.0)
        #expect(waveform.peakToPeak == 12.0)
        #expect(waveform.mean == 0.75)  // (-5-2+3+7)/4
    }
}

