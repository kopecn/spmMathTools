import Foundation
import Testing

@testable import spmMathTools

// MARK: - Basic Peak Detection Suite
@Suite("Basic Peak Detection")
struct BasicPeakDetectionTests {

    /// VALIDATED
    @Test("Simple peak detection")
    func simplePeakDetection() {
        let values = [1.0, 3.0, 2.0, 5.0, 1.0, 4.0, 0.5]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let peaks = waveform.detectPeaks()

        // Expected peaks at indices 1 (value 3.0), 3 (value 5.0), and 5 (value 4.0)
        #expect(peaks.count == 3)
        #expect(peaks[0].index == 1)
        #expect(peaks[0].value == 3.0)
        #expect(peaks[1].index == 3)
        #expect(peaks[1].value == 5.0)
        #expect(peaks[2].index == 5)
        #expect(peaks[2].value == 4.0)
    }

    /// VALIDATED
    @Test("Large peak detection")
    func largePeakDetection() {
        let waveform1 = Waveform1D<Double>.sine(
            frequency: 100,
            amplitude: 1.0,
            phase: 0.0,
            duration: 3,
            samplingRate: 10000
        )
        
        let peaks1 = waveform1.detectPeaks()

        #expect(peaks1.count == 300)

        let waveform2 = Waveform1D<Double>.sine(
            frequency: 100,
            amplitude: 1.0,
            phase: 0.0,
            duration: 3,
            samplingRate: 1000
        )
        
        let peaks2 = waveform2.detectPeaks()

        #expect(peaks2.count == 300)

        /// Push hard against the nyquist frequency.
        let waveform3 = Waveform1D<Double>.sine(
            frequency: 100,
            amplitude: 1.0,
            phase: 0.0,
            duration: 3,
            samplingRate: 210
        )
        
        let peaks3 = waveform3.detectPeaks()

        #expect(peaks3.count == 300)
    }

    @Test("Peak detection with threshold")
    func peakDetectionWithThreshold() {
        let values = [1.0, 3.0, 2.0, 5.0, 1.0, 2.5, 0.5]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let peaks = waveform.detectPeaks(threshold: 4.0)

        // Only peak at index 3 (value 5.0) should be detected
        #expect(peaks.count == 1)
        #expect(peaks[0].index == 3)
        #expect(peaks[0].value == 5.0)
    }

    @Test("Peak detection with minimum distance")
    func peakDetectionWithMinDistance() {
        let values = [1.0, 3.0, 2.8, 5.0, 1.0, 4.0, 0.5, 1.0, 2.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let peaks = waveform.detectPeaks(minDistance: 3)

        // Should detect peaks with at least 3 samples between them
        // Peaks at indices: 1 (3.0), 3 (5.0), 5 (4.0), 8 (2.0)
        // With minDistance=3: picks highest (5.0 at idx 3), then next valid is idx 8 (distance=5)
        #expect(peaks.count >= 1)
        #expect(peaks[0].index == 3)  // Highest peak
        if peaks.count > 1 {
            #expect(peaks[1].index == 8)  // Second peak at least 3 samples away
        }
    }

    @Test("Peak detection with prominence")
    func peakDetectionWithProminence() {
        let values = [0.0, 1.0, 0.5, 4.0, 3.0, 5.0, 2.0, 1.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let peaks = waveform.detectPeaks(prominence: 2.0)

        // Only peaks with prominence >= 2.0 should be detected
        #expect(!peaks.isEmpty)
        #expect(peaks.contains { $0.index == 5 })  // Peak at index 5 should have high prominence
    }

    @Test("Peak detection with all parameters")
    func peakDetectionWithAllParameters() {
        let values = [0.0, 2.0, 1.0, 5.0, 2.0, 1.0, 3.0, 0.5, 4.0, 1.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let peaks = waveform.detectPeaks(
            threshold: 3.0,
            prominence: 1.0,
            minDistance: 2,
            edgePeaks: false
        )

        #expect(!peaks.isEmpty)
        for peak in peaks {
            #expect(peak.value >= 3.0)
        }
    }

    @Test("Peak detection consistency")
    func peakDetectionConsistency() {
        let values = [1.0, 4.0, 2.0, 6.0, 3.0, 5.0, 1.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        // Run peak detection multiple times
        let peaks1 = waveform.detectPeaks()
        let peaks2 = waveform.detectPeaks()

        #expect(peaks1.count == peaks2.count)

        for (peak1, peak2) in zip(peaks1, peaks2) {
            #expect(peak1.index == peak2.index)
            #expect(peak1.value == peak2.value)
            #expect(abs(peak1.timeOffset - peak2.timeOffset) < 1e-10)
        }
    }
}

// MARK: - Edge Peak Detection Suite
@Suite("Edge Peak Detection")
struct EdgePeakDetectionTests {

    @Test("Peak detection with edge peaks enabled")
    func peakDetectionWithEdgePeaks() {
        let values = [5.0, 3.0, 2.0, 4.0, 1.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        // Without edge peaks
        let peaksNoEdge = waveform.detectPeaks(edgePeaks: false)
        #expect(peaksNoEdge.count == 1)  // Only index 3

        // With edge peaks
        let peaksWithEdge = waveform.detectPeaks(edgePeaks: true)
        #expect(peaksWithEdge.count == 2)  // Index 0 and 3
        #expect(peaksWithEdge.contains { $0.index == 0 })
        #expect(peaksWithEdge.contains { $0.index == 3 })
    }

    @Test("Edge peaks behavior")
    func edgePeaksDetection() {
        let values = [5.0, 3.0, 4.0, 2.0, 1.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        // Test with edge peaks enabled
        let peaksWithEdge = waveform.detectPeaks(edgePeaks: true)
        let peaksWithoutEdge = waveform.detectPeaks(edgePeaks: false)

        // Should have more (or equal) peaks when edge peaks are enabled
        #expect(peaksWithEdge.count >= peaksWithoutEdge.count)

        // First element (5.0) should be detected as peak when edge peaks are enabled
        #expect(peaksWithEdge.contains { $0.index == 0 })
    }

    @Test("Peak detection at boundaries")
    func peakDetectionAtBoundaries() {
        let values = [10.0, 1.0, 2.0, 1.0, 15.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let peaksWithEdge = waveform.detectPeaks(edgePeaks: true)
        let peaksWithoutEdge = waveform.detectPeaks(edgePeaks: false)

        // With edge peaks, should detect first and last elements
        #expect(peaksWithEdge.contains { $0.index == 0 && $0.value == 10.0 })
        #expect(peaksWithEdge.contains { $0.index == 4 && $0.value == 15.0 })

        // Without edge peaks, should only detect internal peak at index 2
        #expect(peaksWithoutEdge.contains { $0.index == 2 && $0.value == 2.0 })
        #expect(!peaksWithoutEdge.contains { $0.index == 0 })
        #expect(!peaksWithoutEdge.contains { $0.index == 4 })
    }
}

// MARK: - Valley Detection Suite
@Suite("Valley Detection")
struct ValleyDetectionTests {

    @Test("Simple valley detection")
    func simpleValleyDetection() {
        let values = [3.0, 1.0, 4.0, 0.5, 3.0, 2.0, 5.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let valleys = waveform.detectValleys()

        // Expected valleys at indices 1 (value 1.0) and 3 (value 0.5)
        #expect(valleys.count >= 2)
        #expect(valleys.contains { $0.index == 1 && $0.value == 1.0 })
        #expect(valleys.contains { $0.index == 3 && $0.value == 0.5 })
    }

    /// VALIDATED
    @Test("Large valley detection")
    func largeValleyDetection() {
        let waveform1 = Waveform1D<Double>.sine(
            frequency: 100,
            amplitude: 1.0,
            phase: 0.0,
            duration: 3,
            samplingRate: 10000
        )
        
        let valleys1 = waveform1.detectValleys()

        #expect(valleys1.count == 300)

        let waveform2 = Waveform1D<Double>.sine(
            frequency: 100,
            amplitude: 1.0,
            phase: 0.0,
            duration: 3,
            samplingRate: 1000
        )
        
        let valleys2 = waveform2.detectValleys()

        #expect(valleys2.count == 300)

        /// Push hard against the nyquist frequency.
        let waveform3 = Waveform1D<Double>.sine(
            frequency: 100,
            amplitude: 1.0,
            phase: 0.0,
            duration: 3,
            samplingRate: 250
        )
        
        let valleys3 = waveform3.detectValleys()

        #expect(valleys3.count == 300)
    }

    @Test("Valley detection with threshold")
    func valleyDetectionWithThreshold() {
        let values = [3.0, 1.0, 4.0, 0.5, 3.0, 2.0, 5.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let valleys = waveform.detectValleys(threshold: 0.8)

        // Only valley at index 3 (value 0.5) should be detected
        #expect(valleys.count == 1)
        #expect(valleys[0].index == 3)
        #expect(valleys[0].value == 0.5)
    }

    @Test("Valley detection with prominence")
    func valleyDetectionWithProminence() {
        let values = [5.0, 2.0, 4.0, 0.5, 3.0, 1.0, 4.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let valleys = waveform.detectValleys(prominence: 1.0)

        #expect(!valleys.isEmpty)
        // Valley at index 3 (value 0.5) should have high prominence
        #expect(valleys.contains { $0.index == 3 })
    }

    @Test("Valley detection with minimum distance")
    func valleyDetectionWithMinDistance() {
        let values = [3.0, 1.0, 2.0, 0.8, 3.0, 0.5, 4.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let valleys = waveform.detectValleys(minDistance: 3)

        // Should apply minimum distance filtering
        #expect(!valleys.isEmpty)

        // Check that selected valleys are at least 3 samples apart
        if valleys.count > 1 {
            for i in 1..<valleys.count {
                let distance = abs(valleys[i].index - valleys[i - 1].index)
                #expect(distance >= 3)
            }
        }
    }
}

// MARK: - Peak Prominence Suite
@Suite("Peak Prominence")
struct PeakProminenceTests {

    @Test("Find most prominent peaks")
    func findMostProminentPeaks() {
        let values = [1.0, 6.0, 2.0, 5.0, 1.0, 4.0, 0.5, 7.0, 3.0]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let topPeaks = waveform.findMostProminentPeaks(count: 3, minDistance: 1)

        #expect(topPeaks.count <= 3)  // May be fewer if there aren't enough peaks
        #expect(!topPeaks.isEmpty)

        // Should include the highest peaks by prominence
        #expect(topPeaks.contains { $0.value == 7.0 || $0.value == 6.0 })
    }
}

// MARK: - Time Information Suite
@Suite("Peak Time Information")
struct PeakTimeInformationTests {

    @Test("Peak time offset calculation")
    func peakTimeOffsetCalculation() {
        let values = [1.0, 3.0, 2.0, 5.0, 1.0]
        let dt = 0.2
        let waveform = Waveform1D(values: values, dt: dt)

        let peaks = waveform.detectPeaks()

        for peak in peaks {
            let expectedTimeOffset = TimeInterval(peak.index) * dt
            #expect(abs(peak.timeOffset - expectedTimeOffset) < 1e-10)
        }
    }

    @Test("Peak analysis with absolute timestamps")
    func peakAnalysisWithTimestamps() {
        let startTime = Date()
        let values = [1.0, 3.0, 2.0, 5.0, 1.0]
        let dt = 0.5
        let waveform = Waveform1D(values: values, dt: dt, t0: startTime)

        let peaks = waveform.detectPeaks()

        for peak in peaks {
            // Each peak should have valid time information
            #expect(peak.timeOffset >= 0)
            #expect(peak.timeOffset < waveform.duration)

            // Check that timeOffset matches expected calculation
            let expectedTimeOffset = TimeInterval(peak.index) * dt
            #expect(abs(peak.timeOffset - expectedTimeOffset) < 1e-10)

            // Check absolute time if available
            if let peakTime = peak.time {
                let expectedTime = startTime.addingTimeInterval(expectedTimeOffset)
                #expect(abs(peakTime.timeIntervalSince(expectedTime)) < 1e-10)
            }
        }
    }
}

// MARK: - Edge Cases Suite
@Suite("Peak Detection Edge Cases")
struct PeakDetectionEdgeCasesTests {

    @Test("Peak detection with empty array")
    func peakDetectionEmpty() {
        let waveform = Waveform1D(values: [Double](), dt: 0.1)

        let peaks = waveform.detectPeaks()

        #expect(peaks.isEmpty)
    }

    @Test("Peak detection with single value")
    func peakDetectionSingle() {
        let waveform = Waveform1D(values: [5.0], dt: 0.1)

        let peaks = waveform.detectPeaks()

        #expect(peaks.isEmpty)  // Single value cannot be a peak
    }

    @Test("Peak detection with two values")
    func peakDetectionTwoValues() {
        let waveform = Waveform1D(values: [3.0, 1.0], dt: 0.1)

        let peaks = waveform.detectPeaks()

        #expect(peaks.isEmpty)  // Need at least 3 values for peak detection
    }

    @Test("Peak detection with constant values")
    func peakDetectionConstant() {
        let waveform = Waveform1D(values: [5.0, 5.0, 5.0, 5.0, 5.0], dt: 0.1)

        let peaks = waveform.detectPeaks()

        #expect(peaks.isEmpty)  // No peaks in constant signal
    }

    @Test("Peak detection with all same non-zero values")
    func peakDetectionAllSame() {
        let waveform = Waveform1D(values: Array(repeating: 3.0, count: 10), dt: 0.1)

        let peaks = waveform.detectPeaks()

        #expect(peaks.isEmpty)
    }
}

// MARK: - Integer Type Suite
@Suite("Integer Type Peak Detection")
struct IntegerPeakDetectionTests {

    @Test("Peak detection with integer values")
    func peakDetectionInteger() {
        let values = [1, 5, 3, 8, 2, 6, 1]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let peaks = waveform.detectPeaks()

        // Should detect peaks at indices 1 (5), 3 (8), and 5 (6)
        #expect(peaks.count == 3)
        #expect(peaks.contains { $0.index == 1 && $0.value == 5 })
        #expect(peaks.contains { $0.index == 3 && $0.value == 8 })
        #expect(peaks.contains { $0.index == 5 && $0.value == 6 })
    }

    @Test("Peak detection with negative integer values")
    func peakDetectionNegativeIntegers() {
        let values = [-5, -2, -4, -1, -3, -6]
        let waveform = Waveform1D(values: values, dt: 0.1)

        let peaks = waveform.detectPeaks()

        // Peaks are local maxima, so -2 at index 1 and -1 at index 3 should be peaks
        #expect(peaks.contains { $0.index == 1 && $0.value == -2 })
        #expect(peaks.contains { $0.index == 3 && $0.value == -1 })
    }
}

// MARK: - Performance Tests Suite
@Suite("Performance Tests")
struct PeakDetectionPerformanceTests {

    @Test("Peak detection performance with large dataset", .timeLimit(.minutes(1)))
    func peakDetectionPerformance() {
        // Create a large signal with known peaks
        var values: [Double] = []
        for i in 0..<10000 {
            let x = Double(i) * 0.01
            values.append(sin(x) + 0.5 * sin(3 * x) + 0.1 * Double.random(in: -1...1))
        }

        let waveform = Waveform1D(values: values, dt: 0.01)

        let peaks = waveform.detectPeaks(threshold: 0.5, minDistance: 10)

        // Should complete in reasonable time and find some peaks
        #expect(!peaks.isEmpty)
        #expect(peaks.count < values.count)  // Sanity check
    }
}
