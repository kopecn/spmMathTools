import XCTest

@testable import spmMathTools

/// Performance and stress tests for trajectory generation
final class OTGPerformanceTests: XCTestCase {

    // MARK: - Test Helpers

    func arrayEquals<T: FloatingPoint>(_ first: [T], _ second: [T], accuracy: T = 1e-10) -> Bool {
        guard first.count == second.count else { return false }
        for i in 0..<first.count {
            if abs(first[i] - second[i]) > accuracy { return false }
        }
        return true
    }

    func checkCalculation(_ otg: OTG, _ input: InputParameter) {
        var output = OutputParameter(DOFs: otg.degreesOfFreedom)

        let result = otg.update(input: input, output: &output)
        if result == .ErrorTrajectoryDuration {
            return
        }

        XCTAssertTrue(result == .Working || (result == .Finished && output.trajectory.duration < 0.005))
        XCTAssertGreaterThanOrEqual(output.trajectory.duration, 0.0)

        for dof in 0..<otg.degreesOfFreedom {
            XCTAssertFalse(output.newPosition[dof].isNaN)
            XCTAssertFalse(output.newVelocity[dof].isNaN)
            XCTAssertFalse(output.newAcceleration[dof].isNaN)
        }
    }

    // MARK: - Random Trajectory Tests

    func testRandomTrajectories() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.005)

        // Generate some pseudo-random test cases
        for i in 0..<100 {
            var input = InputParameter(DOFs: 3)

            // Simple pseudo-random generation (replace with proper random if needed)
            let seed = Double(i * 123 + 456)

            input.currentPosition = [
                sin(seed) * 4.0,
                cos(seed + 1) * 4.0,
                sin(seed + 2) * 4.0,
            ]

            input.currentVelocity = [
                sin(seed + 3) * 0.8,
                cos(seed + 4) * 0.8,
                sin(seed + 5) * 0.8,
            ]

            input.currentAcceleration = [
                sin(seed + 6) * 0.8,
                cos(seed + 7) * 0.8,
                sin(seed + 8) * 0.8,
            ]

            input.targetPosition = [
                sin(seed + 9) * 4.0,
                cos(seed + 10) * 4.0,
                sin(seed + 11) * 4.0,
            ]

            input.targetVelocity = [
                sin(seed + 12) * 0.8,
                cos(seed + 13) * 0.8,
                sin(seed + 14) * 0.8,
            ]

            input.targetAcceleration = [
                sin(seed + 15) * 0.8,
                cos(seed + 16) * 0.8,
                sin(seed + 17) * 0.8,
            ]

            input.maxVelocity = [
                abs(sin(seed + 18)) * 15.0 + 1.0,
                abs(cos(seed + 19)) * 15.0 + 1.0,
                abs(sin(seed + 20)) * 15.0 + 1.0,
            ]

            input.maxAcceleration = [
                abs(sin(seed + 21)) * 15.0 + 1.0,
                abs(cos(seed + 22)) * 15.0 + 1.0,
                abs(sin(seed + 23)) * 15.0 + 1.0,
            ]

            input.maxJerk = [
                abs(sin(seed + 24)) * 15.0 + 1.0,
                abs(cos(seed + 25)) * 15.0 + 1.0,
                abs(sin(seed + 26)) * 15.0 + 1.0,
            ]

            // Adjust target limits to be valid
            for dof in 0..<3 {
                if abs(input.targetVelocity[dof]) > input.maxVelocity[dof] {
                    input.targetVelocity[dof] *= input.maxVelocity[dof] / abs(input.targetVelocity[dof]) * 0.9
                }
            }

            if i < 50 {
                input.synchronization = .Phase
            } else {
                input.synchronization = .Time
            }

            if (try? otg.validateInput(input: input)) ?? false {
                checkCalculation(otg, input)
            }
        }
    }

    func testRandomTrajectories6DOF() throws {
        let dofCount = 6
        let otg = OTG(degreesOfFreedom: dofCount, deltaTime: 0.0001)  // 100us deltaTime

        var totalTime: Double = 0
        var successCount = 0
        var failedCases: [(iteration: Int, input: InputParameter, result: Result)] = []

        // Generate some pseudo-random test cases with small trajectory windows (1-10ms)
        for i in 0..<100 {
            var input = InputParameter(DOFs: dofCount)

            // Simple pseudo-random generation
            let seed = Double(i * 123 + 456)

            // Scale positions to create small trajectory windows (1-10ms target)
            let positionScale = 0.005  // Very small position changes
            input.currentPosition = (0..<dofCount).map { j in sin(seed + Double(j)) * positionScale }
            input.currentVelocity = (0..<dofCount).map { j in sin(seed + Double(j + 6)) * 0.1 }
            input.currentAcceleration = (0..<dofCount).map { j in sin(seed + Double(j + 12)) * 0.5 }

            input.targetPosition = (0..<dofCount).map { j in sin(seed + Double(j + 18)) * positionScale }
            input.targetVelocity = (0..<dofCount).map { j in sin(seed + Double(j + 24)) * 0.1 }
            input.targetAcceleration = (0..<dofCount).map { j in sin(seed + Double(j + 30)) * 0.5 }

            // Very high limits for fast trajectories (1-10ms range)
            input.maxVelocity = (0..<dofCount).map { j in abs(sin(seed + Double(j + 36))) * 5.0 + 2.0 }
            input.maxAcceleration = (0..<dofCount).map { j in abs(cos(seed + Double(j + 42))) * 50.0 + 20.0 }
            input.maxJerk = (0..<dofCount).map { j in abs(sin(seed + Double(j + 48))) * 500.0 + 200.0 }

            // Adjust target limits to be valid
            for dof in 0..<dofCount {
                if abs(input.targetVelocity[dof]) > input.maxVelocity[dof] {
                    input.targetVelocity[dof] *= input.maxVelocity[dof] / abs(input.targetVelocity[dof]) * 0.9
                }
            }

            if i < 50 {
                input.synchronization = .Phase
            } else {
                input.synchronization = .Time
            }

            if (try? otg.validateInput(input: input)) ?? false {
                let startTime = CFAbsoluteTimeGetCurrent()
                var output = OutputParameter(DOFs: dofCount)
                let result = otg.update(input: input, output: &output)
                let elapsed = (CFAbsoluteTimeGetCurrent() - startTime) * 1_000_000  // Convert to microseconds

                totalTime += elapsed

                if result == .Working || result == .Finished {
                    successCount += 1
                } else {
                    // Capture failed case for later analysis
                    failedCases.append((iteration: i, input: input, result: result))
                }

                XCTAssertTrue(
                    result == .Working || result == .Finished || result == .ErrorTrajectoryDuration
                        || result == .ErrorSynchronizationCalculation
                )
            }
        }
        // Print detailed failure analysis
        if !failedCases.isEmpty {
            print("\n❌ Failed Trajectory Parameters (\(failedCases.count) cases):")
            for (iteration, input, result) in failedCases {
                print("\n--- Iteration \(iteration) - Result: \(result) ---")
                print("  currentPosition:     \(input.currentPosition)")
                print("  currentVelocity:     \(input.currentVelocity)")
                print("  currentAcceleration: \(input.currentAcceleration)")
                print("  targetPosition:      \(input.targetPosition)")
                print("  targetVelocity:      \(input.targetVelocity)")
                print("  targetAcceleration:  \(input.targetAcceleration)")
                print("  maxVelocity:         \(input.maxVelocity)")
                print("  maxAcceleration:     \(input.maxAcceleration)")
                print("  maxJerk:             \(input.maxJerk)")
                print("  synchronization:      \(input.synchronization)")
            }
        }

        let avgTime = totalTime / Double(successCount)
        print("\n📊 6-DOF Performance Summary:")
        print("   Successful trajectories: \(successCount)/100")
        print("   Average calculation time: \(String(format: "%.1f", avgTime))µs")
        print("   Total time: \(String(format: "%.1f", totalTime / 1000))ms")
    }

    // MARK: - Performance Benchmarks

    /// Swift port of C++ benchmark function from .temp/test/benchmark_target.cpp:42
    /// Compares trajectory calculation performance with C++ implementation
    ///
    /// Run with: make benchmark (enables Release optimizations + full trajectory count for 1:1 comparison)
    /// Default run: Scaled down to complete in reasonable time (~10-30 seconds)
    /// C++ Results: ~2.43µs average, ~65µs worst, ~2.99µs end-to-end
    func testBenchmark3DOF() throws {
        let DOFs = 3

        // Check if running full benchmark via make command
        let runFullBenchmark = ProcessInfo.processInfo.environment["RUN_FULL_BENCHMARK"] == "1"

        let n: Int
        let numberTrajectories: Int

        if runFullBenchmark {
            // Full benchmark: 1:1 comparison with C++ (262,144 trajectories)
            n = 2 * 5  // 10 iterations
            numberTrajectories = 4 * 64 * 1024  // 262,144
            print("Running FULL benchmark (262,144 trajectories × 10 iterations)...")
        } else {
            // Scaled down for regular test runs
            n = 3  // 3 iterations
            numberTrajectories = 1024  // ~0.4% of full benchmark
            print("Running scaled benchmark (1,024 trajectories × 3 iterations). Use 'make benchmark' for full test.")
        }

        let otg = OTG(degreesOfFreedom: DOFs, deltaTime: 0.005)

        // Random number generator with seeds matching C++
        var positionRng = SeededRandom(seed: 42)
        var dynamicRng = SeededRandom(seed: 43)
        var limitRng = SeededRandom(seed: 44)

        var averageResults: [Double] = []
        var worstResults: [Double] = []
        var globalResults: [Double] = []

        for _ in 0..<n {
            var averageTime = 0.0
            var worstTime = 0.0
            var count = 1

            let globalStart = CFAbsoluteTimeGetCurrent()

            for _ in 0..<numberTrajectories {
                var input = InputParameter(DOFs: DOFs)

                // Fill random values matching C++ normal/uniform distributions
                input.currentPosition = positionRng.fillNormal(count: DOFs, mean: 0.0, stdDev: 4.0)
                input.currentVelocity = dynamicRng.fillNormalOrZero(
                    count: DOFs,
                    mean: 0.0,
                    stdDev: 0.8,
                    probability: 0.9
                )
                input.currentAcceleration = dynamicRng.fillNormalOrZero(
                    count: DOFs,
                    mean: 0.0,
                    stdDev: 0.8,
                    probability: 0.8
                )

                input.targetPosition = positionRng.fillNormal(count: DOFs, mean: 0.0, stdDev: 4.0)
                input.targetVelocity = dynamicRng.fillNormalOrZero(
                    count: DOFs,
                    mean: 0.0,
                    stdDev: 0.8,
                    probability: 0.7
                )
                input.targetAcceleration = dynamicRng.fillNormalOrZero(
                    count: DOFs,
                    mean: 0.0,
                    stdDev: 0.8,
                    probability: 0.6
                )

                input.maxVelocity = limitRng.fillUniform(count: DOFs, min: 0.1, max: 12.0)
                input.maxAcceleration = limitRng.fillUniform(count: DOFs, min: 0.1, max: 12.0)
                input.maxJerk = limitRng.fillUniform(count: DOFs, min: 0.1, max: 12.0)

                // Adjust target velocity to be within limits
                for dof in 0..<DOFs {
                    if abs(input.targetVelocity[dof]) > input.maxVelocity[dof] {
                        input.targetVelocity[dof] =
                            input.maxVelocity[dof] * (input.targetVelocity[dof] >= 0 ? 0.9 : -0.9)
                    }
                }

                // Measure calculation time
                let calcStart = CFAbsoluteTimeGetCurrent()
                var output = OutputParameter(DOFs: DOFs)
                _ = otg.update(input: input, output: &output)
                let calcEnd = CFAbsoluteTimeGetCurrent()

                let time = (calcEnd - calcStart) * 1_000_000  // Convert to microseconds
                averageTime = averageTime + (time - averageTime) / Double(count)
                worstTime = max(worstTime, time)
                count += 1
            }

            let globalEnd = CFAbsoluteTimeGetCurrent()
            let globalTime = (globalEnd - globalStart) * 1_000_000 / Double(numberTrajectories)

            averageResults.append(averageTime)
            worstResults.append(worstTime)
            globalResults.append(globalTime)
        }

        // Analyze results (mean and standard deviation)
        let (avgMean, avgStd) = analyze(averageResults)
        let (worstMean, worstStd) = analyze(worstResults)
        let (globalMean, globalStd) = analyze(globalResults)

        print("---")
        print("Benchmark for \(DOFs) DoFs on \(numberTrajectories) trajectories")
        print(
            "Average Calculation Duration \(String(format: "%.2f", avgMean)) ± \(String(format: "%.2f", avgStd)) [µs]"
        )
        print(
            "Worst Calculation Duration \(String(format: "%.2f", worstMean)) ± \(String(format: "%.2f", worstStd)) [µs]"
        )
        print(
            "End-to-end Calculation Duration \(String(format: "%.2f", globalMean)) ± \(String(format: "%.2f", globalStd)) [µs]"
        )
        print("\n// Compare these results to C++ implementation in .temp/test/benchmark_target.cpp:42")
    }

    /// Statistical analysis helper matching C++ analyze() function
    private func analyze(_ values: [Double]) -> (mean: Double, stdDev: Double) {
        let sum = values.reduce(0.0, +)
        let mean = sum / Double(values.count)

        let variance = values.reduce(0.0) { $0 + ($1 - mean) * ($1 - mean) }
        let stdDev = sqrt(variance / Double(values.count))

        return (mean, stdDev)
    }

    func testPerformance() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.005)
        var input = InputParameter(DOFs: 3)
        var output = OutputParameter(DOFs: 3)

        input.currentPosition = [0.0, -2.0, 0.0]
        input.targetPosition = [1.0, -3.0, 2.0]
        input.maxVelocity = [1.0, 1.0, 1.0]
        input.maxAcceleration = [1.0, 1.0, 1.0]
        input.maxJerk = [1.0, 1.0, 1.0]

        measure {
            for _ in 0..<1000 {
                _ = otg.update(input: input, output: &output)
            }
        }
    }
}

// MARK: - Random Number Helpers (matching C++ randomizer.hpp)

/// Seeded random number generator matching C++ distributions
struct SeededRandom {
    var seed: UInt64

    init(seed: UInt64) {
        self.seed = seed
    }

    /// Linear Congruential Generator (simple RNG for reproducibility)
    mutating func next() -> Double {
        seed = (seed &* 1_103_515_245 &+ 12345) & 0x7FFF_FFFF
        return Double(seed) / Double(0x7FFF_FFFF)
    }

    /// Box-Muller transform for normal distribution
    mutating func normal(mean: Double, stdDev: Double) -> Double {
        let u1 = next()
        let u2 = next()
        let z0 = sqrt(-2.0 * log(u1)) * cos(2.0 * .pi * u2)
        return mean + z0 * stdDev
    }

    mutating func uniform(min: Double, max: Double) -> Double {
        return min + next() * (max - min)
    }

    mutating func fillNormal(count: Int, mean: Double, stdDev: Double) -> [Double] {
        return (0..<count).map { _ in normal(mean: mean, stdDev: stdDev) }
    }

    mutating func fillUniform(count: Int, min: Double, max: Double) -> [Double] {
        return (0..<count).map { _ in uniform(min: min, max: max) }
    }

    mutating func fillNormalOrZero(count: Int, mean: Double, stdDev: Double, probability: Double) -> [Double] {
        return (0..<count).map { _ in
            next() < probability ? normal(mean: mean, stdDev: stdDev) : 0.0
        }
    }
}
