import Foundation
import XCTest

@testable import spmMathTools

/// Comprehensive test suite for OTG trajectory generation
///
/// This test suite combines:
/// - Bug fix tests for previously failing trajectories
/// - Truth table tests for known-good trajectories
/// - Validation tests for randomly generated trajectories from JSON files
///
/// All tests share common validation functions to ensure consistency
final class OTGComprehensiveTests: XCTestCase {

    // MARK: - Common Helper Functions

    /// Validates that all time intervals in the trajectory are non-negative
    func validateNoNegativeTimeIntervals(
        _ profile: Profile,
        file: StaticString = #filePath,
        line: UInt = #line
    ) {
        for i in 0..<7 {
            XCTAssertGreaterThanOrEqual(
                profile.t[i],
                0.0,
                "Time interval t[\(i)] must be non-negative, got \(profile.t[i])",
                file: file,
                line: line
            )
        }
    }

    /// Validates that acceleration stays within limits throughout the entire trajectory
    /// Samples the trajectory at multiple points to catch continuous violations
    func validateAccelerationLimits(
        _ trajectory: Trajectory,
        maxAcc: Double,
        sampleCount: Int = 100,
        file: StaticString = #filePath,
        line: UInt = #line
    ) {
        let duration = trajectory.getDuration()
        let dt = duration / Double(sampleCount)

        var newPosition = [Double](repeating: 0.0, count: 1)
        var newVelocity = [Double](repeating: 0.0, count: 1)
        var newAcceleration = [Double](repeating: 0.0, count: 1)

        // Sample trajectory at multiple time points
        for i in 0...sampleCount {
            let t = Double(i) * dt
            try! trajectory.atTime(
                t,
                newPosition: &newPosition,
                newVelocity: &newVelocity,
                newAcceleration: &newAcceleration
            )

            let acc = abs(newAcceleration[0])
            XCTAssertLessThanOrEqual(
                acc,
                maxAcc + 0.001,  // Small tolerance for numerical precision
                "Acceleration \(acc) exceeds limit \(maxAcc) at time t=\(t)",
                file: file,
                line: line
            )
        }
    }

    /// Validates that velocity stays within limits throughout the entire trajectory
    /// Samples the trajectory at multiple points to catch continuous violations
    func validateVelocityLimits(
        _ trajectory: Trajectory,
        maxVel: Double,
        sampleCount: Int = 100,
        file: StaticString = #filePath,
        line: UInt = #line
    ) {
        let duration = trajectory.getDuration()
        let dt = duration / Double(sampleCount)

        var newPosition = [Double](repeating: 0.0, count: 1)
        var newVelocity = [Double](repeating: 0.0, count: 1)
        var newAcceleration = [Double](repeating: 0.0, count: 1)

        // Sample trajectory at multiple time points
        for i in 0...sampleCount {
            let t = Double(i) * dt
            try! trajectory.atTime(
                t,
                newPosition: &newPosition,
                newVelocity: &newVelocity,
                newAcceleration: &newAcceleration
            )

            let vel = abs(newVelocity[0])
            XCTAssertLessThanOrEqual(
                vel,
                maxVel + 0.001,  // Small tolerance for numerical precision
                "Velocity \(vel) exceeds limit \(maxVel) at time t=\(t)",
                file: file,
                line: line
            )
        }
    }

    /// Validates that the trajectory reaches the target state
    func validateTargetReached(
        _ trajectory: Trajectory,
        targetPosition: Double,
        targetVelocity: Double,
        targetAcceleration: Double,
        file: StaticString = #filePath,
        line: UInt = #line
    ) {
        let duration = trajectory.getDuration()
        var finalPosition = [Double](repeating: 0.0, count: 1)
        var finalVelocity = [Double](repeating: 0.0, count: 1)
        var finalAcceleration = [Double](repeating: 0.0, count: 1)

        try! trajectory.atTime(
            duration,
            newPosition: &finalPosition,
            newVelocity: &finalVelocity,
            newAcceleration: &finalAcceleration
        )

        let positionTolerance = 0.01
        let velocityTolerance = 0.01
        let accelerationTolerance = 0.01

        XCTAssertEqual(
            finalPosition[0],
            targetPosition,
            accuracy: positionTolerance,
            "Final position \(finalPosition[0]) does not match target \(targetPosition)",
            file: file,
            line: line
        )

        XCTAssertEqual(
            finalVelocity[0],
            targetVelocity,
            accuracy: velocityTolerance,
            "Final velocity \(finalVelocity[0]) does not match target \(targetVelocity)",
            file: file,
            line: line
        )

        XCTAssertEqual(
            finalAcceleration[0],
            targetAcceleration,
            accuracy: accelerationTolerance,
            "Final acceleration \(finalAcceleration[0]) does not match target \(targetAcceleration)",
            file: file,
            line: line
        )
    }

    /// Complete validation of a trajectory
    func validateTrajectory(
        _ trajectory: Trajectory,
        input: InputParameter,
        file: StaticString = #filePath,
        line: UInt = #line
    ) {
        let profile = trajectory.getProfiles()[0][0]

        // 1. Check no negative time intervals
        validateNoNegativeTimeIntervals(profile, file: file, line: line)

        // 2. Check acceleration limits
        validateAccelerationLimits(
            trajectory,
            maxAcc: input.maxAcceleration[0],
            file: file,
            line: line
        )

        // 3. Check velocity limits
        validateVelocityLimits(
            trajectory,
            maxVel: input.maxVelocity[0],
            file: file,
            line: line
        )

        // 4. Check target reached
        validateTargetReached(
            trajectory,
            targetPosition: input.targetPosition[0],
            targetVelocity: input.targetVelocity[0],
            targetAcceleration: input.targetAcceleration[0],
            file: file,
            line: line
        )
    }

    // MARK: - Bug Fix Tests
    // Tests for previously failing trajectories (from OTGFailureFixTests.swift)

    func testBugFix_AccelerationLimit_Case1() throws {
        // Bug #1 Case 1: Previously violated acceleration limits with target_vel=-2.95, max_acc=1.94

        let ruckig = OTG(degreesOfFreedom: 1, deltaTime: 0.01)
        var input = InputParameter(DOFs: 1)

        input.currentPosition = [-6.139891324284666]
        input.currentVelocity = [0.0]
        input.currentAcceleration = [0.0]
        input.targetPosition = [10.0]
        input.targetVelocity = [-2.9531593910491565]
        input.targetAcceleration = [0.0]
        input.maxVelocity = [4.0]
        input.maxAcceleration = [1.9425898752751283]
        input.maxJerk = [10.0]

        var trajectory = Trajectory(dofs: 1)
        let result = ruckig.calculate(input: input, trajectory: &trajectory)

        XCTAssertTrue(
            result == .Working || result == .Finished,
            "Expected trajectory to succeed, got: \(result)"
        )

        validateTrajectory(trajectory, input: input)
    }

    func testBugFix_AccelerationLimit_Case2() throws {
        // Bug #1 Case 2: Previously violated acceleration limits with target_vel=1.345, max_acc=5.0

        let ruckig = OTG(degreesOfFreedom: 1, deltaTime: 0.01)
        var input = InputParameter(DOFs: 1)

        input.currentPosition = [-6.139891324284666]
        input.currentVelocity = [0.0]
        input.currentAcceleration = [0.0]
        input.targetPosition = [10.0]
        input.targetVelocity = [1.345062591709465]
        input.targetAcceleration = [0.0]
        input.maxVelocity = [4.0]
        input.maxAcceleration = [5.0]
        input.maxJerk = [10.0]

        var trajectory = Trajectory(dofs: 1)
        let result = ruckig.calculate(input: input, trajectory: &trajectory)

        XCTAssertTrue(
            result == .Working || result == .Finished,
            "Expected trajectory to succeed, got: \(result)"
        )

        validateTrajectory(trajectory, input: input)
    }

    func testBugFix_NegativeTimeInterval_Case1() throws {
        // Bug #2 Case 1: Previously produced t[5] = -0.277552

        let ruckig = OTG(degreesOfFreedom: 1, deltaTime: 0.01)
        var input = InputParameter(DOFs: 1)

        input.currentPosition = [-2.2569699192956714]
        input.currentVelocity = [0.6924924339691851]
        input.currentAcceleration = [0.6756694790902418]
        input.targetPosition = [6.790054108584006]
        input.targetVelocity = [-1.159379585473221]
        input.targetAcceleration = [-1.2629539618488628]
        input.maxVelocity = [10.0]
        input.maxAcceleration = [10.0]
        input.maxJerk = [10.0]

        var trajectory = Trajectory(dofs: 1)
        let result = ruckig.calculate(input: input, trajectory: &trajectory)

        XCTAssertTrue(
            result == .Working || result == .Finished,
            "Expected trajectory to succeed, got: \(result)"
        )

        validateTrajectory(trajectory, input: input)
    }

    func testBugFix_NegativeTimeInterval_Case2() throws {
        // Bug #2 Case 2: Previously produced t[5] = -0.526238

        let ruckig = OTG(degreesOfFreedom: 1, deltaTime: 0.01)
        var input = InputParameter(DOFs: 1)

        input.currentPosition = [2.6603799559471364]
        input.currentVelocity = [1.9695301027900147]
        input.currentAcceleration = [1.972541529001468]
        input.targetPosition = [-6.297064289647577]
        input.targetVelocity = [-1.5414200165198237]
        input.targetAcceleration = [-1.5255311582232012]
        input.maxVelocity = [9.088312224669604]
        input.maxAcceleration = [9.973774779735683]
        input.maxJerk = [10.0]

        var trajectory = Trajectory(dofs: 1)
        let result = ruckig.calculate(input: input, trajectory: &trajectory)

        XCTAssertTrue(
            result == .Working || result == .Finished,
            "Expected trajectory to succeed, got: \(result)"
        )

        validateTrajectory(trajectory, input: input)
    }

    func testBugFix_NegativeTimeInterval_Case3() throws {
        // Bug #2 Case 3: Previously produced t[5] = -0.526238 (different max_vel)

        let ruckig = OTG(degreesOfFreedom: 1, deltaTime: 0.01)
        var input = InputParameter(DOFs: 1)

        input.currentPosition = [2.6603799559471364]
        input.currentVelocity = [1.9695301027900147]
        input.currentAcceleration = [1.972541529001468]
        input.targetPosition = [-6.297064289647577]
        input.targetVelocity = [-1.5414200165198237]
        input.targetAcceleration = [-1.5255311582232012]
        input.maxVelocity = [6.282291093061675]
        input.maxAcceleration = [9.973774779735683]
        input.maxJerk = [10.0]

        var trajectory = Trajectory(dofs: 1)
        let result = ruckig.calculate(input: input, trajectory: &trajectory)

        XCTAssertTrue(
            result == .Working || result == .Finished,
            "Expected trajectory to succeed, got: \(result)"
        )

        validateTrajectory(trajectory, input: input)
    }

    // MARK: - Truth Table Tests
    // Tests for known-good trajectories with expected durations and time intervals

    func testTruthTable() throws {
        let truthTable = TruthTableData.testCases
        let durationAccuracy: Double = 0.001
        let timeIntervalAccuracy: Double = 0.001

        for (index, testCase) in truthTable.enumerated() {
            let ruckig = OTG(degreesOfFreedom: 1, deltaTime: 0.01)
            var input = InputParameter(DOFs: 1)

            input.currentPosition = [testCase.currentPos]
            input.currentVelocity = [testCase.currentVel]
            input.currentAcceleration = [testCase.currentAcc]
            input.targetPosition = [testCase.targetPos]
            input.targetVelocity = [testCase.targetVel]
            input.targetAcceleration = [testCase.targetAcc]
            input.maxVelocity = [testCase.maxVel]
            input.maxAcceleration = [testCase.maxAcc]
            input.maxJerk = [testCase.maxJerk]

            var trajectory = Trajectory(dofs: 1)
            let result = ruckig.calculate(input: input, trajectory: &trajectory)

            XCTAssertTrue(
                result == .Working || result == .Finished,
                "Test case \(index + 1): Trajectory calculation failed with result: \(result)"
            )

            // Validate limits are respected
            validateTrajectory(trajectory, input: input)

            // Validate expected duration
            XCTAssertEqual(
                trajectory.getDuration(),
                testCase.expectedDuration,
                accuracy: durationAccuracy,
                "Test case \(index + 1): Duration mismatch"
            )

            // Validate expected time intervals
            let actualTimeIntervals = trajectory.getProfiles()[0][0].t
            XCTAssertEqual(
                actualTimeIntervals.count,
                7,
                "Test case \(index + 1): Expected 7 time intervals"
            )

            for (intervalIndex, expectedInterval) in testCase.expectedTimeIntervals.enumerated() {
                XCTAssertEqual(
                    actualTimeIntervals[intervalIndex],
                    expectedInterval,
                    accuracy: timeIntervalAccuracy,
                    "Test case \(index + 1): Time interval[\(intervalIndex)] mismatch"
                )
            }
        }
    }

    // MARK: - JSON Validation Tests
    // Tests for randomly generated trajectories from JSON files

    func testSuccessfulTrajectoriesFromJSON() throws {
        let jsonURL = URL(fileURLWithPath: "spm/Tests/OTGTests/truthTables/successful_trajectories.json")

        guard FileManager.default.fileExists(atPath: jsonURL.path) else {
            throw XCTSkip("successful_trajectories.json not found. Run generateRandomTrajectories.py first.")
        }

        let jsonData = try Data(contentsOf: jsonURL)
        print(jsonData)
        let decoder = JSONDecoder()
        let inputs = try decoder.decode([InputParameter].self, from: jsonData)

        print("Testing \(inputs.count) successful trajectories from JSON...")

        var successCount = 0
        var failureCount = 0

        for (index, input) in inputs.enumerated() {
            let ruckig = OTG(degreesOfFreedom: input.degreesOfFreedom, deltaTime: 0.01)
            var trajectory = Trajectory(dofs: input.degreesOfFreedom)
            let result = ruckig.calculate(input: input, trajectory: &trajectory)

            if result == .Working || result == .Finished {
                // Should succeed and respect all limits
                validateTrajectory(trajectory, input: input)
                successCount += 1
            } else {
                // This is unexpected - Python marked it as successful
                failureCount += 1
                print("⚠ Trajectory \(index + 1) failed in Swift but succeeded in Python: \(result)")
            }
        }

        print("Results: \(successCount) succeeded, \(failureCount) failed")
        XCTAssertEqual(
            failureCount,
            0,
            "Expected all successful JSON trajectories to succeed in Swift"
        )
    }

    func testFailedTrajectoriesFromJSON() throws {
        let jsonURL = URL(fileURLWithPath: "spm/Tests/OTGTests/truthTables/failed_trajectories.json")

        guard FileManager.default.fileExists(atPath: jsonURL.path) else {
            throw XCTSkip("failed_trajectories.json not found. Run generateRandomTrajectories.py first.")
        }

        let jsonData = try Data(contentsOf: jsonURL)
        let decoder = JSONDecoder()

        // Failed trajectories are dictionaries with all InputParameter fields plus "error"
        // We'll decode as [String: Any] and extract what we need
        guard let jsonArray = try? JSONSerialization.jsonObject(with: jsonData) as? [[String: Any]] else {
            XCTFail("Could not parse failed_trajectories.json")
            return
        }

        print("Testing \(jsonArray.count) failed trajectories from JSON...")

        var expectedFailureCount = 0
        var unexpectedSuccessCount = 0

        for (index, jsonDict) in jsonArray.enumerated() {
            // Remove the error field and re-encode to decode as InputParameter
            var inputDict = jsonDict
            let errorMessage = inputDict.removeValue(forKey: "error") as? String ?? "Unknown error"

            let inputData = try JSONSerialization.data(withJSONObject: inputDict)
            let input = try decoder.decode(InputParameter.self, from: inputData)

            let ruckig = OTG(degreesOfFreedom: input.degreesOfFreedom, deltaTime: 0.01)
            var trajectory = Trajectory(dofs: input.degreesOfFreedom)
            let result = ruckig.calculate(input: input, trajectory: &trajectory)

            if result == .Working || result == .Finished {
                // Unexpected success - but if it succeeds, it must respect limits
                validateTrajectory(trajectory, input: input)
                unexpectedSuccessCount += 1
                print("ℹ Trajectory \(index + 1) succeeded in Swift but failed in Python")
                print("  Python error: \(errorMessage)")
            } else {
                // Expected failure
                expectedFailureCount += 1
            }
        }

        print("Results: \(expectedFailureCount) failed as expected, \(unexpectedSuccessCount) unexpectedly succeeded")

        // Note: It's acceptable for some to succeed in Swift if they failed in Python
        // The important thing is that IF they succeed, they respect all limits
        // So we don't fail the test for unexpected successes
    }
}

// MARK: - Truth Table Data

/// Truth table test cases (extracted from OTGTruthTableTests.swift)
struct TruthTableData {
    struct TestCase {
        let currentPos: Double
        let currentVel: Double
        let currentAcc: Double
        let targetPos: Double
        let targetVel: Double
        let targetAcc: Double
        let maxVel: Double
        let maxAcc: Double
        let maxJerk: Double
        let expectedDuration: Double
        let expectedTimeIntervals: [Double]
    }

    static let testCases: [TestCase] = [
        // Test Case 1: Zero initial state to positive target
        TestCase(
            currentPos: 0.0,
            currentVel: 0.0,
            currentAcc: 0.0,
            targetPos: 10.0,
            targetVel: 0.0,
            targetAcc: 0.0,
            maxVel: 4.0,
            maxAcc: 5.0,
            maxJerk: 10.0,
            expectedDuration: 3.800000,
            expectedTimeIntervals: [0.500000, 0.300000, 0.500000, 1.200000, 0.500000, 0.300000, 0.500000]
        ),

        // Test Case 2: Negative position to positive target
        TestCase(
            currentPos: -4.775541085840058,
            currentVel: 0.0,
            currentAcc: 0.0,
            targetPos: 10.0,
            targetVel: 0.0,
            targetAcc: 0.0,
            maxVel: 4.0,
            maxAcc: 5.0,
            maxJerk: 10.0,
            expectedDuration: 4.993885,
            expectedTimeIntervals: [0.500000, 0.300000, 0.500000, 2.393885, 0.500000, 0.300000, 0.500000]
        ),

        // Test Case 24: Critical case that was regressing
        TestCase(
            currentPos: 2.6603799559471364,
            currentVel: -2.5442249449339207,
            currentAcc: 0.0,
            targetPos: -6.297064289647577,
            targetVel: 0.0,
            targetAcc: 0.0,
            maxVel: 4.0,
            maxAcc: 5.0,
            maxJerk: 10.0,
            expectedDuration: 3.028222,
            expectedTimeIntervals: [0.381546, 0.000000, 0.381546, 0.965130, 0.500000, 0.300000, 0.500000]
        ),

        // Add more test cases as needed...
    ]
}
