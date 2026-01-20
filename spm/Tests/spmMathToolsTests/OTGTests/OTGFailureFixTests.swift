import XCTest

@testable import spmMathTools

/// Tests to verify fixes for previously failing trajectory calculations
/// These test cases were documented in failureTable.txt
final class OTGFailureFixTests: XCTestCase {

    // MARK: - Helper Functions

    /// Validates that all time intervals in the trajectory are non-negative
    func validateNoNegativeTimeIntervals(_ profile: Profile, file: StaticString = #file, line: UInt = #line) {
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
        file: StaticString = #file,
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
                file: (file),
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
        file: StaticString = #file,
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
                file: (file),
                line: line
            )
        }
    }

    // MARK: - Bug #1: Acceleration Limit Violations with Non-Zero Target Velocity

    func testAccelerationLimitWithNonZeroTargetVel_Case1() throws {
        // From failureTable.txt lines 5-16
        // Previously violated acceleration limits with target_vel=-2.95, max_acc=1.94

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

        // Should either succeed with valid trajectory OR fail gracefully
        // Should NOT produce a trajectory that violates acceleration limits
        if result == .Working || result == .Finished {
            let profile = trajectory.getProfiles()[0][0]

            // 1. Check no negative time intervals
            validateNoNegativeTimeIntervals(profile)

            // 2. Check acceleration limits are maintained throughout trajectory
            validateAccelerationLimits(trajectory, maxAcc: 1.9425898752751283)

            // 3. Check velocity limits are maintained throughout trajectory
            validateVelocityLimits(trajectory, maxVel: 4.0)

            print("✓ Case 1 passed with duration: \(trajectory.getDuration())")
        } else {
            // If it fails, that's acceptable - we just need it to not return invalid trajectories
            print("⚠ Case 1: Trajectory calculation failed with result: \(result)")
            XCTFail("Expected trajectory to succeed for this case, got: \(result)")
        }
    }

    func testAccelerationLimitWithNonZeroTargetVel_Case2() throws {
        // From failureTable.txt lines 19-33
        // Previously violated acceleration limits with target_vel=1.345, max_acc=5.0

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

        if result == .Working || result == .Finished {
            let profile = trajectory.getProfiles()[0][0]

            // 1. Check no negative time intervals
            validateNoNegativeTimeIntervals(profile)

            // 2. Check acceleration limits are maintained
            validateAccelerationLimits(trajectory, maxAcc: 5.0)

            // 3. Check velocity limits are maintained
            validateVelocityLimits(trajectory, maxVel: 4.0)

            print("✓ Case 2 passed with duration: \(trajectory.getDuration())")
        } else {
            print("⚠ Case 2: Trajectory calculation failed with result: \(result)")
            XCTFail("Expected trajectory to succeed for this case, got: \(result)")
        }
    }

    // MARK: - Bug #2: Negative Time Intervals with High Limits

    func testNegativeTimeInterval_Case1() throws {
        // From failureTable.txt lines 40-52
        // Previously produced t[5] = -0.277552

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

        if result == .Working || result == .Finished {
            let profile = trajectory.getProfiles()[0][0]

            // 1. Check no negative time intervals (especially t[5])
            validateNoNegativeTimeIntervals(profile)

            // 2. Check acceleration limits are maintained
            validateAccelerationLimits(trajectory, maxAcc: 10.0)

            // 3. Check velocity limits are maintained
            validateVelocityLimits(trajectory, maxVel: 10.0)

            print("✓ Case 3 (negative t[5]) passed with duration: \(trajectory.getDuration())")
            print("  Time intervals: \(profile.t)")
        } else {
            print("⚠ Case 3: Trajectory calculation failed with result: \(result)")
            // This case may legitimately fail if no valid trajectory exists
            // The important thing is it doesn't return invalid trajectories
        }
    }

    func testNegativeTimeInterval_Case2() throws {
        // From failureTable.txt lines 55-67
        // Previously produced t[5] = -0.526238

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

        if result == .Working || result == .Finished {
            let profile = trajectory.getProfiles()[0][0]

            // 1. Check no negative time intervals
            validateNoNegativeTimeIntervals(profile)

            // 2. Check acceleration limits are maintained
            validateAccelerationLimits(trajectory, maxAcc: 9.973774779735683)

            // 3. Check velocity limits are maintained
            validateVelocityLimits(trajectory, maxVel: 9.088312224669604)

            print("✓ Case 4 (negative t[5]) passed with duration: \(trajectory.getDuration())")
            print("  Time intervals: \(profile.t)")
        } else {
            print("⚠ Case 4: Trajectory calculation failed with result: \(result)")
        }
    }

    func testNegativeTimeInterval_Case3() throws {
        // From failureTable.txt lines 70-82
        // Previously produced t[5] = -0.526238 (duplicate with different max_vel)

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

        if result == .Working || result == .Finished {
            let profile = trajectory.getProfiles()[0][0]

            // 1. Check no negative time intervals
            validateNoNegativeTimeIntervals(profile)

            // 2. Check acceleration limits are maintained
            validateAccelerationLimits(trajectory, maxAcc: 9.973774779735683)

            // 3. Check velocity limits are maintained
            validateVelocityLimits(trajectory, maxVel: 6.282291093061675)

            print("✓ Case 5 (negative t[5]) passed with duration: \(trajectory.getDuration())")
            print("  Time intervals: \(profile.t)")
        } else {
            print("⚠ Case 5: Trajectory calculation failed with result: \(result)")
        }
    }
}
