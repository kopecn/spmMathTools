import XCTest

@testable import spmMathTools

/// Tests for trajectory continuity at profile boundaries with non-zero boundary conditions
final class OTGContinuityTests: XCTestCase {

    // MARK: - Test Helpers

    /// Helper to integrate position, velocity, and acceleration over a time interval with constant jerk
    func integrate(_ t: Double, _ p0: Double, _ v0: Double, _ a0: Double, _ j: Double) -> (Double, Double, Double) {
        let p = p0 + t * (v0 + t * (a0 / 2 + t * j / 6))
        let v = v0 + t * (a0 + t * j / 2)
        let a = a0 + t * j
        return (p, v, a)
    }

    /// Helper to check trajectory continuity at profile boundaries
    func checkTrajectoryContinuity(_ trajectory: Trajectory, dof: Int = 0, tolerance: Double = 1e-6) -> Bool {
        guard trajectory.profiles.count > 0, trajectory.profiles[0].count > dof else {
            return false
        }

        let profile = trajectory.profiles[0][dof]

        print("\n🔍 CONTINUITY CHECK for DOF \(dof)")
        print("Profile type: \(profile)")
        print("t = [\(profile.t.map { String(format: "%.6f", $0) }.joined(separator: ", "))]")
        print("j = [\(profile.j.map { String(format: "%.3f", $0) }.joined(separator: ", "))]")
        print("Stored p/v/a values:")
        for i in 0...7 {
            print(
                "  [\(i)] p=\(String(format: "%12.9f", profile.p[i])), v=\(String(format: "%12.9f", profile.v[i])), a=\(String(format: "%12.9f", profile.a[i]))"
            )
        }

        var previousP = profile.p[0]
        var previousV = profile.v[0]
        var previousA = profile.a[0]

        // Check continuity at each profile boundary
        for i in 1...7 {
            let currentP = profile.p[i]
            let currentV = profile.v[i]
            let currentA = profile.a[i]

            // Integrate from previous step
            let t = profile.t[i - 1]
            let j = profile.j[i - 1]
            let (integratedP, integratedV, integratedA) = integrate(t, previousP, previousV, previousA, j)

            print("\nBoundary \(i-1)→\(i):")
            print(
                "  Previous: p=\(String(format: "%.9f", previousP)), v=\(String(format: "%.9f", previousV)), a=\(String(format: "%.9f", previousA))"
            )
            print("  Interval: t=\(String(format: "%.9f", t)), j=\(String(format: "%.3f", j))")
            print(
                "  Integrated: p=\(String(format: "%.9f", integratedP)), v=\(String(format: "%.9f", integratedV)), a=\(String(format: "%.9f", integratedA))"
            )
            print(
                "  Stored:     p=\(String(format: "%.9f", currentP)), v=\(String(format: "%.9f", currentV)), a=\(String(format: "%.9f", currentA))"
            )

            // Check if integration matches stored values
            if abs(integratedP - currentP) > tolerance {
                print(
                    "❌ Position discontinuity at boundary \(i-1)→\(i): integrated=\(integratedP), stored=\(currentP), diff=\(abs(integratedP - currentP))"
                )
                return false
            }
            if abs(integratedV - currentV) > tolerance {
                print(
                    "❌ Velocity discontinuity at boundary \(i-1)→\(i): integrated=\(integratedV), stored=\(currentV), diff=\(abs(integratedV - currentV))"
                )
                return false
            }
            if abs(integratedA - currentA) > tolerance {
                print(
                    "❌ Acceleration discontinuity at boundary \(i-1)→\(i): integrated=\(integratedA), stored=\(currentA), diff=\(abs(integratedA - currentA))"
                )
                return false
            }

            previousP = currentP
            previousV = currentV
            previousA = currentA
        }

        print("✓ Trajectory continuity verified for DOF \(dof)")
        return true
    }

    // MARK: - Continuity Tests

    func testNonZeroVelocityBoundaryConditions() throws {
        // Test case 1: Both initial and target velocities non-zero
        let otg = OTG(degreesOfFreedom: 1, deltaTime: 0.01)
        var input = InputParameter(DOFs: 1)

        input.currentPosition = [2.929012715930903]
        input.currentVelocity = [1.3368222168905952]
        input.currentAcceleration = [0.0]
        input.targetPosition = [10.0]
        input.targetVelocity = [-1.3503178982725528]
        input.targetAcceleration = [0.0]
        input.maxVelocity = [4.0]
        input.maxAcceleration = [5.0]
        input.maxJerk = [10.0]

        var trajectory = Trajectory(dofs: 1)
        let result = otg.calculate(input: input, trajectory: &trajectory)

        print("🟢 Test result: \(result)")
        XCTAssertTrue(result == .Working || result == .Finished, "Calculation should succeed, got: \(result)")
        XCTAssertTrue(checkTrajectoryContinuity(trajectory, dof: 0), "Trajectory must be continuous")

        // Verify boundary conditions
        let (startPos, startVel, startAcc) = trajectory.atTime(0.0)
        XCTAssertEqual(startPos[0], input.currentPosition[0], accuracy: 1e-6, "Initial position mismatch")
        XCTAssertEqual(startVel[0], input.currentVelocity[0], accuracy: 1e-6, "Initial velocity mismatch")
        XCTAssertEqual(startAcc[0], input.currentAcceleration[0], accuracy: 1e-6, "Initial acceleration mismatch")

        let (endPos, endVel, endAcc) = trajectory.atTime(trajectory.getDuration())
        XCTAssertEqual(endPos[0], input.targetPosition[0], accuracy: 1e-6, "Target position mismatch")
        XCTAssertEqual(endVel[0], input.targetVelocity[0], accuracy: 1e-6, "Target velocity mismatch")
        XCTAssertEqual(endAcc[0], input.targetAcceleration[0], accuracy: 1e-6, "Target acceleration mismatch")
    }

    func testNonZeroAccelerationBoundaryConditions() throws {
        // Test case 2: Both initial and target accelerations non-zero
        let otg = OTG(degreesOfFreedom: 1, deltaTime: 0.01)
        var input = InputParameter(DOFs: 1)

        input.currentPosition = [0.0]
        input.currentVelocity = [0.0]
        input.currentAcceleration = [-1.5686480326295587]
        input.targetPosition = [10.0]
        input.targetVelocity = [0.0]
        input.targetAcceleration = [2.3736654270633393]
        input.maxVelocity = [4.0]
        input.maxAcceleration = [5.0]
        input.maxJerk = [10.0]

        var trajectory = Trajectory(dofs: 1)
        let result = otg.calculate(input: input, trajectory: &trajectory)

        XCTAssertTrue(result == .Working || result == .Finished, "Calculation should succeed")
        XCTAssertTrue(checkTrajectoryContinuity(trajectory, dof: 0), "Trajectory must be continuous")

        // Verify boundary conditions
        let (startPos, startVel, startAcc) = trajectory.atTime(0.0)
        XCTAssertEqual(startPos[0], input.currentPosition[0], accuracy: 1e-6)
        XCTAssertEqual(startVel[0], input.currentVelocity[0], accuracy: 1e-6)
        XCTAssertEqual(startAcc[0], input.currentAcceleration[0], accuracy: 1e-6)

        let (endPos, endVel, endAcc) = trajectory.atTime(trajectory.getDuration())
        XCTAssertEqual(endPos[0], input.targetPosition[0], accuracy: 1e-6)
        XCTAssertEqual(endVel[0], input.targetVelocity[0], accuracy: 1e-6)
        XCTAssertEqual(endAcc[0], input.targetAcceleration[0], accuracy: 1e-6)
    }

    func testInvertedBoundaryConditions() throws {
        // Test case 3: Inverted initial and target (should also work)
        let otg = OTG(degreesOfFreedom: 1, deltaTime: 0.01)
        var input = InputParameter(DOFs: 1)

        input.currentPosition = [7.1647072936660265]
        input.currentVelocity = [-1.5912907869481767]
        input.currentAcceleration = [1.2187350047984646]
        input.targetPosition = [-0.1964371401151649]
        input.targetVelocity = [1.8250659788867563]
        input.targetAcceleration = [-1.9871641074856043]
        input.maxVelocity = [4.0]
        input.maxAcceleration = [5.0]
        input.maxJerk = [10.0]

        var trajectory = Trajectory(dofs: 1)
        let result = otg.calculate(input: input, trajectory: &trajectory)

        XCTAssertTrue(result == .Working || result == .Finished, "Calculation should succeed")
        XCTAssertTrue(checkTrajectoryContinuity(trajectory, dof: 0), "Trajectory must be continuous")

        // Verify boundary conditions
        let (startPos, startVel, startAcc) = trajectory.atTime(0.0)
        XCTAssertEqual(startPos[0], input.currentPosition[0], accuracy: 1e-6)
        XCTAssertEqual(startVel[0], input.currentVelocity[0], accuracy: 1e-6)
        XCTAssertEqual(startAcc[0], input.currentAcceleration[0], accuracy: 1e-6)

        let (endPos, endVel, endAcc) = trajectory.atTime(trajectory.getDuration())
        XCTAssertEqual(endPos[0], input.targetPosition[0], accuracy: 1e-6)
        XCTAssertEqual(endVel[0], input.targetVelocity[0], accuracy: 1e-6)
        XCTAssertEqual(endAcc[0], input.targetAcceleration[0], accuracy: 1e-6)
    }

    func testZeroTargetBoundaryConditions() throws {
        // Test case 4: Zero target conditions (known to work)
        let otg = OTG(degreesOfFreedom: 1, deltaTime: 0.01)
        var input = InputParameter(DOFs: 1)

        input.currentPosition = [2.5]
        input.currentVelocity = [1.2]
        input.currentAcceleration = [-0.8]
        input.targetPosition = [10.0]
        input.targetVelocity = [0.0]
        input.targetAcceleration = [0.0]
        input.maxVelocity = [4.0]
        input.maxAcceleration = [5.0]
        input.maxJerk = [10.0]

        var trajectory = Trajectory(dofs: 1)
        let result = otg.calculate(input: input, trajectory: &trajectory)

        XCTAssertTrue(result == .Working || result == .Finished, "Calculation should succeed")
        XCTAssertTrue(checkTrajectoryContinuity(trajectory, dof: 0), "Trajectory must be continuous")
    }

    func testZeroInitialBoundaryConditions() throws {
        // Test case 5: Zero initial conditions (known to work)
        let otg = OTG(degreesOfFreedom: 1, deltaTime: 0.01)
        var input = InputParameter(DOFs: 1)

        input.currentPosition = [0.0]
        input.currentVelocity = [0.0]
        input.currentAcceleration = [0.0]
        input.targetPosition = [10.0]
        input.targetVelocity = [1.5]
        input.targetAcceleration = [-2.0]
        input.maxVelocity = [4.0]
        input.maxAcceleration = [5.0]
        input.maxJerk = [10.0]

        var trajectory = Trajectory(dofs: 1)
        let result = otg.calculate(input: input, trajectory: &trajectory)

        XCTAssertTrue(result == .Working || result == .Finished, "Calculation should succeed")
        XCTAssertTrue(checkTrajectoryContinuity(trajectory, dof: 0), "Trajectory must be continuous")
    }
}
