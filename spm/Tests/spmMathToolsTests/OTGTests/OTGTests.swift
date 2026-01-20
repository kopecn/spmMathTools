import XCTest

@testable import spmMathTools

final class OTGTests: XCTestCase {

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

        if result != .Working && !(result == .Finished && output.trajectory.duration < 0.005) {
            print("[checkCalculation] FAILED: result=\(result), duration=\(output.trajectory.duration)")
        }
        XCTAssertTrue(result == .Working || (result == .Finished && output.trajectory.duration < 0.005))
        XCTAssertGreaterThanOrEqual(output.trajectory.duration, 0.0)

        for dof in 0..<otg.degreesOfFreedom {
            XCTAssertFalse(output.newPosition[dof].isNaN)
            XCTAssertFalse(output.newVelocity[dof].isNaN)
            XCTAssertFalse(output.newAcceleration[dof].isNaN)
        }
    }

    func stepThroughAndCheckCalculation(_ otg: OTG, _ input: inout InputParameter, maxNumberChecks: Int) -> Int {
        var output = OutputParameter(DOFs: 3)
        let otgSecond = OTG(degreesOfFreedom: 3, deltaTime: 0.001)

        checkCalculation(otg, input)

        var numberChecks = 1
        var currentInput = input

        while otg.update(input: currentInput, output: &output) == .Working {
            currentInput.currentPosition = output.newPosition
            currentInput.currentVelocity = output.newVelocity
            currentInput.currentAcceleration = output.newAcceleration

            checkCalculation(otgSecond, currentInput)

            numberChecks += 1
            if numberChecks == maxNumberChecks {
                break
            }
        }

        return numberChecks
    }

    // MARK: - Basic Trajectory Tests

    func testBasicTrajectory() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.005)
        var input = InputParameter(DOFs: 3)
        var output = OutputParameter(DOFs: 3)

        input.currentPosition = [0.0, -2.0, 0.0]
        input.currentVelocity = [0.0, 0.0, 0.0]
        input.currentAcceleration = [0.0, 0.0, 0.0]
        input.targetPosition = [1.0, -3.0, 2.0]
        input.targetVelocity = [0.0, 0.3, 0.0]
        input.targetAcceleration = [0.0, 0.0, 0.0]
        input.maxVelocity = [1.0, 1.0, 1.0]
        input.maxAcceleration = [1.0, 1.0, 1.0]
        input.maxJerk = [1.0, 1.0, 1.0]

        var trajectory = Trajectory(dofs: 3)
        var result = otg.calculate(input: input, trajectory: &trajectory)

        XCTAssertEqual(result, .Working)
        XCTAssertEqual(trajectory.duration, 4.0, accuracy: 1e-10)

        result = otg.update(input: input, output: &output)

        XCTAssertEqual(result, .Working)
        XCTAssertEqual(output.trajectory.duration, 4.0, accuracy: 1e-10)
    }

    func testAtTime() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.005)
        var input = InputParameter(DOFs: 3)
        var output = OutputParameter(DOFs: 3)

        input.currentPosition = [0.0, -2.0, 0.0]
        input.currentVelocity = [0.0, 0.0, 0.0]
        input.currentAcceleration = [0.0, 0.0, 0.0]
        input.targetPosition = [1.0, -3.0, 2.0]
        input.targetVelocity = [0.0, 0.3, 0.0]
        input.targetAcceleration = [0.0, 0.0, 0.0]
        input.maxVelocity = [1.0, 1.0, 1.0]
        input.maxAcceleration = [1.0, 1.0, 1.0]
        input.maxJerk = [1.0, 1.0, 1.0]

        // Save original values before they get modified by update
        let originalPosition = input.currentPosition
        let originalVelocity = input.currentVelocity
        let originalAcceleration = input.currentAcceleration
        let originalTargetPosition = input.targetPosition
        let originalTargetVelocity = input.targetVelocity
        let originalTargetAcceleration = input.targetAcceleration

        let result = otg.update(input: input, output: &output)
        XCTAssertEqual(result, .Working)

        // Test at start time
        var (newPosition, newVelocity, newAcceleration) = output.trajectory.atTime(0.0)

        XCTAssertTrue(arrayEquals(newPosition, originalPosition))
        XCTAssertTrue(arrayEquals(newVelocity, originalVelocity))
        XCTAssertTrue(arrayEquals(newAcceleration, originalAcceleration))

        // Test at end time
        (newPosition, newVelocity, newAcceleration) = output.trajectory.atTime(output.trajectory.duration)

        XCTAssertTrue(arrayEquals(newPosition, originalTargetPosition))
        XCTAssertTrue(arrayEquals(newVelocity, originalTargetVelocity))
        XCTAssertTrue(arrayEquals(newAcceleration, originalTargetAcceleration))

        // Test at specific time
        (newPosition, _, _) = output.trajectory.atTime(2.0)
        let expectedPosition = [0.5, -2.6871268303, 1.0]

        print("=== DEBUGGING TRAJECTORY CALCULATION ERROR ===")
        print("Expected DOF 1 at t=2.0: -2.6871268303")
        print("Actual DOF 1 at t=2.0: \(newPosition[1])")
        print("Error: \(abs(newPosition[1] - expectedPosition[1]))")
        print("Trajectory duration: \(output.trajectory.duration) (C++ expects 4.0)")

        // Check boundary conditions
        let (startPos, startVel, startAcc) = output.trajectory.atTime(0.0)
        let (endPos, endVel, endAcc) = output.trajectory.atTime(output.trajectory.duration)
        print("Boundary check DOF 1:")
        print(
            "  Start: pos=\(startPos[1]) (expect -2.0), vel=\(startVel[1]) (expect 0.0), acc=\(startAcc[1]) (expect 0.0)"
        )
        print("  End: pos=\(endPos[1]) (expect -3.0), vel=\(endVel[1]) (expect 0.3), acc=\(endAcc[1]) (expect 0.0)")

        // Check independent minimum durations
        print("Independent minimum durations: \(output.trajectory.independentMinDurations)")
        print("  C++ test expects DOF 1 independent duration: 3.6860977315")
        print(
            "  Error in DOF 1 independent duration: \(abs(output.trajectory.independentMinDurations[1] - 3.6860977315))"
        )

        // Test DOF 1 in isolation to see if it's a synchronization issue
        print("\n=== Testing DOF 1 in isolation ===")
        let otgSingle = OTG(degreesOfFreedom: 1, deltaTime: 0.005)
        var inputSingle = InputParameter(DOFs: 1)
        var outputSingle = OutputParameter(DOFs: 1)

        inputSingle.currentPosition = [-2.0]
        inputSingle.currentVelocity = [0.0]
        inputSingle.currentAcceleration = [0.0]
        inputSingle.targetPosition = [-3.0]
        inputSingle.targetVelocity = [0.3]
        inputSingle.targetAcceleration = [0.0]
        inputSingle.maxVelocity = [1.0]
        inputSingle.maxAcceleration = [1.0]
        inputSingle.maxJerk = [1.0]

        let resultSingle = otgSingle.update(input: inputSingle, output: &outputSingle)
        print("Single DOF result: \(resultSingle)")
        print("Single DOF duration: \(outputSingle.trajectory.duration)")
        print("Single DOF independent duration: \(outputSingle.trajectory.independentMinDurations[0])")

        let (singlePos, _, _) = outputSingle.trajectory.atTime(2.0)
        print("Single DOF at t=2.0: \(singlePos[0])")

        if output.trajectory.profiles.count > 0 && output.trajectory.profiles[0].count > 1 {
            let profile1 = output.trajectory.profiles[0][1]
            print("\nDOF 1 Profile Details:")
            print("  Initial: p0=\(profile1.p[0]), v0=\(profile1.v[0]), a0=\(profile1.a[0])")
            print("  Final: pf=\(profile1.pf), vf=\(profile1.vf), af=\(profile1.af)")
            print("  Times t: \(profile1.t)")
            print("  Cumulative tSum: \(profile1.tSum)")
            print("  Jerks j: \(profile1.j)")

            // Calculate manually what position should be at t=2.0
            print("\nManual calculation at t=2.0:")
            let t = 2.0
            print("  Target time: \(t)")

            // Find which phase we're in
            for i in 0..<profile1.tSum.count {
                if t <= profile1.tSum[i] {
                    let t_phase = (i == 0) ? t : t - profile1.tSum[i - 1]
                    print("  Phase \(i): t_phase=\(t_phase), tSum[\(i)]=\(profile1.tSum[i])")
                    print(
                        "  Phase \(i): p[\(i)]=\(profile1.p[i]), v[\(i)]=\(profile1.v[i]), a[\(i)]=\(profile1.a[i]), j[\(i)]=\(profile1.j[i])"
                    )

                    // Manual integrate calculation
                    let p0 = profile1.p[i]
                    let v0 = profile1.v[i]
                    let a0 = profile1.a[i]
                    let j = profile1.j[i]

                    let pos = p0 + t_phase * (v0 + t_phase * (a0 / 2 + t_phase * j / 6))
                    print(
                        "  Manual calc: pos = \(p0) + \(t_phase) * (\(v0) + \(t_phase) * (\(a0)/2 + \(t_phase) * \(j)/6)) = \(pos)"
                    )
                    break
                }
            }
        }

        // With the extrema sorting fix, this should now pass with proper tolerance
        XCTAssertTrue(arrayEquals(newPosition, expectedPosition, accuracy: 1e-6))
    }

    func testSingleDOF() throws {
        let otg = OTG(degreesOfFreedom: 1, deltaTime: 0.005)
        var input = InputParameter(DOFs: 1)

        input.currentPosition = [0.0]
        input.targetPosition = [1.0]
        input.maxVelocity = [1.0]
        input.maxAcceleration = [1.0]
        input.maxJerk = [1.0]

        var trajectory = Trajectory(dofs: 1)
        let result = otg.calculate(input: input, trajectory: &trajectory)

        XCTAssertEqual(result, .Working)
        XCTAssertEqual(trajectory.duration, 3.1748, accuracy: 1e-4)

        let (newPosition, _, _) = trajectory.atTime(0.0)
        XCTAssertTrue(arrayEquals(newPosition, input.currentPosition))

        let (midPosition, _, _) = trajectory.atTime(3.1748 / 2)
        XCTAssertEqual(midPosition[0], 0.5, accuracy: 1e-4)
    }

    func testIndependentMinDurations() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.005)
        var input = InputParameter(DOFs: 3)
        var output = OutputParameter(DOFs: 3)

        input.currentPosition = [0.0, -2.0, 0.0]
        input.currentVelocity = [0.0, 0.0, 0.0]
        input.currentAcceleration = [0.0, 0.0, 0.0]
        input.targetPosition = [1.0, -3.0, 2.0]
        input.targetVelocity = [0.0, 0.3, 0.0]
        input.targetAcceleration = [0.0, 0.0, 0.0]
        input.maxVelocity = [1.0, 1.0, 1.0]
        input.maxAcceleration = [1.0, 1.0, 1.0]
        input.maxJerk = [1.0, 1.0, 1.0]

        let result = otg.update(input: input, output: &output)
        XCTAssertEqual(result, .Working)

        let independentMinDurations = output.trajectory.independentMinDurations
        XCTAssertEqual(independentMinDurations[0], 3.1748021039, accuracy: 1e-6)
        XCTAssertEqual(independentMinDurations[1], 3.6860977315, accuracy: 1e-6)
        XCTAssertEqual(independentMinDurations[2], output.trajectory.duration, accuracy: 1e-6)
    }

    // MARK: - Input Validation Tests

    func testInputValidation() throws {
        let otg = OTG(degreesOfFreedom: 2)
        var input = InputParameter(DOFs: 2)

        input.currentPosition = [0.0, -2.0]
        input.currentVelocity = [0.0, 0.0]
        input.currentAcceleration = [0.0, 0.0]
        input.targetPosition = [1.0, -3.0]
        input.targetVelocity = [0.0, 0.3]
        input.targetAcceleration = [0.0, 0.0]
        input.maxVelocity = [1.0, 1.0]
        input.maxAcceleration = [1.0, 1.0]
        input.maxJerk = [1.0, 1.0]

        XCTAssertTrue(try otg.validateInput(input: input))

        // Test with NaN
        input.maxJerk = [1.0, Double.nan]
        XCTAssertNil(try? otg.validateInput(input: input))

        // Test with negative limits
        input.maxJerk = [1.0, 1.0]
        input.maxAcceleration = [1.0, -1.0]
        XCTAssertNil(try? otg.validateInput(input: input))

        // Test velocity exceeding limits
        input.maxAcceleration = [1.0, 1.0]
        input.targetVelocity = [0.0, 1.3]
        XCTAssertNil(try? otg.validateInput(input: input))
    }

    // MARK: - Enabled DOFs Tests

    func testEnabledDOFs() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.005)
        var input = InputParameter(DOFs: 3)
        var output = OutputParameter(DOFs: 3)

        input.enabled = [true, false, false]
        input.currentPosition = [0.0, -2.0, 0.0]
        input.currentVelocity = [0.0, 0.1, 0.0]
        input.currentAcceleration = [0.0, 0.0, -0.2]
        input.targetPosition = [1.0, -3.0, 2.0]
        input.maxVelocity = [1.0, 1.0, 1.0]
        input.maxAcceleration = [1.0, 1.0, 1.0]
        input.maxJerk = [1.0, 1.0, 1.0]

        // Save original values before they get modified by update
        let originalPosition = input.currentPosition
        let originalVelocity = input.currentVelocity
        let originalAcceleration = input.currentAcceleration

        let result = otg.update(input: input, output: &output)
        XCTAssertEqual(result, .Working)
        XCTAssertEqual(output.trajectory.duration, 3.1748021039, accuracy: 1e-6)

        let (startPosition, startVelocity, startAcceleration) = output.trajectory.atTime(0.0)

        XCTAssertTrue(arrayEquals(startPosition, originalPosition))
        XCTAssertTrue(arrayEquals(startVelocity, originalVelocity))
        XCTAssertTrue(arrayEquals(startAcceleration, originalAcceleration))

        let (endPosition, _, _) = output.trajectory.atTime(output.trajectory.duration)
        let expectedEndPosition = [input.targetPosition[0], -1.6825197896, -1.0079368399]
        XCTAssertTrue(arrayEquals(endPosition, expectedEndPosition, accuracy: 1e-6))
    }

    // MARK: - Phase Synchronization Tests

    func testPhaseSynchronization() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.005)
        var input = InputParameter(DOFs: 3)
        var output = OutputParameter(DOFs: 3)

        input.currentPosition = [0.0, -2.0, 0.0]
        input.targetPosition = [1.0, -3.0, 2.0]
        input.maxVelocity = [1.0, 1.0, 1.0]
        input.maxAcceleration = [1.0, 1.0, 1.0]
        input.maxJerk = [1.0, 1.0, 1.0]
        input.synchronization = .Phase

        var trajectory = Trajectory(dofs: 3)
        var result = otg.calculate(input: input, trajectory: &trajectory)
        XCTAssertEqual(result, .Working)
        XCTAssertEqual(trajectory.duration, 4.0, accuracy: 1e-6)

        // Test that all profiles have same timing
        let profiles = trajectory.profiles[0]
        XCTAssertTrue(arrayEquals(profiles[0].t, profiles[1].t))
        XCTAssertTrue(arrayEquals(profiles[0].t, profiles[2].t))

        result = otg.update(input: input, output: &output)
        let (position1s, _, _) = output.trajectory.atTime(1.0)
        let expectedPosition1s = [0.0833333333, -2.0833333333, 0.1666666667]
        XCTAssertTrue(arrayEquals(position1s, expectedPosition1s, accuracy: 1e-6))

        // Test equal start and target state
        input.currentPosition = [1.0, -2.0, 3.0]
        input.targetPosition = [1.0, -2.0, 3.0]
        input.currentVelocity = [0.0, 0.0, 0.0]
        input.currentAcceleration = [0.0, 0.0, 0.0]
        input.targetVelocity = [0.0, 0.0, 0.0]
        input.targetAcceleration = [0.0, 0.0, 0.0]
        otg.reset()  // Force recalculation
        result = otg.update(input: input, output: &output)
        XCTAssertEqual(result, .Finished)
        XCTAssertEqual(output.trajectory.duration, 0.0, accuracy: 1e-6)
    }

    // MARK: - Discretization Tests

    func testDiscretization() throws {
        // First test: DOF 1 in isolation to debug the issue
        print("\n=== Testing DOF 1 with discretization ===")
        let otgSingle = OTG(degreesOfFreedom: 1, deltaTime: 0.01)
        var inputSingle = InputParameter(DOFs: 1)

        inputSingle.currentPosition = [0.0]
        inputSingle.targetPosition = [-3.0]
        inputSingle.targetVelocity = [0.2]
        inputSingle.maxVelocity = [1.0]
        inputSingle.maxAcceleration = [2.0]
        inputSingle.maxJerk = [2.4]

        // Test without discretization first
        inputSingle.durationDiscretization = .Continuous
        var trajSingleCont = Trajectory(dofs: 1)
        let resultSingleCont = otgSingle.calculate(input: inputSingle, trajectory: &trajSingleCont)
        print("Without discretization: result=\(resultSingleCont), duration=\(trajSingleCont.duration)")

        // Now with discretization
        inputSingle.durationDiscretization = .Discrete
        var trajSingleDisc = Trajectory(dofs: 1)
        let resultSingleDisc = otgSingle.calculate(input: inputSingle, trajectory: &trajSingleDisc)
        print("With discretization: result=\(resultSingleDisc), duration=\(trajSingleDisc.duration)")
        print("Expected: result=Working, duration=4.5")

        // Full 3-DOF test
        print("\n=== Testing 3 DOFs with discretization ===")
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.01)
        var input = InputParameter(DOFs: 3)

        input.currentPosition = [0.0, 0.0, 0.0]
        input.targetPosition = [1.0, -3.0, 2.0]
        input.targetVelocity = [0.2, 0.2, 0.2]
        input.maxVelocity = [1.0, 1.0, 1.0]
        input.maxAcceleration = [2.0, 2.0, 2.0]
        input.maxJerk = [1.8, 2.4, 2.0]
        input.durationDiscretization = .Discrete

        var trajectory = Trajectory(dofs: 3)
        let result = otg.calculate(input: input, trajectory: &trajectory)

        print("3-DOF result: \(result), duration: \(trajectory.duration)")
        print("Expected: result=Working, duration=4.5")

        XCTAssertEqual(result, .Working)
        XCTAssertEqual(trajectory.duration, 4.5, accuracy: 1e-6)

        let (endPosition, _, _) = trajectory.atTime(4.5)
        XCTAssertTrue(arrayEquals(endPosition, [1.0, -3.0, 2.0]))
    }

    // MARK: - Per-DOF Settings Tests

    func testPerDOFSettings() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.005)
        var input = InputParameter(DOFs: 3)

        input.currentPosition = [0.0, -2.0, 0.0]
        input.currentVelocity = [0.0, 0.0, 0.0]
        input.currentAcceleration = [0.0, 0.0, 0.0]
        input.targetPosition = [1.0, -3.0, 2.0]
        input.targetVelocity = [0.0, 0.3, 0.0]
        input.targetAcceleration = [0.0, 0.0, 0.0]
        input.maxVelocity = [1.0, 1.0, 1.0]
        input.maxAcceleration = [1.0, 1.0, 1.0]
        input.maxJerk = [1.0, 1.0, 1.0]

        var trajectory = Trajectory(dofs: 3)
        var result = otg.calculate(input: input, trajectory: &trajectory)
        XCTAssertEqual(result, .Working)
        XCTAssertEqual(trajectory.duration, 4.0, accuracy: 1e-6)

        let (position2s, _, _) = trajectory.atTime(2.0)
        let expectedPosition2s = [0.5, -2.6871268303, 1.0]

        XCTAssertTrue(arrayEquals(position2s, expectedPosition2s, accuracy: 1e-6))

        // Test velocity interface
        input.controlInterface = .Velocity
        result = otg.calculate(input: input, trajectory: &trajectory)

        XCTAssertEqual(result, .Working)
        XCTAssertEqual(trajectory.duration, 1.095445115, accuracy: 1e-6)

        // Test per-DOF control interface
        input.perDofControlInterface = [.Position, .Velocity, .Position]
        result = otg.calculate(input: input, trajectory: &trajectory)

        XCTAssertEqual(result, .Working)
        XCTAssertEqual(trajectory.duration, 4.0, accuracy: 1e-6)

        // Test per-DOF synchronization
        input.perDofSynchronization = [.Time, .None, .Time]
        result = otg.calculate(input: input, trajectory: &trajectory)
        XCTAssertEqual(result, .Working)
        XCTAssertEqual(trajectory.duration, 4.0, accuracy: 1e-6)
    }

    // MARK: - Zero Limits Tests

    func testZeroLimits() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.005)
        var input = InputParameter(DOFs: 3)
        var output = OutputParameter(DOFs: 3)

        input.currentPosition = [0.0, -2.0, 0.0]
        input.currentVelocity = [0.2, 0.0, 0.0]
        input.currentAcceleration = [0.0, 0.0, 0.0]
        input.targetPosition = [1.0, -3.0, 0.0]
        input.targetVelocity = [0.2, 0.0, 0.0]
        input.targetAcceleration = [0.0, 0.0, 0.0]
        input.maxVelocity = [1.0, 1.0, 1.0]
        input.maxAcceleration = [0.0, 1.0, 0.0]
        input.maxJerk = [0.0, 1.0, 0.0]

        let result = otg.update(input: input, output: &output)

        XCTAssertEqual(result, .Working)
        XCTAssertEqual(output.trajectory.duration, 5.0, accuracy: 1e-6)
    }

    // MARK: - Velocity Interface Tests

    func testVelocityInterface() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.005)
        var input = InputParameter(DOFs: 3)

        input.controlInterface = .Velocity
        input.currentPosition = [0.0, 0.0, 0.0]
        input.currentVelocity = [-0.2, 0.0, 0.0]
        input.currentAcceleration = [1.0, 0.0, 0.2]
        input.targetVelocity = [0.9, 0.5, 0.4]
        input.targetAcceleration = [1.0, 0.0, 0.2]
        input.maxAcceleration = [1.0, 2.0, 6.0]
        input.maxJerk = [1.0, 2.0, 20.0]

        checkCalculation(otg, input)
    }

    // MARK: - Step Through Tests

    func testStepThrough() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.01)

        // Test a few step-through scenarios
        for i in 0..<10 {
            var input = InputParameter(DOFs: 3)

            let seed = Double(i * 789 + 123)

            input.currentPosition = [
                sin(seed) * 4.0,
                cos(seed + 1) * 4.0,
                sin(seed + 2) * 4.0,
            ]

            input.targetPosition = [
                sin(seed + 9) * 4.0,
                cos(seed + 10) * 4.0,
                sin(seed + 11) * 4.0,
            ]

            input.maxVelocity = [1.0, 1.0, 1.0]
            input.maxAcceleration = [1.0, 1.0, 1.0]
            input.maxJerk = [1.0, 1.0, 1.0]

            if (try? otg.validateInput(input: input)) ?? false {
                let checks = stepThroughAndCheckCalculation(otg, &input, maxNumberChecks: 100)
                XCTAssertGreaterThan(checks, 0)
            }
        }
    }

}

// MARK: - Test Extensions

extension OTGTests {

    func testErrorHandling() throws {
        let otg = OTG(degreesOfFreedom: 3, )
        var input = InputParameter(DOFs: 3)
        var output = OutputParameter(DOFs: 3)

        // Test with impossible target velocity
        input.currentPosition = [0.0, -2.0, 0.0]
        input.targetVelocity = [2.0, 0.3, 0.0]  // Exceeds max velocity
        input.maxVelocity = [1.0, 1.0, 1.0]
        input.maxAcceleration = [1.0, 1.0, 1.0]
        input.maxJerk = [1.0, 1.0, 1.0]

        let result = otg.update(input: input, output: &output)
        // Should handle error gracefully
        XCTAssertNotEqual(result, .Working)
    }

    func testMinimumDuration() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.005)
        var input = InputParameter(DOFs: 3)
        var output = OutputParameter(DOFs: 3)

        input.currentPosition = [0.0, -2.0, 0.0]
        input.targetPosition = [1.0, -3.0, 2.0]
        input.targetVelocity = [0.2, -0.3, 0.8]
        input.maxVelocity = [1.0, 1.0, 1.0]
        input.maxAcceleration = [1.0, 1.0, 1.0]
        input.maxJerk = [1.0, 1.0, 1.0]
        input.minimumDuration = 12.0

        let result = otg.update(input: input, output: &output)

        XCTAssertEqual(result, .Working)
        XCTAssertEqual(output.trajectory.duration, 12.0, accuracy: 1e-6)
    }

    func testHighLimits() throws {
        let otg = OTG(degreesOfFreedom: 3, deltaTime: 0.005)
        var input = InputParameter(DOFs: 3)
        var output = OutputParameter(DOFs: 3)

        input.currentPosition = [1300.0, 0.0, 0.02]
        input.currentVelocity = [1200.0, 0.0, 0.0]
        input.currentAcceleration = [0.0, 0.0, 0.0]
        input.targetPosition = [1400.0, 0.0, 0.02]
        input.targetVelocity = [0.0, 0.0, 0.0]
        input.targetAcceleration = [0.0, 0.0, 0.0]
        input.maxVelocity = [800.0, 1.0, 1.0]
        input.maxAcceleration = [40000.0, 1.0, 1.0]
        input.maxJerk = [200000.0, 1.0, 1.0]

        let result = otg.update(input: input, output: &output)
        XCTAssertEqual(result, .Working)
        XCTAssertEqual(output.trajectory.duration, 0.167347, accuracy: 1e-6)

        let independentMinDurations = output.trajectory.independentMinDurations
        XCTAssertEqual(independentMinDurations[0], output.trajectory.duration, accuracy: 1e-6)
        XCTAssertEqual(independentMinDurations[1], 0.0, accuracy: 1e-6)
        XCTAssertEqual(independentMinDurations[2], 0.0, accuracy: 1e-6)
    }
}
