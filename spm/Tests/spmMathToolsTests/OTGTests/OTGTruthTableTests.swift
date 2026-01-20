import XCTest

@testable import spmMathTools

/// Truth table tests for successful trajectory calculations
/// These tests verify that trajectories can be successfully calculated for various
/// combinations of initial and target states, and kinematic limits
final class OTGTruthTableTests: XCTestCase {

    // MARK: - Truth Table Data

    /// Truth table: (currentPos, currentVel, currentAcc, targetPos, targetVel, targetAcc, maxVel, maxAcc, maxJerk, expectedDuration, expectedTimeIntervals)
    let truthTable:
        [(
            currentPos: Double,
            currentVel: Double,
            currentAcc: Double,
            targetPos: Double,
            targetVel: Double,
            targetAcc: Double,
            maxVel: Double,
            maxAcc: Double,
            maxJerk: Double,
            expectedDuration: Double,
            expectedTimeIntervals: [Double]
        )] = [
            // Test Case 1: Zero initial state to positive target (vel & acc limits hit)
            (
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

            // Test Case 2: Negative position to positive target (vel & acc limits hit)
            (
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

            // Test Case 3: Negative initial velocity (vel & acc limits hit)
            (
                currentPos: 0.0,
                currentVel: -2.933699789068232,
                currentAcc: 0.0,
                targetPos: 10.0,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 4.785261,
                expectedTimeIntervals: [0.500000, 0.886740, 0.500000, 1.598521, 0.500000, 0.300000, 0.500000]
            ),

            // Test Case 4: Negative initial acceleration (vel & acc limits hit)
            (
                currentPos: 0.0,
                currentVel: 0.0,
                currentAcc: -2.2867181768158473,
                targetPos: 10.0,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 4.108977,
                expectedTimeIntervals: [0.728672, 0.352291, 0.500000, 1.228014, 0.500000, 0.300000, 0.500000]
            ),

            // Test Case 5: Higher velocity limit (acc limit hit, vel limit not hit)
            (
                currentPos: -6.139891324284666,
                currentVel: 0.0,
                currentAcc: 0.0,
                targetPos: 10.0,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 9.097991562729273,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 4.127935,
                expectedTimeIntervals: [0.500000, 1.063968, 0.500000, 0.000000, 0.500000, 1.063968, 0.500000]
            ),

            // Test Case 6: Higher acceleration limit (vel limit hit, acc limit not hit)
            (
                currentPos: -6.139891324284666,
                currentVel: 0.0,
                currentAcc: 0.0,
                targetPos: 10.0,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 7.24409333730741,
                maxJerk: 10.0,
                expectedDuration: 5.299884,
                expectedTimeIntervals: [0.632456, 0.000000, 0.632456, 2.770062, 0.632456, 0.000000, 0.632456]
            ),

            // Test Case 7: High velocity and acceleration limits (no limits hit)
            (
                currentPos: -6.139891324284666,
                currentVel: 0.0,
                currentAcc: 0.0,
                targetPos: 10.0,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 10.0,
                maxAcc: 10.0,
                maxJerk: 10.0,
                expectedDuration: 3.724062,
                expectedTimeIntervals: [0.931015, 0.000000, 1.862031, 0.000000, 0.000000, 0.000000, 0.931015]
            ),

            // Test Case 8: Negative position and velocity (vel & acc limits hit)
            (
                currentPos: -6.139891324284666,
                currentVel: -2.9076772285399857,
                currentAcc: 0.0,
                targetPos: 10.0,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 6.309603,
                expectedTimeIntervals: [0.500000, 0.881535, 0.500000, 3.128067, 0.500000, 0.300000, 0.500000]
            ),

            // Test Case 9: Negative position and acceleration (vel & acc limits hit)
            (
                currentPos: -6.139891324284666,
                currentVel: 0.0,
                currentAcc: -2.5220102714600143,
                targetPos: 10.0,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 5.686552,
                expectedTimeIntervals: [0.752201, 0.363605, 0.500000, 2.770746, 0.500000, 0.300000, 0.500000]
            ),

            // Test Case 10: All negative initial states (vel & acc limits hit)
            (
                currentPos: -6.139891324284666,
                currentVel: -2.4332813646368305,
                currentAcc: -2.3549270909757887,
                targetPos: 10.0,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 6.619806,
                expectedTimeIntervals: [0.735493, 0.842113, 0.500000, 3.242200, 0.500000, 0.300000, 0.500000]
            ),

            // Test Case 11: Lower acceleration limit with symmetric profile (vel & acc limits hit)
            (
                currentPos: -6.139891324284666,
                currentVel: 0.0,
                currentAcc: 0.0,
                targetPos: 10.0,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 3.058381671863536,
                maxJerk: 10.0,
                expectedDuration: 5.648692,
                expectedTimeIntervals: [0.305838, 1.002043, 0.305838, 2.421253, 0.305838, 1.002043, 0.305838]
            ),

            // Test Case 12: Complex state transitions (negative to positive)
            (
                currentPos: -2.2569699192956714,
                currentVel: -2.2910170579603815,
                currentAcc: -2.2093383162142333,
                targetPos: 6.790054108584006,
                targetVel: 0.5880296221570069,
                targetAcc: 0.57461711298606,
                maxVel: 3.9501329787234045,
                maxAcc: 10.0,
                maxJerk: 10.0,
                expectedDuration: 4.516839,
                expectedTimeIntervals: [1.026242, 0.000000, 0.805308, 1.465311, 0.581258, 0.000000, 0.638720]
            ),

            // Test Case 13: Complex state transitions (positive to negative) A
            (
                currentPos: -2.2569699192956714,
                currentVel: 0.6924924339691851,
                currentAcc: 0.6756694790902418,
                targetPos: 6.790054108584006,
                targetVel: -1.159379585473221,
                targetAcc: -1.2629539618488628,
                maxVel: 3.9501329787234045,
                maxAcc: 10.0,
                maxJerk: 10.0,
                expectedDuration: 3.491257,
                expectedTimeIntervals: [0.505187, 0.000000, 0.572754, 1.098882, 0.720366, 0.000000, 0.594070]
            ),

            // Test Case 14: Complex state transitions with higher velocity limit
            (
                currentPos: -2.2569699192956714,
                currentVel: 0.6924924339691851,
                currentAcc: 0.6756694790902418,
                targetPos: 6.790054108584006,
                targetVel: -1.159379585473221,
                targetAcc: -1.2629539618488628,
                maxVel: 5.783597189104916,
                maxAcc: 10.0,
                maxJerk: 10.0,
                expectedDuration: 3.001705,
                expectedTimeIntervals: [0.647551, 0.000000, 0.715118, 0.089298, 0.838017, 0.000000, 0.711722]
            ),

            // Test Case 15: Mixed velocity and acceleration signs (negative vel, positive acc)
            (
                currentPos: -2.2569699192956714,
                currentVel: -0.9487630685986792,
                currentAcc: 0.6756694790902418,
                targetPos: 6.790054108584006,
                targetVel: 0.6247420671313275,
                targetAcc: -1.2629539618488628,
                maxVel: 4.0,
                maxAcc: 10.0,
                maxJerk: 10.0,
                expectedDuration: 3.453726,
                expectedTimeIntervals: [0.637528, 0.000000, 0.705095, 1.061812, 0.587793, 0.000000, 0.461498]
            ),

            // Test Case 16: Mixed velocity and acceleration signs (positive vel, negative acc) A
            (
                currentPos: -2.2569699192956714,
                currentVel: 1.0843211206896552,
                currentAcc: -0.908898110785033,
                targetPos: 6.790054108584006,
                targetVel: -0.8658519809244314,
                targetAcc: 0.8345847853998534,
                maxVel: 4.0,
                maxAcc: 10.0,
                maxJerk: 10.0,
                expectedDuration: 3.690308,
                expectedTimeIntervals: [0.634671, 0.000000, 0.543782, 1.028300, 0.700048, 0.000000, 0.783507]
            ),

            // Test Case 17: Mixed velocity and acceleration signs (positive vel, negative acc) B
            (
                currentPos: -2.2569699192956714,
                currentVel: 1.0843211206896552,
                currentAcc: -0.908898110785033,
                targetPos: 6.790054108584006,
                targetVel: 0.7905642424798243,
                targetAcc: -1.0706220194424065,
                maxVel: 4.0,
                maxAcc: 10.0,
                maxJerk: 10.0,
                expectedDuration: 3.110478,
                expectedTimeIntervals: [0.634671, 0.000000, 0.543782, 0.895978, 0.571555, 0.000000, 0.464492]
            ),

            // Test Case 18: Positive to negative target (vel & acc limits hit)
            (
                currentPos: 2.6603799559471364,
                currentVel: 0.0,
                currentAcc: 0.0,
                targetPos: -6.297064289647577,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 3.539361,
                expectedTimeIntervals: [0.500000, 0.300000, 0.500000, 0.939361, 0.500000, 0.300000, 0.500000]
            ),

            // Test Case 19: Positive to negative with negative initial velocity (vel & acc limits hit)
            (
                currentPos: 2.6603799559471364,
                currentVel: -1.1073731185756241,
                currentAcc: 0.0,
                targetPos: -6.297064289647577,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 3.279332,
                expectedTimeIntervals: [0.500000, 0.078525, 0.500000, 0.900807, 0.500000, 0.300000, 0.500000]
            ),

            // Test Case 20: Positive to negative with negative initial acceleration (vel & acc limits hit)
            (
                currentPos: 2.6603799559471364,
                currentVel: 0.0,
                currentAcc: -1.3502374724669601,
                targetPos: -6.297064289647577,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 3.426422,
                expectedTimeIntervals: [0.364976, 0.318231, 0.500000, 0.943215, 0.500000, 0.300000, 0.500000]
            ),

            // Test Case 21: Positive to negative with higher velocity limit (acc limit hit, vel limit not hit)
            (
                currentPos: 2.6603799559471364,
                currentVel: 0.0,
                currentAcc: 0.0,
                targetPos: -6.297064289647577,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 9.363831910792952,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 3.223225,
                expectedTimeIntervals: [0.500000, 0.611613, 0.500000, 0.000000, 0.500000, 0.611613, 0.500000]
            ),

            // Test Case 22: Positive to negative with higher acceleration limit (vel limit hit, acc limit not hit)
            (
                currentPos: 2.6603799559471364,
                currentVel: 0.0,
                currentAcc: 0.0,
                targetPos: -6.297064289647577,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 7.814659622797357,
                maxJerk: 10.0,
                expectedDuration: 3.504272,
                expectedTimeIntervals: [0.632456, 0.000000, 0.632456, 0.974450, 0.632456, 0.000000, 0.632456]
            ),

            // Test Case 23: Positive to negative with high velocity and acceleration limits (no limits hit)
            (
                currentPos: 2.6603799559471364,
                currentVel: 0.0,
                currentAcc: 0.0,
                targetPos: -6.297064289647577,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 9.045025123898679,
                maxAcc: 9.155501445484582,
                maxJerk: 10.0,
                expectedDuration: 3.060399,
                expectedTimeIntervals: [0.765100, 0.000000, 1.530199, 0.000000, 0.000000, 0.000000, 0.765100]
            ),

            // Test Case 24: Positive to negative with negative velocity (vel & acc limits hit)
            (
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

            // Test Case 25: Positive to negative with negative acceleration (vel & acc limits hit)
            (
                currentPos: 2.6603799559471364,
                currentVel: 0.0,
                currentAcc: -3.347616097650514,
                targetPos: -6.297064289647577,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 3.328272,
                expectedTimeIntervals: [0.165238, 0.412065, 0.500000, 0.950968, 0.500000, 0.300000, 0.500000]
            ),

            // Test Case 26: Positive to negative with negative velocity and acceleration (vel & acc limits hit)
            (
                currentPos: 2.6603799559471364,
                currentVel: -2.5039578744493394,
                currentAcc: -3.347616097650514,
                targetPos: -6.297064289647577,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 2.966020,
                expectedTimeIntervals: [0.118710, 0.000000, 0.453472, 1.093838, 0.500000, 0.300000, 0.500000]
            ),

            // Test Case 27: Positive to negative with negative velocity and acceleration variant (vel & acc limits hit)
            (
                currentPos: 2.6603799559471364,
                currentVel: -2.5536320668135097,
                currentAcc: -3.8897875367107195,
                targetPos: -6.297064289647577,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 2.958146,
                expectedTimeIntervals: [0.080371, 0.000000, 0.469350, 1.108425, 0.500000, 0.300000, 0.500000]
            ),

            // Test Case 28: Positive to negative at max velocity with negative acceleration (vel & acc limits hit)
            (
                currentPos: 2.6603799559471364,
                currentVel: -4.0,
                currentAcc: -3.0555937958883996,
                targetPos: -6.297064289647577,
                targetVel: 0.0,
                targetAcc: 0.0,
                maxVel: 4.0,
                maxAcc: 5.0,
                maxJerk: 10.0,
                expectedDuration: 2.890803,
                expectedTimeIntervals: [0.521622, 0.000000, 0.216063, 0.241999, 0.500000, 0.300000, 0.500000]
            ),

            // Test Case 29: Positive to negative with complex state transitions (negative to positive targets)
            (
                currentPos: 2.6603799559471364,
                currentVel: -2.9630139500734214,
                currentAcc: -3.0555937958883996,
                targetPos: -6.297064289647577,
                targetVel: 1.0698019915565349,
                targetAcc: 1.8217407764317184,
                maxVel: 4.0,
                maxAcc: 9.973774779735683,
                maxJerk: 10.0,
                expectedDuration: 2.993353,
                expectedTimeIntervals: [0.082232, 0.000000, 0.387791, 1.258336, 0.723584, 0.000000, 0.541410]
            ),

            // Test Case 30: Positive to negative with complex state transitions (positive to negative targets)
            (
                currentPos: 2.6603799559471364,
                currentVel: 1.9695301027900147,
                currentAcc: 1.972541529001468,
                targetPos: -6.297064289647577,
                targetVel: -1.5414200165198237,
                targetAcc: -1.5255311582232012,
                maxVel: 4.0,
                maxAcc: 9.973774779735683,
                maxJerk: 10.0,
                expectedDuration: 4.173396,
                expectedTimeIntervals: [0.982370, 0.000000, 0.785116, 1.238478, 0.507439, 0.000000, 0.659992]
            ),

            // Test Case 31: Positive to negative with complex state transitions and higher velocity limit
            (
                currentPos: 2.6603799559471364,
                currentVel: 1.9695301027900147,
                currentAcc: 1.972541529001468,
                targetPos: -6.297064289647577,
                targetVel: -1.5414200165198237,
                targetAcc: -1.5255311582232012,
                maxVel: 6.155966065528634,
                maxAcc: 9.973774779735683,
                maxJerk: 10.0,
                expectedDuration: 3.597270,
                expectedTimeIntervals: [1.109397, 0.000000, 0.912143, 0.047546, 0.687816, 0.000000, 0.840369]
            ),
        ]

    // MARK: - Tests

    func testTruthTable() throws {
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

            XCTAssertEqual(
                trajectory.getDuration(),
                testCase.expectedDuration,
                accuracy: durationAccuracy,
                "Test case \(index + 1): Duration mismatch"
            )

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
}
