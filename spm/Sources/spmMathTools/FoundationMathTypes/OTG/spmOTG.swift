//
//  OTG.swift
//
//
//  Created by Nicholas Bergantz on 3/23/24.
//

import Foundation

public class OTG {
    var currentInput: InputParameter
    var currentInputInitialized: Bool = false
    let calculator: TargetCalculator
    let maxNumberOfWaypoints: Int
    let degreesOfFreedom: Int
    var deltaTime: Double = 0.0

    public init(
        degreesOfFreedom: Int, 
        deltaTime: Double = -1.0, 
        maxNumberOfWaypoints: Int = 0
    ) {
        self.currentInput = InputParameter(DOFs: degreesOfFreedom)
        self.calculator = TargetCalculator(dofs: degreesOfFreedom)
        self.maxNumberOfWaypoints = maxNumberOfWaypoints
        self.degreesOfFreedom = degreesOfFreedom
        self.deltaTime = deltaTime
    }

    func reset() {
        currentInputInitialized = false
    }

    func filterIntermediatePositions(input: InputParameter, thresholdDistance: [Double]) -> [[Double]] {
        if input.intermediatePositions.isEmpty {
            return input.intermediatePositions
        }

        let nWaypoints = input.intermediatePositions.count
        var isActive = Array(repeating: true, count: nWaypoints)

        var start = 0
        var end = start + 2

        while end < nWaypoints + 2 {
            let posStart = (start == 0) ? input.currentPosition : input.intermediatePositions[start - 1]
            let posEnd = (end == nWaypoints + 1) ? input.targetPosition : input.intermediatePositions[end - 1]

            var areAllBelow = true

            for current in (start + 1)..<end {
                let posCurrent = input.intermediatePositions[current - 1]
                var tStartMax = 0.0
                var tEndMin = 1.0

                for dof in 0..<degreesOfFreedom {
                    let delta = posEnd[dof] - posStart[dof]
                    if delta == 0.0 { continue }
                    let h0 = (posCurrent[dof] - posStart[dof]) / delta
                    let tStart = h0 - thresholdDistance[dof] / abs(delta)
                    let tEnd = h0 + thresholdDistance[dof] / abs(delta)

                    tStartMax = max(tStart, tStartMax)
                    tEndMin = min(tEnd, tEndMin)

                    if tStartMax > tEndMin {
                        areAllBelow = false
                        break
                    }
                }

                if !areAllBelow {
                    break
                }
            }

            isActive[end - 2] = !areAllBelow
            if !areAllBelow {
                start = end - 1
            }

            end += 1
        }

        var filteredPositions: [[Double]] = []
        for i in 0..<nWaypoints {
            if isActive[i] {
                filteredPositions.append(input.intermediatePositions[i])
            }
        }

        return filteredPositions
    }

    func validateInput(
        input: InputParameter,
        checkCurrentStateWithinLimits: Bool = false,
        checkTargetStateWithinLimits: Bool = true,
    ) throws -> Bool {
        // Additional OTG-specific validations before calling InputParameter validation
        if !input.intermediatePositions.isEmpty && input.controlInterface == .Position {
            if input.intermediatePositions.count > maxNumberOfWaypoints {
                throw RuckigError(
                    "The number of intermediate positions \(input.intermediatePositions.count) exceeds the maximum number of waypoints \(maxNumberOfWaypoints)."
                )
            }
        }

        if deltaTime <= 0.0 && input.durationDiscretization != .Continuous {
            throw RuckigError("delta time (control rate) parameter \(deltaTime) should be larger than zero.")
        }

        return try input.validate(
            checkCurrentStateWithinLimits: checkCurrentStateWithinLimits,
            checkTargetStateWithinLimits: checkTargetStateWithinLimits
        )
    }

    public func calculate(
        input: InputParameter, 
        trajectory: inout Trajectory, 
    ) -> Result {
        do {
            let isValid = try validateInput(input: input)
            if !isValid {
                return .ErrorInvalidInput
            }
        } catch {
            return .ErrorInvalidInput
        }

        return calculator.calculate(
            input: input,
            trajectory: &trajectory,
            deltaTime: deltaTime,
        )
    }

    func update(
        input: InputParameter, 
        output: inout OutputParameter
    ) -> Result {
        let start = DispatchTime.now()

        output.newCalculation = false

        var result: Result = .Working
        if !currentInputInitialized || input != currentInput {
            result = calculate(
                input: input,
                trajectory: &output.trajectory,
            )
            if result != .Working && result != .ErrorPositionalLimits {
                return result
            }

            currentInput = input
            currentInputInitialized = true
            output.time = 0.0
            output.newCalculation = true
        }

        let oldSection = output.newSection
        output.time += deltaTime
        try! output.trajectory.atTime(
            output.time,
            newPosition: &output.newPosition,
            newVelocity: &output.newVelocity,
            newAcceleration: &output.newAcceleration,
            newJerk: &output.newJerk,
            newSection: &output.newSection
        )
        output.didSectionChange = (output.newSection > oldSection)

        let stop = DispatchTime.now()
        output.calculationDuration = Double(stop.uptimeNanoseconds - start.uptimeNanoseconds) / 1_000.0

        output.passToInput(&currentInput)

        if output.time > output.trajectory.getDuration() {
            return .Finished
        }

        return result
    }
}
