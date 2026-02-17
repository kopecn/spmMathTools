//
//  CalculatorTarget.swift
//
//
//  Created by Nicholas Bergantz on 4/25/24.
//

import Foundation

/// Calculation class for a state-to-state trajectory.
public class TargetCalculator {
    let degreesOfFreedom: Int

    private static let eps: Double = Double.ulpOfOne
    private static let returnErrorAtMaximalDuration = true

    private var newPhaseControl: [Double]
    private var pd: [Double]  // For phase synchronization
    private var possibleTSyncs: [Double]
    private var idx: [Int]

    private var blocks: [Block]
    private var inpMinVelocity: [Double]
    private var inpMinAcceleration: [Double]
    private var inpPerDofControlInterface: [ControlInterface]
    private var inpPerDofSynchronization: [Synchronization]

    public init(dofs: Int) {
        self.degreesOfFreedom = dofs
        self.blocks = (0..<dofs).map { _ in Block() }
        self.inpMinVelocity = Array(repeating: 0.0, count: dofs)
        self.inpMinAcceleration = Array(repeating: 0.0, count: dofs)
        self.inpPerDofControlInterface = Array(repeating: .Position, count: dofs)
        self.inpPerDofSynchronization = Array(repeating: .None, count: dofs)
        self.newPhaseControl = Array(repeating: 0.0, count: dofs)
        self.pd = Array(repeating: 0.0, count: dofs)
        self.possibleTSyncs = Array(repeating: Double.infinity, count: 3 * dofs + 1)
        self.idx = Array(repeating: 0, count: 3 * dofs + 1)
    }

    /// Is the trajectory (in principle) phase synchronizable?
    private func isInputCollinear(
        _ inp: InputParameter,
        _ limitingDirection: Direction,
        _ limitingDoF: Int
    ) -> Bool {
        // Check that vectors pd, v0, a0, vf, af are collinear
        for dof in 0..<degreesOfFreedom {
            pd[dof] = inp.targetPosition[dof] - inp.currentPosition[dof]
        }

        var scaleVector: [Double]? = nil
        var scaleDoF: Int? = nil  // Need to find a scale DOF because limiting DOF might not be phase synchronized

        for dof in 0..<degreesOfFreedom {
            if inpPerDofSynchronization[dof] != .Phase {
                continue
            }

            if inpPerDofControlInterface[dof] == .Position && abs(pd[dof]) > Self.eps {
                scaleVector = pd
                scaleDoF = dof
                break
            } else if abs(inp.currentVelocity[dof]) > Self.eps {
                scaleVector = inp.currentVelocity
                scaleDoF = dof
                break
            } else if abs(inp.currentAcceleration[dof]) > Self.eps {
                scaleVector = inp.currentAcceleration
                scaleDoF = dof
                break
            } else if abs(inp.targetVelocity[dof]) > Self.eps {
                scaleVector = inp.targetVelocity
                scaleDoF = dof
                break
            } else if abs(inp.targetAcceleration[dof]) > Self.eps {
                scaleVector = inp.targetAcceleration
                scaleDoF = dof
                break
            }
        }

        guard let scaleDoF = scaleDoF, let scaleVector = scaleVector else {
            return false  // Zero everywhere is in theory collinear, but that trivial case is better handled elsewhere
        }

        let scale = scaleVector[scaleDoF]
        let pdScale = pd[scaleDoF] / scale
        let v0Scale = inp.currentVelocity[scaleDoF] / scale
        let vfScale = inp.targetVelocity[scaleDoF] / scale
        let a0Scale = inp.currentAcceleration[scaleDoF] / scale
        let afScale = inp.targetAcceleration[scaleDoF] / scale

        let scaleLimiting = scaleVector[limitingDoF]
        var controlLimiting = (limitingDirection == .UP) ? inp.maxJerk[limitingDoF] : -inp.maxJerk[limitingDoF]
        if inp.maxJerk[limitingDoF].isInfinite {
            controlLimiting =
                (limitingDirection == .UP) ? inp.maxAcceleration[limitingDoF] : inpMinAcceleration[limitingDoF]
        }

        for dof in 0..<degreesOfFreedom {
            if inpPerDofSynchronization[dof] != .Phase {
                continue
            }

            let currentScale = scaleVector[dof]
            if (inpPerDofControlInterface[dof] == .Position && abs(pd[dof] - pdScale * currentScale) > Self.eps)
                || abs(inp.currentVelocity[dof] - v0Scale * currentScale) > Self.eps
                || abs(inp.currentAcceleration[dof] - a0Scale * currentScale) > Self.eps
                || abs(inp.targetVelocity[dof] - vfScale * currentScale) > Self.eps
                || abs(inp.targetAcceleration[dof] - afScale * currentScale) > Self.eps
            {
                return false
            }

            newPhaseControl[dof] = controlLimiting * currentScale / scaleLimiting
        }

        return true
    }

    private func synchronize(
        _ tMin: Double?,
        _ tSync: inout Double,
        _ limitingDoF: inout Int?,
        _ profiles: inout [Profile],
        _ discreteDuration: Bool,
        _ deltaTime: Double
    ) -> Bool {
        // Possible tSyncs are the start times of the intervals and optional tMin
        var anyInterval = false
        for dof in 0..<degreesOfFreedom {
            // Ignore DoFs without synchronization here
            if inpPerDofSynchronization[dof] == .None {
                possibleTSyncs[dof] = 0.0
                possibleTSyncs[degreesOfFreedom + dof] = Double.infinity
                possibleTSyncs[2 * degreesOfFreedom + dof] = Double.infinity
                continue
            }

            possibleTSyncs[dof] = blocks[dof].tMin
            possibleTSyncs[degreesOfFreedom + dof] = blocks[dof].a?.right ?? Double.infinity
            possibleTSyncs[2 * degreesOfFreedom + dof] = blocks[dof].b?.right ?? Double.infinity
            anyInterval = anyInterval || blocks[dof].a != nil || blocks[dof].b != nil
        }
        possibleTSyncs[3 * degreesOfFreedom] = tMin ?? Double.infinity
        anyInterval = anyInterval || tMin != nil

        if discreteDuration {
            for i in 0..<possibleTSyncs.count {
                if possibleTSyncs[i].isInfinite {
                    continue
                }

                let remainder = possibleTSyncs[i].truncatingRemainder(dividingBy: deltaTime)  // in [0, deltaTime)
                if remainder > Self.eps {
                    possibleTSyncs[i] += deltaTime - remainder
                }
            }
        }

        // Test them in sorted order
        let idxEnd = anyInterval ? idx.count : degreesOfFreedom
        for i in 0..<idxEnd {
            idx[i] = i
        }
        idx[0..<idxEnd].sort { i, j in possibleTSyncs[i] < possibleTSyncs[j] }

        // Start at last tmin (or worse)
        for i in (degreesOfFreedom - 1)..<idxEnd {
            let possibleTSync = possibleTSyncs[idx[i]]
            var isBlocked = false
            for dof in 0..<degreesOfFreedom {
                if inpPerDofSynchronization[dof] == .None {
                    continue  // inner dof loop
                }
                if blocks[dof].isBlocked(t: possibleTSync) {
                    isBlocked = true
                    break  // inner dof loop
                }
            }
            if isBlocked || possibleTSync < (tMin ?? 0.0) || possibleTSync.isInfinite {
                continue
            }

            tSync = possibleTSync
            if idx[i] == 3 * degreesOfFreedom {  // Optional tMin
                limitingDoF = nil
                return true
            }

            let divQuot = idx[i] / degreesOfFreedom
            let divRem = idx[i] % degreesOfFreedom
            limitingDoF = divRem
            switch divQuot {
            case 0:
                profiles[divRem] = blocks[divRem].pMin
            case 1:
                profiles[divRem] = blocks[divRem].a!.profile!
            case 2:
                profiles[divRem] = blocks[divRem].b!.profile!
            default:
                break
            }
            return true
        }

        print("[DEBUG] synchronize() returning false - no valid tSync found")
        print("[DEBUG] anyInterval=\(anyInterval), idxEnd=\(idxEnd), degreesOfFreedom=\(degreesOfFreedom)")
        return false
    }

    /// Calculate the time-optimal waypoint-based trajectory
    public func calculate(
        input: InputParameter,
        trajectory: inout Trajectory,
        deltaTime: Double,
    ) -> Result {

        // Check for trivial case: all positions equal and all velocities/accelerations zero
        var allPositionsEqual = true
        var allVelocitiesZero = true
        var allAccelerationsZero = true

        for dof in 0..<degreesOfFreedom {
            if abs(input.targetPosition[dof] - input.currentPosition[dof]) > Self.eps {
                allPositionsEqual = false
            }
            if abs(input.currentVelocity[dof]) > Self.eps || abs(input.targetVelocity[dof]) > Self.eps {
                allVelocitiesZero = false
            }
            if abs(input.currentAcceleration[dof]) > Self.eps || abs(input.targetAcceleration[dof]) > Self.eps {
                allAccelerationsZero = false
            }
        }

        if allPositionsEqual && allVelocitiesZero && allAccelerationsZero {
            // Trivial case: already at target with zero motion - initialize trajectory properly
            trajectory.duration = 0.0
            trajectory.cumulativeTimes[0] = 0.0

            for dof in 0..<degreesOfFreedom {
                // Initialize all profile arrays to represent the current state with zero duration
                var p = trajectory.profiles[0][dof]

                // Clear existing profile data
                p.t = Array(repeating: 0.0, count: 7)
                p.tSum = Array(repeating: 0.0, count: 7)
                p.j = Array(repeating: 0.0, count: 7)
                p.a = Array(repeating: input.currentAcceleration[dof], count: 7)
                p.v = Array(repeating: input.currentVelocity[dof], count: 7)
                p.p = Array(repeating: input.currentPosition[dof], count: 7)

                trajectory.profiles[0][dof] = p
                trajectory.independentMinDurations[dof] = 0.0
            }

            return .Working
        }

        for dof in 0..<degreesOfFreedom {
            var p = trajectory.profiles[0][dof]

            inpMinVelocity[dof] = input.minVelocity?[dof] ?? -input.maxVelocity[dof]
            inpMinAcceleration[dof] = input.minAcceleration?[dof] ?? -input.maxAcceleration[dof]
            inpPerDofControlInterface[dof] = input.perDofControlInterface?[dof] ?? input.controlInterface
            inpPerDofSynchronization[dof] = input.perDofSynchronization?[dof] ?? input.synchronization

            if !input.enabled[dof] {
                p.p[p.p.count - 1] = input.currentPosition[dof]
                p.v[p.v.count - 1] = input.currentVelocity[dof]
                p.a[p.a.count - 1] = input.currentAcceleration[dof]
                p.tSum[p.tSum.count - 1] = 0.0
                blocks[dof].tMin = 0.0
                blocks[dof].a = nil
                blocks[dof].b = nil
                trajectory.profiles[0][dof] = p
                continue
            }

            // Calculate brake (if input exceeds or will exceed limits)
            switch inpPerDofControlInterface[dof] {
            case .Position:
                if !input.maxJerk[dof].isInfinite {
                    p.brake.getPositionBrakeTrajectory(
                        input.currentVelocity[dof],
                        input.currentAcceleration[dof],
                        input.maxVelocity[dof],
                        inpMinVelocity[dof],
                        input.maxAcceleration[dof],
                        inpMinAcceleration[dof],
                        input.maxJerk[dof]
                    )
                } else if !input.maxAcceleration[dof].isInfinite {
                    p.brake.getSecondOrderPositionBrakeTrajectory(
                        input.currentVelocity[dof],
                        input.maxVelocity[dof],
                        inpMinVelocity[dof],
                        input.maxAcceleration[dof],
                        inpMinAcceleration[dof]
                    )
                }
                p.setBoundary(
                    input.currentPosition[dof],
                    input.currentVelocity[dof],
                    input.currentAcceleration[dof],
                    input.targetPosition[dof],
                    input.targetVelocity[dof],
                    input.targetAcceleration[dof]
                )
            case .Velocity:
                if !input.maxJerk[dof].isInfinite {
                    p.brake.getVelocityBrakeTrajectory(
                        input.currentAcceleration[dof],
                        input.maxAcceleration[dof],
                        inpMinAcceleration[dof],
                        input.maxJerk[dof]
                    )
                } else {
                    p.brake.getSecondOrderVelocityBrakeTrajectory()
                }
                p.setBoundaryForVelocity(
                    input.currentPosition[dof],
                    input.currentVelocity[dof],
                    input.currentAcceleration[dof],
                    input.targetVelocity[dof],
                    input.targetAcceleration[dof]
                )
            }

            // Finalize pre & post-trajectories
            // Extract to temp vars to avoid simultaneous access to p
            if !input.maxJerk[dof].isInfinite {
                var p0 = p.p[0]; var v0 = p.v[0]; var a0 = p.a[0]
                p.brake.finalize(&p0, &v0, &a0)
                p.p[0] = p0; p.v[0] = v0; p.a[0] = a0
            } else if !input.maxAcceleration[dof].isInfinite {
                var p0 = p.p[0]; var v0 = p.v[0]; var a0 = p.a[0]
                p.brake.finalizeSecondOrder(&p0, &v0, &a0)
                p.p[0] = p0; p.v[0] = v0; p.a[0] = a0
            }

            var foundProfile = false
            switch inpPerDofControlInterface[dof] {
            case .Position:
                if !input.maxJerk[dof].isInfinite {
                    // MARK: - P-ThirdOrderStep1
                    let step1 = PositionThirdOrderStep1(
                        p0: p.p[0],
                        v0: p.v[0],
                        a0: p.a[0],
                        pf: p.pf,
                        vf: p.vf,
                        af: p.af,
                        vMax: input.maxVelocity[dof],
                        vMin: inpMinVelocity[dof],
                        aMax: input.maxAcceleration[dof],
                        aMin: inpMinAcceleration[dof],
                        jMax: input.maxJerk[dof]
                    )
                    foundProfile = step1.getProfile(&p, &blocks[dof])
                } else if !input.maxAcceleration[dof].isInfinite {
                    // MARK: - P-SecondOrderStep1
                    let step1 = PositionSecondOrderStep1(
                        p0: p.p[0],
                        v0: p.v[0],
                        pf: p.pf,
                        vf: p.vf,
                        vMax: input.maxVelocity[dof],
                        vMin: inpMinVelocity[dof],
                        aMax: input.maxAcceleration[dof],
                        aMin: inpMinAcceleration[dof]
                    )
                    foundProfile = step1.getProfile(input: &p, block: &blocks[dof])
                } else {
                    // MARK: - P-FirstOrderStep1
                    let step1 = PositionFirstOrderStep1(
                        p0: p.p[0],
                        pf: p.pf,
                        vMax: input.maxVelocity[dof],
                        vMin: inpMinVelocity[dof]
                    )
                    foundProfile = step1.getProfile(input: &p, block: &blocks[dof])
                }
            case .Velocity:
                if !input.maxJerk[dof].isInfinite {
                    // MARK: - V-ThirdOrderStep1
                    let step1 = VelocityThirdOrderStep1(
                        v0: p.v[0],
                        a0: p.a[0],
                        vf: p.vf,
                        af: p.af,
                        aMax: input.maxAcceleration[dof],
                        aMin: inpMinAcceleration[dof],
                        jMax: input.maxJerk[dof]
                    )
                    foundProfile = step1.getProfile(&p, &blocks[dof])
                } else {
                    // MARK: - V-SecondOrderStep1
                    let step1 = VelocitySecondOrderStep1(
                        v0: p.v[0],
                        vf: p.vf,
                        aMax: input.maxAcceleration[dof],
                        aMin: inpMinAcceleration[dof]
                    )
                    foundProfile = step1.getProfile(&p, &blocks[dof])
                }
            }

            if !foundProfile {
                let hasZeroLimits =
                    (input.maxAcceleration[dof] == 0.0 || inpMinAcceleration[dof] == 0.0
                        || input.maxJerk[dof] == 0.0)
                if hasZeroLimits {
                    return .ErrorZeroLimits
                } else {
                    return .ErrorExecutionTimeCalculation
                }
            }

            trajectory.independentMinDurations[dof] = blocks[dof].tMin
            trajectory.profiles[0][dof] = p
        }

        let discreteDuration = (input.durationDiscretization == .Discrete)
        if degreesOfFreedom == 1 && input.minimumDuration == nil && !discreteDuration {
            trajectory.duration = blocks[0].tMin
            trajectory.profiles[0][0] = blocks[0].pMin
            trajectory.cumulativeTimes[0] = trajectory.duration
            return .Working
        }

        var limitingDoF: Int? = nil  // The DoF that doesn't need step 2
        let foundSynchronization = synchronize(
            input.minimumDuration,
            &trajectory.duration,
            &limitingDoF,
            &trajectory.profiles[0],
            discreteDuration,
            deltaTime
        )
        if !foundSynchronization {
            print("[DEBUG] Synchronization FAILED for controlInterface: \(input.controlInterface)")
            var hasZeroLimits = false
            for dof in 0..<degreesOfFreedom {
                if input.maxAcceleration[dof] == 0.0 || inpMinAcceleration[dof] == 0.0 || input.maxJerk[dof] == 0.0
                {
                    hasZeroLimits = true
                    break
                }
            }

            if hasZeroLimits {
                return .ErrorZeroLimits
            } else {
                return .ErrorSynchronizationCalculation
            }
        }

        // None Synchronization
        for dof in 0..<degreesOfFreedom {
            if input.enabled[dof] && inpPerDofSynchronization[dof] == .None {
                trajectory.profiles[0][dof] = blocks[dof].pMin
                if blocks[dof].tMin > trajectory.duration {
                    trajectory.duration = blocks[dof].tMin
                    limitingDoF = dof
                }
            }
        }
        trajectory.cumulativeTimes[0] = trajectory.duration

        if Self.returnErrorAtMaximalDuration {
            if trajectory.duration > 7.6e3 {
                return .ErrorTrajectoryDuration
            }
        }

        if trajectory.duration == 0.0 {
            // Copy all profiles for end state
            for dof in 0..<degreesOfFreedom {
                trajectory.profiles[0][dof] = blocks[dof].pMin
            }
            return .Working
        }

        if !discreteDuration && inpPerDofSynchronization.allSatisfy({ $0 == .None }) {
            return .Working
        }

        // Phase Synchronization
        if let limitingDoF = limitingDoF, inpPerDofSynchronization.contains(.Phase) {
            let pLimiting = trajectory.profiles[0][limitingDoF]
            if isInputCollinear(input, pLimiting.direction, limitingDoF) {
                var foundTimeSynchronization = true
                for dof in 0..<degreesOfFreedom {
                    if !input.enabled[dof] || dof == limitingDoF || inpPerDofSynchronization[dof] != .Phase {
                        continue
                    }

                    var p = trajectory.profiles[0][dof]
                    let tProfile = trajectory.duration - p.brake.duration - p.accel.duration

                    p.t = pLimiting.t  // Copy timing information from limiting DoF
                    p.controlSigns = pLimiting.controlSigns

                    // Profile.ReachedLimits.NONE is a small hack, as there is no specialization for that in the check function
                    switch inpPerDofControlInterface[dof] {
                    case .Position:
                        switch p.controlSigns {
                        case .UDDU:
                            if !input.maxJerk[dof].isInfinite {
                                foundTimeSynchronization =
                                    foundTimeSynchronization
                                    && p.checkWithTiming(
                                        tProfile,
                                        newPhaseControl[dof],
                                        input.maxVelocity[dof],
                                        inpMinVelocity[dof],
                                        input.maxAcceleration[dof],
                                        inpMinAcceleration[dof],
                                        input.maxJerk[dof],
                                        .UDDU,
                                        .NONE
                                    )
                            } else if !input.maxAcceleration[dof].isInfinite {
                                foundTimeSynchronization =
                                    foundTimeSynchronization
                                    && p.checkForSecondOrderWithTiming(
                                        tProfile,
                                        newPhaseControl[dof],
                                        -newPhaseControl[dof],
                                        input.maxVelocity[dof],
                                        inpMinVelocity[dof],
                                        input.maxAcceleration[dof],
                                        inpMinAcceleration[dof],
                                        .UDDU,
                                        .NONE
                                    )
                            } else {
                                foundTimeSynchronization =
                                    foundTimeSynchronization
                                    && p.checkForFirstOrderWithTiming(
                                        tProfile,
                                        newPhaseControl[dof],
                                        input.maxVelocity[dof],
                                        inpMinVelocity[dof],
                                        .UDDU,
                                        .NONE
                                    )
                            }
                        case .UDUD:
                            if !input.maxJerk[dof].isInfinite {
                                foundTimeSynchronization =
                                    foundTimeSynchronization
                                    && p.checkWithTiming(
                                        tProfile,
                                        newPhaseControl[dof],
                                        input.maxVelocity[dof],
                                        inpMinVelocity[dof],
                                        input.maxAcceleration[dof],
                                        inpMinAcceleration[dof],
                                        input.maxJerk[dof],
                                        .UDUD,
                                        .NONE
                                    )
                            } else {
                                foundTimeSynchronization =
                                    foundTimeSynchronization
                                    && p.checkForSecondOrderWithTiming(
                                        tProfile,
                                        newPhaseControl[dof],
                                        -newPhaseControl[dof],
                                        input.maxVelocity[dof],
                                        inpMinVelocity[dof],
                                        input.maxAcceleration[dof],
                                        inpMinAcceleration[dof],
                                        .UDUD,
                                        .NONE
                                    )
                            }
                        }
                    case .Velocity:
                        switch p.controlSigns {
                        case .UDDU:
                            if !input.maxJerk[dof].isInfinite {
                                foundTimeSynchronization =
                                    foundTimeSynchronization
                                    && p.checkForVelocityWithTiming(
                                        tProfile,
                                        newPhaseControl[dof],
                                        input.maxAcceleration[dof],
                                        inpMinAcceleration[dof],
                                        input.maxJerk[dof],
                                        .UDDU,
                                        .NONE
                                    )
                            } else {
                                foundTimeSynchronization =
                                    foundTimeSynchronization
                                    && p.checkForSecondOrderVelocityWithTiming(
                                        tProfile,
                                        newPhaseControl[dof],
                                        input.maxAcceleration[dof],
                                        inpMinAcceleration[dof],
                                        .UDDU,
                                        .NONE
                                    )
                            }
                        case .UDUD:
                            if !input.maxJerk[dof].isInfinite {
                                foundTimeSynchronization =
                                    foundTimeSynchronization
                                    && p.checkForVelocityWithTiming(
                                        tProfile,
                                        newPhaseControl[dof],
                                        input.maxAcceleration[dof],
                                        inpMinAcceleration[dof],
                                        input.maxJerk[dof],
                                        .UDUD,
                                        .NONE
                                    )
                            } else {
                                foundTimeSynchronization =
                                    foundTimeSynchronization
                                    && p.checkForSecondOrderVelocityWithTiming(
                                        tProfile,
                                        newPhaseControl[dof],
                                        input.maxAcceleration[dof],
                                        inpMinAcceleration[dof],
                                        .UDUD,
                                        .NONE
                                    )
                            }
                        }
                    }

                    p.limits = pLimiting.limits  // After check method call to set correct limits
                    trajectory.profiles[0][dof] = p
                }

                if foundTimeSynchronization && inpPerDofSynchronization.allSatisfy({ $0 == .Phase || $0 == .None })
                {
                    return .Working
                }
            }
        }

        // Time Synchronization
        for dof in 0..<degreesOfFreedom {
            let skipSynchronization =
                (dof == limitingDoF || inpPerDofSynchronization[dof] == .None) && !discreteDuration
            if !input.enabled[dof] || skipSynchronization {
                continue
            }

            var p = trajectory.profiles[0][dof]
            let tProfile = trajectory.duration - p.brake.duration - p.accel.duration

            if inpPerDofSynchronization[dof] == .TimeIfNecessary && abs(input.targetVelocity[dof]) < Self.eps
                && abs(input.targetAcceleration[dof]) < Self.eps
            {
                p = blocks[dof].pMin
                trajectory.profiles[0][dof] = p
                continue
            }

            // Check if the final time corresponds to an extremal profile calculated in step 1
            // Use 2*eps because of numerical robustness in duration discretization
            if abs(tProfile - blocks[dof].tMin) < 2 * Self.eps {
                p = blocks[dof].pMin
                trajectory.profiles[0][dof] = p
                continue
            } else if let a = blocks[dof].a, abs(tProfile - a.right) < 2 * Self.eps {
                p = a.profile!
                trajectory.profiles[0][dof] = p
                continue
            } else if let b = blocks[dof].b, abs(tProfile - b.right) < 2 * Self.eps {
                p = b.profile!
                trajectory.profiles[0][dof] = p
                continue
            }

            var foundTimeSynchronization = false
            switch inpPerDofControlInterface[dof] {
            case .Position:
                if !input.maxJerk[dof].isInfinite {

                    // MARK: - P-ThirdOrderStep2
                    let step2 = PositionThirdOrderStep2(
                        tProfile,
                        p.p[0],
                        p.v[0],
                        p.a[0],
                        p.pf,
                        p.vf,
                        p.af,
                        input.maxVelocity[dof],
                        inpMinVelocity[dof],
                        input.maxAcceleration[dof],
                        inpMinAcceleration[dof],
                        input.maxJerk[dof]
                    )
                    foundTimeSynchronization = step2.getProfile(&p)
                } else if !input.maxAcceleration[dof].isInfinite {
                    // MARK: - P-SecondOrderStep2
                    let step2 = PositionSecondOrderStep2(
                        tProfile,
                        p.p[0],
                        p.v[0],
                        p.pf,
                        p.vf,
                        input.maxVelocity[dof],
                        inpMinVelocity[dof],
                        input.maxAcceleration[dof],
                        inpMinAcceleration[dof]
                    )
                    foundTimeSynchronization = step2.getProfile(&p)
                } else {
                    // MARK: - P-FirstOrderStep2
                    let step2 = PositionFirstOrderStep2(
                        tf: tProfile,
                        p0: p.p[0],
                        pf: p.pf,
                        vMax: input.maxVelocity[dof],
                        vMin: inpMinVelocity[dof]
                    )
                    foundTimeSynchronization = step2.getProfile(profile: &p)
                }
            case .Velocity:
                if !input.maxJerk[dof].isInfinite {
                    // MARK: - V-ThirdOrderStep2
                    let step2 = VelocityThirdOrderStep2(
                        tf: tProfile,
                        v0: p.v[0],
                        a0: p.a[0],
                        vf: p.vf,
                        af: p.af,
                        aMax: input.maxAcceleration[dof],
                        aMin: inpMinAcceleration[dof],
                        jMax: input.maxJerk[dof]
                    )
                    foundTimeSynchronization = step2.getProfile(&p)
                } else {
                    // MARK: - V-SecondOrderStep2
                    let step2 = VelocitySecondOrderStep2(
                        tf: tProfile,
                        v0: p.v[0],
                        vf: p.vf,
                        aMax: input.maxAcceleration[dof],
                        aMin: inpMinAcceleration[dof]
                    )
                    foundTimeSynchronization = step2.getProfile(&p)
                }
            }

            if !foundTimeSynchronization {
                print("[DEBUG] ErrorSynchronizationCalculation at DOF \(dof)")
                print("  controlInterface: \(inpPerDofControlInterface[dof])")
                print("  tProfile: \(tProfile)")
                print("  p.v[0]=\(p.v[0]), p.vf=\(p.vf)")
                return .ErrorSynchronizationCalculation
            }

            trajectory.profiles[0][dof] = p
        }

        return .Working
    }
}
