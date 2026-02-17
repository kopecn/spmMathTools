//
//  PositionSecondOrderStep1.swift
//
//
//  Created by Nicholas Bergantz on 4/25/24.
//

import Foundation

/// Mathematical equations for Step 1 in second-order position interface: Extremal profiles
class PositionSecondOrderStep1 {

    var v0: Double
    var vf: Double
    var _vMax: Double
    var _vMin: Double
    var _aMax: Double
    var _aMin: Double

    /// Pre-calculated expressions
    let pd: Double

    // Max 3 valid profiles
    // CRITICAL: Create separate Profile instances, not shared references
    // Using `repeating:` with a class creates multiple references to the SAME object!
    var validProfiles: [Profile] = (0..<3).map { _ in Profile() }
    private var profileCount: Int = 0

    init(
        p0: Double,
        v0: Double,
        pf: Double,
        vf: Double,
        vMax: Double,
        vMin: Double,
        aMax: Double,
        aMin: Double
    ) {
        self.v0 = v0
        self.vf = vf
        self._vMax = vMax
        self._vMin = vMin
        self._aMax = aMax
        self._aMin = aMin
        self.pd = pf - p0
    }

    private func addProfile() {
        profileCount += 1
        if profileCount < validProfiles.count {
            validProfiles[profileCount].setBoundary(validProfiles[profileCount - 1])
        }
    }

    private func resetProfiles() {
        profileCount = 0
    }

    private func hasProfiles() -> Bool {
        profileCount > 0
    }

    /// Generates a valid motion profile for a given input and block.
    ///
    /// This function handles the special case where the maximum and minimum velocities are both zero, as well as the case where the final velocity is zero. It iterates through different combinations of maximum and minimum velocities and accelerations to find a valid motion profile.
    ///
    /// - Parameters:
    ///   - input: The input profile to be set.
    ///   - block: The block to be calculated.
    /// - Returns: `true` if a valid motion profile is found, `false` otherwise.
    func getProfile(input: inout Profile, block: inout Block) -> Bool {
        /// Zero-limits special case
        if _vMax == 0.0 && _vMin == 0.0 {
            var p = block.pMin
            p.setBoundary(input)

            if timeAllSingleStep(&p, _vMax, _vMin, _aMax, _aMin) {
                block.pMin = p
                block.tMin = p.tSum.last! + p.brake.duration + p.accel.duration
                if abs(v0) > Double.leastNormalMagnitude {
                    block.a = Interval(block.tMin, Double.infinity)
                }
                return true
            }
            return false
        }

        resetProfiles()
        validProfiles[profileCount].setBoundary(input)

        if abs(vf) < Double.leastNormalMagnitude {
            // There is no blocked interval when vf==0, so return after first found profile
            let vMax = (pd >= 0) ? _vMax : _vMin
            let vMin = (pd >= 0) ? _vMin : _vMax
            let aMax = (pd >= 0) ? _aMax : _aMin
            let aMin = (pd >= 0) ? _aMin : _aMax

            timeNone(vMax, vMin, aMax, aMin, true)
            if hasProfiles() {
                var count = profileCount
                return Block.calculateBlock(&block, &validProfiles, &count)
            }
            timeAcc0(vMax, vMin, aMax, aMin, true)
            if hasProfiles() {
                var count = profileCount
                return Block.calculateBlock(&block, &validProfiles, &count)
            }

            timeNone(vMin, vMax, aMin, aMax, true)
            if hasProfiles() {
                var count = profileCount
                return Block.calculateBlock(&block, &validProfiles, &count)
            }
            timeAcc0(vMin, vMax, aMin, aMax, true)

        } else {
            timeNone(_vMax, _vMin, _aMax, _aMin, false)
            timeNone(_vMin, _vMax, _aMin, _aMax, false)
            timeAcc0(_vMax, _vMin, _aMax, _aMin, false)
            timeAcc0(_vMin, _vMax, _aMin, _aMax, false)
        }

        var count = profileCount
        return Block.calculateBlock(&block, &validProfiles, &count)
    }

    /// Attempts to find a valid time profile for a single step motion, given the maximum and minimum velocity and acceleration constraints.
    ///
    /// - Parameters:
    ///   - vMax: The maximum allowed velocity.
    ///   - vMin: The minimum allowed velocity.
    ///   - aMax: The maximum allowed acceleration.
    ///   - aMin: The minimum allowed acceleration.
    ///   - returnAfterFound: If true, the function will return after the first valid profile is found.
    ///
    /// This function generates a single time profile for the given constraints and adds it to the `validProfiles` array if it is valid. The function returns once the first valid profile is found.
    func timeAcc0(
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ returnAfterFound: Bool
    ) {
        validProfiles[profileCount].t[0] = (-v0 + vMax) / aMax
        validProfiles[profileCount].t[1] =
            (aMin * v0 * v0 - aMax * vf * vf) / (2 * aMax * aMin * vMax) + vMax * (aMax - aMin)
            / (2 * aMax * aMin) + pd / vMax
        validProfiles[profileCount].t[2] = (vf - vMax) / aMin
        validProfiles[profileCount].t[3] = 0
        validProfiles[profileCount].t[4] = 0
        validProfiles[profileCount].t[5] = 0
        validProfiles[profileCount].t[6] = 0

        if validProfiles[profileCount]
            .checkForSecondOrder(aMax, aMin, vMax, vMin, .UDDU, .ACC0)
        {
            addProfile()
        }
    }

    /// Attempts to find a valid time profile for a single step motion, given the maximum and minimum velocity and acceleration constraints.
    ///
    /// - Parameters:
    ///   - vMax: The maximum allowed velocity.
    ///   - vMin: The minimum allowed velocity.
    ///   - aMax: The maximum allowed acceleration.
    ///   - aMin: The minimum allowed acceleration.
    ///   - returnAfterFound: If true, the function will return after the first valid profile is found.
    ///
    /// This function generates two possible time profiles for the given constraints and adds them to the `validProfiles` array if they are valid. The function returns once the first valid profile is found, or after both profiles have been checked.
    func timeNone(
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ returnAfterFound: Bool
    ) {
        var h1 = (aMax * vf * vf - aMin * v0 * v0 - 2 * aMax * aMin * pd) / (aMax - aMin)
        if h1 >= 0.0 {
            h1 = sqrt(h1)

            // Solution 1
            validProfiles[profileCount].t[0] = -(v0 + h1) / aMax
            validProfiles[profileCount].t[1] = 0
            validProfiles[profileCount].t[2] = (vf + h1) / aMin
            validProfiles[profileCount].t[3] = 0
            validProfiles[profileCount].t[4] = 0
            validProfiles[profileCount].t[5] = 0
            validProfiles[profileCount].t[6] = 0

            if validProfiles[profileCount]
                .checkForSecondOrder(aMax, aMin, vMax, vMin, .UDDU, .NONE)
            {
                addProfile()
                if returnAfterFound { return }
            }

            // Solution 2
            validProfiles[profileCount].t[0] = (-v0 + h1) / aMax
            validProfiles[profileCount].t[1] = 0
            validProfiles[profileCount].t[2] = (vf - h1) / aMin
            validProfiles[profileCount].t[3] = 0
            validProfiles[profileCount].t[4] = 0
            validProfiles[profileCount].t[5] = 0
            validProfiles[profileCount].t[6] = 0

            if validProfiles[profileCount]
                .checkForSecondOrder(aMax, aMin, vMax, vMin, .UDDU, .NONE)
            {
                addProfile()
            }
        }
    }

    /// Attempts to time a single step of a motion profile, given the maximum and minimum velocity and acceleration constraints.
    ///
    /// - Parameters:
    ///   - profile: An inout `Profile` object that will be updated with the timing information.
    ///   - vMax: The maximum allowed velocity.
    ///   - vMin: The minimum allowed velocity.
    ///   - aMax: The maximum allowed acceleration.
    ///   - aMin: The minimum allowed acceleration.
    /// - Returns: `true` if the single step can be timed within the constraints, `false` otherwise.
    func timeAllSingleStep(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double
    ) -> Bool {
        if abs(vf - v0) > Double.leastNormalMagnitude {
            return false
        }

        profile.t[0] = 0
        profile.t[1] = 0
        profile.t[2] = 0
        profile.t[3] = 0
        profile.t[4] = 0
        profile.t[5] = 0
        profile.t[6] = 0

        if abs(v0) > Double.leastNormalMagnitude {
            profile.t[3] = pd / v0
            if profile.checkForSecondOrder(0.0, 0.0, vMax, vMin, .UDDU, .NONE) {
                return true
            }

        } else if abs(pd) < Double.leastNormalMagnitude {
            if profile.checkForSecondOrder(0.0, 0.0, vMax, vMin, .UDDU, .NONE) {
                return true
            }
        }

        return false
    }
}
