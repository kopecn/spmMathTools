//
//  PositionSecondOrderStep2.swift
//
//
//  Created by Nicholas Bergantz on 4/25/24.
//

import Foundation

//! Mathematical equations for Step 2 in second-order position interface: Time synchronization
class PositionSecondOrderStep2 {

    var v0: Double
    var tf: Double
    var vf: Double
    var _vMax: Double
    var _vMin: Double
    var _aMax: Double
    var _aMin: Double

    // Pre-calculated expressions
    let pd: Double
    let vd: Double

    init(
        _ tf: Double,
        _ p0: Double,
        _ v0: Double,
        _ pf: Double,
        _ vf: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double
    ) {
        self.v0 = v0
        self.tf = tf
        self.vf = vf
        self._vMax = vMax
        self._vMin = vMin
        self._aMax = aMax
        self._aMin = aMin
        self.pd = pf - p0
        self.vd = vf - v0
    }

    /// Checks if the given profile can be executed within the specified velocity and acceleration limits using a time-none solution.
    ///
    /// This function first checks if the initial and final velocities and position are all close to zero. If so, it sets the profile times accordingly and checks if the profile is valid. Otherwise, it calculates a time-none solution and checks if the resulting acceleration is within the specified limits. If a valid profile is found, the function returns `true`, otherwise it returns `false`.
    ///
    /// - Parameters:
    ///   - profile: The profile to check.
    ///   - vMax: The maximum allowed velocity.
    ///   - vMin: The minimum allowed velocity.
    ///   - aMax: The maximum allowed acceleration.
    ///   - aMin: The minimum allowed acceleration.
    /// - Returns: `true` if the profile can be executed within the limits, `false` otherwise.
    func timeAcc0(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double
    ) -> Bool {
        // UD Solution 1/2

        let h1_1 = sqrt(
            (2 * aMax * (pd - tf * vf) - 2 * aMin * (pd - tf * v0) + vd * vd) / (aMax * aMin) + tf * tf
        )

        profile.t[0] = (aMax * vd - aMax * aMin * (tf - h1_1)) / (aMax * (aMax - aMin))
        profile.t[1] = h1_1
        profile.t[2] = tf - (profile.t[0] + h1_1)
        profile.t[3] = 0
        profile.t[4] = 0
        profile.t[5] = 0
        profile.t[6] = 0

        if profile.checkForSecondOrderWithTiming(tf, aMax, aMin, vMax, vMin, .UDDU, .ACC0) {
            profile.pf = profile.p.last!
            return true
        }

        // UU Solution

        let h1_2 = (-vd + aMax * tf)

        profile.t[0] = -vd * vd / (2 * aMax * h1_2) + (pd - v0 * tf) / h1_2
        profile.t[1] = -vd / aMax + tf
        profile.t[2] = 0
        profile.t[3] = 0
        profile.t[4] = 0
        profile.t[5] = 0
        profile.t[6] = tf - (profile.t[0] + profile.t[1])

        if profile.checkForSecondOrderWithTiming(tf, aMax, aMin, vMax, vMin, .UDDU, .ACC0) {
            profile.pf = profile.p.last!
            return true
        }

        // UU Solution - 2 step

        profile.t[0] = 0
        profile.t[1] = -vd / aMax + tf
        profile.t[2] = 0
        profile.t[3] = 0
        profile.t[4] = 0
        profile.t[5] = 0
        profile.t[6] = vd / aMax

        if profile.checkForSecondOrderWithTiming(tf, aMax, aMin, vMax, vMin, .UDDU, .ACC0) {
            profile.pf = profile.p.last!
            return true
        }

        return false
    }

    /// Checks if the given profile can be executed within the specified velocity and acceleration limits using a time-none solution.
    ///
    /// This function first checks if the initial and final velocities and position are all close to zero. If so, it sets the profile times accordingly and checks if the profile is valid. Otherwise, it calculates a time-none solution and checks if the resulting acceleration is within the specified limits. If a valid profile is found, the function returns `true`, otherwise it returns `false`.
    ///
    /// - Parameters:
    ///   - profile: The profile to check.
    ///   - vMax: The maximum allowed velocity.
    ///   - vMin: The minimum allowed velocity.
    ///   - aMax: The maximum allowed acceleration.
    ///   - aMin: The minimum allowed acceleration.
    /// - Returns: `true` if the profile can be executed within the limits, `false` otherwise.
    func timeNone(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double
    ) -> Bool {

        if abs(v0) < Double.leastNormalMagnitude && abs(vf) < Double.leastNormalMagnitude
            && abs(pd) < Double.leastNormalMagnitude
        {
            profile.t[0] = 0
            profile.t[1] = tf
            profile.t[2] = 0
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = 0

            if profile.checkForSecondOrderWithTiming(tf, aMax, aMin, vMax, vMin, .UDDU, .NONE) {
                profile.pf = profile.p.last!
                return true
            }
        }

        /// UD Solution 1/2

        let h1 = 2 * (vf * tf - pd)

        profile.t[0] = h1 / vd
        profile.t[1] = tf - profile.t[0]
        profile.t[2] = 0
        profile.t[3] = 0
        profile.t[4] = 0
        profile.t[5] = 0
        profile.t[6] = 0

        let af = vd * vd / h1

        if (aMin - 1e-12 < af) && (af < aMax + 1e-12)
            && profile.checkForSecondOrderWithTiming(tf, af, -af, vMax, vMin, .UDDU, .NONE)
        {
            profile.pf = profile.p.last!
            return true
        }

        return false
    }

    /// Checks if the given profile can be executed within the specified velocity and acceleration limits.
    ///
    /// This function tests all possible cases to find a profile that matches the given velocity and acceleration limits. It first checks the case where the desired position is positive, and then checks the case where the desired position is negative. The function returns `true` if a valid profile is found, and `false` otherwise.
    ///
    /// - Parameters:
    ///   - profile: The profile to check.
    ///   - vMax: The maximum allowed velocity.
    ///   - vMin: The minimum allowed velocity.
    ///   - aMax: The maximum allowed acceleration.
    ///   - aMin: The minimum allowed acceleration.
    /// - Returns: `true` if the profile can be executed within the limits, `false` otherwise.
    func getProfile(
        _ profile: inout Profile
    ) -> Bool {
        /// Test all cases to get ones that match
        /// However we should guess which one is correct and try them first...
        if pd > 0 {
            return checkAll(&profile, _vMax, _vMin, _aMax, _aMin)
                || checkAll(&profile, _vMin, _vMax, _aMin, _aMax)
        }

        return checkAll(&profile, _vMin, _vMax, _aMin, _aMax)
            || checkAll(&profile, _vMax, _vMin, _aMax, _aMin)
    }

    /// Checks if the given profile can be executed within the specified velocity and acceleration limits.
    ///
    /// - Parameters:
    ///   - profile: The profile to check.
    ///   - vMax: The maximum allowed velocity.
    ///   - vMin: The minimum allowed velocity.
    ///   - aMax: The maximum allowed acceleration.
    ///   - aMin: The minimum allowed acceleration.
    /// - Returns: `true` if the profile can be executed within the limits, `false` otherwise.
    func checkAll(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double
    ) -> Bool {
        timeAcc0(&profile, vMax, vMin, aMax, aMin)
            || timeNone(&profile, vMax, vMin, aMax, aMin)
    }
}
