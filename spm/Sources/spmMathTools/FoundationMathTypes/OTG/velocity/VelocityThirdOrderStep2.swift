//
//  VelocityThirdOrderStep2.swift
//
//
//  Created by Nicholas Bergantz on 4/25/24.
//

import Foundation

/// Mathematical equations for Step 2 in third-order velocity interface: Time synchronization
class VelocityThirdOrderStep2 {
    let a0: Double
    let tf: Double
    let af: Double
    let _aMax: Double
    let _aMin: Double
    let _jMax: Double

    /// Pre-calculated expressions
    let vd: Double
    let ad: Double

    init(tf: Double, v0: Double, a0: Double, vf: Double, af: Double, aMax: Double, aMin: Double, jMax: Double) {
        self.tf = tf
        self.a0 = a0
        self.af = af
        self._aMax = aMax
        self._aMin = aMin
        self._jMax = jMax
        self.vd = vf - v0
        self.ad = af - a0
    }

    private func timeAcc0(_ profile: inout Profile, _ aMax: Double, _ aMin: Double, _ jMax: Double) -> Bool {
        // UD Solution 1/2
        do {
            let h1 = sqrt((-ad * ad + 2 * jMax * ((a0 + af) * tf - 2 * vd)) / (jMax * jMax) + tf * tf)

            profile.t[0] = ad / (2 * jMax) + (tf - h1) / 2
            profile.t[1] = h1
            profile.t[2] = tf - (profile.t[0] + h1)
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = 0

            if profile.checkForVelocityWithTiming(tf, jMax, aMax, aMin, .UDDU, .ACC0) {
                profile.pf = profile.p[7]
                return true
            }
        }

        // UU Solution
        do {
            let h1 = (-ad + jMax * tf)

            profile.t[0] = -ad * ad / (2 * jMax * h1) + (vd - a0 * tf) / h1
            profile.t[1] = -ad / jMax + tf
            profile.t[2] = 0
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = tf - (profile.t[0] + profile.t[1])

            if profile.checkForVelocityWithTiming(tf, jMax, aMax, aMin, .UDDU, .ACC0) {
                profile.pf = profile.p[7]
                return true
            }
        }

        // UU Solution - 2 step
        do {
            profile.t[0] = 0
            profile.t[1] = -ad / jMax + tf
            profile.t[2] = 0
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = ad / jMax

            if profile.checkForVelocityWithTiming(tf, jMax, aMax, aMin, .UDDU, .ACC0) {
                profile.pf = profile.p[7]
                return true
            }
        }

        return false
    }

    private func timeNone(_ profile: inout Profile, _ aMax: Double, _ aMin: Double, _ jMax: Double) -> Bool {
        // Special case: all zero
        if abs(a0) < Double.ulpOfOne && abs(af) < Double.ulpOfOne && abs(vd) < Double.ulpOfOne {
            profile.t[0] = 0
            profile.t[1] = tf
            profile.t[2] = 0
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = 0

            if profile.checkForVelocityWithTiming(tf, jMax, aMax, aMin, .UDDU, .NONE) {
                profile.pf = profile.p[7]
                return true
            }
        }

        // UD Solution 1/2
        do {
            let h1 = 2 * (af * tf - vd)

            profile.t[0] = h1 / ad
            profile.t[1] = tf - profile.t[0]
            profile.t[2] = 0
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = 0

            let jf = ad * ad / h1

            if abs(jf) < abs(jMax) + 1e-12 && profile.checkForVelocityWithTiming(tf, jf, aMax, aMin, .UDDU, .NONE) {
                profile.pf = profile.p[7]
                return true
            }
        }

        return false
    }

    private func checkAll(_ profile: inout Profile, _ aMax: Double, _ aMin: Double, _ jMax: Double) -> Bool {
        timeAcc0(&profile, aMax, aMin, jMax) || timeNone(&profile, aMax, aMin, jMax)
    }

    func getProfile(_ profile: inout Profile) -> Bool {
        if vd > 0 {
            return checkAll(&profile, _aMax, _aMin, _jMax) || checkAll(&profile, _aMin, _aMax, -_jMax)
        }

        return checkAll(&profile, _aMin, _aMax, -_jMax) || checkAll(&profile, _aMax, _aMin, _jMax)
    }
}
