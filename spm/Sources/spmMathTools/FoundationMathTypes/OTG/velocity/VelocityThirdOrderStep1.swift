//
//  VelocityThirdOrderStep1.swift
//
//
//  Created by Nicholas Bergantz on 4/25/24.
//

import Foundation

/// Mathematical equations for Step 1 in third-order velocity interface: Extremal profiles
class VelocityThirdOrderStep1 {
    let a0: Double
    let af: Double
    let _aMax: Double
    let _aMin: Double
    let _jMax: Double

    /// Pre-calculated expressions
    let vd: Double

    // Max 3 valid profiles
    var validProfiles: [Profile] = (0..<3).map { _ in Profile() }
    private var profileCount: Int = 0

    init(v0: Double, a0: Double, vf: Double, af: Double, aMax: Double, aMin: Double, jMax: Double) {
        self.a0 = a0
        self.af = af
        self._aMax = aMax
        self._aMin = aMin
        self._jMax = jMax
        self.vd = vf - v0
    }

    private func addProfile() {
        profileCount += 1
        if profileCount < validProfiles.count {
            validProfiles[profileCount].setBoundary(&validProfiles[profileCount - 1])
        }
    }

    private func resetProfiles() {
        profileCount = 0
    }

    private func hasProfiles() -> Bool {
        profileCount > 0
    }

    private func timeAcc0(_ aMax: Double, _ aMin: Double, _ jMax: Double, _ returnAfterFound: Bool) {
        let profile = validProfiles[profileCount]

        profile.t[0] = (-a0 + aMax) / jMax
        profile.t[1] = (a0 * a0 + af * af) / (2 * aMax * jMax) - aMax / jMax + vd / aMax
        profile.t[2] = (-af + aMax) / jMax
        profile.t[3] = 0
        profile.t[4] = 0
        profile.t[5] = 0
        profile.t[6] = 0

        if profile.checkForVelocity(jMax, aMax, aMin, .UDDU, .ACC0) {
            validProfiles[profileCount] = profile
            addProfile()
            if returnAfterFound { return }
        }
    }

    private func timeNone(_ aMax: Double, _ aMin: Double, _ jMax: Double, _ returnAfterFound: Bool) {
        let profile = validProfiles[profileCount]

        var h1 = (a0 * a0 + af * af) / 2 + jMax * vd
        if h1 < 0.0 {
            return
        }
        h1 = sqrt(h1)

        // Solution 1
        profile.t[0] = -(a0 + h1) / jMax
        profile.t[1] = 0
        profile.t[2] = -(af + h1) / jMax
        profile.t[3] = 0
        profile.t[4] = 0
        profile.t[5] = 0
        profile.t[6] = 0

        if profile.checkForVelocity(jMax, aMax, aMin, .UDDU, .NONE) {
            validProfiles[profileCount] = profile
            addProfile()
            if returnAfterFound { return }
        }

        // Solution 2
        profile.t[0] = (-a0 + h1) / jMax
        profile.t[1] = 0
        profile.t[2] = (-af + h1) / jMax
        profile.t[3] = 0
        profile.t[4] = 0
        profile.t[5] = 0
        profile.t[6] = 0

        if profile.checkForVelocity(jMax, aMax, aMin, .UDDU, .NONE) {
            validProfiles[profileCount] = profile
            addProfile()
        }
    }

    /// Only for zero-limits case
    private func timeAllSingleStep(_ profile: inout Profile, _ aMax: Double, _ aMin: Double, _ jMax: Double) -> Bool
    {
        if abs(af - a0) > Double.leastNormalMagnitude {
            return false
        }

        profile.t[0] = 0
        profile.t[1] = 0
        profile.t[2] = 0
        profile.t[3] = 0
        profile.t[4] = 0
        profile.t[5] = 0
        profile.t[6] = 0

        if abs(a0) > Double.leastNormalMagnitude {
            profile.t[3] = vd / a0
            if profile.checkForVelocity(0.0, aMax, aMin, .UDDU, .NONE) {
                return true
            }
        } else if abs(vd) < Double.leastNormalMagnitude {
            if profile.checkForVelocity(0.0, aMax, aMin, .UDDU, .NONE) {
                return true
            }
        }

        return false
    }

    func getProfile(_ input: inout Profile, _ block: inout Block) -> Bool {
        // Zero-limits special case
        if _jMax == 0.0 {
            var p = block.pMin
            p.setBoundary(&input)

            if timeAllSingleStep(&p, _aMax, _aMin, _jMax) {
                block.tMin = p.tSum.last! + p.brake.duration + p.accel.duration
                if abs(a0) > Double.leastNormalMagnitude {
                    block.a = Interval(block.tMin, .infinity)
                }
                return true
            }
            return false
        }

        resetProfiles()
        validProfiles[profileCount].setBoundary(&input)

        if abs(af) < Double.leastNormalMagnitude {
            // There is no blocked interval when af==0, so return after first found profile
            let aMax = (vd >= 0) ? _aMax : _aMin
            let aMin = (vd >= 0) ? _aMin : _aMax
            let jMax = (vd >= 0) ? _jMax : -_jMax

            timeNone(aMax, aMin, jMax, true)
            if hasProfiles() {
                var count = profileCount
                return Block.calculateBlock(&block, &validProfiles, &count)
            }

            timeAcc0(aMax, aMin, jMax, true)
            if hasProfiles() {
                var count = profileCount
                return Block.calculateBlock(&block, &validProfiles, &count)
            }

            timeNone(aMin, aMax, -jMax, true)
            if hasProfiles() {
                var count = profileCount
                return Block.calculateBlock(&block, &validProfiles, &count)
            }

            timeAcc0(aMin, aMax, -jMax, true)

        } else {
            timeNone(_aMax, _aMin, _jMax, false)
            timeNone(_aMin, _aMax, -_jMax, false)
            timeAcc0(_aMax, _aMin, _jMax, false)
            timeAcc0(_aMin, _aMax, -_jMax, false)
        }

        var count = profileCount
        return Block.calculateBlock(&block, &validProfiles, &count)
    }
}
