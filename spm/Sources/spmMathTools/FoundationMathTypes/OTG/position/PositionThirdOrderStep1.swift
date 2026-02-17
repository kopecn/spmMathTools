//
//  PositionThirdOrderStep1.swift
//
//  Created by Nicholas Bergantz on 4/25/24.
//

import Foundation

// Assuming Profile and Block are defined somewhere else in Swift.
typealias ProfileIter = Int  // Swift doesn't use iterators the same way

class PositionThirdOrderStep1 {

    /// initial velocity
    let v0: Double
    /// initial acceleration
    let a0: Double
    /// final velocity
    let vf: Double
    /// final acceleration
    let af: Double
    /// max velocity limit
    let _vMax: Double
    /// min velocity limit
    let _vMin: Double
    /// max acceleration limit
    let _aMax: Double
    /// min acceleration limit
    let _aMin: Double
    /// max Jerk limit
    let _jMax: Double

    // Pre-calculated expressions
    /// position delta
    let pd: Double
    /// initial velocity raised to the power of 2
    let v0P2: Double
    /// final velocity raised to the power of 2
    let vfP2: Double
    /// initial acceleration raised to the power of 2
    let a0P2: Double
    /// initial acceleration raised to the power of 3
    let a0P3: Double
    /// initial acceleration raised to the power of 4
    let a0P4: Double
    /// final acceleration raised to the power of 2
    let afP2: Double
    /// final acceleration raised to the power of 3
    let afP3: Double
    /// final acceleration raised to the power of 4
    let afP4: Double
    /// max jerck raised to the power of 2
    let jMaxP2: Double

    // Max 5 valid profiles + 1 spare for numerical issues
    // CRITICAL: Create separate Profile instances, not shared references
    // Using `repeating:` with a class creates multiple references to the SAME object!
    var validProfiles: [Profile] = (0..<6).map { _ in Profile() }
    private var profileCount: Int = 0

    init(
        p0: Double,
        v0: Double,
        a0: Double,
        pf: Double,
        vf: Double,
        af: Double,
        vMax: Double,
        vMin: Double,
        aMax: Double,
        aMin: Double,
        jMax: Double
    ) {

        self.v0 = v0
        self.a0 = a0
        self.vf = vf
        self.af = af
        self._vMax = vMax
        self._vMin = vMin
        self._aMax = aMax
        self._aMin = aMin
        self._jMax = jMax

        self.pd = pf - p0
        self.v0P2 = v0 * v0
        self.vfP2 = vf * vf
        self.a0P2 = a0 * a0
        self.afP2 = af * af
        self.a0P3 = a0 * a0P2
        self.a0P4 = a0P2 * a0P2
        self.afP3 = af * afP2
        self.afP4 = afP2 * afP2
        self.jMaxP2 = jMax * jMax
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

    private func timeAllVel(
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double,
        _ returnAfterFound: Bool
    ) {
        var profile = validProfiles[profileCount]
        // ACC0_ACC1_VEL
        // NOTE: This profile type explicitly uses aMax/aMin in formulas (lines 110-125),
        // so it inherently respects acceleration limits. No additional checking needed.

        profile.t[0] = (-a0 + aMax) / jMax
        profile.t[1] = (a0P2 / 2 - aMax * aMax - jMax * (v0 - vMax)) / (aMax * jMax)
        profile.t[2] = aMax / jMax
        profile.t[3] =
            (3 * (a0P4 * aMin - afP4 * aMax)
                + 8 * aMax * aMin * (afP3 - a0P3 + 3 * jMax * (a0 * v0 - af * vf))
                + 6 * a0P2 * aMin * (aMax * aMax - 2 * jMax * v0)
                - 6 * afP2 * aMax * (aMin * aMin - 2 * jMax * vf)
                - 12 * jMax
                * (aMax * aMin * (aMax * (v0 + vMax) - aMin * (vf + vMax) - 2 * jMax * pd)
                    + (aMin - aMax) * jMax * vMax * vMax
                    + jMax * (aMax * vfP2 - aMin * v0P2)))
            / (24 * aMax * aMin * jMaxP2 * vMax)
        profile.t[4] = -aMin / jMax
        profile.t[5] = -(afP2 / 2 - aMin * aMin - jMax * (vf - vMax)) / (aMin * jMax)
        profile.t[6] = profile.t[4] + af / jMax

        if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC0_ACC1_VEL) {
            validProfiles[profileCount] = profile
            addProfile()
            if returnAfterFound { return }
        }

        // ACC1_VEL
        // BUG FIX: This profile type calculates times based only on vMax and jMax for the
        // acceleration phase, without explicitly constraining to aMax/aMin. We need to check
        // if the resulting peak acceleration in the acceleration phase would violate limits.
        // The deceleration phase explicitly uses aMin, so it's automatically constrained.

        let tAcc0 = sqrt(a0P2 / (2 * jMaxP2) + (vMax - v0) / jMax)
        let aPeakAcc0 = jMax * tAcc0  // Peak acceleration in acceleration phase

        // Check if peak is within [aMin, aMax] range (handles both UP and DOWN directions)
        let aLimitMin = min(aMin, aMax)
        let aLimitMax = max(aMin, aMax)

        if aPeakAcc0 >= aLimitMin && aPeakAcc0 <= aLimitMax {
            // Acceleration phase peak is within limits. Deceleration phase uses aMin directly,
            // so it's automatically constrained. Generate the profile.
            profile.t[0] = tAcc0 - a0 / jMax
            profile.t[1] = 0
            profile.t[2] = tAcc0
            profile.t[3] =
                (-3 * afP4
                    + 8 * aMin * (afP3 - a0P3)
                    + 24 * aMin * jMax * (a0 * v0 - af * vf)
                    - 6 * afP2 * (aMin * aMin - 2 * jMax * vf)
                    + 12 * jMax
                    * (2 * aMin * jMax * pd
                        + aMin * aMin * (vf + vMax)
                        + jMax * (vMax * vMax - vfP2)
                        + aMin * tAcc0 * (a0P2 - 2 * jMax * (v0 + vMax))))
                / (24 * aMin * jMaxP2 * vMax)
            profile.t[4] = -aMin / jMax
            profile.t[5] = -(afP2 / 2 - aMin * aMin - jMax * (vf - vMax)) / (aMin * jMax)
            profile.t[6] = profile.t[4] + af / jMax

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC1_VEL) {
                validProfiles[profileCount] = profile
                addProfile()
                if returnAfterFound { return }
            }
        }

        // ACC0_VEL
        let tAcc1 = sqrt(afP2 / (2 * jMaxP2) + (vMax - vf) / jMax)

        // BUG FIX: Check if peak acceleration would exceed limits (preserve sign!)
        let aPeakAcc1Acc0vel = jMax * tAcc1

        // Check if peak is within [aMin, aMax] range
        if aPeakAcc1Acc0vel >= aLimitMin && aPeakAcc1Acc0vel <= aLimitMax {
            profile.t[0] = (-a0 + aMax) / jMax
            profile.t[1] = (a0P2 / 2 - aMax * aMax - jMax * (v0 - vMax)) / (aMax * jMax)
            profile.t[2] = aMax / jMax
            profile.t[3] =
                (3 * a0P4
                    + 8 * aMax * (afP3 - a0P3)
                    + 24 * aMax * jMax * (a0 * v0 - af * vf)
                    + 6 * a0P2 * (aMax * aMax - 2 * jMax * v0)
                    - 12 * jMax
                    * (-2 * aMax * jMax * pd
                        + aMax * aMax * (v0 + vMax)
                        + jMax * (vMax * vMax - v0P2)
                        + aMax * tAcc1 * (-afP2 + 2 * (vf + vMax) * jMax)))
                / (24 * aMax * jMaxP2 * vMax)
            profile.t[4] = tAcc1
            profile.t[5] = 0
            profile.t[6] = tAcc1 + af / jMax

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC0_VEL) {
                validProfiles[profileCount] = profile
                addProfile()
                if returnAfterFound { return }
            }
        }

        // VEL
        // BUG FIX: Re-check tAcc0 and tAcc1 for VEL profile (they were calculated earlier)
        // Peak accelerations: aPeakAcc0 (already calculated), aPeakAcc1 (need tAcc1)
        let tAcc1Vel = sqrt(afP2 / (2 * jMaxP2) + (vMax - vf) / jMax)
        let aPeakAcc1Vel = jMax * tAcc1Vel

        // Check if both phases respect limits (directional check, not abs!)
        // NOTE: aLimitMin/Max already calculated above for ACC1_VEL
        if aPeakAcc0 >= aLimitMin && aPeakAcc0 <= aLimitMax
            && aPeakAcc1Vel >= aLimitMin && aPeakAcc1Vel <= aLimitMax
        {
            profile.t[0] = tAcc0 - a0 / jMax
            profile.t[1] = 0
            profile.t[2] = tAcc0
            profile.t[3] =
                (afP3 - a0P3) / (3 * jMaxP2 * vMax)
                + (a0 * v0 - af * vf + (afP2 * tAcc1Vel + a0P2 * tAcc0) / 2) / (jMax * vMax)
                - (v0 / vMax + 1.0) * tAcc0
                - (vf / vMax + 1.0) * tAcc1Vel
                + pd / vMax

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .VEL) {
                validProfiles[profileCount] = profile
                addProfile()
                if returnAfterFound { return }
            }
        }
    }

    private func timeAcc0Acc1(
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double,
        _ returnAfterFound: Bool
    ) {
        var profile = validProfiles[profileCount]

        var h1 =
            (3 * (afP4 * aMax - a0P4 * aMin)
                + aMax * aMin
                * (8 * (a0P3 - afP3)
                    + 3 * aMax * aMin * (aMax - aMin)
                    + 6 * aMin * afP2
                    - 6 * aMax * a0P2)
                + 12 * jMax
                * (aMax * aMin * ((aMax - 2 * a0) * v0 - (aMin - 2 * af) * vf)
                    + aMin * a0P2 * v0
                    - aMax * afP2 * vf))
            / (3 * (aMax - aMin) * jMaxP2)

        h1 += 4 * (aMax * vfP2 - aMin * v0P2 - 2 * aMin * aMax * pd) / (aMax - aMin)

        if h1 >= 0 {
            h1 = sqrt(h1) / 2

            let h2 = a0P2 / (2 * aMax * jMax) + (aMin - 2 * aMax) / (2 * jMax) - v0 / aMax
            let h3 = -afP2 / (2 * aMin * jMax) - (aMax - 2 * aMin) / (2 * jMax) + vf / aMin

            // UDDU: Solution 2
            if h2 > h1 / aMax && h3 > -h1 / aMin {
                profile.t[0] = (-a0 + aMax) / jMax
                profile.t[1] = h2 - h1 / aMax
                profile.t[2] = aMax / jMax
                profile.t[3] = 0
                profile.t[4] = -aMin / jMax
                profile.t[5] = h3 + h1 / aMin
                profile.t[6] = profile.t[4] + af / jMax

                if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC0_ACC1, true) {
                    validProfiles[profileCount] = profile
                    addProfile()
                    if returnAfterFound {
                        return
                    }
                }
            }

            // UDDU: Solution 1
            if h2 > -h1 / aMax && h3 > h1 / aMin {
                profile.t[0] = (-a0 + aMax) / jMax
                profile.t[1] = h2 + h1 / aMax
                profile.t[2] = aMax / jMax
                profile.t[3] = 0
                profile.t[4] = -aMin / jMax
                profile.t[5] = h3 - h1 / aMin
                profile.t[6] = profile.t[4] + af / jMax

                if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC0_ACC1, true) {
                    validProfiles[profileCount] = profile
                    addProfile()
                }
            }
        }
    }

    private func timeAllNoneAcc0Acc1(
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double,
        _ returnAfterFound: Bool
    ) {

        // BUG FIX: Use 'var' instead of 'let' to allow getting fresh references after addProfile()
        var profile: Profile = validProfiles[profileCount]

        // NONE UDDU / UDUD Strategy
        let h2None = (a0P2 - afP2) / (2 * jMax) + (vf - v0)
        let h2P2 = h2None * h2None
        let tMinNone = (a0 - af) / jMax
        let tMaxNone = (aMax - aMin) / jMax

        // Construct polynomials
        var polynomNone: [Double] = [
            0,
            -2 * (a0P2 + afP2 - 2 * jMax * (v0 + vf)) / jMaxP2,
            4 * (a0P3 - afP3 + 3 * jMax * (af * vf - a0 * v0)) / (3 * jMax * jMaxP2) - 4 * pd / jMax,
            -h2P2 / jMaxP2,
        ]

        let h3Acc0 = (a0P2 - afP2) / (2 * aMax * jMax) + (vf - v0) / aMax
        let tMinAcc0 = (aMax - af) / jMax
        let tMaxAcc0 = (aMax - aMin) / jMax

        let h0Acc0 =
            3 * (afP4 - a0P4) + 8 * (a0P3 - afP3) * aMax + 24 * aMax * jMax * (af * vf - a0 * v0) - 6 * a0P2
            * (aMax * aMax - 2 * jMax * v0) + 6 * afP2 * (aMax * aMax - 2 * jMax * vf) + 12 * jMax
            * (jMax * (vfP2 - v0P2 - 2 * aMax * pd) - aMax * aMax * (vf - v0))
        let h2Acc0 = -afP2 + aMax * aMax + 2 * jMax * vf

        var polynomAcc0: [Double] = [
            -2 * aMax / jMax,
            h2Acc0 / jMaxP2,
            0,
            h0Acc0 / (12 * jMaxP2 * jMaxP2),
        ]

        let h3Acc1 = -(a0P2 + afP2) / (2 * jMax * aMin) + aMin / jMax + (vf - v0) / aMin
        let tMinAcc1 = (aMin - a0) / jMax
        let tMaxAcc1 = (aMax - a0) / jMax

        let h0Acc1 =
            (a0P4 - afP4) / 4 + 2 * (afP3 - a0P3) * aMin / 3 + (a0P2 - afP2) * aMin * aMin / 2 + jMax
            * (afP2 * vf + a0P2 * v0 + 2 * aMin * (jMax * pd - a0 * v0 - af * vf) + aMin * aMin * (v0 + vf) + jMax
                * (v0P2 - vfP2))
        let h2Acc1 = a0P2 - a0 * aMin + 2 * jMax * v0

        var polynomAcc1: [Double] = [
            2 * (2 * a0 - aMin) / jMax,
            (5 * a0P2 + aMin * (aMin - 6 * a0) + 2 * jMax * v0) / jMaxP2,
            2 * (a0 - aMin) * h2Acc1 / (jMaxP2 * jMax),
            h0Acc1 / (jMaxP2 * jMaxP2),
        ]

        let polynomAcc0Min = [
            polynomAcc0[0] + 4 * tMinAcc0,
            polynomAcc0[1] + (3 * polynomAcc0[0] + 6 * tMinAcc0) * tMinAcc0,
            polynomAcc0[2] + (2 * polynomAcc0[1] + (3 * polynomAcc0[0] + 4 * tMinAcc0) * tMinAcc0) * tMinAcc0,
            polynomAcc0[3]
                + (polynomAcc0[2] + (polynomAcc0[1] + (polynomAcc0[0] + tMinAcc0) * tMinAcc0) * tMinAcc0)
                    * tMinAcc0,
        ]

        let polynomAcc0HasSolution =
            polynomAcc0Min[0] < 0 || polynomAcc0Min[1] < 0 || polynomAcc0Min[2] < 0 || polynomAcc0Min[3] <= 0

        let polynomAcc1_has_solution = polynomAcc1.contains(where: { $0 < 0 }) || polynomAcc1[3] <= 0

        // Solve quartic polynomial roots
        let rootsNone = solveQuarticMonic(&polynomNone).sorted()
        let rootsAcc0 = polynomAcc0HasSolution ? solveQuarticMonic(&polynomAcc0).sorted() : []
        let rootsAcc1 = polynomAcc1_has_solution ? solveQuarticMonic(&polynomAcc1).sorted() : []

        // Process roots for NONE
        for var t in rootsNone {
            guard t >= tMinNone && t <= tMaxNone else { continue }
            if t > Double.leastNormalMagnitude {
                let h1 = jMax * t * t
                let orig =
                    -h2P2 / (4 * jMax * t) + h2None * (af / jMax + t)
                    + (4 * a0P3 + 2 * afP3 - 6 * a0P2 * (af + 2 * jMax * t) + 12 * (af - a0) * jMax * v0 + 3
                        * jMaxP2 * (-4 * pd + (h1 + 8 * v0) * t)) / (12 * jMaxP2)
                let deriv = h2None + 2 * v0 - a0P2 / jMax + h2P2 / (4 * h1) + (3 * h1) / 4
                t -= orig / deriv
            }

            let h0 = h2None / (2 * jMax * t)
            profile.t[0] = h0 + t / 2 - a0 / jMax
            profile.t[1] = 0
            profile.t[2] = t
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = -h0 + t / 2 + af / jMax

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                validProfiles[profileCount] = profile
                addProfile()
                if returnAfterFound { return }
                // BUG FIX: Get fresh profile reference to prevent corrupting previously added profile
                profile = validProfiles[profileCount]
            }
        }

        // Process ACC0 roots
        for var t in rootsAcc0 {
            guard t >= tMinAcc0 && t <= tMaxAcc0 else { continue }
            if t > Double.leastNormalMagnitude {
                let h1 = jMax * t
                let orig = h0Acc0 / (12 * jMaxP2 * t) + t * (h2Acc0 + h1 * (h1 - 2 * aMax))
                let deriv = 2 * (h2Acc0 + h1 * (2 * h1 - 3 * aMax))
                t -= orig / deriv
            }

            profile.t[0] = (-a0 + aMax) / jMax
            profile.t[1] = h3Acc0 - 2 * t + jMax / aMax * t * t
            profile.t[2] = t
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = (af - aMax) / jMax + t

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC0) {
                validProfiles[profileCount] = profile
                addProfile()
                if returnAfterFound { return }
                // BUG FIX: Get fresh profile reference to prevent corrupting previously added profile
                profile = validProfiles[profileCount]
            }
        }

        // Process ACC1 roots
        for var t in rootsAcc1 {
            guard t >= tMinAcc1 && t <= tMaxAcc1 else { continue }
            if t > Double.leastNormalMagnitude {
                let h5 = a0P3 + 2 * jMax * a0 * v0
                var h1 = jMax * t
                var orig =
                    -(h0Acc1 / 2 + h1
                    * (h5 + a0 * (aMin - 2 * h1) * (aMin - h1) + a0P2 * (5 * h1 / 2 - 2 * aMin) + aMin * aMin * h1 / 2
                        + jMax * (h1 / 2 - aMin) * (h1 * t + 2 * v0)))
                    / jMax
                var deriv = (aMin - a0 - h1) * (h2Acc1 + h1 * (4 * a0 - aMin + 2 * h1))
                t -= min(orig / deriv, t)

                h1 = jMax * t
                orig =
                    -(h0Acc1 / 2 + h1
                    * (h5 + a0 * (aMin - 2 * h1) * (aMin - h1) + a0P2 * (5 * h1 / 2 - 2 * aMin) + aMin * aMin * h1 / 2
                        + jMax * (h1 / 2 - aMin) * (h1 * t + 2 * v0)))
                    / jMax
                if abs(orig) > 1e-9 {
                    deriv = (aMin - a0 - h1) * (h2Acc1 + h1 * (4 * a0 - aMin + 2 * h1))
                    t -= orig / deriv

                    h1 = jMax * t
                    orig =
                        -(h0Acc1 / 2 + h1
                        * (h5 + a0 * (aMin - 2 * h1) * (aMin - h1) + a0P2 * (5 * h1 / 2 - 2 * aMin) + aMin * aMin * h1
                            / 2 + jMax * (h1 / 2 - aMin) * (h1 * t + 2 * v0)))
                        / jMax
                    if abs(orig) > 1e-9 {
                        deriv = (aMin - a0 - h1) * (h2Acc1 + h1 * (4 * a0 - aMin + 2 * h1))
                        t -= orig / deriv
                    }
                }
            }

            profile.t[0] = t
            profile.t[1] = 0
            profile.t[2] = (a0 - aMin) / jMax + t
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = h3Acc1 - (2 * a0 + jMax * t) * t / aMin
            profile.t[6] = (af - aMin) / jMax

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC1, true) {
                validProfiles[profileCount] = profile
                addProfile()
                if returnAfterFound { return }
                // BUG FIX: Get fresh profile reference to prevent corrupting previously added profile
                profile = validProfiles[profileCount]
            }
        }
    }

    // Only for numerical issues
    private func timeAcc1VelTwoStep(_ vMax: Double, _ vMin: Double, _ aMax: Double, _ aMin: Double, _ jMax: Double)
    {
        var profile: Profile = validProfiles[profileCount]

        profile.t[0] = 0
        profile.t[1] = 0
        profile.t[2] = a0 / jMax
        profile.t[3] =
            -(3 * afP4 - 8 * aMin * (afP3 - a0P3) - 24 * aMin * jMax * (a0 * v0 - af * vf) + 6 * afP2
            * (aMin * aMin - 2 * jMax * vf) - 12 * jMax
            * (2 * aMin * jMax * pd + aMin * aMin * (vf + vMax) + jMax * (vMax * vMax - vfP2) + aMin * a0
                * (a0P2 - 2 * jMax * (v0 + vMax)) / jMax))
            / (24 * aMin * jMaxP2 * vMax)
        profile.t[4] = -aMin / jMax
        profile.t[5] = -(afP2 / 2 - aMin * aMin + jMax * (vMax - vf)) / (aMin * jMax)
        profile.t[6] = profile.t[4] + af / jMax

        if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC1_VEL) {
            validProfiles[profileCount] = profile
            addProfile()
        }
    }

    private func timeAcc0TwoStep(_ vMax: Double, _ vMin: Double, _ aMax: Double, _ aMin: Double, _ jMax: Double) {
        var profile: Profile = validProfiles[profileCount]

        // Two-step profile
        do {
            profile.t[0] = 0
            profile.t[1] = (afP2 - a0P2 + 2 * jMax * (vf - v0)) / (2 * a0 * jMax)
            profile.t[2] = (a0 - af) / jMax
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = 0

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC0) {
                validProfiles[profileCount] = profile
                addProfile()
                return
            }
        }

        // Three-step profile - removed pf
        do {
            profile.t[0] = (-a0 + aMax) / jMax
            profile.t[1] = (a0P2 + afP2 - 2 * aMax * aMax + 2 * jMax * (vf - v0)) / (2 * aMax * jMax)
            profile.t[2] = (-af + aMax) / jMax
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = 0

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC0) {
                validProfiles[profileCount] = profile
                addProfile()
                return
            }
        }

        // Three-step profile - removed aMax
        do {
            let h0 = 3 * (afP2 - a0P2 + 2 * jMax * (v0 + vf))
            let h2 = a0P3 + 2 * afP3 + 6 * jMaxP2 * pd + 6 * (af - a0) * jMax * vf - 3 * a0 * afP2

            let discriminantNumerator =
                2
                * (2 * h2 * h2 + h0
                    * (a0P4 - 6 * a0P2 * (afP2 + 2 * jMax * vf) + 8 * a0
                        * (afP3 + 3 * jMaxP2 * pd + 3 * af * jMax * vf) - 3
                        * (afP4 + 4 * afP2 * jMax * vf + 4 * jMaxP2 * (vfP2 - v0P2))))

            let h1 = sqrt(discriminantNumerator) * abs(jMax) / jMax

            profile.t[0] =
                (4 * afP3 + 2 * a0P3 - 6 * a0 * afP2 + 12 * jMaxP2 * pd + 12 * (af - a0) * jMax * vf + h1)
                / (2 * jMax * h0)

            profile.t[1] = -h1 / (jMax * h0)

            profile.t[2] =
                (-4 * a0P3 - 2 * afP3 + 6 * a0P2 * af + 12 * jMaxP2 * pd - 12 * (af - a0) * jMax * v0 + h1)
                / (2 * jMax * h0)

            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = 0

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC0) {
                validProfiles[profileCount] = profile
                addProfile()
                return
            }
        }

        // Three-step profile - t = (aMax - aMin)/jMax
        do {
            let t = (aMax - aMin) / jMax

            profile.t[0] = (-a0 + aMax) / jMax
            profile.t[1] = (a0P2 - afP2) / (2 * aMax * jMax) + (vf - v0 + jMax * t * t) / aMax - 2 * t
            profile.t[2] = t
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = (af - aMin) / jMax

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC0) {
                validProfiles[profileCount] = profile
                addProfile()
                return
            }
        }
    }

    private func timeVelTwoStep(_ vMax: Double, _ vMin: Double, _ aMax: Double, _ aMin: Double, _ jMax: Double) {
        var profile = validProfiles[profileCount]

        let h1 = sqrt(afP2 / (2 * jMaxP2) + (vMax - vf) / jMax)

        // Four-step profile: Solution 3/4
        do {
            profile.t[0] = -a0 / jMax
            profile.t[1] = 0
            profile.t[2] = 0
            profile.t[3] =
                (afP3 - a0P3) / (3 * jMaxP2 * vMax) + (a0 * v0 - af * vf + (afP2 * h1) / 2) / (jMax * vMax)
                - (vf / vMax + 1.0) * h1 + pd / vMax
            profile.t[4] = h1
            profile.t[5] = 0
            profile.t[6] = h1 + af / jMax

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .VEL) {
                validProfiles[profileCount] = profile
                addProfile()
                return
            }
        }

        // Four-step profile: Alternative
        do {
            profile.t[0] = 0
            profile.t[1] = 0
            profile.t[2] = a0 / jMax
            profile.t[3] =
                (afP3 - a0P3) / (3 * jMaxP2 * vMax) + (a0 * v0 - af * vf + (afP2 * h1 + a0P3 / jMax) / 2)
                / (jMax * vMax) - (v0 / vMax + 1.0) * a0 / jMax - (vf / vMax + 1.0) * h1 + pd / vMax
            profile.t[4] = h1
            profile.t[5] = 0
            profile.t[6] = h1 + af / jMax

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .VEL) {
                validProfiles[profileCount] = profile
                addProfile()
                return
            }
        }
    }

    private func timeNoneTwoStep(_ vMax: Double, _ vMin: Double, _ aMax: Double, _ aMin: Double, _ jMax: Double) {
        var profile = validProfiles[profileCount]

        // Two step
        do {
            let h0 = sqrt((a0P2 + afP2) / 2 + jMax * (vf - v0)) * abs(jMax) / jMax
            profile.t[0] = (h0 - a0) / jMax
            profile.t[1] = 0
            profile.t[2] = (h0 - af) / jMax
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = 0

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                validProfiles[profileCount] = profile
                addProfile()
                return
            }
        }

        // Single step
        do {
            profile.t[0] = (af - a0) / jMax
            profile.t[1] = 0
            profile.t[2] = 0
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = 0

            if profile.check(jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                validProfiles[profileCount] = profile
                addProfile()
                return
            }
        }

    }

    // Only for zero-limits case
    private func timeAllSingleStep(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) -> Bool {

        if abs(af - a0) > Double.leastNormalMagnitude {
            return false
        }

        profile.t = [Double](repeating: 0.0, count: 7)

        if abs(a0) > Double.leastNormalMagnitude {
            let q = sqrt(2 * a0 * pd + v0P2)

            // Solution 1
            profile.t[3] = (-v0 + q) / a0
            if profile.t[3] >= 0.0 && profile.check(0.0, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                return true
            }

            // Solution 2
            profile.t[3] = -(v0 + q) / a0
            if profile.t[3] >= 0.0 && profile.check(0.0, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                return true
            }

        } else if abs(v0) > Double.ulpOfOne {
            profile.t[3] = pd / v0
            if profile.check(0.0, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                return true
            }

        } else if abs(pd) < Double.ulpOfOne {
            if profile.check(0.0, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                return true
            }
        }

        return false
    }

    func getProfile(_ input: inout Profile, _ block: inout Block) -> Bool {
        // Zero-limits special case
        if _jMax == 0.0 || _aMax == 0.0 || _aMin == 0.0 {
            var p = block.pMin
            p.setBoundary(input)

            if timeAllSingleStep(&p, _vMax, _vMin, _aMax, _aMin, _jMax) {
                block.pMin = p
                block.tMin = p.tSum.last! + p.brake.duration + p.accel.duration
                if abs(v0) > Double.leastNormalMagnitude || abs(a0) > Double.leastNormalMagnitude {
                    block.a = Interval(block.tMin, .infinity)
                }
                return true
            }
            return false
        }

        resetProfiles()
        validProfiles[profileCount].setBoundary(input)

        if abs(vf) < Double.ulpOfOne && abs(af) < Double.ulpOfOne {
            let vMax = (pd >= 0) ? _vMax : _vMin
            let vMin = (pd >= 0) ? _vMin : _vMax
            let aMax = (pd >= 0) ? _aMax : _aMin
            let aMin = (pd >= 0) ? _aMin : _aMax
            let jMax = (pd >= 0) ? _jMax : -_jMax

            // Special case: when v0, a0, and pd are all essentially zero
            if abs(v0) < Double.ulpOfOne && abs(a0) < Double.ulpOfOne && abs(pd) < Double.ulpOfOne {
                // For the trivial case where we're already at the target, try timeAllSingleStep first
                var p = validProfiles[profileCount]
                if timeAllSingleStep(&p, vMax, vMin, aMax, aMin, jMax) {
                    profileCount += 1
                } else {
                    timeAllNoneAcc0Acc1(vMax, vMin, aMax, aMin, jMax, true)
                }
            } else {
                // There is no blocked interval when vf==0 && af==0, so return after first found profile
                timeAllVel(vMax, vMin, aMax, aMin, jMax, true)
                if hasProfiles() {
                    var count = profileCount
                    return Block.calculateBlock(&block, &validProfiles, &count)
                }

                timeAllNoneAcc0Acc1(vMax, vMin, aMax, aMin, jMax, true)
            }
            if hasProfiles() {
                var count = profileCount
                return Block.calculateBlock(&block, &validProfiles, &count)
            }

            // Try profiles one by one and exit early if found
            let directions: [(Double, Double, Double, Double, Double)] = [
                (vMax, vMin, aMax, aMin, jMax),
                (vMax, vMin, aMax, aMin, jMax),
                (vMax, vMin, aMax, aMin, jMax),
                (vMin, vMax, aMin, aMax, -jMax),
                (vMin, vMax, aMin, aMax, -jMax),
                (vMin, vMax, aMin, aMax, -jMax),
            ]

            for (vMax, vMin, aMax, aMin, jMax) in directions {
                timeAllVel(vMax, vMin, aMax, aMin, jMax, true)
                if hasProfiles() {
                    var count = profileCount
                    return Block.calculateBlock(&block, &validProfiles, &count)
                }

                timeAllNoneAcc0Acc1(vMax, vMin, aMax, aMin, jMax, true)
                if hasProfiles() {
                    var count = profileCount
                    return Block.calculateBlock(&block, &validProfiles, &count)
                }

                timeAcc0Acc1(vMax, vMin, aMax, aMin, jMax, true)
                if hasProfiles() {
                    var count = profileCount
                    return Block.calculateBlock(&block, &validProfiles, &count)
                }
            }

        } else {
            // General search for profiles
            // Prioritize profiles with proper limit handling to avoid acceleration violations
            timeAllVel(_vMax, _vMin, _aMax, _aMin, _jMax, false)
            timeAllVel(_vMin, _vMax, _aMin, _aMax, -_jMax, false)
            timeAcc0Acc1(_vMax, _vMin, _aMax, _aMin, _jMax, false)
            timeAcc0Acc1(_vMin, _vMax, _aMin, _aMax, -_jMax, false)
            timeAllNoneAcc0Acc1(_vMax, _vMin, _aMax, _aMin, _jMax, false)
            timeAllNoneAcc0Acc1(_vMin, _vMax, _aMin, _aMax, -_jMax, false)
        }

        var count = profileCount
        return Block.calculateBlock(&block, &validProfiles, &count)
    }
}
