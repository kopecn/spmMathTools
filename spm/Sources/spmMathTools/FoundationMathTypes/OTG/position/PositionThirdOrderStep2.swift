//
//  PositionThirdOrderStep2.swift
//
//
//  Created by Nicholas Bergantz on 4/25/24.
//

import Foundation

private let tolerance: Double = 1e-14

// Helper function for squaring
private func pow2(_ x: Double) -> Double {
    x * x
}

/// Mathematical equations for Step 2 in third-order position interface: Time synchronization
class PositionThirdOrderStep2 {

    let p0: Double  // Store p0 to reconstruct profile.p[0]
    let v0: Double
    let a0: Double
    let tf: Double
    let pf: Double  // Store pf for profile.pf
    let vf: Double
    let af: Double
    let _vMax: Double
    let _vMin: Double
    let _aMax: Double
    let _aMin: Double
    let _jMax: Double

    // Pre-calculated expressions
    let pd: Double
    let tfP2: Double
    let tfP3: Double
    let tfP4: Double
    let vd: Double
    let vdP2: Double
    let ad: Double
    let adP2: Double
    let v0P2: Double
    let vfP2: Double
    let a0P2: Double
    let a0P3: Double
    let a0P4: Double
    let a0P5: Double
    let a0P6: Double
    let afP2: Double
    let afP3: Double
    let afP4: Double
    let afP5: Double
    let afP6: Double
    let jMaxP2: Double
    let g1: Double
    let g2: Double

    // Configuration flag for jerk minimization
    let minimize_jerk: Bool = false

    init(
        _ tf: Double,
        _ p0: Double,
        _ v0: Double,
        _ a0: Double,
        _ pf: Double,
        _ vf: Double,
        _ af: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) {
        self.p0 = p0
        self.v0 = v0
        self.a0 = a0
        self.tf = tf
        self.pf = pf
        self.vf = vf
        self.af = af
        self._vMax = vMax
        self._vMin = vMin
        self._aMax = aMax
        self._aMin = aMin
        self._jMax = jMax
        self.pd = pf - p0
        self.tfP2 = tf * tf
        self.tfP3 = tfP2 * tf
        self.tfP4 = tfP2 * tfP2

        self.vd = vf - v0
        self.vdP2 = vd * vd
        self.v0P2 = v0 * v0
        self.vfP2 = vf * vf

        self.ad = af - a0
        self.adP2 = ad * ad
        self.a0P2 = a0 * a0
        self.afP2 = af * af

        self.a0P3 = a0 * a0P2
        self.a0P4 = a0P2 * a0P2
        self.a0P5 = a0P3 * a0P2
        self.a0P6 = a0P4 * a0P2
        self.afP3 = af * afP2
        self.afP4 = afP2 * afP2
        self.afP5 = afP3 * afP2
        self.afP6 = afP4 * afP2

        self.jMaxP2 = jMax * jMax

        self.g1 = -pd + tf * v0
        self.g2 = -2 * pd + tf * (v0 + vf)
    }

    func timeAcc0Acc1Vel(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) -> Bool {
        /// Profile UDDU, Solution 1

        if (2 * (aMax - aMin) + ad) / jMax < tf {
            let h1 =
                sqrt(
                    (a0P4 + afP4 - 4 * a0P3 * (2 * aMax + aMin) / 3 - 4 * afP3 * (aMax + 2 * aMin) / 3 + 2
                        * (a0P2 - afP2) * aMax * aMax + (4 * a0 * aMax - 2 * a0P2)
                        * (afP2 - 2 * af * aMin + (aMin - aMax) * aMin + 2 * jMax * (aMin * tf - vd)) + 2
                        * afP2 * (aMin * aMin + 2 * jMax * (aMax * tf - vd)) + 4 * jMax
                        * (2 * aMin * (af * vd + jMax * g1) + (aMax * aMax - aMin * aMin) * vd + jMax * vdP2)
                        + 8 * aMax * jMaxP2 * (pd - tf * vf)) / (aMax * aMin) + 4 * afP2 + 2 * a0P2
                        + (4 * af + aMax - aMin) * (aMax - aMin) + 4 * jMax * (aMin - aMax + jMax * tf - 2 * af)
                        * tf
                ) * abs(jMax) / jMax

            profile.t[0] = (-a0 + aMax) / jMax
            profile.t[1] =
                (-(afP2 - a0P2 + 2 * aMax * aMax + aMin * (aMin - 2 * ad - 3 * aMax) + 2 * jMax
                    * (aMin * tf - vd)) + aMin * h1) / (2 * (aMax - aMin) * jMax)
            profile.t[2] = aMax / jMax
            profile.t[3] = (aMin - aMax + h1) / (2 * jMax)
            profile.t[4] = -aMin / jMax
            profile.t[5] =
                tf
                - (profile.t[0] + profile.t[1] + profile.t[2] + profile.t[3] + 2 * profile.t[4] + af / jMax)
            profile.t[6] = profile.t[4] + af / jMax

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC0_ACC1_VEL) {
                return true
            }
        }

        // Profile UDUD
        if (-a0 + 4 * aMax - af) / jMax < tf {
            profile.t[0] = (-a0 + aMax) / jMax
            profile.t[1] =
                (3 * (a0P4 + afP4) - 4 * (a0P3 + afP3) * aMax - 4 * afP3 * aMax + 24 * (a0 + af) * aMax
                    * aMax * aMax - 6 * (afP2 + a0P2) * (aMax * aMax - 2 * jMax * vd) + 6 * a0P2
                    * (afP2 - 2 * af * aMax - 2 * aMax * jMax * tf) - 12 * aMax * aMax
                    * (2 * aMax * aMax - 2 * aMax * jMax * tf + jMax * vd) - 24 * af * aMax * jMax * vd + 12
                    * jMaxP2 * (2 * aMax * g1 + vdP2))
                / (12 * aMax * jMax
                    * (a0P2 + afP2 - 2 * (a0 + af) * aMax + 2 * (aMax * aMax - aMax * jMax * tf + jMax * vd)))
            profile.t[2] = aMax / jMax
            profile.t[3] =
                (-a0P2 - afP2 + 2 * aMax * (a0 + af - 2 * aMax) - 2 * jMax * vd) / (2 * aMax * jMax) + tf
            profile.t[4] = profile.t[2]
            profile.t[5] =
                tf
                - (profile.t[0] + profile.t[1] + profile.t[2] + profile.t[3] + 2 * profile.t[4] - af / jMax)
            profile.t[6] = profile.t[4] - af / jMax

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDUD, .ACC0_ACC1_VEL) {
                return true
            }
        }

        return false
    }

    func timeAcc1Vel(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) -> Bool {

        /// Profile UDDU

        let ph1uddu = a0P2 + afP2 - aMin * (a0 + 2 * af - aMin) - 2 * jMax * (vd - aMin * tf)
        let ph2uddu = 2 * aMin * (jMax * g1 + af * vd) - aMin * aMin * vd + jMax * vdP2
        let ph3uddu = afP2 + aMin * (aMin - 2 * af) - 2 * jMax * (vd - aMin * tf)

        var polynomuddu: [Double] = [0, 0, 0, 0]
        polynomuddu[0] = (2 * (2 * a0 - aMin)) / jMax
        polynomuddu[1] = (4 * a0P2 + ph1uddu - 3 * a0 * aMin) / jMaxP2
        polynomuddu[2] = (2 * a0 * ph1uddu) / (jMaxP2 * jMax)
        polynomuddu[3] =
            (3 * (a0P4 + afP4) - 4 * (a0P3 + 2 * afP3) * aMin + 6 * afP2
                * (aMin * aMin - 2 * jMax * vd) + 12 * jMax * ph2uddu + 6 * a0P2 * ph3uddu)
            / (12 * jMaxP2 * jMaxP2)

        let tMinuddu = -a0 / jMax
        let tMaxuddu = min((tf + 2 * aMin / jMax - (a0 + af) / jMax) / 2, (aMax - a0) / jMax)

        let rootsuddu = solveQuarticMonic(&polynomuddu).sorted()
        for t0 in rootsuddu {
            if t0 < tMinuddu || t0 > tMaxuddu {
                continue
            }
            var t = t0

            /// Single Newton step (regarding pd)
            if abs(a0 + jMax * t) > 16 * .ulpOfOne {
                let h0 = jMax * t * t
                let orig =
                    -pd
                    + (3 * (a0P4 + afP4) - 8 * afP3 * aMin - 4 * a0P3 * aMin + 6 * afP2
                        * (aMin * aMin + 2 * jMax * (h0 - vd)) + 6 * a0P2
                        * (afP2 - 2 * af * aMin + aMin * aMin + 2 * aMin * jMax * (-2 * t + tf) + 2 * jMax
                            * (5 * h0 - vd)) + 24 * a0 * jMax * t
                        * (a0P2 + afP2 - 2 * af * aMin + aMin * aMin + 2 * jMax * (aMin * (-t + tf) + h0 - vd))
                        - 24 * af * aMin * jMax * (h0 - vd) + 12 * jMax
                        * (aMin * aMin * (h0 - vd) + jMax * (h0 - vd) * (h0 - vd))) / (24 * aMin * jMaxP2)
                    + h0 * (tf - t) + tf * v0
                let deriv =
                    (a0 + jMax * t)
                    * ((a0P2 + afP2) / (aMin * jMax) + (aMin - a0 - 2 * af) / jMax
                        + (4 * a0 * t + 2 * h0 - 2 * vd) / aMin + 2 * tf - 3 * t)

                t -= orig / deriv
            }

            let h1 = -((a0P2 + afP2) / 2 + jMax * (-vd + 2 * a0 * t + jMax * t * t)) / aMin

            profile.t[0] = t
            profile.t[1] = 0
            profile.t[2] = a0 / jMax + t
            profile.t[3] = tf - (h1 - aMin + a0 + af) / jMax - 2 * t
            profile.t[4] = -aMin / jMax
            profile.t[5] = (h1 + aMin) / jMax
            profile.t[6] = profile.t[4] + af / jMax

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC1_VEL) {
                return true
            }
        }

        /// Profile UDUD

        let ph1udud = a0P2 - afP2 + (2 * af - a0) * aMax - aMax * aMax - 2 * jMax * (vd - aMax * tf)
        let ph2udud = aMax * aMax + 2 * jMax * vd
        let ph3udud = afP2 + ph2udud - 2 * aMax * (af + jMax * tf)
        let ph4udud = 2 * aMax * jMax * g1 + aMax * aMax * vd + jMax * vdP2

        var polynom: [Double] = [0, 0, 0, 0]
        polynom[0] = (4 * a0 - 2 * aMax) / jMax
        polynom[1] = (4 * a0P2 - 3 * a0 * aMax + ph1udud) / jMaxP2
        polynom[2] = (2 * a0 * ph1udud) / (jMaxP2 * jMax)
        polynom[3] =
            (3 * (a0P4 + afP4) - 4 * (a0P3 + 2 * afP3) * aMax - 24 * af * aMax * jMax * vd + 12 * jMax
                * ph4udud - 6 * a0P2 * ph3udud + 6 * afP2 * ph2udud) / (12 * jMaxP2 * jMaxP2)

        let tMin = -a0 / jMax
        let tMax = min((tf + ad / jMax - 2 * aMax / jMax) / 2, (aMax - a0) / jMax)

        let roots = solveQuarticMonic(&polynom).sorted()
        for t in roots {
            if t > tMax || t < tMin {
                continue
            }

            let h1 = ((a0P2 - afP2) / 2 + jMaxP2 * t * t - jMax * (vd - 2 * a0 * t)) / aMax

            profile.t[0] = t
            profile.t[1] = 0
            profile.t[2] = t + a0 / jMax
            profile.t[3] = tf + (h1 + ad - aMax) / jMax - 2 * t
            profile.t[4] = aMax / jMax
            profile.t[5] = -(h1 + aMax) / jMax
            profile.t[6] = profile.t[4] - af / jMax

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDUD, .ACC1_VEL) {
                return true
            }
        }

        return false
    }

    func timeAcc0Vel(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) -> Bool {

        if tf < max((-a0 + aMax) / jMax, 0.0) + max(aMax / jMax, 0.0) {
            return false
        }

        let ph1 = 12 * jMax * (-aMax * aMax * vd - jMax * vdP2 + 2 * aMax * jMax * (-pd + tf * vf))

        /// Profile UDDU
        var polynomuddu: [Double] = [0, 0, 0, 0]
        polynomuddu[0] = (2 * aMax) / jMax
        polynomuddu[1] =
            (a0P2 - afP2 + 2 * ad * aMax + aMax * aMax + 2 * jMax * (vd - aMax * tf)) / jMaxP2
        polynomuddu[2] = 0
        polynomuddu[3] =
            -(-3 * (a0P4 + afP4) + 4 * (afP3 + 2 * a0P3) * aMax - 12 * a0 * aMax
            * (afP2 - 2 * jMax * vd) + 6 * a0P2 * (afP2 - aMax * aMax - 2 * jMax * vd) + 6 * afP2
            * (aMax * aMax - 2 * aMax * jMax * tf + 2 * jMax * vd) + ph1) / (12 * jMaxP2 * jMaxP2)

        let tMinuddu = -af / jMax
        let tMaxuddu = min(tf - (2 * aMax - a0) / jMax, -aMin / jMax)

        let rootsuddu = solveQuarticMonic(&polynomuddu).sorted()
        for t0 in rootsuddu {
            if t0 < tMinuddu || t0 > tMaxuddu {
                continue
            }
            var t = t0

            // Single Newton step (regarding pd)
            if t > .ulpOfOne {
                let h1 = jMax * t * t + vd
                let orig =
                    (-3 * (a0P4 + afP4) + 4 * (afP3 + 2 * a0P3) * aMax - 24 * af * aMax * jMaxP2 * t
                        * t - 12 * a0 * aMax * (afP2 - 2 * jMax * h1) + 6 * a0P2
                        * (afP2 - aMax * aMax - 2 * jMax * h1) + 6 * afP2
                        * (aMax * aMax - 2 * aMax * jMax * tf + 2 * jMax * h1) - 12 * jMax
                        * (aMax * aMax * h1 + jMax * h1 * h1 + 2 * aMax * jMax
                            * (pd + jMax * t * t * (t - tf) - tf * vf)))
                    / (24 * aMax * jMaxP2)
                let deriv =
                    -t
                    * (a0P2 - afP2 + 2 * aMax * (ad - jMax * tf) + aMax * aMax + 3 * aMax * jMax * t + 2
                        * jMax * h1) / aMax

                t -= orig / deriv
            }

            let h1uddu = ((a0P2 - afP2) / 2 + jMax * (jMax * t * t + vd)) / aMax

            profile.t[0] = (-a0 + aMax) / jMax
            profile.t[1] = (h1uddu - aMax) / jMax
            profile.t[2] = aMax / jMax
            profile.t[3] = tf - (h1uddu + ad + aMax) / jMax - 2 * t
            profile.t[4] = t
            profile.t[5] = 0
            profile.t[6] = af / jMax + t

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC0_VEL) {
                return true
            }
        }

        // Profile UDUD

        var polynom: [Double] = [0, 0, 0, 0]
        polynom[0] = (-2 * aMax) / jMax
        polynom[1] =
            -(a0P2 + afP2 - 2 * (a0 + af) * aMax + aMax * aMax + 2 * jMax * (vd - aMax * tf))
            / jMaxP2
        polynom[2] = 0
        polynom[3] =
            (3 * (a0P4 + afP4) - 4 * (afP3 + 2 * a0P3) * aMax + 6 * a0P2
                * (afP2 + aMax * aMax + 2 * jMax * vd) - 12 * a0 * aMax * (afP2 + 2 * jMax * vd) + 6
                * afP2 * (aMax * aMax - 2 * aMax * jMax * tf + 2 * jMax * vd) - ph1)
            / (12 * jMaxP2 * jMaxP2)

        let tMin = af / jMax
        let tMax = min(tf - aMax / jMax, aMax / jMax)

        let roots = solveQuarticMonic(&polynom).sorted()
        for t0 in roots {
            if t0 < tMin || t0 > tMax {
                continue
            }
            var t = t0

            /// Single Newton step (regarding pd)
            let h1ududpd = jMax * t * t - vd
            let orig =
                -(3 * (a0P4 + afP4) - 4 * (2 * a0P3 + afP3) * aMax + 24 * af * aMax * jMaxP2 * t * t
                - 12 * a0 * aMax * (afP2 - 2 * jMax * h1ududpd) + 6 * a0P2
                * (afP2 + aMax * aMax - 2 * jMax * h1ududpd) + 6 * afP2
                * (aMax * aMax - 2 * jMax * (tf * aMax + h1ududpd)) + 12 * jMax
                * (-aMax * aMax * h1ududpd + jMax * h1ududpd * h1ududpd - 2 * aMax * jMax
                    * (-pd + jMax * t * t * (t - tf) + tf * vf)))
                / (24 * aMax * jMaxP2)
            let deriv =
                t
                * (a0P2 + afP2 - 2 * jMax * h1ududpd - 2 * (a0 + af + jMax * tf) * aMax + aMax * aMax + 3
                    * aMax * jMax * t) / aMax

            t -= orig / deriv

            let h1udud = ((a0P2 + afP2) / 2 + jMax * (vd - jMax * t * t)) / aMax

            profile.t[0] = (-a0 + aMax) / jMax
            profile.t[1] = (h1udud - aMax) / jMax
            profile.t[2] = aMax / jMax
            profile.t[3] = tf - (h1udud - a0 - af + aMax) / jMax - 2 * t
            profile.t[4] = t
            profile.t[5] = 0
            profile.t[6] = -(af / jMax) + t

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDUD, .ACC0_VEL) {
                return true
            }
        }

        return false
    }

    func timeVel(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) -> Bool {

        // Debug logging disabled - see UNIT_TEST_DEBUG.md for Session 3 findings
        let debugTimeVel = false

        let tzMin = max(0.0, -a0 / jMax)
        let tzMax = min((tf - a0 / jMax) / 2, (aMax - a0) / jMax)

        /// Profile UDDU
        if abs(v0) < .ulpOfOne && abs(a0) < .ulpOfOne && abs(vf) < .ulpOfOne && abs(af) < .ulpOfOne {

            var polynom: [Double] = [0, 0, 0, 0]
            polynom[0] = 1
            polynom[1] = -tf / 2
            polynom[2] = 0
            polynom[3] = pd / (2 * jMax)

            let roots = solveCubic(polynom[0], polynom[1], polynom[2], polynom[3]).sorted()

            for t0 in roots {
                if t0 > tf / 4 {
                    continue
                }
                var t = t0

                /// Single Newton step (regarding pd)
                if t > .ulpOfOne {
                    let orig = -pd + jMax * t * t * (tf - 2 * t)
                    let deriv = 2 * jMax * t * (tf - 3 * t)
                    t -= orig / deriv
                }

                profile.t[0] = t
                profile.t[1] = 0
                profile.t[2] = t
                profile.t[3] = tf - 4 * t
                profile.t[4] = t
                profile.t[5] = 0
                profile.t[6] = t

                if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .VEL) {
                    return true
                }
            }

        } else {
            if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                print("  [timeVel] Taking ELSE path (v0=\(v0), vf=\(vf), a0=\(a0), af=\(af))")
                print("  [timeVel] tzMin=\(tzMin), tzMax=\(tzMax)")
            }

            let p1 = afP2 - 2 * jMax * (-2 * af * tf + jMax * tfP2 + 3 * vd)
            let ph1 = afP3 - 3 * jMaxP2 * g1 - 3 * af * jMax * vd
            let ph2 =
                afP4 + 8 * afP3 * jMax * tf + 12 * jMax
                * (3 * jMax * vdP2 - afP2 * vd + 2 * af * jMax * (g1 - tf * vd) - 2 * jMaxP2 * tf * g1)
            let ph3 = a0 * (af - jMax * tf)
            let ph4 = jMax * (-ad + jMax * tf)

            // Find root of 5th order polynom
            var polynom: [Double] = [0, 0, 0, 0, 0, 0]
            polynom[0] = 1.0
            polynom[1] =
                (15 * a0P2 + afP2 + 4 * af * jMax * tf - 16 * ph3 - 2 * jMax * (jMax * tfP2 + 3 * vd))
                / (4 * ph4)

            polynom[2] =
                (29 * a0P3 - 2 * afP3 - 33 * a0 * ph3 + 6 * jMaxP2 * g1 + 6 * af * jMax * vd + 6 * a0
                    * p1) / (6 * jMax * ph4)
            polynom[3] =
                (61 * a0P4 - 76 * a0P2 * ph3 - 16 * a0 * ph1 + 30 * a0P2 * p1 + ph2)
                / (24 * jMaxP2 * ph4)
            polynom[4] =
                (a0 * (7 * a0P4 - 10 * a0P2 * ph3 - 4 * a0 * ph1 + 6 * a0P2 * p1 + ph2))
                / (12 * jMaxP2 * jMax * ph4)
            polynom[5] =
                (7 * a0P6 + afP6 - 12 * a0P4 * ph3 + 48 * afP3 * jMaxP2 * g1 - 8 * a0P3 * ph1 - 72
                    * jMaxP2 * jMax * (jMax * g1 * g1 + vdP2 * vd + 2 * af * g1 * vd) - 6 * afP4 * jMax
                    * vd + 36 * afP2 * jMaxP2 * vdP2 + 9 * a0P4 * p1 + 3 * a0P2 * ph2)
                / (144 * jMaxP2 * jMaxP2 * ph4)

            var deriv = polynomialMonicDerivative(&polynom)
            let dderiv = polynomialDerivative(&deriv)

            if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                print("  [timeVel] UDDU 5th-order: p1=\(p1), ph1=\(ph1), ph2=\(ph2), ph3=\(ph3), ph4=\(ph4)")
                print("  [timeVel] UDDU 5th-order: polynom=\(polynom)")
                print("  [timeVel] UDDU 5th-order: deriv=\(deriv)")
            }

            // Solve 4th order derivative analytically
            let dExtremas = solveQuarticMonic(deriv[1], deriv[2], deriv[3], deriv[4]).sorted()

            if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                print("  [timeVel] UDDU 5th-order: Found \(dExtremas.count) extremas (sorted): \(dExtremas)")
                print("  [timeVel] UDDU 5th-order: tzMin=\(tzMin), tzMax=\(tzMax)")
            }

            if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                print("  [timeVel] Found \(dExtremas.count) extremas: \(dExtremas)")
            }

            var tzCurrent = tzMin

            func checkRootUddu(_ t0: Double) -> Bool {
                // Single Newton step (regarding pd)

                var t = t0

                if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                    print("    [checkRootUddu] t0=\(t0)")
                }

                let h1pd = sqrt((a0P2 + afP2) / (2 * jMaxP2) + (2 * a0 * t + jMax * t * t - vd) / jMax)
                let orig =
                    -pd
                    - (2 * a0P3 + 4 * afP3 + 24 * a0 * jMax * t * (af + jMax * (h1pd + t - tf)) + 6 * a0P2
                        * (af + jMax * (2 * t - tf)) + 6 * (a0P2 + afP2) * jMax * h1pd + 12 * af * jMax
                        * (jMax * t * t - vd) + 12 * jMaxP2
                        * (jMax * t * t * (h1pd + t - tf) - tf * v0 - h1pd * vd)) / (12 * jMaxP2)
                let derivNewton = -(a0 + jMax * t) * (3 * (h1pd + t) - 2 * tf + (a0 + 2 * af) / jMax)
                if !orig.isNaN && !derivNewton.isNaN && abs(derivNewton) > .ulpOfOne {
                    t -= orig / derivNewton
                }

                if t > tf || t.isNaN {
                    return false
                }

                let h1 = sqrt((a0P2 + afP2) / (2 * jMaxP2) + (t * (2 * a0 + jMax * t) - vd) / jMax)

                profile.t[0] = t
                profile.t[1] = 0
                profile.t[2] = t + a0 / jMax
                profile.t[3] = tf - 2 * (t + h1) - (a0 + af) / jMax
                profile.t[4] = h1
                profile.t[5] = 0
                profile.t[6] = h1 + af / jMax

                return profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .VEL)
            }

            for tz0 in dExtremas {
                /// for (double tz: dExtremas) {
                if tz0 >= tzMax {
                    if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                        print("  [timeVel] Skipping extrema tz=\(tz0) (>= tzMax=\(tzMax))")
                    }
                    continue
                }
                var tz = tz0
                let orig = evaluatePolynomial(deriv, tz)
                if abs(orig) > tolerance {
                    tz -= orig / evaluatePolynomial(dderiv, tz)
                }

                let valNew = evaluatePolynomial(polynom, tz)
                let dderivVal = evaluatePolynomial(dderiv, tz)

                if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                    let valCurrent = evaluatePolynomial(polynom, tzCurrent)
                    print(
                        "  [timeVel] Processing extrema tz=\(tz): valCurrent=\(valCurrent), valNew=\(valNew), product=\(valCurrent * valNew)"
                    )
                }

                if abs(valNew) < 64 * abs(dderivVal) * tolerance {
                    if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                        print("  [timeVel] Checking root at tz=\(tz) (valNew=\(valNew))")
                    }
                    if checkRootUddu(tz) {
                        if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                            print("  [timeVel] ✓ Root at tz=\(tz) PASSED validation!")
                        }
                        return true
                    } else if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                        print("  [timeVel] ✗ Root at tz=\(tz) FAILED validation")
                    }
                } else if evaluatePolynomial(polynom, tzCurrent) * valNew < 0 {
                    if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                        print("  [timeVel] Shrinking interval [\(tzCurrent), \(tz)]")
                    }
                    if checkRootUddu(shrinkInterval(&polynom, tzCurrent, tz)) {
                        if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                            print("  [timeVel] ✓ Shrunk interval PASSED validation!")
                        }
                        return true
                    } else if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                        print("  [timeVel] ✗ Shrunk interval FAILED validation")
                    }
                }
                tzCurrent = tz
            }
            let valMax = evaluatePolynomial(polynom, tzMax)
            let valCurrent = evaluatePolynomial(polynom, tzCurrent)

            if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                print("  [timeVel] Final check: tzCurrent=\(tzCurrent), tzMax=\(tzMax)")
                print(
                    "  [timeVel] Final check: valCurrent=\(valCurrent), valMax=\(valMax), product=\(valCurrent * valMax)"
                )
            }

            if valCurrent * valMax < 0 {
                if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                    print("  [timeVel] Final check: Trying shrinkInterval [\(tzCurrent), \(tzMax)]")
                }
                if checkRootUddu(shrinkInterval(&polynom, tzCurrent, tzMax)) {
                    if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                        print("  [timeVel] ✓ Final shrinkInterval SUCCEEDED!")
                    }
                    return true
                } else if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                    print("  [timeVel] ✗ Final shrinkInterval FAILED validation")
                }
            } else if abs(valMax) < 8 * .ulpOfOne {
                if checkRootUddu(tzMax) {
                    return true
                }
            }
        }

        /// Profile UDUD
        if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
            print("  [timeVel] Trying UDUD profile...")
        }

        let ph1 = afP2 - 2 * jMax * (2 * af * tf + jMax * tfP2 - 3 * vd)
        let ph2 = afP3 - 3 * jMaxP2 * g1 + 3 * af * jMax * vd
        let ph3 = 2 * jMax * tf * g1 + 3 * vdP2
        let ph4 =
            afP4 - 8 * afP3 * jMax * tf + 12 * jMax
            * (jMax * ph3 + afP2 * vd + 2 * af * jMax * (g1 - tf * vd))
        let ph5 = af + jMax * tf

        if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
            print("  [timeVel] UDUD: ph1=\(ph1), ph2=\(ph2), ph3=\(ph3), ph4=\(ph4), ph5=\(ph5)")
        }

        // Find root of 6th order polynom
        var polynom: [Double] = [0, 0, 0, 0, 0, 0, 0]
        polynom[0] = 1.0
        polynom[1] = (5 * a0 - ph5) / jMax
        polynom[2] = (39 * a0P2 - ph1 - 16 * a0 * ph5) / (4 * jMaxP2)
        polynom[3] = (55 * a0P3 - 33 * a0P2 * ph5 - 6 * a0 * ph1 + 2 * ph2) / (6 * jMaxP2 * jMax)
        polynom[4] =
            (101 * a0P4 + ph4 - 76 * a0P3 * ph5 - 30 * a0P2 * ph1 + 16 * a0 * ph2)
            / (24 * jMaxP2 * jMaxP2)
        polynom[5] =
            (a0 * (11 * a0P4 + ph4 - 10 * a0P3 * ph5 - 6 * a0P2 * ph1 + 4 * a0 * ph2))
            / (12 * jMaxP2 * jMaxP2 * jMax)
        polynom[6] =
            (11 * a0P6 - afP6 - 12 * a0P5 * ph5 - 48 * afP3 * jMaxP2 * g1 - 9 * a0P4 * ph1 + 72
                * jMaxP2 * jMax * (jMax * g1 * g1 - vdP2 * vd - 2 * af * g1 * vd) - 6 * afP4 * jMax
                * vd - 36 * afP2 * jMaxP2 * vdP2 + 8 * a0P3 * ph2 + 3 * a0P2 * ph4)
            / (144 * jMaxP2 * jMaxP2 * jMaxP2)

        if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
            print("  [timeVel] UDUD: polynom=\(polynom)")
        }

        var deriv = polynomialMonicDerivative(&polynom)
        var dderiv = polynomialMonicDerivative(&deriv)

        if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
            print("  [timeVel] UDUD: deriv=\(deriv)")
            print("  [timeVel] UDUD: dderiv=\(dderiv)")
        }

        var ddTzCurrent = tzMin
        var ddTzIntervals: [(Double, Double)] = []

        let ddExtremas = solveQuarticMonic(dderiv[1], dderiv[2], dderiv[3], dderiv[4]).sorted()

        if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
            print("  [timeVel] UDUD: Found \(ddExtremas.count) ddExtremas: \(ddExtremas)")
        }

        for tz0 in ddExtremas {
            if tz0 >= tzMax {
                continue
            }
            var tz = tz0

            let orig = evaluatePolynomial(dderiv, tz)
            if abs(orig) > tolerance {
                tz -= orig / evaluatePolynomial(polynomialDerivative(&dderiv), tz)
            }

            let evalCurrent = evaluatePolynomial(deriv, ddTzCurrent)
            let evalTz = evaluatePolynomial(deriv, tz)
            let product = evalCurrent * evalTz

            if product < 0 {
                ddTzIntervals.append((ddTzCurrent, tz))
            }
            ddTzCurrent = tz
        }

        let finalEvalCurrent = evaluatePolynomial(deriv, ddTzCurrent)
        let finalEvalMax = evaluatePolynomial(deriv, tzMax)
        let finalDdProduct = finalEvalCurrent * finalEvalMax

        if finalDdProduct < 0 {
            ddTzIntervals.append((ddTzCurrent, tzMax))
        }

        if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
            print("  [timeVel] UDUD: Found \(ddTzIntervals.count) intervals to check: \(ddTzIntervals)")
        }

        var tzCurrent = tzMin

        func checkRootUdud(_ t0: Double) -> Bool {
            /// Double Newton step (regarding pd)
            var t = t0
            var h1pd = sqrt((afP2 - a0P2) / (2 * jMaxP2) - ((2 * a0 + jMax * t) * t - vd) / jMax)
            var orig =
                -pd + (afP3 - a0P3 + 3 * a0P2 * jMax * (tf - 2 * t)) / (6 * jMaxP2)
                + (2 * a0 + jMax * t) * t * (tf - t) + (jMax * h1pd - af) * h1pd * h1pd + tf * v0
            var derivNewton =
                (a0 + jMax * t) * (2 * (af + jMax * tf) - 3 * jMax * (h1pd + t) - a0) / jMax

            t -= orig / derivNewton

            h1pd = sqrt((afP2 - a0P2) / (2 * jMaxP2) - ((2 * a0 + jMax * t) * t - vd) / jMax)
            orig =
                -pd + (afP3 - a0P3 + 3 * a0P2 * jMax * (tf - 2 * t)) / (6 * jMaxP2)
                + (2 * a0 + jMax * t) * t * (tf - t) + (jMax * h1pd - af) * h1pd * h1pd + tf * v0
            if abs(orig) > 1e-9 {
                derivNewton = (a0 + jMax * t) * (2 * (af + jMax * tf) - 3 * jMax * (h1pd + t) - a0) / jMax

                t -= orig / derivNewton
            }

            let h1 = sqrt((afP2 - a0P2) / (2 * jMaxP2) - ((2 * a0 + jMax * t) * t - vd) / jMax)

            // CRITICAL: Calculate times first, validate BEFORE assigning to prevent corrupting shared profile
            let t0 = t
            let t1 = 0.0
            let t2 = t + a0 / jMax
            let t3 = tf - 2 * (t + h1) + ad / jMax
            let t4 = h1
            let t5 = 0.0
            let t6 = h1 - af / jMax

            // Validate all times are non-negative before modifying the profile
            if t0 < 0 || t2 < 0 || t3 < 0 || t4 < 0 || t6 < 0 {
                return false
            }

            profile.t[0] = t0
            profile.t[1] = t1
            profile.t[2] = t2
            profile.t[3] = t3
            profile.t[4] = t4
            profile.t[5] = t5
            profile.t[6] = t6

            return profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDUD, .VEL)
        }

        for interval in ddTzIntervals {

            if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                print("  [timeVel] UDUD: Processing interval [\(interval.0), \(interval.1)]")
            }

            let tz = shrinkInterval(&deriv, interval.0, interval.1)

            if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                print("  [timeVel] UDUD: Shrunk to tz=\(tz), tzMax=\(tzMax)")
            }

            if tz >= tzMax {
                if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                    print(
                        "  [timeVel] UDUD: Skipping interval [\(interval.0), \(interval.1)] -> tz=\(tz) >= tzMax=\(tzMax)"
                    )
                }
                continue
            }

            let pVal = evaluatePolynomial(polynom, tz)
            let dderivVal = evaluatePolynomial(dderiv, tz)
            let threshold = 64 * abs(dderivVal) * tolerance

            if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                print("  [timeVel] UDUD: pVal=\(pVal), threshold=\(threshold), abs(pVal)=\(abs(pVal))")
                print("  [timeVel] UDUD: First condition: abs(pVal) < threshold? \(abs(pVal) < threshold)")
            }

            if abs(pVal) < threshold {
                if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                    print("  [timeVel] UDUD: Checking root at tz=\(tz) (pVal=\(pVal))")
                }
                if checkRootUdud(tz) {
                    if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                        print("  [timeVel] UDUD: ✓ Root at tz=\(tz) PASSED validation!")
                    }
                    return true
                } else if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                    print("  [timeVel] UDUD: ✗ Root at tz=\(tz) FAILED validation")
                }

            } else {
                let valCurrent = evaluatePolynomial(polynom, tzCurrent)
                let product = valCurrent * pVal

                if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                    print("  [timeVel] UDUD: valCurrent=\(valCurrent), product=\(product)")
                    print("  [timeVel] UDUD: Second condition: product < 0? \(product < 0)")
                }

                if product < 0 {
                    if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                        print("  [timeVel] UDUD: Shrinking interval [\(tzCurrent), \(tz)]")
                    }
                    if checkRootUdud(shrinkInterval(&polynom, tzCurrent, tz)) {
                        if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                            print("  [timeVel] UDUD: ✓ Shrunk interval PASSED validation!")
                        }
                        return true
                    } else if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                        print("  [timeVel] UDUD: ✗ Shrunk interval FAILED validation")
                    }
                }
            }
            tzCurrent = tz
        }

        let finalUdudValCurrent = evaluatePolynomial(polynom, tzCurrent)
        let finalUdudValMax = evaluatePolynomial(polynom, tzMax)
        let finalUdudProduct = finalUdudValCurrent * finalUdudValMax

        if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
            print("  [timeVel] UDUD: Final check - tzCurrent=\(tzCurrent), tzMax=\(tzMax)")
            print(
                "  [timeVel] UDUD: Final check - valCurrent=\(finalUdudValCurrent), valMax=\(finalUdudValMax), product=\(finalUdudProduct)"
            )
            print("  [timeVel] UDUD: Final condition: product < 0? \(finalUdudProduct < 0)")
        }

        if finalUdudProduct < 0 {
            if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                print("  [timeVel] UDUD: Checking final interval [\(tzCurrent), \(tzMax)]")
            }
            if checkRootUdud(shrinkInterval(&polynom, tzCurrent, tzMax)) {
                if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                    print("  [timeVel] UDUD: ✓ Final interval PASSED validation!")
                }
                return true
            } else if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
                print("  [timeVel] UDUD: ✗ Final interval FAILED validation")
            }
        }

        if debugTimeVel && abs(pd - (-3.0)) < 0.01 && abs(tf - 4.5) < 0.01 {
            print("  [timeVel] UDUD section failed, returning false")
        }

        return false
    }

    func timeAcc0Acc1(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) -> Bool {

        if abs(a0) < .ulpOfOne && abs(af) < .ulpOfOne {
            let h1 = 2 * aMin * g1 + vdP2 + aMax * (2 * pd + aMin * tfP2 - 2 * tf * vf)
            let h2 = ((aMax - aMin) * (-aMin * vd + aMax * (aMin * tf - vd)))

            let jf = h2 / h1
            profile.t[0] = aMax / jf
            profile.t[1] = (-2 * aMax * h1 + aMin * aMin * g2) / h2
            profile.t[2] = profile.t[0]
            profile.t[3] = 0
            profile.t[4] = -aMin / jf
            profile.t[5] = tf - (2 * profile.t[0] + profile.t[1] + 2 * profile.t[4])
            profile.t[6] = profile.t[4]

            return profile.checkWithTiming(tf, jf, vMax, vMin, aMax, aMin, jMax, .UDDU, .ACC0_ACC1)
        }

        /// UDDU

        let h1 = sqrt(
            144
                * pow(
                    (aMax - aMin) * (-aMin * vd + aMax * (aMin * tf - vd)) - afP2 * (aMax * tf - vd) + 2 * af
                        * aMin * (aMax * tf - vd) + a0P2 * (aMin * tf + v0 - vf) - 2 * a0 * aMax
                        * (aMin * tf - vd),
                    2
                ) + 48 * ad
                * (3 * a0P3 - 3 * afP3 + 12 * aMax * aMin * (-aMax + aMin) + 4 * afP2 * (aMax + 2 * aMin)
                    + a0
                    * (-3 * afP2 + 8 * af * (aMin - aMax) + 6 * (aMax * aMax + 2 * aMax * aMin - aMin * aMin))
                    + 6 * af * (aMax * aMax - 2 * aMax * aMin - aMin * aMin) + a0P2
                    * (3 * af - 4 * (2 * aMax + aMin)))
                * (2 * aMin * g1 + vd * vd + aMax * (2 * pd + aMin * tf * tf - 2 * tf * vf))
        )

        let jf =
            -(3 * afP2 * aMax * tf - 3 * a0P2 * aMin * tf - 6 * ad * aMax * aMin * tf + 3 * aMax * aMin
            * (aMin - aMax) * tf + 3 * (a0P2 - afP2) * vd + 6 * vd * (af * aMin - a0 * aMax) + 3
            * (aMax * aMax - aMin * aMin) * vd + h1 / 4)
            / (6 * (2 * aMin * g1 + vd * vd + aMax * (2 * pd + aMin * tfP2 - 2 * tf * vf)))
        profile.t[0] = (aMax - a0) / jf
        profile.t[1] =
            (a0P2 - afP2 + 2 * ad * aMin - 2
                * (aMax * aMax - 2 * aMax * aMin + aMin * aMin + aMin * jf * tf - jf * vd))
            / (2 * (aMax - aMin) * jf)
        profile.t[2] = aMax / jf
        profile.t[3] = 0
        profile.t[4] = -aMin / jf
        profile.t[5] = tf - (profile.t[0] + profile.t[1] + profile.t[2] + 2 * profile.t[4] + af / jf)
        profile.t[6] = profile.t[4] + af / jf

        if profile.checkWithTiming(tf, jf, vMax, vMin, aMax, aMin, jMax, .UDDU, .ACC0_ACC1) {
            return true
        }

        return false
    }

    func timeAcc1(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) -> Bool {
        // a3 != 0
        /// Case UDDU

        let h0uddu =
            sqrt(
                jMaxP2
                    * (a0P4 + afP4 - 4 * afP3 * jMax * tf + 6 * afP2 * jMaxP2 * tfP2 - 4 * a0P3
                        * (af - jMax * tf) + 6 * a0P2 * (af - jMax * tf) * (af - jMax * tf) + 24 * af
                        * jMaxP2 * g1 - 4 * a0
                        * (afP3 - 3 * afP2 * jMax * tf + 6 * jMaxP2 * (-pd + tf * vf)) - 12 * jMaxP2
                        * (-vdP2 + jMax * tf * g2)) / 3
            ) / jMax
        let h1UDDUSol1 = sqrt(
            (a0P2 + afP2 - 2 * a0 * af - 2 * ad * jMax * tf + 2 * h0uddu) / jMaxP2 + tfP2
        )

        profile.t[0] =
            -(a0P2 + afP2 + 2 * a0 * (jMax * tf - af) - 2 * jMax * vd + h0uddu)
            / (2 * jMax * (-ad + jMax * tf))
        profile.t[1] = 0
        profile.t[2] = (tf - h1UDDUSol1) / 2 - ad / (2 * jMax)
        profile.t[3] = 0
        profile.t[4] = 0
        profile.t[5] = h1UDDUSol1
        profile.t[6] = tf - (profile.t[0] + profile.t[2] + profile.t[5])

        if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC1) {
            return true
        }

        // Case UDUD

        let h0udud =
            sqrt(
                jMaxP2
                    * (a0P4 + afP4 + 4 * (afP3 - a0P3) * jMax * tf + 6 * afP2 * jMaxP2 * tfP2 + 6
                        * a0P2 * (af + jMax * tf) * (af + jMax * tf) + 24 * af * jMaxP2 * g1 - 4 * a0
                        * (a0P2 * af + afP3 + 3 * afP2 * jMax * tf + 6 * jMaxP2 * (-pd + tf * vf)) + 12
                        * jMaxP2 * (vdP2 + jMax * tf * g2)) / 3
            ) / jMax
        let h1udud = sqrt(
            (a0P2 + afP2 - 2 * a0 * af + 2 * ad * jMax * tf + 2 * h0udud) / jMaxP2 + tfP2
        )

        profile.t[0] = 0
        profile.t[1] = 0
        profile.t[2] =
            -(a0P2 + afP2 - 2 * a0 * af + 2 * jMax * (vd - a0 * tf) + h0udud)
            / (2 * jMax * (ad + jMax * tf))
        profile.t[3] = 0
        profile.t[4] = ad / (2 * jMax) + (tf - h1udud) / 2
        profile.t[5] = h1udud
        profile.t[6] = tf - (profile.t[5] + profile.t[4] + profile.t[2])

        if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDUD, .ACC1) {
            return true
        }

        // Case UDDU, Solution 2

        let h0a =
            a0P3 - afP3 - 3 * a0P2 * aMin + 3 * aMin * aMin * (a0 + jMax * tf) + 3 * af * aMin
            * (-aMin - 2 * jMax * tf) - 3 * afP2 * (-aMin - jMax * tf) - 3 * jMaxP2
            * (-2 * pd - aMin * tfP2 + 2 * tf * vf)
        let h0b = a0P2 + afP2 - 2 * (a0 + af) * aMin + 2 * (aMin * aMin - jMax * (-aMin * tf + vd))
        let h0c =
            a0P4 + 3 * afP4 - 4 * (a0P3 + 2 * afP3) * aMin + 6 * a0P2 * aMin * aMin + 6 * afP2
            * (aMin * aMin - 2 * jMax * vd) + 12 * jMax
            * (2 * aMin * jMax * g1 - aMin * aMin * vd + jMax * vdP2) + 24 * af * aMin * jMax * vd - 4
            * a0
            * (afP3 - 3 * af * aMin * (-aMin - 2 * jMax * tf) + 3 * afP2 * (-aMin - jMax * tf) + 3
                * jMax * (-aMin * aMin * tf + jMax * (-2 * pd - aMin * tfP2 + 2 * tf * vf)))
        let h1UDDUSol2: Double = abs(jMax) / jMax * sqrt(4 * h0a * h0a - 6 * h0b * h0c)
        let h2 = 6 * jMax * h0b

        profile.t[0] = 0
        profile.t[1] = 0
        profile.t[2] = (2 * h0a + h1UDDUSol2) / h2
        profile.t[3] =
            -(a0P2 + afP2 - 2 * (a0 + af) * aMin + 2 * (aMin * aMin + aMin * jMax * tf - jMax * vd))
            / (2 * jMax * (a0 - aMin - jMax * profile.t[2]))
        profile.t[4] = (a0 - aMin) / jMax - profile.t[2]
        profile.t[5] = tf - (profile.t[2] + profile.t[3] + profile.t[4] + (af - aMin) / jMax)
        profile.t[6] = (af - aMin) / jMax

        if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC1) {
            return true
        }

        // Case UDUD, Solution 1

        let h0aUdud =
            -a0P3 + afP3 + 3 * (a0P2 - afP2) * aMax - 3 * ad * aMax * aMax - 6 * af * aMax * jMax * tf
            + 3 * afP2 * jMax * tf + 3 * jMax
            * (aMax * aMax * tf + jMax * (-2 * pd - aMax * tfP2 + 2 * tf * vf))
        let h0bUdud = a0P2 - afP2 + 2 * ad * aMax + 2 * jMax * (aMax * tf - vd)
        let h0cUdud =
            a0P4 + 3 * afP4 - 4 * (a0P3 + 2 * afP3) * aMax + 6 * a0P2 * aMax * aMax - 24 * af * aMax
            * jMax * vd + 12 * jMax * (2 * aMax * jMax * g1 + jMax * vdP2 + aMax * aMax * vd) + 6 * afP2
            * (aMax * aMax + 2 * jMax * vd) - 4 * a0
            * (afP3 + 3 * af * aMax * (aMax - 2 * jMax * tf) - 3 * afP2 * (aMax - jMax * tf) + 3 * jMax
                * (aMax * aMax * tf + jMax * (-2 * pd - aMax * tfP2 + 2 * tf * vf)))
        let h1 = abs(jMax) / jMax * sqrt(4 * h0aUdud * h0aUdud - 6 * h0bUdud * h0cUdud)
        let h2Udud = 6 * jMax * h0bUdud

        profile.t[0] = 0
        profile.t[1] = 0
        profile.t[2] = -(2 * h0aUdud + h1) / h2Udud
        profile.t[3] = 2 * h1 / h2Udud
        profile.t[4] = (aMax - a0) / jMax + profile.t[2]
        profile.t[5] = tf - (profile.t[2] + profile.t[3] + profile.t[4] + (-af + aMax) / jMax)
        profile.t[6] = (-af + aMax) / jMax

        if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDUD, .ACC1) {
            return true
        }

        return false
    }

    func timeAcc0(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) -> Bool {
        // UDUD
        do {
            let h1 = sqrt(
                self.adP2 / (2 * self.jMaxP2) - self.ad * (aMax - self.a0) / (self.jMaxP2)
                    + (aMax * self.tf - self.vd) / jMax
            )

            profile.t[0] = (aMax - self.a0) / jMax
            profile.t[1] = self.tf - self.ad / jMax - 2 * h1
            profile.t[2] = h1
            profile.t[3] = 0
            profile.t[4] = (self.af - aMax) / jMax + h1
            profile.t[5] = 0
            profile.t[6] = 0

            if profile.checkWithTiming(self.tf, jMax, vMax, vMin, aMax, aMin, .UDUD, .NONE) {
                return true
            }
        }

        // UDUD
        do {
            let h0a = -a0P2 + afP2 - 2 * ad * aMax + 2 * jMax * (aMax * tf - vd)
            let h0b =
                a0P3 + 2 * afP3 - 6 * afP2 * aMax - 3 * a0P2 * (af - jMax * tf) - 3 * a0 * aMax
                * (aMax - 2 * af + 2 * jMax * tf) - 3 * jMax
                * (jMax * (-2 * pd + aMax * tfP2 + 2 * tf * v0) + aMax * (aMax * tf - 2 * vd)) + 3 * af
                * (aMax * aMax + 2 * aMax * jMax * tf - 2 * jMax * vd)
            let h0 = abs(jMax) * sqrt(4 * h0b * h0b - 18 * h0a * h0a * h0a)
            let h1 = 3 * jMax * h0a

            profile.t[0] = (-a0 + aMax) / jMax
            profile.t[1] =
                (-a0P3 + afP3 + afP2 * (-6 * aMax + 3 * jMax * tf) + a0P2
                    * (-3 * af + 6 * aMax + 3 * jMax * tf) + 6 * af * (aMax * aMax - jMax * vd) + 3 * a0
                    * (afP2 - 2 * (aMax * aMax + jMax * vd)) - 6 * jMax
                    * (aMax * (aMax * tf - 2 * vd) + jMax * g2)) / h1
            profile.t[2] = -(ad + h0 / h1) / (2 * jMax) + tf / 2 - profile.t[1] / 2
            profile.t[3] = h0 / (jMax * h1)
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = tf - (profile.t[0] + profile.t[1] + profile.t[2] + profile.t[3])

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                return true
            }
        }

        // a3 != 0

        // UDDU Solution 1
        do {
            let h0a =
                a0P3 + 2 * afP3 - 6 * (afP2 + aMax * aMax) * aMax - 6 * (a0 + af) * aMax * jMax * tf + 9
                * aMax * aMax * (af + jMax * tf) + 3 * a0 * aMax * (-2 * af + 3 * aMax) + 3 * a0P2
                * (af - 2 * aMax + jMax * tf) - 6 * jMaxP2 * g1 + 6 * (af - aMax) * jMax * vd - 3 * aMax
                * jMaxP2 * tfP2
            let h0b = a0P2 + afP2 + 2 * (aMax * aMax - (a0 + af) * aMax + jMax * (vd - aMax * tf))
            let h1 = abs(jMax) / jMax * sqrt(4 * h0a * h0a - 18 * h0b * h0b * h0b)
            let h2 = 6 * jMax * h0b

            profile.t[0] = (-a0 + aMax) / jMax
            profile.t[1] = ad / jMax - 2 * profile.t[0] - (2 * h0a - h1) / h2 + tf
            profile.t[2] = -(2 * h0a + h1) / h2
            profile.t[3] = (2 * h0a - h1) / h2
            profile.t[4] = tf - (profile.t[0] + profile.t[1] + profile.t[2] + profile.t[3])
            profile.t[5] = 0
            profile.t[6] = 0

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .ACC0) {
                return true
            }
        }

        return false
    }

    func timeNone(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) -> Bool {
        if abs(v0) < .ulpOfOne && abs(a0) < .ulpOfOne && abs(af) < .ulpOfOne {
            let h1 = sqrt(tfP2 * vfP2 + pow2(4 * pd - tf * vf))
            let jf = 4 * (4 * pd - 2 * tf * vf + h1) / tfP3

            profile.t[0] = tf / 4
            profile.t[1] = 0
            profile.t[2] = 2 * profile.t[0]
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = profile.t[0]

            if profile.checkWithTiming(tf, jf, vMax, vMin, aMax, aMin, jMax, .UDDU, .NONE) {
                return true
            }
        }

        if abs(a0) < .ulpOfOne && abs(af) < .ulpOfOne {
            // Solution 1
            // {
            //     let h1 = sqrt(16 * pd * (pd - tf * (v0 + vf)) + tfP2 * (5 * v0P2 + 6 * v0 * vf + 5 * vfP2))
            //     let jf = 4 * (4 * pd - 2 * tf * (v0 + vf) - h1)/tfP3

            //     profile.t[0] = (tf * (v0 + 3 * vf) - 4 * pd)/(4 * vd)
            //     profile.t[1] = 0
            //     profile.t[2] = tf/2
            //     profile.t[3] = 0
            //     profile.t[4] = 0
            //     profile.t[5] = 0
            //     profile.t[6] = profile.t[4]

            //     if (profile.checkWithTiming(tf, jf, vMax, vMin, aMax, aMin, jMax, .UDDU, .NONE)) {
            //         std::cout << "i2" << std::endl
            //         return true
            //     }
            // }

            // Is that really needed?
            // Profiles with a3 != 0, Solution UDDU
            do {
                // First acc, then constant
                do {
                    var polynom = [
                        -2 * self.tf,
                        2 * self.vd / jMax + self.tfP2,
                        4 * (self.pd - self.tf * self.vf) / jMax,
                        (self.vdP2 + jMax * self.tf * self.g2) / (self.jMaxP2),
                    ]
                    let roots = solveQuarticMonic(&polynom).sorted()
                    for var t in roots {
                        if t > self.tf / 2 || t > (aMax - self.a0) / jMax {
                            continue
                        }

                        // Single Newton step (regarding pd)
                        do {
                            let h1 = (jMax * t * (t - self.tf) + self.vd) / (jMax * (2 * t - self.tf))
                            let h2 =
                                (2 * jMax * t * (t - self.tf) + jMax * self.tfP2 - 2 * self.vd)
                                / (jMax * (2 * t - self.tf) * (2 * t - self.tf))
                            let orig =
                                (-2 * pd + 2 * tf * v0 + h1 * h1 * jMax * (tf - 2 * t) + jMax * tf
                                    * (2 * h1 * t - t * t - (h1 - t) * tf)) / 2
                            let deriv =
                                (jMax * tf * (2 * t - tf) * (h2 - 1)) / 2 + h1 * jMax
                                * (tf - (2 * t - tf) * h2 - h1)

                            t -= orig / deriv
                        }

                        profile.t[0] = t
                        profile.t[1] = 0
                        profile.t[2] = (jMax * t * (t - tf) + vd) / (jMax * (2 * t - tf))
                        profile.t[3] = tf - 2 * t
                        profile.t[4] = t - profile.t[2]
                        profile.t[5] = 0
                        profile.t[6] = 0

                        if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                            return true
                        }
                    }
                }
            }
        }

        // UDUD T 0246
        do {
            let h0 =
                sqrt(
                    2 * jMaxP2
                        * (2
                            * pow2(
                                a0P3 - afP3 - 3 * afP2 * jMax * tf + 9 * af * jMaxP2 * tfP2 - 3 * a0P2
                                    * (af + jMax * tf) + 3 * a0 * pow2(af + jMax * tf) + 3 * jMaxP2
                                    * (8 * pd + jMax * tfP2 * tf - 8 * tf * vf)
                            ) - 3
                            * (a0P2 + afP2 - 2 * af * jMax * tf - 2 * a0 * (af + jMax * tf) - jMax
                                * (jMax * tfP2 + 4 * v0 - 4 * vf))
                            * (a0P4 + afP4 + 4 * afP3 * jMax * tf + 6 * afP2 * jMaxP2 * tfP2 - 3
                                * jMaxP2 * jMaxP2 * tfP2 * tfP2 - 4 * a0P3 * (af + jMax * tf) + 6 * a0P2
                                * pow2(af + jMax * tf) - 12 * af * jMaxP2
                                * (8 * pd + jMax * tfP2 * tf - 8 * tf * v0) + 48 * jMaxP2 * vdP2 + 48
                                * jMaxP2 * jMax * tf * g2 - 4 * a0
                                * (afP3 + 3 * afP2 * jMax * tf - 9 * af * jMaxP2 * tfP2 - 3 * jMaxP2
                                    * (8 * pd + jMax * tfP2 * tf - 8 * tf * vf))))
                ) / jMax
            let h1 =
                12 * jMax
                * (-a0P2 - afP2 + 2 * af * jMax * tf + 2 * a0 * (af + jMax * tf) + jMax
                    * (jMax * tfP2 + 4 * v0 - 4 * vf))
            let h2 =
                -4 * a0P3 + 4 * afP3 + 12 * a0P2 * af - 12 * a0 * afP2 + 48 * jMaxP2 * pd + 12
                * (a0P2 - afP2) * jMax * tf - 24 * jMaxP2 * tf * (v0 + vf) + 24 * ad * jMax * vd
            let h3 = 2 * a0P3 - 2 * afP3 - 6 * a0P2 * af + 6 * a0 * afP2

            profile.t[0] =
                (h3 - 48 * jMaxP2 * (tf * vf - pd) - 6 * (a0P2 + afP2) * jMax * tf + 12 * a0 * af
                    * jMax * tf + 6 * (a0 + 3 * af + jMax * tf) * tfP2 * jMaxP2 - h0) / h1
            profile.t[1] = 0
            profile.t[2] = (h2 + h0) / h1
            profile.t[3] = 0
            profile.t[4] = (-h2 + h0) / h1
            profile.t[5] = 0
            profile.t[6] =
                (-h3 + 48 * jMaxP2 * (tf * v0 - pd) - 6 * (a0P2 + afP2) * jMax * tf + 12 * a0 * af
                    * jMax * tf + 6 * (af + 3 * a0 + jMax * tf) * tfP2 * jMaxP2 - h0) / h1

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDUD, .NONE) {
                return true
            }
        }

        // Profiles with a3 != 0, Solution UDDU
        do {
            // T 0234
            do {
                let ph1 = af + jMax * tf

                var polynom = [
                    -2 * (ad + jMax * tf) / jMax,
                    2 * (a0P2 + afP2 + jMax * (af * tf + vd) - 2 * a0 * ph1) / jMaxP2 + tfP2,
                    2
                        * (a0P3 - afP3 - 3 * afP2 * jMax * tf + 3 * a0 * ph1 * (ph1 - a0) - 6 * jMaxP2
                            * (-pd + tf * vf)) / (3 * jMaxP2 * jMax),
                    (a0P4 + afP4 + 4 * afP3 * jMax * tf - 4 * a0P3 * ph1 + 6 * a0P2 * ph1 * ph1 + 24
                        * jMaxP2 * af * g1 - 4 * a0
                        * (afP3 + 3 * afP2 * jMax * tf + 6 * jMaxP2 * (-pd + tf * vf)) + 6 * jMaxP2
                        * afP2 * tfP2 + 12 * jMaxP2 * (vdP2 + jMax * tf * g2))
                        / (12 * jMaxP2 * jMaxP2),
                ]
                let tMin = ad / jMax
                let tMax = min((aMax - a0) / jMax, (ad / jMax + tf) / 2)

                let roots = solveQuarticMonic(&polynom).sorted()
                for _t in roots {
                    var tM = _t
                    if tM < tMin || tM > tMax {
                        continue
                    }

                    // Single Newton step (regarding pd)
                    do {
                        let h0 = jMax * (2 * tM - tf) - ad
                        let h1 =
                            (adP2 - 2 * af * jMax * tM + 2 * a0 * jMax * (tM - tf) + 2 * jMax
                                * (jMax * tM * (tM - tf) + vd)) / (2 * jMax * h0)
                        let h2 =
                            (-adP2 + 2 * jMaxP2 * (tfP2 + tM * (tM - tf)) + (a0 + af) * jMax * tf - ad * h0 - 2
                                * jMax * vd) / (h0 * h0)

                        // MARK: - FIXME
                        let orig =
                            (-a0P3 + afP3 + 3 * adP2 * jMax * (h1 - tM) + 3 * ad * jMaxP2 * (h1 - tM)
                                * (h1 - tM) - 3 * a0 * af * ad + 3 * jMaxP2
                                * (a0 * tfP2 - 2 * pd + 2 * tf * v0 + h1 * h1 * jMax * (tf - 2 * tM) + jMax * tf
                                    * (2 * h1 * tM - tM * tM - (h1 - tM) * tf)))
                            / (6 * jMaxP2)
                        // MARK: - FIXME

                        let deriv =
                            (h0 * (-ad + jMax * tf) * (h2 - 1)) / (2 * jMax) + h1
                            * (-ad + jMax * (tf - h1) - h0 * h2)

                        tM -= Double(orig) / deriv
                    }

                    profile.t[0] = tM
                    profile.t[1] = 0
                    profile.t[2] =
                        (adP2 + 2 * jMax * (-a0 * tf - ad * tM + jMax * tM * (tM - tf) + vd))
                        / (2 * jMax * (-ad + jMax * (2 * tM - tf)))
                    profile.t[3] = ad / jMax + tf - 2 * tM
                    profile.t[4] = tf - (tM + profile.t[2] + profile.t[3])
                    profile.t[5] = 0
                    profile.t[6] = 0

                    if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                        return true
                    }
                }
            }

            // T 3456
            do {
                let h1 = 3 * jMax * (adP2 + 2 * jMax * (a0 * tf - vd))
                let h2 = adP2 + 2 * jMax * (a0 * tf - vd)
                let h0 =
                    sqrt(
                        4
                            * pow2(
                                2 * (a0P3 - afP3) - 6 * a0P2 * (af - jMax * tf) + 6 * jMaxP2 * g1 + 3 * a0
                                    * (2 * afP2 - 2 * jMax * af * tf + jMaxP2 * tfP2) + 6 * ad * jMax * vd
                            ) - 18
                            * h2 * h2 * h2
                    ) / h1 * abs(jMax) / jMax

                profile.t[0] = 0
                profile.t[1] = 0
                profile.t[2] = 0
                profile.t[3] =
                    (afP3 - a0P3 + 3 * (afP2 - a0P2) * jMax * tf - 3 * ad * (a0 * af + 2 * jMax * vd) - 6
                        * jMaxP2 * g2) / h1
                profile.t[4] = (tf - profile.t[3] - h0) / 2 - ad / (2 * jMax)
                profile.t[5] = h0
                profile.t[6] = (tf - profile.t[3] + ad / jMax - h0) / 2

                if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                    return true
                }
            }

            // T 2346
            do {
                let ph1 = adP2 + 2 * (af + a0) * jMax * tf - jMax * (jMax * tfP2 + 4 * vd)
                let ph2 = jMax * tfP2 * g1 - vd * (-2 * pd - tf * v0 + 3 * tf * vf)
                let ph3 = 5 * afP2 - 8 * af * jMax * tf + 2 * jMax * (2 * jMax * tfP2 - vd)
                let ph4 = jMaxP2 * tfP4 - 2 * vdP2 + 8 * jMax * tf * (-pd + tf * vf)
                let ph5 =
                    (5 * afP4 - 8 * afP3 * jMax * tf - 12 * afP2 * jMax * (jMax * tfP2 + vd) + 24 * af
                        * jMaxP2 * (-2 * pd + jMax * tfP3 + 2 * tf * vf) - 6 * jMaxP2 * ph4)
                let ph6 = -vdP2 + jMax * tf * (-2 * pd + 3 * tf * v0 - tf * vf) - af * g2

                let poly1 =
                    -(4 * (a0P3 - afP3) - 12 * a0P2 * (af - jMax * tf) + 6 * a0
                    * (2 * afP2 - 2 * af * jMax * tf + jMax * (jMax * tfP2 - 2 * vd)) + 6 * af * jMax
                    * (3 * jMax * tfP2 + 2 * vd) - 6 * jMaxP2
                    * (-4 * pd + jMax * tfP3 - 2 * tf * v0 + 6 * tf * vf)) / (3 * jMax * ph1)

                let poly2 =
                    -(-a0P4 - afP4 + 4 * a0P3 * (af - jMax * tf) + a0P2
                    * (-6 * afP2 + 8 * af * jMax * tf - 4 * jMax * (jMax * tfP2 - vd)) + 2 * afP2 * jMax
                    * (jMax * tfP2 + 2 * vd) - 4 * af * jMaxP2
                    * (-3 * pd + jMax * tfP3 + 2 * tf * v0 + tf * vf) + jMaxP2
                    * (jMaxP2 * tfP4 - 8 * vdP2 + 4 * jMax * tf * (-3 * pd + tf * v0 + 2 * tf * vf))
                    + 2 * a0
                    * (2 * afP3 - 2 * afP2 * jMax * tf + af * jMax * (-3 * jMax * tfP2 - 4 * vd)
                        + jMaxP2 * (-6 * pd + jMax * tfP3 - 4 * tf * v0 + 10 * tf * vf)))
                    / (jMaxP2 * ph1)

                let poly3 =
                    -(a0P5 - afP5 + afP4 * jMax * tf - 5 * a0P4 * (af - jMax * tf) + 2 * a0P3 * ph3 + 4
                    * afP3 * jMax * (jMax * tfP2 + vd) + 12 * jMaxP2 * af * ph6 - 2 * a0P2
                    * (5 * afP3 - 9 * afP2 * jMax * tf - 6 * af * jMax * vd + 6 * jMaxP2
                        * (-2 * pd - tf * v0 + 3 * tf * vf)) - 12 * jMaxP2 * jMax * ph2 + a0 * ph5)
                    / (3 * jMaxP2 * jMax * ph1)

                let poly4 =
                    -(-a0P6 - afP6 + 6 * a0P5 * (af - jMax * tf) - 48 * afP3 * jMaxP2 * g1 + 72
                    * jMaxP2 * jMax * (jMax * g1 * g1 + vdP2 * vd + 2 * af * g1 * vd) - 3 * a0P4 * ph3
                    - 36 * afP2 * jMaxP2 * vdP2 + 6 * afP4 * jMax * vd + 4 * a0P3
                    * (5 * afP3 - 9 * afP2 * jMax * tf - 6 * af * jMax * vd + 6 * jMaxP2
                        * (-2 * pd - tf * v0 + 3 * tf * vf)) - 3 * a0P2 * ph5 + 6 * a0
                    * (afP5 - afP4 * jMax * tf - 4 * afP3 * jMax * (jMax * tfP2 + vd) + 12 * jMaxP2
                        * (-af * ph6 + jMax * ph2)))
                    / (18 * jMaxP2 * jMaxP2 * ph1)

                var polynom = [poly1, poly2, poly3, poly4]
                let tMax = (a0 - aMin) / jMax

                let roots = solveQuarticMonic(&polynom).sorted()
                for _t in roots {
                    var tM = _t
                    if tM > tMax {
                        continue
                    }

                    // Single Newton step (regarding pd)
                    do {
                        let h1 = adP2 / 2 + jMax * (af * tM + (jMax * tM - a0) * (tM - tf) - vd)
                        let h2 = -ad + jMax * (tf - 2 * tM)
                        let h3 = sqrt(h1)
                        let orig =
                            (afP3 - a0P3 + 3 * af * jMax * tM * (af + jMax * tM) + 3 * a0P2 * (af + jMax * tM) - 3
                                * a0 * (afP2 + 2 * af * jMax * tM + jMaxP2 * (tM * tM - tfP2)) + 3 * jMaxP2
                                * (-2 * pd + jMax * tM * (tM - tf) * tf + 2 * tf * v0)) / (6 * jMaxP2) - h3 * h3
                            * h3 / (jMax * abs(jMax)) + ((-ad - jMax * tM) * h1) / (jMaxP2)
                        let deriv =
                            (6 * jMax * h2 * h3 / abs(jMax) + 2 * (-ad - jMax * tf) * h2 - 2
                                * (3 * adP2 + af * jMax * (8 * tM - 2 * tf) + 4 * a0 * jMax * (-2 * tM + tf) + 2
                                    * jMax * (jMax * tM * (3 * tM - 2 * tf) - vd)))
                            / (4 * jMax)

                        tM -= orig / deriv
                    }

                    let h1 =
                        sqrt(2 * adP2 + 4 * jMax * (ad * tM + a0 * tf + jMax * tM * (tM - tf) - vd)) / abs(jMax)

                    // Solution 2 with aPlat
                    profile.t[0] = 0
                    profile.t[1] = 0
                    profile.t[2] = tM
                    profile.t[3] = tf - 2 * tM - ad / jMax - h1
                    profile.t[4] = h1 / 2
                    profile.t[5] = 0
                    profile.t[6] = tf - (tM + profile.t[3] + profile.t[4])

                    if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                        return true
                    }
                }
            }
        }

        // Profiles with a3 != 0, Solution UDUD
        do {
            // T 0124
            do {
                let ph0 = -2 * pd - tf * v0 + 3 * tf * vf
                let ph1 = -ad + jMax * tf
                let ph2 = jMax * tfP2 * g1 - vd * ph0
                let ph3 = 5 * afP2 + 2 * jMax * (2 * jMax * tfP2 - vd - 4 * af * tf)
                let ph4 = jMaxP2 * tfP4 - 2 * vdP2 + 8 * jMax * tf * (-pd + tf * vf)
                let ph5 =
                    (5 * afP4 - 8 * afP3 * jMax * tf - 12 * afP2 * jMax * (jMax * tfP2 + vd) + 24 * af
                        * jMaxP2 * (-2 * pd + jMax * tfP3 + 2 * tf * vf) - 6 * jMaxP2 * ph4)
                let ph6 = -vdP2 + jMax * tf * (-2 * pd + 3 * tf * v0 - tf * vf)
                let ph7 = 3 * jMaxP2 * ph1 * ph1

                let poly1 = (4 * af * tf - 2 * jMax * tfP2 - 4 * vd) / ph1
                let poly2 =
                    (-2 * (a0P4 + afP4) + 8 * afP3 * jMax * tf + 6 * afP2 * jMaxP2 * tfP2 + 8 * a0P3
                        * (af - jMax * tf) - 12 * a0P2 * (af - jMax * tf) * (af - jMax * tf) - 12 * af
                        * jMaxP2 * (-pd + jMax * tfP3 - 2 * tf * v0 + 3 * tf * vf) + 2 * a0
                        * (4 * afP3 - 12 * afP2 * jMax * tf + 9 * af * jMaxP2 * tfP2 - 3 * jMaxP2
                            * (2 * pd + jMax * tfP3 - 2 * tf * vf)) + 3 * jMaxP2
                        * (jMaxP2 * tfP4 + 4 * vdP2 - 4 * jMax * tf * (pd + tf * v0 - 2 * tf * vf))) / ph7
                let poly3 =
                    (-a0P5 + afP5 - afP4 * jMax * tf + 5 * a0P4 * (af - jMax * tf) - 2 * a0P3 * ph3 - 4
                        * afP3 * jMax * (jMax * tfP2 + vd) + 12 * afP2 * jMaxP2 * g2 - 12 * af * jMaxP2
                        * ph6 + 2 * a0P2
                        * (5 * afP3 - 9 * afP2 * jMax * tf - 6 * af * jMax * vd + 6 * jMaxP2 * ph0) + 12
                        * jMaxP2 * jMax * ph2 + a0
                        * (-5 * afP4 + 8 * afP3 * jMax * tf + 12 * afP2 * jMax * (jMax * tfP2 + vd) - 24
                            * af * jMaxP2 * (-2 * pd + jMax * tfP3 + 2 * tf * vf) + 6 * jMaxP2 * ph4))
                    / (jMax * ph7)
                let poly4 =
                    -(a0P6 + afP6 - 6 * a0P5 * (af - jMax * tf) + 48 * afP3 * jMaxP2 * g1 - 72
                    * jMaxP2 * jMax * (jMax * g1 * g1 + vdP2 * vd + 2 * af * g1 * vd) + 3 * a0P4 * ph3
                    - 6 * afP4 * jMax * vd + 36 * afP2 * jMaxP2 * vdP2 - 4 * a0P3
                    * (5 * afP3 - 9 * afP2 * jMax * tf - 6 * af * jMax * vd + 6 * jMaxP2 * ph0) + 3
                    * a0P2 * ph5 - 6 * a0
                    * (afP5 - afP4 * jMax * tf - 4 * afP3 * jMax * (jMax * tfP2 + vd) + 12 * jMaxP2
                        * (afP2 * g2 - af * ph6 + jMax * ph2)))
                    / (6 * jMaxP2 * ph7)

                var polynom = [
                    poly1,
                    poly2,
                    poly3,
                    poly4,
                ]

                let roots = solveQuarticMonic(&polynom).sorted()
                for t in roots {
                    if t > tf || t > (aMax - a0) / jMax { continue }

                    let h1 = sqrt(
                        adP2 / (2 * jMaxP2) + (a0 * (t + tf) - af * t + jMax * t * tf - vd) / jMax
                    )

                    profile.t[0] = t
                    profile.t[1] = tf - ad / jMax - 2 * h1
                    profile.t[2] = h1
                    profile.t[3] = 0
                    profile.t[4] = ad / jMax + h1 - t
                    profile.t[5] = 0
                    profile.t[6] = 0

                    if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDUD, .NONE) {
                        return true
                    }
                }
            }
        }

        // 3 step profile (ak. UZD), sometimes missed because of numerical errors T 012
        do {
            let h1 = sqrt(-adP2 + jMax * (2 * (a0 + af) * tf - 4 * vd + jMax * tfP2)) / abs(jMax)

            profile.t[0] = (tf - h1 + ad / jMax) / 2
            profile.t[1] = h1
            profile.t[2] = (tf - h1 - ad / jMax) / 2
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = 0

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                return true
            }
        }

        // 3 step profile (ak. UZU), sometimes missed because of numerical errors
        do {
            let polynom = [
                adP2,
                adP2 * tf,
                (a0P2 + afP2 + 10 * a0 * af) * tfP2 + 24 * (tf * (af * v0 - a0 * vf) - pd * ad) + 12
                    * vdP2,
                -3 * tf * ((a0P2 + afP2 + 2 * a0 * af) * tfP2 - 4 * vd * (a0 + af) * tf + 4 * vdP2),
            ]
            let roots = solveCubic(polynom[0], polynom[1], polynom[2], polynom[3]).sorted()
            for t in roots {
                if t > tf { continue }

                let jf = ad / (tf - t)

                profile.t[0] = (2 * (vd - a0 * tf) + ad * (t - tf)) / (2 * jf * t)
                profile.t[1] = t
                profile.t[2] = 0
                profile.t[3] = 0
                profile.t[4] = 0
                profile.t[5] = 0
                profile.t[6] = tf - (profile.t[0] + profile.t[1])

                if profile.checkWithTiming(tf, jf, vMax, vMin, aMax, aMin, jMax, .UDDU, .NONE) {
                    return true
                }
            }
        }

        // 3 step profile (ak. UDU), sometimes missed because of numerical errors
        do {
            profile.t[0] =
                (adP2 / jMax + 2 * (a0 + af) * tf - jMax * tfP2 - 4 * vd) / (4 * (ad - jMax * tf))
            profile.t[1] = 0
            profile.t[2] = -ad / (2 * jMax) + tf / 2
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = tf - (profile.t[0] + profile.t[2])

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                return true
            }
        }

        return false
    }

    func timeNoneSmooth(
        _ profile: inout Profile,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) -> Bool {
        do {
            let h0 = adP2 + 2 * jMax * (a0 * tf - vd)
            let h1a =
                2 * (a0P3 - afP3) - 6 * a0P2 * (af - jMax * tf) + 6 * jMaxP2 * (-pd + tf * v0) + 6
                * a0 * afP2 + 3 * a0 * jMax * (jMax * tfP2 - 2 * vd) + 6 * af * jMax * (vd - tf * a0)
            let h1 = sqrt(4 * h1a * h1a - 18 * h0 * h0 * h0) * abs(jMax) / jMax

            profile.t[0] = 0
            profile.t[1] =
                (-a0P3 + afP3 + 3 * (afP2 - a0P2) * jMax * tf - 3 * a0 * af * ad - 6 * jMax * ad * vd
                    - 6 * jMaxP2 * (-2 * pd + tf * (v0 + vf))) / (3 * jMax * h0)
            profile.t[2] =
                (4 * (a0P3 - afP3) + 6 * jMaxP2 * a0 * tfP2 + 12 * a0 * af * ad + 12 * jMax
                    * (jMax * (tf * v0 - pd) + ad * (vd - a0 * tf)) - h1) / (6 * jMax * h0)
            profile.t[3] = h1 / (3 * jMax * h0)
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = tf - (profile.t[1] + profile.t[2] + profile.t[3])

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                return true
            }
        }

        do {
            let h0 = adP2 + 2 * jMax * (vd - af * tf)
            let h0b = afP3 - 3 * jMaxP2 * (af * tfP2 + 2 * (pd - tf * vf))
            let h1a = a0P3 + 3 * a0 * af * ad - h0b
            let h1 =
                sqrt(
                    4 * h1a * h1a - 6 * h0
                        * (a0P4 + afP4 - 4 * a0P3 * af + 6 * a0P2 * afP2 + 12 * jMaxP2
                            * (vdP2 - 2 * af * (pd - tf * v0)) - 4 * a0 * h0b)
                ) * abs(jMax) / jMax

            profile.t[0] = -(2 * h1a + h1) / (6 * jMax * h0)
            profile.t[1] = h1 / (3 * jMax * h0)
            profile.t[2] = profile.t[0] - (af - a0) / jMax
            profile.t[3] = 0
            profile.t[4] = 0
            profile.t[5] = tf - (profile.t[0] + profile.t[1] + profile.t[2])
            profile.t[6] = 0

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                return true
            }
        }

        // Solution 3
        do {
            let h0 =
                sqrt(
                    3
                        * (a0P4 + afP4 - 4 * afP3 * jMax * tf + 6 * afP2 * jMaxP2 * tfP2 - 4 * a0P3
                            * (af - jMax * tf) + 6 * a0P2 * (af - jMax * tf) * (af - jMax * tf) + 24 * af
                            * jMaxP2 * (-pd + tf * v0) - 4 * a0
                            * (afP3 - 3 * afP2 * jMax * tf + 6 * jMaxP2 * (-pd + tf * vf)) - 12 * jMaxP2
                            * (-vdP2 + jMax * tf * (-2 * pd + tf * (v0 + vf))))
                ) * abs(jMax) / jMax
            let h1 =
                sqrt(
                    3
                        * (3 * a0P2 + 3 * afP2 - 6 * a0 * af - 6 * ad * jMax * tf + 3 * jMaxP2 * tfP2 - 2
                            * h0)
                ) * abs(jMax) / jMax

            profile.t[0] =
                (-3 * (a0P2 + afP2) + 6 * a0 * af + 6 * jMax * (vd - a0 * tf) + h0)
                / (6 * jMax * (-ad + jMax * tf))
            profile.t[1] = 0
            profile.t[2] = (3 * jMax * tf - 3 * ad - h1) / (6 * jMax)
            profile.t[3] = h1 / (3 * jMax)
            profile.t[4] = 0
            profile.t[5] = 0
            profile.t[6] = tf - (profile.t[0] + profile.t[2] + profile.t[3])

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                return true
            }
        }

        // Solution 2
        do {
            let h0 = 6 * (adP2 + 2 * af * jMax * tf - 2 * jMax * vd)
            let h1a =
                2
                * (a0P3 - afP3 + 3 * a0 * af * ad + 6 * jMaxP2 * (pd - tf * vf) + 3 * jMaxP2 * af
                    * tfP2)
            let h1 =
                sqrt(
                    h1a * h1a - h0
                        * (a0P4 - 4 * a0P3 * af + 6 * a0P2 * afP2 + afP4 + 24 * af * jMaxP2
                            * (-pd + tf * v0) + 12 * jMaxP2 * vdP2 - 4 * a0
                            * (afP3 - 3 * af * jMaxP2 * tfP2 + 6 * jMaxP2 * (-pd + tf * vf)))
                ) * abs(jMax)
                / jMax
            let h2 =
                4 * a0P3 - 4 * afP3 + 12 * a0 * af * ad - 12 * jMaxP2 * (pd - tf * vf) - 6 * jMaxP2
                * af * tfP2 + 12 * ad * jMax * (vd - af * tf)
            let h3 = jMax * h0

            profile.t[0] = 0
            profile.t[1] = 0
            profile.t[2] = (h1a + h1) / h3
            profile.t[3] = -(h2 + h1) / h3
            profile.t[4] = (h2 - h1) / h3
            profile.t[5] = tf - (profile.t[2] + profile.t[3] + profile.t[4])
            profile.t[6] = 0

            if profile.checkWithTiming(tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                return true
            }
        }

        // Solution 1
        do {
            // Break down complex expression for h0
            let term1 = self.a0P4 + self.afP4 - 4 * self.afP3 * jMax * self.tf
            let term2 = 6 * self.afP2 * self.jMaxP2 * self.tfP2
            let term3 = -4 * self.a0P3 * (self.af - jMax * self.tf)
            let term4 = 6 * self.a0P2 * (self.af - jMax * self.tf) * (self.af - jMax * self.tf)
            let term5 = 24 * self.af * self.jMaxP2 * (-self.pd + self.tf * self.v0)
            let term6 =
                -4 * self.a0
                * (self.afP3 - 3 * self.afP2 * jMax * self.tf + 6 * self.jMaxP2
                    * (-self.pd + self.tf * self.vf))
            let term7 =
                -12 * self.jMaxP2
                * (-self.vdP2 + jMax * self.tf * (-2 * self.pd + self.tf * (self.v0 + self.vf)))

            let h0 = sqrt((term1 + term2 + term3 + term4 + term5 + term6 + term7) / 3) * abs(jMax) / jMax
            let h1 =
                sqrt(self.adP2 - 2 * self.ad * jMax * self.tf + self.jMaxP2 * self.tfP2 + 2 * h0)
                * abs(jMax) / jMax

            profile.t[0] =
                -(self.adP2 + 2 * jMax * (self.a0 * self.tf - self.vd) + h0)
                / (2 * jMax * (-self.ad + jMax * self.tf))
            profile.t[1] = 0
            profile.t[2] = 0
            profile.t[3] = 0
            profile.t[4] = (-self.ad + jMax * self.tf - h1) / (2 * jMax)
            profile.t[5] = h1 / jMax
            profile.t[6] = self.tf - (profile.t[0] + profile.t[4] + profile.t[5])

            if profile.checkWithTiming(self.tf, jMax, vMax, vMin, aMax, aMin, .UDDU, .NONE) {
                return true
            }
        }

        return false
    }

    func getProfile(
        _ profile: inout Profile
    ) -> Bool {
        // Set the profile target values and initial state
        // These are used by check() to validate the generated trajectory
        profile.pf = pf
        profile.vf = vf
        profile.af = af
        profile.p[0] = p0
        profile.v[0] = v0
        profile.a[0] = a0

        // Test all cases to get ones that match
        // However we should guess which one is correct and try them first...
        let upFirst = (pd > tf * v0)
        let vMax = upFirst ? _vMax : _vMin
        let vMin = upFirst ? _vMin : _vMax
        let aMax = upFirst ? _aMax : _aMin
        let aMin = upFirst ? _aMin : _aMax
        let jMax = upFirst ? _jMax : -_jMax

        if timeAcc0Acc1Vel(&profile, vMax, vMin, aMax, aMin, jMax) {
            return true
        }

        if timeVel(&profile, vMax, vMin, aMax, aMin, jMax) {
            return true
        }

        if timeAcc0Vel(&profile, vMax, vMin, aMax, aMin, jMax) {
            return true
        }

        if timeAcc1Vel(&profile, vMax, vMin, aMax, aMin, jMax) {
            return true
        }
        
        if timeAcc0Acc1Vel(&profile, vMin, vMax, aMin, aMax, -jMax) {
            return true
        }

        if timeVel(&profile, vMin, vMax, aMin, aMax, -jMax) {
            return true
        }

        if timeAcc0Vel(&profile, vMin, vMax, aMin, aMax, -jMax) {
            return true
        }

        if timeAcc1Vel(&profile, vMin, vMax, aMin, aMax, -jMax) {
            return true
        }

        if timeAcc0Acc1(&profile, vMax, vMin, aMax, aMin, jMax) {
            return true
        }

        if timeAcc0(&profile, vMax, vMin, aMax, aMin, jMax) {
            return true
        }

        if timeAcc1(&profile, vMax, vMin, aMax, aMin, jMax) {
            return true
        }

        if timeNone(&profile, vMax, vMin, aMax, aMin, jMax) {
            return true
        }
        
        if timeAcc0Acc1(&profile, vMin, vMax, aMin, aMax, -jMax) {
            return true
        }

        if timeAcc0(&profile, vMin, vMax, aMin, aMax, -jMax) {
            return true
        }

        if timeAcc1(&profile, vMin, vMax, aMin, aMax, -jMax) {
            return true
        }

        if timeNone(&profile, vMin, vMax, aMin, aMax, -jMax) {
            return true
        }
        return false
    }
}
