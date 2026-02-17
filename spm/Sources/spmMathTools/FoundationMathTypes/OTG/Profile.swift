//
//  Profile.swift
//
//
//  Created by Nicholas Bergantz on 3/25/24.
//

import Foundation

private let vEps: Double = 1e-12
private let aEps: Double = 1e-12
private let jEps: Double = 1e-12
private let pPrecision: Double = 1e-8
private let vPrecision: Double = 1e-8
private let aPrecision: Double = 1e-10
private let tPrecision: Double = 1e-12
private let tMax: Double = 1e12

//! @brief A single-dof kinematic profile with position, velocity, acceleration and jerk
//!
//! The class members are only available in the Ruckig Community Version.
public struct Profile {

    public var t: [Double] = Array(repeating: 0.0, count: 7)
    var tSum: [Double] = Array(repeating: 0.0, count: 7)
    var j: [Double] = Array(repeating: 0.0, count: 7)
    var a: [Double] = Array(repeating: 0.0, count: 8)
    var v: [Double] = Array(repeating: 0.0, count: 8)
    var p: [Double] = Array(repeating: 0.0, count: 8)

    //! Brake sub-profiles
    var brake: BrakeProfile = BrakeProfile()
    var accel: BrakeProfile = BrakeProfile()

    //! Target (final) kinematic state
    var pf: Double = 0
    var vf: Double = 0
    var af: Double = 0

    var limits: ReachedLimits = .NONE
    var direction: Direction = .DOWN
    var controlSigns: ControlSigns = .UDDU

    init() {
        // Default initialization - arrays are already initialized above
    }

    // For third-order velocity interface
    mutating func checkForVelocity(
        _ jf: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {

        if t[0] < 0 { return false }

        tSum[0] = t[0]

        for i in 0..<6 {
            if t[i + 1] < 0 { return false }
            tSum[i + 1] = tSum[i] + t[i + 1]
        }

        if limits == .ACC0 && t[1] < Double.ulpOfOne {
            return false
        }

        if tSum.last! > tMax {
            return false
        }

        if controlSigns == .UDDU {
            j = [
                (t[0] > 0 ? jf : 0),
                0,
                (t[2] > 0 ? -jf : 0),
                0,
                (t[4] > 0 ? -jf : 0),
                0,
                (t[6] > 0 ? jf : 0),
            ]
        } else {
            j = [
                (t[0] > 0 ? jf : 0),
                0,
                (t[2] > 0 ? -jf : 0),
                0,
                (t[4] > 0 ? jf : 0),
                0,
                (t[6] > 0 ? -jf : 0),
            ]
        }

        for i in 0..<7 {
            a[i + 1] = a[i] + t[i] * j[i]
            v[i + 1] = v[i] + t[i] * (a[i] + t[i] * j[i] / 2)
            p[i + 1] = p[i] + t[i] * (v[i] + t[i] * (a[i] / 2 + t[i] * j[i] / 6))
        }

        self.controlSigns = controlSigns
        self.limits = limits

        direction = aMax > 0 ? .UP : .DOWN
        let aUppLim = (direction == .UP ? aMax : aMin) + aEps
        let aLowLim = (direction == .UP ? aMin : aMax) - aEps

        return abs(v.last! - vf) < vPrecision && abs(a.last! - af) < aPrecision
            && a[1] >= aLowLim && a[3] >= aLowLim && a[5] >= aLowLim
            && a[1] <= aUppLim && a[3] <= aUppLim && a[5] <= aUppLim
    }

    mutating func checkForVelocityWithTiming(
        _ tf: Double,
        _ jf: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {
        checkForVelocity(jf, aMax, aMin, controlSigns, limits)
    }

    mutating func checkForVelocityWithTiming(
        _ tf: Double,
        _ jf: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {
        (abs(jf) < abs(jMax) + jEps)
            && checkForVelocityWithTiming(tf, jf, aMax, aMin, controlSigns, limits)
    }

    mutating func setBoundaryForVelocity(
        _ p0New: Double,
        _ v0New: Double,
        _ a0New: Double,
        _ vfNew: Double,
        _ afNew: Double
    ) {
        a[0] = a0New
        v[0] = v0New
        p[0] = p0New
        af = afNew
        vf = vfNew
    }

    // For second-order velocity interface
    mutating func checkForSecondOrderVelocity(
        _ aUp: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {
        // .ACC0
        if t[1] < 0.0 { return false }

        tSum = [0, t[1], t[1], t[1], t[1], t[1], t[1]]
        if tSum.last! > tMax {
            return false
        }

        j = [0, 0, 0, 0, 0, 0, 0]
        a = [0, (t[1] > 0) ? aUp : 0, 0, 0, 0, 0, 0, af]
        for i in 0..<7 {
            v[i + 1] = v[i] + t[i] * a[i]
            p[i + 1] = p[i] + t[i] * (v[i] + t[i] * a[i] / 2)
        }

        self.controlSigns = controlSigns
        self.limits = limits

        direction = (aUp > 0) ? .UP : .DOWN

        return abs(v.last! - vf) < vPrecision
    }

    mutating func checkForSecondOrderVelocityWithTiming(
        _ tf: Double,
        _ aUp: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {
        checkForSecondOrderVelocity(aUp, controlSigns, limits)
    }

    mutating func checkForSecondOrderVelocityWithTiming(
        _ tf: Double,
        _ aUp: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {
        (aMin - aEps < aUp) && (aUp < aMax + aEps)
            && checkForSecondOrderVelocityWithTiming(tf, aUp, controlSigns, limits)
    }

    // For third-order position interface
    mutating func check(
        _ jf: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits,
        _ setLimits: Bool = false
    ) -> Bool {

        if t[0] < 0 { return false }

        tSum[0] = t[0]
        for i in 0..<6 {
            if t[i + 1] < 0 {
                return false
            }

            tSum[i + 1] = tSum[i] + t[i + 1]
        }

        if limits == .ACC0_ACC1_VEL || limits == .ACC0_VEL || limits == .ACC1_VEL || limits == .VEL {
            if t[3] < Double.ulpOfOne { return false }
        }

        if limits == .ACC0 || limits == .ACC0_ACC1 {
            if t[1] < Double.ulpOfOne { return false }
        }

        if limits == .ACC1 || limits == .ACC0_ACC1 {
            if t[5] < Double.ulpOfOne { return false }
        }

        if tSum.last! > tMax {
            return false
        }

        if controlSigns == .UDDU {
            j = [
                (t[0] > 0 ? jf : 0),
                0,
                (t[2] > 0 ? -jf : 0),
                0,
                (t[4] > 0 ? -jf : 0),
                0,
                (t[6] > 0 ? jf : 0),
            ]
        } else {
            j = [
                (t[0] > 0 ? jf : 0),
                0,
                (t[2] > 0 ? -jf : 0),
                0,
                (t[4] > 0 ? jf : 0),
                0,
                (t[6] > 0 ? -jf : 0),
            ]
        }

        direction = vMax > 0 ? .UP : .DOWN
        let vUppLim = (direction == .UP ? vMax : vMin) + vEps
        let vLowLim = (direction == .UP ? vMin : vMax) - vEps

        for i in 0..<7 {
            a[i + 1] = a[i] + t[i] * j[i]
            v[i + 1] = v[i] + t[i] * (a[i] + t[i] * j[i] / 2)
            p[i + 1] = p[i] + t[i] * (v[i] + t[i] * (a[i] / 2 + t[i] * j[i] / 6))

            if limits == .ACC0_ACC1_VEL || limits == .ACC0_ACC1 || limits == .ACC0_VEL
                || limits == .ACC1_VEL || limits == .VEL
            {
                if i == 2 {
                    a[3] = 0.0
                    if t[2] > Double.ulpOfOne {
                        j[2] = (a[3] - a[2]) / t[2]
                    }
                }
            }

            if setLimits {
                if limits == .ACC1 {
                    if i == 2 {
                        a[3] = aMin
                        if t[2] > Double.ulpOfOne {
                            j[2] = (a[3] - a[2]) / t[2]
                        }
                    }
                }

                if limits == .ACC0_ACC1 {
                    if i == 0 {
                        a[1] = aMax
                        if t[0] > Double.ulpOfOne {
                            j[0] = (a[1] - a[0]) / t[0]
                        }
                    }

                    if i == 4 {
                        a[5] = aMin
                        if t[4] > Double.ulpOfOne {
                            j[4] = (a[5] - a[4]) / t[4]
                        }
                    }
                }
            }

            if i > 1 && a[i + 1] * a[i] < -Double.ulpOfOne {
                let vAZero = v[i] - (a[i] * a[i]) / (2 * j[i])
                if vAZero > vUppLim || vAZero < vLowLim { return false }
            }
        }

        self.controlSigns = controlSigns
        self.limits = limits

        let aUppLim = (direction == .UP ? aMax : aMin) + aEps
        let aLowLim = (direction == .UP ? aMin : aMax) - aEps

        let pCheck = abs(p.last! - pf) < pPrecision
        let vCheck = abs(v.last! - vf) < vPrecision
        let aCheck = abs(a.last! - af) < aPrecision

        let aLimCheck =
            a[1] >= aLowLim && a[3] >= aLowLim && a[5] >= aLowLim
            && a[1] <= aUppLim && a[3] <= aUppLim && a[5] <= aUppLim

        let vLimCheck =
            v[3] <= vUppLim && v[4] <= vUppLim && v[5] <= vUppLim && v[6] <= vUppLim
            && v[3] >= vLowLim && v[4] >= vLowLim && v[5] >= vLowLim && v[6] >= vLowLim

        return pCheck && vCheck && aCheck && aLimCheck && vLimCheck
    }

    mutating func checkWithTiming(
        _ tf: Double,
        _ jf: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {
        check(jf, vMax, vMin, aMax, aMin, controlSigns, limits)
    }

    mutating func checkWithTiming(
        _ tf: Double,
        _ jf: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {

        if !(abs(jf) < abs(jMax) + jEps) {
            return false
        }

        return checkWithTiming(tf, jf, vMax, vMin, aMax, aMin, controlSigns, limits)
    }

    mutating func setBoundary(
        _ profile: Profile
    ) {
        a[0] = profile.a[0]
        v[0] = profile.v[0]
        p[0] = profile.p[0]
        af = profile.af
        vf = profile.vf
        pf = profile.pf
        brake = profile.brake
        accel = profile.accel
    }

    mutating func setBoundary(
        _ p0New: Double,
        _ v0New: Double,
        _ a0New: Double,
        _ pfNew: Double,
        _ vfNew: Double,
        _ afNew: Double
    ) {
        a[0] = a0New
        v[0] = v0New
        p[0] = p0New
        af = afNew
        vf = vfNew
        pf = pfNew
    }

    // For second-order position interface
    mutating func checkForSecondOrder(
        _ aUp: Double,
        _ aDown: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {
        if t[0] < 0 { return false }

        tSum[0] = t[0]
        for i in 0..<6 {
            if t[i + 1] < 0 { return false }

            tSum[i + 1] = tSum[i] + t[i + 1]
        }

        if tSum.last! > tMax {
            return false
        }

        j = [0, 0, 0, 0, 0, 0, 0]
        if controlSigns == .UDDU {
            a = [
                (t[0] > 0 ? aUp : 0), 0, (t[2] > 0 ? aDown : 0), 0, (t[4] > 0 ? aDown : 0), 0,
                (t[6] > 0 ? aUp : 0), af,
            ]
        } else {
            a = [
                (t[0] > 0 ? aUp : 0), 0, (t[2] > 0 ? aDown : 0), 0, (t[4] > 0 ? aUp : 0), 0,
                (t[6] > 0 ? aDown : 0), af,
            ]
        }

        direction = (vMax > 0) ? .UP : .DOWN
        let vUppLim = (direction == .UP ? vMax : vMin) + vEps
        let vLowLim = (direction == .UP ? vMin : vMax) - vEps

        for i in 0..<7 {
            v[i + 1] = v[i] + t[i] * a[i]
            p[i + 1] = p[i] + t[i] * (v[i] + t[i] * a[i] / 2)
        }

        self.controlSigns = controlSigns
        self.limits = limits

        return abs(p.last! - pf) < pPrecision && abs(v.last! - vf) < vPrecision
            && v[2] <= vUppLim && v[3] <= vUppLim && v[4] <= vUppLim && v[5] <= vUppLim && v[6] <= vUppLim
            && v[2] >= vLowLim && v[3] >= vLowLim && v[4] >= vLowLim && v[5] >= vLowLim && v[6] >= vLowLim
    }

    mutating func checkForSecondOrderWithTiming(
        _ tf: Double,
        _ aUp: Double,
        _ aDown: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {
        checkForSecondOrder(aUp, aDown, vMax, vMin, controlSigns, limits)
    }

    mutating func checkForSecondOrderWithTiming(
        _ tf: Double,
        _ aUp: Double,
        _ aDown: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {
        (aMin - aEps < aUp) && (aUp < aMax + aEps) && (aMin - aEps < aDown)
            && (aDown < aMax + aEps)
            && checkForSecondOrderWithTiming(tf, aUp, aDown, vMax, vMin, controlSigns, limits)
    }

    // For first-order position interface
    mutating func checkForFirstOrder(
        _ vUp: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {
        // .VEL
        if t[3] < 0.0 { return false }

        tSum = [0, 0, 0, t[3], t[3], t[3], t[3]]
        if tSum.last! > tMax {
            return false
        }

        j = [0, 0, 0, 0, 0, 0, 0]
        a = [0, 0, 0, 0, 0, 0, 0, af]
        v = [0, 0, 0, t[3] > 0 ? vUp : 0, 0, 0, 0, vf]
        for i in 0..<7 {
            p[i + 1] = p[i] + t[i] * (v[i] + t[i] * a[i] / 2)
        }

        self.controlSigns = controlSigns
        self.limits = limits

        direction = (vUp > 0) ? .UP : .DOWN

        return abs(p.last! - pf) < pPrecision
    }

    mutating func checkForFirstOrderWithTiming(
        _ tf: Double,
        _ vUp: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {

        checkForFirstOrder(vUp, controlSigns, limits)
    }

    mutating func checkForFirstOrderWithTiming(
        _ tf: Double,
        _ vUp: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ controlSigns: ControlSigns,
        _ limits: ReachedLimits
    ) -> Bool {

        (vMin - vEps < vUp) && (vUp < vMax + vEps)
            && checkForFirstOrderWithTiming(tf, vUp, controlSigns, limits)
    }

    // Secondary features
    static func checkPositionExtremum(
        _ tExt: Double,
        _ tSum: Double,
        _ t: Double,
        _ p: Double,
        _ v: Double,
        _ a: Double,
        _ j: Double,
        _ ext: inout Bound
    ) {
        if 0 < tExt && tExt < t {
            var pExt: Double
            var aExt: Double
            (pExt, _, aExt) = integrate(tExt, p, v, a, j)
            if aExt > 0 && pExt < ext.min {
                ext.min = pExt
                ext.tMin = tSum + tExt
            } else if aExt < 0 && pExt > ext.max {
                ext.max = pExt
                ext.tMax = tSum + tExt
            }
        }
    }

    static func checkStepForPositionExtremum(
        _ tSum: Double,
        _ t: Double,
        _ p: Double,
        _ v: Double,
        _ a: Double,
        _ j: Double,
        _ ext: inout Bound
    ) {
        if p < ext.min {
            ext.min = p
            ext.tMin = tSum
        }
        if p > ext.max {
            ext.max = p
            ext.tMax = tSum
        }

        if j != 0 {
            let D = a * a - 2 * j * v
            if abs(D) < Double.ulpOfOne {
                checkPositionExtremum(-a / j, tSum, t, p, v, a, j, &ext)

            } else if D > 0.0 {
                let D_sqrt = sqrt(D)
                checkPositionExtremum((-a - D_sqrt) / j, tSum, t, p, v, a, j, &ext)
                checkPositionExtremum((-a + D_sqrt) / j, tSum, t, p, v, a, j, &ext)
            }
        }
    }

    func getPositionExtrema() -> Bound {
        var extrema = Bound(
            min: Double.infinity,
            max: -Double.infinity,
            tMin: -Double.infinity,
            tMax: -Double.infinity
        )

        if brake.duration > 0.0 {
            if brake.t[0] > 0.0 {
                Profile.checkStepForPositionExtremum(
                    0.0,
                    brake.t[0],
                    brake.p[0],
                    brake.v[0],
                    brake.a[0],
                    brake.j[0],
                    &extrema
                )

                if brake.t[1] > 0.0 {
                    Profile.checkStepForPositionExtremum(
                        brake.t[0],
                        brake.t[1],
                        brake.p[1],
                        brake.v[1],
                        brake.a[1],
                        brake.j[1],
                        &extrema
                    )
                }
            }
        }

        var t_currentSum: Double = 0.0
        for i in 0..<7 {
            if i > 0 {
                t_currentSum = tSum[i - 1]
            }
            Profile.checkStepForPositionExtremum(
                t_currentSum + brake.duration,
                t[i],
                p[i],
                v[i],
                a[i],
                j[i],
                &extrema
            )
        }

        if pf < extrema.min {
            extrema.min = pf
            extrema.tMin = tSum.last! + brake.duration
        }
        if pf > extrema.max {
            extrema.max = pf
            extrema.tMax = tSum.last! + brake.duration
        }

        return extrema
    }

    func getFirstStateAtPosition(
        pt: Double,
        time: inout Double,
        timeAfter: Double = 0.0
    ) -> Bool {
        var tCum: Double = 0.0

        for i in 0..<7 {
            if t[i] == 0.0 { continue }

            if abs(p[i] - pt) < Double.ulpOfOne && tCum >= timeAfter {
                time = tCum
                return true
            }

            for _t in CubicUnivariatePolynomial(j[i] / 6, a[i] / 2, v[i], p[i] - pt).roots {
                if 0 < _t && timeAfter - tCum <= _t && _t <= t[i] {
                    time = _t + tCum
                    return true
                }
            }

            tCum += t[i]
        }

        if (t[6] > 0.0 || tSum.last! == 0.0) && abs(pf - pt) < 1e-9 && tSum.last! >= timeAfter {
            time = tSum.last!
            return true
        }

        return false
    }
}

extension Profile: CustomStringConvertible {

    public var description: String {
        var result: String = ""
        switch direction {
        case .UP: result += "UP_"
        case .DOWN: result += "DOWN_"
        }
        switch limits {
        case .ACC0_ACC1_VEL: result += "ACC0_ACC1_VEL"
        case .VEL: result += "VEL"
        case .ACC0: result += "ACC0"
        case .ACC1: result += "ACC1"
        case .ACC0_ACC1: result += "ACC0_ACC1"
        case .ACC0_VEL: result += "ACC0_VEL"
        case .ACC1_VEL: result += "ACC1_VEL"
        case .NONE: result += "NONE"
        }
        switch controlSigns {
        case .UDDU: result += "_UDDU"
        case .UDUD: result += "_UDUD"
        }
        return result
    }
}

extension Profile: Equatable {
    public static func == (lhs: Profile, rhs: Profile) -> Bool {
        lhs.t == rhs.t && lhs.tSum == rhs.tSum && lhs.j == rhs.j && lhs.a == rhs.a && lhs.v == rhs.v && lhs.p == rhs.p
            && lhs.brake == rhs.brake && lhs.accel == rhs.accel && lhs.pf == rhs.pf && lhs.vf == rhs.vf
            && lhs.af == rhs.af && lhs.limits == rhs.limits && lhs.direction == rhs.direction
            && lhs.controlSigns == rhs.controlSigns
    }
}
