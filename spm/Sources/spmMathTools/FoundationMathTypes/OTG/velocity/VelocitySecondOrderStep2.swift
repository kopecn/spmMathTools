//
//  VelocitySecondOrderStep2.swift
//
//
//  Created by Nicholas Bergantz on 4/25/24.
//

import Foundation

/// Mathematical equations for Step 2 in second-order velocity interface: Time synchronization
class VelocitySecondOrderStep2 {
    let tf: Double
    let _aMax: Double
    let _aMin: Double

    /// Pre-calculated expressions
    let vd: Double

    init(tf: Double, v0: Double, vf: Double, aMax: Double, aMin: Double) {
        self.tf = tf
        self._aMax = aMax
        self._aMin = aMin
        self.vd = vf - v0
    }

    func getProfile(_ profile: inout Profile) -> Bool {
        let af = vd / tf

        let debugThis = true  // (abs(vd) < 0.01)
        if debugThis {
            print("[VelocitySecondOrderStep2] vd=\(vd), tf=\(tf), af=\(af)")
            print("  aMin=\(_aMin), aMax=\(_aMax)")
        }

        profile.t[0] = 0
        profile.t[1] = tf
        profile.t[2] = 0
        profile.t[3] = 0
        profile.t[4] = 0
        profile.t[5] = 0
        profile.t[6] = 0

        let checkResult = profile.checkForSecondOrderVelocityWithTiming(tf, af, _aMax, _aMin, .UDDU, .NONE)
        if debugThis {
            print("  checkForSecondOrderVelocityWithTiming returned: \(checkResult)")
        }

        if checkResult {
            profile.pf = profile.p.last!
            return true
        }

        return false
    }
}
