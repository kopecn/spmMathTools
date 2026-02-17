//
//  VelocitySecondOrderStep1.swift
//
//
//  Created by Nicholas Bergantz on 4/25/24.
//

import Foundation

/// Mathematical equations for Step 1 in second-order velocity interface: Extremal profiles
class VelocitySecondOrderStep1 {
    let _aMax: Double
    let _aMin: Double

    /// Pre-calculated expressions
    let vd: Double

    init(v0: Double, vf: Double, aMax: Double, aMin: Double) {
        self._aMax = aMax
        self._aMin = aMin
        self.vd = vf - v0
    }

    func getProfile(_ input: inout Profile, _ block: inout Block) -> Bool {
        var p = block.pMin
        p.setBoundary(input)

        let af = (vd > 0) ? _aMax : _aMin
        p.t[0] = 0
        p.t[1] = vd / af
        p.t[2] = 0
        p.t[3] = 0
        p.t[4] = 0
        p.t[5] = 0
        p.t[6] = 0

        if p.checkForSecondOrderVelocity(af, .UDDU, .ACC0) {
            block.tMin = p.tSum.last! + p.brake.duration + p.accel.duration
            block.pMin = p
            return true
        }
        return false
    }
}
