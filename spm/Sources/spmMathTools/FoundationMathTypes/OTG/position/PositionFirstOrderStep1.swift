//
//  PositionFirstOrderStep1.swift
//
//
//  Created by Nicholas Bergantz on 4/25/24.
//

import Foundation

/// Mathematical equations for Step 1 in first-order position interface: Extremal profiles
class PositionFirstOrderStep1 {

    let _vMax: Double
    let _vMin: Double
    let pd: Double
    /// Pre-calculated expressions

    init(
        p0: Double,
        pf: Double,
        vMax: Double,
        vMin: Double
    ) {
        _vMax = vMax
        _vMin = vMin
        pd = pf - p0
    }

    /// Calculates the profile for a given block, setting the boundary conditions and
    /// determining the time parameters for the first-order motion.
    ///
    /// - Parameters:
    ///   - input: The input `Profile` object to be updated.
    ///   - block: The `Block` object containing the motion parameters.
    /// - Returns: `true` if the first-order motion profile is valid, `false` otherwise.
    func getProfile(
        input: inout Profile,
        block: inout Block
    ) -> Bool {

        var p = block.pMin
        p.setBoundary(input)

        let vf = (pd > 0) ? _vMax : _vMin
        p.t[0] = 0
        p.t[1] = 0
        p.t[2] = 0
        p.t[3] = pd / vf
        p.t[4] = 0
        p.t[5] = 0
        p.t[6] = 0

        if p.checkForFirstOrder(vf, ControlSigns.UDDU, ReachedLimits.VEL) {
            assert(!p.tSum.isEmpty, "tSum should not be empty when checkForFirstOrder returns true")
            block.pMin = p
            block.tMin = p.tSum.last! + p.brake.duration + p.accel.duration
            return true
        }
        return false
    }
}
