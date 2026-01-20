//
//  PositionFirstOrderStep2.swift
//
//
//  Created by Nicholas Bergantz on 4/25/24.
//

import Foundation

/// Mathematical equations for Step 2 in first-order position interface: Time synchronization
class PositionFirstOrderStep2 {
    let tf: Double
    let _vMax: Double
    let _vMin: Double
    let pd: Double  // Pre-calculated expressions

    init(
        tf: Double,
        p0: Double,
        pf: Double,
        vMax: Double,
        vMin: Double
    ) {
        self.tf = tf
        self._vMax = vMax
        self._vMin = vMin
        self.pd = pf - p0
    }

    /// Calculates the time profile for a first-order position interface with timing constraints.
    ///
    /// - Parameter profile: An inout `Profile` object to store the calculated time profile.
    /// - Returns: `true` if the time profile can be calculated within the given constraints, `false` otherwise.
    ///
    /// This function calculates the time profile for a first-order position interface, where the position changes linearly over time. The function takes into account the maximum and minimum velocity constraints, as well as the total time `tf` for the motion. The calculated time profile is stored in the provided `Profile` object.
    func getProfile(profile: inout Profile) -> Bool {
        let vf = pd / tf

        profile.t[0] = 0
        profile.t[1] = 0
        profile.t[2] = 0
        profile.t[3] = tf
        profile.t[4] = 0
        profile.t[5] = 0
        profile.t[6] = 0

        return profile.checkForFirstOrderWithTiming(
            tf,
            vf,
            _vMax,
            _vMin,
            ControlSigns.UDDU,
            ReachedLimits.NONE
        )
    }
}
