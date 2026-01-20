//
//  Interval.swift
//
//
//  Created by Nicholas Bergantz on 5/4/24.
//

import Foundation

struct Interval {

    /// seconds
    var left: Double
    /// seconds
    var right: Double
    /// Profile corresponding to right (end) time
    var profile: Profile?

    init(
        _ left: Double,
        _ right: Double
    ) {
        self.left = left
        self.right = right
    }

    init(
        profileLeft: Profile,
        profileRight: Profile
    ) {
        let leftDuration =
            profileLeft.tSum.last! + profileLeft.brake.duration + profileLeft.accel.duration
        let rightDuration =
            profileRight.tSum.last! + profileRight.brake.duration + profileRight.accel.duration
        if leftDuration < rightDuration {
            left = leftDuration
            right = rightDuration
            profile = profileRight
        } else {
            left = rightDuration
            right = leftDuration
            profile = profileLeft
        }
    }

    func isBlocked(_ t: Double) -> Bool {
        t > left && t < right
    }
}
