import Foundation

private let eps: Double = 2.2e-14

/// `BrakeProfile` models a pre-trajectory braking motion to bring a system
/// within acceptable kinematic limits (position, velocity, acceleration)
/// *before* initiating the main synchronized trajectory.
///
/// The braking motion is computed based on initial conditions and system limits,
/// using either second-order (acceleration-limited) or third-order (jerk-limited) profiles.
///
/// Common use cases:
/// - Avoiding kinematic constraint violations at trajectory start
/// - Stabilizing initial overshoots in trajectory blending
/// - Ensuring valid input to a trajectory generator
struct BrakeProfile {
    private(set) var duration: Double = 0.0
    public var t: [Double] = [0, 0]  // Durations of braking segments
    var j: [Double] = [0, 0]  // Jerk values for braking segments
    var a: [Double] = [0, 0]  // Acceleration values for second-order braking
    var v: [Double] = [0, 0]  // Velocity at segment start
    var p: [Double] = [0, 0]  // Position at segment start

    /// Finalizes the full (third-order) braking profile.
    ///
    /// Applies braking segments to the input kinematic state, updating
    /// the internal trajectory arrays with intermediate values.
    ///
    /// - Parameters:
    ///   - pS: Initial position (in-out)
    ///   - vS: Initial velocity (in-out)
    ///   - aS: Initial acceleration (in-out)
    mutating func finalize(
        _ pS: inout Double,
        _ vS: inout Double,
        _ aS: inout Double
    ) {
        if t[0] <= 0.0 && t[1] <= 0.0 {
            duration = 0.0
            return
        }

        duration = t[0]
        p[0] = pS
        v[0] = vS
        a[0] = aS
        (pS, vS, aS) = integrate(t[0], pS, vS, aS, j[0])

        if t[1] > 0.0 {
            duration += t[1]
            p[1] = pS
            v[1] = vS
            a[1] = aS
            (pS, vS, aS) = integrate(t[1], pS, vS, aS, j[1])
        }
    }

    /// Finalizes a second-order (acceleration-only) braking profile.
    ///
    /// This integrates the state with constant acceleration (no jerk).
    ///
    /// - Parameters:
    ///   - pS: Initial position (in-out)
    ///   - vS: Initial velocity (in-out)
    ///   - aS: Initial acceleration (in-out)
    mutating func finalizeSecondOrder(
        _ pS: inout Double,
        _ vS: inout Double,
        _ aS: inout Double
    ) {
        if t[0] <= 0.0 {
            duration = 0.0
            return
        }

        duration = t[0]
        p[0] = pS
        v[0] = vS
        (pS, vS, aS) = integrate(t[0], pS, vS, a[0], 0.0)
    }

    /// Computes velocity at time `t`, given initial velocity `v0`, acceleration `a0`, and jerk `j`.
    func vAtT(_ v0: Double, _ a0: Double, _ j: Double, _ t: Double) -> Double {
        v0 + t * (a0 + j * t / 2)
    }

    /// Computes velocity at the moment when acceleration becomes zero,
    /// based on initial velocity `v0`, acceleration `a0`, and jerk `j`.
    func vAtAZero(_ v0: Double, _ a0: Double, _ j: Double) -> Double {
        v0 + (a0 * a0) / (2 * j)
    }

    /// Builds an acceleration-based braking profile using jerk to bring the
    /// system into velocity bounds by first adjusting acceleration.
    ///
    /// - Parameters:
    ///   - v0: Initial velocity
    ///   - a0: Initial acceleration
    ///   - vMax/vMin: Velocity bounds
    ///   - aMax/aMin: Acceleration bounds
    ///   - jMax: Maximum jerk (always positive; direction handled internally)
    mutating func accelerationBrake(
        _ v0: Double,
        _ a0: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) {
        j[0] = -jMax

        let tToAMax = (a0 - aMax) / jMax
        let tToAZero = a0 / jMax

        let vAtAMax = vAtT(v0, a0, -jMax, tToAMax)
        let vAtAZero = vAtT(v0, a0, -jMax, tToAZero)

        if (vAtAZero > vMax && jMax > 0) || (vAtAZero < vMax && jMax < 0) {
            velocityBrake(v0, a0, vMax, vMin, aMax, aMin, jMax)

        } else if (vAtAMax < vMin && jMax > 0) || (vAtAMax > vMin && jMax < 0) {
            let tToVMin = -(vAtAMax - vMin) / aMax
            let tToVMax = -aMax / (2 * jMax) - (vAtAMax - vMax) / aMax

            t[0] = tToAMax + eps
            t[1] = max(min(tToVMin, tToVMax - eps), 0.0)

        } else {
            t[0] = tToAMax + eps
        }
    }

    /// Builds a braking profile focused on reducing velocity directly.
    ///
    /// Handles the case where velocity must be brought into bounds without
    /// exceeding acceleration or jerk limits.
    ///
    /// - Parameters:
    ///   - v0, a0: Initial velocity and acceleration
    ///   - vMax/vMin: Velocity bounds
    ///   - aMax/aMin: Acceleration bounds
    ///   - jMax: Maximum jerk (always positive)
    mutating func velocityBrake(
        _ v0: Double,
        _ a0: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) {
        j[0] = -jMax

        let tToAMin = (a0 - aMin) / jMax
        let tToVMax = a0 / jMax + sqrt(a0 * a0 + 2 * jMax * (v0 - vMax)) / abs(jMax)
        let tToVMin = a0 / jMax + sqrt(a0 * a0 / 2 + jMax * (v0 - vMin)) / abs(jMax)
        let tMinToV = min(tToVMax, tToVMin)

        if tToAMin < tMinToV {
            let vAtAMin = vAtT(v0, a0, -jMax, tToAMin)
            let tToVMaxWithConstant = -(vAtAMin - vMax) / aMin
            let tToVMinWithConstant = aMin / (2 * jMax) - (vAtAMin - vMin) / aMin

            t[0] = max(tToAMin - eps, 0.0)
            t[1] = max(min(tToVMaxWithConstant, tToVMinWithConstant), 0.0)

        } else {
            t[0] = max(tMinToV - eps, 0.0)
        }
    }

    /// Computes the full braking trajectory from position interface.
    ///
    /// Selects acceleration-based or velocity-based profiles depending on whether
    /// acceleration or velocity are out of bounds.
    ///
    /// - Parameters:
    ///   - v0, a0: Initial velocity and acceleration
    ///   - vMax/vMin: Velocity limits
    ///   - aMax/aMin: Acceleration limits
    ///   - jMax: Maximum jerk
    mutating func getPositionBrakeTrajectory(
        _ v0: Double,
        _ a0: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) {
        t = [0.0, 0.0]
        j = [0.0, 0.0]

        guard jMax != 0.0 && aMax != 0.0 && aMin != 0.0 else { return }

        if a0 > aMax {
            accelerationBrake(v0, a0, vMax, vMin, aMax, aMin, jMax)
        } else if a0 < aMin {
            accelerationBrake(v0, a0, vMin, vMax, aMin, aMax, -jMax)
        } else if (v0 > vMax && vAtAZero(v0, a0, -jMax) > vMin)
            || (a0 > 0 && vAtAZero(v0, a0, jMax) > vMax)
        {
            velocityBrake(v0, a0, vMax, vMin, aMax, aMin, jMax)
        } else if (v0 < vMin && vAtAZero(v0, a0, jMax) < vMax)
            || (a0 < 0 && vAtAZero(v0, a0, -jMax) < vMin)
        {
            velocityBrake(v0, a0, vMin, vMax, aMin, aMax, -jMax)
        }
    }

    /// Computes second-order (acceleration-only) braking from a position interface.
    ///
    /// - Brakes to `vMax` if initial velocity exceeds it.
    /// - Brakes to `vMin` if initial velocity is below it.
    mutating func getSecondOrderPositionBrakeTrajectory(
        _ v0: Double,
        _ vMax: Double,
        _ vMin: Double,
        _ aMax: Double,
        _ aMin: Double
    ) {
        t = [0.0, 0.0]
        j = [0.0, 0.0]
        a = [0.0, 0.0]

        guard aMax != 0.0 && aMin != 0.0 else { return }

        if v0 > vMax {
            a[0] = aMin
            t[0] = (vMax - v0) / aMin + eps
        } else if v0 < vMin {
            a[0] = aMax
            t[0] = (vMin - v0) / aMax + eps
        }
    }

    /// Computes braking trajectory for second-order velocity interface.
    ///
    /// Brings acceleration into valid range using jerk-limited braking.
    mutating func getVelocityBrakeTrajectory(
        _ a0: Double,
        _ aMax: Double,
        _ aMin: Double,
        _ jMax: Double
    ) {
        t = [0.0, 0.0]
        j = [0.0, 0.0]

        guard jMax != 0.0 else { return }

        if a0 > aMax {
            j[0] = -jMax
            t[0] = (a0 - aMax) / jMax + eps
        } else if a0 < aMin {
            j[0] = jMax
            t[0] = -(a0 - aMin) / jMax + eps
        }
    }

    /// Resets the braking profile for second-order velocity interface.
    ///
    /// This is a no-op placeholder for future implementation or interface consistency.
    mutating func getSecondOrderVelocityBrakeTrajectory() {
        t = [0.0, 0.0]
        j = [0.0, 0.0]
    }
}

extension BrakeProfile: Equatable {
    static func == (lhs: BrakeProfile, rhs: BrakeProfile) -> Bool {
        lhs.duration == rhs.duration && lhs.t == rhs.t && lhs.j == rhs.j && lhs.a == rhs.a && lhs.v == rhs.v
            && lhs.p == rhs.p
    }
}
