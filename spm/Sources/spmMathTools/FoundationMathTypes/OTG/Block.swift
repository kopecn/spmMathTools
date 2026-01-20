import Foundation

/// Represents a temporal block used in motion profile synchronization.
/// A Block encapsulates the fastest motion profile (`pMin`) and up to two optional time intervals (`a` and `b`)
/// during which other, slower profiles must not be selected due to constraints (i.e. these intervals are "blocked").
///
/// This class supports constructing blocked intervals from a set of valid profiles and determining
/// whether a given time `t` is blocked or what motion profile should be selected at that time.
class Block {

    // MARK: - Properties

    /// The fastest profile found among the set of candidates (based on total time).
    /// It is used as a default when no blocking applies.
    var pMin: Profile

    /// Total time duration of the fastest profile (in seconds).
    /// Includes time for braking and acceleration.
    var tMin: Double

    /// Optional first blocked time interval `a`, with an associated profile.
    /// This interval is used to exclude a range of times where synchronization is not valid.
    var a: Interval?

    /// Optional second blocked time interval `b`, also used to exclude invalid synchronization times.
    var b: Interval?

    // MARK: - Initialization

    /// Creates an empty block with default values.
    init() {
        self.pMin = Profile()
        self.tMin = 0.0
    }

    /// Creates a new block with a given minimal profile and its duration.
    ///
    /// - Parameters:
    ///   - pMin: The fastest profile among all candidates.
    ///   - tMin: Duration of the fastest profile, including any acceleration/braking.
    init(pMin: Profile, tMin: Double) {
        self.pMin = pMin
        self.tMin = tMin
    }

    // MARK: - Static Methods

    /// Removes a profile from a list of valid profiles in-place by shifting remaining elements left.
    ///
    /// - Parameters:
    ///   - validProfiles: Array of valid profiles to operate on.
    ///   - validProfileCounter: Number of profiles currently considered valid (may be < `validProfiles.count`).
    ///   - index: The index of the profile to remove.
    static func removeProfile(
        validProfiles: inout [Profile],
        validProfileCounter: inout Int,
        index: Int
    ) {
        for i in index..<(validProfileCounter - 1) {
            validProfiles[i] = validProfiles[i + 1]
        }
        validProfileCounter -= 1
    }

    // MARK: - Instance Methods

    /// Sets the given profile as the new minimal profile, and updates `tMin`.
    /// Also clears any previously defined blocked intervals.
    ///
    /// - Parameter profile: The profile to set as the new minimal profile.
    func setMinProfile(profile: Profile) {

        // CRITICAL: Copy the profile to prevent shared reference corruption
        // Profile is a class, so without copy(), pMin would share the same object
        // that could be modified later, causing discontinuities
        pMin = profile.copy()
        tMin = pMin.tSum.last! + pMin.brake.duration + pMin.accel.duration
        a = nil
        b = nil
    }

    /// Calculates the blocking intervals (`a` and optionally `b`) from a list of valid profiles.
    ///
    /// This function applies different heuristics based on the number of valid profiles:
    /// - For 1 profile: Sets it as the minimum and returns.
    /// - For 2 profiles: Checks if durations are nearly equal, or else determines one blocking interval.
    /// - For 4 profiles: Tries to detect and remove duplicate profiles caused by numerical errors.
    /// - For 3 or 5 profiles: Constructs one or two blocking intervals accordingly.
    ///
    /// - Parameters:
    ///   - block: The block instance to modify.
    ///   - validProfiles: List of candidate motion profiles.
    ///   - validProfileCounter: Number of profiles currently considered valid.
    ///   - numericalRobust: If true, uses tolerance-based filtering for near-equality.
    /// - Returns: `true` if a valid blocking configuration could be determined, `false` otherwise.
    static func calculateBlock(
        _ block: inout Block,
        _ validProfiles: inout [Profile],
        _ validProfileCounter: inout Int,
        _ numericalRobust: Bool = true
    ) -> Bool {

        // Case 1: Only one profile is valid
        if validProfileCounter == 1 {
            block.setMinProfile(profile: validProfiles[0])
            return true

            // Case 2: Two profiles - check for near-equality
        } else if validProfileCounter == 2 {
            if abs(validProfiles[0].tSum.last! - validProfiles[1].tSum.last!) < 8.0 * Double.ulpOfOne {
                block.setMinProfile(profile: validProfiles[0])
                return true
            }

            if numericalRobust {
                let idxMin = (validProfiles[0].tSum.last! < validProfiles[1].tSum.last!) ? 0 : 1
                let idxElse1 = (idxMin + 1) % 2

                block.setMinProfile(profile: validProfiles[idxMin])
                block.a = Interval(
                    profileLeft: validProfiles[idxMin],
                    profileRight: validProfiles[idxElse1]
                )
                return true
            }

            // Case 4: Heuristic to remove duplicate profiles due to numerical issues
        } else if validProfileCounter == 4 {
            if abs(validProfiles[0].tSum.last! - validProfiles[1].tSum.last!) < 32 * Double.ulpOfOne
                && validProfiles[0].direction != validProfiles[1].direction
            {
                Block.removeProfile(
                    validProfiles: &validProfiles,
                    validProfileCounter: &validProfileCounter,
                    index: 1
                )
            } else if abs(validProfiles[2].tSum.last! - validProfiles[3].tSum.last!) < 256
                * Double.ulpOfOne
                && validProfiles[2].direction != validProfiles[3].direction
            {
                Block.removeProfile(
                    validProfiles: &validProfiles,
                    validProfileCounter: &validProfileCounter,
                    index: 3
                )
            } else if abs(validProfiles[0].tSum.last! - validProfiles[3].tSum.last!) < 256
                * Double.ulpOfOne
                && validProfiles[0].direction != validProfiles[3].direction
            {
                Block.removeProfile(
                    validProfiles: &validProfiles,
                    validProfileCounter: &validProfileCounter,
                    index: 3
                )
            } else {
                return false
            }

            // Even number of profiles not covered above is considered invalid
        } else if validProfileCounter % 2 == 0 {
            return false
        }

        // Find the profile with the shortest total time
        let idxMinIt = validProfiles.prefix(validProfileCounter).enumerated().min {
            $0.element.tSum.last! < $1.element.tSum.last!
        }
        let idxMin = idxMinIt!.offset

        block.setMinProfile(profile: validProfiles[idxMin])

        // Case 3: One blocking interval from remaining two profiles
        if validProfileCounter == 3 {
            let idxElse1 = (idxMin + 1) % 3
            let idxElse2 = (idxMin + 2) % 3

            block.a = Interval(
                profileLeft: validProfiles[idxElse1],
                profileRight: validProfiles[idxElse2]
            )
            return true

            // Case 5: Two blocking intervals formed from four other profiles
        } else if validProfileCounter == 5 {
            let idxElse1 = (idxMin + 1) % 5
            let idxElse2 = (idxMin + 2) % 5
            let idxElse3 = (idxMin + 3) % 5
            let idxElse4 = (idxMin + 4) % 5

            if validProfiles[idxElse1].direction == validProfiles[idxElse2].direction {
                block.a = Interval(
                    profileLeft: validProfiles[idxElse1],
                    profileRight: validProfiles[idxElse2]
                )
                block.b = Interval(
                    profileLeft: validProfiles[idxElse3],
                    profileRight: validProfiles[idxElse4]
                )
            } else {
                block.a = Interval(
                    profileLeft: validProfiles[idxElse1],
                    profileRight: validProfiles[idxElse4]
                )
                block.b = Interval(
                    profileLeft: validProfiles[idxElse2],
                    profileRight: validProfiles[idxElse3]
                )
            }
            return true
        }

        return false
    }

    // MARK: - Querying Methods

    /// Determines whether a given time `t` is within a blocked interval or before the minimum profile time.
    ///
    /// - Parameter t: The time to check.
    /// - Returns: `true` if the time is not allowed for synchronization, `false` otherwise.
    func isBlocked(t: Double) -> Bool {
        (t < tMin) || (a?.isBlocked(t) ?? false) || (b?.isBlocked(t) ?? false)
    }

    /// Retrieves the profile that is active at a given time `t`.
    ///
    /// - Parameter t: Time at which to retrieve the corresponding profile.
    /// - Returns: The profile valid at that time, according to the block's intervals.
    func getProfile(t: Double) -> Profile {
        if let b = b, t >= b.right {
            return b.profile!
        }
        if let a = a, t >= a.right {
            return a.profile!
        }
        return pMin
    }
}

// MARK: - Printable

/// Allows the block to be printed directly using `print()` by conforming to `CustomStringConvertible`.
extension Block: CustomStringConvertible {
    var description: String {
        let intervals = [a, b].compactMap { interval in
            interval.map { "\($0.left)] [\($0.right) " }
        }.joined()
        return "[\(tMin) \(intervals)-"
    }
}
