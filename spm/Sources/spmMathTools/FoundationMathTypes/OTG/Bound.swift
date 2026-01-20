import Foundation

/// Represents information about position extrema (minimum and maximum positions)
/// within a motion profile trajectory.
///
/// This is typically used to determine when a system reaches its furthest
/// extents in position, which can be helpful for:
/// - Validating if position limits are exceeded
/// - Collision detection
/// - Calculating bounds for trajectory planning
///
/// Values should be derived from a computed motion `Profile`.
struct Bound {

    /// The minimum position value reached during the motion.
    ///
    /// This corresponds to the lowest point on the position trajectory curve.
    var min: Double

    /// The maximum position value reached during the motion.
    ///
    /// This corresponds to the highest point on the position trajectory curve.
    var max: Double

    /// The time (in seconds) at which the minimum position `min` is reached.
    ///
    /// Useful for time-aligned checks or dynamic limit analysis.
    var tMin: Double

    /// The time (in seconds) at which the maximum position `max` is reached.
    ///
    /// This may be used for trajectory peak detection or enforcing
    /// duration-based constraints.
    var tMax: Double

    /// Default initializer with sensible defaults
    init() {
        self.min = Double.infinity
        self.max = -Double.infinity
        self.tMin = 0.0
        self.tMax = 0.0
    }

    /// Initialize with specific values
    init(min: Double, max: Double, tMin: Double, tMax: Double) {
        self.min = min
        self.max = max
        self.tMin = tMin
        self.tMax = tMax
    }
}
