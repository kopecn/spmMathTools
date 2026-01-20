import Foundation

// MARK: - Helper Functions

/// Helper function to format arrays for printing, mimicking ruckig/utils.hpp join
@inline(__always)
private func join<T>(_ values: [T], _ precise: Bool = false) -> String {
    if precise, let arr = values as? [Double] {
        return arr.map { (x: Double) -> String in
            if x.isFinite {
                return String(format: "%.15g", x)
            } else if x.isNaN {
                return "nan"
            } else if x == Double.infinity {
                return "inf"
            } else {
                return "-inf"
            }
        }.joined(separator: ", ")
    } else {
        return values.map { "\($0)" }.joined(separator: ", ")
    }
}

@inline(__always)
private func join<T>(_ values: [[T]], _ precise: Bool = false) -> String {
    values.map { "[" + join($0, precise) + "]" }.joined(separator: ", ")
}

// MARK: - InputParameter

/// Input parameter for trajectory generation.
///
/// This struct contains all the configuration and constraints needed to generate
/// a trajectory, including current and target kinematic states, limits on velocity,
/// acceleration, and jerk, as well as optional per-section constraints and intermediate
/// waypoints.
///
/// The input parameter supports multiple degrees of freedom (DOFs) and can be configured
/// for position or velocity control with various synchronization modes.
public struct InputParameter: Equatable {

    // MARK: - Properties

    /// Number of degrees of freedom for this input parameter
    public var degreesOfFreedom: Int

    // MARK: Control Configuration

    /// Control interface mode (Position or Velocity)
    public var controlInterface: ControlInterface = .Position

    /// Synchronization mode across degrees of freedom
    public var synchronization: Synchronization = .Time

    /// Duration discretization mode (Continuous or Discrete)
    public var durationDiscretization: DurationDiscretization = .Continuous

    // MARK: Current Kinematic State

    /// Current position for each degree of freedom
    public var currentPosition: [Double]

    /// Current velocity for each degree of freedom
    public var currentVelocity: [Double]

    /// Current acceleration for each degree of freedom
    public var currentAcceleration: [Double]

    // MARK: Target Kinematic State

    /// Target position for each degree of freedom
    public var targetPosition: [Double]

    /// Target velocity for each degree of freedom
    public var targetVelocity: [Double]

    /// Target acceleration for each degree of freedom
    public var targetAcceleration: [Double]

    // MARK: Global Limits

    /// Maximum velocity limit for each degree of freedom
    public var maxVelocity: [Double]

    /// Maximum acceleration limit for each degree of freedom
    public var maxAcceleration: [Double]

    /// Maximum jerk limit for each degree of freedom
    public var maxJerk: [Double]

    /// Optional minimum velocity limit (defaults to -maxVelocity)
    public var minVelocity: [Double]? = nil

    /// Optional minimum acceleration limit (defaults to -maxAcceleration)
    public var minAcceleration: [Double]? = nil

    /// Optional maximum position limit
    public var maxPosition: [Double]? = nil

    /// Optional minimum position limit
    public var minPosition: [Double]? = nil

    // MARK: Per-Section Constraints

    /// Intermediate waypoint positions to pass through
    public var intermediatePositions: [[Double]] = []

    /// Per-section maximum velocity limits
    public var perSectionMaxVelocity: [[Double]]? = nil

    /// Per-section maximum acceleration limits
    public var perSectionMaxAcceleration: [[Double]]? = nil

    /// Per-section maximum jerk limits
    public var perSectionMaxJerk: [[Double]]? = nil

    /// Per-section minimum velocity limits
    public var perSectionMinVelocity: [[Double]]? = nil

    /// Per-section minimum acceleration limits
    public var perSectionMinAcceleration: [[Double]]? = nil

    /// Per-section maximum position limits
    public var perSectionMaxPosition: [[Double]]? = nil

    /// Per-section minimum position limits
    public var perSectionMinPosition: [[Double]]? = nil

    // MARK: DOF-Specific Configuration

    /// Enabled/disabled flag for each degree of freedom
    public var enabled: [Bool]

    /// Optional per-DOF control interface override
    public var perDofControlInterface: [ControlInterface]? = nil

    /// Optional per-DOF synchronization override
    public var perDofSynchronization: [Synchronization]? = nil

    // MARK: Duration Constraints

    /// Optional global minimum trajectory duration
    public var minimumDuration: Double? = nil

    /// Optional per-section minimum durations
    public var perSectionMinimumDuration: [Double]? = nil

    /// Optional duration threshold for interrupting calculation
    public var interruptCalculationDuration: Double? = nil

    // MARK: - Initialization

    /// Creates an input parameter with the specified degrees of freedom.
    ///
    /// - Parameter DOFs: Number of degrees of freedom (must be positive)
    public init(DOFs: Int) {
        precondition(DOFs > 0, "Degrees of freedom must be positive")

        self.degreesOfFreedom = DOFs

        // Initialize all arrays with appropriate sizes and defaults
        self.currentPosition = Array(repeating: 0.0, count: DOFs)
        self.currentVelocity = Array(repeating: 0.0, count: DOFs)
        self.currentAcceleration = Array(repeating: 0.0, count: DOFs)
        self.targetPosition = Array(repeating: 0.0, count: DOFs)
        self.targetVelocity = Array(repeating: 0.0, count: DOFs)
        self.targetAcceleration = Array(repeating: 0.0, count: DOFs)
        self.maxVelocity = Array(repeating: 0.0, count: DOFs)
        self.maxAcceleration = Array(repeating: Double.infinity, count: DOFs)
        self.maxJerk = Array(repeating: Double.infinity, count: DOFs)
        self.enabled = Array(repeating: true, count: DOFs)
    }

    // MARK: - Validation

    /// Validates the input parameters and throws an error if invalid.
    ///
    /// - Parameters:
    ///   - checkCurrentStateWithinLimits: Whether to validate current state is within limits
    ///   - checkTargetStateWithinLimits: Whether to validate target state is within limits
    /// - Throws: `RuckigError` if validation fails
    @discardableResult
    public func validate(
        checkCurrentStateWithinLimits: Bool = false,
        checkTargetStateWithinLimits: Bool = true
    ) throws -> Bool {
        try validateJerkLimits()
        try validateAccelerationLimits(
            checkCurrent: checkCurrentStateWithinLimits,
            checkTarget: checkTargetStateWithinLimits
        )
        try validateVelocityLimits(
            checkCurrent: checkCurrentStateWithinLimits,
            checkTarget: checkTargetStateWithinLimits
        )
        try validatePositionLimits()
        try validateIntermediatePositions()

        return true
    }

    // MARK: - Private Validation Methods

    private func validateJerkLimits() throws {
        for dof in 0..<degreesOfFreedom {
            let jMax = maxJerk[dof]
            if jMax.isNaN || jMax < 0.0 {
                throw RuckigError("maximum jerk limit \(jMax) of DoF \(dof) should be larger than or equal to zero.")
            }
        }
    }

    private func validateAccelerationLimits(checkCurrent: Bool, checkTarget: Bool) throws {
        for dof in 0..<degreesOfFreedom {
            let aMax = maxAcceleration[dof]
            if aMax.isNaN || aMax < 0.0 {
                throw RuckigError("maximum acceleration limit \(aMax) of DoF \(dof) should be larger than or equal to zero.")
            }

            let aMin = minAcceleration?[dof] ?? (-aMax)
            if aMin.isNaN || aMin > 0.0 {
                throw RuckigError("minimum acceleration limit \(aMin) of DoF \(dof) should be smaller than or equal to zero.")
            }

            let a0 = currentAcceleration[dof]
            if a0.isNaN {
                throw RuckigError("current acceleration \(a0) of DoF \(dof) should be a valid number.")
            }

            let af = targetAcceleration[dof]
            if af.isNaN {
                throw RuckigError("target acceleration \(af) of DoF \(dof) should be a valid number.")
            }

            if checkCurrent {
                if a0 > aMax {
                    throw RuckigError("current acceleration \(a0) of DoF \(dof) exceeds its maximum acceleration limit \(aMax).")
                }
                if a0 < aMin {
                    throw RuckigError("current acceleration \(a0) of DoF \(dof) undercuts its minimum acceleration limit \(aMin).")
                }
            }

            if checkTarget {
                if af > aMax {
                    throw RuckigError("target acceleration \(af) of DoF \(dof) exceeds its maximum acceleration limit \(aMax).")
                }
                if af < aMin {
                    throw RuckigError("target acceleration \(af) of DoF \(dof) undercuts its minimum acceleration limit \(aMin).")
                }
            }
        }
    }

    private func validateVelocityLimits(checkCurrent: Bool, checkTarget: Bool) throws {
        for dof in 0..<degreesOfFreedom {
            let v0 = currentVelocity[dof]
            if v0.isNaN {
                throw RuckigError("current velocity \(v0) of DoF \(dof) should be a valid number.")
            }

            let vf = targetVelocity[dof]
            if vf.isNaN {
                throw RuckigError("target velocity \(vf) of DoF \(dof) should be a valid number.")
            }

            let controlInterface_ = perDofControlInterface?[dof] ?? controlInterface
            if controlInterface_ == .Position {
                let vMax = maxVelocity[dof]
                if vMax.isNaN || vMax < 0.0 {
                    throw RuckigError("maximum velocity limit \(vMax) of DoF \(dof) should be larger than or equal to zero.")
                }

                let vMin = minVelocity?[dof] ?? (-vMax)
                if vMin.isNaN || vMin > 0.0 {
                    throw RuckigError("minimum velocity limit \(vMin) of DoF \(dof) should be smaller than or equal to zero.")
                }

                if checkCurrent {
                    if v0 > vMax {
                        throw RuckigError("current velocity \(v0) of DoF \(dof) exceeds its maximum velocity limit \(vMax).")
                    }
                    if v0 < vMin {
                        throw RuckigError("current velocity \(v0) of DoF \(dof) undercuts its minimum velocity limit \(vMin).")
                    }

                    // Check inevitable velocity violations
                    let a0 = currentAcceleration[dof]
                    let jMax = maxJerk[dof]
                    if a0 > 0 && jMax > 0 {
                        let inevitable = velocityAtAccelerationZero(v0, a0, jMax)
                        if inevitable > vMax {
                            throw RuckigError("DoF \(dof) will inevitably reach a velocity \(inevitable) from the current kinematic state that will exceed its maximum velocity limit \(vMax).")
                        }
                    }
                    if a0 < 0 && jMax > 0 {
                        let inevitable = velocityAtAccelerationZero(v0, a0, -jMax)
                        if inevitable < vMin {
                            throw RuckigError("DoF \(dof) will inevitably reach a velocity \(inevitable) from the current kinematic state that will undercut its minimum velocity limit \(vMin).")
                        }
                    }
                }

                if checkTarget {
                    if vf > vMax {
                        throw RuckigError("target velocity \(vf) of DoF \(dof) exceeds its maximum velocity limit \(vMax).")
                    }
                    if vf < vMin {
                        throw RuckigError("target velocity \(vf) of DoF \(dof) undercuts its minimum velocity limit \(vMin).")
                    }

                    // Check inevitable velocity violations for target state
                    let af = targetAcceleration[dof]
                    let jMax = maxJerk[dof]
                    if af < 0 && jMax > 0 {
                        let inevitable = velocityAtAccelerationZero(vf, af, jMax)
                        if inevitable > vMax {
                            throw RuckigError("DoF \(dof) will inevitably have reached a velocity \(inevitable) from the target kinematic state that will exceed its maximum velocity limit \(vMax).")
                        }
                    }
                    if af > 0 && jMax > 0 {
                        let inevitable = velocityAtAccelerationZero(vf, af, -jMax)
                        if inevitable < vMin {
                            throw RuckigError("DoF \(dof) will inevitably have reached a velocity \(inevitable) from the target kinematic state that will undercut its minimum velocity limit \(vMin).")
                        }
                    }
                }
            }
        }
    }

    private func validatePositionLimits() throws {
        for dof in 0..<degreesOfFreedom {
            let controlInterface_ = perDofControlInterface?[dof] ?? controlInterface
            if controlInterface_ == .Position {
                let p0 = currentPosition[dof]
                if p0.isNaN {
                    throw RuckigError("current position \(p0) of DoF \(dof) should be a valid number.")
                }

                let pf = targetPosition[dof]
                if pf.isNaN {
                    throw RuckigError("target position \(pf) of DoF \(dof) should be a valid number.")
                }
            }
        }
    }

    private func validateIntermediatePositions() throws {
        if !intermediatePositions.isEmpty && controlInterface == .Position {
            if minimumDuration != nil || durationDiscretization != .Continuous {
                throw RuckigError("Intermediate position can not be used together with a global minimum or discrete duration.")
            }

            if perDofControlInterface != nil || perDofSynchronization != nil {
                throw RuckigError("Intermediate positions can only be used together with the position control interface and a global synchronization.")
            }

            for dof in 0..<degreesOfFreedom {
                let jMax = maxJerk[dof]
                if jMax.isInfinite {
                    throw RuckigError("infinite jerk limit of DoF \(dof) is currently not supported with intermediate positions.")
                }
            }
        }
    }

    /// Calculates the velocity when acceleration reaches zero given initial conditions.
    ///
    /// - Parameters:
    ///   - v0: Initial velocity
    ///   - a0: Initial acceleration
    ///   - j: Jerk value
    /// - Returns: Velocity when acceleration reaches zero
    @inline(__always)
    private func velocityAtAccelerationZero(_ v0: Double, _ a0: Double, _ j: Double) -> Double {
        v0 + (a0 * a0) / (2 * j)
    }

    // MARK: - Equatable

    public static func == (lhs: InputParameter, rhs: InputParameter) -> Bool {
        // Compare core kinematic state first (most likely to differ)
        guard lhs.currentPosition == rhs.currentPosition,
              lhs.currentVelocity == rhs.currentVelocity,
              lhs.currentAcceleration == rhs.currentAcceleration,
              lhs.targetPosition == rhs.targetPosition,
              lhs.targetVelocity == rhs.targetVelocity,
              lhs.targetAcceleration == rhs.targetAcceleration else {
            return false
        }

        // Compare limits
        guard lhs.maxVelocity == rhs.maxVelocity,
              lhs.maxAcceleration == rhs.maxAcceleration,
              lhs.maxJerk == rhs.maxJerk,
              lhs.minVelocity == rhs.minVelocity,
              lhs.minAcceleration == rhs.minAcceleration,
              lhs.maxPosition == rhs.maxPosition,
              lhs.minPosition == rhs.minPosition else {
            return false
        }

        // Compare intermediate positions
        guard lhs.intermediatePositions == rhs.intermediatePositions else {
            return false
        }

        // Compare per-section constraints
        guard lhs.perSectionMaxVelocity == rhs.perSectionMaxVelocity,
              lhs.perSectionMaxAcceleration == rhs.perSectionMaxAcceleration,
              lhs.perSectionMaxJerk == rhs.perSectionMaxJerk,
              lhs.perSectionMinVelocity == rhs.perSectionMinVelocity,
              lhs.perSectionMinAcceleration == rhs.perSectionMinAcceleration,
              lhs.perSectionMaxPosition == rhs.perSectionMaxPosition,
              lhs.perSectionMinPosition == rhs.perSectionMinPosition else {
            return false
        }

        // Compare configuration
        guard lhs.enabled == rhs.enabled,
              lhs.controlInterface == rhs.controlInterface,
              lhs.synchronization == rhs.synchronization,
              lhs.durationDiscretization == rhs.durationDiscretization,
              lhs.perDofControlInterface == rhs.perDofControlInterface,
              lhs.perDofSynchronization == rhs.perDofSynchronization else {
            return false
        }

        // Compare duration constraints
        guard lhs.minimumDuration == rhs.minimumDuration,
              lhs.perSectionMinimumDuration == rhs.perSectionMinimumDuration else {
            return false
        }

        return true
    }
}

// MARK: - CustomStringConvertible

extension InputParameter: CustomStringConvertible {
    public var description: String {
        var ss = "\n"

        // Control configuration
        if controlInterface == .Velocity {
            ss += "inp.controlInterface = ControlInterface.Velocity\n"
        }
        if synchronization == .Phase {
            ss += "inp.synchronization = Synchronization.Phase\n"
        } else if synchronization == .None {
            ss += "inp.synchronization = Synchronization.No\n"
        }
        if durationDiscretization == .Discrete {
            ss += "inp.durationDiscretization = DurationDiscretization.Discrete\n"
        }

        // Kinematic state
        ss += "inp.currentPosition = [\(join(currentPosition, true))]\n"
        ss += "inp.currentVelocity = [\(join(currentVelocity, true))]\n"
        ss += "inp.currentAcceleration = [\(join(currentAcceleration, true))]\n"
        ss += "inp.targetPosition = [\(join(targetPosition, true))]\n"
        ss += "inp.targetVelocity = [\(join(targetVelocity, true))]\n"
        ss += "inp.targetAcceleration = [\(join(targetAcceleration, true))]\n"

        // Limits
        ss += "inp.maxVelocity = [\(join(maxVelocity, true))]\n"
        ss += "inp.maxAcceleration = [\(join(maxAcceleration, true))]\n"
        ss += "inp.maxJerk = [\(join(maxJerk, true))]\n"

        if let minVelocity {
            ss += "inp.minVelocity = [\(join(minVelocity, true))]\n"
        }
        if let minAcceleration {
            ss += "inp.minAcceleration = [\(join(minAcceleration, true))]\n"
        }
        if let minimumDuration {
            ss += "inp.minimumDuration = \(minimumDuration)\n"
        }

        // Intermediate positions
        if !intermediatePositions.isEmpty {
            ss += "inp.intermediatePositions = [\n"
            for p in intermediatePositions {
                ss += "    [\(join(p, true))],\n"
            }
            ss += "]\n"
        }

        // Position limits
        if let minPosition {
            ss += "inp.minPosition = [\(join(minPosition, true))]\n"
        }
        if let maxPosition {
            ss += "inp.maxPosition = [\(join(maxPosition, true))]\n"
        }

        return ss
    }
}
