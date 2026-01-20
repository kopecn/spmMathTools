import Foundation

/// Output parameter containing the trajectory result and current kinematic state.
///
/// This struct holds the calculated trajectory and maintains the current position,
/// velocity, acceleration, and jerk values. It also tracks timing information,
/// section changes, and computational metrics.
public struct OutputParameter {
    // MARK: - Properties

    /// Number of degrees of freedom for this output parameter
    public var degreesOfFreedom: Int

    /// The calculated trajectory containing all profile information
    public var trajectory: Trajectory

    // MARK: Current kinematic state

    /// Current position values for each degree of freedom
    public var newPosition: [Double]

    /// Current velocity values for each degree of freedom
    public var newVelocity: [Double]

    /// Current acceleration values for each degree of freedom
    public var newAcceleration: [Double]

    /// Current jerk values for each degree of freedom
    public var newJerk: [Double]

    // MARK: Trajectory state

    /// Current time on the trajectory
    public var time: Double = 0.0

    /// Index of the current section between intermediate positions
    public var newSection: Int = 0

    // MARK: State flags

    /// Indicates whether the trajectory section changed during the last update
    public var didSectionChange: Bool = false

    /// Indicates whether a new trajectory was calculated
    public var newCalculation: Bool = false

    /// Computational duration of the last update call [µs]
    public var calculationDuration: Double = 0.0

    // MARK: - Initialization

    /// Creates an output parameter with the specified degrees of freedom.
    ///
    /// - Parameter DOFs: Number of degrees of freedom (must be positive)
    public init(DOFs: Int) {
        precondition(DOFs > 0, "Degrees of freedom must be positive")

        self.degreesOfFreedom = DOFs
        self.trajectory = Trajectory(dofs: DOFs)
        self.newPosition = Array(repeating: 0.0, count: DOFs)
        self.newVelocity = Array(repeating: 0.0, count: DOFs)
        self.newAcceleration = Array(repeating: 0.0, count: DOFs)
        self.newJerk = Array(repeating: 0.0, count: DOFs)
    }

    /// Creates an output parameter with the specified degrees of freedom.
    ///
    /// - Parameters:
    ///   - dofs: Number of degrees of freedom (must be positive)
    ///   - maxNumberOfWaypoints: Optional maximum number of intermediate waypoints.
    ///                           If provided, pre-allocates trajectory storage for waypoints.
    public init(dofs: Int, maxNumberOfWaypoints: Int? = nil) {
        precondition(dofs > 0, "Degrees of freedom must be positive")

        self.degreesOfFreedom = dofs

        // Initialize trajectory with waypoint support if specified
        if let maxWaypoints = maxNumberOfWaypoints {
            self.trajectory = Trajectory(dofs: dofs, maxNumberOfWaypoints: maxWaypoints)
        } else {
            self.trajectory = Trajectory(dofs: dofs)
        }

        // Initialize kinematic state arrays
        self.newPosition = Array(repeating: 0.0, count: dofs)
        self.newVelocity = Array(repeating: 0.0, count: dofs)
        self.newAcceleration = Array(repeating: 0.0, count: dofs)
        self.newJerk = Array(repeating: 0.0, count: dofs)
    }

    // MARK: - Public Methods

    /// Passes the current kinematic state to an input parameter for trajectory continuation.
    ///
    /// This method updates the input parameter's current state with this output's state,
    /// enabling trajectory chaining. If a section change occurred and intermediate positions
    /// exist, the first intermediate position is removed.
    ///
    /// - Parameter input: The input parameter to update
    /// - Throws: `preconditionFailure` if array sizes don't match degrees of freedom
    public func passToInput(_ input: inout InputParameter) {
        // Validate array sizes match
        precondition(
            newPosition.count == degreesOfFreedom &&
            newVelocity.count == degreesOfFreedom &&
            newAcceleration.count == degreesOfFreedom,
            "Output parameter arrays must match degreesOfFreedom"
        )

        input.currentPosition = self.newPosition
        input.currentVelocity = self.newVelocity
        input.currentAcceleration = self.newAcceleration

        if didSectionChange && !input.intermediatePositions.isEmpty {
            input.intermediatePositions.removeFirst()
        }
    }
}

// MARK: - CustomStringConvertible

extension OutputParameter: CustomStringConvertible {
    public var description: String {
        func fmt(_ values: [Double]) -> String {
            values.map { String(format: "%.6f", $0) }.joined(separator: ", ")
        }

        var out = "\nout.newPosition = [\(fmt(newPosition))]\n"
        out += "out.newVelocity = [\(fmt(newVelocity))]\n"
        out += "out.newAcceleration = [\(fmt(newAcceleration))]\n"
        out += "out.newJerk = [\(fmt(newJerk))]\n"
        out += "out.time = [\(String(format: "%.16f", time))]\n"
        out += "out.calculationDuration = [\(String(format: "%.16f", calculationDuration))]\n"
        return out
    }
}
