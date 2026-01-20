public enum ControlInterface: String, Codable {
    /// Position-control: Full control over the entire kinematic state (Default)
    case Position = "Position"
    /// Velocity-control: Ignores the current position, target position, and velocity limits
    case Velocity = "Velocity"
}
