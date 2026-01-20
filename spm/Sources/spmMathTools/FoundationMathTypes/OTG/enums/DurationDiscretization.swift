public enum DurationDiscretization: String, Codable {
    /// Every trajectory synchronization duration is allowed (Default)
    case Continuous = "Continuous"
    /// The trajectory synchronization duration must be a multiple of the control cycle
    case Discrete = "Discrete"
}
