public enum Synchronization: String, Codable {

    ///< Always synchronize the DoFs to reach the target at the same time (Default)
    case Time = "Time"

    ///< Synchronize only when necessary (e.g. for non-zero target velocity or acceleration)
    case TimeIfNecessary = "TimeIfNecessary"

    ///< Phase synchronize the DoFs when this is possible, else fallback to "Time" strategy. Phase synchronization will result a straight-line trajectory
    case Phase = "Phase"

    // ///< Always phase synchronize the DoFs (even when this is not time-optimal), else returns "ErrorNoPhaseSynchronization". Ruckig will then guarantee a straight-line trajectory
    // case PhaseOnly = "PhaseOnly",

    ///< Calculate every DoF independently
    case None = "None"
}
