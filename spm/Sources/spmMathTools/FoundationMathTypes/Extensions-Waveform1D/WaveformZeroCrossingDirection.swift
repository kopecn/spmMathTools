/// Direction filter for zero crossing detection
public enum WaveformZeroCrossingDirection {
    /// Detect all crossings
    case all
    /// Detect only rising crossings
    case rising
    /// Detect only falling crossings
    case falling

    func includes(_ type: WaveformZeroCrossingType) -> Bool {
        switch self {
        case .all:
            return true
        case .rising:
            return type == .rising
        case .falling:
            return type == .falling
        }
    }
}
