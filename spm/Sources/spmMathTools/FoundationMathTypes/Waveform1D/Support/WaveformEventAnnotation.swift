/// Event annotation types
public enum WaveformEventAnnotation {
    /// Simple marker
    case marker
    /// Peak marker
    case peak
    /// Valley marker
    case valley
    /// Custom annotation with label
    case custom(String)
}
