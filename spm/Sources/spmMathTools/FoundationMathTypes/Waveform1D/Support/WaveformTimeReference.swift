/// Time reference options for time-based value access
public enum WaveformTimeReference {
    /// Time is relative to the start of the waveform (t=0 is first sample)
    case waveformStart
    /// Time is Unix epoch seconds (requires t0 to be set)
    case epoch
}
