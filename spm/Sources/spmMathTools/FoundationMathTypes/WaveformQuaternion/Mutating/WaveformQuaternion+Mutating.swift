import Foundation
import simd
import FoundationTypes

// MARK: - Mutating Operations
extension WaveformQuaternion {
    
    /// Extend  another quaternion waveform to this one
    /// Both waveforms must have the same sampling rate
    public mutating func extend(_ other: WaveformQuaternion<T>) throws {
        // Compare dt with relative tolerance
        let dtEqual: Bool
        if self.dt == 0 && other.dt == 0 {
            dtEqual = true
        } else if self.dt == 0 || other.dt == 0 {
            dtEqual = abs(self.dt - other.dt) < 1e-10
        } else {
            let relativeDifference = abs(self.dt - other.dt) / max(abs(self.dt), abs(other.dt))
            dtEqual = relativeDifference < 1e-10
        }

        guard dtEqual else {
            throw WaveformError.incompatibleSamplingRates
        }

        self.values.append(contentsOf: other.values)
    }

    /// Create a new waveform by concatenating this one with another
    public func concatenated(with other: WaveformQuaternion<T>) throws -> WaveformQuaternion<T> {
        var result = self
        try result.append(other)
        return result
    }

    /// Append a single Quaternion to the end of the waveform
    /// - Parameter quaternion: The Quaternion to append
    /// - Note: This adds one sample to the waveform
    public mutating func append(_ quaternion: Quaternion<T>) {
        self.values.append(quaternion)
    }

    /// Append multiple Quaternion values to the end of the waveform
    /// - Parameter quaternions: Array of Quaternion values to append
    public mutating func append(contentsOf quaternions: [Quaternion<T>]) {
        self.values.append(contentsOf: quaternions)
    }

    /// Prepend a single Quaternion to the beginning of the waveform
    /// - Parameter quaternion: The Quaternion to prepend
    /// - Note: This adds one sample to the beginning and may shift t0 if it's set
    public mutating func prepend(_ quaternion: Quaternion<T>) {
        self.values.insert(quaternion, at: 0)

        // Adjust t0 if it exists (shift back by dt)
        if let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(-dt)
        }
    }

    /// Prepend multiple Quaternion values to the beginning of the waveform
    /// - Parameter quaternions: Array of Quaternion values to prepend
    /// - Note: The quaternions are prepended in order, so quaternions[0] becomes the first sample
    public mutating func prepend(contentsOf quaternions: [Quaternion<T>]) {
        self.values.insert(contentsOf: quaternions, at: 0)

        // Adjust t0 if it exists (shift back by dt * count)
        if let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(-dt * TimeInterval(quaternions.count))
        }
    }

    /// Insert a single Quaternion at the specified index
    /// - Parameters:
    ///   - quaternion: The Quaternion to insert
    ///   - index: The index at which to insert the quaternion
    /// - Note: If inserting at index 0, t0 is adjusted. Otherwise, this creates a temporal discontinuity
    ///         in the uniformly sampled data.
    public mutating func insert(_ quaternion: Quaternion<T>, at index: Int) {
        precondition(index >= 0 && index <= values.count, "Index out of bounds")

        self.values.insert(quaternion, at: index)

        // Adjust t0 if inserting at the beginning
        if index == 0, let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(-dt)
        }
    }

    /// Insert multiple Quaternion values at the specified index
    /// - Parameters:
    ///   - quaternions: Array of Quaternion values to insert
    ///   - index: The index at which to insert the quaternions
    /// - Note: If inserting at index 0, t0 is adjusted. Otherwise, this creates a temporal discontinuity
    ///         in the uniformly sampled data.
    public mutating func insert(contentsOf quaternions: [Quaternion<T>], at index: Int) {
        precondition(index >= 0 && index <= values.count, "Index out of bounds")

        self.values.insert(contentsOf: quaternions, at: index)

        // Adjust t0 if inserting at the beginning
        if index == 0, let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(-dt * TimeInterval(quaternions.count))
        }
    }

    /// Replace the Quaternion at the specified index
    /// - Parameters:
    ///   - index: The index of the quaternion to replace
    ///   - quaternion: The new Quaternion value
    /// - Precondition: index must be within bounds [0, sampleCount)
    public mutating func replace(at index: Int, with quaternion: Quaternion<T>) {
        precondition(index >= 0 && index < values.count, "Index out of bounds")
        self.values[index] = quaternion
    }

    /// Replace a range of Quaternion values
    /// - Parameters:
    ///   - range: The range of indices to replace
    ///   - quaternions: The new Quaternion values
    public mutating func replaceSubrange<C>(_ range: Range<Int>, with quaternions: C)
    where C: Collection, C.Element == Quaternion<T> {
        self.values.replaceSubrange(range, with: quaternions)
    }

    /// Remove the Quaternion at the specified index
    /// - Parameter index: The index of the quaternion to remove
    /// - Returns: The removed Quaternion
    /// - Note: If removing at index 0, t0 is adjusted forward by dt
    @discardableResult
    public mutating func remove(at index: Int) -> Quaternion<T> {
        precondition(index >= 0 && index < values.count, "Index out of bounds")

        let removedQuaternion = self.values.remove(at: index)

        // Adjust t0 if removing the first element
        if index == 0, let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(dt)
        }

        return removedQuaternion
    }

    /// Remove all quaternions
    public mutating func removeAll() {
        self.values.removeAll()
    }

    /// Remove all quaternions and optionally keep capacity
    /// - Parameter keepingCapacity: If true, keeps the underlying storage capacity
    public mutating func removeAll(keepingCapacity: Bool) {
        self.values.removeAll(keepingCapacity: keepingCapacity)
    }

    /// Remove the first quaternion
    /// - Returns: The removed quaternion, or nil if the waveform is empty
    @discardableResult
    public mutating func removeFirst() -> Quaternion<T>? {
        guard !values.isEmpty else { return nil }
        return remove(at: 0)
    }

    /// Remove the last quaternion
    /// - Returns: The removed quaternion, or nil if the waveform is empty
    @discardableResult
    public mutating func removeLast() -> Quaternion<T>? {
        guard !values.isEmpty else { return nil }
        return values.removeLast()
    }
}

// MARK: - Subscript Access
extension WaveformQuaternion {

    /// Access a Quaternion at the specified index
    /// - Parameter index: The index of the quaternion to access
    public subscript(index: Int) -> Quaternion<T> {
        get {
            precondition(index >= 0 && index < values.count, "Index out of bounds")
            return values[index]
        }
        set {
            precondition(index >= 0 && index < values.count, "Index out of bounds")
            values[index] = newValue
        }
    }
}
