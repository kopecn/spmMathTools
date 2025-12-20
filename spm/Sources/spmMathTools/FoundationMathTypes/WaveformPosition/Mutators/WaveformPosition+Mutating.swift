import Foundation
import FoundationTypes
import simd

// MARK: - Mutating Operations
extension WaveformPosition {

    /// Extend another position waveform to this one
    /// Both waveforms must have the same sampling rate
    public mutating func extend(_ other: WaveformPosition<T>) throws {
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
    public func concatenated(with other: WaveformPosition<T>) throws -> WaveformPosition<T> {
        var result = self
        try result.extend(other)
        return result
    }

    /// Append a single Position to the end of the waveform
    /// - Parameter position: The Position to append
    /// - Note: This adds one sample to the waveform
    public mutating func append(_ position: Position<T>) {
        self.values.append(position)
    }

    /// Append multiple Position values to the end of the waveform
    /// - Parameter positions: Array of Position values to append
    public mutating func append(contentsOf positions: [Position<T>]) {
        self.values.append(contentsOf: positions)
    }

    /// Prepend a single Position to the beginning of the waveform
    /// - Parameter position: The Position to prepend
    /// - Note: This adds one sample to the beginning and may shift t0 if it's set
    public mutating func prepend(_ position: Position<T>) {
        self.values.insert(position, at: 0)

        // Adjust t0 if it exists (shift back by dt)
        if let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(-dt)
        }
    }

    /// Prepend multiple Position values to the beginning of the waveform
    /// - Parameter positions: Array of Position values to prepend
    /// - Note: The positions are prepended in order, so positions[0] becomes the first sample
    public mutating func prepend(contentsOf positions: [Position<T>]) {
        self.values.insert(contentsOf: positions, at: 0)

        // Adjust t0 if it exists (shift back by dt * count)
        if let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(-dt * T(positions.count))
        }
    }

    /// Insert a single Position at the specified index
    /// - Parameters:
    ///   - position: The Position to insert
    ///   - index: The index at which to insert the position
    /// - Note: If inserting at index 0, t0 is adjusted. Otherwise, this creates a temporal discontinuity
    ///         in the uniformly sampled data.
    public mutating func insert(_ position: Position<T>, at index: Int) {
        precondition(index >= 0 && index <= values.count, "Index out of bounds")

        self.values.insert(position, at: index)

        // Adjust t0 if inserting at the beginning
        if index == 0, let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(-dt)
        }
    }

    /// Insert multiple Position values at the specified index
    /// - Parameters:
    ///   - positions: Array of Position values to insert
    ///   - index: The index at which to insert the positions
    /// - Note: If inserting at index 0, t0 is adjusted. Otherwise, this creates a temporal discontinuity
    ///         in the uniformly sampled data.
    public mutating func insert(contentsOf positions: [Position<T>], at index: Int) {
        precondition(index >= 0 && index <= values.count, "Index out of bounds")

        self.values.insert(contentsOf: positions, at: index)

        // Adjust t0 if inserting at the beginning
        if index == 0, let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(-dt * T(positions.count))
        }
    }

    /// Replace the Position at the specified index
    /// - Parameters:
    ///   - index: The index of the position to replace
    ///   - position: The new Position value
    /// - Precondition: index must be within bounds [0, sampleCount)
    public mutating func replace(at index: Int, with position: Position<T>) {
        precondition(index >= 0 && index < values.count, "Index out of bounds")
        self.values[index] = position
    }

    /// Replace a range of Position values
    /// - Parameters:
    ///   - range: The range of indices to replace
    ///   - positions: The new Position values
    public mutating func replaceSubrange<C>(_ range: Range<Int>, with positions: C)
    where C: Collection, C.Element == Position<T> {
        self.values.replaceSubrange(range, with: positions)
    }

    /// Remove the Position at the specified index
    /// - Parameter index: The index of the position to remove
    /// - Returns: The removed Position
    /// - Note: If removing at index 0, t0 is adjusted forward by dt
    @discardableResult
    public mutating func remove(at index: Int) -> Position<T> {
        precondition(index >= 0 && index < values.count, "Index out of bounds")

        let removedPosition = self.values.remove(at: index)

        // Adjust t0 if removing the first element
        if index == 0, let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(dt)
        }

        return removedPosition
    }

    /// Remove all positions
    public mutating func removeAll() {
        self.values.removeAll()
    }

    /// Remove all positions and optionally keep capacity
    /// - Parameter keepingCapacity: If true, keeps the underlying storage capacity
    public mutating func removeAll(keepingCapacity: Bool) {
        self.values.removeAll(keepingCapacity: keepingCapacity)
    }

    /// Remove the first position
    /// - Returns: The removed position, or nil if the waveform is empty
    @discardableResult
    public mutating func removeFirst() -> Position<T>? {
        guard !values.isEmpty else { return nil }
        return remove(at: 0)
    }

    /// Remove the last position
    /// - Returns: The removed position, or nil if the waveform is empty
    @discardableResult
    public mutating func removeLast() -> Position<T>? {
        guard !values.isEmpty else { return nil }
        return values.removeLast()
    }
}

// MARK: - Subscript Access
extension WaveformPosition {

    /// Access a Position at the specified index
    /// - Parameter index: The index of the position to access
    public subscript(index: Int) -> Position<T> {
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
