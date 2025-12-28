import Foundation
import FoundationTypes

// MARK: - Mutating Operations
extension Waveform1D {

    /// Append a single value to the end of the waveform
    /// - Parameter value: The value to append
    /// - Note: This adds one sample to the waveform
    public mutating func append(_ value: T) {
        self.values.append(value)
    }

    /// Append multiple values to the end of the waveform
    /// - Parameter values: Array of values to append
    public mutating func append(contentsOf values: [T]) {
        self.values.append(contentsOf: values)
    }

    /// Prepend a single value to the beginning of the waveform
    /// - Parameter value: The value to prepend
    /// - Note: This adds one sample to the beginning and may shift t0 if it's set
    public mutating func prepend(_ value: T) {
        self.values.insert(value, at: 0)

        // Adjust t0 if it exists (shift back by dt)
        if let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(add: -dt)
        }
    }

    /// Prepend multiple values to the beginning of the waveform
    /// - Parameter values: Array of values to prepend
    /// - Note: The values are prepended in order, so values[0] becomes the first sample
    public mutating func prepend(contentsOf values: [T]) {
        self.values.insert(contentsOf: values, at: 0)

        // Adjust t0 if it exists (shift back by dt * count)
        if let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(add: -dt * U(values.count))
        }
    }

    /// Insert a single value at the specified index
    /// - Parameters:
    ///   - value: The value to insert
    ///   - index: The index at which to insert the value
    /// - Note: If inserting at index 0, t0 is adjusted. Otherwise, this creates a temporal discontinuity
    ///         in the uniformly sampled data.
    public mutating func insert(_ value: T, at index: Int) {
        precondition(index >= 0 && index <= values.count, "Index out of bounds")

        self.values.insert(value, at: index)

        // Adjust t0 if inserting at the beginning
        if index == 0, let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(add: -dt)
        }
    }

    /// Insert multiple values at the specified index
    /// - Parameters:
    ///   - values: Array of values to insert
    ///   - index: The index at which to insert the values
    /// - Note: If inserting at index 0, t0 is adjusted. Otherwise, this creates a temporal discontinuity
    ///         in the uniformly sampled data.
    public mutating func insert(contentsOf values: [T], at index: Int) {
        precondition(index >= 0 && index <= self.values.count, "Index out of bounds")

        self.values.insert(contentsOf: values, at: index)

        // Adjust t0 if inserting at the beginning
        if index == 0, let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(add: -dt * U(values.count))
        }
    }

    /// Replace the value at the specified index
    /// - Parameters:
    ///   - index: The index of the value to replace
    ///   - value: The new value
    /// - Precondition: index must be within bounds [0, sampleCount)
    public mutating func replace(at index: Int, with value: T) {
        precondition(index >= 0 && index < values.count, "Index out of bounds")
        self.values[index] = value
    }

    /// Replace a range of values
    /// - Parameters:
    ///   - range: The range of indices to replace
    ///   - values: The new values
    public mutating func replaceSubrange<C>(_ range: Range<Int>, with values: C)
    where C: Collection, C.Element == T {
        self.values.replaceSubrange(range, with: values)
    }

    /// Remove the value at the specified index
    /// - Parameter index: The index of the value to remove
    /// - Returns: The removed value
    /// - Note: If removing at index 0, t0 is adjusted forward by dt
    @discardableResult
    public mutating func remove(at index: Int) -> T {
        precondition(index >= 0 && index < values.count, "Index out of bounds")

        let removedValue = self.values.remove(at: index)

        // Adjust t0 if removing the first element
        if index == 0, let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(add: dt)
        }

        return removedValue
    }

    /// Remove all values
    public mutating func removeAll() {
        self.values.removeAll()
    }

    /// Remove all values and optionally keep capacity
    /// - Parameter keepingCapacity: If true, keeps the underlying storage capacity
    public mutating func removeAll(keepingCapacity: Bool) {
        self.values.removeAll(keepingCapacity: keepingCapacity)
    }

    /// Remove the first value
    /// - Returns: The removed value, or nil if the waveform is empty
    @discardableResult
    public mutating func removeFirst() -> T? {
        guard !values.isEmpty else { return nil }
        return remove(at: 0)
    }

    /// Remove the last value
    /// - Returns: The removed value, or nil if the waveform is empty
    @discardableResult
    public mutating func removeLast() -> T? {
        guard !values.isEmpty else { return nil }
        return values.removeLast()
    }
}

// MARK: - Subscript Access
extension Waveform1D {

    /// Access a value at the specified index
    /// - Parameter index: The index of the value to access
    public subscript(index: Int) -> T {
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
