import Foundation
import FoundationTypes

// MARK: - Mutating Operations
extension WaveformSpatialPose {

    /// Extend another pose waveform to this one
    /// Both waveforms must have the same sampling rate
    public mutating func extend(_ other: WaveformSpatialPose<T>) throws {
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

        self.positions.append(contentsOf: other.positions)
        self.quaternions.append(contentsOf: other.quaternions)
    }

    /// Create a new waveform by concatenating this one with another
    public func concatenated(with other: WaveformSpatialPose<T>) throws -> WaveformSpatialPose<T> {
        var result = self
        try result.extend(other)
        return result
    }

    /// Append a single SpatialPose to the end of the waveform
    /// - Parameter pose: The SpatialPose to append
    /// - Note: This adds one sample to the waveform
    public mutating func append(_ pose: SpatialPose<T>) {
        self.positions.append(pose.position)
        self.quaternions.append(pose.quaternion)
    }

    /// Append multiple SpatialPose instances to the end of the waveform
    /// - Parameter poses: Array of SpatialPose instances to append
    public mutating func append(contentsOf poses: [SpatialPose<T>]) {
        let newPositions = poses.map { $0.position }
        let newQuaternions = poses.map { $0.quaternion }
        self.positions.append(contentsOf: newPositions)
        self.quaternions.append(contentsOf: newQuaternions)
    }

    /// Prepend a single SpatialPose to the beginning of the waveform
    /// - Parameter pose: The SpatialPose to prepend
    /// - Note: This adds one sample to the beginning and may shift t0 if it's set
    public mutating func prepend(_ pose: SpatialPose<T>) {
        self.positions.insert(pose.position, at: 0)
        self.quaternions.insert(pose.quaternion, at: 0)

        // Adjust t0 if it exists (shift back by dt)
        if let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(add: -dt)
        }
    }

    /// Prepend multiple SpatialPose instances to the beginning of the waveform
    /// - Parameter poses: Array of SpatialPose instances to prepend
    /// - Note: The poses are prepended in order, so poses[0] becomes the first sample
    public mutating func prepend(contentsOf poses: [SpatialPose<T>]) {
        let newPositions = poses.map { $0.position }
        let newQuaternions = poses.map { $0.quaternion }
        self.positions.insert(contentsOf: newPositions, at: 0)
        self.quaternions.insert(contentsOf: newQuaternions, at: 0)

        // Adjust t0 if it exists (shift back by dt * count)
        if let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(add: -dt * T(poses.count))
        }
    }

    /// Insert a single SpatialPose at the specified index
    /// - Parameters:
    ///   - pose: The SpatialPose to insert
    ///   - index: The index at which to insert the pose
    /// - Note: If inserting at index 0, t0 is adjusted. Otherwise, this creates a temporal discontinuity
    ///         in the uniformly sampled data.
    public mutating func insert(_ pose: SpatialPose<T>, at index: Int) {
        precondition(index >= 0 && index <= positions.count, "Index out of bounds")
        precondition(index >= 0 && index <= quaternions.count, "Index out of bounds")

        self.positions.insert(pose.position, at: index)
        self.quaternions.insert(pose.quaternion, at: index)

        // Adjust t0 if inserting at the beginning
        if index == 0, let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(add: -dt)
        }
    }

    /// Insert multiple SpatialPose instances at the specified index
    /// - Parameters:
    ///   - poses: Array of SpatialPose instances to insert
    ///   - index: The index at which to insert the poses
    /// - Note: If inserting at index 0, t0 is adjusted. Otherwise, this creates a temporal discontinuity
    ///         in the uniformly sampled data.
    public mutating func insert(contentsOf poses: [SpatialPose<T>], at index: Int) {
        precondition(index >= 0 && index <= positions.count, "Index out of bounds")
        precondition(index >= 0 && index <= quaternions.count, "Index out of bounds")

        let newPositions = poses.map { $0.position }
        let newQuaternions = poses.map { $0.quaternion }
        self.positions.insert(contentsOf: newPositions, at: index)
        self.quaternions.insert(contentsOf: newQuaternions, at: index)

        // Adjust t0 if inserting at the beginning
        if index == 0, let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(add: -dt * T(poses.count))
        }
    }

    /// Replace the SpatialPose at the specified index
    /// - Parameters:
    ///   - index: The index of the pose to replace
    ///   - pose: The new SpatialPose value
    /// - Precondition: index must be within bounds [0, sampleCount)
    public mutating func replace(at index: Int, with pose: SpatialPose<T>) {
        precondition(index >= 0 && index < positions.count, "Index out of bounds for positions")
        precondition(index >= 0 && index < quaternions.count, "Index out of bounds for quaternions")

        self.positions[index] = pose.position
        self.quaternions[index] = pose.quaternion
    }

    /// Replace a range of SpatialPose values
    /// - Parameters:
    ///   - range: The range of indices to replace
    ///   - poses: The new SpatialPose values
    /// - Note: The range count must match the poses array count
    public mutating func replaceSubrange<C>(_ range: Range<Int>, with poses: C)
    where C: Collection, C.Element == SpatialPose<T> {
        let newPositions = poses.map { $0.position }
        let newQuaternions = poses.map { $0.quaternion }
        self.positions.replaceSubrange(range, with: newPositions)
        self.quaternions.replaceSubrange(range, with: newQuaternions)
    }

    /// Remove the SpatialPose at the specified index
    /// - Parameter index: The index of the pose to remove
    /// - Returns: The removed SpatialPose
    /// - Note: If removing at index 0, t0 is adjusted forward by dt
    @discardableResult
    public mutating func remove(at index: Int) -> (position: Position<T>, quaternion: Quaternion<T>) {
        precondition(index >= 0 && index < positions.count, "Index out of bounds for positions")
        precondition(index >= 0 && index < quaternions.count, "Index out of bounds for quaternions")

        let removedPosition = self.positions.remove(at: index)
        let removedQuaternion = self.quaternions.remove(at: index)

        // Adjust t0 if removing the first element
        if index == 0, let currentT0 = self.t0 {
            self.t0 = currentT0.addingTimeInterval(add: dt)
        }

        return (position: removedPosition, quaternion: removedQuaternion)
    }

    /// Remove all poses
    public mutating func removeAll() {
        self.positions.removeAll()
        self.quaternions.removeAll()
    }

    /// Remove all poses and optionally keep capacity
    /// - Parameter keepingCapacity: If true, keeps the underlying storage capacity
    public mutating func removeAll(keepingCapacity: Bool) {
        self.positions.removeAll(keepingCapacity: keepingCapacity)
        self.quaternions.removeAll(keepingCapacity: keepingCapacity)
    }
}

// MARK: - Subscript Access
extension WaveformSpatialPose {

    /// Access a SpatialPose at the specified index
    /// - Parameter index: The index of the pose to access
    /// - Returns: A SpatialPose constructed from the position and quaternion at the index
    /// - Note: This creates a new SpatialPose on each access. For read-only access to components,
    ///         consider accessing positions[index] or quaternions[index] directly.
    public subscript(index: Int) -> SpatialPose<T> {
        get {
            precondition(index >= 0 && index < sampleCount, "Index out of bounds")
            return SpatialPose<T>(
                position: positions[index],
                rotation: quaternions[index]
            )
        }
        set {
            precondition(index >= 0 && index < positions.count, "Index out of bounds for positions")
            precondition(index >= 0 && index < quaternions.count, "Index out of bounds for quaternions")
            positions[index] = newValue.position
            quaternions[index] = newValue.quaternion
        }
    }
}
