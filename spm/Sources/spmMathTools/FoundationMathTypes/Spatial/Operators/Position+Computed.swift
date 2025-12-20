import Foundation
import FoundationTypes
import simd

// MARK: - Computed Properties for Double Position

extension Position where T == Double {

    /// The magnitude (distance from origin) of the position
    @inlinable
    public var magnitude: T {
        return simd_length(vector)
    }

    /// The squared magnitude of the position
    @inlinable
    public var magnitudeSquared: T {
        return simd_length_squared(vector)
    }

    /// Get a normalized copy of the position (unit vector)
    @inlinable
    public var normalized: Position<T> {
        let mag = magnitude
        if mag > T.ulpOfOne {
            return Position<T>(vector: simd_normalize(vector))
        } else {
            return .origin
        }
    }

    /// Normalize the position in place to make it a unit vector
    @inlinable
    public mutating func normalize() {
        let mag = magnitude
        if mag > T.ulpOfOne {
            vector = simd_normalize(vector)
        } else {
            self = .origin
        }
    }

    /// Check if this is a unit position (magnitude ≈ 1)
    @inlinable
    public var isUnit: Bool {
        let mag = magnitude
        return abs(mag - 1) < T.ulpOfOne * 10
    }

    /// Convert to cylindrical coordinates (radius, angle, height)
    public var cylindrical: (radius: T, angle: T, height: T) {
        let radius = sqrt(x * x + y * y)
        let angle = simd.atan2(y, x)
        return (radius: radius, angle: angle, height: z)
    }

    /// Convert to spherical coordinates (mathematical/geographic convention)
    /// - Returns: A tuple with (radius, azimuth, elevation) where:
    ///   - radius: The radial distance from the origin
    ///   - azimuth: The azimuthal angle in radians (longitude)
    ///   - elevation: The elevation angle in radians (latitude), measured from the equator
    public var spherical: (radius: T, azimuth: T, elevation: T) {
        let radius = magnitude
        let azimuth = simd.atan2(y, x)
        let elevation = simd.asin(z / radius)
        return (radius: radius, azimuth: azimuth, elevation: elevation)
    }

    /// Convert to spherical coordinates (ISO 80000-2:2019 physics convention)
    /// - Returns: A tuple with (radius, azimuth, polar) where:
    ///   - radius: The radial distance from the origin
    ///   - azimuth: Azimuthal angle in radians (0 to 2π), measured from the positive x-axis
    ///   - polar: Polar angle (colatitude/zenith angle) in radians (0 to π), measured from +z axis
    public var sphericalISO: (radius: T, azimuth: T, polar: T) {
        let radius = magnitude
        let azimuth = simd.atan2(y, x)
        let polar = simd.acos(z / radius)
        return (radius: radius, azimuth: azimuth, polar: polar)
    }
}

// MARK: - Computed Properties for Float Position

extension Position where T == Float {

    /// The magnitude (distance from origin) of the position
    @inlinable
    public var magnitude: T {
        return simd_length(vector)
    }

    /// The squared magnitude of the position
    @inlinable
    public var magnitudeSquared: T {
        return simd_length_squared(vector)
    }

    /// Get a normalized copy of the position (unit vector)
    @inlinable
    public var normalized: Position<T> {
        let mag = magnitude
        if mag > T.ulpOfOne {
            return Position<T>(vector: simd_normalize(vector))
        } else {
            return .origin
        }
    }

    /// Normalize the position in place to make it a unit vector
    @inlinable
    public mutating func normalize() {
        let mag = magnitude
        if mag > T.ulpOfOne {
            vector = simd_normalize(vector)
        } else {
            self = .origin
        }
    }

    /// Check if this is a unit position (magnitude ≈ 1)
    @inlinable
    public var isUnit: Bool {
        let mag = magnitude
        return abs(mag - 1) < T.ulpOfOne * 10
    }

    /// Convert to cylindrical coordinates (radius, angle, height)
    public var cylindrical: (radius: T, angle: T, height: T) {
        let radius = sqrt(x * x + y * y)
        let angle = simd.atan2(y, x)
        return (radius: radius, angle: angle, height: z)
    }

    /// Convert to spherical coordinates (mathematical/geographic convention)
    /// - Returns: A tuple with (radius, azimuth, elevation) where:
    ///   - radius: The radial distance from the origin
    ///   - azimuth: The azimuthal angle in radians (longitude)
    ///   - elevation: The elevation angle in radians (latitude), measured from the equator
    public var spherical: (radius: T, azimuth: T, elevation: T) {
        let radius = magnitude
        let azimuth = simd.atan2(y, x)
        let elevation = simd.asin(z / radius)
        return (radius: radius, azimuth: azimuth, elevation: elevation)
    }

    /// Convert to spherical coordinates (ISO 80000-2:2019 physics convention)
    /// - Returns: A tuple with (radius, azimuth, polar) where:
    ///   - radius: The radial distance from the origin
    ///   - azimuth: Azimuthal angle in radians (0 to 2π), measured from the positive x-axis
    ///   - polar: Polar angle (colatitude/zenith angle) in radians (0 to π), measured from +z axis
    public var sphericalISO: (radius: T, azimuth: T, polar: T) {
        let radius = magnitude
        let azimuth = simd.atan2(y, x)
        let polar = simd.acos(z / radius)
        return (radius: radius, azimuth: azimuth, polar: polar)
    }
}
