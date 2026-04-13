import Foundation
import FoundationTypes
import simd

// MARK: - Computed Properties for Double Position

extension Position where T == Double {

    /// Convert to cylindrical coordinates (radius, angle, height)
    @inlinable
    public var cylindrical: (radius: T, angle: T, height: T) {
        let radius = simd_length(SIMD2<Double>(vector.x, vector.y))
        let angle = simd.atan2(y, x)
        return (radius: radius, angle: angle, height: z)
    }

    /// Convert to spherical coordinates (mathematical/geographic convention)
    /// - Returns: A tuple with (radius, azimuth, elevation) where:
    ///   - radius: The radial distance from the origin
    ///   - azimuth: The azimuthal angle in radians (longitude)
    ///   - elevation: The elevation angle in radians (latitude), measured from the equator
    @inlinable
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
    @inlinable
    public var sphericalISO: (radius: T, azimuth: T, polar: T) {
        let radius = magnitude
        let azimuth = simd.atan2(y, x)
        let polar = simd.acos(z / radius)
        return (radius: radius, azimuth: azimuth, polar: polar)
    }
}

// MARK: - Computed Properties for Float Position

extension Position where T == Float {

    /// Convert to cylindrical coordinates (radius, angle, height)
    @inlinable
    public var cylindrical: (radius: T, angle: T, height: T) {
        let radius = simd_length(SIMD2<Float>(vector.x, vector.y))
        let angle = simd.atan2(y, x)
        return (radius: radius, angle: angle, height: z)
    }

    /// Convert to spherical coordinates (mathematical/geographic convention)
    /// - Returns: A tuple with (radius, azimuth, elevation) where:
    ///   - radius: The radial distance from the origin
    ///   - azimuth: The azimuthal angle in radians (longitude)
    ///   - elevation: The elevation angle in radians (latitude), measured from the equator
    @inlinable
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
    @inlinable
    public var sphericalISO: (radius: T, azimuth: T, polar: T) {
        let radius = magnitude
        let azimuth = simd.atan2(y, x)
        let polar = simd.acos(z / radius)
        return (radius: radius, azimuth: azimuth, polar: polar)
    }
}
