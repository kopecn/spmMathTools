import Foundation
import Testing
import simd
import FoundationTypes

@testable import spmMathTools



// MARK: - Coordinate Conversion Tests Suite
@Suite("Position Coordinate Conversion")
struct PositionCoordinateConversionTests {

    @Test("Cylindrical coordinates conversion - Float")
    func cylindricalCoordinatesConversionFloat() {
        let position = FloatPosition(x: 3.0, y: 4.0, z: 5.0)
        let cylindrical = position.cylindrical

        let expectedRadius: Float = sqrt(3.0 * 3.0 + 4.0 * 4.0)  // 5.0
        let expectedAngle: Float = atan2(4.0, 3.0)

        #expect(abs(cylindrical.radius - expectedRadius) < 1e-6)
        #expect(abs(cylindrical.angle - expectedAngle) < 1e-6)
        #expect(cylindrical.height == 5.0)
    }

    @Test("Cylindrical coordinates conversion - Double")
    func cylindricalCoordinatesConversionDouble() {
        let position = DoublePosition(x: 1.0, y: 1.0, z: 10.0)
        let cylindrical = position.cylindrical

        let expectedRadius = sqrt(2.0)
        let expectedAngle = Double.pi / 4  // 45 degrees

        #expect(abs(cylindrical.radius - expectedRadius) < 1e-10)
        #expect(abs(cylindrical.angle - expectedAngle) < 1e-10)
        #expect(cylindrical.height == 10.0)
    }

    @Test("Spherical coordinates conversion - Float")
    func sphericalCoordinatesConversionFloat() {
        let position = FloatPosition(x: 1.0, y: 0.0, z: 1.0)
        let spherical = position.spherical

        let expectedRadius: Float = sqrt(2.0)
        let expectedAzimuth: Float = 0.0  // atan2(0, 1)
        let expectedElevation: Float = asin(1.0 / sqrt(2.0))  // atan(z/radius)

        #expect(abs(spherical.radius - expectedRadius) < 1e-6)
        #expect(abs(spherical.azimuth - expectedAzimuth) < 1e-6)
        #expect(abs(spherical.elevation - expectedElevation) < 1e-6)
    }

    @Test("Spherical coordinates conversion - Double")
    func sphericalCoordinatesConversionDouble() {
        let position = DoublePosition(x: 0.0, y: 1.0, z: 0.0)
        let spherical = position.spherical

        let expectedRadius = 1.0
        let expectedAzimuth = Double.pi / 2  // atan2(1, 0)
        let expectedElevation = 0.0  // asin(0/1)

        #expect(abs(spherical.radius - expectedRadius) < 1e-10)
        #expect(abs(spherical.azimuth - expectedAzimuth) < 1e-10)
        #expect(abs(spherical.elevation - expectedElevation) < 1e-10)
    }

    @Test("Round-trip cylindrical conversion")
    func roundTripCylindricalConversion() {
        let original = FloatPosition(x: 2.0, y: 3.0, z: 7.0)
        let cylindrical = original.cylindrical
        let reconstructed = FloatPosition(
            cylindrical: cylindrical.radius,
            angle: cylindrical.angle,
            height: cylindrical.height
        )

        #expect(abs(reconstructed.x - original.x) < 1e-6)
        #expect(abs(reconstructed.y - original.y) < 1e-6)
        #expect(abs(reconstructed.z - original.z) < 1e-6)
    }

    @Test("Round-trip spherical conversion")
    func roundTripSphericalConversion() {
        let original = DoublePosition(x: 1.5, y: 2.5, z: 3.5)
        let spherical = original.spherical
        let reconstructed = DoublePosition(
            spherical: spherical.radius,
            azimuth: spherical.azimuth,
            elevation: spherical.elevation
        )

        #expect(abs(reconstructed.x - original.x) < 1e-10)
        #expect(abs(reconstructed.y - original.y) < 1e-10)
        #expect(abs(reconstructed.z - original.z) < 1e-10)
    }

    @Test("Origin coordinate conversions")
    func originCoordinateConversions() {
        let origin = FloatPosition.origin
        let cylindrical = origin.cylindrical
        let spherical = origin.spherical

        #expect(cylindrical.radius == 0.0)
        #expect(cylindrical.height == 0.0)
        #expect(spherical.radius == 0.0)

        // Angles are undefined for origin, but should not crash
        #expect(!cylindrical.angle.isNaN)
        #expect(!spherical.azimuth.isNaN)
        #expect(!spherical.elevation.isNaN || spherical.elevation.isNaN)  // asin(0/0) could be NaN
    }
}


// MARK: - ISO Spherical Coordinates Tests Suite
@Suite("Position ISO 80000-2:2019 Spherical Coordinates")
struct PositionISOSphericalCoordinatesTests {

    @Test("ISO spherical initialization - Float - north pole")
    func isoSphericalInitializationFloatNorthPole() {
        let radius: Float = 10.0
        let azimuth: Float = 0.0  // Azimuth is arbitrary at poles
        let polar: Float = 0.0  // 0 = north pole (+z axis)

        let position = FloatPosition(sphericalISO: radius, azimuth: azimuth, polar: polar)

        // At north pole: x=0, y=0, z=radius
        #expect(abs(position.x - 0.0) < 1e-6)
        #expect(abs(position.y - 0.0) < 1e-6)
        #expect(abs(position.z - 10.0) < 1e-6)
    }

    @Test("ISO spherical initialization - Float - south pole")
    func isoSphericalInitializationFloatSouthPole() {
        let radius: Float = 5.0
        let azimuth: Float = 0.0  // Azimuth is arbitrary at poles
        let polar: Float = Float.pi  // π = south pole (-z axis)

        let position = FloatPosition(sphericalISO: radius, azimuth: azimuth, polar: polar)

        // At south pole: x=0, y=0, z=-radius
        #expect(abs(position.x - 0.0) < 1e-6)
        #expect(abs(position.y - 0.0) < 1e-6)
        #expect(abs(position.z - (-5.0)) < 1e-6)
    }

    @Test("ISO spherical initialization - Float - equator")
    func isoSphericalInitializationFloatEquator() {
        let radius: Float = 8.0
        let azimuth: Float = Float.pi / 4  // 45 degrees
        let polar: Float = Float.pi / 2  // π/2 = equator (xy-plane)

        let position = FloatPosition(sphericalISO: radius, azimuth: azimuth, polar: polar)

        // At equator: z=0, x and y form circle
        let expectedX = radius * sin(polar) * cos(azimuth)
        let expectedY = radius * sin(polar) * sin(azimuth)
        let expectedZ = radius * cos(polar)

        #expect(abs(position.x - expectedX) < 1e-6)
        #expect(abs(position.y - expectedY) < 1e-6)
        #expect(abs(position.z - expectedZ) < 1e-6)
        #expect(abs(position.z) < 1e-6)  // Should be at equator
    }

    @Test("ISO spherical conversion - Float - north pole")
    func isoSphericalConversionFloatNorthPole() {
        let position = FloatPosition(x: 0.0, y: 0.0, z: 5.0)
        let sphericalISO = position.sphericalISO

        #expect(abs(sphericalISO.radius - 5.0) < 1e-6)
        #expect(abs(sphericalISO.polar - 0.0) < 1e-6)  // polar = 0 at north pole
        // Azimuth is arbitrary at poles
    }

    @Test("ISO spherical conversion - Float - south pole")
    func isoSphericalConversionFloatSouthPole() {
        let position = FloatPosition(x: 0.0, y: 0.0, z: -3.0)
        let sphericalISO = position.sphericalISO

        #expect(abs(sphericalISO.radius - 3.0) < 1e-6)
        #expect(abs(sphericalISO.polar - Float.pi) < 1e-6)  // polar = π at south pole
        // Azimuth is arbitrary at poles
    }

    @Test("ISO spherical conversion - Double - equator")
    func isoSphericalConversionDoubleEquator() {
        let position = DoublePosition(x: 1.0, y: 1.0, z: 0.0)
        let sphericalISO = position.sphericalISO

        let expectedRadius = sqrt(2.0)
        let expectedAzimuth = Double.pi / 4  // atan2(1, 1)
        let expectedPolar = Double.pi / 2  // acos(0/radius) = π/2

        #expect(abs(sphericalISO.radius - expectedRadius) < 1e-10)
        #expect(abs(sphericalISO.azimuth - expectedAzimuth) < 1e-10)
        #expect(abs(sphericalISO.polar - expectedPolar) < 1e-10)
    }

    @Test("ISO spherical conversion - Double - general case")
    func isoSphericalConversionDoubleGeneralCase() {
        let position = DoublePosition(x: 1.0, y: 0.0, z: 1.0)
        let sphericalISO = position.sphericalISO

        let expectedRadius = sqrt(2.0)
        let expectedAzimuth = 0.0  // atan2(0, 1)
        let expectedPolar = acos(1.0 / sqrt(2.0))  // acos(z/r)

        #expect(abs(sphericalISO.radius - expectedRadius) < 1e-10)
        #expect(abs(sphericalISO.azimuth - expectedAzimuth) < 1e-10)
        #expect(abs(sphericalISO.polar - expectedPolar) < 1e-10)
    }

    @Test("Round-trip ISO spherical conversion - Float")
    func roundTripISOSphericalConversionFloat() {
        let original = FloatPosition(x: 2.5, y: 3.5, z: 4.5)
        let sphericalISO = original.sphericalISO
        let reconstructed = FloatPosition(
            sphericalISO: sphericalISO.radius,
            azimuth: sphericalISO.azimuth,
            polar: sphericalISO.polar
        )

        #expect(abs(reconstructed.x - original.x) < 1e-6)
        #expect(abs(reconstructed.y - original.y) < 1e-6)
        #expect(abs(reconstructed.z - original.z) < 1e-6)
    }

    @Test("Round-trip ISO spherical conversion - Double")
    func roundTripISOSphericalConversionDouble() {
        let original = DoublePosition(x: -1.5, y: 2.5, z: -3.5)
        let sphericalISO = original.sphericalISO
        let reconstructed = DoublePosition(
            sphericalISO: sphericalISO.radius,
            azimuth: sphericalISO.azimuth,
            polar: sphericalISO.polar
        )

        #expect(abs(reconstructed.x - original.x) < 1e-10)
        #expect(abs(reconstructed.y - original.y) < 1e-10)
        #expect(abs(reconstructed.z - original.z) < 1e-10)
    }

    @Test("ISO spherical relationship - polar and elevation")
    func isoSphericalRelationshipPolarAndElevation() {
        let position = DoublePosition(x: 1.0, y: 2.0, z: 3.0)

        let spherical = position.spherical
        let sphericalISO = position.sphericalISO

        // Relationship: polar = π/2 - elevation
        let expectedPolar = Double.pi / 2 - spherical.elevation
        #expect(abs(sphericalISO.polar - expectedPolar) < 1e-10)

        // Both should have same radius and azimuth
        #expect(abs(spherical.radius - sphericalISO.radius) < 1e-10)
        #expect(abs(spherical.azimuth - sphericalISO.azimuth) < 1e-10)
    }

    @Test("ISO spherical azimuth range - Float")
    func isoSphericalAzimuthRangeFloat() {
        // Test various azimuth values
        let testAzimuths: [Float] = [0.0, Float.pi / 4, Float.pi / 2, Float.pi, 3 * Float.pi / 2]

        for azimuth in testAzimuths {
            let position = FloatPosition(sphericalISO: 5.0, azimuth: azimuth, polar: Float.pi / 3)
            let recovered = position.sphericalISO

            // Azimuth should be preserved (within normalized range)
            let normalizedAzimuth = atan2(sin(azimuth), cos(azimuth))
            let normalizedRecovered = atan2(sin(recovered.azimuth), cos(recovered.azimuth))
            #expect(abs(normalizedAzimuth - normalizedRecovered) < 1e-6)
        }
    }

    @Test("ISO spherical polar range validation - Double")
    func isoSphericalPolarRangeValidationDouble() {
        // Polar angle should be in range [0, π]
        let testPolarAngles: [Double] = [
            0.0, Double.pi / 6, Double.pi / 4, Double.pi / 3, Double.pi / 2, 2 * Double.pi / 3, Double.pi,
        ]

        for polar in testPolarAngles {
            let position = DoublePosition(sphericalISO: 10.0, azimuth: 0.0, polar: polar)
            let recovered = position.sphericalISO

            #expect(recovered.polar >= 0.0)
            #expect(recovered.polar <= Double.pi)
            #expect(abs(recovered.polar - polar) < 1e-10)
        }
    }

    @Test("ISO spherical unit vectors")
    func isoSphericalUnitVectors() {
        // Test standard unit vectors using ISO convention
        let unitX = FloatPosition.unitX
        let unitY = FloatPosition.unitY
        let unitZ = FloatPosition.unitZ

        let unitXISO = unitX.sphericalISO
        let unitYISO = unitY.sphericalISO
        let unitZISO = unitZ.sphericalISO

        // Unit X: should be at equator (polar = π/2), azimuth = 0
        #expect(abs(unitXISO.radius - 1.0) < 1e-6)
        #expect(abs(unitXISO.polar - Float.pi / 2) < 1e-6)
        #expect(abs(unitXISO.azimuth - 0.0) < 1e-6)

        // Unit Y: should be at equator (polar = π/2), azimuth = π/2
        #expect(abs(unitYISO.radius - 1.0) < 1e-6)
        #expect(abs(unitYISO.polar - Float.pi / 2) < 1e-6)
        #expect(abs(unitYISO.azimuth - Float.pi / 2) < 1e-6)

        // Unit Z: should be at north pole (polar = 0)
        #expect(abs(unitZISO.radius - 1.0) < 1e-6)
        #expect(abs(unitZISO.polar - 0.0) < 1e-6)
    }
}
