import XCTest
import FoundationTypes
import simd
@testable import spmMathTools

final class SpatialPoseDenavitHartenbergTests: XCTestCase {

    // MARK: - Double DH Tests

    func testDoubleDHInitializerMatchesMatrix() {
        // Test parameters
        let a = 0.5
        let alpha = Double.pi / 4
        let d = 1.0
        let theta = Double.pi / 3

        // Create using DH initializer
        let poseFromDH = SpatialPose<Double>(
            denavitHartenberg: a,
            alpha: alpha,
            d: d,
            theta: theta
        )

        // Create using matrix approach
        let dhMatrix = simd_double4x4(
            denavitHartenberg: a,
            alpha: alpha,
            d: d,
            theta: theta
        )
        let poseFromMatrix = SpatialPose<Double>(homogeneousTransform: dhMatrix)

        // Verify position matches
        XCTAssertEqual(poseFromDH.x, poseFromMatrix.x, accuracy: 1e-10)
        XCTAssertEqual(poseFromDH.y, poseFromMatrix.y, accuracy: 1e-10)
        XCTAssertEqual(poseFromDH.z, poseFromMatrix.z, accuracy: 1e-10)

        // Verify quaternion matches (or represents same rotation)
        // Quaternions q and -q represent the same rotation
        let dotProduct = abs(
            poseFromDH.qx * poseFromMatrix.qx +
            poseFromDH.qy * poseFromMatrix.qy +
            poseFromDH.qz * poseFromMatrix.qz +
            poseFromDH.qw * poseFromMatrix.qw
        )

        XCTAssertEqual(dotProduct, 1.0, accuracy: 1e-10, "Quaternions should represent the same rotation")
    }

    func testDoubleDHWithPrecomputedAlpha() {
        let a = 0.5
        let alpha = Double.pi / 4
        let ca = cos(alpha)
        let sa = sin(alpha)
        let d = 1.0
        let theta = Double.pi / 3

        // Standard DH
        let pose1 = SpatialPose<Double>(
            denavitHartenberg: a,
            alpha: alpha,
            d: d,
            theta: theta
        )

        // Precomputed alpha DH
        let pose2 = SpatialPose<Double>(
            denavitHartenberg: a,
            ca: ca,
            sa: sa,
            d: d,
            theta: theta
        )

        // Should produce identical results
        XCTAssertEqual(pose1.x, pose2.x, accuracy: 1e-10)
        XCTAssertEqual(pose1.y, pose2.y, accuracy: 1e-10)
        XCTAssertEqual(pose1.z, pose2.z, accuracy: 1e-10)
        XCTAssertEqual(pose1.qx, pose2.qx, accuracy: 1e-10)
        XCTAssertEqual(pose1.qy, pose2.qy, accuracy: 1e-10)
        XCTAssertEqual(pose1.qz, pose2.qz, accuracy: 1e-10)
        XCTAssertEqual(pose1.qw, pose2.qw, accuracy: 1e-10)
    }

    func testDoubleDHZeroTheta() {
        let a = 1.0
        let alpha = Double.pi / 2
        let d = 0.5
        let theta = 0.0

        let pose = SpatialPose<Double>(
            denavitHartenberg: a,
            alpha: alpha,
            d: d,
            theta: theta
        )

        // When theta = 0: position should be (a, 0, d)
        XCTAssertEqual(pose.x, a, accuracy: 1e-10)
        XCTAssertEqual(pose.y, 0.0, accuracy: 1e-10)
        XCTAssertEqual(pose.z, d, accuracy: 1e-10)
    }

    func testDoubleDHIdentityParameters() {
        // All zero parameters should give identity-like pose
        let pose = SpatialPose<Double>(
            denavitHartenberg: 0,
            alpha: 0,
            d: 0,
            theta: 0
        )

        XCTAssertEqual(pose.x, 0.0, accuracy: 1e-10)
        XCTAssertEqual(pose.y, 0.0, accuracy: 1e-10)
        XCTAssertEqual(pose.z, 0.0, accuracy: 1e-10)
        // Rotation should be identity quaternion
        XCTAssertEqual(pose.qw, 1.0, accuracy: 1e-10)
    }

    // MARK: - Float DH Tests

    func testFloatDHInitializerMatchesMatrix() {
        let a: Float = 0.5
        let alpha: Float = .pi / 4
        let d: Float = 1.0
        let theta: Float = .pi / 3

        let poseFromDH = SpatialPose<Float>(
            denavitHartenberg: a,
            alpha: alpha,
            d: d,
            theta: theta
        )

        let dhMatrix = simd_float4x4(
            denavitHartenberg: a,
            alpha: alpha,
            d: d,
            theta: theta
        )
        let poseFromMatrix = SpatialPose<Float>(homogeneousTransform: dhMatrix)

        XCTAssertEqual(poseFromDH.x, poseFromMatrix.x, accuracy: 1e-5)
        XCTAssertEqual(poseFromDH.y, poseFromMatrix.y, accuracy: 1e-5)
        XCTAssertEqual(poseFromDH.z, poseFromMatrix.z, accuracy: 1e-5)

        let dotProduct = abs(
            poseFromDH.qx * poseFromMatrix.qx +
            poseFromDH.qy * poseFromMatrix.qy +
            poseFromDH.qz * poseFromMatrix.qz +
            poseFromDH.qw * poseFromMatrix.qw
        )

        XCTAssertEqual(dotProduct, 1.0, accuracy: 1e-5)
    }

    func testFloatDHWithPrecomputedAlpha() {
        let a: Float = 0.5
        let alpha: Float = .pi / 4
        let ca = cos(alpha)
        let sa = sin(alpha)
        let d: Float = 1.0
        let theta: Float = .pi / 3

        let pose1 = SpatialPose<Float>(
            denavitHartenberg: a,
            alpha: alpha,
            d: d,
            theta: theta
        )

        let pose2 = SpatialPose<Float>(
            denavitHartenberg: a,
            ca: ca,
            sa: sa,
            d: d,
            theta: theta
        )

        XCTAssertEqual(pose1.x, pose2.x, accuracy: 1e-5)
        XCTAssertEqual(pose1.y, pose2.y, accuracy: 1e-5)
        XCTAssertEqual(pose1.z, pose2.z, accuracy: 1e-5)
        XCTAssertEqual(pose1.qx, pose2.qx, accuracy: 1e-5)
        XCTAssertEqual(pose1.qy, pose2.qy, accuracy: 1e-5)
        XCTAssertEqual(pose1.qz, pose2.qz, accuracy: 1e-5)
        XCTAssertEqual(pose1.qw, pose2.qw, accuracy: 1e-5)
    }

    // MARK: - Chain Tests

    func testDHChainComposition() {
        // Create a simple 2-joint chain
        let joint1 = SpatialPose<Double>(
            denavitHartenberg: 1.0,
            alpha: Double.pi / 2,
            d: 0.0,
            theta: Double.pi / 4
        )

        let joint2 = SpatialPose<Double>(
            denavitHartenberg: 0.5,
            alpha: 0.0,
            d: 0.0,
            theta: Double.pi / 6
        )

        // Compose using pose multiplication
        let composed = joint1 * joint2

        // Verify composition produces valid normalized quaternion
        let qMag = sqrt(
            composed.qx * composed.qx +
            composed.qy * composed.qy +
            composed.qz * composed.qz +
            composed.qw * composed.qw
        )

        XCTAssertEqual(qMag, 1.0, accuracy: 1e-10, "Composed quaternion should be normalized")
    }

    func testDHVariousThetaValues() {
        let testThetas: [Double] = [0, .pi / 6, .pi / 4, .pi / 3, .pi / 2, .pi, 3 * .pi / 2, 2 * .pi]
        let a = 1.0
        let alpha = Double.pi / 4
        let d = 0.5

        for theta in testThetas {
            let pose = SpatialPose<Double>(
                denavitHartenberg: a,
                alpha: alpha,
                d: d,
                theta: theta
            )

            // Position should match: (a*cos(θ), a*sin(θ), d)
            XCTAssertEqual(pose.x, a * cos(theta), accuracy: 1e-10)
            XCTAssertEqual(pose.y, a * sin(theta), accuracy: 1e-10)
            XCTAssertEqual(pose.z, d, accuracy: 1e-10)

            // Quaternion should be normalized
            let qMag = sqrt(pose.qx * pose.qx + pose.qy * pose.qy + pose.qz * pose.qz + pose.qw * pose.qw)
            XCTAssertEqual(qMag, 1.0, accuracy: 1e-10)
        }
    }
}
