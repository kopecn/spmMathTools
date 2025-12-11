import XCTest
import FoundationTypes
import simd
@testable import spmMathTools

final class QuaternionVectorTests: XCTestCase {

    // MARK: - Conjugate Tests

    func testDoubleQuaternionConjugate() {
        let q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let conjugate = q.conjugate

        XCTAssertEqual(conjugate.x, -1.0, accuracy: 1e-10)
        XCTAssertEqual(conjugate.y, -2.0, accuracy: 1e-10)
        XCTAssertEqual(conjugate.z, -3.0, accuracy: 1e-10)
        XCTAssertEqual(conjugate.w, 4.0, accuracy: 1e-10)
    }

    func testFloatQuaternionConjugate() {
        let q = Quaternion<Float>(x: 1, y: 2, z: 3, w: 4)
        let conjugate = q.conjugate

        XCTAssertEqual(conjugate.x, -1.0, accuracy: 1e-5)
        XCTAssertEqual(conjugate.y, -2.0, accuracy: 1e-5)
        XCTAssertEqual(conjugate.z, -3.0, accuracy: 1e-5)
        XCTAssertEqual(conjugate.w, 4.0, accuracy: 1e-5)
    }

    func testConjugateOfConjugateIsOriginal() {
        let q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let doubleConjugate = q.conjugate.conjugate

        XCTAssertEqual(doubleConjugate.x, q.x, accuracy: 1e-10)
        XCTAssertEqual(doubleConjugate.y, q.y, accuracy: 1e-10)
        XCTAssertEqual(doubleConjugate.z, q.z, accuracy: 1e-10)
        XCTAssertEqual(doubleConjugate.w, q.w, accuracy: 1e-10)
    }

    // MARK: - Inverse Tests

    func testDoubleQuaternionInverse() {
        let q = Quaternion<Double>(x: 1, y: 0, z: 0, w: 1).normalized
        let inverse = q.inverse

        // q * q.inverse should equal identity
        let result = q * inverse

        XCTAssertEqual(result.x, 0.0, accuracy: 1e-10)
        XCTAssertEqual(result.y, 0.0, accuracy: 1e-10)
        XCTAssertEqual(result.z, 0.0, accuracy: 1e-10)
        XCTAssertEqual(result.w, 1.0, accuracy: 1e-10)
    }

    func testFloatQuaternionInverse() {
        let q = Quaternion<Float>(x: 0, y: 1, z: 0, w: 1).normalized
        let inverse = q.inverse

        let result = q * inverse

        XCTAssertEqual(result.x, 0.0, accuracy: 1e-5)
        XCTAssertEqual(result.y, 0.0, accuracy: 1e-5)
        XCTAssertEqual(result.z, 0.0, accuracy: 1e-5)
        XCTAssertEqual(result.w, 1.0, accuracy: 1e-5)
    }

    func testInverseOfUnitQuaternionIsConjugate() {
        let q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4).normalized
        let inverse = q.inverse
        let conjugate = q.conjugate

        XCTAssertEqual(inverse.x, conjugate.x, accuracy: 1e-10)
        XCTAssertEqual(inverse.y, conjugate.y, accuracy: 1e-10)
        XCTAssertEqual(inverse.z, conjugate.z, accuracy: 1e-10)
        XCTAssertEqual(inverse.w, conjugate.w, accuracy: 1e-10)
    }

    // MARK: - Dot Product Tests

    func testDoubleQuaternionDotProduct() {
        let q1 = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let q2 = Quaternion<Double>(x: 5, y: 6, z: 7, w: 8)

        let dot = q1.dot(q2)
        let expected = 1.0 * 5.0 + 2.0 * 6.0 + 3.0 * 7.0 + 4.0 * 8.0

        XCTAssertEqual(dot, expected, accuracy: 1e-10)
    }

    func testFloatQuaternionDotProduct() {
        let q1 = Quaternion<Float>(x: 1, y: 2, z: 3, w: 4)
        let q2 = Quaternion<Float>(x: 5, y: 6, z: 7, w: 8)

        let dot = q1.dot(q2)
        // 1*5 + 2*6 + 3*7 + 4*8 = 5 + 12 + 21 + 32 = 70
        let expected: Float = 70.0

        XCTAssertEqual(dot, expected, accuracy: 1e-5)
    }

    func testDotProductOperator() {
        let q1 = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let q2 = Quaternion<Double>(x: 5, y: 6, z: 7, w: 8)

        let dotMethod = q1.dot(q2)
        let dotOperator = q1 • q2

        XCTAssertEqual(dotMethod, dotOperator, accuracy: 1e-10)
    }

    func testDotProductWithSelf() {
        let q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let dot = q.dot(q)
        let magnitudeSquared = q.magnitudeSquared

        XCTAssertEqual(dot, magnitudeSquared, accuracy: 1e-10)
    }

    // MARK: - Rotation Tests

    func testDoubleQuaternionRotatePosition() {
        // 90-degree rotation around Z-axis
        let q = Quaternion<Double>(axis: SIMD3<Double>(0, 0, 1), angle: .pi / 2)
        let position = Position<Double>(x: 1, y: 0, z: 0)

        let rotated = q.rotate(position)

        // After 90-degree rotation around Z, (1,0,0) -> (0,1,0)
        XCTAssertEqual(rotated.x, 0.0, accuracy: 1e-10)
        XCTAssertEqual(rotated.y, 1.0, accuracy: 1e-10)
        XCTAssertEqual(rotated.z, 0.0, accuracy: 1e-10)
    }

    func testFloatQuaternionRotatePosition() {
        // 90-degree rotation around Y-axis
        let q = Quaternion<Float>(axis: SIMD3<Float>(0, 1, 0), angle: .pi / 2)
        let position = Position<Float>(x: 1, y: 0, z: 0)

        let rotated = q.rotate(position)

        // After 90-degree rotation around Y, (1,0,0) -> (0,0,-1)
        XCTAssertEqual(rotated.x, 0.0, accuracy: 1e-5)
        XCTAssertEqual(rotated.y, 0.0, accuracy: 1e-5)
        XCTAssertEqual(rotated.z, -1.0, accuracy: 1e-5)
    }

    func testRotationOperator() {
        let q = Quaternion<Double>(axis: SIMD3<Double>(0, 0, 1), angle: .pi / 2)
        let position = Position<Double>(x: 1, y: 0, z: 0)

        let rotatedMethod = q.rotate(position)
        let rotatedOperator = q * position

        XCTAssertEqual(rotatedMethod.x, rotatedOperator.x, accuracy: 1e-10)
        XCTAssertEqual(rotatedMethod.y, rotatedOperator.y, accuracy: 1e-10)
        XCTAssertEqual(rotatedMethod.z, rotatedOperator.z, accuracy: 1e-10)
    }

    func testIdentityQuaternionDoesNotRotate() {
        let q = Quaternion<Double>.identity
        let position = Position<Double>(x: 1, y: 2, z: 3)

        let rotated = q.rotate(position)

        XCTAssertEqual(rotated.x, position.x, accuracy: 1e-10)
        XCTAssertEqual(rotated.y, position.y, accuracy: 1e-10)
        XCTAssertEqual(rotated.z, position.z, accuracy: 1e-10)
    }

    func test180DegreeRotation() {
        // 180-degree rotation around X-axis
        let q = Quaternion<Double>(axis: SIMD3<Double>(1, 0, 0), angle: .pi)
        let position = Position<Double>(x: 0, y: 1, z: 0)

        let rotated = q.rotate(position)

        // After 180-degree rotation around X, (0,1,0) -> (0,-1,0)
        XCTAssertEqual(rotated.x, 0.0, accuracy: 1e-10)
        XCTAssertEqual(rotated.y, -1.0, accuracy: 1e-10)
        XCTAssertEqual(rotated.z, 0.0, accuracy: 1e-10)
    }

    func testRotationPreservesDistance() {
        let axis = simd_normalize(SIMD3<Double>(1, 1, 1))
        let q = Quaternion<Double>(axis: axis, angle: .pi / 3)
        let position = Position<Double>(x: 1, y: 2, z: 3)

        let originalDistanceSquared = (position.x * position.x) + (position.y * position.y) + (position.z * position.z)
        let originalDistance = sqrt(originalDistanceSquared)

        let rotated = q.rotate(position)
        let rotatedDistanceSquared = (rotated.x * rotated.x) + (rotated.y * rotated.y) + (rotated.z * rotated.z)
        let rotatedDistance = sqrt(rotatedDistanceSquared)

        XCTAssertEqual(originalDistance, rotatedDistance, accuracy: 1e-10)
    }

    func testCombinedRotations() {
        // Rotate 90° around Z, then 90° around X
        let q1 = Quaternion<Double>(axis: SIMD3<Double>(0, 0, 1), angle: .pi / 2)
        let q2 = Quaternion<Double>(axis: SIMD3<Double>(1, 0, 0), angle: .pi / 2)
        let combined = q2 * q1

        let position = Position<Double>(x: 1, y: 0, z: 0)

        // Apply rotations separately
        let step1 = q1.rotate(position)
        let step2 = q2.rotate(step1)

        // Apply combined rotation
        let result = combined.rotate(position)

        XCTAssertEqual(result.x, step2.x, accuracy: 1e-10)
        XCTAssertEqual(result.y, step2.y, accuracy: 1e-10)
        XCTAssertEqual(result.z, step2.z, accuracy: 1e-10)
    }
}
