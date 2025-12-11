import XCTest
import FoundationTypes
@testable import spmMathTools

final class QuaternionArithmeticTests: XCTestCase {

    // MARK: - Double Quaternion Tests

    func testDoubleQuaternionMultiplication() {
        let q1 = Quaternion<Double>(x: 1, y: 0, z: 0, w: 1)
        let q2 = Quaternion<Double>(x: 0, y: 1, z: 0, w: 1)

        let result = q1 * q2

        // Verify quaternion multiplication formula
        XCTAssertEqual(result.w, q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z, accuracy: 1e-10)
        XCTAssertEqual(result.x, q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y, accuracy: 1e-10)
        XCTAssertEqual(result.y, q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x, accuracy: 1e-10)
        XCTAssertEqual(result.z, q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w, accuracy: 1e-10)
    }

    func testDoubleQuaternionMultiplicationNotCommutative() {
        let q1 = Quaternion<Double>(x: 1, y: 0, z: 0, w: 1)
        let q2 = Quaternion<Double>(x: 0, y: 1, z: 0, w: 1)

        let result1 = q1 * q2
        let result2 = q2 * q1

        // At least one component should be different
        let isDifferent = result1.x != result2.x ||
                         result1.y != result2.y ||
                         result1.z != result2.z ||
                         result1.w != result2.w

        XCTAssertTrue(isDifferent, "Quaternion multiplication should not be commutative")
    }

    func testDoubleQuaternionScalarMultiplication() {
        let q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let scalar = 2.5

        let result1 = q * scalar
        let result2 = scalar * q

        XCTAssertEqual(result1.x, q.x * scalar, accuracy: 1e-10)
        XCTAssertEqual(result1.y, q.y * scalar, accuracy: 1e-10)
        XCTAssertEqual(result1.z, q.z * scalar, accuracy: 1e-10)
        XCTAssertEqual(result1.w, q.w * scalar, accuracy: 1e-10)

        // Scalar multiplication is commutative
        XCTAssertEqual(result1.x, result2.x, accuracy: 1e-10)
        XCTAssertEqual(result1.y, result2.y, accuracy: 1e-10)
        XCTAssertEqual(result1.z, result2.z, accuracy: 1e-10)
        XCTAssertEqual(result1.w, result2.w, accuracy: 1e-10)
    }

    func testDoubleQuaternionDivision() {
        let q = Quaternion<Double>(x: 2, y: 4, z: 6, w: 8)
        let scalar = 2.0

        let result = q / scalar

        XCTAssertEqual(result.x, 1.0, accuracy: 1e-10)
        XCTAssertEqual(result.y, 2.0, accuracy: 1e-10)
        XCTAssertEqual(result.z, 3.0, accuracy: 1e-10)
        XCTAssertEqual(result.w, 4.0, accuracy: 1e-10)
    }

    func testDoubleQuaternionAddition() {
        let q1 = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let q2 = Quaternion<Double>(x: 5, y: 6, z: 7, w: 8)

        let result = q1 + q2

        XCTAssertEqual(result.x, 6.0, accuracy: 1e-10)
        XCTAssertEqual(result.y, 8.0, accuracy: 1e-10)
        XCTAssertEqual(result.z, 10.0, accuracy: 1e-10)
        XCTAssertEqual(result.w, 12.0, accuracy: 1e-10)
    }

    func testDoubleQuaternionScalarAddition() {
        let q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let scalar = 5.0

        let result1 = q + scalar
        let result2 = scalar + q

        XCTAssertEqual(result1.x, 6.0, accuracy: 1e-10)
        XCTAssertEqual(result1.y, 7.0, accuracy: 1e-10)
        XCTAssertEqual(result1.z, 8.0, accuracy: 1e-10)
        XCTAssertEqual(result1.w, 9.0, accuracy: 1e-10)

        // Addition is commutative
        XCTAssertEqual(result1.x, result2.x, accuracy: 1e-10)
        XCTAssertEqual(result1.y, result2.y, accuracy: 1e-10)
        XCTAssertEqual(result1.z, result2.z, accuracy: 1e-10)
        XCTAssertEqual(result1.w, result2.w, accuracy: 1e-10)
    }

    func testDoubleQuaternionSubtraction() {
        let q1 = Quaternion<Double>(x: 5, y: 6, z: 7, w: 8)
        let q2 = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)

        let result = q1 - q2

        XCTAssertEqual(result.x, 4.0, accuracy: 1e-10)
        XCTAssertEqual(result.y, 4.0, accuracy: 1e-10)
        XCTAssertEqual(result.z, 4.0, accuracy: 1e-10)
        XCTAssertEqual(result.w, 4.0, accuracy: 1e-10)
    }

    func testDoubleQuaternionScalarSubtraction() {
        let q = Quaternion<Double>(x: 5, y: 6, z: 7, w: 8)
        let scalar = 2.0

        let result1 = q - scalar
        let result2 = scalar - q

        XCTAssertEqual(result1.x, 3.0, accuracy: 1e-10)
        XCTAssertEqual(result1.y, 4.0, accuracy: 1e-10)
        XCTAssertEqual(result1.z, 5.0, accuracy: 1e-10)
        XCTAssertEqual(result1.w, 6.0, accuracy: 1e-10)

        XCTAssertEqual(result2.x, -3.0, accuracy: 1e-10)
        XCTAssertEqual(result2.y, -4.0, accuracy: 1e-10)
        XCTAssertEqual(result2.z, -5.0, accuracy: 1e-10)
        XCTAssertEqual(result2.w, -6.0, accuracy: 1e-10)
    }

    // MARK: - Float Quaternion Tests

    func testFloatQuaternionMultiplication() {
        let q1 = Quaternion<Float>(x: 1, y: 0, z: 0, w: 1)
        let q2 = Quaternion<Float>(x: 0, y: 1, z: 0, w: 1)

        let result = q1 * q2

        XCTAssertEqual(result.w, q1.w * q2.w - q1.x * q2.x - q1.y * q2.y - q1.z * q2.z, accuracy: 1e-5)
        XCTAssertEqual(result.x, q1.w * q2.x + q1.x * q2.w + q1.y * q2.z - q1.z * q2.y, accuracy: 1e-5)
        XCTAssertEqual(result.y, q1.w * q2.y - q1.x * q2.z + q1.y * q2.w + q1.z * q2.x, accuracy: 1e-5)
        XCTAssertEqual(result.z, q1.w * q2.z + q1.x * q2.y - q1.y * q2.x + q1.z * q2.w, accuracy: 1e-5)
    }

    func testFloatQuaternionScalarMultiplication() {
        let q = Quaternion<Float>(x: 1, y: 2, z: 3, w: 4)
        let scalar: Float = 2.5

        let result1 = q * scalar
        let result2 = scalar * q

        XCTAssertEqual(result1.x, q.x * scalar, accuracy: 1e-5)
        XCTAssertEqual(result1.y, q.y * scalar, accuracy: 1e-5)
        XCTAssertEqual(result1.z, q.z * scalar, accuracy: 1e-5)
        XCTAssertEqual(result1.w, q.w * scalar, accuracy: 1e-5)

        XCTAssertEqual(result1.x, result2.x, accuracy: 1e-5)
        XCTAssertEqual(result1.y, result2.y, accuracy: 1e-5)
        XCTAssertEqual(result1.z, result2.z, accuracy: 1e-5)
        XCTAssertEqual(result1.w, result2.w, accuracy: 1e-5)
    }

    func testFloatQuaternionDivision() {
        let q = Quaternion<Float>(x: 2, y: 4, z: 6, w: 8)
        let scalar: Float = 2.0

        let result = q / scalar

        XCTAssertEqual(result.x, 1.0, accuracy: 1e-5)
        XCTAssertEqual(result.y, 2.0, accuracy: 1e-5)
        XCTAssertEqual(result.z, 3.0, accuracy: 1e-5)
        XCTAssertEqual(result.w, 4.0, accuracy: 1e-5)
    }

    func testFloatQuaternionAddition() {
        let q1 = Quaternion<Float>(x: 1, y: 2, z: 3, w: 4)
        let q2 = Quaternion<Float>(x: 5, y: 6, z: 7, w: 8)

        let result = q1 + q2

        XCTAssertEqual(result.x, 6.0, accuracy: 1e-5)
        XCTAssertEqual(result.y, 8.0, accuracy: 1e-5)
        XCTAssertEqual(result.z, 10.0, accuracy: 1e-5)
        XCTAssertEqual(result.w, 12.0, accuracy: 1e-5)
    }

    func testFloatQuaternionSubtraction() {
        let q1 = Quaternion<Float>(x: 5, y: 6, z: 7, w: 8)
        let q2 = Quaternion<Float>(x: 1, y: 2, z: 3, w: 4)

        let result = q1 - q2

        XCTAssertEqual(result.x, 4.0, accuracy: 1e-5)
        XCTAssertEqual(result.y, 4.0, accuracy: 1e-5)
        XCTAssertEqual(result.z, 4.0, accuracy: 1e-5)
        XCTAssertEqual(result.w, 4.0, accuracy: 1e-5)
    }

    // MARK: - Identity Tests

    func testQuaternionMultiplicationWithIdentity() {
        let q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let identity = Quaternion<Double>.identity

        let result1 = q * identity
        let result2 = identity * q

        XCTAssertEqual(result1.x, q.x, accuracy: 1e-10)
        XCTAssertEqual(result1.y, q.y, accuracy: 1e-10)
        XCTAssertEqual(result1.z, q.z, accuracy: 1e-10)
        XCTAssertEqual(result1.w, q.w, accuracy: 1e-10)

        XCTAssertEqual(result2.x, q.x, accuracy: 1e-10)
        XCTAssertEqual(result2.y, q.y, accuracy: 1e-10)
        XCTAssertEqual(result2.z, q.z, accuracy: 1e-10)
        XCTAssertEqual(result2.w, q.w, accuracy: 1e-10)
    }
}
