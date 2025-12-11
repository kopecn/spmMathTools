import XCTest
import FoundationTypes
@testable import spmMathTools

final class QuaternionUnaryTests: XCTestCase {

    // MARK: - Double Quaternion Unary Tests

    func testDoubleQuaternionNegation() {
        let q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let negated = -q

        XCTAssertEqual(negated.x, -1.0, accuracy: 1e-10)
        XCTAssertEqual(negated.y, -2.0, accuracy: 1e-10)
        XCTAssertEqual(negated.z, -3.0, accuracy: 1e-10)
        XCTAssertEqual(negated.w, -4.0, accuracy: 1e-10)
    }

    func testDoubleQuaternionUnaryPlus() {
        let q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let plus = +q

        XCTAssertEqual(plus.x, q.x, accuracy: 1e-10)
        XCTAssertEqual(plus.y, q.y, accuracy: 1e-10)
        XCTAssertEqual(plus.z, q.z, accuracy: 1e-10)
        XCTAssertEqual(plus.w, q.w, accuracy: 1e-10)
    }

    func testDoubleNegationIsInvolutive() {
        let q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let doubleNegated = -(-q)

        XCTAssertEqual(doubleNegated.x, q.x, accuracy: 1e-10)
        XCTAssertEqual(doubleNegated.y, q.y, accuracy: 1e-10)
        XCTAssertEqual(doubleNegated.z, q.z, accuracy: 1e-10)
        XCTAssertEqual(doubleNegated.w, q.w, accuracy: 1e-10)
    }

    // MARK: - Float Quaternion Unary Tests

    func testFloatQuaternionNegation() {
        let q = Quaternion<Float>(x: 1, y: 2, z: 3, w: 4)
        let negated = -q

        XCTAssertEqual(negated.x, -1.0, accuracy: 1e-5)
        XCTAssertEqual(negated.y, -2.0, accuracy: 1e-5)
        XCTAssertEqual(negated.z, -3.0, accuracy: 1e-5)
        XCTAssertEqual(negated.w, -4.0, accuracy: 1e-5)
    }

    func testFloatQuaternionUnaryPlus() {
        let q = Quaternion<Float>(x: 1, y: 2, z: 3, w: 4)
        let plus = +q

        XCTAssertEqual(plus.x, q.x, accuracy: 1e-5)
        XCTAssertEqual(plus.y, q.y, accuracy: 1e-5)
        XCTAssertEqual(plus.z, q.z, accuracy: 1e-5)
        XCTAssertEqual(plus.w, q.w, accuracy: 1e-5)
    }

    // MARK: - Negation Properties

    func testNegatedQuaternionAdditionCancels() {
        let q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let result = q + (-q)

        XCTAssertEqual(result.x, 0.0, accuracy: 1e-10)
        XCTAssertEqual(result.y, 0.0, accuracy: 1e-10)
        XCTAssertEqual(result.z, 0.0, accuracy: 1e-10)
        XCTAssertEqual(result.w, 0.0, accuracy: 1e-10)
    }

    func testNegatedQuaternionRepresentsOppositeRotation() {
        // For unit quaternions, -q represents the same rotation as q
        // (quaternions q and -q represent the same rotation)
        let q = Quaternion<Double>(x: 1, y: 0, z: 0, w: 1).normalized
        let negated = -q

        let position = Position<Double>(x: 1, y: 2, z: 3)
        let rotated1 = q.rotate(position)
        let rotated2 = negated.rotate(position)

        XCTAssertEqual(rotated1.x, rotated2.x, accuracy: 1e-10)
        XCTAssertEqual(rotated1.y, rotated2.y, accuracy: 1e-10)
        XCTAssertEqual(rotated1.z, rotated2.z, accuracy: 1e-10)
    }

    func testNegationWithScalarMultiplication() {
        let q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let negated1 = -q
        let negated2 = q * -1.0

        XCTAssertEqual(negated1.x, negated2.x, accuracy: 1e-10)
        XCTAssertEqual(negated1.y, negated2.y, accuracy: 1e-10)
        XCTAssertEqual(negated1.z, negated2.z, accuracy: 1e-10)
        XCTAssertEqual(negated1.w, negated2.w, accuracy: 1e-10)
    }

    func testNegationDistributesOverAddition() {
        let q1 = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let q2 = Quaternion<Double>(x: 5, y: 6, z: 7, w: 8)

        let result1 = -(q1 + q2)
        let result2 = (-q1) + (-q2)

        XCTAssertEqual(result1.x, result2.x, accuracy: 1e-10)
        XCTAssertEqual(result1.y, result2.y, accuracy: 1e-10)
        XCTAssertEqual(result1.z, result2.z, accuracy: 1e-10)
        XCTAssertEqual(result1.w, result2.w, accuracy: 1e-10)
    }

    func testIdentityNegation() {
        let identity = Quaternion<Double>.identity
        let negated = -identity

        XCTAssertEqual(negated.x, 0.0, accuracy: 1e-10)
        XCTAssertEqual(negated.y, 0.0, accuracy: 1e-10)
        XCTAssertEqual(negated.z, 0.0, accuracy: 1e-10)
        XCTAssertEqual(negated.w, -1.0, accuracy: 1e-10)
    }

    func testZeroQuaternionNegation() {
        let zero = Quaternion<Double>.zero
        let negated = -zero

        XCTAssertEqual(negated.x, 0.0, accuracy: 1e-10)
        XCTAssertEqual(negated.y, 0.0, accuracy: 1e-10)
        XCTAssertEqual(negated.z, 0.0, accuracy: 1e-10)
        XCTAssertEqual(negated.w, 0.0, accuracy: 1e-10)
    }
}
