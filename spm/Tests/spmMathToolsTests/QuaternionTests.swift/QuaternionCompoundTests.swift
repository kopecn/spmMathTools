import XCTest
import FoundationTypes
@testable import spmMathTools

final class QuaternionCompoundTests: XCTestCase {

    // MARK: - Double Quaternion Compound Tests

    func testDoubleQuaternionMultiplyAssign() {
        var q1 = Quaternion<Double>(x: 1, y: 0, z: 0, w: 1)
        let q2 = Quaternion<Double>(x: 0, y: 1, z: 0, w: 1)
        let expected = q1 * q2

        q1 *= q2

        XCTAssertEqual(q1.x, expected.x, accuracy: 1e-10)
        XCTAssertEqual(q1.y, expected.y, accuracy: 1e-10)
        XCTAssertEqual(q1.z, expected.z, accuracy: 1e-10)
        XCTAssertEqual(q1.w, expected.w, accuracy: 1e-10)
    }

    func testDoubleQuaternionScalarMultiplyAssign() {
        var q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let scalar = 2.5
        let expected = q * scalar

        q *= scalar

        XCTAssertEqual(q.x, expected.x, accuracy: 1e-10)
        XCTAssertEqual(q.y, expected.y, accuracy: 1e-10)
        XCTAssertEqual(q.z, expected.z, accuracy: 1e-10)
        XCTAssertEqual(q.w, expected.w, accuracy: 1e-10)
    }

    func testDoubleQuaternionDivideAssign() {
        var q = Quaternion<Double>(x: 2, y: 4, z: 6, w: 8)
        let scalar = 2.0

        q /= scalar

        XCTAssertEqual(q.x, 1.0, accuracy: 1e-10)
        XCTAssertEqual(q.y, 2.0, accuracy: 1e-10)
        XCTAssertEqual(q.z, 3.0, accuracy: 1e-10)
        XCTAssertEqual(q.w, 4.0, accuracy: 1e-10)
    }

    func testDoubleQuaternionAddAssign() {
        var q1 = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let q2 = Quaternion<Double>(x: 5, y: 6, z: 7, w: 8)

        q1 += q2

        XCTAssertEqual(q1.x, 6.0, accuracy: 1e-10)
        XCTAssertEqual(q1.y, 8.0, accuracy: 1e-10)
        XCTAssertEqual(q1.z, 10.0, accuracy: 1e-10)
        XCTAssertEqual(q1.w, 12.0, accuracy: 1e-10)
    }

    func testDoubleQuaternionScalarAddAssign() {
        var q = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)
        let scalar = 5.0

        q += scalar

        XCTAssertEqual(q.x, 6.0, accuracy: 1e-10)
        XCTAssertEqual(q.y, 7.0, accuracy: 1e-10)
        XCTAssertEqual(q.z, 8.0, accuracy: 1e-10)
        XCTAssertEqual(q.w, 9.0, accuracy: 1e-10)
    }

    func testDoubleQuaternionSubtractAssign() {
        var q1 = Quaternion<Double>(x: 5, y: 6, z: 7, w: 8)
        let q2 = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4)

        q1 -= q2

        XCTAssertEqual(q1.x, 4.0, accuracy: 1e-10)
        XCTAssertEqual(q1.y, 4.0, accuracy: 1e-10)
        XCTAssertEqual(q1.z, 4.0, accuracy: 1e-10)
        XCTAssertEqual(q1.w, 4.0, accuracy: 1e-10)
    }

    func testDoubleQuaternionScalarSubtractAssign() {
        var q = Quaternion<Double>(x: 5, y: 6, z: 7, w: 8)
        let scalar = 2.0

        q -= scalar

        XCTAssertEqual(q.x, 3.0, accuracy: 1e-10)
        XCTAssertEqual(q.y, 4.0, accuracy: 1e-10)
        XCTAssertEqual(q.z, 5.0, accuracy: 1e-10)
        XCTAssertEqual(q.w, 6.0, accuracy: 1e-10)
    }

    // MARK: - Float Quaternion Compound Tests

    func testFloatQuaternionMultiplyAssign() {
        var q1 = Quaternion<Float>(x: 1, y: 0, z: 0, w: 1)
        let q2 = Quaternion<Float>(x: 0, y: 1, z: 0, w: 1)
        let expected = q1 * q2

        q1 *= q2

        XCTAssertEqual(q1.x, expected.x, accuracy: 1e-5)
        XCTAssertEqual(q1.y, expected.y, accuracy: 1e-5)
        XCTAssertEqual(q1.z, expected.z, accuracy: 1e-5)
        XCTAssertEqual(q1.w, expected.w, accuracy: 1e-5)
    }

    func testFloatQuaternionScalarMultiplyAssign() {
        var q = Quaternion<Float>(x: 1, y: 2, z: 3, w: 4)
        let scalar: Float = 2.5
        let expected = q * scalar

        q *= scalar

        XCTAssertEqual(q.x, expected.x, accuracy: 1e-5)
        XCTAssertEqual(q.y, expected.y, accuracy: 1e-5)
        XCTAssertEqual(q.z, expected.z, accuracy: 1e-5)
        XCTAssertEqual(q.w, expected.w, accuracy: 1e-5)
    }

    func testFloatQuaternionDivideAssign() {
        var q = Quaternion<Float>(x: 2, y: 4, z: 6, w: 8)
        let scalar: Float = 2.0

        q /= scalar

        XCTAssertEqual(q.x, 1.0, accuracy: 1e-5)
        XCTAssertEqual(q.y, 2.0, accuracy: 1e-5)
        XCTAssertEqual(q.z, 3.0, accuracy: 1e-5)
        XCTAssertEqual(q.w, 4.0, accuracy: 1e-5)
    }

    func testFloatQuaternionAddAssign() {
        var q1 = Quaternion<Float>(x: 1, y: 2, z: 3, w: 4)
        let q2 = Quaternion<Float>(x: 5, y: 6, z: 7, w: 8)

        q1 += q2

        XCTAssertEqual(q1.x, 6.0, accuracy: 1e-5)
        XCTAssertEqual(q1.y, 8.0, accuracy: 1e-5)
        XCTAssertEqual(q1.z, 10.0, accuracy: 1e-5)
        XCTAssertEqual(q1.w, 12.0, accuracy: 1e-5)
    }

    func testFloatQuaternionSubtractAssign() {
        var q1 = Quaternion<Float>(x: 5, y: 6, z: 7, w: 8)
        let q2 = Quaternion<Float>(x: 1, y: 2, z: 3, w: 4)

        q1 -= q2

        XCTAssertEqual(q1.x, 4.0, accuracy: 1e-5)
        XCTAssertEqual(q1.y, 4.0, accuracy: 1e-5)
        XCTAssertEqual(q1.z, 4.0, accuracy: 1e-5)
        XCTAssertEqual(q1.w, 4.0, accuracy: 1e-5)
    }

    // MARK: - Chained Compound Operations

    func testChainedCompoundOperations() {
        var q = Quaternion<Double>(x: 1, y: 1, z: 1, w: 1)

        q *= 2.0
        q += Quaternion<Double>(x: 1, y: 1, z: 1, w: 1)
        q -= 1.0

        XCTAssertEqual(q.x, 2.0, accuracy: 1e-10)
        XCTAssertEqual(q.y, 2.0, accuracy: 1e-10)
        XCTAssertEqual(q.z, 2.0, accuracy: 1e-10)
        XCTAssertEqual(q.w, 2.0, accuracy: 1e-10)
    }

    func testMultiplyAssignEquivalence() {
        var q1 = Quaternion<Double>(x: 1, y: 2, z: 3, w: 4).normalized
        var q2 = q1
        let q3 = Quaternion<Double>(x: 5, y: 6, z: 7, w: 8).normalized

        // Test that *= and = * are equivalent
        q1 *= q3
        q2 = q2 * q3

        XCTAssertEqual(q1.x, q2.x, accuracy: 1e-10)
        XCTAssertEqual(q1.y, q2.y, accuracy: 1e-10)
        XCTAssertEqual(q1.z, q2.z, accuracy: 1e-10)
        XCTAssertEqual(q1.w, q2.w, accuracy: 1e-10)
    }
}
