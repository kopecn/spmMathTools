import XCTest
import spmMathTools

@testable import spmMathTools

/// Testing Roots
let quarticRootTruthTable:
    [(
        coeff: [Double],
        roots: [Double]
    )] = [
        ([1, 2, 3, 4, -1], [0.21129, -1.76565])
    ]

let cubicRootTruthTable:
    [(
        coeff: [Double],
        roots: [Double]
    )] = [
        ([1, 2, 3, 4], [-1.65063, -0.17469, -0.17469]),
        ([3, 6, 2, 1], [-1.72560]),
        ([1 / 3, 6, 2, 1], [-17.67005]),
        ([1 / 3, 6, 22, 1], [-0.04603, -12.90288, -5.05108]),
        ([1 / 3, 6, 22, 1.21443], [-0.05605, -12.90924, -5.03469]),
        ([1 / 3, 6.2342, 22, 1.21443], [-0.05609, -14.01036, -4.63614]),
        ([1 / 3, 6.2342, -22, 1.21443], [2.98673, -21.74542, 0.05609]),
    ]

let quadraticRootTruthTable:
    [(
        coeff: [Double],
        roots: [Double]
    )] = [
        ([1, 2, 3], [-1]),
        ([-1, 2, 3], [-1, 3]),
        ([4, 9, 3], [-0.40693, -1.84307]),
        ([0.25, 9, 3], [-0.33647, -35.66352]),
        ([0.25, 1 / 3, 8 / 7], [-2 / 3]),
        ([0.25, 18.703, 1 / 7], [-0.00763, -74.80436]),
    ]

final class PolynomialUnivariateFunctionTests: XCTestCase {
    func testPolynomialFunctionDerivatives() throws {
        let cubic = CubicUnivariatePolynomial(coefficients: [1, 2, 3, 4])
        XCTAssertEqual(cubic.degree, .cubic)
        XCTAssertEqual(cubic.evaluate(2), 26.0)

        let quadratic = cubic.derivative
        XCTAssertEqual(quadratic.degree, .quadratic)
        XCTAssert(quadratic is QuadraticUnivariatePolynomial)
        XCTAssertEqual(quadratic.evaluate(2), 23.0)

        let linear = quadratic.derivative
        XCTAssertEqual(linear.degree, .linear)
        XCTAssert(linear is LinearUnivariateFunction)
        XCTAssertEqual(linear.evaluate(2), 16.0)
        print(linear)

        let constant = linear.derivative
        XCTAssertEqual(constant.degree, .constant)
        XCTAssert(constant is ConstantUnivariateFunction)
        XCTAssertEqual(constant.evaluate(2), 6.0)

        let zeroth = constant.derivative
        XCTAssertEqual(zeroth.degree, .zeroth)
        XCTAssert(zeroth is ZeroUnivariatePolynomial)
        XCTAssertEqual(zeroth.evaluate(2), 0)

        let zerothN = zeroth.derivative
        XCTAssertEqual(zerothN.degree, .zeroth)
        XCTAssert(zerothN is ZeroUnivariatePolynomial)
    }

    func testPolynomialFunctionIntegration() throws {

        let zeroth = ZeroUnivariatePolynomial()

        let constant = zeroth.integrate(6)
        XCTAssertEqual(constant.degree, .constant)
        XCTAssert(constant is ConstantUnivariateFunction)
        XCTAssertEqual(constant.evaluate(2), 6.0)

        let linear = constant.integrate(4.0)
        XCTAssertEqual(linear.degree, .linear)
        XCTAssert(linear is LinearUnivariateFunction)
        XCTAssertEqual(linear.evaluate(2), 16.0)

        let quadratic = linear.integrate(3.0)
        XCTAssertEqual(quadratic.degree, .quadratic)
        XCTAssert(quadratic is QuadraticUnivariatePolynomial)
        XCTAssertEqual(quadratic.evaluate(2), 23.0)
        print(quadratic)

        let cubic = quadratic.integrate(4)
        XCTAssertEqual(cubic.degree, .cubic)
        XCTAssert(cubic is CubicUnivariatePolynomial)
        XCTAssertEqual(cubic.evaluate(2), 26.0)
    }

    func testPolynomialFunctionSolveRoots() throws {

        for r in quadraticRootTruthTable {
            let quad = QuadraticUnivariatePolynomial(coefficients: r.coeff)
            for (qr, rr) in zip(quad.roots, r.roots) {
                XCTAssertEqual(qr, rr, accuracy: 0.0001)
            }
        }

        for r in cubicRootTruthTable {
            let cubic = CubicUnivariatePolynomial(coefficients: r.coeff)
            for (qr, rr) in zip(cubic.roots, r.roots) {
                XCTAssertEqual(qr, rr, accuracy: 0.0001)
            }
        }

        for r in quarticRootTruthTable {
            let quartic = QuarticUnivariatePolynomial(coefficients: r.coeff)
            if let roots = try? quartic.roots {
                for (qr, rr) in zip(roots, r.roots) {
                    XCTAssertEqual(qr, rr, accuracy: 0.0001)
                }
            }
        }
    }
}
