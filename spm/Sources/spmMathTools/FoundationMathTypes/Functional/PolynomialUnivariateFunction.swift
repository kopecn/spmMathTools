//
//  PolynomialUnivariateFunction.swift
//
//
//  Created by Nicholas Bergantz on 4/21/24.
//
import Foundation

/// if the number is within this range round it to zero
///    - e.g. `if polyUniThresh ~= val { return 0 }`
private let polyUniThresh: Double = 1e-9
private let polyUniThreshRange: Range = -polyUniThresh..<polyUniThresh
private let eps16: Double = 16 * .ulpOfOne

extension Array where Element == Double {

    /// Gate allowing only positive values to be appended.
    ///
    /// - Parameters:
    ///        - val: Double, If true, elements will be joined with higher precision formatting.
    fileprivate mutating func insertIfPositive(_ val: Double) {
        if val >= 0 {
            self.append(val)
        }
    }
}

// MARK: - PolynomialUnivariateOrder

/// the various [degrees](https://en.wikipedia.org/wiki/Degree_of_a_polynomial) of a polynomial
///     Special case – zero (see § Degree of the zero polynomial, below)
///     Degree ⁰ – non-zero constant[5]
///     Degree ¹ – linear
///     Degree ² – quadratic
///     Degree ³ – cubic
///     Degree ⁴ – quartic (or, if all terms have even degree, biquadratic)
///     Degree ⁵ – quintic
///     Degree ⁶ – sextic (or, less commonly, hexic)
///     Degree ⁷ – septic (or, less commonly, heptic)
///     Degree ⁸ – octic
///     Degree ⁹ – nonic
///     Degree ¹⁰ – decic
public enum PolynomialUnivariateOrder: Equatable, Sendable {

    case zeroth, constant, linear, quadratic, cubic, quartic, quintic, sextic, septic, octic, nonic,
        decic
    case higher(Int)

    public init(fromIndex: Int) {
        switch fromIndex {
        case ...0: self = .zeroth
        case 1: self = .constant
        case 2: self = .linear
        case 3: self = .quadratic
        case 4: self = .cubic
        case 5: self = .quartic
        case 6: self = .quintic
        case 7: self = .sextic
        case 8: self = .septic
        case 9: self = .octic
        case 10: self = .nonic
        case 11: self = .decic
        case _: self = .higher(fromIndex)
        }
    }

    public var description: String {
        switch self {
        case .zeroth: return "zero univariate polynomial"
        case .constant: return "constant univariate polynomial"
        case .linear: return "linear univariate polynomial"
        case .quadratic: return "quadratic univariate polynomial"
        case .cubic: return "cubic univariate polynomial"
        case .quartic: return "quartic univariate polynomial"
        case .quintic: return "quintic univariate polynomial"
        case .sextic: return "sextic univariate polynomial"
        case .septic: return "septic univariate polynomial"
        case .octic: return "octic univariate polynomial"
        case .nonic: return "nonic univariate polynomial"
        case .decic: return "decic univariate polynomial"
        case .higher(let degree): return "\(degree)-degree univariate polynomial"
        }
    }

    public var degree: Int {
        switch self {
        case .zeroth: return 0
        case .constant: return 0
        case .linear: return 1
        case .quadratic: return 2
        case .cubic: return 3
        case .quartic: return 4
        case .quintic: return 5
        case .sextic: return 6
        case .septic: return 7
        case .octic: return 8
        case .nonic: return 9
        case .decic: return 10
        case .higher(let degree): return degree
        }
    }

    public func construct(coefficients: [Double]) -> PolynomialUnivariateFunction {
        switch self {
        case .zeroth: return ZeroUnivariatePolynomial(coefficients: coefficients)
        case .constant: return ConstantUnivariateFunction(coefficients: coefficients)
        case .linear: return LinearUnivariateFunction(coefficients: coefficients)
        case .quadratic: return QuadraticUnivariatePolynomial(coefficients: coefficients)
        case .cubic: return CubicUnivariatePolynomial(coefficients: coefficients)
        case .quartic: return QuarticUnivariatePolynomial(coefficients: coefficients)
        case .quintic: return QuinticUnivariatePolynomial(coefficients: coefficients)
        case .sextic: return SexticUnivariatePolynomial(coefficients: coefficients)
        case .septic: return SepticUnivariatePolynomial(coefficients: coefficients)
        case .octic: return OcticUnivariatePolynomial(coefficients: coefficients)
        case .nonic: return NonicUnivariatePolynomial(coefficients: coefficients)
        case .decic: return DecicUnivariatePolynomial(coefficients: coefficients)
        case .higher(_): return PolynomialUnivariateFunction(coefficients: coefficients)
        }
    }
}

// MARK: - PolynomialUnivariateFunction

public class PolynomialUnivariateFunction: CustomStringConvertible {

    static let tolerance = 1e-14

    // MARK: - Properties
    let coefficients: [Double]

    /// CustomStringConvertible
    public var description: String { "\(type(of: self)) with coeffs: \(self.coefficients)" }

    /// the [degree](https://en.wikipedia.org/wiki/Degree_of_a_polynomial) of this polynomial
    public let degree: PolynomialUnivariateOrder

    /// returns the derivative of this function, yet another polynomical derivative one order less
    public var derivative: PolynomialUnivariateFunction {
        /// the zeroth' derivative is zero
        if self.degree == .zeroth { return ZeroUnivariatePolynomial() }
        let c = self.coefficients.count - 1

        // 4      3      2      1      0
        // a x³ + b x² + c x¹ + d x⁰ + 0  --> Cubic
        //        3ax² + 2bx¹ + 1cx⁰ + 0  --> Quadratic
        return PolynomialUnivariateOrder(fromIndex: self.coefficients.count - 1)
            .construct(
                coefficients: self.coefficients.enumerated().compactMap {
                    c == $0 ? nil : Double(c - $0) * $1
                }
            )
    }

    /// roots, get the roots of this polynomial
    var roots: [Double] {
        get throws {
            throw MathErrors.rootsUnsolveableFor(self.degree, self.coefficients)
        }
    }

    // MARK: - Initializers

    /// constructs using an array of coeff's
    public required init(coefficients: [Double]) {
        self.coefficients = coefficients
        self.degree = PolynomialUnivariateOrder(fromIndex: coefficients.count)
    }

    // MARK: - Methods

    /// returns the coefficient for said degree,
    public func coeffAt(_ deg: Int) -> Double {
        self.coefficients[deg]
    }

    /// This will evaluate the result of the polynomial function for f(x) ...
    /// Uses Horner's method for efficient polynomial evaluation
    public func evaluate(_ x: Double) -> Double {
        if self.degree == .zeroth { return 0 }

        // Horner's method: evaluates polynomial by factoring as a₀ + x(a₁ + x(a₂ + ...))
        // Coefficients are stored as [aₙ, aₙ₋₁, ..., a₁, a₀] (highest degree first)
        guard !coefficients.isEmpty else { return 0 }

        var result = coefficients[0]
        for i in 1..<coefficients.count {
            result = result * x + coefficients[i]
        }
        return result
    }

    /// returns the integral of this function, yet another polynomical derivative one order more
    public func integrate(_ c: Double) -> PolynomialUnivariateFunction {
        /// the zeroth' derivative is zero
        if self.degree == .zeroth { return ConstantUnivariateFunction(c) }
        let co = self.coefficients.count
        var coeffs = self.coefficients.enumerated().map {
            $1 / Double(co - $0)
        }
        coeffs.append(c)
        //           4         3        2       1      0
        //           a x³   +  b x²   +  cx¹ +  dx⁰ +  0  --> Cubic
        // ax⁴/4  +  bx³/3  +  cx²/2  +  dx¹  + ex⁰ +  0  --> Quartic
        return PolynomialUnivariateOrder(fromIndex: coeffs.count)
            .construct(coefficients: coeffs)
    }

}

// MARK: - ZeroUnivariatePolynomial

/// Zeroth Polynomial of the form: 0
public class ZeroUnivariatePolynomial: PolynomialUnivariateFunction {
    // MARK: - Properties

    /// the roots for the zeroth polynomial is none
    override var roots: [Double] { [] }

    // MARK: - Initializers

    /// constructs from the form: 0
    public init() {
        super.init(coefficients: [])
    }

    public required init(coefficients: [Double]) {
        super.init(coefficients: coefficients)
    }
}

// MARK: - ConstantUnivariateFunction

/// Constant Polynomial of the form: a x⁰
public class ConstantUnivariateFunction: PolynomialUnivariateFunction {
    // MARK: - Properties

    /// the real roots for a constant function is zero if the constant is zero
    override var roots: [Double] {
        if a == 0 { return [0] }
        return []
    }

    /// The constant coeff for the polynomial
    var a: Double { self.coefficients[0] }

    // MARK: - Initializers

    /// constructs from the form: a x⁰
    public init(
        _ a: Double
    ) {
        super.init(coefficients: [a])
    }

    public required init(coefficients: [Double]) {
        super.init(coefficients: coefficients)
    }
}

// MARK: - LinearUnivariateFunction

/// Linear Polynomial of the form: a x¹ + b x⁰
public class LinearUnivariateFunction: PolynomialUnivariateFunction {
    // MARK: - Properties

    /// the roots for the linear equation is b/m = x
    override var roots: [Double] {
        if a == 0 {
            return ConstantUnivariateFunction(b).roots
        }
        return [-b / a]
    }

    /// The constant coeff for the polynomial
    var b: Double { self.coefficients[1] }

    /// The linear coeff for the polynomial
    var a: Double { self.coefficients[0] }

    // MARK: - Initializers

    /// constructs from the form: a x¹ + b x⁰
    public init(
        _ a: Double,
        _ b: Double
    ) {
        super.init(coefficients: [a, b])
    }

    public required init(coefficients: [Double]) {
        super.init(coefficients: coefficients)
    }
}

// MARK: - QuadraticUnivariatePolynomial

/// Quadratic Polynomial of the form: a x² + b x¹ + c x⁰
public class QuadraticUnivariatePolynomial: PolynomialUnivariateFunction {
    // MARK: - Properties

    /// the real roots for the quadratic equation is b/m = x
    ///     - https://kieranb662.github.io/blog/2020/04/21/Swift-Equation-Solvers
    override var roots: [Double] {
        if a == 0 {
            return LinearUnivariateFunction(b, c).roots
        }

        let d = pow(b, 2) - 4 * a * c

        switch d {
        case ...polyUniThresh:
            return [-b / (2 * a)]

        case 0...:
            return [
                (-b + sqrt(d)) / (2 * a),
                (-b - sqrt(d)) / (2 * a),
            ]

        default:
            return []
        }
    }

    /// The constant coeff for the polynomial
    var c: Double { self.coefficients[2] }

    /// The linear coeff for the polynomial
    var b: Double { self.coefficients[1] }

    /// The quadratic coeff for the polynomial
    var a: Double { self.coefficients[0] }

    // MARK: - Initializers

    /// constructs from the form: a x² + b x¹ + c x⁰
    public init(
        _ a: Double,
        _ b: Double,
        _ c: Double
    ) {
        super.init(coefficients: [a, b, c])
    }

    public required init(coefficients: [Double]) {
        super.init(coefficients: coefficients)
    }
}

// MARK: - CubicUnivariatePolynomial

/// Cubic Polynomial of the form: a x³ + b x² + c x¹ + d x⁰
public class CubicUnivariatePolynomial: PolynomialUnivariateFunction {
    // MARK: - Properties

    /// the real roots for the quadratic equation is b/m = x
    ///     - https://kieranb662.github.io/blog/2020/04/21/Swift-Equation-Solvers
    override var roots: [Double] {
        if a == 0 {
            return QuadraticUnivariatePolynomial(b, c, d).roots
        }

        let a_1 = b / a
        let a_2 = c / a
        let a_3 = d / a

        let q = (3 * a_2 - pow(a_1, 2)) / 9
        let r = (9 * a_1 * a_2 - 27 * a_3 - 2 * pow(a_1, 3)) / 54

        let s = cbrt(r + sqrt(pow(q, 3) + pow(r, 2)))
        let t = cbrt(r - sqrt(pow(q, 3) + pow(r, 2)))

        let d = pow(q, 3) + pow(r, 2)  // discriminant

        switch d {
        case ...polyUniThresh:
            let theta = acos(r / sqrt(-pow(q, 3)))
            let twoSqrtQ = 2 * sqrt(-q)
            return [
                twoSqrtQ * cos((1 / 3) * theta) - (1 / 3) * a_1,
                twoSqrtQ * cos((1 / 3) * theta + 2 * Double.pi / 3) - (1 / 3) * a_1,
                twoSqrtQ * cos((1 / 3) * theta + 4 * Double.pi / 3) - (1 / 3) * a_1,
            ]

        case 0...:
            return [
                s + t - (1 / 3) * a_1,
                -(1 / 2) * (s + t) - (1 / 3) * a_1,
                -(1 / 2) * (s + t) - (1 / 3) * a_1,
            ]

        default:
            return []
        }
    }

    /// The constant coeff for the polynomial
    var d: Double { self.coefficients[3] }

    /// The linear coeff for the polynomial
    var c: Double { self.coefficients[2] }

    /// The quadratic coeff for the polynomial
    var b: Double { self.coefficients[1] }

    /// The cubic coeff for the polynomial
    var a: Double { self.coefficients[0] }

    // MARK: - Initializers

    /// constructs from the form: a x³ + b x² + c x¹ + d x⁰
    public init(
        _ a: Double,
        _ b: Double,
        _ c: Double,
        _ d: Double
    ) {
        super.init(coefficients: [a, b, c, d])
    }

    public required init(coefficients: [Double]) {
        super.init(coefficients: coefficients)
    }

    func cubicSolve(
        a: Double,
        b: Double,
        c: Double,
        d: Double,
        threshold: Double = 0.0001
    ) -> [Double] {

        var roots = [Double]()

        let a_1 = b / a
        let a_2 = c / a
        let a_3 = d / a

        let q = (3 * a_2 - pow(a_1, 2)) / 9
        let r = (9 * a_1 * a_2 - 27 * a_3 - 2 * pow(a_1, 3)) / 54

        let s = cbrt(r + sqrt(pow(q, 3) + pow(r, 2)))
        let t = cbrt(r - sqrt(pow(q, 3) + pow(r, 2)))

        var d = pow(q, 3) + pow(r, 2)  // discriminant

        // Check if d is within the zero threshold
        if -threshold < d && d < threshold { d = 0 }
        if d > 0 {
            let x_1 = s + t - (1 / 3) * a_1
            let x_2 = -(1 / 2) * (s + t) - (1 / 3) * a_1
            let x_3 = -(1 / 2) * (s + t) - (1 / 3) * a_1
            roots = [x_1, x_2, x_3]

        } else if d <= 0 {
            let theta = acos(r / sqrt(-pow(q, 3)))
            let x_1 = 2 * sqrt(-q) * cos((1 / 3) * theta) - (1 / 3) * a_1
            let x_2 = 2 * sqrt(-q) * cos((1 / 3) * theta + 2 * Double.pi / 3) - (1 / 3) * a_1
            let x_3 = 2 * sqrt(-q) * cos((1 / 3) * theta + 4 * Double.pi / 3) - (1 / 3) * a_1
            roots = [x_1, x_2, x_3]
        }
        return roots
    }
}

// MARK: - QuarticUnivariatePolynomial

/// Quartic Polynomial of the form: a x⁴ + b x³ + c x² + d x¹ + e x⁰
public class QuarticUnivariatePolynomial: PolynomialUnivariateFunction {
    // MARK: - Properties

    /// the real roots for the quadratic equation is b/m = x
    ///     - https://kieranb662.github.io/blog/2020/04/21/Swift-Equation-Solvers
    override var roots: [Double] {
        get throws {
            if a == 0 {
                return CubicUnivariatePolynomial(b, c, d, e).roots
            }
            if e != 0 {
                throw MathErrors.rootsUnsolveableFor(self.degree, self.coefficients)
            }
            return self.solveQuartMonic
        }
    }

    /// The constant coeff for the polynomial
    var e: Double { self.coefficients[4] }

    /// The linear coeff for the polynomial
    var d: Double { self.coefficients[3] }

    /// The quadratic coeff for the polynomial
    var c: Double { self.coefficients[2] }

    /// The cubic coeff for the polynomial
    var b: Double { self.coefficients[1] }

    /// The quartic coeff for the polynomial
    var a: Double { self.coefficients[0] }

    // MARK: - Initializers

    /// constructs from the form: a x⁴ + b x³ + c x² + d x¹ + e x⁰
    public init(
        _ a: Double,
        _ b: Double,
        _ c: Double,
        _ d: Double,
        _ e: Double
    ) {
        super.init(coefficients: [a, b, c, d, e])
    }

    public required init(coefficients: [Double]) {
        super.init(coefficients: coefficients)
    }

    /// Solves a (quartic polynomial)[https://en.wikipedia.org/wiki/Quartic_function] with real coefficients for its real roots.
    ///
    /// Calculate all roots of the [monic quartic equation](https://en.wikipedia.org/wiki/Quartic_equation): x^4 + a*x^3 + b*x^2 + c*x + d = 0
    ///
    /// - Parameters:
    ///     - aS: Coefficient of x^3 term in original quartic polynomial.
    ///     - bS: Coefficient of x^2 term in original quartic polynomial.
    ///     - cS: Coefficient of x term in original quartic polynomial.
    ///     - dS: Constant term in original quartic polynomial.
    /// - Returns: Array containing all real roots of the original quartic polynomial.
    ///
    private var solveQuartMonic: [Double] {

        let aM = b / a
        let bM = c / a
        let cM = d / a
        let dM = e / a

        var roots: [Double] = []

        if abs(dM) < .ulpOfOne {
            if abs(cM) < .ulpOfOne {
                roots.insertIfPositive(0.0)

                let D = aM * aM - 4 * bM
                if abs(D) < .ulpOfOne {
                    roots.insertIfPositive(-aM / 2)
                } else if D > 0.0 {
                    let sqrtD = sqrt(D)
                    roots.insertIfPositive((-aM - sqrtD) / 2)
                    roots.insertIfPositive((-aM + sqrtD) / 2)
                }
                return roots
            }

            if abs(aM) < .ulpOfOne && abs(bM) < .ulpOfOne {
                roots.insertIfPositive(0.0)
                roots.insertIfPositive(-cbrt(cM))
                return roots
            }
        }

        let a3 = -bM
        let b3 = aM * cM - 4 * dM
        let c3 = -aM * aM * dM - cM * cM + 4 * bM * dM

        var x3: [Double] = [0, 0, 0]
        let number_zeroes = solveResolventQuartic(&x3, a3, b3, c3)

        var y = x3[0]
        // Choosing Y with maximal absolute value.
        if number_zeroes != 1 {
            if abs(x3[1]) > abs(y) {
                y = x3[1]
            }
            if abs(x3[2]) > abs(y) {
                y = x3[2]
            }
        }

        var q1: Double
        var q2: Double
        var p1: Double
        var p2: Double

        var D = y * y - 4 * dM
        if abs(D) < .ulpOfOne {
            q1 = y / 2
            q2 = q1
            D = aM * aM - 4 * (bM - y)
            if abs(D) < .ulpOfOne {
                p1 = aM / 2
                p2 = p1
            } else {
                let sqrtD = sqrt(D)
                p1 = (aM + sqrtD) / 2
                p2 = (aM - sqrtD) / 2
            }
        } else {
            let sqrtD = sqrt(D)
            q1 = (y + sqrtD) / 2
            q2 = (y - sqrtD) / 2
            p1 = (aM * q1 - cM) / (q1 - q2)
            p2 = (cM - aM * q2) / (q1 - q2)
        }

        D = p1 * p1 - 4 * q1
        if abs(D) < eps16 {
            roots.insertIfPositive(-p1 / 2)
        } else if D > 0.0 {
            let sqrtD = sqrt(D)
            roots.insertIfPositive((-p1 - sqrtD) / 2)
            roots.insertIfPositive((-p1 + sqrtD) / 2)
        }

        D = p2 * p2 - 4 * q2
        if abs(D) < eps16 {
            roots.insertIfPositive(-p2 / 2)
        } else if D > 0.0 {
            let sqrtD = sqrt(D)
            roots.insertIfPositive((-p2 - sqrtD) / 2)
            roots.insertIfPositive((-p2 + sqrtD) / 2)
        }

        return roots
    }

    /// Solve the resolvent cubic equation corresponding to a quartic polynomial.
    ///
    /// This function solves the resolvent cubic equation corresponding to the input quartic polynomial, whose coefficients are passed as parameters. It stores the three roots of the resolvent cubic in the passed x array, and returns the number of real roots.
    ///
    /// - Parameters:
    ///   - x: Array to store the three roots of the resolvent cubic. Passed by reference to allow modification.
    ///   - aS: Coefficient of x^3 term in original quartic polynomial.
    ///   - bS: Coefficient of x^2 term in original quartic polynomial.
    ///   - cS: Coefficient of x term in original quartic polynomial.
    /// - Returns: The number of real roots (1, 2 or 3).
    private func solveResolventQuartic(
        _ x: inout [Double],
        _ aS: Double,
        _ bS: Double,
        _ cS: Double
    ) -> Int {

        let aR = aS / 3
        let bR = bS
        let cR = cS

        let a2 = aR * aR
        let r = (aR * (2 * a2 - bR) + cR) / 2
        let r2 = r * r
        var q = a2 - bR / 3
        let q3 = q * q * q

        if r2 < q3 {
            let qsqrt = sqrt(q)
            let t = min(max(r / (q * qsqrt), -1.0), 1.0)
            q = -2 * qsqrt

            let theta = acos(t) / 3
            let ux = cos(theta) * q
            let uyi = sin(theta) * q
            x[0] = ux - aR
            x[1] = ux * cos120 - uyi * sin120 - aR
            x[2] = ux * cos120 + uyi * sin120 - aR
            return 3

        }

        var A = -cbrt(abs(r) + sqrt(r2 - q3))
        if r < 0.0 {
            A = -A
        }
        let B = (0.0 == A ? 0.0 : q / A)

        x[0] = (A + B) - aR
        x[1] = -(A + B) / 2 - aR
        x[2] = sqrt(3) * (A - B) / 2

        if abs(x[2]) < .ulpOfOne {
            x[2] = x[1]
            return 2
        }

        return 1
    }

}

// MARK: - QuinticUnivariatePolynomial

/// Quintic Polynomial of the form: a x⁵ + b x⁴ + c x³ + d x² + e x¹ + f x⁰
///     Notes on the roots solution for the quintic form and higher:
///         - https://en.wikipedia.org/wiki/Quintic_function
///         - https://mathworld.wolfram.com/QuinticEquation.html
///         - https://en.wikipedia.org/wiki/Bring_radical
public class QuinticUnivariatePolynomial: PolynomialUnivariateFunction {
    // MARK: - Initializers

    /// constructs from the form: a x⁵ + b x⁴ + c x³ + d x² + e x¹ + f x⁰
    public init(
        _ a: Double,
        _ b: Double,
        _ c: Double,
        _ d: Double,
        _ e: Double,
        _ f: Double
    ) {
        super.init(coefficients: [a, b, c, d, e, f])
    }

    public required init(coefficients: [Double]) {
        super.init(coefficients: coefficients)
    }
}

// MARK: - SexticUnivariatePolynomial

/// Sextic Polynomial of the form: a x⁶ + b x⁵ + c x⁴ + d x³ + e x² + f x¹ + g x⁰
public class SexticUnivariatePolynomial: PolynomialUnivariateFunction {
    // MARK: - Initializers

    /// constructs from the form: a x⁶ + b x⁵ + c x⁴ + d x³ + e x² + f x¹ + g x⁰
    public init(
        _ a: Double,
        _ b: Double,
        _ c: Double,
        _ d: Double,
        _ e: Double,
        _ f: Double,
        _ g: Double
    ) {
        super.init(coefficients: [a, b, c, d, e, f, g])
    }

    public required init(coefficients: [Double]) {
        super.init(coefficients: coefficients)
    }
}

// MARK: - SepticUnivariatePolynomial

/// Septic Polynomial of the form: a x⁷ + b x⁶ + c x⁵ + d x⁴ + e x³ + f x² + g x¹ + h x⁰
public class SepticUnivariatePolynomial: PolynomialUnivariateFunction {
    // MARK: - Initializers

    /// constructs from the form: a x⁷ + b x⁶ + c x⁵ + d x⁴ + e x³ + f x² + g x¹ + h x⁰
    public init(
        _ a: Double,
        _ b: Double,
        _ c: Double,
        _ d: Double,
        _ e: Double,
        _ f: Double,
        _ g: Double,
        _ h: Double
    ) {
        super.init(coefficients: [a, b, c, d, e, f, g, h])
    }

    public required init(coefficients: [Double]) {
        super.init(coefficients: coefficients)
    }
}

// MARK: - OcticUnivariatePolynomial

/// Octic Polynomial of the form: a x⁸ + b x⁷ + c x⁶ + d x⁵ + e x⁴ + f x³ + g x² + h x¹ + i x⁰
public class OcticUnivariatePolynomial: PolynomialUnivariateFunction {
    // MARK: - Initializers

    /// constructs from the form: a x⁸ + b x⁷ + c x⁶ + d x⁵ + e x⁴ + f x³ + g x² + h x¹ + i x⁰
    public init(
        _ a: Double,
        _ b: Double,
        _ c: Double,
        _ d: Double,
        _ e: Double,
        _ f: Double,
        _ g: Double,
        _ h: Double,
        _ i: Double
    ) {
        super.init(coefficients: [a, b, c, d, e, f, g, h, i])
    }

    public required init(coefficients: [Double]) {
        super.init(coefficients: coefficients)
    }
}

// MARK: - NonicUnivariatePolynomial

/// Nonic Polynomial of the form: a x⁹ + b x⁸ + c x⁷ + d x⁶ + e x⁵ + f x⁴ + g x³ + h x² + i x¹ + j x⁰
public class NonicUnivariatePolynomial: PolynomialUnivariateFunction {
    // MARK: - Initializers

    /// constructs from the form: a x⁹ + b x⁸ + c x⁷ + d x⁶ + e x⁵ + f x⁴ + g x³ + h x² + i x¹ + j x⁰
    public init(
        _ a: Double,
        _ b: Double,
        _ c: Double,
        _ d: Double,
        _ e: Double,
        _ f: Double,
        _ g: Double,
        _ h: Double,
        _ i: Double,
        _ j: Double
    ) {
        super.init(coefficients: [a, b, c, d, e, f, g, h, i, j])
    }

    public required init(coefficients: [Double]) {
        super.init(coefficients: coefficients)
    }
}

// MARK: - DecicUnivariatePolynomial

/// Decic Polynomial of the form: a x¹⁰ + b x⁹ + c x⁸ + d x⁷ + e x⁶ + f x⁵ + g x⁴ + h x³ + i x² + j x¹ + k x⁰
public class DecicUnivariatePolynomial: PolynomialUnivariateFunction {
    // MARK: - Initializers

    /// constructs from the form: a x¹⁰ + b x⁹ + c x⁸ + d x⁷ + e x⁶ + f x⁵ + g x⁴ + h x³ + i x² + j x¹ + k x⁰
    public init(
        _ a: Double,
        _ b: Double,
        _ c: Double,
        _ d: Double,
        _ e: Double,
        _ f: Double,
        _ g: Double,
        _ h: Double,
        _ i: Double,
        _ j: Double,
        _ k: Double
    ) {
        super.init(coefficients: [a, b, c, d, e, f, g, h, i, j, k])
    }

    public required init(coefficients: [Double]) {
        super.init(coefficients: coefficients)
    }
}

///⁰ ₀
///¹ ₁
///² ₂
///³ ₃
///⁴ ₄
///⁵ ₅
///⁶ ₆
///⁷ ₇
///⁸ ₈
///⁹ ₉
