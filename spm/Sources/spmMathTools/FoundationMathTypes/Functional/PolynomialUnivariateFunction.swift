//
//  PolynomialUnivariateFunction.swift
//
//
//  Created by Nicholas Bergantz on 4/21/24.
//
import Foundation

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
        case ...(-1): self = .zeroth
        case 0: self = .constant
        case 1: self = .linear
        case 2: self = .quadratic
        case 3: self = .cubic
        case 4: self = .quartic
        case 5: self = .quintic
        case 6: self = .sextic
        case 7: self = .septic
        case 8: self = .octic
        case 9: self = .nonic
        case 10: self = .decic
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

    // MARK: - Properties
    let coefficients: [Double]

    /// CustomStringConvertible
    public var description: String { "\(type(of: self)) with coeffs: \(self.coefficients)" }

    /// the [degree](https://en.wikipedia.org/wiki/Degree_of_a_polynomial) of this polynomial
    public let degree: PolynomialUnivariateOrder

    /// returns the derivative of this function, yet another polynomical derivative one order less
    public var derivative: PolynomialUnivariateFunction {
        /// the zeroth derivative is zero
        if self.degree == .zeroth { return ZeroUnivariatePolynomial() }
        let c = self.coefficients.count - 1

        // Calculate derivative coefficients
        let derivCoeffs = self.coefficients.enumerated().compactMap { (index, coeff) in
            c == index ? nil : Double(c - index) * coeff
        }

        // Fix: The derivative degree should be one less than the original
        return PolynomialUnivariateOrder(fromIndex: derivCoeffs.count - 1)
            .construct(coefficients: derivCoeffs)
    }

    /// roots, get the roots of this polynomial
    public var roots: [Double] {
        get throws {
            throw MathErrors.rootsUnsolveableFor(self.degree, self.coefficients)
        }
    }

    // MARK: - Initializers

    /// constructs using an array of coeff's
    public required init(coefficients: [Double]) {
        self.coefficients = coefficients
        // Fix: degree should be coefficients.count - 1, not coefficients.count
        self.degree = PolynomialUnivariateOrder(fromIndex: coefficients.count - 1)
    }

    // MARK: - Methods

    /// Returns the coefficient for the specified degree with bounds checking
    /// - Parameter deg: The degree index (0 = constant term, 1 = linear, etc.)
    /// - Returns: The coefficient at the specified degree, or 0 if out of bounds
    public func coeffAt(_ deg: Int) -> Double {
        guard deg >= 0 && deg < coefficients.count else { return 0 }
        return coefficients[deg]
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
        return PolynomialUnivariateOrder(fromIndex: coeffs.count - 1)
            .construct(coefficients: coeffs)
    }

}

// MARK: - ZeroUnivariatePolynomial

/// Zeroth Polynomial of the form: 0
public class ZeroUnivariatePolynomial: PolynomialUnivariateFunction {
    // MARK: - Properties

    /// the roots for the zeroth polynomial is none
    public override var roots: [Double] { [] }

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
    public override var roots: [Double] {
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
    public override var roots: [Double] {
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
    public override var roots: [Double] {
        if a == 0 {
            return LinearUnivariateFunction(b, c).roots
        }

        let d = pow(b, 2) - 4 * a * c

        switch d {
        case ...polynomialZeroThreshold:
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
    public override var roots: [Double] {
        if a == 0 {
            return QuadraticUnivariatePolynomial(b, c, d).roots
        }

        // Use the existing solveCubic function from Roots.swift
        return solveCubic(a, b, c, d)
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
}

// MARK: - QuarticUnivariatePolynomial

/// Quartic Polynomial of the form: a x⁴ + b x³ + c x² + d x¹ + e x⁰
public class QuarticUnivariatePolynomial: PolynomialUnivariateFunction {
    // MARK: - Properties

    /// the real roots for the quartic equation
    ///     - https://kieranb662.github.io/blog/2020/04/21/Swift-Equation-Solvers
    public override var roots: [Double] {
        if a == 0 {
            return CubicUnivariatePolynomial(b, c, d, e).roots
        }

        // Use the existing solveQuarticMonic function from Roots.swift
        let aM = b / a
        let bM = c / a
        let cM = d / a
        let dM = e / a

        return solveQuarticMonic(aM, bM, cM, dM)
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
