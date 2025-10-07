//
//  Roots.swift
//
//
//  Created by Nicholas Bergantz on 3/23/24.
//

import Foundation

public let eps16: Double = 16 * .ulpOfOne
private let tolerance: Double = 1e-14
private let maxIts: Int = 128

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

/// Calculate all real roots of a cubic polynomial equation a*x^3 + b*x^2 + c*x + d = 0.
///
/// - Parameters:
///   - aS: Coefficient of x^3 term.
///   - bS: Coefficient of x^2 term.
///   - cS: Coefficient of x term.
///   - dS: Constant term.
/// - Returns: Array containing all real roots.
public func solveCubic(
    _ aS: Double,
    _ bS: Double,
    _ cS: Double,
    _ dS: Double
) -> [Double] {
    var a = aS
    var b = bS
    var c = cS
    var d = dS

    var roots: [Double] = []

    if abs(d) < .ulpOfOne {
        // First solution is x = 0
        roots.insertIfPositive(0.0)

        // Converting to a quadratic equation
        d = c
        c = b
        b = a
        a = 0.0
    }

    if abs(a) < .ulpOfOne {
        if abs(b) < .ulpOfOne {
            // Linear equation
            if abs(c) > .ulpOfOne {
                roots.insertIfPositive(-d / c)
            }

        } else {
            // Quadratic equation
            let discriminant = c * c - 4 * b * d
            if discriminant >= 0 {
                let inv2b = 1.0 / (2 * b)
                let y = sqrt(discriminant)
                roots.insertIfPositive((-c + y) * inv2b)
                roots.insertIfPositive((-c - y) * inv2b)
            }
        }

    } else {
        // Cubic equation
        let inva = 1.0 / a
        let invaa = inva * inva
        let bb = b * b
        let bover3a = b * inva / 3
        let p = (a * c - bb / 3) * invaa
        let halfq = (2 * bb * b - 9 * a * b * c + 27 * a * a * d) / 54 * invaa * invaa
        let yy = p * p * p / 27 + halfq * halfq

        if yy > .ulpOfOne {
            // Sqrt is positive: one real solution
            let y = sqrt(yy)
            let uuu = -halfq + y
            let vvv = -halfq - y
            let www = abs(uuu) > abs(vvv) ? uuu : vvv
            let w = cbrt(www)
            roots.insertIfPositive(w - p / (3 * w) - bover3a)
        } else if yy < -.ulpOfOne {
            // Sqrt is negative: three real solutions
            let x = -halfq
            let y = sqrt(-yy)
            var theta: Double
            var r: Double

            // Convert to polar form
            if abs(x) > .ulpOfOne {
                theta = (x > 0.0) ? atan(y / x) : (atan(y / x) + .pi)
                r = sqrt(x * x - yy)
            } else {
                // Vertical line
                theta = .pi / 2
                r = y
            }
            // Calculate cube root
            theta /= 3
            r = 2 * cbrt(r)
            // Convert to complex coordinate
            let ux = cos(theta) * r
            let uyi = sin(theta) * r

            roots.insertIfPositive(ux - bover3a)
            roots.insertIfPositive(ux * cos120 - uyi * sin120 - bover3a)
            roots.insertIfPositive(ux * cos120 + uyi * sin120 - bover3a)
        } else {
            // Sqrt is zero: two real solutions
            let www = -halfq
            let w = 2 * cbrt(www)

            roots.insertIfPositive(w - bover3a)
            roots.insertIfPositive(w * cos120 - bover3a)
        }
    }
    return roots
}

/// Solve the resolvent cubic equation corresponding to a quartic polynomial.
///
/// This function solves the resolvent cubic equation corresponding to the input quartic polynomial, whose coefficients are passed as parameters. It stores the three roots of the resolvent cubic in the passed x array, and returns the number of real roots.
///
/// - Parameters:
///   - x: Array to store the three roots of the resolvent cubic. Must have at least 3 elements.
///   - aS: Coefficient of x^3 term in original quartic polynomial.
///   - bS: Coefficient of x^2 term in original quartic polynomial.
///   - cS: Coefficient of x term in original quartic polynomial.
/// - Returns: The number of real roots (1, 2 or 3).
public func solveResolvent(
    _ x: inout [Double],
    _ aS: Double,
    _ bS: Double,
    _ cS: Double
) -> Int {

    // Ensure array has at least 3 elements
    while x.count < 3 {
        x.append(0.0)
    }

    var a = aS
    let b = bS
    let c = cS

    a /= 3
    let a2 = a * a
    var q = a2 - b / 3
    let r = (a * (2 * a2 - b) + c) / 2
    let r2 = r * r
    let q3 = q * q * q

    if r2 < q3 {
        let qsqrt = sqrt(q)
        let t = min(max(r / (q * qsqrt), -1.0), 1.0)
        q = -2 * qsqrt

        let theta = acos(t) / 3
        let ux = cos(theta) * q
        let uyi = sin(theta) * q
        x[0] = ux - a
        x[1] = ux * cos120 - uyi * sin120 - a
        x[2] = ux * cos120 + uyi * sin120 - a
        return 3

    } else {
        var A = -cbrt(abs(r) + sqrt(r2 - q3))
        if r < 0.0 {
            A = -A
        }
        let B = (0.0 == A ? 0.0 : q / A)

        x[0] = (A + B) - a
        x[1] = -(A + B) / 2 - a
        x[2] = sqrt(3) * (A - B) / 2
        if abs(x[2]) < .ulpOfOne {
            x[2] = x[1]
            return 2
        }

        return 1
    }
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
public func solveQuarticMonic(
    _ a: Double,
    _ b: Double,
    _ c: Double,
    _ d: Double
) -> [Double] {
    var roots: [Double] = []

    if abs(d) < .ulpOfOne {
        if abs(c) < .ulpOfOne {
            roots.insertIfPositive(0.0)

            let D = a * a - 4 * b
            if abs(D) < .ulpOfOne {
                roots.insertIfPositive(-a / 2)
            } else if D > 0.0 {
                let sqrtD = sqrt(D)
                roots.insertIfPositive((-a - sqrtD) / 2)
                roots.insertIfPositive((-a + sqrtD) / 2)
            }
            return roots
        }

        if abs(a) < .ulpOfOne && abs(b) < .ulpOfOne {
            roots.insertIfPositive(0.0)
            roots.insertIfPositive(-cbrt(c))
            return roots
        }
    }

    let a3 = -b
    let b3 = a * c - 4 * d
    let c3 = -a * a * d - c * c + 4 * b * d

    // Fix: Initialize with 3 zeros instead of empty array
    var x3: [Double] = [0.0, 0.0, 0.0]
    let numberOfZeroes = solveResolvent(&x3, a3, b3, c3)
    
    var y = x3[0]
    // Choosing Y with maximal absolute value.
    if numberOfZeroes != 1 {
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

    var D = y * y - 4 * d
    if abs(D) < .ulpOfOne {
        q1 = y / 2
        q2 = q1
        D = a * a - 4 * (b - y)
        if abs(D) < .ulpOfOne {
            p1 = a / 2
            p2 = p1
        } else {
            let sqrtD = sqrt(D)
            p1 = (a + sqrtD) / 2
            p2 = (a - sqrtD) / 2
        }
    } else {
        let sqrtD = sqrt(D)
        q1 = (y + sqrtD) / 2
        q2 = (y - sqrtD) / 2
        p1 = (a * q1 - c) / (q1 - q2)
        p2 = (c - a * q2) / (q1 - q2)
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

/// Solves a quartic polynomial with monic coefficients by calling the overloaded function that accepts the individual coefficients as parameters.
public func solveQuarticMonic(
    _ polynom: inout [Double]
) -> [Double] {
    solveQuarticMonic(polynom[0], polynom[1], polynom[2], polynom[3])
}

/// Evaluates a polynomial with coefficients in p at the point x.
///
/// - Parameters:
///   - p: Array of polynomial coefficients from highest to lowest degree.
///   - x: Point to evaluate the polynomial at.
/// - Returns: Value of the polynomial evaluated at x.
public func evaluatePolynomial(
    _ p: [Double],
    _ x: Double
) -> Double {

    var retVal: Double = 0.0

    if p.count == 0 {
        return retVal
    }

    if abs(x) < .ulpOfOne {
        retVal = p[p.count - 1]
    } else if x == 1.0 {
        for val in p.reversed() {
            retVal += val
        }
    } else {
        var xn: Double = 1.0

        for val in p.reversed() {
            retVal += val * xn
            xn *= x
        }
    }

    return retVal
}

/// Returns the derivative of a polynomial with the given coefficients.
///
/// - Parameters:
///   - coeffs: Array of polynomial coefficients from highest to lowest degree.
/// - Returns: Array containing the derivative polynomial's coefficients.
public func polynomialDerivative(_ coeffs: inout [Double]) -> [Double] {
    guard coeffs.count > 1 else { return [] }

    var deriv: [Double] = Array(repeating: 0.0, count: coeffs.count - 1)
    let count = coeffs.count - 1
    for i in 0..<count {
        deriv[i] = Double(count - i) * coeffs[i]
    }
    return deriv
}

/// Returns the derivative of a monic polynomial with the given coefficients.
///
/// - Parameters:
///   - monicCoeffs: Array of monic polynomial coefficients from highest to lowest degree.
/// - Returns: Array containing the derivative polynomial's coefficients.
public func polynomialMonicDerivative(_ monicCoeffs: inout [Double]) -> [Double] {
    guard monicCoeffs.count > 1 else { return [] }

    var deriv: [Double] = Array(repeating: 0.0, count: monicCoeffs.count - 1)
    deriv[0] = 1.0
    let N = monicCoeffs.count
    for i in 1..<(N - 1) {
        deriv[i] = Double(N - 1 - i) * monicCoeffs[i] / Double(N - 1)
    }
    return deriv
}

/// Shrinks the interval [li, hi] to contain a root of the polynomial p.
///
/// - Parameters:
///   - p: Array of polynomial coefficients
///   - li: Lower bound of initial interval
///   - hi: Upper bound of initial interval
/// - Returns: Approximation of a root contained in the shrunken interval
public func shrinkInterval(
    _ p: inout [Double],
    _ li: Double,
    _ hi: Double
) -> Double {
    var l = li
    var h = hi

    let fl = evaluatePolynomial(p, l)
    if fl == 0.0 { return l }

    let fh = evaluatePolynomial(p, h)
    if fh == 0.0 { return h }

    if fl > 0.0 { swap(&l, &h) }

    var rts = (l + h) / 2
    var dxold = abs(h - l)
    var dx = dxold
    let deriv = polynomialDerivative(&p)
    var f = evaluatePolynomial(p, rts)
    var df = evaluatePolynomial(deriv, rts)
    var temp: Double

    for _ in 0..<maxIts {
        if (((rts - h) * df - f) * ((rts - l) * df - f) > 0.0) || (abs(2 * f) > abs(dxold * df)) {
            dxold = dx
            dx = (h - l) / 2
            rts = l + dx
            if l == rts { break }

        } else {
            dxold = dx
            dx = f / df
            temp = rts
            rts -= dx
            if temp == rts { break }
        }

        if abs(dx) < tolerance { break }

        f = evaluatePolynomial(p, rts)
        df = evaluatePolynomial(deriv, rts)

        if f < 0.0 {
            l = rts
        } else {
            h = rts
        }
    }

    return rts
}

/// Custom cube root function that handles negative values correctly
private func cbrt(_ x: Double) -> Double {
    if x >= 0 {
        return pow(x, 1.0 / 3.0)
    } else {
        return -pow(-x, 1.0 / 3.0)
    }
}


private func solveWithTrigonometric(p: Double, q: Double, a: Double, b: Double) -> [Double] {
    // For three real roots case: discriminant > 0
    let m = 2.0 * sqrt(-p / 3.0)
    let theta = (1.0 / 3.0) * acos((3.0 * q) / (p * m))
    
    // Transform back from depressed cubic
    let shift = -b / (3.0 * a)
    
    return [
        m * cos(theta) + shift,
        m * cos(theta - 2.0 * .pi / 3.0) + shift,
        m * cos(theta - 4.0 * .pi / 3.0) + shift
    ]
}

private func solveWithCardano(p: Double, q: Double, a: Double, b: Double) -> [Double] {
    // For one real root case: discriminant <= 0
    let discriminant = (q * q / 4.0) + (p * p * p / 27.0)
    
    if discriminant >= 0 {
        let sqrtDisc = sqrt(discriminant)
        let u = cbrt(-q / 2.0 + sqrtDisc)
        let v = cbrt(-q / 2.0 - sqrtDisc)
        
        // Transform back from depressed cubic
        let shift = -b / (3.0 * a)
        return [u + v + shift]
    } else {
        // Fallback to trigonometric method for edge case
        return solveWithTrigonometric(p: p, q: q, a: a, b: b)
    }
}