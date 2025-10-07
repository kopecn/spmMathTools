//
//  Utils.swift
//
//
//  Created by Nicholas Bergantz on 3/23/24.
//
import Foundation

// MARK: - Trigonometric Constants

/// Cosine of 120 degrees (-0.5)
let cos120 = -0.50
/// Sine of 120 degrees (√3/2)
let sin120 = 0.866025403784438646764

// MARK: - Numerical Tolerance Constants

/// Machine epsilon scaled by 16 for numerical comparisons
public let eps16: Double = 16 * .ulpOfOne

/// Standard tolerance for polynomial root finding and general numerical operations
public let polynomialTolerance: Double = 1e-14

/// Threshold for values considered effectively zero in polynomial operations
/// Values within ±polynomialZeroThreshold are rounded to zero
public let polynomialZeroThreshold: Double = 1e-9

extension Array where Element == Double {

    /// Gate allowing only positive values to be appended.
    ///
    /// - Parameter val: Double value to potentially append. Only appended if >= 0
    mutating func insertIfPositive(_ val: Double) {
        if val >= 0 {
            self.append(val)
        }
    }

    /// Joins the array into a string such that it can be displayed
    ///
    /// Joins the elements of the array into a string, with each element separated by a comma and space.
    ///
    /// - Parameters:
    ///        - highPrecision: boolean, If true, elements will be joined with higher precision formatting.
    ///- Returns: A string containing the joined elements of the array.
    func join(_ highPrecision: Bool = false) -> String {
        let precise: String = highPrecision ? "%.016f" : "%.06f"
        return self.reduce("") { partialResult, next in
            partialResult + ", " + String(format: precise, next)
        }
    }
}

/// performs an integration on duration t, p0,v0,a0 for a given constant jerk
///
/// - Parameters:
///        - t: Double, time in seconds
///        - p0: Double, current position
///        - v0: Double, current velocity
///        - a0: Double, current acceleration
///        - j: Double, applied jerk for region
///- Returns: Tuple(postion, velocity and accleration)
func integrate(
    _ t: Double,
    _ p0: Double,
    _ v0: Double,
    _ a0: Double,
    _ j: Double
) -> (Double, Double, Double) {
    (
        p0 + t * (v0 + t * (a0 / 2 + t * j / 6)),
        v0 + t * (a0 + t * j / 2),
        a0 + t * j
    )
}
