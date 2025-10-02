//
//  Utils.swift
//
//
//  Created by Nicholas Bergantz on 3/23/24.
//
import Foundation

let cos120 = -0.50
let sin120 = 0.866025403784438646764

extension Array where Element == Double {

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
