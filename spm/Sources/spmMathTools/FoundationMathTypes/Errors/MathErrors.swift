//
//  File.swift
//
//
//  Created by Nicholas Bergantz on 4/21/24.
//

import Foundation

public struct MathError: Error, CustomStringConvertible {
    public let message: String
    public var description: String { message }
    public init(_ message: String) { self.message = message }
}

public enum MathErrors: Error {
    /// General Runtime Errors
    case runtimeError(_ msg: String)

    /// Size of polynomial is invalid
    case invalidPolynomialSize(_ poly: PolynomialUnivariateOrder, _ gotSize: Int)
    /// the roots for this polynomial are unsolvable.
    case rootsUnsolveableFor(_ poly: PolynomialUnivariateOrder, _ withCoeffs: [Double])
}

extension MathErrors: CustomStringConvertible {
    public var description: String {
        switch self {
        case .invalidPolynomialSize(let poly, let size):
            return "The provided polynomial size \(size) was invalid for \(poly)."

        case .rootsUnsolveableFor(let poly, let coeffs):
            return "Polynomial \(poly) has unsolvable roots for coeffs: \(coeffs)."

        case .runtimeError(let msg):
            return "OTG Runtime Error, \(msg)"
        }
    }
}
