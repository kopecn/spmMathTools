//
//  File.swift
//
//
//  Created by Nicholas Bergantz on 4/21/24.
//

import Foundation

public struct RuckigError: Error, CustomStringConvertible, LocalizedError {
    public let message: String
    public var description: String { message }
    public var errorDescription: String? { message }
    public init(_ message: String) { self.message = message }
}

enum OTGErrors: Error, LocalizedError {
    /// General Runtime Errors
    case runtimeError(_ msg: String)
}

extension OTGErrors: CustomStringConvertible {
    public var description: String {
        switch self {
        case .runtimeError(let msg):
            return "OTG Runtime Error, \(msg)"
        }
    }

    public var errorDescription: String? { description }
}
