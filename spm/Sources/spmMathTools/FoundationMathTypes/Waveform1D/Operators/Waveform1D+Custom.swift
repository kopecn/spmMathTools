import Foundation
import FoundationTypes

// MARK: - Custom Operators for All Numeric Types

extension Waveform1D {

    // MARK: Custom Prefix Operators

    // TODO: Add custom prefix operators here
    // Example structure:
    // public static prefix func ~~(waveform: Waveform1D<T>) -> Waveform1D<T> {
    //     // Implementation
    // }

    // MARK: Custom Infix Operators

    // TODO: Add custom infix operators here
    // Example structure:
    // public static func <*>(lhs: Waveform1D<T>, rhs: Waveform1D<T>) -> Waveform1D<T> {
    //     // Implementation
    // }

    // MARK: Custom Postfix Operators

    // TODO: Add custom postfix operators here
    // Example structure:
    // public static postfix func !!(waveform: Waveform1D<T>) -> Waveform1D<T> {
    //     // Implementation
    // }
}

// MARK: - Custom Operators for Comparable Types

extension Waveform1D where T: Comparable {

    // MARK: Custom Prefix Operators

    // TODO: Add custom prefix operators for Comparable types

    // MARK: Custom Infix Operators

    // TODO: Add custom infix operators for Comparable types

    // MARK: Custom Postfix Operators

    // TODO: Add custom postfix operators for Comparable types
}

// MARK: - Custom Operators for Floating Point Types

extension Waveform1D where T: BinaryFloatingPoint {

    // MARK: Custom Prefix Operators

    // TODO: Add custom prefix operators for BinaryFloatingPoint types

    // MARK: Custom Infix Operators

    // TODO: Add custom infix operators for BinaryFloatingPoint types

    // MARK: Custom Postfix Operators

    // TODO: Add custom postfix operators for BinaryFloatingPoint types
}

// MARK: - Custom Operators for Integer Types

extension Waveform1D where T: BinaryInteger {

    // MARK: Custom Prefix Operators

    // TODO: Add custom prefix operators for BinaryInteger types

    // MARK: Custom Infix Operators

    // TODO: Add custom infix operators for BinaryInteger types

    // MARK: Custom Postfix Operators

    // TODO: Add custom postfix operators for BinaryInteger types
}

// MARK: - Custom Operators for Signed Integer Types

extension Waveform1D where T: SignedInteger {

    // MARK: Custom Prefix Operators

    // TODO: Add custom prefix operators for SignedInteger types

    // MARK: Custom Infix Operators

    // TODO: Add custom infix operators for SignedInteger types

    // MARK: Custom Postfix Operators

    // TODO: Add custom postfix operators for SignedInteger types
}

// MARK: - Operator Declarations

// NOTE: Custom operators must be declared at the global scope
// Declare your custom operators here before implementing them above

// Example operator declarations:
// prefix operator ~~
// infix operator <*>: MultiplicationPrecedence
// postfix operator !!

// Available precedence groups:
// - AssignmentPrecedence
// - TernaryPrecedence (ternary conditional operator)
// - DefaultPrecedence
// - LogicalDisjunctionPrecedence (||)
// - LogicalConjunctionPrecedence (&&)
// - ComparisonPrecedence (<, <=, >, >=, ==, !=, ===, !==, ~=)
// - NilCoalescencePrecedence (??)
// - CastingPrecedence (is, as, as?, as!)
// - RangeFormationPrecedence (..., ..<)
// - AdditionPrecedence (+, -, |, ^, etc.)
// - MultiplicationPrecedence (*, /, %, &, etc.)
// - BitwiseShiftPrecedence (<<, >>, etc.)
