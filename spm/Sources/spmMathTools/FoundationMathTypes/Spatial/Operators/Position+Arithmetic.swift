import Foundation
import FoundationTypes
import simd

// MARK: - Arithmetic Operators for Double Position

extension Position where T == Double {

    // MARK: Addition

    /// Add two positions (vector addition)
    public static func + (lhs: Position<Double>, rhs: Position<Double>) -> Position<Double> {
        return Position(vector: lhs.vector + rhs.vector)
    }

    /// Add a scalar to all components
    public static func + (lhs: Position<Double>, rhs: Double) -> Position<Double> {
        return Position(vector: lhs.vector + SIMD3<Double>(repeating: rhs))
    }

    /// Add a position to a scalar
    public static func + (lhs: Double, rhs: Position<Double>) -> Position<Double> {
        return rhs + lhs
    }

    // MARK: Subtraction

    /// Subtract two positions (vector subtraction)
    public static func - (lhs: Position<Double>, rhs: Position<Double>) -> Position<Double> {
        return Position(vector: lhs.vector - rhs.vector)
    }

    /// Subtract a scalar from all components
    public static func - (lhs: Position<Double>, rhs: Double) -> Position<Double> {
        return Position(vector: lhs.vector - SIMD3<Double>(repeating: rhs))
    }

    /// Subtract a position from a scalar
    public static func - (lhs: Double, rhs: Position<Double>) -> Position<Double> {
        return Position(vector: SIMD3<Double>(repeating: lhs) - rhs.vector)
    }

    // MARK: Multiplication (Scaling)

    /// Multiply position by a scalar (scaling)
    public static func * (lhs: Position<Double>, rhs: Double) -> Position<Double> {
        return Position(vector: lhs.vector * rhs)
    }

    /// Multiply scalar by a position (scaling)
    public static func * (lhs: Double, rhs: Position<Double>) -> Position<Double> {
        return Position(vector: lhs * rhs.vector)
    }

    // MARK: Division (Scaling)

    /// Divide position by a scalar
    public static func / (lhs: Position<Double>, rhs: Double) -> Position<Double> {
        return Position(vector: lhs.vector / rhs)
    }
}


// MARK: - Arithmetic Operators for Float Position

extension Position where T == Float {

    // MARK: Addition

    /// Add two positions (vector addition)
    public static func + (lhs: Position<Float>, rhs: Position<Float>) -> Position<Float> {
        return Position(vector: lhs.vector + rhs.vector)
    }

    /// Add a scalar to all components
    public static func + (lhs: Position<Float>, rhs: Float) -> Position<Float> {
        return Position(vector: lhs.vector + SIMD3<Float>(repeating: rhs))
    }

    /// Add a position to a scalar
    public static func + (lhs: Float, rhs: Position<Float>) -> Position<Float> {
        return rhs + lhs
    }

    // MARK: Subtraction

    /// Subtract two positions (vector subtraction)
    public static func - (lhs: Position<Float>, rhs: Position<Float>) -> Position<Float> {
        return Position(vector: lhs.vector - rhs.vector)
    }

    /// Subtract a scalar from all components
    public static func - (lhs: Position<Float>, rhs: Float) -> Position<Float> {
        return Position(vector: lhs.vector - SIMD3<Float>(repeating: rhs))
    }

    /// Subtract a position from a scalar
    public static func - (lhs: Float, rhs: Position<Float>) -> Position<Float> {
        return Position(vector: SIMD3<Float>(repeating: lhs) - rhs.vector)
    }

    // MARK: Multiplication (Scaling)

    /// Multiply position by a scalar (scaling)
    public static func * (lhs: Position<Float>, rhs: Float) -> Position<Float> {
        return Position(vector: lhs.vector * rhs)
    }

    /// Multiply scalar by a position (scaling)
    public static func * (lhs: Float, rhs: Position<Float>) -> Position<Float> {
        return Position(vector: lhs * rhs.vector)
    }

    // MARK: Division (Scaling)

    /// Divide position by a scalar
    public static func / (lhs: Position<Float>, rhs: Float) -> Position<Float> {
        return Position(vector: lhs.vector / rhs)
    }
}
