import Foundation
import FoundationTypes
import simd

// MARK: - Arithmetic Operators for Double Position

extension Position where T == Double {

    // MARK: Addition

    /// Add two positions (vector addition)
    @inlinable
    public static func + (lhs: Position<Double>, rhs: Position<Double>) -> Position<Double> {
        Position(vector: lhs.vector + rhs.vector)
    }

    /// Add a scalar to all components
    @inlinable
    public static func + (lhs: Position<Double>, rhs: Double) -> Position<Double> {
        Position(vector: lhs.vector + SIMD3<Double>(repeating: rhs))
    }

    /// Add a position to a scalar
    @inlinable
    public static func + (lhs: Double, rhs: Position<Double>) -> Position<Double> {
        rhs + lhs
    }

    // MARK: Subtraction

    /// Subtract two positions (vector subtraction)
    @inlinable
    public static func - (lhs: Position<Double>, rhs: Position<Double>) -> Position<Double> {
        Position(vector: lhs.vector - rhs.vector)
    }

    /// Subtract a scalar from all components
    @inlinable
    public static func - (lhs: Position<Double>, rhs: Double) -> Position<Double> {
        Position(vector: lhs.vector - SIMD3<Double>(repeating: rhs))
    }

    /// Subtract a position from a scalar
    @inlinable
    public static func - (lhs: Double, rhs: Position<Double>) -> Position<Double> {
        Position(vector: SIMD3<Double>(repeating: lhs) - rhs.vector)
    }

    // MARK: Multiplication (Scaling)

    /// Multiply position by a scalar (scaling)
    @inlinable
    public static func * (lhs: Position<Double>, rhs: Double) -> Position<Double> {
        Position(vector: lhs.vector * rhs)
    }

    /// Multiply scalar by a position (scaling)
    @inlinable
    public static func * (lhs: Double, rhs: Position<Double>) -> Position<Double> {
        Position(vector: lhs * rhs.vector)
    }

    // MARK: Division (Scaling)

    /// Divide position by a scalar
    @inlinable
    public static func / (lhs: Position<Double>, rhs: Double) -> Position<Double> {
        Position(vector: lhs.vector / rhs)
    }
}

// MARK: - Arithmetic Operators for Float Position

extension Position where T == Float {

    // MARK: Addition

    /// Add two positions (vector addition)
    @inlinable
    public static func + (lhs: Position<Float>, rhs: Position<Float>) -> Position<Float> {
        Position(vector: lhs.vector + rhs.vector)
    }

    /// Add a scalar to all components
    @inlinable
    public static func + (lhs: Position<Float>, rhs: Float) -> Position<Float> {
        Position(vector: lhs.vector + SIMD3<Float>(repeating: rhs))
    }

    /// Add a position to a scalar
    @inlinable
    public static func + (lhs: Float, rhs: Position<Float>) -> Position<Float> {
        rhs + lhs
    }

    // MARK: Subtraction

    /// Subtract two positions (vector subtraction)
    @inlinable
    public static func - (lhs: Position<Float>, rhs: Position<Float>) -> Position<Float> {
        Position(vector: lhs.vector - rhs.vector)
    }

    /// Subtract a scalar from all components
    @inlinable
    public static func - (lhs: Position<Float>, rhs: Float) -> Position<Float> {
        Position(vector: lhs.vector - SIMD3<Float>(repeating: rhs))
    }

    /// Subtract a position from a scalar
    @inlinable
    public static func - (lhs: Float, rhs: Position<Float>) -> Position<Float> {
        Position(vector: SIMD3<Float>(repeating: lhs) - rhs.vector)
    }

    // MARK: Multiplication (Scaling)

    /// Multiply position by a scalar (scaling)
    @inlinable
    public static func * (lhs: Position<Float>, rhs: Float) -> Position<Float> {
        Position(vector: lhs.vector * rhs)
    }

    /// Multiply scalar by a position (scaling)
    @inlinable
    public static func * (lhs: Float, rhs: Position<Float>) -> Position<Float> {
        Position(vector: lhs * rhs.vector)
    }

    // MARK: Division (Scaling)

    /// Divide position by a scalar
    @inlinable
    public static func / (lhs: Position<Float>, rhs: Float) -> Position<Float> {
        Position(vector: lhs.vector / rhs)
    }
}
