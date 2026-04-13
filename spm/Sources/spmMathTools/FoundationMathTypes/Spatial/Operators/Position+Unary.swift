import Foundation
import FoundationTypes
import simd

// MARK: - Unary Operators for Double Position

extension Position where T == Double {

    /// Negate a position (flip direction)
    @inlinable
    public static prefix func - (position: Position<Double>) -> Position<Double> {
        Position(vector: -position.vector)
    }

    /// Unary plus (returns copy)
    @inlinable
    public static prefix func + (position: Position<Double>) -> Position<Double> {
        position
    }
}

// MARK: - Unary Operators for Float Position

extension Position where T == Float {

    /// Negate a position (flip direction)
    @inlinable
    public static prefix func - (position: Position<Float>) -> Position<Float> {
        Position(vector: -position.vector)
    }

    /// Unary plus (returns copy)
    @inlinable
    public static prefix func + (position: Position<Float>) -> Position<Float> {
        position
    }
}
