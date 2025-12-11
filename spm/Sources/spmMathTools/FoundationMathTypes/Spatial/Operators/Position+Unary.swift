import Foundation
import FoundationTypes
import simd

// MARK: - Unary Operators for Double Position

extension Position where T == Double {

    /// Negate a position (flip direction)
    public static prefix func - (position: Position<Double>) -> Position<Double> {
        return Position(vector: -position.vector)
    }

    /// Unary plus (returns copy)
    public static prefix func + (position: Position<Double>) -> Position<Double> {
        return position
    }
}

// MARK: - Unary Operators for Float Position

extension Position where T == Float {

    /// Negate a position (flip direction)
    public static prefix func - (position: Position<Float>) -> Position<Float> {
        return Position(vector: -position.vector)
    }

    /// Unary plus (returns copy)
    public static prefix func + (position: Position<Float>) -> Position<Float> {
        return position
    }
}
