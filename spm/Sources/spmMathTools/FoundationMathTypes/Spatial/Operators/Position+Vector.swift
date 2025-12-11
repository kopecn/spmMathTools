import Foundation
import FoundationTypes
import simd

// MARK: - Vector Operators for Double Position

extension Position where T == Double {

    // MARK: Dot Product

    /// Compute dot product of two positions
    public func dot(_ other: Position<Double>) -> Double {
        return simd_dot(self.vector, other.vector)
    }

    /// Compute dot product using operator
    /// macOS: Option + 8, windows: Alt + 7
    public static func • (lhs: Position<Double>, rhs: Position<Double>) -> Double {
        return simd_dot(lhs.vector, rhs.vector)
    }

    // MARK: Cross Product

    /// Compute cross product of two positions
    public func cross(_ other: Position<Double>) -> Position<Double> {
        return Position(vector: simd_cross(self.vector, other.vector))
    }

    /// Compute cross product using operator
    /// macOS: Option + 00D7, windows: Alt + 0215
    public static func × (lhs: Position<Double>, rhs: Position<Double>) -> Position<Double> {
        return Position(vector: simd_cross(lhs.vector, rhs.vector))
    }

    // MARK: Distance

    /// Compute distance to another position
    public func distance(to other: Position<Double>) -> Double {
        return simd_distance(self.vector, other.vector)
    }

    /// Compute squared distance to another position (more efficient than distance)
    public func distanceSquared(to other: Position<Double>) -> Double {
        return simd_distance_squared(self.vector, other.vector)
    }
}

// MARK: - Vector Operators for Float Position

extension Position where T == Float {

    // MARK: Dot Product

    /// Compute dot product of two positions
    public func dot(_ other: Position<Float>) -> Float {
        return simd_dot(self.vector, other.vector)
    }

    /// Compute dot product using operator
    public static func • (lhs: Position<Float>, rhs: Position<Float>) -> Float {
        return simd_dot(lhs.vector, rhs.vector)
    }

    // MARK: Cross Product

    /// Compute cross product of two positions
    public func cross(_ other: Position<Float>) -> Position<Float> {
        return Position(vector: simd_cross(self.vector, other.vector))
    }

    /// Compute cross product using operator
    public static func × (lhs: Position<Float>, rhs: Position<Float>) -> Position<Float> {
        return Position(vector: simd_cross(lhs.vector, rhs.vector))
    }

    // MARK: Distance

    /// Compute distance to another position
    public func distance(to other: Position<Float>) -> Float {
        return simd_distance(self.vector, other.vector)
    }

    /// Compute squared distance to another position (more efficient than distance)
    public func distanceSquared(to other: Position<Float>) -> Float {
        return simd_distance_squared(self.vector, other.vector)
    }
}

// MARK: - Operator Declarations

infix operator •: MultiplicationPrecedence  // Dot product
infix operator ×: MultiplicationPrecedence  // Cross product
