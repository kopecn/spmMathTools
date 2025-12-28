import Foundation
import FoundationTypes

// MARK: - Comparable
extension PrecisionTimeInterval: @retroactive Comparable {
    public static func < (lhs: PrecisionTimeInterval, rhs: PrecisionTimeInterval) -> Bool {
        // Handle different signs
        if lhs.sign != rhs.sign {
            // Negative is less than positive
            return lhs.sign == .negative && rhs.sign == .positive
        }

        // Same sign - compare magnitudes
        if lhs.seconds != rhs.seconds {
            let secondsLess = lhs.seconds < rhs.seconds
            // If positive, normal comparison; if negative, reverse
            return lhs.sign == .positive ? secondsLess : !secondsLess
        }

        // Seconds are equal, compare attoseconds
        let attosecondsLess = lhs.attoseconds < rhs.attoseconds
        // If positive, normal comparison; if negative, reverse
        return lhs.sign == .positive ? attosecondsLess : !attosecondsLess
    }
}
