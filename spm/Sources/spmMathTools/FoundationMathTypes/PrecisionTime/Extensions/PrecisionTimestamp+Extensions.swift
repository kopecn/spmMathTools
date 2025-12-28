import Foundation
import FoundationTypes

// MARK: - Comparable
extension PrecisionTimestamp: @retroactive Comparable {
    public static func < (lhs: PrecisionTimestamp, rhs: PrecisionTimestamp) -> Bool {
        // Compare based on the interval from epoch
        return lhs.interval < rhs.interval
    }
}