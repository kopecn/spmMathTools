import Foundation


extension Date {
    func addingTimeInterval<T: BinaryFloatingPoint>(_ interval: T) -> Date {
        return Date(timeIntervalSinceReferenceDate: self.timeIntervalSinceReferenceDate + Double(interval))
    }
    func addingTimeInterval<T: BinaryInteger>(_ interval: T) -> Date {
        return Date(timeIntervalSinceReferenceDate: self.timeIntervalSinceReferenceDate + Double(interval))
    }
}