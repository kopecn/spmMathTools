//
//  Result.swift
//
//
//  Created by Nicholas Bergantz on 4/25/24.
//

import Foundation

/// Result type of Ruckig's update function
public enum Result: Int {

    /// The trajectory is calculated normally
    case Working = 0
    /// The trajectory has reached its final position
    case Finished = 1
    /// Unclassified error
    case Error = -1
    /// Error in the input parameter
    case ErrorInvalidInput = -100
    /// The trajectory duration exceeds its numerical limits
    case ErrorTrajectoryDuration = -101
    /// The trajectory exceeds the given positional limits (only in Ruckig Pro)
    case ErrorPositionalLimits = -102

    /// The trajectory cannot be phase synchronized
    // ErrorNoPhaseSynchronization = -103

    /// The trajectory is not valid due to a conflict with zero limits
    case ErrorZeroLimits = -104
    /// Error during the extremel time calculation (Step 1)
    case ErrorExecutionTimeCalculation = -110
    /// Error during the synchronization calculation (Step 2)
    case ErrorSynchronizationCalculation = -111
}
