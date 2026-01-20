import Foundation

extension InputParameter: Codable {

    enum CodingKeys: String, CodingKey {
        case degreesOfFreedom
        case controlInterface
        case synchronization
        case durationDiscretization
        case currentPosition
        case currentVelocity
        case currentAcceleration
        case targetPosition
        case targetVelocity
        case targetAcceleration
        case maxVelocity
        case maxAcceleration
        case maxJerk
        case minVelocity
        case minAcceleration
        case intermediatePositions
        case perSectionMaxVelocity
        case perSectionMaxAcceleration
        case perSectionMaxJerk
        case perSectionMinVelocity
        case perSectionMinAcceleration
        case perSectionMaxPosition
        case perSectionMinPosition
        case maxPosition
        case minPosition
        case enabled
        case perDofControlInterface
        case perDofSynchronization
        case minimumDuration
        case perSectionMinimumDuration
        case interruptCalculationDuration
    }

    public init(from decoder: Decoder) throws {

        let container = try decoder.container(keyedBy: CodingKeys.self)

        let degreesOfFreedom = try container.decode(Int.self, forKey: .degreesOfFreedom)

        self.degreesOfFreedom = degreesOfFreedom

        self.controlInterface = try container.decode(ControlInterface.self, forKey: .controlInterface)
        self.synchronization = try container.decode(Synchronization.self, forKey: .synchronization)
        self.durationDiscretization = try container.decode(
            DurationDiscretization.self,
            forKey: .durationDiscretization
        )

        self.currentPosition = try container.decode([Double].self, forKey: .currentPosition)
        self.currentVelocity = try container.decode([Double].self, forKey: .currentVelocity)
        self.currentAcceleration = try container.decode([Double].self, forKey: .currentAcceleration)

        self.targetPosition = try container.decode([Double].self, forKey: .targetPosition)
        self.targetVelocity = try container.decode([Double].self, forKey: .targetVelocity)
        self.targetAcceleration = try container.decode([Double].self, forKey: .targetAcceleration)

        self.intermediatePositions = try container.decode([[Double]].self, forKey: .intermediatePositions)

        self.enabled = try container.decode([Bool].self, forKey: .enabled)

        self.maxVelocity = try container.decode([Double].self, forKey: .maxVelocity)
        self.maxAcceleration = try container.decode([Double].self, forKey: .maxAcceleration)
        self.maxJerk = try container.decode([Double].self, forKey: .maxJerk)

        // Decode Optionals
        self.minVelocity = try container.decodeIfPresent([Double].self, forKey: .minVelocity)
        self.minAcceleration = try container.decodeIfPresent([Double].self, forKey: .minAcceleration)

        self.perSectionMaxVelocity = try container.decodeIfPresent(
            [[Double]].self,
            forKey: .perSectionMaxVelocity
        )
        self.perSectionMaxAcceleration = try container.decodeIfPresent(
            [[Double]].self,
            forKey: .perSectionMaxAcceleration
        )
        self.perSectionMaxJerk = try container.decodeIfPresent([[Double]].self, forKey: .perSectionMaxJerk)
        self.perSectionMinVelocity = try container.decodeIfPresent(
            [[Double]].self,
            forKey: .perSectionMinVelocity
        )
        self.perSectionMinAcceleration = try container.decodeIfPresent(
            [[Double]].self,
            forKey: .perSectionMinAcceleration
        )
        self.perSectionMaxPosition = try container.decodeIfPresent(
            [[Double]].self,
            forKey: .perSectionMaxPosition
        )
        self.perSectionMinPosition = try container.decodeIfPresent(
            [[Double]].self,
            forKey: .perSectionMinPosition
        )

        self.maxPosition = try container.decodeIfPresent([Double].self, forKey: .maxPosition)
        self.minPosition = try container.decodeIfPresent([Double].self, forKey: .minPosition)

        self.perDofControlInterface = try container.decodeIfPresent(
            [ControlInterface].self,
            forKey: .perDofControlInterface
        )
        self.perDofSynchronization = try container.decodeIfPresent(
            [Synchronization].self,
            forKey: .perDofSynchronization
        )

        self.minimumDuration = try container.decodeIfPresent(Double.self, forKey: .minimumDuration)
        self.perSectionMinimumDuration = try container.decodeIfPresent(
            [Double].self,
            forKey: .perSectionMinimumDuration
        )

        self.interruptCalculationDuration = try container.decodeIfPresent(
            Double.self,
            forKey: .interruptCalculationDuration
        )

    }

    public func encode(to encoder: Encoder) throws {
        var container = encoder.container(keyedBy: CodingKeys.self)

        try container.encode(degreesOfFreedom, forKey: .degreesOfFreedom)
        try container.encode(controlInterface, forKey: .controlInterface)
        try container.encode(synchronization, forKey: .synchronization)
        try container.encode(durationDiscretization, forKey: .durationDiscretization)

        try container.encode(currentPosition, forKey: .currentPosition)
        try container.encode(currentVelocity, forKey: .currentVelocity)
        try container.encode(currentAcceleration, forKey: .currentAcceleration)

        try container.encode(targetPosition, forKey: .targetPosition)
        try container.encode(targetVelocity, forKey: .targetVelocity)
        try container.encode(targetAcceleration, forKey: .targetAcceleration)

        try container.encode(intermediatePositions, forKey: .intermediatePositions)

        try container.encode(enabled, forKey: .enabled)

        try container.encode(maxVelocity, forKey: .maxVelocity)
        try container.encode(maxAcceleration, forKey: .maxAcceleration)
        try container.encode(maxJerk, forKey: .maxJerk)

        // Encode Optionals
        try container.encodeIfPresent(minVelocity, forKey: .minVelocity)
        try container.encodeIfPresent(minAcceleration, forKey: .minAcceleration)

        try container.encodeIfPresent(perSectionMaxVelocity, forKey: .perSectionMaxVelocity)
        try container.encodeIfPresent(perSectionMaxAcceleration, forKey: .perSectionMaxAcceleration)
        try container.encodeIfPresent(perSectionMaxJerk, forKey: .perSectionMaxJerk)
        try container.encodeIfPresent(perSectionMinVelocity, forKey: .perSectionMinVelocity)
        try container.encodeIfPresent(perSectionMinAcceleration, forKey: .perSectionMinAcceleration)
        try container.encodeIfPresent(perSectionMaxPosition, forKey: .perSectionMaxPosition)
        try container.encodeIfPresent(perSectionMinPosition, forKey: .perSectionMinPosition)

        try container.encodeIfPresent(maxPosition, forKey: .maxPosition)
        try container.encodeIfPresent(minPosition, forKey: .minPosition)

        try container.encodeIfPresent(perDofControlInterface, forKey: .perDofControlInterface)
        try container.encodeIfPresent(perDofSynchronization, forKey: .perDofSynchronization)

        try container.encodeIfPresent(minimumDuration, forKey: .minimumDuration)
        try container.encodeIfPresent(perSectionMinimumDuration, forKey: .perSectionMinimumDuration)

        try container.encodeIfPresent(interruptCalculationDuration, forKey: .interruptCalculationDuration)
    }
}
