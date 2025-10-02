import Foundation

// MARK: - Digital Filtering
extension Waveform1D where T: BinaryFloatingPoint {

    /// Apply a digital filter to the waveform
    /// - Parameters:
    ///   - filterType: Type of filter to apply
    ///   - order: Filter order (higher order = steeper rolloff)
    /// - Returns: Filtered waveform
    public func filtered(with filterType: WaveformFilterType, order: Int = 4) -> Waveform1D<T> {
        guard !values.isEmpty && order > 0 else { return self }

        switch filterType {
        case .lowPass(let cutoffFrequency):
            return applyButterworthFilter(type: .lowPass, cutoffFrequency: cutoffFrequency, order: order)
        case .highPass(let cutoffFrequency):
            return applyButterworthFilter(type: .highPass, cutoffFrequency: cutoffFrequency, order: order)
        case .bandPass(let lowFrequency, let highFrequency):
            return applyBandPassFilter(lowFrequency: lowFrequency, highFrequency: highFrequency, order: order)
        case .bandStop(let lowFrequency, let highFrequency):
            return applyBandStopFilter(lowFrequency: lowFrequency, highFrequency: highFrequency, order: order)
        case .movingAverage(let windowSize):
            return applyMovingAverage(windowSize: windowSize)
        case .exponential(let alpha):
            return applyExponentialFilter(alpha: alpha)
        }
    }

    /// Apply a simple moving average filter
    /// - Parameter windowSize: Size of the moving average window
    /// - Returns: Filtered waveform
    public func movingAverageFilter(windowSize: Int) -> Waveform1D<T> {
        return filtered(with: .movingAverage(windowSize: windowSize))
    }

    /// Apply an exponential smoothing filter
    /// - Parameter alpha: Smoothing factor (0 < alpha <= 1, higher = less smoothing)
    /// - Returns: Filtered waveform
    public func exponentialFilter(alpha: Double) -> Waveform1D<T> {
        return filtered(with: .exponential(alpha: alpha))
    }

    /// Apply a Butterworth low-pass filter
    /// - Parameters:
    ///   - cutoffFrequency: Cutoff frequency in Hz
    ///   - order: Filter order
    /// - Returns: Filtered waveform
    public func lowPassFilter(cutoffFrequency: Double, order: Int = 4) -> Waveform1D<T> {
        return filtered(with: .lowPass(cutoffFrequency: cutoffFrequency), order: order)
    }

    /// Apply a Butterworth high-pass filter
    /// - Parameters:
    ///   - cutoffFrequency: Cutoff frequency in Hz
    ///   - order: Filter order
    /// - Returns: Filtered waveform
    public func highPassFilter(cutoffFrequency: Double, order: Int = 4) -> Waveform1D<T> {
        return filtered(with: .highPass(cutoffFrequency: cutoffFrequency), order: order)
    }

    /// Apply a band-pass filter
    /// - Parameters:
    ///   - lowFrequency: Lower cutoff frequency in Hz
    ///   - highFrequency: Upper cutoff frequency in Hz
    ///   - order: Filter order
    /// - Returns: Filtered waveform
    public func bandPassFilter(lowFrequency: Double, highFrequency: Double, order: Int = 4) -> Waveform1D<T> {
        return filtered(with: .bandPass(lowFrequency: lowFrequency, highFrequency: highFrequency), order: order)
    }

    // MARK: - Private Implementation Methods

    private func applyButterworthFilter(
        type: WaveformFilterButterworthType,
        cutoffFrequency: Double,
        order: Int
    ) -> Waveform1D<T> {
        let nyquist = samplingFrequency / 2.0
        let normalizedCutoff = cutoffFrequency / nyquist

        // Clamp normalized frequency to valid range
        let clampedCutoff = max(0.001, min(0.999, normalizedCutoff))

        // Calculate filter coefficients
        let coefficients = calculateButterworthCoefficients(type: type, normalizedCutoff: clampedCutoff, order: order)

        // Apply IIR filter
        return applyIIRFilter(coefficients: coefficients)
    }

    private func applyBandPassFilter(lowFrequency: Double, highFrequency: Double, order: Int) -> Waveform1D<T> {
        // Apply high-pass first, then low-pass
        let highPassed = highPassFilter(cutoffFrequency: lowFrequency, order: order)
        return highPassed.lowPassFilter(cutoffFrequency: highFrequency, order: order)
    }

    private func applyBandStopFilter(lowFrequency: Double, highFrequency: Double, order: Int) -> Waveform1D<T> {
        // Create band-pass filter and subtract from original
        let bandPass = bandPassFilter(lowFrequency: lowFrequency, highFrequency: highFrequency, order: order)

        let filtered = zip(values, bandPass.values).map { original, bandPassValue in
            original - bandPassValue
        }

        return Waveform1D(values: filtered, dt: dt, t0: t0)
    }

    private func applyMovingAverage(windowSize: Int) -> Waveform1D<T> {
        guard windowSize > 0 && windowSize <= values.count else { return self }

        var filtered: [T] = []
        let halfWindow = windowSize / 2

        for i in 0..<values.count {
            let start = max(0, i - halfWindow)
            let end = min(values.count, i + halfWindow + 1)

            let windowSum = values[start..<end].reduce(T.zero, +)
            let average = windowSum / T(end - start)
            filtered.append(average)
        }

        return Waveform1D(values: filtered, dt: dt, t0: t0)
    }

    private func applyExponentialFilter(alpha: Double) -> Waveform1D<T> {
        guard !values.isEmpty && alpha > 0 && alpha <= 1 else { return self }

        var filtered: [T] = []
        let alphaT = T(alpha)
        let oneMinusAlpha = T(1.0 - alpha)

        // Initialize with first value
        var smoothed = values[0]
        filtered.append(smoothed)

        // Apply exponential smoothing
        for i in 1..<values.count {
            smoothed = alphaT * values[i] + oneMinusAlpha * smoothed
            filtered.append(smoothed)
        }

        return Waveform1D(values: filtered, dt: dt, t0: t0)
    }

    private func calculateButterworthCoefficients(
        type: WaveformFilterButterworthType,
        normalizedCutoff: Double,
        order: Int
    ) -> WaveformFilterCoefficients<T> {
        // Simplified Butterworth coefficient calculation
        // For a more complete implementation, you'd use proper pole-zero placement

        let wc = tan(.pi * normalizedCutoff / 2.0)
        let wc2 = wc * wc

        switch type {
        case .lowPass:
            return calculateLowPassCoefficients(wc: wc, wc2: wc2, order: order)
        case .highPass:
            return calculateHighPassCoefficients(wc: wc, wc2: wc2, order: order)
        }
    }

    private func calculateLowPassCoefficients(wc: Double, wc2: Double, order: Int) -> WaveformFilterCoefficients<T> {
        // Simplified 2nd order Butterworth low-pass
        let k1 = sqrt(2.0) * wc
        let k2 = wc2
        let denominator = k2 + k1 + 1.0

        let b = [
            T(k2 / denominator),
            T(2.0 * k2 / denominator),
            T(k2 / denominator),
        ]

        let a = [
            T(1.0),
            T(2.0 * (k2 - 1.0) / denominator),
            T((k2 - k1 + 1.0) / denominator),
        ]

        return WaveformFilterCoefficients(b: b, a: a)
    }

    private func calculateHighPassCoefficients(wc: Double, wc2: Double, order: Int) -> WaveformFilterCoefficients<T> {
        // Simplified 2nd order Butterworth high-pass
        let k1 = sqrt(2.0) * wc
        let k2 = wc2
        let denominator = k2 + k1 + 1.0

        let b = [
            T(1.0 / denominator),
            T(-2.0 / denominator),
            T(1.0 / denominator),
        ]

        let a = [
            T(1.0),
            T(2.0 * (k2 - 1.0) / denominator),
            T((k2 - k1 + 1.0) / denominator),
        ]

        return WaveformFilterCoefficients(b: b, a: a)
    }

    private func applyIIRFilter(coefficients: WaveformFilterCoefficients<T>) -> Waveform1D<T> {
        guard !values.isEmpty else { return self }

        let b = coefficients.b
        let a = coefficients.a
        let orderB = b.count
        let orderA = a.count

        var filtered: [T] = []
        var inputHistory: [T] = Array(repeating: T.zero, count: orderB)
        var outputHistory: [T] = Array(repeating: T.zero, count: orderA)

        for input in values {
            // Shift input history
            for i in stride(from: orderB - 1, to: 0, by: -1) {
                inputHistory[i] = inputHistory[i - 1]
            }
            inputHistory[0] = input

            // Calculate output
            var output = T.zero

            // Feedforward (numerator)
            for i in 0..<orderB {
                output += b[i] * inputHistory[i]
            }

            // Feedback (denominator, excluding a[0])
            for i in 1..<orderA {
                output -= a[i] * outputHistory[i]
            }

            // Shift output history
            for i in stride(from: orderA - 1, to: 0, by: -1) {
                outputHistory[i] = outputHistory[i - 1]
            }
            outputHistory[0] = output

            filtered.append(output)
        }

        return Waveform1D(values: filtered, dt: dt, t0: t0)
    }
}

// MARK: - Frequency Response Analysis
extension Waveform1D where T: BinaryFloatingPoint {

    /// Calculate the frequency response of the signal
    /// - Parameters:
    ///   - frequencies: Array of frequencies to evaluate
    ///   - filterType: Filter type to analyze
    ///   - order: Filter order
    /// - Returns: Tuple containing magnitude and phase response
    public func frequencyResponse(
        at frequencies: [Double],
        filterType: WaveformFilterType,
        order: Int = 4
    ) -> (magnitude: [T], phase: [T]) {

        let nyquist = samplingFrequency / 2.0
        var magnitudes: [T] = []
        var phases: [T] = []

        for freq in frequencies {
            let normalizedFreq = freq / nyquist
            let omega = .pi * normalizedFreq

            let response = calculateFilterResponse(omega: omega, filterType: filterType, order: order)
            magnitudes.append(T(response.magnitude))
            phases.append(T(response.phase))
        }

        return (magnitudes, phases)
    }

    private func calculateFilterResponse(
        omega: Double,
        filterType: WaveformFilterType,
        order: Int
    ) -> (magnitude: Double, phase: Double) {
        // Simplified frequency response calculation
        // For demonstration purposes - a complete implementation would use proper transfer functions

        switch filterType {
        case .lowPass(let cutoffFreq):
            let wc = 2.0 * .pi * cutoffFreq / samplingFrequency
            let ratio = omega / wc
            let magnitude = 1.0 / sqrt(1.0 + pow(ratio, 2.0 * Double(order)))
            let phase = -Double(order) * atan(ratio)
            return (magnitude, phase)

        case .highPass(let cutoffFreq):
            let wc = 2.0 * .pi * cutoffFreq / samplingFrequency
            let ratio = wc / omega
            let magnitude = 1.0 / sqrt(1.0 + pow(ratio, 2.0 * Double(order)))
            let phase = Double(order) * atan(ratio)
            return (magnitude, phase)

        default:
            return (1.0, 0.0)  // Placeholder for other filter types
        }
    }
}

// MARK: - Savitzky-Golay and Polynomial Filters
extension Waveform1D where T: BinaryFloatingPoint {

    /// Apply Savitzky-Golay filter for smoothing and derivative estimation
    /// - Parameters:
    ///   - windowSize: Size of the filter window (must be odd and >= 3)
    ///   - polynomialOrder: Order of the fitting polynomial (must be < windowSize)
    ///   - derivative: Derivative order (0 = smoothing, 1 = first derivative, etc.)
    /// - Returns: Filtered waveform
    public func savitzkyGolayFilter(
        windowSize: Int,
        polynomialOrder: Int,
        derivative: Int = 0
    ) -> Waveform1D<T> {
        guard
            windowSize >= 3 && windowSize % 2 == 1 && polynomialOrder >= 0 && polynomialOrder < windowSize
                && derivative >= 0 && derivative <= polynomialOrder
        else {
            return self
        }

        let coefficients = calculateSavitzkyGolayCoefficients(
            windowSize: windowSize,
            polynomialOrder: polynomialOrder,
            derivative: derivative
        )

        return applySavitzkyGolayFilter(coefficients: coefficients, windowSize: windowSize, derivative: derivative)
    }

    /// Apply Savitzky-Golay smoothing filter (derivative = 0)
    /// - Parameters:
    ///   - windowSize: Size of the filter window (must be odd)
    ///   - polynomialOrder: Order of the fitting polynomial
    /// - Returns: Smoothed waveform
    public func savitzkyGolaySmooth(windowSize: Int = 5, polynomialOrder: Int = 2) -> Waveform1D<T> {
        return savitzkyGolayFilter(windowSize: windowSize, polynomialOrder: polynomialOrder, derivative: 0)
    }

    /// Compute Savitzky-Golay first derivative
    /// - Parameters:
    ///   - windowSize: Size of the filter window (must be odd)
    ///   - polynomialOrder: Order of the fitting polynomial
    /// - Returns: First derivative waveform
    public func savitzkyGolayDerivative(windowSize: Int = 5, polynomialOrder: Int = 2) -> Waveform1D<T> {
        let result = savitzkyGolayFilter(windowSize: windowSize, polynomialOrder: polynomialOrder, derivative: 1)
        // Scale by 1/dt to get proper derivative units
        let scaledValues = result.values.map { $0 / T(dt) }
        return Waveform1D(values: scaledValues, dt: dt, t0: t0)
    }

    /// Compute Savitzky-Golay second derivative
    /// - Parameters:
    ///   - windowSize: Size of the filter window (must be odd)
    ///   - polynomialOrder: Order of the fitting polynomial (must be >= 2)
    /// - Returns: Second derivative waveform
    public func savitzkyGolaySecondDerivative(windowSize: Int = 5, polynomialOrder: Int = 3) -> Waveform1D<T> {
        guard polynomialOrder >= 2 else { return self }
        let result = savitzkyGolayFilter(windowSize: windowSize, polynomialOrder: polynomialOrder, derivative: 2)
        // Scale by 1/dt² to get proper derivative units
        let scaledValues = result.values.map { $0 / T(dt * dt) }
        return Waveform1D(values: scaledValues, dt: dt, t0: t0)
    }

    /// Apply local polynomial regression filter (LOESS-style)
    /// - Parameters:
    ///   - windowSize: Size of the local fitting window
    ///   - polynomialOrder: Order of the fitting polynomial
    ///   - weights: Optional weight function for local regression
    /// - Returns: Filtered waveform
    public func localPolynomialFilter(
        windowSize: Int,
        polynomialOrder: Int,
        weights: WaveformLocalWeightFunction = .uniform
    ) -> Waveform1D<T> {
        guard windowSize >= polynomialOrder + 1 && windowSize % 2 == 1 else { return self }

        var filtered: [T] = []
        let halfWindow = windowSize / 2

        for i in 0..<values.count {
            let startIdx = max(0, i - halfWindow)
            let endIdx = min(values.count, i + halfWindow + 1)
            let windowValues = Array(values[startIdx..<endIdx])

            // Create local coordinate system centered at current point
            let localX = (startIdx..<endIdx).map { T($0 - i) }

            // Generate weights
            let weightValues = generateLocalWeights(
                distances: localX.map { abs($0) },
                maxDistance: T(halfWindow),
                function: weights
            )

            // Fit polynomial and evaluate at center point
            if let fittedValue = fitLocalPolynomial(
                x: localX,
                y: windowValues,
                weights: weightValues,
                order: polynomialOrder,
                evaluateAt: T.zero
            ) {
                filtered.append(fittedValue)
            } else {
                filtered.append(values[i])  // Fallback to original value
            }
        }

        return Waveform1D(values: filtered, dt: dt, t0: t0)
    }

    /// Apply Whittaker-Henderson filter for smoothing
    /// - Parameters:
    ///   - lambda: Smoothing parameter (higher = more smoothing)
    ///   - order: Difference order (typically 2)
    /// - Returns: Smoothed waveform
    public func whittakerHendersonFilter(lambda: Double, order: Int = 2) -> Waveform1D<T> {
        guard order > 0 && lambda > 0 && values.count > order else { return self }

        let smoothed = solveWhittakerSystem(
            values: values.map { Double($0) },
            lambda: lambda,
            order: order
        )

        let result = smoothed.map { T($0) }
        return Waveform1D(values: result, dt: dt, t0: t0)
    }

    // MARK: - Private Implementation Methods

    private func calculateSavitzkyGolayCoefficients(
        windowSize: Int,
        polynomialOrder: Int,
        derivative: Int
    ) -> [T] {
        let halfWindow = windowSize / 2
        let n = windowSize
        let order = polynomialOrder

        // Create matrix for least squares polynomial fitting
        var A: [[Double]] = Array(repeating: Array(repeating: 0.0, count: order + 1), count: n)

        // Fill design matrix
        for i in 0..<n {
            let x = Double(i - halfWindow)
            for j in 0...order {
                A[i][j] = pow(x, Double(j))
            }
        }

        // Solve normal equations: (A^T A) c = A^T e_k
        // where e_k is the k-th unit vector (for derivative order)
        let coefficients = solveSavitzkyGolaySystem(A: A, derivative: derivative)

        // Apply factorial scaling for derivatives
        let derivativeScale = factorial(derivative)
        return coefficients.map { T($0 * derivativeScale) }
    }

    private func solveSavitzkyGolaySystem(A: [[Double]], derivative: Int) -> [Double] {
        let n = A.count
        let m = A[0].count

        // Compute A^T A
        var ATA: [[Double]] = Array(repeating: Array(repeating: 0.0, count: m), count: m)
        for i in 0..<m {
            for j in 0..<m {
                for k in 0..<n {
                    ATA[i][j] += A[k][i] * A[k][j]
                }
            }
        }

        // Compute A^T e_derivative (right-hand side)
        var rhs: [Double] = Array(repeating: 0.0, count: m)
        for i in 0..<n {
            if derivative < m {
                rhs[derivative] += A[i][derivative]
            }
        }

        // Solve linear system using Gaussian elimination
        let polynomialCoeffs = gaussianElimination(matrix: ATA, rhs: rhs)

        // Convert polynomial coefficients to filter coefficients
        var filterCoeffs: [Double] = Array(repeating: 0.0, count: n)
        let halfWindow = n / 2

        for i in 0..<n {
            let x = Double(i - halfWindow)
            for j in 0..<polynomialCoeffs.count {
                filterCoeffs[i] += polynomialCoeffs[j] * pow(x, Double(j))
            }
        }

        return filterCoeffs
    }

    private func applySavitzkyGolayFilter(coefficients: [T], windowSize: Int, derivative: Int) -> Waveform1D<T> {
        guard coefficients.count == windowSize else { return self }

        var filtered: [T] = []
        let halfWindow = windowSize / 2

        for i in 0..<values.count {
            var sum = T.zero

            for j in 0..<windowSize {
                let sampleIndex = i - halfWindow + j
                let clampedIndex = max(0, min(values.count - 1, sampleIndex))
                sum += coefficients[j] * values[clampedIndex]
            }

            filtered.append(sum)
        }

        return Waveform1D(values: filtered, dt: dt, t0: t0)
    }

    private func generateLocalWeights(
        distances: [T],
        maxDistance: T,
        function: WaveformLocalWeightFunction
    ) -> [T] {
        switch function {
        case .uniform:
            return Array(repeating: T(1.0), count: distances.count)

        case .tricube:
            return distances.map { d in
                let normalized = d / maxDistance
                if normalized <= T(1.0) {
                    let u = T(1.0) - normalized * normalized * normalized
                    return u * u * u
                } else {
                    return T.zero
                }
            }

        case .gaussian(let sigma):
            let sigmaT = T(sigma)
            return distances.map { d in
                let normalized = d / sigmaT
                return T(exp(-0.5 * Double(normalized * normalized)))
            }

        case .epanechnikov:
            return distances.map { d in
                let normalized = d / maxDistance
                if normalized <= T(1.0) {
                    return T(0.75) * (T(1.0) - normalized * normalized)
                } else {
                    return T.zero
                }
            }
        }
    }

    private func fitLocalPolynomial(
        x: [T],
        y: [T],
        weights: [T],
        order: Int,
        evaluateAt: T
    ) -> T? {
        guard x.count == y.count && y.count == weights.count && x.count > order else {
            return nil
        }

        let n = x.count
        let m = order + 1

        // Create weighted design matrix
        var A: [[Double]] = Array(repeating: Array(repeating: 0.0, count: m), count: n)
        var weightedY: [Double] = []

        for i in 0..<n {
            let weight = sqrt(Double(weights[i]))
            weightedY.append(Double(y[i]) * weight)

            for j in 0..<m {
                A[i][j] = pow(Double(x[i]), Double(j)) * weight
            }
        }

        // Solve weighted least squares
        let coefficients = solveLeastSquares(A: A, b: weightedY)

        // Evaluate polynomial at target point
        var result = 0.0
        for j in 0..<coefficients.count {
            result += coefficients[j] * pow(Double(evaluateAt), Double(j))
        }

        return T(result)
    }

    private func solveWhittakerSystem(values: [Double], lambda: Double, order: Int) -> [Double] {
        let n = values.count

        // Create difference matrix D
        let D = createDifferenceMatrix(size: n, order: order)

        // Create weight matrix W (identity for simple case)
        var W = Array(repeating: Array(repeating: 0.0, count: n), count: n)
        for i in 0..<n {
            W[i][i] = 1.0
        }

        // Solve (W + λD^T D) z = W y
        var system = W
        let DTD = matrixMultiply(transpose(D), D)

        // Add λD^T D to W
        for i in 0..<n {
            for j in 0..<n {
                system[i][j] += lambda * DTD[i][j]
            }
        }

        return gaussianElimination(matrix: system, rhs: values)
    }

    // MARK: - Linear Algebra Helpers

    private func gaussianElimination(matrix: [[Double]], rhs: [Double]) -> [Double] {
        let n = matrix.count
        guard n == rhs.count && n > 0 else { return [] }

        var A = matrix
        var b = rhs

        // Forward elimination
        for i in 0..<(n - 1) {
            // Find pivot
            var maxRow = i
            for k in (i + 1)..<n {
                if abs(A[k][i]) > abs(A[maxRow][i]) {
                    maxRow = k
                }
            }

            // Swap rows
            A.swapAt(i, maxRow)
            b.swapAt(i, maxRow)

            // Make all rows below this one 0 in current column
            for k in (i + 1)..<n {
                if A[i][i] != 0 {
                    let factor = A[k][i] / A[i][i]
                    for j in i..<n {
                        A[k][j] -= factor * A[i][j]
                    }
                    b[k] -= factor * b[i]
                }
            }
        }

        // Back substitution
        var x = Array(repeating: 0.0, count: n)
        for i in stride(from: n - 1, through: 0, by: -1) {
            x[i] = b[i]
            for j in (i + 1)..<n {
                x[i] -= A[i][j] * x[j]
            }
            if A[i][i] != 0 {
                x[i] /= A[i][i]
            }
        }

        return x
    }

    private func solveLeastSquares(A: [[Double]], b: [Double]) -> [Double] {
        // Solve normal equations A^T A x = A^T b
        let AT = transpose(A)
        let ATA = matrixMultiply(AT, A)
        let ATb = matrixVectorMultiply(AT, b)

        return gaussianElimination(matrix: ATA, rhs: ATb)
    }

    private func transpose(_ matrix: [[Double]]) -> [[Double]] {
        guard !matrix.isEmpty else { return [] }
        let rows = matrix.count
        let cols = matrix[0].count

        var result: [[Double]] = Array(repeating: Array(repeating: 0.0, count: rows), count: cols)
        for i in 0..<rows {
            for j in 0..<cols {
                result[j][i] = matrix[i][j]
            }
        }
        return result
    }

    private func matrixMultiply(_ A: [[Double]], _ B: [[Double]]) -> [[Double]] {
        guard !A.isEmpty && !B.isEmpty && A[0].count == B.count else { return [] }

        let rows = A.count
        let cols = B[0].count
        let inner = A[0].count

        var result: [[Double]] = Array(repeating: Array(repeating: 0.0, count: cols), count: rows)

        for i in 0..<rows {
            for j in 0..<cols {
                for k in 0..<inner {
                    result[i][j] += A[i][k] * B[k][j]
                }
            }
        }

        return result
    }

    private func matrixVectorMultiply(_ A: [[Double]], _ b: [Double]) -> [Double] {
        guard !A.isEmpty && A[0].count == b.count else { return [] }

        let rows = A.count
        var result: [Double] = Array(repeating: 0.0, count: rows)

        for i in 0..<rows {
            for j in 0..<b.count {
                result[i] += A[i][j] * b[j]
            }
        }

        return result
    }

    private func createDifferenceMatrix(size: Int, order: Int) -> [[Double]] {
        guard order > 0 && size > order else { return [] }

        // Start with first-order differences
        var currentD: [[Double]] = Array(repeating: Array(repeating: 0.0, count: size), count: size - 1)
        for i in 0..<(size - 1) {
            currentD[i][i] = -1.0
            currentD[i][i + 1] = 1.0
        }

        // Apply higher-order differences
        for _ in 1..<order {
            let newRows = currentD.count - 1
            var nextD: [[Double]] = Array(repeating: Array(repeating: 0.0, count: size), count: newRows)

            for i in 0..<newRows {
                for j in 0..<size {
                    nextD[i][j] = currentD[i + 1][j] - currentD[i][j]
                }
            }
            currentD = nextD
        }

        return currentD
    }

    private func factorial(_ n: Int) -> Double {
        guard n >= 0 else { return 0.0 }
        guard n > 1 else { return 1.0 }

        var result = 1.0
        for i in 2...n {
            result *= Double(i)
        }
        return result
    }
}
