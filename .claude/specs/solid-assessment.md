# spmMathTools — SOLID Principles Assessment

## Single Responsibility (SRP)

| Area | Rating | Notes |
|------|--------|-------|
| Waveform1D extensions | Good | Each file owns one concern (calc, generators, operators) |
| Spatial operators | Good | Clean separation by operator category |
| Polynomial hierarchy | Good | Each degree is its own class with focused root logic |
| **OTG Profile** | **Poor** | One 735-line struct handles: boundary conditions, state validation, timing checks across 3 interface types (position/velocity, 1st/2nd/3rd order), extrema calculation, brake sub-profiles. Should split into `ProfileValidator`, `ProfileComputer`, `ProfileBoundary` |
| **OTG CalculatorTarget** | **Fair** | 765 lines with clear structure but mixes trivial-case detection, per-DOF stepping, synchronization, and phase sync into one method |
| **OTG Trajectory** | **Fair** | Mixes storage, querying (6 `atTime()` overloads), and extrema computation |

## Open/Closed (OCP)

| Area | Rating | Notes |
|------|--------|-------|
| PolynomialUnivariateOrder | Good | `construct()` factory + `.higher(Int)` case allows extension |
| Waveform1D operators | Good | Extension-based design is inherently open for extension |
| **OTG Block.calculateBlock()** | **Poor** | 5 hardcoded conditional branches; adding a new heuristic means modifying the method |
| **OTG control interface dispatch** | **Fair** | switch on `ControlInterface` in multiple places; adding a 3rd interface type requires touching many files |

## Liskov Substitution (LSP)

| Area | Rating | Notes |
|------|--------|-------|
| Polynomial hierarchy | Good | Each subclass correctly overrides `roots` with appropriate fallback to lower degree |
| OTG types | N/A | No inheritance hierarchy (flat struct design) |
| Spatial types | Good | Generic constraints ensure substitutability |

## Interface Segregation (ISP)

| Area | Rating | Notes |
|------|--------|-------|
| Waveform1D | Excellent | 40+ focused support types, each modeling one concept |
| Spatial types | Good | Operators cleanly separated by category |
| **OTG Profile** | **Poor** | No protocol decomposition; clients must depend on the entire 735-line struct even if they only need timing checks |
| **OTG InputParameter** | **Fair** | Large struct (30+ fields) but logically cohesive for trajectory input |

## Dependency Inversion (DIP)

| Area | Rating | Notes |
|------|--------|-------|
| spmMathTools → FoundationTypes | Good | Depends on foundational abstractions, not concretions |
| Waveform1D extensions | Good | Extend generic type, don't depend on concrete implementations |
| **OTG internal** | **Poor** | `TargetCalculator`, `Block`, `Trajectory` all depend on concrete `Profile` struct. No protocol abstraction for profile behavior. |
| **Root solvers** | **Fair** | Free functions (`solveCubic`, `solveQuarticMonic`) are usable independently but `insertIfPositive` is an implicit dependency via Array extension |

## Overall SOLID Grade: C+

**Strengths:** The extension-based architecture for Waveform1D, Spatial, and Polynomial types is well-structured and follows SOLID naturally. The separation between `spmFoundationTools` (types) and `spmMathTools` (operations) is a sound DIP pattern.

**Weaknesses concentrated in OTG:** The Ruckig port carries over C++ structural patterns (switch-based dispatch, large types mixing concerns) that don't fully align with Swift's protocol-oriented idioms. OTG types (Profile, Block, Trajectory) have been converted to value types (structs), eliminating manual copy semantics — but further protocol decomposition would yield additional SOLID and performance benefits toward closing the 30x performance gap with C++.
