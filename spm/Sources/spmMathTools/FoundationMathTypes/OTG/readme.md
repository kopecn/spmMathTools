# Swift based Multi DoF Jerk Limited Profiles

- Based upon [Ruckig](https://github.com/pantor/ruckig)

### Major Control Flow

```mermaid
flowchart TD
    A[Start: calculate] --> B[Reset was_interrupted, maybe clear traj]
    B --> C[Loop over DoFs]

    C -->|if DoF disabled| D[Copy current state to trajectory]
    D --> C

    C -->|else| E[Set limits and control settings]
    E --> F[Compute brake trajectory]
    F --> G[Set boundary conditions]
    G --> H[Finalize brake profile]
    H --> I[Step 1: Try finding profile]

    I -->|if not found| J[Check zero limits?]
    J -->|yes| K[Return ErrorZeroLimits]
    J -->|no| L[Return ErrorExecutionTimeCalculation]

    I -->|found| M[Store t_min to traj.independent_min_durations]
    M --> C

    C -->|All DoFs done| N[Check if 1 DoF & no min duration & continuous]
    N -->|yes| O[Assign duration & Return Result::Working]

    N -->|no| P[Try synchronize across DoFs]

    P -->|failed| Q[Check zero limits?]
    Q -->|yes| R[Return ErrorZeroLimits]
    Q -->|no| S[Return ErrorSynchronizationCalculation]

    P -->|success| T[Handle Synchronization::None]
    T --> U[Update traj.duration if needed]

    U --> V[Check if duration == 0]
    V -->|yes| W[Copy p_min to all profiles → Return Result::Working]

    V -->|no| X[Check if all sync == None and continuous]
    X -->|yes| Y[Return Result::Working]

    X -->|no| Z[Check if Phase Sync is needed]
    Z -->|true| AA[Try to align all DoFs to limiting DoF]

    AA -->|success| AB[Return Result::Working]

    AA -->|fail or not needed| AC[Loop over DoFs for Time Sync]
    AC -->|if already synced| AD[Skip or copy profile]
    AD --> AC

    AC -->|else| AE[Try Step 2 to match timing]
    AE -->|fail| AF[Return ErrorSynchronizationCalculation]
    AE -->|success| AG[Store profile]

    AG --> AC

    AC -->|All done| AH[Return Result::Working]
```

## Citation

If you use this project or build upon it, please cite the following paper:

**Jerk-limited Real-time Trajectory Generation with Arbitrary Target States**  
Lars Berscheid, Torsten Kröger  
_Robotics: Science and Systems XVII_, 2021.

[BibTeX](#bibtex)

---

### BibTeX

```bibtex
@article{berscheid2021jerk,
  title={Jerk-limited Real-time Trajectory Generation with Arbitrary Target States},
  author={Berscheid, Lars and Kr{\"o}ger, Torsten},
  journal={Robotics: Science and Systems XVII},
  year={2021}
}
```

### Notes on Benchmarking:

Cpp Ruckig on Apple M1 Max, 64Gb ~ 2021 provided:

```bash
  Benchmark for 3 DoFs on 262144 trajectories
  Average Calculation Duration 2.427 pm 0.0122862 [µs]
  Worst Calculation Duration 65.1749 pm 19.3034 [µs]
  End-to-end Calculation Duration 2.98597 pm 0.0150864 [µs]
```

Swift Ruckig on Apple M1 Max, 64Gb ~ 2021 provided:

```bash
  Benchmark for 3 DoFs on 262144 trajectories
  Average Calculation Duration 73.65 ± 1.18 [µs]
  Worst Calculation Duration 2397.72 ± 4981.77 [µs]
  End-to-end Calculation Duration 76.79 ± 1.25 [µs]
```

Suspect there are failing conditions which causes the Worst Calculation to jump significantly into the 2.4ms range.

to recreate the Swift benchmarks run:

```bash
make benchmark
```

and to recreate the cPP benchmarks run:

```bash
mkdir -p build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DBUILD_BENCHMARK=ON -DBUILD_TESTS=ON -DBUILD_EXAMPLES=ON ..
make
./benchmark-target
```


--- 

# Class-Functional Hierarchy


InputParameter + OutputParameter (inout) --> OTG

OTG --> OutputParameter/Trajectory --> Calculator/TargetCalculator (primary calculator)

TargetCalculator --> PositionFirstOrderStep1 ... PositionThirdOrderStep2



