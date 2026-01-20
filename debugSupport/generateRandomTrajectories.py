"""
Generate random trajectories and save to JSON files matching Swift Codable schema.

This script generates random trajectory input parameters, tests them with Ruckig,
and saves the results to JSON files that match the InputParameter+Codable.swift schema.

Usage:
    # Generate 1000 random 1-DOF trajectories (default)
    python3 debugSupport/generateRandomTrajectories.py

    # Generate custom number of trajectories with reproducible seed
    # Edit main() call at bottom of file with desired parameters:
    # main(dof=1, num_trajectories=5000, seed=42)

Output files:
    - debugSupport/successful_trajectories.json: Trajectories that Ruckig could solve
    - debugSupport/failed_trajectories.json: Trajectories that Ruckig could not solve (includes error message)

JSON Schema (matches InputParameter+Codable.swift):
    - All required fields: degrees_of_freedom, control_interface, synchronization, etc.
    - Arrays match DOF count
    - Enums as strings: "Position"/"Velocity", "Time"/"TimeIfNecessary"/"Phase"/"None", etc.
"""

import json
import random
from typing import Dict, List, Any

from ruckig import InputParameter, OutputParameter, Result, Ruckig


def generate_random_input(dof: int = 1) -> InputParameter:
    """Generate a random InputParameter for testing."""
    inp = InputParameter(dof)

    # Generate kinematic limits FIRST
    inp.max_velocity = [random.uniform(2, 8) for _ in range(dof)]
    inp.max_acceleration = [random.uniform(3, 10) for _ in range(dof)]
    inp.max_jerk = [random.uniform(5, 20) for _ in range(dof)]

    # Random current state - respecting limits
    inp.current_position = [random.uniform(-10, 10) for _ in range(dof)]
    inp.current_velocity = [random.uniform(-inp.max_velocity[i], inp.max_velocity[i]) for i in range(dof)]
    inp.current_acceleration = [random.uniform(-inp.max_acceleration[i], inp.max_acceleration[i]) for i in range(dof)]

    # Random target state - respecting limits
    inp.target_position = [random.uniform(-10, 10) for _ in range(dof)]
    inp.target_velocity = [random.uniform(-inp.max_velocity[i], inp.max_velocity[i]) for i in range(dof)]
    inp.target_acceleration = [random.uniform(-inp.max_acceleration[i], inp.max_acceleration[i]) for i in range(dof)]

    return inp


def input_to_dict(inp: InputParameter, dof: int) -> Dict[str, Any]:
    """Convert InputParameter to dictionary matching Swift Codable schema."""
    return {
        "degrees_of_freedom": dof,
        "control_interface": "Position",  # Default value
        "synchronization": "Time",  # Default value
        "duration_discretization": "Continuous",  # Default value

        "current_position": inp.current_position,
        "current_velocity": inp.current_velocity,
        "current_acceleration": inp.current_acceleration,

        "target_position": inp.target_position,
        "target_velocity": inp.target_velocity,
        "target_acceleration": inp.target_acceleration,

        "max_velocity": inp.max_velocity,
        "max_acceleration": inp.max_acceleration,
        "max_jerk": inp.max_jerk,

        "intermediate_positions": [],  # Empty for simple trajectories
        "enabled": [True] * dof,  # All DOFs enabled
    }


def test_trajectory(inp: InputParameter, dof: int, otg: Ruckig) -> tuple[bool, Dict[str, Any], str]:
    """
    Test a trajectory with Ruckig and validate it respects limits.

    Returns:
        (success, input_dict, error_message)
    """
    out = OutputParameter(dof)
    input_dict = input_to_dict(inp, dof)

    try:
        result = otg.update(inp, out)

        if result == Result.Working or result == Result.Finished:
            # Validate that trajectory actually respects limits
            duration = out.trajectory.duration
            max_vel = inp.max_velocity[0]
            max_acc = inp.max_acceleration[0]

            # Sample trajectory at 100 points
            for i in range(101):
                t = duration * i / 100.0
                pos, vel, acc = out.trajectory.at_time(t)

                if abs(vel[0]) > max_vel + 0.001:
                    return (False, input_dict, f"Velocity limit violated: {abs(vel[0]):.4f} > {max_vel:.4f} at t={t:.4f}")
                if abs(acc[0]) > max_acc + 0.001:
                    return (False, input_dict, f"Acceleration limit violated: {abs(acc[0]):.4f} > {max_acc:.4f} at t={t:.4f}")

            # Success - trajectory calculated and respects all limits
            return (True, input_dict, "")
        else:
            # Failed - could not calculate trajectory
            error_msg = str(result)
            return (False, input_dict, error_msg)
    except Exception as e:
        # Exception during calculation
        return (False, input_dict, str(e))


def main(dof: int = 1, num_trajectories: int = 1000, seed: int | None = None):
    """
    Generate random trajectories and save to JSON files.

    Args:
        dof: Degrees of freedom (default: 1)
        num_trajectories: Number of random trajectories to generate (default: 1000)
        seed: Random seed for reproducibility (default: None)
    """
    if seed is not None:
        random.seed(seed)

    # Configuration
    control_cycle = 0.01

    otg = Ruckig(dof, control_cycle)

    successful_trajectories: List[Dict[str, Any]] = []
    failed_trajectories: List[Dict[str, Any]] = []

    print(f"Generating {num_trajectories} random trajectories...")

    for i in range(num_trajectories):
        if (i + 1) % 100 == 0:
            print(f"  Progress: {i + 1}/{num_trajectories}")

        inp = generate_random_input(dof)
        success, input_dict, error_msg = test_trajectory(inp, dof, otg)

        if success:
            successful_trajectories.append(input_dict)
        else:
            # Add error message to failed trajectories
            input_dict["error"] = error_msg
            failed_trajectories.append(input_dict)

    # Save to JSON files
    success_file = "debugSupport/successful_trajectories.json"
    failed_file = "debugSupport/failed_trajectories.json"

    with open(success_file, 'w') as f:
        json.dump(successful_trajectories, f, indent=2)

    with open(failed_file, 'w') as f:
        json.dump(failed_trajectories, f, indent=2)

    # Print summary
    print(f"\nResults:")
    print(f"  Successful: {len(successful_trajectories)} ({len(successful_trajectories)/num_trajectories*100:.1f}%)")
    print(f"  Failed: {len(failed_trajectories)} ({len(failed_trajectories)/num_trajectories*100:.1f}%)")
    print(f"\nFiles saved:")
    print(f"  {success_file}")
    print(f"  {failed_file}")

    # Show some failed trajectory reasons if any
    if failed_trajectories:
        print(f"\nSample failure reasons:")
        error_counts: Dict[str, int] = {}
        for traj in failed_trajectories:
            error = traj.get("error", "Unknown")
            error_counts[error] = error_counts.get(error, 0) + 1

        for error, count in sorted(error_counts.items(), key=lambda x: -x[1]):
            print(f"  {error}: {count} occurrences")


if __name__ == "__main__":
    main()
