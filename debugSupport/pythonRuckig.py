"""Main module."""

from copy import copy

from ruckig import InputParameter, OutputParameter, Result, Ruckig


if __name__ == "__main__":
    # Create instances: the Ruckig OTG as well as input and output parameters
    otg = Ruckig(1, 0.01)  # DoFs, control cycle
    inp = InputParameter(1)
    out = OutputParameter(1)

    inp.current_position = [2.6603799559471364]
    inp.current_velocity = [-2.5442249449339207]
    inp.current_acceleration = [0.0]

    inp.target_position = [-6.297064289647577]
    inp.target_velocity = [0.0]
    inp.target_acceleration = [0.0]

    inp.max_velocity = [4.0]
    inp.max_acceleration = [5.0]
    inp.max_jerk = [10.0]

    print("\t".join(["t"] + [str(i) for i in range(otg.degrees_of_freedom)]))

    # Generate the trajectory within the control loop
    first_output, out_list = None, []
    res = Result.Working
    while res == Result.Working:
        res = otg.update(inp, out)

        print("\t".join([f"{out.time:0.3f}"] + [f"{p:0.3f}" for p in out.new_position]))
        out_list.append(copy(out))

        out.pass_to_input(inp)

        if not first_output:
            first_output = copy(out)

    print(f"\nCalculation duration: {first_output.calculation_duration:0.1f} [µs]")
    print(f"Trajectory duration: {first_output.trajectory.duration:0.4f} [s]")

    # Print profile details
    print("\nProfile for DoF 0:")
    profile = first_output.trajectory.profiles[0]
    print(f"  Profile type: {profile}")

    # Check what intermediate_durations actually contains
    print("\nIntermediate durations structure:")
    print(f"  Type: {type(first_output.trajectory.intermediate_durations)}")
    print(f"  Value: {first_output.trajectory.intermediate_durations}")

    # Sample the trajectory at key points to extract acceleration values
    print("\nSampling trajectory at key times:")
    sample_times = [0.0, 0.01, 0.381546, 0.763092, 1.728222, 2.228222, 2.528222, 3.028222]
    for t in sample_times:
        if t <= first_output.trajectory.duration:
            new_pos, new_vel, new_acc = first_output.trajectory.at_time(t)
            print(f"  t={t:0.6f}: p={new_pos[0]:8.4f}, v={new_vel[0]:8.4f}, a={new_acc[0]:8.4f}")

    print(f"\nAcceleration limits: aMin={-inp.max_acceleration[0]}, aMax={inp.max_acceleration[0]}")
