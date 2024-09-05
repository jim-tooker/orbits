"""
Tests for Orbit Simulator
"""
import sys
import argparse
import pytest
from orbits.motion_tracker import MotionType
from orbits.orbit_simulator import OrbitSimulator
from orbits import config


@pytest.mark.earth
def test_earth_rotation_period() -> None:
    """Verify the Earth's rotational period matches the defined value."""
    # Relative tolerance for test results
    rel_tolerance = 0.001

    # Initialize simulator
    time_scale_factor: float = 10_000
    config.time_scale_factor = time_scale_factor
    sim: OrbitSimulator = OrbitSimulator()

    # Run long enough for 3 revolutions plus a margin
    num_of_revolutions: float = 3
    time_margin: float = 0.5
    run_time: float = ((sim.mode.earth.params.rotation_period / time_scale_factor) * num_of_revolutions) + time_margin

    # Run it
    sim.run(run_time)

    # Check if Earth's rotation time matches given rotation period parameters
    assert sim.mode.tracker.full_angle_times[MotionType.EARTH_ROTATION] == \
        pytest.approx(sim.mode.earth.params.rotation_period, rel=rel_tolerance)

@pytest.mark.moon
def test_moon_orbital_and_rotation_period() -> None:
    """
    Verify the Moon's orbital period matches the defined value.  
    Verify the Moon's rotational period matches the defined value.  
    Verify the number of Earth's rotations within one Moon's orbit matches the defined value.  
    """
    # Relative tolerance for test results
    rel_tolerance = 0.1

    # Initialize simulator
    time_scale_factor: float = 500_000
    config.time_scale_factor = time_scale_factor
    sim: OrbitSimulator = OrbitSimulator()

    # Run long enough for one orbit
    time_margin: float = 0.5
    run_time: float = (sim.mode.moon.orbit.params.period / time_scale_factor) + time_margin

    # Run it
    sim.run(run_time)

    # Check if Moon's orbit time matches given orbit time parameters
    assert sim.mode.tracker.full_angle_times[MotionType.MOON_ORBIT] == \
        pytest.approx(sim.mode.moon.orbit.params.period, rel=rel_tolerance)

    # Check if Moon's rotation time matches given rotation period parameters
    assert sim.mode.tracker.full_angle_times[MotionType.MOON_ROTATION] == \
        pytest.approx(sim.mode.moon.params.rotation_period, rel=rel_tolerance)

    # Check if the number of Earth's revolutions during one Moon's orbit matches given parameters
    assert sim.mode.tracker.totals[MotionType.EARTH_ROTATION] == \
        pytest.approx(sim.mode.moon.orbit.params.period_days, rel=rel_tolerance)

@pytest.mark.earth
def test_earth_orbital_period() -> None:
    """
    Verify the Earth's orbital period matches the defined value.  
    """
    # Relative tolerance for test results
    rel_tolerance = 0.001

    # Initialize simulator
    time_scale_factor: float = 50_000
    config.time_scale_factor = time_scale_factor
    sim: OrbitSimulator = OrbitSimulator()

    # Run long enough for one orbit
    time_margin: float = 0.5
    run_time: float = (sim.mode.earth.orbit.params.period / time_scale_factor) + time_margin

    # Run it
    sim.run(run_time)

    # Check if Earth's orbit time matches given orbit time parameters
    assert sim.mode.tracker.full_angle_times[MotionType.EARTH_ORBIT] == \
        pytest.approx(sim.mode.earth.orbit.params.period, rel=rel_tolerance)

    # Check if the number of Earth's revolutions during one Earth's orbit matches given parameters
    assert sim.mode.tracker.totals[MotionType.EARTH_ROTATION] == \
        pytest.approx(sim.mode.earth.orbit.params.period_days, rel=rel_tolerance)

@pytest.mark.sun
def test_sun_rotation_period() -> None:
    """Verify the Sun's rotational period matches the defined value."""
    # Relative tolerance for test results
    rel_tolerance = 0.001

    # Initialize simulator
    time_scale_factor: float = 100_000
    config.time_scale_factor = time_scale_factor
    sim: OrbitSimulator = OrbitSimulator()

    # Run long enough for 3 revolutions plus a margin
    num_of_revolutions: float = 3
    time_margin: float = 0.5
    run_time: float = ((sim.mode.sun.params.rotation_period / time_scale_factor) * num_of_revolutions) + time_margin

    # Run it
    sim.run(run_time)

    # Check if Sun's rotation time matches given rotation period parameters
    assert sim.mode.tracker.full_angle_times[MotionType.SUN_ROTATION] == \
        pytest.approx(sim.mode.sun.params.rotation_period, rel=rel_tolerance)


def main() -> None:
    """Main entry point for test."""
    parser = argparse.ArgumentParser(description='Orbit Simulator Tester')
    parser.add_argument('--no_gui', action='store_true', help='Run without GUI')
    parser.add_argument('--test', help='Specify a test to run')
    parser.add_argument('-k', help='Only run tests which match the given substring expression')
    parser.add_argument('-m', help='Only run tests matching given mark expression')
    args = parser.parse_args()

    if args.no_gui:
        OrbitSimulator.disable_gui(True)

    pytest_args = ["tests/test_orbit_simulator.py"]
    if args.test:
        pytest_args.append(f"::test_{args.test}")
    if args.k:
        pytest_args.extend(["-k", args.k])
    if args.m:
        pytest_args.extend(["-m", args.m])

    result = pytest.main(pytest_args)

    if args.no_gui:
        sys.exit(result)
    else:
        OrbitSimulator.quit_simulation()

if __name__ == '__main__':
    main()
