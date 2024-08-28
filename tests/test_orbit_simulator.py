"""
Tests for Orbit Simulator
"""
import sys
import argparse
import pytest
from orbits.celestial_body import MotionType
from orbits.orbit_simulator import OrbitSimulator


class TestOrbitSimulation:
    """Class for `pytest` testing of Orbit Simulator."""

    def test_earth_rotation_period(self) -> None:
        """Verify the Earth's rotational period matches the defined value."""
        # Relative tolerance for test results
        rel_tolerance = 0.001

        # Initialize simulator
        time_scale_factor: float = 10_000
        sim: OrbitSimulator = OrbitSimulator(time_scale_factor=time_scale_factor)

        # Run long enough for 3 revolutions
        num_of_revolutions: float = 3
        time_margin: float = 0.5
        run_time: float = ((sim.earth.params.rotation_period / time_scale_factor) * num_of_revolutions) + time_margin

        # Run it
        sim.run(run_time)

        # Check if Earth's rotation time matches given rotation period parameters
        assert sim.tracker.full_angle_times[MotionType.EARTH_ROTATION] == \
            pytest.approx(sim.earth.params.rotation_period, rel=rel_tolerance)

    def test_moon_orbital_and_rotation_period(self) -> None:
        """
        Verify the Moon's orbital period matches the defined value.  
        Verify the Moon's rotational period matches the defined value.  
        Verify the number of Earth's rotations within one Moon's orbit matches the defined value.  
        """
        # Relative tolerance for test results
        rel_tolerance = 0.1

        # Initialize simulator
        time_scale_factor: float = OrbitSimulator.DEFAULT_TIME_SCALE_FACTOR
        sim: OrbitSimulator = OrbitSimulator(time_scale_factor=time_scale_factor)

        # Run long enough for one orbit
        time_margin: float = 0.5
        run_time: float = (sim.moon.orbit.params.period / time_scale_factor) + time_margin

        # Run it
        sim.run(run_time)

        # Check if Moon's orbit time matches given orbit time parameters
        assert sim.tracker.full_angle_times[MotionType.MOON_ORBIT] == \
            pytest.approx(sim.moon.orbit.params.period, rel=rel_tolerance)

        # Check if Moon's rotation time matches given rotation period parameters
        assert sim.tracker.full_angle_times[MotionType.MOON_ROTATION] == \
            pytest.approx(sim.moon.params.rotation_period, rel=rel_tolerance)

        # Check if the number of Earth's revolutions during one Moon's orbit matches given parameters
        assert sim.tracker.totals[MotionType.EARTH_ROTATION] == \
            pytest.approx(sim.moon.orbit.params.period_days, rel=rel_tolerance)


def main() -> None:
    """Main entry point for test."""
    parser = argparse.ArgumentParser(description='Orbit Simulator Tester')
    parser.add_argument('--no_gui', action='store_true', help='Run without GUI')
    args = parser.parse_args()

    if args.no_gui:
        OrbitSimulator.disable_gui(True)

    result = pytest.main(["tests/test_orbit_simulator.py"])

    if args.no_gui:
        sys.exit(result)
    else:
        OrbitSimulator.quit_simulation()

if __name__ == '__main__':
    main()
