"""
Tests for Orbit Simulator
"""
import sys
import pytest
from orbits.orbit_simulator import OrbitSimulator


class TestOrbitSimulation:
    """
    Class for `pytest` testing of Orbit Simulator.
    """
    tolerance = 0.001


    # def test_earth_rotation_period(self, simulator: OrbitSimulator):
    #     """Verify the Earth's rotation period matches the defined value."""
    #     simulator.run()
    #     assert simulator._earth.sim_rotation_time in pytest.approx(simulator._earth.params.rotation_period, abs=1)

    # def test_moon_rotation_period(self, simulator: OrbitSimulator):
    #     """Verify the Moon's rotation period matches the defined value."""
    #     simulator.run()
    #     assert simulator._moon.sim_rotation_time in pytest.approx(simulator._moon.params.rotation_period, abs=1)

    def test_moon_orbital_period(self) -> None:
        """Verify the Moon's orbital period matches the defined value."""
        # Initialize simulator
        time_scale_factor: float = 100_000
        dist_scale_factor: float = 0.1
        simulator: OrbitSimulator = OrbitSimulator(time_scale_factor=time_scale_factor,
                                                   dist_scale_factor=dist_scale_factor)
        # Run long enough for one orbit
        time_margin: float = 1.1
        run_time: float = (simulator.moon.orbit.params.period / time_scale_factor) * time_margin

        # Run it
        simulator.run(run_time)

        assert simulator.sim_moon_orbit_time == pytest.approx(simulator.moon.orbit.params.period, rel=self.tolerance)


def main():
    #result = pytest.main(["tests/test_orbit_simulator.py"])
    result = pytest.main()
    print(result)
    OrbitSimulator._stop_vp_server()
    sys.exit(result)

if __name__ == '__main__':
    main()
