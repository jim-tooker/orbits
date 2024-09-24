#!/usr/bin/env python
"""
Welcome to the Orbit Simulator.

This program simulates the rotation and orbital physics of celestial bodies.  For each celestial
body, the following are visualized:  
- Rotation period around its axis  
- Axial tilt with respect to its orbit  
- Orbital inclination  
- Orbital period  
- Orbital radius  
- Path of the orbit, including direction and semi-major, semi-minor axis, and eccentricity  
"""
import os
import sys
import argparse
from typing import Final
import vpython as vp
from constants import SECS_IN_HR
import config
from simulation_mode import SimulationMode, SunEarthMoonMode, EarthMoonMode

__author__ = "Jim Tooker"


class OrbitSimulator:
    """
    Main class responsible for setting up and running the orbit simulation.

    Attributes:
      mode (orbits.simulation_mode.SimulationMode): Reference to the simulation mode instance
    """

    @staticmethod
    def quit_simulation() -> None:
        """Stops the VPython server."""
        if config.no_gui is False:
            # We don't import vp_services until needed, because importing it will start
            # the server, if not started already.
            import vpython.no_notebook as vp_services  # type: ignore[import-untyped]
            vp_services.stop_server()

    def __init__(self) -> None:
        if not 0 <= config.time_scale_factor <= config.MAX_TIME_SCALE_FACTOR:
            raise ValueError(f'time_scale_factor must be between 0 and {config.MAX_TIME_SCALE_FACTOR}')

        self.mode: SimulationMode
        if config.sim_mode == config.SimMode.EARTH_MOON:
            self.mode = EarthMoonMode()
        else:
            self.mode = SunEarthMoonMode()

    def run(self, runtime: float = 0) -> None:
        """
        Execute the main simulation loop, updating positions of celestial bodies.

        This method runs as long as `runtime` specifies, or indefinitely if `runtime=0`.
        It updates the positions and rotations of the celestial bodies based on the elapsed simulation time.

        Args:
            runtime (float): How long to run the simulation (s).  Defaults to 0 (indefinite runtime)
        """
        # The max time for a simulation to run, if a time is specified
        max_runtime: Final[float] = 1 * SECS_IN_HR  # 1 hr

        if not 0 <= runtime <= max_runtime:
            raise ValueError(f'Runtime must be between 0 and {max_runtime}')

        # Initialize simulation loop time variables
        t: float = 0
        t_prime: float = 0
        dt: Final[float] = 0.01
        dt_prime: Final[float] = dt * config.time_scale_factor

        # If GUI enabled and this is an indefinite runtime, add a quit button
        if config.no_gui is False and runtime == 0:
            self.mode.add_quit_button()

        print()

        while self.mode.quit_sim is False and True if runtime == 0 else t < runtime:
            vp.rate(100)

            if config.no_gui is False:
                self.mode.update_celestial_bodies(t_prime, dt_prime)

                # If we're on a timed runtime, display how much time we have left
                if runtime != 0:
                    self.mode.runtime_left_label.text = f'Simulation time left: {(runtime - t):.1f}'


            # Check the objects we're tracking for full angles
            self.mode.check_for_full_angle(t_prime)

            t += dt
            t_prime += dt_prime


class CustomArgparseFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawDescriptionHelpFormatter):
    """Custom class for `argparse` that combines two `formatter_class` classes."""


def main() -> None:
    """
    Parse command-line arguments and run the orbit simulation.

    This function creates an OrbitSimulator instance with the specified
    parameters, and runs the simulation.

    Raises:
        SystemExit: If there's a ValueError during OrbitSimulator initialization.
    """
    parser: argparse.ArgumentParser = argparse.ArgumentParser(prog=os.path.basename(__file__),
                                     formatter_class=CustomArgparseFormatter,
                                     description=__doc__)
    parser.add_argument('-t', '--time_scale_factor',
                        type=float,
                        default=0,
                        help='How much to scale up the sense of time. 0 means use best-fit time scaling.')
    parser.add_argument('-r', '--runtime',
                        type=float,
                        default=0,
                        help='How many secs to run the simulation. 0 means run indefinitely.')
    parser.add_argument('-m', '--mode',
                        choices=['sun_earth_moon', 'earth_moon'],
                        default='sun_earth_moon',
                        help='Simulation mode to use.')
    parser.add_argument('--no_gui', action='store_true', help='Run without GUI')
    args = parser.parse_args()

    if args.no_gui is True:
        config.no_gui = True
        print('Hit Ctrl-C to exit.')

    if args.mode == 'earth_moon':
        config.sim_mode = config.SimMode.EARTH_MOON

    if args.time_scale_factor != 0:
        config.time_scale_factor = args.time_scale_factor

    try:
        simulation: OrbitSimulator = OrbitSimulator()
        simulation.run(args.runtime)
    except ValueError as e:
        print(f'Error: {e}')
        sys.exit(1)

    # Quit simulation
    simulation.quit_simulation()


if __name__ == '__main__':
    main()
