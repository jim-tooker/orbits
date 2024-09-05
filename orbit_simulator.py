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

To ease visualizing the simulation, there is a time scaling factor that can *optionally* be used
when running this program:

- `time_scale_factor`: This factor increases the simulation's time reference vs. real-time.
                       This allows the simulation to progress faster than reality so
                       that observing rotations and orbits is possible.
- FIXME (add others)

The simulation runs indefinitely unless you specify a `runtime` with the call to `run()`.
"""
import os
import sys
import math
import argparse
from typing import Final
import vpython as vp
from orbits.constants import SECS_IN_HR
from orbits.simulation_mode import SunEarthMoonMode, EarthMoonMode
from orbits import config
from orbits.config import SimMode

class OrbitSimulator:
    """
    Main class responsible for setting up and running the orbit simulation.

    Attributes:
        sun (Sun): Representation of the Sun.
        earth (Earth): Representation of the Earth.
        moon (Moon): Representation of the Moon.
        tracker (MotionTracker): Tracker to keep track of simulation object's angles and times.
    """

    MAX_RUNTIME: Final[float] = SECS_IN_HR  # 1 hr
    """The max time for a simulation to run, if a time is specified"""

    def __init__(self):
        if not 1 <= config.time_scale_factor <= config.MAX_TIME_SCALE_FACTOR:
            raise ValueError(f'time_scale_factor must be between 1 and {config.MAX_TIME_SCALE_FACTOR}')

        if config.sim_mode == SimMode.EARTH_MOON:
            self.mode = EarthMoonMode()
        else:
            self.mode = SunEarthMoonMode()

        if config.no_gui is False:
            self._canvas: vp.canvas = vp.canvas(title='Orbit Simulator',
                                                width=1600,
                                                height=1000,
                                                align='left')

        self.mode.create_celestial_bodies()

        self._exit_sim: bool = False

        if config.no_gui is False:
            self._create_orientation_figure()
            self._runtime_left_label: vp.label
            self._info_canvas: vp.canvas = self.mode.setup_info_canvas()

            # Set default canvas back to normal canvas
            self._canvas.select()

            # rotate camera around the x axis to see the orbits better (not straight on)
            camera_rotate_angle: Final[float] = 30  # degrees
            self._canvas.camera.rotate(angle=-math.radians(camera_rotate_angle), axis=vp.vector(1, 0, 0))

    def __del__(self) -> None:
        """
        Deletes the canvases and sets the reference to None to allow canvases to disappear from GUI
        """
        if hasattr(self, '_canvas') and self._canvas:
            self._canvas.delete()

        if hasattr(self, '_info_canvas') and self._info_canvas:
            self._info_canvas.delete()

    def _exit_sim_loop(self) -> None:
        """Causes the sim loop to exit."""
        self._exit_sim = True

    def _create_orientation_figure(self) -> None:
        """Create orientation arrows and labels to show the x, y, and z axes."""
        # Define orientation figure size and position using the largest object on the canvas.
        # Currently this is the Earth's orbit around the sun.
        #   (This will need to change if a larger orbit is added)
        canvas_mag: Final[float] = self.mode.earth.orbit.orbit_mag     #FIXME (get the biggest??)
        orient_size: Final[float] = canvas_mag/8
        orient_pos: vp.vector = vp.vector(-canvas_mag, canvas_mag/2, 0)

        # Create the arrows in each direction of x,y,z axis
        vp.arrow(pos=orient_pos,
                 axis=vp.vector(orient_size, 0, 0),
                 color=vp.color.red,
                 round=True,
                 emissive=True)
        vp.arrow(pos=orient_pos,
                 axis=vp.vector(0, orient_size, 0),
                 color=vp.color.green,
                 round=True,
                 emissive=True)
        vp.arrow(pos=orient_pos,
                 axis=vp.vector(0, 0, orient_size),
                 color=vp.color.yellow,
                 round=True,
                 emissive=True)

        # Label the arrows
        label_distance: float = orient_size * 1.1
        vp.label(pos=orient_pos + vp.vector(label_distance, 0, 0),
                 text='X',
                 color=vp.color.red,
                 opacity=0,
                 box=False)
        vp.label(pos=orient_pos + vp.vector(0, label_distance, 0),
                 text='Y',
                 color=vp.color.green,
                 opacity=0,
                 box=False)
        vp.label(pos=orient_pos + vp.vector(0, 0, label_distance),
                 text='Z',
                 color=vp.color.yellow,
                 opacity=0,
                 box=False)

    def _handle_quit_button(self, button: vp.button) -> None:
        """Handles quit simulation button"""
        self._exit_sim_loop()

    @staticmethod
    def quit_simulation() -> None:
        """Stops the VPython server."""
        if config.no_gui is False:
            # We don't import vp_services until needed, because importing it will start
            # the server, if not started already.
            import vpython.no_notebook as vp_services  # type: ignore[import-untyped]
            vp_services.stop_server()

    def run(self, runtime: float = 0) -> None:
        """
        Execute the main simulation loop, updating positions of celestial bodies.

        This method runs as long as `runtime` specifies, or indefinitely if `runtime=0`.
        It updates the positions and rotations of the celestial bodies based on the elapsed simulation time.

        Args:
            runtime (float): How long to run the simulation (s).  Defaults to 0 (indefinite runtime)
        """
        if not 0 <= runtime <= self.MAX_RUNTIME:
            raise ValueError(f'Runtime must be between 0 and {self.MAX_RUNTIME}')

        # Initialize simulation loop time variables
        t: float = 0
        t_prime: float = 0
        dt: Final[float] = 0.01
        dt_prime: Final[float] = dt * config.time_scale_factor

        # If GUI enabled and this is an indefinite runtime, add a quit button
        if config.no_gui is False and runtime == 0:
            # Reduce the height of the info canvas by the height of the quit button
            height_of_quit_button: Final[int] = 25
            self._info_canvas.height -= height_of_quit_button

            # Add the Quit simulation button, and bind it to _handle_quit_button()
            vp.button(pos=self._info_canvas.title_anchor,
                    text='                                Quit Simulation                                ',
                    color=vp.color.red,
                    background=vp.color.gray(0.8),
                    bind=self._handle_quit_button)

        print()

        while self._exit_sim is False and True if runtime == 0 else t < runtime:
            vp.rate(100)

            if config.no_gui is False:
                self.mode.update_celestial_bodies(t_prime, dt_prime)

                # If we're on a timed runtime, display how much time we have left
                if runtime != 0:
                    self._runtime_left_label.text = f'Simulation time left: {(runtime - t):.1f}'


            # Check the objects we're tracking for full angles
            self.mode.check_for_full_angle(t_prime)

            t += dt
            t_prime += dt_prime


class CustomArgparseFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawDescriptionHelpFormatter):
    """Custom class for `argparse` that combines two `formatter_class` classes."""
    ...


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
    parser.add_argument('-t', '--time-scale-factor',
                        type=float,
                        default=0,
                        help='How much to scale up the sense of time.')
    parser.add_argument('-r', '--runtime',
                        type=float,
                        default=0,
                        help='How many secs to run the simulation.')
    parser.add_argument('-m', '--mode',
                        choices=['sun_earth_moon', 'earth_moon'],
                        default='sun_earth_moon',
                        help='Simulation mode to use')
    parser.add_argument('--no_gui', action='store_true', help='Run without GUI')
    args = parser.parse_args()

    if args.no_gui is True:
        config.no_gui = True
        print('Hit Ctrl-C to exit.')

    if args.mode == 'earth_moon':
        config.sim_mode = SimMode.EARTH_MOON

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
