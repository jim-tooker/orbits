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

The simulation runs indefinitely unless you specify a `runtime` with the call to `run()`.
"""
import os
import sys
import math
import argparse
from typing import List, Final
from dataclasses import dataclass, field
import vpython as vp
from orbits.celestial_body import Earth, Moon, MotionType, Sun
from orbits.orbit import Orbit, EarthOrbit, MoonOrbit
from orbits.constants import FULL_ANGLE, HRS_IN_DAY, SECS_IN_HR

@dataclass
class MotionTracker:
    """Data class to hold the times and angles for objects we track."""
    last_event_times: dict[MotionType, float] = field(default_factory=lambda: {body: 0.0 for body in MotionType})
    full_angle_times: dict[MotionType, float] = field(default_factory=lambda: {body: 0.0 for body in MotionType})
    angles: dict[MotionType, float] = field(default_factory=lambda: {body: 0.0 for body in MotionType})
    totals: dict[MotionType, float] = field(default_factory=lambda: {body: 0.0 for body in MotionType})


class OrbitSimulator:
    """
    Main class responsible for setting up and running the orbit simulation.

    Attributes:
        sun (Sun): Representation of the Sun.
        earth (Earth): Representation of the Earth.
        moon (Moon): Representation of the Moon.
        tracker (MotionTracker): Tracker to keep track of simulation object's angles and times.
    """

    DEFAULT_TIME_SCALE_FACTOR: Final[float] = 1_000_000
    """The default value for the time_scale_factor."""

    MAX_TIME_SCALE_FACTOR: Final[float] = 2_000_000
    """The max value allowed for the time_scale_factor."""

    MAX_RUNTIME: Final[float] = SECS_IN_HR  # 1 hr
    """The max time for a simulation to run, if a time is specified"""

   # Flag to indicate whether the GUI should be disabled (True = no GUI)
    _no_gui: bool = False

    def __init__(self,
                 time_scale_factor: float = DEFAULT_TIME_SCALE_FACTOR):
        """
        Args:
            time_scale_factor (float): How much to scale up the sense of time.
        """
        if not 1 <= time_scale_factor <= self.MAX_TIME_SCALE_FACTOR:
            raise ValueError(f'time_scale_factor must be between 1 and {self.MAX_TIME_SCALE_FACTOR}')

        self._exit_sim: bool = False
        self._time_scale_factor: float = time_scale_factor

        # Create self.tracker so we can monitor times and angles of objects
        self.tracker: MotionTracker = MotionTracker()

        if OrbitSimulator._no_gui is False:
            self._canvas: vp.canvas = vp.canvas(title='Orbit Simulator',
                                                width=1600,
                                                height=1000,
                                                align='left')

        # Create Sun
        self.sun: Sun = Sun(no_gui=OrbitSimulator._no_gui)

        # Create orbits List
        orbits: List[Orbit] = []

        # Create Earth and its orbit
        earth_orbit: EarthOrbit = EarthOrbit()
        orbits.append(earth_orbit)
        self.earth: Earth = Earth(orbits=orbits,
                                  no_gui=OrbitSimulator._no_gui)

        # Create Moon and its orbit
        moon_orbit: MoonOrbit = MoonOrbit()
        orbits.append(moon_orbit)
        self.moon: Moon = Moon(orbits=orbits,
                               earth=self.earth,
                               no_gui=OrbitSimulator._no_gui)

        if OrbitSimulator._no_gui is False:
            self._create_orientation_figure()
            self._runtime_left_label: vp.label
            self._info_canvas: vp.canvas = self._setup_info_canvas()

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

    @classmethod
    def disable_gui(cls, no_gui: bool) -> None:
        """
        Enables or disables the GUI.

        Args:
            no_gui (bool): Flag to indicate where GUI should be disabled (True = disable GUI).
        """
        cls._no_gui = no_gui

    def _exit_sim_loop(self) -> None:
        """Causes the sim loop to exit."""
        self._exit_sim = True

    def _create_orientation_figure(self) -> None:
        """Create orientation arrows and labels to show the x, y, and z axes."""
        # Define orientation figure size and position using the largest object on the canvas.
        # Currently this is the Earth's orbit around the sun.
        #   (This will need to change if a larger orbit is added)
        canvas_mag: Final[float] = self.earth.orbit.orbit_mag
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

    def _create_info_label(self,
                           text: str,
                           left_margin: int,
                           line_number: int,
                           color: vp.vector=vp.color.white) -> vp.label:
        """
        Create a label with given text at specified position on the info canvas.

        Args:
            text (str): The text to display on the label.
            left_margin (int): The left margin position for the label.
            line_number (int): The vertical position (line number) for the label.
            color (vp.vector): The color of the label text. Defaults to white.

        Returns:
            vp.label: The created label object.
        """
        return vp.label(pos=vp.vector(left_margin, line_number, 0),
                                      text=text,
                                      color=color,
                                      height=16,
                                      align='left',
                                      box=False)

    def _setup_info_canvas(self) -> vp.canvas:
        """
        Set up the information canvas with labels displaying simulation details.

        This method creates a separate canvas for displaying information about
        the simulation, including time and distance scales, and details about
        the Earth and Moon.

        Returns:
            vp.canvas: The create info canvas
        """
        # Create canvas and configure
        width: Final[int] = 400
        height: Final[int] = 1000
        info_canvas: vp.canvas = vp.canvas(width=width, height=height, align='left')
        info_canvas.range = 10
        info_canvas.userzoom = False
        info_canvas.userspin = False
        info_canvas.userpan = False

        # Start labels 1 over from left margin (range -10 to 10)
        left_margin: int = -info_canvas.range + 1

        # Start lines at the top minus 2
        line_number: int = int(info_canvas.range * height/width - 2)

        # Create time scale label
        self._create_info_label(f'Time scale: {self._time_scale_factor:,.0f}x. 1 sec = {
              self._time_scale_factor/(SECS_IN_HR*HRS_IN_DAY):.1f} day(s).',
              left_margin, line_number)
        line_number -= 2

        # Create Sun info label
        self._create_info_label('Sun Info:\n' +
                                f'  Radius: {self.sun.params.radius:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> {self.sun.scale_factor:.0f}x\n' +
                                f'  Tilt: {self.sun.params.tilt_degrees:.1f}°\n' +
                                f'  Rotation period: {self.sun.params.rotation_period_days:.2f} days',
                                left_margin, line_number)
        line_number -= 6

        # Create Earth info label
        self._create_info_label('Earth Info:\n' +
                                f'  Radius: {self.earth.params.radius:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> {self.earth.scale_factor:.0f}x\n' +
                                f'  Tilt: {self.earth.params.tilt_degrees:.1f}°\n' +
                                f'  Sidereal day: {self.earth.SIDEREAL_DAY:.2f} hrs',
                                left_margin, line_number)
        line_number -= 6

        # Create Moon info label
        self._create_info_label('Moon Info:\n' +
                                f'  Radius: {self.moon.params.radius:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> {self.moon.scale_factor:.0f}x\n' +
                                f'  Tilt: {self.moon.params.tilt_degrees:.1f}°\n' +
                                f'  Sidereal month: {self.moon.SIDEREAL_MONTH:.2f} days',
                                left_margin, line_number)
        line_number -= 6

        # Create Earth's Orbit info label
        self._create_info_label("Earth's Orbit Info:\n" +
                                f'  Semi-major axis: {self.earth.orbit.params.semi_major_axis:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> 1/{1/self.earth.orbit.scale_factor:.0f}x\n'
                                f'  Eccentricity: {self.earth.orbit.params.eccentricity:.3f}\n' +
                                f'  Inclination: {self.earth.orbit.params.inclination_degs:.2f}°\n' +
                                f'  Period: {self.earth.orbit.params.period_days:.2f} days',
                                left_margin, line_number)
        line_number -= 7

        # Create Moon's Orbit info label
        self._create_info_label("Moon's Orbit Info:\n" +
                                f'  Semi-major axis: {self.moon.orbit.params.semi_major_axis:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> 1/{1/self.moon.orbit.scale_factor:.0f}x\n'
                                f'  Eccentricity: {self.moon.orbit.params.eccentricity:.3f}\n' +
                                f'  Inclination: {self.moon.orbit.params.inclination_degs:.2f}°\n' +
                                f'  Period: {self.moon.orbit.params.period_days:.2f} days',
                                left_margin, line_number)
        line_number -= 7

        # Create Camera view info label
        self._create_info_label(text='To change view:\n' +
                                     'Rotate view: Drag with right mouse button.\n' +
                                     'Zoom: Scroll wheel or drag left/right mouse buttons.\n' +
                                     'Pan: Shift drag with left mouse button.',
                                left_margin=left_margin,
                                line_number=line_number,
                                color=vp.color.yellow)
        line_number -= 5

        # Create Sim view scale info label
        self._create_info_label(text='Note:  Celestial body sizes and orbits have\n' +
                                     '           been scaled to fit screen.\n' +
                                     'See "Simulation visual scale factor" notes above.',
                                left_margin=left_margin,
                                line_number=line_number,
                                color=vp.color.orange)
        line_number -= 4

        # Create Time remaining label
        self._runtime_left_label = self._create_info_label(text='',
                                                           left_margin=left_margin,
                                                           line_number=line_number)
        line_number -= 2

        # Set default canvas back to normal canvas
        self._canvas.select()

        return info_canvas

    def _check_for_full_angle(self, t: float) -> None:
        # If the Earth has rotated 360°, store rotation time
        if abs(self.earth.angle(t) - self.tracker.angles[MotionType.EARTH_ROTATION]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.EARTH_ROTATION] = \
                t - self.tracker.last_event_times[MotionType.EARTH_ROTATION]
            self.tracker.last_event_times[MotionType.EARTH_ROTATION] = t
            self.tracker.angles[MotionType.EARTH_ROTATION] = self.earth.angle(t)
            self.tracker.totals[MotionType.EARTH_ROTATION] += 1
            print(f'Earth Rotation {self.tracker.totals[MotionType.EARTH_ROTATION]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.EARTH_ROTATION]/(SECS_IN_HR):.2f} hours')

        # If the Moon has orbited 360°, store the orbit time
        if abs(self.moon.orbit.angle(t) - self.tracker.angles[MotionType.MOON_ORBIT]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.MOON_ORBIT] = \
                t - self.tracker.last_event_times[MotionType.MOON_ORBIT]
            self.tracker.last_event_times[MotionType.MOON_ORBIT] = t
            self.tracker.angles[MotionType.MOON_ORBIT] = self.moon.orbit.angle(t)
            self.tracker.totals[MotionType.MOON_ORBIT] += 1
            print(f'Moon Orbit {self.tracker.totals[MotionType.MOON_ORBIT]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.MOON_ORBIT]/(SECS_IN_HR*HRS_IN_DAY):.2f} days')
            print(f'There were {self.tracker.totals[MotionType.EARTH_ROTATION] / \
                self.tracker.totals[MotionType.MOON_ORBIT]:.2f} Earth rotations during the last Moon orbit.')

        # If the Moon has rotated 360°, store rotation time
        if abs(self.moon.angle(t) - self.tracker.angles[MotionType.MOON_ROTATION]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.MOON_ROTATION] = \
                t - self.tracker.last_event_times[MotionType.MOON_ROTATION]
            self.tracker.last_event_times[MotionType.MOON_ROTATION] = t
            self.tracker.angles[MotionType.MOON_ROTATION] = self.moon.angle(t)
            self.tracker.totals[MotionType.MOON_ROTATION] += 1
            print(f'Moon Rotation {self.tracker.totals[MotionType.MOON_ROTATION]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.MOON_ROTATION]/(SECS_IN_HR*HRS_IN_DAY):.2f} days')

        # If the Earth has orbited 360° around the Sun, store the orbit time
        if abs(self.earth.orbit.angle(t) - self.tracker.angles[MotionType.EARTH_ORBIT]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.EARTH_ORBIT] = \
                t - self.tracker.last_event_times[MotionType.EARTH_ORBIT]
            self.tracker.last_event_times[MotionType.EARTH_ORBIT] = t
            self.tracker.angles[MotionType.EARTH_ORBIT] = self.earth.orbit.angle(t)
            self.tracker.totals[MotionType.EARTH_ORBIT] += 1
            print(f'Earth Orbit {self.tracker.totals[MotionType.EARTH_ORBIT]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.EARTH_ORBIT]/(SECS_IN_HR*HRS_IN_DAY):.2f} days')

        # If the Sun has rotated 360°, store rotation time
        if abs(self.sun.angle(t) - self.tracker.angles[MotionType.SUN_ROTATION]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.SUN_ROTATION] = \
                t - self.tracker.last_event_times[MotionType.SUN_ROTATION]
            self.tracker.last_event_times[MotionType.SUN_ROTATION] = t
            self.tracker.angles[MotionType.SUN_ROTATION] = self.sun.angle(t)
            self.tracker.totals[MotionType.SUN_ROTATION] += 1
            print(f'Sun Rotation {self.tracker.totals[MotionType.SUN_ROTATION]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.SUN_ROTATION]/(SECS_IN_HR*HRS_IN_DAY):.2f} days')

    def _handle_quit_button(self, button: vp.button) -> None:
        """Handles quit simulation button"""
        self._exit_sim_loop()

    @staticmethod
    def quit_simulation() -> None:
        """Stops the VPython server."""
        if OrbitSimulator._no_gui is False:
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
        dt_prime: Final[float] = dt * self._time_scale_factor

        # If GUI enabled and this is an indefinite runtime, add a quit button
        if OrbitSimulator._no_gui is False and runtime == 0:
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

            if OrbitSimulator._no_gui is False:
                # Update celestial bodies orbit position based on total time elapsed
                self.earth.update_position(t_prime)
                self.moon.update_position(t_prime)

                # Rotate celestial bodies based on the small time step and their respective angular velocities
                self.sun.rotate(dt_prime)
                self.earth.rotate(dt_prime)
                self.moon.rotate(dt_prime)

                # If we're on a timed runtime, display how much time we have left
                if runtime != 0:
                    self._runtime_left_label.text = f'Simulation time left: {(runtime - t):.1f}'


            # Check the objects we're tracking for full angles
            self._check_for_full_angle(t_prime)

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
                        default=int(OrbitSimulator.DEFAULT_TIME_SCALE_FACTOR),
                        help='How much to scale up the sense of time.')
    parser.add_argument('-r', '--runtime',
                        type=float,
                        default=0,
                        help='How many secs to run the simulation.')
    parser.add_argument('--no_gui', action='store_true', help='Run without GUI')
    args = parser.parse_args()

    if args.no_gui is True:
        OrbitSimulator.disable_gui(True)
        print('Hit Ctrl-C to exit.')

    try:
        simulation: OrbitSimulator = OrbitSimulator(time_scale_factor=args.time_scale_factor)
        simulation.run(args.runtime)
    except ValueError as e:
        print(f'Error: {e}')
        sys.exit(1)

    # Quit simulation
    simulation.quit_simulation()


if __name__ == '__main__':
    main()
