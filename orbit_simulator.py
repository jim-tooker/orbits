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

To ease visualizing the simulation, there are two scaling factors that can *optionally* be used
when running this program:

- `time_scale_factor`: This factor increases the simulation's time reference vs. real-time.  
                   This allows the simulation to progress faster than reality so  
                   that observing rotations and orbits is possible.  

- `dist_scale_factor`: This factor decreases the orbital distance vs. the actual distance.  
                   This allows easier viewing of the planets. Without this scaling, planets  
                   are generally too small to view because orbital distances are relatively  
                   much larger than the planet sizes.  
"""
import os
import sys
import math
import argparse
import vpython as vp
from orbits.celestial_body import Earth, Moon
from orbits.constants import HRS_IN_DAY, SECS_IN_HR


class OrbitSimulator:
    """Main class responsible for setting up and running the orbit simulation"""

    DEFAULT_TIME_SCALE_FACTOR: float = SECS_IN_HR * HRS_IN_DAY
    """The default value for the time_scale_factor."""

    DEFAULT_DIST_SCALE_FACTOR: float = 0.1
    """The default value for the dist_scale_factor."""

    MAX_TIME_SCALE_FACTOR: float = 1_000_000
    """The max value allowed for the time_scale_factor."""

    MIN_DIST_SCALE_FACTOR: float = 0.1
    """The min value allowed for the dist_scale_factor."""

    def __init__(self,
                 time_scale_factor: float = DEFAULT_TIME_SCALE_FACTOR,
                 dist_scale_factor: float = DEFAULT_DIST_SCALE_FACTOR):
        """
        Args:
            time_scale_factor (float): How much to scale up the sense of time.
            dist_scale_factor (float): How much to scale down the orbital distances.
        """
        if not 1 <= time_scale_factor <= self.MAX_TIME_SCALE_FACTOR:
            raise ValueError(f'time_scale_factor must be between 1 and {self.MAX_TIME_SCALE_FACTOR}')
        if not self.MIN_DIST_SCALE_FACTOR <= dist_scale_factor <= 1:
            raise ValueError(f'dist_scale_factor must be between {self.MIN_DIST_SCALE_FACTOR} and 1')

        self._quit_simulation = False

        self._time_scale_factor: float = time_scale_factor
        self._dist_scale_factor: float = dist_scale_factor
        self._canvas: vp.canvas = vp.canvas(title='Orbit Simulator',
                                            width=1600,
                                            height=1000,
                                            align='left')

        self._earth = Earth()
        self._moon = Moon(dist_scale_factor)

        self._create_orientation_figure()
        self._info_canvas: vp.canvas
        self._setup_info_canvas()

        # rotate camera around the x axis to see the orbits better (not straight on)
        camera_rotate_angle = 2  # degrees
        self._canvas.camera.rotate(angle=-math.radians(camera_rotate_angle), axis=vp.vector(1, 0, 0))

    def __del__(self) -> None:
        """
        Deletes the canvases and sets the reference to None to allow canvases to disappear from GUI
        """
        if self._canvas:
            self._canvas.delete()

        if self._info_canvas:
            self._info_canvas.delete()

    def quit_simulation(self) -> None:
        """Stops the VPython server."""
        self._quit_simulation = True

    def _create_orientation_figure(self) -> None:
        """Create orientation arrows and labels to show the x, y, and z axes."""
        # Define orientation figure size and position using the largest object on the canvas.
        # Currently this is the Moon's orbit.  This will need to change when adding more orbits.
        canvas_mag: float = self._moon.orbit.orbit_mag
        orient_size: float = canvas_mag/8
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
        label_distance = orient_size * 1.1
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

    def _create_info_label(self, text: str, left_margin: int, line_number: int) -> vp.label:
        """
        Create a label with given text at specified position on the info canvas.

        Args:
            text (str): The text to display on the label.
            left_margin (int): The left margin position for the label.
            line_number (int): The vertical position (line number) for the label.

        Returns:
            vp.label: The created label object.
        """
        return vp.label(pos=vp.vector(left_margin, line_number, 0),
                                      text=text,
                                      height=16,
                                      align='left',
                                      box=False)

    def _handle_quit_button(self, button: vp.button) -> None:
        """Handles quit simulation button"""
        self.quit_simulation()

    def _stop_vp_server(self) -> None:
        """Stops the VPython server."""
        # We don't import vp_services until needed, because importing it will start
        # the server, if not started already.
        import vpython.no_notebook as vp_services  # type: ignore[import-untyped]
        vp_services.stop_server()
    
    def _setup_info_canvas(self) -> None:
        """
        Set up the information canvas with labels displaying simulation details.

        This method creates a separate canvas for displaying information about
        the simulation, including time and distance scales, and details about
        the Earth and Moon.
        """
        # Create canvas and configure
        width = 400
        height_of_quit_button = 25
        height = 1000 - height_of_quit_button
        self._info_canvas = vp.canvas(width=width, height=height, align='left')
        self._info_canvas.range = 10
        self._info_canvas.userzoom = False
        self._info_canvas.userspin = False
        self._info_canvas.userpan = False

        # Start labels 1 over from left margin (range -10 to 10)
        left_margin: int = -self._info_canvas.range + 1

        # Start lines at the top minus 1
        line_number: int = int(self._info_canvas.range * height/width - 1)

        # Create time scale label
        self._create_info_label(f'Time scale: {self._time_scale_factor:,.0f}x. 1 sec = {
              self._time_scale_factor/(SECS_IN_HR*HRS_IN_DAY):.1f} day(s).',
              left_margin, line_number)
        line_number -= 1

        # Create distance scale label
        self._create_info_label(f'Orbital distance scale: 1/{1/self._dist_scale_factor:.0f}.',
                                left_margin, line_number)
        line_number -= 2

        # Create Earth info label
        self._create_info_label(f'Earth Info:\n  Radius: {self._earth.params.radius:,.0f} km\n  Tilt: {
            self._earth.params.tilt_degrees:.1f}°\n  Sidereal day: {self._earth.sidereal_day:.2f} hrs',
                                left_margin, line_number)
        line_number -= 5

        # Create Moon info label
        self._create_info_label(f'Moon Info:\n  Radius: {self._moon.params.radius:,.0f} km\n  Tilt: {
            self._moon.params.tilt_degrees:.1f}°\n  Sidereal month: {self._moon.sidereal_month:.2f} days',
                                left_margin, line_number)
        line_number -= 5

        # Create Moon's Orbit info label
        self._create_info_label(f"Moon's Orbit Info:\n  Semi-major axis: {
            self._moon.orbit.params.semi_major_axis:,.0f} km\n  Eccentricity: {
            self._moon.orbit.params.eccentricity:.3f}\n  Inclination: {
            self._moon.orbit.params.inclination_degs:.2f}°\n  Period: {
            self._moon.orbit.params.period_days:.2f} days",
                                left_margin, line_number)
        line_number -= 6

        # Add the Quit simulation button, and bind it to _handle_quit_button()
        vp.button(pos=self._info_canvas.title_anchor,
                  text='                                Quit Simulation                                ',
                  color=vp.color.red,
                  bind=self._handle_quit_button)

        # Set default canvas back to normal canvas
        self._canvas.select()

    def run(self) -> None:
        """
        Execute the main simulation loop, updating positions of celestial bodies.

        This method runs indefinitely, updating the positions and rotations of
        the Earth and Moon based on the elapsed simulation time.
        """
        # Initialize simulation loop time variables
        t: float = 0
        dt: float = 0.01 * self._time_scale_factor

        while self._quit_simulation is False:
            vp.rate(100)

            # Update Moon's position based on total time elapsed
            self._moon.update_position(t)

            # Rotate the moon and earth based on the small time step and their respective angular velocities
            self._earth.rotate(dt)
            self._moon.rotate(dt)

            t += dt

        # Exit simulation
        self._stop_vp_server()


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
    parser.add_argument('-d', '--dist-scale-factor',
                        type=float,
                        default=OrbitSimulator.DEFAULT_DIST_SCALE_FACTOR,
                        help='How much to scale down the orbital distances.')
    args = parser.parse_args()

    try:
        simulation: OrbitSimulator = OrbitSimulator(time_scale_factor=args.time_scale_factor,
                                    dist_scale_factor=args.dist_scale_factor)
        simulation.run()
    except ValueError as e:
        print(f'Error: {e}')
        sys.exit(1)


if __name__ == '__main__':
    main()
