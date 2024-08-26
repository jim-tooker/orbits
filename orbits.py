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

time_scale_factor: This factor increases the simulation's time reference vs. real-time.
                   This allows the simulation to progress faster than reality so
                   that observing rotations and orbits is possible.

dist_scale_factor: This factor decreases the orbital distance vs. the actual distance.
                   This allows easier viewing of the planets. Without this scaling, planets
                   are generally too small to view because orbital distances are relatively
                   much larger than the planet sizes.
"""
import os
import sys
import math
from abc import ABC, abstractmethod
from enum import Enum, auto
import argparse
import vpython as vp

# Constants
SECS_IN_HR: int = 3600
HRS_IN_DAY: float = 23.9344696


class OrbitDirection(Enum):
    """
    Enum for orbit direction
    """
    CLOCKWISE = auto()
    COUNTER_CLOCKWISE = auto()


class Orbit(ABC):
    """
    FIXME
    """
    def __init__(self, dist_scale_factor: float = 1):
        self.dist_scale_factor: float = dist_scale_factor
        self._orbit_mag: float = 1

        self._create_path()

    @property
    @abstractmethod
    def semi_major_axis(self) -> float:
        ...

    @property
    @abstractmethod
    def eccentricity(self) -> float:
        ...

    @property
    @abstractmethod
    def inclination_degrees(self) -> float:
        ...

    @property
    @abstractmethod
    def period(self) -> float:
        ...

    @property
    @abstractmethod
    def direction(self) -> OrbitDirection:
        ...

    @property
    def inclination(self) -> float:
        return math.radians(self.inclination_degrees)

    @property
    def orbit_mag(self) -> float:
        return self._orbit_mag

    @property
    def a(self) -> float:
        """
        Semi-major axis of the Orbit
        """
        return self.semi_major_axis * self.dist_scale_factor

    @property
    def b(self) -> float:
        """
        Semi-minor axis of the Orbit
        """
        return self.a * math.sqrt(1 - self.eccentricity**2)

    @property
    def angular_velocity(self) -> float:
        av: float = 2 * math.pi / self.period
        if self.direction == OrbitDirection.COUNTER_CLOCKWISE:
            av = -av

        return av

    def calculate_next_point_on_path(self, angle: float) -> vp.vector:
        # Calculate the radial distance for this angle
        r = self.a * (1 - self.eccentricity**2) / (1 + self.eccentricity * math.cos(angle))

        # Calculate position in the x-y-z plane
        x_zero_inclination = r * math.cos(angle)  # length of max x with no inclination
        x = x_zero_inclination * math.cos(self.inclination)
        z = r * math.sin(angle)
        y = x_zero_inclination * math.sin(self.inclination)

        next_point = vp.vector(x,y,z)
        self._orbit_mag = next_point.mag

        return next_point

    def _create_path(self) -> None:
        orbit_ellipse = vp.curve(color=vp.color.gray(0.5))

        for theta in range(0, 360+1):
            theta_rad = math.radians(theta)
            next_point = self.calculate_next_point_on_path(theta_rad)
            orbit_ellipse.append(next_point)


class MoonOrbit(Orbit):
    semi_major_axis: float = 384405  # km
    eccentricity: float = 0.0549
    inclination_degrees: float = 5.145  # degrees
    period_days: float = 27.321661  # sidereal month
    period: float = period_days * HRS_IN_DAY * SECS_IN_HR
    direction: OrbitDirection = OrbitDirection.COUNTER_CLOCKWISE


class CelestialBody(ABC):
    def __init__(self, position=vp.vector(0, 0, 0)):
        self.sphere = vp.sphere(pos=position, radius=self.radius, texture=self.texture)
        self.axis = self._calculate_axis()
        self.axis_line = self._create_axis_line()

    @property
    @abstractmethod
    def radius(self) -> float:
        """
        FIXME
        """
        ...

    @property
    @abstractmethod
    def rotation_period(self) -> float:
        ...

    @property
    @abstractmethod
    def texture(self) -> object:
        ...

    @property
    @abstractmethod
    def tilt_degrees(self) -> float:
        ...

    @property
    def tilt(self) -> float:
        return math.radians(self.tilt_degrees)

    @property
    def angular_velocity(self) -> float:
        return 2 * math.pi / self.rotation_period

    def _calculate_axis(self):
        assert self.tilt
        return vp.vector(math.sin(self.tilt), math.cos(self.tilt), 0)

    def _create_axis_line(self):
        assert self.radius
        axis_length = self.radius * 3
        return vp.cylinder(pos=self.sphere.pos - axis_length/2 * self.axis,
                           axis=axis_length * self.axis,
                           radius=self.radius/50,
                           color=vp.color.white)

    def rotate(self, dt):
        self.sphere.rotate(angle=self.angular_velocity * dt, axis=self.axis, origin=self.sphere.pos)
        self.axis_line.rotate(angle=self.angular_velocity * dt, axis=self.axis, origin=self.sphere.pos)


class Earth(CelestialBody):
    radius: float = 6378  # km
    tilt_degrees: float = 23.44  # degrees
    sidereal_day: float = HRS_IN_DAY  # hours
    rotation_period: float = sidereal_day * SECS_IN_HR  # seconds
    texture: object = vp.textures.earth


class Moon(CelestialBody):
    radius: float = 1738  # km
    tilt_degrees: float  = 6.68  # degrees
    sidereal_month: float = 27.321661  # days
    rotation_period: float = sidereal_month * HRS_IN_DAY * SECS_IN_HR  # seconds
    texture: object = 'moon_texture.jpg'

    def __init__(self, dist_scale_factor: float = 1):
        self.orbit = MoonOrbit(dist_scale_factor)
        super().__init__(position=vp.vector(self.orbit.a, 0, 0))
        self.arrow = self._create_arrow()

    def _create_arrow(self):
        arrow_length = self.radius * 2
        return vp.arrow(pos=self.sphere.pos,
                        axis=arrow_length*vp.vector(1,0,0),
                        color=vp.color.yellow,
                        round=True)

    def update_position(self, t):
        '''Compute the Moon's position in the orbital plane'''
        theta = self.orbit.angular_velocity * t

        # Update the position with the rotated values
        next_point = self.orbit.calculate_next_point_on_path(theta)
        self.sphere.pos = next_point

        # Update axis line and arrow positions
        self.axis_line.pos = self.sphere.pos - self.axis_line.axis/2
        self.arrow.pos = self.sphere.pos

    def rotate(self, dt):
        super().rotate(dt)

        # Axis of arrow should be a normalized vector pointing to Earth multiplied by
        # the magnitude of the arrow's axis vector
        self.arrow.axis = vp.norm(-self.sphere.pos) * self.arrow.axis.mag


class OrbitSimulator:
    """
    This class is the primary class for setting up and running the Orbit Simulator.
    """
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

    def _create_orientation_figure(self) -> None:
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
        return vp.label(pos=vp.vector(left_margin, line_number, 0),
                                      text=text,
                                      height=16,
                                      align='left',
                                      box=False)

    def _setup_info_canvas(self) -> None:
        # Create canvas and configure
        width = 400
        height = 1000
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
        self._create_info_label(f'Earth Info:\n  Radius: {self._earth.radius:,.0f} km\n  Tilt: {
            self._earth.tilt_degrees:.1f}°\n  Sidereal day: {self._earth.sidereal_day:.2f} hrs',
                                left_margin, line_number)
        line_number -= 5

        # Create Moon info label
        self._create_info_label(f'Moon Info:\n  Radius: {self._moon.radius:,.0f} km\n  Tilt: {
            self._moon.tilt_degrees:.1f}°\n  Sidereal month: {self._moon.sidereal_month:.2f} days',
                                left_margin, line_number)
        line_number -= 5

        # Create Moon's Orbit info label
        self._create_info_label(f"Moon's Orbit Info:\n  Semi-major axis: {
            self._moon.orbit.semi_major_axis:,.0f} km\n  Eccentricity: {
            self._moon.orbit.eccentricity:.3f}\n  Inclination: {
            self._moon.orbit.inclination_degrees:.2f}°\n  Period: {self._moon.orbit.period_days:.2f} days",
                                left_margin, line_number)
        line_number -= 6

        # Set default canvas back to normal canvas
        self._canvas.select()

    def run(self):
        """
        FIXME
        """

        # Initialize simulation loop time variables
        t: float = 0
        dt: float = 0.01 * self._time_scale_factor

        while True:
            vp.rate(100)

            # Update Moon's position based on total time elapsed
            self._moon.update_position(t)

            # Rotate the moon and earth based on the small time step and their respective angular velocities
            self._earth.rotate(dt)
            self._moon.rotate(dt)

            t += dt


class CustomArgparseFormatter(argparse.ArgumentDefaultsHelpFormatter,
                              argparse.RawDescriptionHelpFormatter):
    """
    Custom class for `argparse` that combines two `formatter_class` classes.
    """
    ...


def main():
    parser = argparse.ArgumentParser(prog=os.path.basename(__file__),
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
        simulation = OrbitSimulator(time_scale_factor=args.time_scale_factor,
                                    dist_scale_factor=args.dist_scale_factor)
        simulation.run()
    except ValueError as e:
        print(f'Error: {e}')
        sys.exit(1)

if __name__ == '__main__':
    main()
