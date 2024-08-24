"""
FIXME
"""
import math
from abc import ABC, abstractmethod
from enum import Enum, auto
import vpython as vp

# Constants
SECS_IN_HR: int = 3600
HRS_IN_DAY: float = 23.9344696


class OrbitDirection(Enum):
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
    def direction(self) -> OrbitDirection:
        ...

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
    def inclination(self) -> float:
        ...

    @property
    @abstractmethod
    def period(self) -> float:
        ...

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
    inclination: float = math.radians(inclination_degrees)  # radians
    period_days: float = 27.321661  # sidereal month
    period: float = period_days * HRS_IN_DAY * SECS_IN_HR
    direction: OrbitDirection = OrbitDirection.COUNTER_CLOCKWISE

    def __init__(self, dist_scale_factor: float = 1):
        super().__init__(dist_scale_factor)


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
    def texture(self) -> float:
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

    def __init__(self):
        super().__init__()


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


class Simulation:
    DEFAULT_TIME_SCALE_FACTOR: float = SECS_IN_HR * HRS_IN_DAY
    """
    Determines how sped up the simulation is vs. real-time.
    """

    DEFAULT_DIST_SCALE_FACTOR: float = 0.1
    """
    This will reduce the orbit distance by this scale
    """

    def __init__(self,
                 time_scale_factor: float = DEFAULT_TIME_SCALE_FACTOR,
                 dist_scale_factor: float = DEFAULT_DIST_SCALE_FACTOR):
        self._time_scale_factor: float = time_scale_factor
        self._dist_scale_factor: float = dist_scale_factor

        self._canvas: vp.canvas = vp.canvas(title='Moon Orbiting Earth',
                                            width=1600,
                                            height=1000,
                                            align='left')

        self.earth = Earth()
        self.moon = Moon(dist_scale_factor)

        self._create_orientation_figure()
        self._info_canvas: vp.canvas
        self._time_scale_label: vp.label
        self._distance_scale_label: vp.label
        self._setup_info_canvas()

    def _create_orientation_figure(self) -> None:
        # Define orientation figure size and position using the largest object on the scene
        # Currently this is the Moon's orbit, be will change. FIXME
        canvas_mag: float = self.moon.orbit.orbit_mag
        orient_size: float = canvas_mag/8
        orient_pos: vp.vector = vp.vector(-canvas_mag, canvas_mag/2, 0)

        # Create the arrows in each direction of x,y,z axis
        vp.arrow(pos=orient_pos, axis=vp.vector(orient_size, 0, 0), color=vp.color.red, round=True, emissive=True)
        vp.arrow(pos=orient_pos, axis=vp.vector(0, orient_size, 0), color=vp.color.green, round=True, emissive=True)
        vp.arrow(pos=orient_pos, axis=vp.vector(0, 0, orient_size), color=vp.color.yellow, round=True, emissive=True)

        # Label the arrows
        label_distance = orient_size * 1.1
        vp.label(pos=orient_pos + vp.vector(label_distance, 0, 0), text='X', color=vp.color.red, box=False)
        vp.label(pos=orient_pos + vp.vector(0, label_distance, 0), text='Y', color=vp.color.green, box=False)
        vp.label(pos=orient_pos + vp.vector(0, 0, label_distance), text='Z', color=vp.color.yellow, box=False)

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
        self._time_scale_label = vp.label(pos=vp.vector(left_margin, line_number, 0),
                                          text=f'Time scale: {self._time_scale_factor:,.0f}x. 1 sec = {
                                               self._time_scale_factor/(SECS_IN_HR*HRS_IN_DAY):.1f} day(s).',
                                          height=16,
                                          align='left',
                                          box=False)

        line_number -= 1

        # Create distance scale label
        self._distance_scale_label = vp.label(pos=vp.vector(left_margin, line_number, 0),
                                              text=f'Orbital distance scale: 1/{1/self._dist_scale_factor:.0f}.',
                                              height=16,
                                              align='left',
                                              box=False)

        # Set default canvas back to normal canvas
        self._canvas.select()

    def run(self):
        """
        FIXME
        """
        t = 0
        dt = 0.01 * self._time_scale_factor

        while True:
            vp.rate(100)
            
            # Update moon's position based on total time elapsed
            self.moon.update_position(t)
            
            # Rotate the moon and earth based on the small time step and their respective angular velocities
            self.earth.rotate(dt)
            self.moon.rotate(dt)
            
            t += dt


def main():
    simulation = Simulation(time_scale_factor=50000, dist_scale_factor=0.1)
    simulation.run()

if __name__ == '__main__':
    main()