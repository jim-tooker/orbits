import math
from abc import ABC, abstractmethod
from typing import List
from enum import Enum, auto
import vpython as vp

TIME_SCALE_FACTOR = 100000
DIST_SCALE_FACTOR = 0.1  # This will reduce the orbit distance by this scale
SECS_IN_HR = 3600
HRS_IN_DAY = 23.9344696

class OrbitDirection(Enum):
    CLOCKWISE = auto()
    COUNTER_CLOCKWISE = auto()

class Orbit(ABC):
    def __init__(self):
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
    def a(self) -> float:
        """
        Semi-major axis of the Orbit
        """
        return self.semi_major_axis * DIST_SCALE_FACTOR

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
    
    def _calculate_next_point_on_path(self, angle: float) -> vp.vector:
        # Calculate the radial distance for this angle
        r = self.a * (1 - self.eccentricity**2) / (1 + self.eccentricity * math.cos(angle))

        # Calculate position in the x-y-z plane
        c = r * math.cos(angle)  # length of max x with no inclination
        x = c * math.cos(self.inclination)
        z = r * math.sin(angle)
        y = c * math.sin(self.inclination)

        return vp.vector(x, y, z)


    def _create_path(self) -> None:
        orbit_ellipse = vp.curve(color=vp.color.red)

        for theta in range(0, 360+1):
            theta_rad = math.radians(theta)
            next_point = self._calculate_next_point_on_path(theta_rad)
            orbit_ellipse.append(next_point)

class MoonOrbit(Orbit):
    semi_major_axis: float = 384405  # km
    eccentricity: float = 0.0549
    inclination_degrees: float = 5.145  # degrees
    inclination: float = math.radians(inclination_degrees)  # radians
    period_days: float = 27.321661  # sidereal month
    period: float = period_days * HRS_IN_DAY * SECS_IN_HR
    direction: OrbitDirection = OrbitDirection.COUNTER_CLOCKWISE

    def __init__(self):
        super().__init__()

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
    texture: object = vp.textures.rock

    def __init__(self):
        self.orbit = MoonOrbit()
        super().__init__(position=vp.vector(self.orbit.a, 0, 0))
        self.arrow = self._create_arrow()

    def _create_arrow(self):
        arrow_length = self.radius * 2
        return vp.arrow(pos=self.sphere.pos,
                        axis=arrow_length*vp.vector(1,0,0),
                        color=vp.color.yellow,
                        round=True)

    def update_position(self, t):
        # Compute the Moon's position in the orbital plane
        theta = self.orbit.angular_velocity * t

        # Update the position with the rotated values
        next_point = self.orbit._calculate_next_point_on_path(theta)
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
    def __init__(self):
        self.scene = self._setup_scene()
        self.earth = Earth()
        self.moon = Moon()
        self.orientation = self._create_orientation()

    def _setup_scene(self):
        return vp.canvas(title="Moon Orbiting Earth", width=1600, height=1000, background=vp.color.black)

    def _create_orientation(self):
        orient_size = self.earth.radius
        orient_ps = vp.vector(-self.moon.orbit.a, self.moon.orbit.a/2, 0)  #FIXME
        orient_object = vp.compound([
            vp.arrow(pos=orient_ps, axis=vp.vector(orient_size, 0, 0), color=vp.color.red, round=True),
            vp.arrow(pos=orient_ps, axis=vp.vector(0, orient_size, 0), color=vp.color.green, round=True),
            vp.arrow(pos=orient_ps, axis=vp.vector(0, 0, orient_size), color=vp.color.blue, round=True)
        ])
        label_distance = orient_size * 1.2
        vp.label(pos=orient_ps + vp.vector(label_distance, 0, 0), text="X", color=vp.color.red, box=False)
        vp.label(pos=orient_ps + vp.vector(0, label_distance, 0), text="Y", color=vp.color.green, box=False)
        vp.label(pos=orient_ps + vp.vector(0, 0, label_distance), text="Z", color=vp.color.blue, box=False)
        return orient_object

    def run(self):
        t = 0
        dt = 0.01  # Small time step for smooth animation

        while True:
            vp.rate(100)
            
            # Update moon's position based on total time elapsed
            self.moon.update_position(t * TIME_SCALE_FACTOR)
            
            # Rotate the moon and earth based on the small time step and their respective angular velocities
            self.earth.rotate(dt * TIME_SCALE_FACTOR)
            self.moon.rotate(dt * TIME_SCALE_FACTOR)
            
            t += dt


def main():
    simulation = Simulation()
    simulation.run()

if __name__ == "__main__":
    main()