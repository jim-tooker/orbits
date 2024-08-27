"""
Celestial Body Module

This module defines classes and data structures for representing celestial bodies
in an orbital simulation. It includes abstract base classes and specific implementations
for Earth and Moon, along with their physical parameters and visualization properties.
"""
from abc import ABC
from dataclasses import dataclass
from typing import Any
import math
import vpython as vp
from orbits.orbit import MoonOrbit
from orbits.constants import HRS_IN_DAY, SECS_IN_HR


@dataclass
class CelestialBodyParams:
    """
    Data class to encapsulate the key parameters of a celestial body
    
    Attributes:
        radius (float): Radius of the celestial body in km
        tilt (float): Axial tilt of the celestial body in radians
        rotation_period (float): Rotation period of the body in seconds
        texture (Any): Texture object for the celestial body
    """
    radius: float
    tilt: float
    rotation_period: float
    texture: Any

    @property
    def tilt_degrees(self) -> float:
        """
        Convert the axial tilt from radians to degrees.

        Returns:
            float: The axial tilt in degrees.
        """
        return math.degrees(self.tilt)


class CelestialBody(ABC):
    """Abstract base class for visualizing and animating celestial bodies"""
    def __init__(self, params: CelestialBodyParams, position: vp.vector = vp.vector(0, 0, 0)):
        self.params = params
        self.sphere = vp.sphere(pos=position, radius=self.params.radius, texture=self.params.texture)
        self.axis = self._calculate_axis()
        self.axis_line = self._create_axis_line()

    @property
    def angular_velocity(self) -> float:
        """
        Calculate the angular velocity of the celestial body.

        Returns:
            float: Angular velocity in radians per second.
        """
        return 2 * math.pi / self.params.rotation_period

    def _calculate_axis(self) -> vp.vector:
        """
        Calculate the axis vector based on the body's tilt.

        Returns:
            vp.vector: The axis vector of the celestial body.
        """
        return vp.vector(math.sin(self.params.tilt), math.cos(self.params.tilt), 0)

    def _create_axis_line(self) -> vp.cylinder:
        """
        Create a visual representation of the body's axis.

        Returns:
            vp.cylinder: A cylinder object representing the axis.
        """
        axis_length = self.params.radius * 3
        return vp.cylinder(pos=self.sphere.pos - axis_length/2 * self.axis,
                           axis=axis_length * self.axis,
                           radius=self.params.radius/50,
                           color=vp.color.white)

    def rotate(self, dt: float) -> None:
        """
        Rotate the celestial body around its axis by the given time step.
        
        Args:
            dt (float): The small time increment (seconds) to rotate by.
        """
        self.sphere.rotate(angle=self.angular_velocity * dt, axis=self.axis, origin=self.sphere.pos)
        self.axis_line.rotate(angle=self.angular_velocity * dt, axis=self.axis, origin=self.sphere.pos)


class Earth(CelestialBody):
    """
    Representation of Earth
    
    Attributes:
        sidereal_day (float): The sidereal day duration in hours.
    """
    sidereal_day: float = 23.9344696  # hours

    params = CelestialBodyParams(
        radius = 6378,  # km
        tilt = math.radians(23.44),  # radians
        rotation_period = sidereal_day * SECS_IN_HR,  # seconds
        texture=vp.textures.earth
    )

    def __init__(self, position: vp.vector = vp.vector(0, 0, 0)):
        """
        Args:
            position (Optional[vp.vector]): Initial position of Earth. Defaults to (0, 0, 0).
        """
        super().__init__(params=self.params, position=position)


class Moon(CelestialBody):
    """
    Represents the Earth's moon with its own orbit and visualization.
    
    Attributes:
        sidereal_month (float): The sidereal month duration in days.
        params (CelestialBodyParams): Parameters defining the Moon's physical properties.
        orbit (MoonOrbit): Object handling the moon's orbital mechanics.
        arrow (vp.arrow): Visual indicator of moon's orientation.
    """
    sidereal_month: float = 27.321661  # days

    params = CelestialBodyParams(
        radius = 1738,  # km
        tilt = math.radians(6.68),  # radians
        rotation_period = sidereal_month * HRS_IN_DAY * SECS_IN_HR,  # seconds
        texture = 'images/moon_texture.jpg')

    def __init__(self, dist_scale_factor: float = 1):
        """
        Initialize the Moon object with optional distance scaling factor.
        
        Args:
            dist_scale_factor (float): Factor to scale down the orbital distances
        """
        self.orbit = MoonOrbit(dist_scale_factor)
        super().__init__(params=self.params, position=vp.vector(self.orbit.a, 0, 0))
        self.arrow = self._create_arrow()

    def _create_arrow(self) -> vp.arrow:
        """
        Create an arrow to indicate the Moon's orientation.

        Returns:
            vp.arrow: An arrow object pointing towards Earth.
        """
        arrow_length = self.params.radius * 2
        return vp.arrow(pos=self.sphere.pos,
                        axis=arrow_length*vp.vector(1,0,0),
                        color=vp.color.yellow,
                        round=True)

    def update_position(self, t: float) -> None:
        """
        Compute the Moon's position in the orbital plane based on the given time.

        Args:
            t (float): The current simulation time.
        """
        theta = self.orbit.angular_velocity * t

        # Update the position with the rotated values
        next_point = self.orbit.calculate_next_point_on_path(theta)
        self.sphere.pos = next_point

        # Update axis line and arrow positions
        self.axis_line.pos = self.sphere.pos - self.axis_line.axis/2
        self.arrow.pos = self.sphere.pos

    def rotate(self, dt: float) -> None:
        """
        Rotate the Moon and update its arrow orientation.

        Args:
            dt (float): The time step for rotation.
        """
        super().rotate(dt)

        # Axis of arrow should be a normalized vector pointing to Earth multiplied by
        # the magnitude of the arrow's axis vector
        self.arrow.axis = vp.norm(-self.sphere.pos) * self.arrow.axis.mag
