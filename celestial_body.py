"""
Celestial Body Module

This module defines classes and data structures for representing celestial bodies
in an orbital simulation. It includes abstract base classes and subclasses
for the Sun, Earth and Moon, along with their physical parameters and visualization properties.
"""
from abc import ABC
from dataclasses import dataclass, field
from copy import copy
from typing import List, Any, Final, Optional
import math
import vpython as vp
from constants import FULL_ANGLE, HRS_IN_DAY, SECS_IN_HR
import config
from orbit import Orbit

__author__ = "Jim Tooker"


@dataclass
class CelestialBodyParams:
    """
    Data class to encapsulate the key parameters of a celestial body.
    
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

    @property
    def rotation_period_hrs(self) -> float:
        """
        Convert the period from seconds to hours.

        Returns:
            float: Period in hours
        """
        return self.rotation_period / SECS_IN_HR

    @property
    def rotation_period_days(self) -> float:
        """
        Convert the period from seconds to days.

        Returns:
            float: Period in days
        """
        return self.rotation_period_hrs / HRS_IN_DAY


@dataclass
class TrailParams:
    """Parameters for the trail made by a VPython sphere."""
    trail_radius: float = 0
    trail_color: vp.vector = field(default_factory=lambda: vp.color.white)
    trail_retain: int = 1000


class CelestialBody(ABC):
    """
    Abstract base class for visualizing and animating celestial bodies
    
    Attributes:
        params (CelestialBodyParams): Parameters defining the celestial body.
        scale_factor (float): How much to scale the size of the celestial body.
    """
    def __init__(self,
                 params: CelestialBodyParams,
                 orbits: Optional[List[Orbit]] = None,
                 scale_factor: float = 1):
        """
        Args:
            params (CelestialBodyParams): Parameters defining the celestial body.
            orbits (Optional[List[orbits.orbit.Orbit]]): A list of orbits relevant to the celestial body.
            scale_factor (float): How much to scale the size of the celestial body.
        """
        self.params: CelestialBodyParams = params
        self.scale_factor = scale_factor

        # Copy the orbits, if they exist
        self._orbits: List[Orbit]
        if orbits is None:
            self._orbits = []
        else:
            self._orbits = copy(orbits)

        if config.no_gui is False:
            self._sphere: vp.sphere = vp.sphere(radius=self.radius,
                                                texture=self.params.texture,
                                                shininess=0,
                                                make_trail=True)

            self._axis: vp.vector = self._calculate_axis()
            self._axis_line: vp.cylinder = self._create_axis_line()

    @property
    def radius(self) -> float:
        """
        The radius of the celestial body.

        Returns:
            float: The radius of the celestial body.
        """
        return self.params.radius * self.scale_factor

    @property
    def orbit(self) -> Orbit:
        """
        The primary orbit of the celestial body.

        Returns:
            orbits.orbit.Orbit: The primary orbit of the celestial body.
        """
        assert self._orbits
        return self._orbits[-1]

    @property
    def position(self) -> vp.vector:
        """
        The position of the celestial body.

        Returns:
            vp.vector: The position of the celestial body.
        """
        return self._sphere.pos

    @property
    def angular_velocity(self) -> float:
        """
        Calculate the angular velocity of the celestial body.

        Returns:
            float: Angular velocity in radians per second.
        """
        return FULL_ANGLE / self.params.rotation_period

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
        axis_length: Final[float] = self.radius * 3
        return vp.cylinder(pos=self._sphere.pos - axis_length/2 * self._axis,
                           axis=axis_length * self._axis,
                           radius=self.radius/50,
                           color=vp.color.white)

    def angle(self, t: float) -> float:
        """
        Calculates the current rotational angle of the celestial body based on the given time.

        Args:
            t (float): The current time.

        Returns:
            float: Current orbit angle
        """
        return self.angular_velocity * t

    def update_position(self, t: float) -> None:
        """
        Compute position in the orbital plane based on the given time.

        Args:
            t (float): The current simulation time.
        """
        # Add up all the orbits
        new_pos: vp.vector = vp.vector(0, 0, 0)
        for orbit in self._orbits:
            new_pos += orbit.position(t)

        # Update positions of celestial body and its axis line
        self._sphere.pos = new_pos
        self._axis_line.pos = new_pos - self._axis_line.axis/2

    def rotate(self, dt: float) -> None:
        """
        Rotate the celestial body around its axis by the given time step.
        
        Args:
            dt (float): The small time increment (seconds) to rotate by.
        """
        angle: float = self.angle(dt)
        self._sphere.rotate(angle=angle, axis=self._axis, origin=self._sphere.pos)
        self._axis_line.rotate(angle=angle, axis=self._axis, origin=self._sphere.pos)

    def set_trail_params(self, params: TrailParams) -> None:
        """
        Sets the trail parameters for a VPython sphere.
        
        Args:
          params (TrailParams): The parameters for the sphere trail.
        """
        if config.no_gui is False:
            self._sphere.trail_radius = params.trail_radius
            self._sphere.trail_color = params.trail_color
            self._sphere.retain = params.trail_retain


class Sun(CelestialBody):
    """
    Represents the Sun
    
    Attributes:
        params (CelestialBodyParams): Parameters defining the Sun.
    """
    params = CelestialBodyParams(
        radius = 695_700,  # km
        tilt = 0,  # radians
        rotation_period = 27 * HRS_IN_DAY * SECS_IN_HR,  # seconds
        texture = 'images/sun_texture.jpg')

    def __init__(self,
                 scale_factor: float = 1):
        """
        Args:
            scale_factor (float): How much to scale the size of the Sun.
        """
        super().__init__(params=self.params,
                         scale_factor=scale_factor)

        if config.no_gui is False:
            # Make Sun glow
            self._sphere.emissive = True


class Earth(CelestialBody):
    """
    Representation of Earth
    
    Attributes:
        SIDEREAL_DAY (float): The sidereal day duration in hours.
        params (CelestialBodyParams): Parameters defining the Earth.
    """
    SIDEREAL_DAY: Final[float] = 23.9344696  # hours

    params = CelestialBodyParams(
        radius = (6378.137 + 6356.752) / 2,  # km
        tilt = math.radians(23.44),  # radians
        rotation_period = SIDEREAL_DAY * SECS_IN_HR,  # seconds
        texture=vp.textures.earth)

    def __init__(self,
                 orbits: Optional[List[Orbit]] = None,
                 scale_factor: float = 1):
        """
        Args:
            orbits (Optional[List[orbits.orbit.Orbit]]): A list of orbits relevant to the Earth.
            scale_factor (float): How much to scale the size of the Earth.
        """
        super().__init__(params=self.params,
                         orbits=orbits,
                         scale_factor=scale_factor)


class Moon(CelestialBody):
    """
    Represents the Earth's moon with its own orbit and visualization.
    
    Attributes:
        SIDEREAL_MONTH (float): The sidereal month duration in days.
        params (CelestialBodyParams): Parameters defining the Moon.
    """
    SIDEREAL_MONTH: Final[float] = 27.321661  # days

    params = CelestialBodyParams(
        radius = 1738,  # km
        tilt = math.radians(6.68),  # radians
        rotation_period = SIDEREAL_MONTH * HRS_IN_DAY * SECS_IN_HR,  # seconds
        texture = 'images/moon_texture.jpg')

    def __init__(self,
                 orbits: List[Orbit],
                 earth: Earth,
                 scale_factor: float = 1):
        """
        Args:
            orbits (List[orbits.orbit.Orbit]): A list of orbits relevant to the Moon.
            earth (Earth): A reference to the Earth object
            scale_factor (float): How much to scale the size of the Moon.
        """
        self.earth: Earth = earth
        super().__init__(params=self.params,
                         orbits=orbits,
                         scale_factor=scale_factor)

        if config.no_gui is False:
            # Create arrow for Moon that points to Earth
            self._arrow: vp.arrow = self._create_arrow()

    def _create_arrow(self) -> vp.arrow:
        """
        Create an arrow to indicate the Moon's orientation.

        Returns:
            vp.arrow: An arrow object pointing towards Earth.
        """
        arrow_length: Final[float] = self.radius * 2
        return vp.arrow(pos=self._sphere.pos,
                        axis=arrow_length*vp.vector(1,0,0),
                        color=vp.color.yellow,
                        round=True)

    def update_position(self, t: float) -> None:
        """
        Compute the Moon's position in the orbital plane based on the given time.

        Args:
            t (float): The current simulation time.
        """
        super().update_position(t)

        # Update arrow position
        self._arrow.pos = self._sphere.pos

    def rotate(self, dt: float) -> None:
        """
        Rotate the Moon and update its arrow orientation.

        Args:
            dt (float): The time step for rotation.
        """
        super().rotate(dt)

        # Axis of arrow should be a normalized vector pointing to Earth multiplied by
        # the magnitude of the arrow's axis vector
        self._arrow.axis = vp.norm(self.earth.position - self._sphere.pos) * self._arrow.axis.mag
