"""
Celestial Body Module

This module defines classes and data structures for representing celestial bodies
in an orbital simulation. It includes abstract base classes and specific implementations
for Earth and Moon, along with their physical parameters and visualization properties.
"""
from abc import ABC
from dataclasses import dataclass
from enum import Enum, auto
from copy import copy
from typing import List, Any, Final, Optional
import math
import vpython as vp
from orbits.orbit import Orbit
from orbits.constants import FULL_ANGLE, HRS_IN_DAY, SECS_IN_HR


class MotionType(Enum):
    """Enum for different types of motions we track."""
    EARTH_ROTATION = auto()
    EARTH_ORBIT = auto()
    MOON_ROTATION = auto()
    MOON_ORBIT = auto()


@dataclass
class CelestialBodyParams:
    """
    Data class to encapsulate the key parameters of a celestial body.
    
    Attributes:
        radius (float): Radius of the celestial body in km
        tilt (float): Axial tilt of the celestial body in radians
        rotation_period (float): Rotation period of the body in seconds
        texture (Any): Texture object for the celestial body
        no_gui (bool): Whether to display a GUI (True = no GUI). Defaults to False
    """
    radius: float
    tilt: float
    rotation_period: float
    texture: Any
    no_gui: bool = False

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

class CelestialBody(ABC):
    """
    Abstract base class for visualizing and animating celestial bodies
    
    Attributes:
        params (CelestialBodyParams): Parameters defining the celestial body.
        orbits (List[Orbit]): Orbit of the celestial body. FIXME because not  
                                 all bodies will have orbits.  FIXME Must be overridden in  
                                 subclass.
    """
    def __init__(self,
                 params: CelestialBodyParams,
                 orbits: Optional[List[Orbit]] = None):
        """
        Args:
            params (CelestialBodyParams): Parameters defining the celestial body.
            FIXME
            position (vp.vector): Initial position of celestial body. Defaults to (0, 0, 0).
        """
        self.params: CelestialBodyParams = params

        # Copy the orbits, if they exist
        self._orbits: List[Orbit]
        if orbits is None:
            self._orbits = []
        else:
            self._orbits = copy(orbits)

        if self.params.no_gui is False:
            self._sphere: vp.sphere = vp.sphere(radius=self.params.radius,
                                                texture=self.params.texture,
                                                make_trail=True)
            self._axis: vp.vector = self._calculate_axis()
            self._axis_line: vp.cylinder = self._create_axis_line()

    @property
    def orbit(self) -> Orbit:
        """
        The primary orbit of the celestial body.

        Returns:
            Orbit: The primary orbit of the celestial body.
        """
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

    def angle(self, t: float) -> float:
        """
        Calculates the current rotational angle of the celestial body based on the given time.

        Args:
            t (float): The current time.

        Returns:
            float: Current orbit angle
        """
        return self.angular_velocity * t

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
        axis_length: float = self.params.radius * 3
        return vp.cylinder(pos=self._sphere.pos - axis_length/2 * self._axis,
                           axis=axis_length * self._axis,
                           radius=self.params.radius/50,
                           color=vp.color.white)

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
                 no_gui: bool = False):
        """
        Args:
            no_gui (bool): Whether to display a GUI (True = no GUI). Defaults to False.
        """
        self.params.no_gui = no_gui
        super().__init__(params=self.params)

        # Make Sun glow
        self._sphere.emissive = True
        self._sphere.shininess = 1

        
class Earth(CelestialBody):
    """
    Representation of Earth
    
    Attributes:
        sidereal_day (float): The sidereal day duration in hours.
    """
    sidereal_day: Final[float] = 23.9344696  # hours

    params = CelestialBodyParams(
        radius = 200_000,  # km
        #radius = (6378.137 + 6356.752) / 2,  # km
        tilt = math.radians(23.44),  # radians
        rotation_period = sidereal_day * SECS_IN_HR,  # seconds
        texture=vp.textures.earth)

    def __init__(self,
                 orbits: List[Orbit],
                 no_gui: bool = False):
        """
        Args:
            no_gui (bool): Whether to display a GUI (True = no GUI). Defaults to False.
        """
        self.params.no_gui = no_gui
        super().__init__(params=self.params,
                         orbits=orbits)


class Moon(CelestialBody):
    """
    Represents the Earth's moon with its own orbit and visualization.
    
    Attributes:
        sidereal_month (float): The sidereal month duration in days.
        orbit (MoonOrbit): Object handling the moon's orbital mechanics.
        arrow (vp.arrow): Visual indicator of moon's orientation.
    """
    sidereal_month: Final[float] = 27.321661  # days

    params = CelestialBodyParams(
        radius = 50_000,  # km
        #radius = 1738,  # km
        tilt = math.radians(6.68),  # radians
        rotation_period = sidereal_month * HRS_IN_DAY * SECS_IN_HR,  # seconds
        texture = 'images/moon_texture.jpg')

    def __init__(self,
                 orbits: List[Orbit],
                 no_gui: bool = False):
        """
        Args:
            orbits: FIXME
            no_gui (bool): Whether to display a GUI (True = no GUI). Defaults to False.
        """
        self.params.no_gui = no_gui
        super().__init__(params=self.params,
                         orbits=orbits)

        if self.params.no_gui is False:
            # Create arrow for Moon that points to Earth
            self.arrow: vp.arrow = self._create_arrow()

    def _create_arrow(self) -> vp.arrow:
        """
        Create an arrow to indicate the Moon's orientation.

        Returns:
            vp.arrow: An arrow object pointing towards Earth.
        """
        arrow_length: float = self.params.radius * 2
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
        self.arrow.pos = self._sphere.pos

    def rotate(self, dt: float) -> None:
        """
        Rotate the Moon and update its arrow orientation.

        Args:
            dt (float): The time step for rotation.
        """
        super().rotate(dt)

        # Axis of arrow should be a normalized vector pointing to Earth multiplied by
        # the magnitude of the arrow's axis vector
        self.arrow.axis = vp.norm(-self._sphere.pos) * self.arrow.axis.mag
