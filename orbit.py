"""
Orbit Module

This module defines classes and data structures for representing orbits
in the celestial body simulation. It includes abstract base classes and
specific implementations for orbits, along with their associated parameters
and visualization elements.
"""
from enum import Enum, auto
from abc import ABC
from dataclasses import dataclass
import math
import vpython as vp
from orbits.constants import HRS_IN_DAY, SECS_IN_HR

class OrbitDirection(Enum):
    """Enum to represent the direction of an orbit (clockwise or counter-clockwise)."""
    CLOCKWISE = auto()
    COUNTER_CLOCKWISE = auto()


@dataclass
class OrbitParams:
    """
    Data class to encapsulate the key parameters of a celestial body's orbit
    
    Attributes:
        semi_major_axis (float): The semi-major axis of the orbit in km
        eccentricity (float): The eccentricity of the orbit (0 <= e < 1)
        inclination (float): The inclination of the orbit in radians
        period (float): The orbital period in seconds
        direction (OrbitDirection): The direction of the orbit (clockwise or counter-clockwise)
        no_gui (bool): Whether to display a GUI (True = no GUI). Defaults to False
"""
    semi_major_axis: float
    eccentricity: float
    inclination: float
    period: float
    direction: OrbitDirection
    no_gui: bool = False

    @property
    def inclination_degs(self) -> float:
        """
        Convert the inclination from radians to degrees.

        Returns:
            float: Inclination in degrees
        """
        return math.degrees(self.inclination)

    @property
    def period_hrs(self) -> float:
        """
        Convert the period from seconds to hours.

        Returns:
            float: Period in hours
        """
        return self.period / SECS_IN_HR

    @property
    def period_days(self) -> float:
        """
        Convert the period from seconds to days.

        Returns:
            float: Period in days
        """
        return self.period_hrs / HRS_IN_DAY


class Orbit(ABC):
    """
    Abstract base class that defines common properties and methods for all orbits.
    
    Attributes:
        params (OrbitParams): Parameters defining the orbit.
    """
    def __init__(self, params: OrbitParams, dist_scale_factor: float = 1):
        """
        Args:
            params (OrbitParams): Parameters defining the orbit
            dist_scale_factor (float): Scaling factor for orbital distance. Defaults to 1.
        """
        self.params: OrbitParams = params
        self.__dist_scale_factor: float = dist_scale_factor
        self.__orbit_mag: float = 1

        if self.params.no_gui is False:
            self.__create_path()

    @property
    def a(self) -> float:
        """
        Calculate the scaled semi-major axis of the orbit.

        Returns:
            float: Scaled semi-major axis
        """
        return self.params.semi_major_axis * self.__dist_scale_factor

    @property
    def b(self) -> float:
        """
        Calculate the scaled semi-minor axis of the orbit.

        Returns:
            float: Scaled semi-minor axis
        """
        return self.a * math.sqrt(1 - self.params.eccentricity**2)

    @property
    def angular_velocity(self) -> float:
        """
        Calculate the angular velocity of the orbit.

        Returns:
            float: Angular velocity in radians per second
        """
        av: float = 2 * math.pi / self.params.period
        if self.params.direction == OrbitDirection.COUNTER_CLOCKWISE:
            av = -av

        return av

    @property
    def orbit_mag(self) -> float:
        """
        Get the current magnitude of the orbit.

        Returns:
            float: Orbit magnitude
        """
        return self.__orbit_mag

    def angle(self, t: float) -> float:
        """
        Calculates the current angle of the orbit based on the given time.

        Args:
            t (float): The current time.

        Returns:
            float: Current orbit angle
        """
        return self.angular_velocity * t

    def _calculate_next_point_on_path(self, angle: float) -> vp.vector:
        """
        Calculate the next point on the orbit path at the given angle.
        
        Args:
            angle (float): The current angle in radians along the orbit
            
        Returns:
            vp.vector: The 3D position of the body on the orbit
        """
        # Calculate the radial distance for this angle
        r: float = self.a * (1 - self.params.eccentricity**2) / (1 + self.params.eccentricity * math.cos(angle))

        # Calculate position in the x-y-z plane
        x_zero_inclination: float = r * math.cos(angle)  # length of max x with no inclination
        x: float = x_zero_inclination * math.cos(self.params.inclination)
        z: float = r * math.sin(angle)
        y: float = x_zero_inclination * math.sin(self.params.inclination)

        next_point: vp.vector = vp.vector(x,y,z)
        self.__orbit_mag = next_point.mag

        return next_point

    def __create_path(self) -> None:
        """Create the visual representation of the orbit path."""
        orbit_ellipse: vp.curve = vp.curve(color=vp.color.gray(0.5))

        for theta in range(0, 360+1):
            theta_rad: float = math.radians(theta)
            next_point: vp.vector = self._calculate_next_point_on_path(theta_rad)
            orbit_ellipse.append(next_point)

    def update_position(self, t: float) -> vp.vector:
        """
        Compute position in the orbital plane based on the given time.

        Args:
            t (float): The current simulation time.
        """
        return self._calculate_next_point_on_path(self.angle(t))


class MoonOrbit(Orbit):
    """
    Represents the orbit of the Moon around Earth.

    This class inherits from the Orbit base class and defines
    specific parameters for the Moon's orbit.
    """
    sidereal_month: float = 27.321661  # days

    params = OrbitParams(
        semi_major_axis = 384405,  # km
        eccentricity = 0.0549,
        inclination = math.radians(5.145),  # radians
        period = sidereal_month * HRS_IN_DAY * SECS_IN_HR,  # secs
        direction = OrbitDirection.COUNTER_CLOCKWISE)

    def __init__(self, dist_scale_factor: float = 1, no_gui: bool = False):
        """
        Args:
            dist_scale_factor (float): Scaling factor for orbital distance. Defaults to 1.
        """
        self.params.no_gui = no_gui
        super().__init__(params=self.params, dist_scale_factor=dist_scale_factor)
