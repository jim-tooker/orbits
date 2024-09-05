"""
Orbit Module

This module defines classes and data structures for representing orbits
in the celestial body simulation. It includes abstract base classes and
specific implementations for orbits, along with their associated parameters
and visualization elements.
"""
from __future__ import annotations
from enum import Enum, auto
from abc import ABC
from typing import Final
from dataclasses import dataclass
import math
import vpython as vp
from orbits.constants import FULL_ANGLE, HRS_IN_DAY, SECS_IN_HR

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
"""
    semi_major_axis: float
    eccentricity: float
    inclination: float
    period: float
    direction: OrbitDirection

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
        scale_factor (float): How much to scale the size of orbit.
    """
    def __init__(self,
                 params: OrbitParams,
                 scale_factor: float = 1):
        """
        Args:
            params (OrbitParams): Parameters defining the orbit.
            scale_factor (float): How much to scale the size of orbit.
        """
        self.params: OrbitParams = params
        self.scale_factor: float = scale_factor

    @property
    def a(self) -> float:
        """
        Calculate the scaled semi-major axis of the orbit.

        Returns:
            float: Scaled semi-major axis
        """
        return self.params.semi_major_axis * self.scale_factor

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
        av: float = FULL_ANGLE / self.params.period
        if self.params.direction == OrbitDirection.COUNTER_CLOCKWISE:
            av = -av

        return av

    @property
    def orbit_mag(self) -> float:
        """
        The magnitude of the orbit.

        Returns:
            float: The magnitude of the orbit.
        """
        return self.calculate_next_point_on_path(0).mag

    def angle(self, t: float) -> float:
        """
        Calculates the current angle of the orbit based on the given time.

        Args:
            t (float): The current time.

        Returns:
            float: Current orbit angle
        """
        return self.angular_velocity * t

    def calculate_next_point_on_path(self, angle: float) -> vp.vector:
        """
        Calculate the next point on the orbit path at the given angle.
        
        Args:
            angle (float): The current angle in radians along the orbit
            
        Returns:
            vp.vector: The 3D position of the body on the orbit
        """
        # Shift the angle by pi to start with the orbit coming towards the user at input angle=0
        angle = angle + math.pi

        # Calculate the radial distance for this angle
        r: float = self.a * (1 - self.params.eccentricity**2) / (1 + self.params.eccentricity * math.cos(angle))

        # Calculate position in the x-y-z plane
        x_zero_inclination: float = r * math.cos(angle)  # length of max x with no inclination
        x: float = x_zero_inclination * math.cos(self.params.inclination)
        z: float = r * math.sin(angle)
        y: float = x_zero_inclination * math.sin(self.params.inclination)

        next_point: vp.vector = vp.vector(x,y,z)

        return next_point

    def position(self, t: float) -> vp.vector:
        """
        Compute position in the orbital plane based on the given time.

        Args:
            t (float): The current simulation time.

        Returns:
            vp.vector: The position at time t.
        """
        return self.calculate_next_point_on_path(self.angle(t))


class EarthOrbit(Orbit):
    """
    Represents the orbit of the Earth around the Sun.

    This class inherits from the Orbit base class and defines
    specific parameters for the Earth's orbit.
    
    Attributes:
        SIDEREAL_YEAR (float): The sidereal year duration in days.
        SCALE_FACTOR (float): How much to scale the size of the Earth orbit.
        params (OrbitParams): Parameters defining the Earth's orbit.
    """
    SIDEREAL_YEAR: Final[float] = 365.256 # days
    SCALE_FACTOR: Final[float] = 1/30

    params = OrbitParams(
        semi_major_axis = 149_597_870,  # km
        eccentricity = 0.0167,
        inclination = math.radians(0),  # radians
        period = SIDEREAL_YEAR * HRS_IN_DAY * SECS_IN_HR,  # secs
        direction = OrbitDirection.COUNTER_CLOCKWISE)

    def __init__(self) -> None:
        """
        """
        super().__init__(params=self.params,
                         scale_factor=self.SCALE_FACTOR)


class MoonOrbit(Orbit):
    """
    Represents the orbit of the Moon around Earth.

    This class inherits from the Orbit base class and defines
    specific parameters for the Moon's orbit.

    Attributes:
        SIDEREAL_MONTH (float): The sidereal month duration in days.
        SCALE_FACTOR (float): How much to scale the size of the Moon's orbit.
        params (OrbitParams): Parameters defining the Moon's orbit.
    """
    SIDEREAL_MONTH: Final[float] = 27.321661  # days
    SCALE_FACTOR: Final[float] = 1/4

    params = OrbitParams(
        semi_major_axis = 384_405,  # km
        eccentricity = 0.0549,
        inclination = math.radians(5.145),  # radians
        period = SIDEREAL_MONTH * HRS_IN_DAY * SECS_IN_HR,  # secs
        direction = OrbitDirection.COUNTER_CLOCKWISE)

    def __init__(self) -> None:
        """
        """
        super().__init__(params=self.params,
                         scale_factor=self.SCALE_FACTOR)
