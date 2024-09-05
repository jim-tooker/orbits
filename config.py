"""
This file holds all the configuration items for the Orbit Simulator.
"""
from typing import Final
from enum import Enum, auto

no_gui: bool = False
"""Flag to enable or disable the GUI"""

DEFAULT_TIME_SCALE_FACTOR: Final[float] = 1_000_000
"""The default value for the time_scale_factor."""

MAX_TIME_SCALE_FACTOR: Final[float] = 2_000_000
"""The max value allowed for the time_scale_factor."""

time_scale_factor: float = DEFAULT_TIME_SCALE_FACTOR
"""
This factor increases the simulation's time reference vs. real-time.
This allows the simulation to progress faster than reality so that
observing rotations and orbits is possible.
"""

class SimMode(Enum):
    """
    Enum for different Simulation modes.
    """
    SUN_EARTH_MOON = auto()
    EARTH_MOON = auto()

sim_mode: SimMode = SimMode.SUN_EARTH_MOON
"""
What simulation mode is selected.  Defaults to SUN_EARTH_MOON
"""
