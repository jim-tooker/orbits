"""
This file holds all the configuration items for the Orbit Simulator.
"""
from typing import Final
from enum import Enum, auto
from orbits.constants import SECS_IN_HR, HRS_IN_DAY

no_gui: bool = False
"""Flag to enable or disable the GUI"""

MAX_TIME_SCALE_FACTOR: Final[float] = 2_000_000
"""The max value allowed for the time_scale_factor."""

time_scale_factor: float = 0
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
