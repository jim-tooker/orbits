"""
FIXME
"""
from typing import Final
from enum import Enum, auto


no_gui: bool = False
"""FIXME"""

DEFAULT_TIME_SCALE_FACTOR: Final[float] = 1_000_000
"""The default value for the time_scale_factor."""

MAX_TIME_SCALE_FACTOR: Final[float] = 2_000_000
"""The max value allowed for the time_scale_factor."""

time_scale_factor: float = DEFAULT_TIME_SCALE_FACTOR

class SimMode(Enum):
    SUN_EARTH_MOON = auto()
    EARTH_MOON = auto()


sim_mode: SimMode = SimMode.SUN_EARTH_MOON