"""
FIXME
"""
from dataclasses import dataclass, field
from enum import Enum, auto

class MotionType(Enum):
    """Enum for different types of motions we track."""
    EARTH_ROTATION = auto()
    EARTH_ORBIT = auto()
    MOON_ROTATION = auto()
    MOON_ORBIT = auto()
    SUN_ROTATION = auto()

@dataclass
class MotionTracker:
    """Data class to hold the times and angles for objects we track."""
    last_event_times: dict[MotionType, float] = field(default_factory=lambda: {body: 0.0 for body in MotionType})
    full_angle_times: dict[MotionType, float] = field(default_factory=lambda: {body: 0.0 for body in MotionType})
    angles: dict[MotionType, float] = field(default_factory=lambda: {body: 0.0 for body in MotionType})
    totals: dict[MotionType, float] = field(default_factory=lambda: {body: 0.0 for body in MotionType})
