"""
FIXME
"""
from abc import ABC, abstractmethod
from typing import List, Final
import vpython as vp
from orbits import config
from orbits.celestial_body import Earth, Moon, Sun
from orbits.orbit import Orbit, EarthOrbit, MoonOrbit
from orbits.motion_tracker import MotionTracker, MotionType
from orbits.constants import FULL_ANGLE, HRS_IN_DAY, SECS_IN_HR

class SimulationMode(ABC):
    @abstractmethod
    def create_celestial_bodies(self):
        pass

    @abstractmethod
    def update_celestial_bodies(self, t_prime, dt_prime):
        pass

    def _create_info_label(self,
                           text: str,
                           left_margin: int,
                           line_number: int,
                           color: vp.vector=vp.color.white) -> vp.label:
        """
        Create a label with given text at specified position on the info canvas.

        Args:
            text (str): The text to display on the label.
            left_margin (int): The left margin position for the label.
            line_number (int): The vertical position (line number) for the label.
            color (vp.vector): The color of the label text. Defaults to white.

        Returns:
            vp.label: The created label object.
        """
        return vp.label(pos=vp.vector(left_margin, line_number, 0),
                                      text=text,
                                      color=color,
                                      height=16,
                                      align='left',
                                      box=False)



class SunEarthMoonMode(SimulationMode):
    def __init__(self):
        self.sun: Sun
        self.earth: Earth
        self.moon: Moon

        # Create self.tracker so we can monitor times and angles of objects
        self.tracker: MotionTracker = MotionTracker()

    def create_celestial_bodies(self):
        # Create Sun
        self.sun = Sun()

        # Create orbits List
        orbits: List[Orbit] = []

        # Create Earth and its orbit
        earth_orbit: EarthOrbit = EarthOrbit()
        orbits.append(earth_orbit)
        self.earth = Earth(orbits=orbits)

        # Create Moon and its orbit
        moon_orbit: MoonOrbit = MoonOrbit()
        orbits.append(moon_orbit)
        self.moon = Moon(orbits=orbits,
                         earth=self.earth)

    def update_celestial_bodies(self, t_prime, dt_prime):
        # Update celestial bodies orbit position based on total time elapsed
        self.earth.update_position(t_prime)
        self.moon.update_position(t_prime)

        # Rotate celestial bodies based on the small time step and their respective angular velocities
        self.sun.rotate(dt_prime)
        self.earth.rotate(dt_prime)
        self.moon.rotate(dt_prime)

    def check_for_full_angle(self, t: float) -> None:
        # If the Earth has rotated 360°, store rotation time
        if abs(self.earth.angle(t) - self.tracker.angles[MotionType.EARTH_ROTATION]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.EARTH_ROTATION] = \
                t - self.tracker.last_event_times[MotionType.EARTH_ROTATION]
            self.tracker.last_event_times[MotionType.EARTH_ROTATION] = t
            self.tracker.angles[MotionType.EARTH_ROTATION] = self.earth.angle(t)
            self.tracker.totals[MotionType.EARTH_ROTATION] += 1
            print(f'Earth Rotation {self.tracker.totals[MotionType.EARTH_ROTATION]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.EARTH_ROTATION]/(SECS_IN_HR):.2f} hours')

        # If the Moon has orbited 360°, store the orbit time
        if abs(self.moon.orbit.angle(t) - self.tracker.angles[MotionType.MOON_ORBIT]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.MOON_ORBIT] = \
                t - self.tracker.last_event_times[MotionType.MOON_ORBIT]
            self.tracker.last_event_times[MotionType.MOON_ORBIT] = t
            self.tracker.angles[MotionType.MOON_ORBIT] = self.moon.orbit.angle(t)
            self.tracker.totals[MotionType.MOON_ORBIT] += 1
            print(f'Moon Orbit {self.tracker.totals[MotionType.MOON_ORBIT]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.MOON_ORBIT]/(SECS_IN_HR*HRS_IN_DAY):.2f} days')
            print(f'There were {self.tracker.totals[MotionType.EARTH_ROTATION] / \
                self.tracker.totals[MotionType.MOON_ORBIT]:.2f} Earth rotations during the last Moon orbit.')

        # If the Moon has rotated 360°, store rotation time
        if abs(self.moon.angle(t) - self.tracker.angles[MotionType.MOON_ROTATION]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.MOON_ROTATION] = \
                t - self.tracker.last_event_times[MotionType.MOON_ROTATION]
            self.tracker.last_event_times[MotionType.MOON_ROTATION] = t
            self.tracker.angles[MotionType.MOON_ROTATION] = self.moon.angle(t)
            self.tracker.totals[MotionType.MOON_ROTATION] += 1
            print(f'Moon Rotation {self.tracker.totals[MotionType.MOON_ROTATION]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.MOON_ROTATION]/(SECS_IN_HR*HRS_IN_DAY):.2f} days')

        # If the Earth has orbited 360° around the Sun, store the orbit time
        if abs(self.earth.orbit.angle(t) - self.tracker.angles[MotionType.EARTH_ORBIT]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.EARTH_ORBIT] = \
                t - self.tracker.last_event_times[MotionType.EARTH_ORBIT]
            self.tracker.last_event_times[MotionType.EARTH_ORBIT] = t
            self.tracker.angles[MotionType.EARTH_ORBIT] = self.earth.orbit.angle(t)
            self.tracker.totals[MotionType.EARTH_ORBIT] += 1
            print(f'Earth Orbit {self.tracker.totals[MotionType.EARTH_ORBIT]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.EARTH_ORBIT]/(SECS_IN_HR*HRS_IN_DAY):.2f} days')

        # If the Sun has rotated 360°, store rotation time
        if abs(self.sun.angle(t) - self.tracker.angles[MotionType.SUN_ROTATION]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.SUN_ROTATION] = \
                t - self.tracker.last_event_times[MotionType.SUN_ROTATION]
            self.tracker.last_event_times[MotionType.SUN_ROTATION] = t
            self.tracker.angles[MotionType.SUN_ROTATION] = self.sun.angle(t)
            self.tracker.totals[MotionType.SUN_ROTATION] += 1
            print(f'Sun Rotation {self.tracker.totals[MotionType.SUN_ROTATION]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.SUN_ROTATION]/(SECS_IN_HR*HRS_IN_DAY):.2f} days')

    def setup_info_canvas(self) -> vp.canvas:
        """
        Set up the information canvas with labels displaying simulation details.

        This method creates a separate canvas for displaying information about
        the simulation, including time and distance scales, and details about
        the Earth and Moon.

        Returns:
            vp.canvas: The create info canvas
        """
        # Create canvas and configure
        width: Final[int] = 400
        height: Final[int] = 1000
        info_canvas: vp.canvas = vp.canvas(width=width, height=height, align='left')
        info_canvas.range = 10
        info_canvas.userzoom = False
        info_canvas.userspin = False
        info_canvas.userpan = False

        # Start labels 1 over from left margin (range -10 to 10)
        left_margin: int = -info_canvas.range + 1

        # Start lines at the top minus 2
        line_number: int = int(info_canvas.range * height/width - 2)

        # Create time scale label
        self._create_info_label(f'Time scale: {config.time_scale_factor:,.0f}x. 1 sec = {
              config.time_scale_factor/(SECS_IN_HR*HRS_IN_DAY):.1f} day(s).',
              left_margin, line_number)
        line_number -= 2

        # Create Sun info label
        self._create_info_label('Sun Info:\n' +
                                f'  Radius: {self.sun.params.radius:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> {self.sun.scale_factor:.0f}x\n' +
                                f'  Tilt: {self.sun.params.tilt_degrees:.1f}°\n' +
                                f'  Rotation period: {self.sun.params.rotation_period_days:.2f} days',
                                left_margin, line_number)
        line_number -= 6

        # Create Earth info label
        self._create_info_label('Earth Info:\n' +
                                f'  Radius: {self.earth.params.radius:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> {self.earth.scale_factor:.0f}x\n' +
                                f'  Tilt: {self.earth.params.tilt_degrees:.1f}°\n' +
                                f'  Sidereal day: {self.earth.SIDEREAL_DAY:.2f} hrs',
                                left_margin, line_number)
        line_number -= 6

        # Create Moon info label
        self._create_info_label('Moon Info:\n' +
                                f'  Radius: {self.moon.params.radius:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> {self.moon.scale_factor:.0f}x\n' +
                                f'  Tilt: {self.moon.params.tilt_degrees:.1f}°\n' +
                                f'  Sidereal month: {self.moon.SIDEREAL_MONTH:.2f} days',
                                left_margin, line_number)
        line_number -= 6

        # Create Earth's Orbit info label
        self._create_info_label("Earth's Orbit Info:\n" +
                                f'  Semi-major axis: {self.earth.orbit.params.semi_major_axis:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> 1/{1/self.earth.orbit.scale_factor:.0f}x\n'
                                f'  Eccentricity: {self.earth.orbit.params.eccentricity:.3f}\n' +
                                f'  Inclination: {self.earth.orbit.params.inclination_degs:.2f}°\n' +
                                f'  Period: {self.earth.orbit.params.period_days:.2f} days',
                                left_margin, line_number)
        line_number -= 7

        # Create Moon's Orbit info label
        self._create_info_label("Moon's Orbit Info:\n" +
                                f'  Semi-major axis: {self.moon.orbit.params.semi_major_axis:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> 1/{1/self.moon.orbit.scale_factor:.0f}x\n'
                                f'  Eccentricity: {self.moon.orbit.params.eccentricity:.3f}\n' +
                                f'  Inclination: {self.moon.orbit.params.inclination_degs:.2f}°\n' +
                                f'  Period: {self.moon.orbit.params.period_days:.2f} days',
                                left_margin, line_number)
        line_number -= 7

        # Create Camera view info label
        self._create_info_label(text='To change view:\n' +
                                     'Rotate view: Drag with right mouse button.\n' +
                                     'Zoom: Scroll wheel or drag left/right mouse buttons.\n' +
                                     'Pan: Shift drag with left mouse button.',
                                left_margin=left_margin,
                                line_number=line_number,
                                color=vp.color.yellow)
        line_number -= 5

        # Create Sim view scale info label
        self._create_info_label(text='Note:  Celestial body sizes and orbits have\n' +
                                     '           been scaled to fit screen.\n' +
                                     'See "Simulation visual scale factor" notes above.',
                                left_margin=left_margin,
                                line_number=line_number,
                                color=vp.color.orange)
        line_number -= 4

        # Create Time remaining label
        self._runtime_left_label = self._create_info_label(text='',
                                                           left_margin=left_margin,
                                                           line_number=line_number)
        line_number -= 2

        return info_canvas


class EarthMoonMode(SimulationMode):
    def __init__(self):
        self.earth: Earth
        self.moon: Moon

    def create_celestial_bodies(self):
        # Create orbits List
        orbits: List[Orbit] = []

        # Create Earth
        self.earth = Earth()

        # Create Moon and its orbit
        moon_orbit: MoonOrbit = MoonOrbit()
        orbits.append(moon_orbit)
        self.moon = Moon(orbits=orbits,
                         earth=self.earth)

    def update_celestial_bodies(self, t_prime, dt_prime):
        # Update celestial bodies orbit position based on total time elapsed
        self.moon.update_position(t_prime)

        # Rotate celestial bodies based on the small time step and their respective angular velocities
        self.earth.rotate(dt_prime)
        self.moon.rotate(dt_prime)

    def check_for_full_angle(self, t: float) -> None:
        """
        FIXME
        """
        pass

    def setup_info_canvas(self) -> vp.canvas:
        return vp.canvas()
