"""
FIXME
"""
from abc import ABC, abstractmethod
from typing import List, Final
import math
import vpython as vp
from orbits import config
from orbits.celestial_body import Earth, Moon, Sun
from orbits.orbit import Orbit, EarthOrbit, MoonOrbit
from orbits.motion_tracker import MotionTracker, MotionType
from orbits.constants import FULL_ANGLE, HRS_IN_DAY, SECS_IN_HR

class SimulationMode(ABC):
    """FIXME"""
    def __init__(self):
        # Create celestial bodies
        self.sun: Sun
        self.earth: Earth
        self.moon: Moon

        # Create self.tracker so we can monitor times and angles of objects
        self.tracker: MotionTracker = MotionTracker()

        self.runtime_left_label: vp.label
        self._info_canvas_left_margin: float = 0
        self._info_canvas_line_number: float = 0

        if config.no_gui is False:
            self._canvas: vp.canvas = vp.canvas(title='Orbit Simulator',
                                                width=1600,
                                                height=1000,
                                                align='left')

        self._create_celestial_bodies()

        if config.no_gui is False:
            self._create_orientation_figure()
            self._info_canvas: vp.canvas

    @abstractmethod
    def update_celestial_bodies(self, t_prime: float, dt_prime: float) -> None:
        """FIXME"""
        pass

    @abstractmethod
    def _create_celestial_bodies(self) -> None:
        pass

    @abstractmethod
    def _largest_orbit_mag(self) -> float:
        pass

    def add_quit_button(self, sim):
        """FIXME"""
        # Reduce the height of the info canvas by the height of the quit button
        height_of_quit_button: Final[int] = 25
        self._info_canvas.height -= height_of_quit_button

        # Add the Quit simulation button, and bind it to handle_quit_button()
        vp.button(pos=self._info_canvas.title_anchor,
                text='                                Quit Simulation                                ',
                color=vp.color.red,
                background=vp.color.gray(0.8),
                bind=sim.handle_quit_button)

    def _create_orientation_figure(self) -> None:
        """Create orientation arrows and labels to show the x, y, and z axes."""
        # Define orientation figure size and position using the largest object on the canvas.
        canvas_mag: Final[float] = self._largest_orbit_mag()
        orient_size: Final[float] = canvas_mag/8
        orient_pos: vp.vector = vp.vector(-canvas_mag, canvas_mag/2, 0)

        # Create the arrows in each direction of x,y,z axis
        vp.arrow(pos=orient_pos,
                 axis=vp.vector(orient_size, 0, 0),
                 color=vp.color.red,
                 round=True,
                 emissive=True)
        vp.arrow(pos=orient_pos,
                 axis=vp.vector(0, orient_size, 0),
                 color=vp.color.green,
                 round=True,
                 emissive=True)
        vp.arrow(pos=orient_pos,
                 axis=vp.vector(0, 0, orient_size),
                 color=vp.color.yellow,
                 round=True,
                 emissive=True)

        # Label the arrows
        label_distance: float = orient_size * 1.1
        vp.label(pos=orient_pos + vp.vector(label_distance, 0, 0),
                 text='X',
                 color=vp.color.red,
                 opacity=0,
                 box=False)
        vp.label(pos=orient_pos + vp.vector(0, label_distance, 0),
                 text='Y',
                 color=vp.color.green,
                 opacity=0,
                 box=False)
        vp.label(pos=orient_pos + vp.vector(0, 0, label_distance),
                 text='Z',
                 color=vp.color.yellow,
                 opacity=0,
                 box=False)

    def _create_info_label(self,
                           text: str,
                           color: vp.vector=vp.color.white) -> vp.label:
        """
        Create a label with given text at specified position on the info canvas.

        Args:
            text (str): The text to display on the label.
            color (vp.vector): The color of the label text. Defaults to white.

        Returns:
            vp.label: The created label object.
        """
        return vp.label(pos=vp.vector(self._info_canvas_left_margin,
                                      self._info_canvas_line_number,
                                      0),
                        text=text,
                        color=color,
                        height=16,
                        align='left',
                        box=False)

    def _setup_info_canvas(self) -> None:
        """
        Set up the information canvas with labels displaying simulation details.

        This method creates a separate canvas for displaying information about
        the simulation, including time and distance scales, and details about
        the Earth and Moon.
        """
        # Create canvas and configure
        width: Final[int] = 400
        height: Final[int] = 1000
        self._info_canvas = vp.canvas(width=width, height=height, align='left')
        self._info_canvas.range = 10
        self._info_canvas.userzoom = False
        self._info_canvas.userspin = False
        self._info_canvas.userpan = False

        # Start labels 1 over from left margin (range -10 to 10)
        self._info_canvas_left_margin = -self._info_canvas.range + 1

        # Start lines at the top minus 2
        self._info_canvas_line_number = int(self._info_canvas.range * height/width - 2)

        # Create Sim view scale info label
        self._create_info_label(text='Note:  Celestial body sizes and orbits have\n' +
                                     '           been scaled to fit screen.\n' +
                                     'See "Simulation visual scale factor" notes below.',
                                color=vp.color.orange)
        self._info_canvas_line_number -= 4

        # Create Camera view info label
        self._create_info_label(text='To change view:\n' +
                                     'Rotate view: Drag with right mouse button.\n' +
                                     'Zoom: Scroll wheel or drag left/right mouse buttons.\n' +
                                     'Pan: Shift drag with left mouse button.',
                                color=vp.color.yellow)
        self._info_canvas_line_number -= 5

        # Create Time remaining label
        self.runtime_left_label = self._create_info_label(text='')
        self._info_canvas_line_number -= 2

    def _rotate_camera_angle(self, angle: float) -> None:
        self._canvas.camera.rotate(angle=-math.radians(angle), axis=vp.vector(1, 0, 0))

    def _check_earth_rotation(self, t: float) -> None:
        # If the Earth has rotated 360°, store rotation time
        if abs(self.earth.angle(t) - self.tracker.angles[MotionType.EARTH_ROTATION]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.EARTH_ROTATION] = \
                t - self.tracker.last_event_times[MotionType.EARTH_ROTATION]
            self.tracker.last_event_times[MotionType.EARTH_ROTATION] = t
            self.tracker.angles[MotionType.EARTH_ROTATION] = self.earth.angle(t)
            self.tracker.totals[MotionType.EARTH_ROTATION] += 1
            print(f'Earth Rotation {self.tracker.totals[MotionType.EARTH_ROTATION]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.EARTH_ROTATION]/(SECS_IN_HR):.2f} hours')

    def _check_moon_rotation(self, t: float) -> None:
        # If the Moon has rotated 360°, store rotation time
        if abs(self.moon.angle(t) - self.tracker.angles[MotionType.MOON_ROTATION]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.MOON_ROTATION] = \
                t - self.tracker.last_event_times[MotionType.MOON_ROTATION]
            self.tracker.last_event_times[MotionType.MOON_ROTATION] = t
            self.tracker.angles[MotionType.MOON_ROTATION] = self.moon.angle(t)
            self.tracker.totals[MotionType.MOON_ROTATION] += 1
            print(f'Moon Rotation {self.tracker.totals[MotionType.MOON_ROTATION]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.MOON_ROTATION]/(SECS_IN_HR*HRS_IN_DAY):.2f} days')

    def _check_moon_orbit(self, t: float) -> None:
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

    def _check_earth_orbit(self, t: float) -> None:
        # If the Earth has orbited 360° around the Sun, store the orbit time
        if abs(self.earth.orbit.angle(t) - self.tracker.angles[MotionType.EARTH_ORBIT]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.EARTH_ORBIT] = \
                t - self.tracker.last_event_times[MotionType.EARTH_ORBIT]
            self.tracker.last_event_times[MotionType.EARTH_ORBIT] = t
            self.tracker.angles[MotionType.EARTH_ORBIT] = self.earth.orbit.angle(t)
            self.tracker.totals[MotionType.EARTH_ORBIT] += 1
            print(f'Earth Orbit {self.tracker.totals[MotionType.EARTH_ORBIT]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.EARTH_ORBIT]/(SECS_IN_HR*HRS_IN_DAY):.2f} days')

    def _check_sun_rotation(self, t: float) -> None:
        # If the Sun has rotated 360°, store rotation time
        if abs(self.sun.angle(t) - self.tracker.angles[MotionType.SUN_ROTATION]) >= FULL_ANGLE:
            self.tracker.full_angle_times[MotionType.SUN_ROTATION] = \
                t - self.tracker.last_event_times[MotionType.SUN_ROTATION]
            self.tracker.last_event_times[MotionType.SUN_ROTATION] = t
            self.tracker.angles[MotionType.SUN_ROTATION] = self.sun.angle(t)
            self.tracker.totals[MotionType.SUN_ROTATION] += 1
            print(f'Sun Rotation {self.tracker.totals[MotionType.SUN_ROTATION]:.0f}: Time: {
                self.tracker.full_angle_times[MotionType.SUN_ROTATION]/(SECS_IN_HR*HRS_IN_DAY):.2f} days')

    def _create_sun_info_label(self) -> None:
        """
        FIXME
        """
        # Create Sun info label
        self._create_info_label('Sun Info:\n' +
                                f'  Radius: {self.sun.params.radius:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> {self.sun.scale_factor:.0f}x\n' +
                                f'  Tilt: {self.sun.params.tilt_degrees:.1f}°\n' +
                                f'  Rotation period: {self.sun.params.rotation_period_days:.2f} days')
        self._info_canvas_line_number -= 6

    def _create_earth_info_label(self) -> None:
        """
        FIXME
        """
        # Create Earth info label
        self._create_info_label('Earth Info:\n' +
                                f'  Radius: {self.earth.params.radius:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> {self.earth.scale_factor:.0f}x\n' +
                                f'  Tilt: {self.earth.params.tilt_degrees:.1f}°\n' +
                                f'  Sidereal day: {self.earth.SIDEREAL_DAY:.2f} hrs')
        self._info_canvas_line_number -= 6

    def _create_moon_info_label(self) -> None:
        """
        FIXME
        """
        # Create Moon info label
        self._create_info_label('Moon Info:\n' +
                                f'  Radius: {self.moon.params.radius:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> {self.moon.scale_factor:.0f}x\n' +
                                f'  Tilt: {self.moon.params.tilt_degrees:.1f}°\n' +
                                f'  Sidereal month: {self.moon.SIDEREAL_MONTH:.2f} days')
        self._info_canvas_line_number -= 6

    def _create_earth_orbit_info_label(self) -> None:
        """
        FIXME
        """
        # Create Earth's Orbit info label
        self._create_info_label("Earth's Orbit Info:\n" +
                                f'  Semi-major axis: {self.earth.orbit.params.semi_major_axis:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> 1/{1/self.earth.orbit.scale_factor:.0f}x\n'
                                f'  Eccentricity: {self.earth.orbit.params.eccentricity:.3f}\n' +
                                f'  Inclination: {self.earth.orbit.params.inclination_degs:.2f}°\n' +
                                f'  Period: {self.earth.orbit.params.period_days:.2f} days')
        self._info_canvas_line_number -= 7

    def _create_moon_orbit_info_label(self) -> None:
        """
        FIXME
        """
        # Create Moon's Orbit info label
        self._create_info_label("Moon's Orbit Info:\n" +
                                f'  Semi-major axis: {self.moon.orbit.params.semi_major_axis:,.0f} km\n' +
                                f'    <i>Simulation visual scale factor:</i> 1/{1/self.moon.orbit.scale_factor:.0f}x\n'
                                f'  Eccentricity: {self.moon.orbit.params.eccentricity:.3f}\n' +
                                f'  Inclination: {self.moon.orbit.params.inclination_degs:.2f}°\n' +
                                f'  Period: {self.moon.orbit.params.period_days:.2f} days')
        self._info_canvas_line_number -= 7

    def _create_time_scale_info_label(self )-> None:
        """
        FIXME
        """
        # Create time scale label
        self._create_info_label(f'Time scale: {config.time_scale_factor:,.0f}x. 1 sec = {
              config.time_scale_factor/(SECS_IN_HR*HRS_IN_DAY):.1f} day(s).')
        self._info_canvas_line_number -= 2


class SunEarthMoonMode(SimulationMode):
    """
    FIXME
    """
    def __init__(self):
        super().__init__()

        # Set time scale, if not set
        if config.time_scale_factor == 0:
            default_time_scale_factor: Final[float] = 10 * SECS_IN_HR * HRS_IN_DAY  # 10 day = 1 sec
            config.time_scale_factor = default_time_scale_factor

        if config.no_gui is False:
            self._setup_info_canvas()

        # rotate camera around the x axis to see the orbits better (not straight on)
        camera_rotate_angle: Final[float] = 30  # degrees
        self._rotate_camera_angle(camera_rotate_angle)

    def update_celestial_bodies(self, t_prime: float, dt_prime: float) -> None:
        """
        FIXME
        """
        # Update celestial bodies orbit position based on total time elapsed
        self.earth.update_position(t_prime)
        self.moon.update_position(t_prime)

        # Rotate celestial bodies based on the small time step and their respective angular velocities
        self.sun.rotate(dt_prime)
        self.earth.rotate(dt_prime)
        self.moon.rotate(dt_prime)

    def check_for_full_angle(self, t: float) -> None:
        """
        FIXME
        """
        self._check_earth_rotation(t)
        self._check_moon_rotation(t)
        self._check_moon_orbit(t)
        self._check_earth_orbit(t)
        self._check_sun_rotation(t)

    def _create_celestial_bodies(self) -> None:
        """
        FIXME
        """
        # Create Sun
        self.sun = Sun()

        # Create orbits List
        orbits: List[Orbit] = []

        # Create Earth and its orbit
        earth_orbit_scale_factor: Final[float] = 1/30
        earth_orbit: EarthOrbit = EarthOrbit(scale_factor=earth_orbit_scale_factor)
        orbits.append(earth_orbit)
        earth_scale_factor: Final[float] = 8
        self.earth = Earth(scale_factor=earth_scale_factor, orbits=orbits)

        # Create Moon and its orbit
        moon_orbit_scale_factor: Final[float] = 1/4
        moon_orbit: MoonOrbit = MoonOrbit(scale_factor=moon_orbit_scale_factor)
        orbits.append(moon_orbit)
        moon_scale_factor: Final[float] = 8
        self.moon = Moon(scale_factor=moon_scale_factor, orbits=orbits, earth=self.earth)

    def _setup_info_canvas(self) -> None:
        super()._setup_info_canvas()

        self._create_time_scale_info_label()
        self._create_sun_info_label()
        self._create_earth_info_label()
        self._create_moon_info_label()
        self._create_earth_orbit_info_label()
        self._create_moon_orbit_info_label()

        # Set default canvas back to normal canvas
        self._canvas.select()

    def _largest_orbit_mag(self) -> float:
        """
        FIXME
        """
        return self.earth.orbit.orbit_mag


class EarthMoonMode(SimulationMode):
    """FIXME"""
    def __init__(self):
        super().__init__()

        # Set time scale, if not set
        if config.time_scale_factor == 0:
            default_time_scale_factor: Final[float] = SECS_IN_HR * HRS_IN_DAY  # 1 day = 1 sec
            config.time_scale_factor = default_time_scale_factor

        if config.no_gui is False:
            self._setup_info_canvas()

        # rotate camera around the x axis to see the orbits better (not straight on)
        camera_rotate_angle: Final[float] = 2  # degrees
        self._rotate_camera_angle(camera_rotate_angle)

    def update_celestial_bodies(self, t_prime: float, dt_prime: float) -> None:
        """
        FIXME
        """
        # Update celestial bodies orbit position based on total time elapsed
        self.moon.update_position(t_prime)

        # Rotate celestial bodies based on the small time step and their respective angular velocities
        self.earth.rotate(dt_prime)
        self.moon.rotate(dt_prime)

    def check_for_full_angle(self, t: float) -> None:
        """
        FIXME
        """
        self._check_earth_rotation(t)
        self._check_moon_rotation(t)
        self._check_moon_orbit(t)

    def _create_celestial_bodies(self) -> None:
        """
        FIXME
        """
        # Create orbits List
        orbits: List[Orbit] = []

        # Create Earth
        self.earth = Earth()

        # Create Moon and its orbit
        moon_orbit_scale_factor: Final[float] = 1/10
        moon_orbit: MoonOrbit = MoonOrbit(scale_factor=moon_orbit_scale_factor)
        orbits.append(moon_orbit)
        self.moon = Moon(orbits=orbits, earth=self.earth)

    def _setup_info_canvas(self) -> None:
        super()._setup_info_canvas()

        self._create_time_scale_info_label()
        self._create_earth_info_label()
        self._create_moon_info_label()
        self._create_moon_orbit_info_label()

        # Set default canvas back to normal canvas
        self._canvas.select()

    def _largest_orbit_mag(self) -> float:
        """
        FIXME
        """
        return self.moon.orbit.orbit_mag
