"""
Simulation Mode Module

This module defines classes that handle the simulation mode we are in.
Currently, there are two possible modes:  
- Sun, Earth, and Moon:  The simulator shows the Sun with the Earth's and Moon's orbits around it.  
- Earth and Moon: The simulator shows the Earth and Moon's orbit around it.  

This module also handles all the VPython GUI and canvas updates.
"""
from abc import ABC, abstractmethod
from typing import List, Final
import math
import vpython as vp
from constants import FULL_ANGLE, HRS_IN_DAY, SECS_IN_HR
import config
from celestial_body import Earth, Moon, Sun, TrailParams
from orbit import Orbit, EarthOrbit, MoonOrbit
from motion_tracker import MotionTracker, MotionType

__author__ = "Jim Tooker"


class SimulationMode(ABC):
    """Abstract class to handle all the common simulation mode methods and states."""
    def __init__(self) -> None:
        # Create celestial bodies
        self.sun: Sun
        self.earth: Earth
        self.moon: Moon

        # Create tracker so we can monitor times and angles of objects
        self.tracker: MotionTracker = MotionTracker()

        # Label to update the runtime left
        self.runtime_left_label: vp.label

        # Flag to indicate when the user has pressed the quit simulation button
        self.quit_sim: bool = False

        # Variables for the info canvas left margin and current line number
        self._info_canvas_left_margin: float = 0
        self._info_canvas_line_number: float = 0

        # Create main canvas if GUI enabled
        if config.no_gui is False:
            self._canvas: vp.canvas = vp.canvas(title='Orbit Simulator',
                                                width=1600,
                                                height=1000,
                                                align='left')

        # Create all the celestial bodies (handled in the sub classes)
        self._create_celestial_bodies()

        if config.no_gui is False:
            self._create_orientation_figure()

            # Create info canvas (will be populated by sub classes)
            self._info_canvas: vp.canvas

    def __del__(self) -> None:
        """
        Deletes the canvases to allow canvases to disappear from GUI
        """
        if hasattr(self, '_canvas') and self._canvas:
            self._canvas.delete()

        if hasattr(self, '_info_canvas') and self._info_canvas:
            self._info_canvas.delete()

    @abstractmethod
    def update_celestial_bodies(self, t: float, dt: float) -> None:
        """
        Abstract method to update the positions of celestial bodies based on time.

        Args:
          t (float): The current time.
          dt (float): The current delta time.
        """

    @abstractmethod
    def check_for_full_angle(self, t: float) -> None:
        """
        Abstract method to check if any celestial bodies have made a full rotation or full orbit.

        Args:
          t (float): The current time.
        """

    @abstractmethod
    def _create_celestial_bodies(self) -> None:
        """Abstract method to create the celestial bodies."""

    @abstractmethod
    def _largest_orbit_mag(self) -> float:
        """
        Abstract method to determine the largest orbit on the canvas.
        
        Returns:
          float:  The magnitude of the largest orbit.
        """

    def _handle_quit_button(self, _: vp.button) -> None:
        """Handles quit simulation button"""
        self.quit_sim = True

    def add_quit_button(self) -> None:
        """Adds a "Quit" button to the info canvas."""
        # Reduce the height of the info canvas by the height of the quit button
        height_of_quit_button: Final[int] = 25
        self._info_canvas.height -= height_of_quit_button

        # Add the Quit simulation button, and bind it to _handle_quit_button()
        vp.button(pos=self._info_canvas.title_anchor,
                text='                                Quit Simulation                                ',
                color=vp.color.red,
                background=vp.color.gray(0.8),
                bind=self._handle_quit_button)

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
        the celestial bodies.
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
        """
        Rotates the camera angle for a better initial viewing angle
        
        Args:
          angle (float): The angle to rotate the camera (in degrees).
        """
        self._canvas.camera.rotate(angle=-math.radians(angle), axis=vp.vector(1, 0, 0))

    def _check_earth_rotation(self, t: float) -> None:
        """
        Check if the Earth has done a full rotation. If so, store the information in the tracker.
        
        Args:
          t (float): The current time.
        """
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
        """
        Check if the Moon has done a full rotation. If so, store the information in the tracker.
        
        Args:
          t (float): The current time.
        """
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
        """
        Check if the Moon has done a full orbit. If so, store the information in the tracker.
        
        Args:
          t (float): The current time.
        """
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
        """
        Check if the Earth has done a full orbit. If so, store the information in the tracker.
        
        Args:
          t (float): The current time.
        """
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
        """
        Check if the Sun has done a full rotation. If so, store the information in the tracker.
        
        Args:
          t (float): The current time.
        """
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
        Create a label with the Sun's info.
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
        Create a label with the Earth's info.
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
        Create a label with the Moon's info.
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
        Create a label with the Earth's orbit info.
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
        Create a label with the Moon's orbit info.
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
        Create a label with time scale info.
        """
        # Create time scale label
        self._create_info_label(f'Time scale: {config.time_scale_factor:,.0f}x. 1 sec = {
              config.time_scale_factor/(SECS_IN_HR*HRS_IN_DAY):.1f} day(s).')
        self._info_canvas_line_number -= 2


class SunEarthMoonMode(SimulationMode):
    """
    Subclass to handle the Simulation mode of Sun, Earth, and Moon.
    """
    def __init__(self) -> None:
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

    def update_celestial_bodies(self, t: float, dt: float) -> None:
        """
        Updates the positions of celestial bodies based on time.

        Args:
          t (float): The current time.
          dt (float): The current delta time.
        """
        # Update celestial bodies orbit position based on total time elapsed
        self.earth.update_position(t)
        self.moon.update_position(t)

        # Rotate celestial bodies based on the small time step and their respective angular velocities
        self.sun.rotate(dt)
        self.earth.rotate(dt)
        self.moon.rotate(dt)

    def check_for_full_angle(self, t: float) -> None:
        """
        Check if any celestial bodies have made a full rotation or full orbit.

        Args:
          t (float): The current time.
        """
        self._check_earth_rotation(t)
        self._check_moon_rotation(t)
        self._check_moon_orbit(t)
        self._check_earth_orbit(t)
        self._check_sun_rotation(t)

    def _create_celestial_bodies(self) -> None:
        """Creates the celestial bodies."""
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
        trail_params: TrailParams = TrailParams(
            trail_radius=0.1*self.earth.radius,
            trail_retain=1000)
        self.earth.set_trail_params(trail_params)

        # Create Moon and its orbit
        moon_orbit_scale_factor: Final[float] = 1/4
        moon_orbit: MoonOrbit = MoonOrbit(scale_factor=moon_orbit_scale_factor)
        orbits.append(moon_orbit)
        moon_scale_factor: Final[float] = 8
        self.moon = Moon(scale_factor=moon_scale_factor, orbits=orbits, earth=self.earth)
        trail_params = TrailParams(
            trail_radius=0.25*self.moon.radius,
            trail_retain=1000)
        self.moon.set_trail_params(trail_params)

    def _setup_info_canvas(self) -> None:
        """
        Set up the information canvas with labels displaying simulation details.

        This method creates a separate canvas for displaying information about
        the simulation, including time and distance scales, and details about
        the celestial bodies.
        """
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
        Determines the largest orbit on the canvas.
        
        Returns:
          float:  The magnitude of the largest orbit.
        """
        return self.earth.orbit.orbit_mag


class EarthMoonMode(SimulationMode):
    """
    Subclass to handle the Simulation mode of Earth and Moon.
    """
    def __init__(self) -> None:
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

    def update_celestial_bodies(self, t: float, dt: float) -> None:
        """
        Updates the positions of celestial bodies based on time.

        Args:
          t (float): The current time.
          dt (float): The current delta time.
        """
        # Update celestial bodies orbit position based on total time elapsed
        self.moon.update_position(t)

        # Rotate celestial bodies based on the small time step and their respective angular velocities
        self.earth.rotate(dt)
        self.moon.rotate(dt)

    def check_for_full_angle(self, t: float) -> None:
        """
        Check if any celestial bodies have made a full rotation or full orbit.

        Args:
          t (float): The current time.
        """
        self._check_earth_rotation(t)
        self._check_moon_rotation(t)
        self._check_moon_orbit(t)

    def _create_celestial_bodies(self) -> None:
        """Creates the celestial bodies."""
        # Create orbits List
        orbits: List[Orbit] = []

        # Create Earth
        self.earth = Earth()

        # Create Moon and its orbit
        moon_orbit_scale_factor: Final[float] = 1/10
        moon_orbit: MoonOrbit = MoonOrbit(scale_factor=moon_orbit_scale_factor)
        orbits.append(moon_orbit)
        self.moon = Moon(orbits=orbits, earth=self.earth)
        trail_params: TrailParams = TrailParams(
            trail_radius=0,
            trail_retain=1000)
        self.moon.set_trail_params(trail_params)

    def _setup_info_canvas(self) -> None:
        """
        Set up the information canvas with labels displaying simulation details.

        This method creates a separate canvas for displaying information about
        the simulation, including time and distance scales, and details about
        the celestial bodies.
        """
        super()._setup_info_canvas()

        self._create_time_scale_info_label()
        self._create_earth_info_label()
        self._create_moon_info_label()
        self._create_moon_orbit_info_label()

        # Set default canvas back to normal canvas
        self._canvas.select()

    def _largest_orbit_mag(self) -> float:
        """
        Determines the largest orbit on the canvas.
        
        Returns:
          float:  The magnitude of the largest orbit.
        """
        return self.moon.orbit.orbit_mag
