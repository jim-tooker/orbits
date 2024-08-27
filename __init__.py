"""
Welcome to the Orbit Simulator.

This program simulates the rotation and orbital physics of celestial bodies.  For each celestial
body, the following are visualized:
    - Rotation period around its axis
    - Axial tilt with respect to its orbit
    - Orbital inclination
    - Orbital period
    - Orbital radius
    - Path of the orbit, including direction and semi-major, semi-minor axis, and eccentricity

To ease visualizing the simulation, there are two scaling factors that can *optionally* be used
when running this program:

time_scale_factor: This factor increases the simulation's time reference vs. real-time.  
                   This allows the simulation to progress faster than reality so  
                   that observing rotations and orbits is possible.  

dist_scale_factor: This factor decreases the orbital distance vs. the actual distance.  
                   This allows easier viewing of the planets. Without this scaling, planets  
                   are generally too small to view because orbital distances are relatively  
                   much larger than the planet sizes.
"""

__author__ = "Jim Tooker"
