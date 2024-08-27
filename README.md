## Orbit Simulator
(Package `orbits`)

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

- `time_scale_factor`: This factor increases the simulation's time reference vs. real-time.  
                   This allows the simulation to progress faster than reality so  
                   that observing rotations and orbits is possible.  

- `dist_scale_factor`: This factor decreases the orbital distance vs. the actual distance.  
                   This allows easier viewing of the planets. Without this scaling, planets  
                   are generally too small to view because orbital distances are relatively  
                   much larger than the planet sizes.  

## Modules:
- `orbit_simulator`: Main module responsible for setting up and running the orbit simulation
- `celestial_body`: Defines classes and data structures for representing celestial bodies
- `orbit`: Defines classes and data structures for representing orbits
- `constants`: Constants for the orbits package

## Documentation
For detailed API documentation, see:
[Orbit Simulator API Documentation](https://jim-tooker.github.io/orbits/docs/orbits/index.html)

## Sample Screenshot
<img width="1509" alt="orbits-screenshot" src="https://github.com/user-attachments/assets/abca1b0c-a3c6-4a7e-9778-f43632b80673">


