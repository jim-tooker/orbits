"""
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

Here are the runtime options: `python orbit_simulator.py -h`:
```
options:
  -h, --help            show this help message and exit
  -t TIME_SCALE_FACTOR, --time_scale_factor TIME_SCALE_FACTOR
                        How much to scale up the sense of time. 0 means use best-fit time scaling. (default: 0)
  -r RUNTIME, --runtime RUNTIME
                        How many secs to run the simulation. 0 means run indefinitely. (default: 0)
  -m {sun_earth_moon,earth_moon}, --mode {sun_earth_moon,earth_moon}
                        Simulation mode to use. (default: sun_earth_moon)
  --no_gui              Run without GUI (default: False)
```

"""
__author__ = "Jim Tooker"
