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

## Modules:
- `orbit_simulator`: Main module responsible for setting up and running the orbit simulation.
- `simulation_mode`: Handles what simulation mode we are in and all the corresponding GUI components.
- `motion_tracker`: Defines classes that keep track of the motion of the celestial bodies.
- `celestial_body`: Defines classes and data structures for representing celestial bodies.
- `orbit`: Defines classes and data structures for representing orbits.
- `config`: Stores the configuration information shared among the modules.
- `constants`: Constants for the orbits package.

## Documentation
For detailed API documentation, see:
[Orbit Simulator API Documentation](https://jim-tooker.github.io/orbits/docs/orbits/index.html)

## Sample Screenshots
<img width="1509" alt="orbits-screenshot" src="https://github.com/user-attachments/assets/abca1b0c-a3c6-4a7e-9778-f43632b80673">
<img width="1497" alt="Orbits with Sun-Screenshot" src="https://github.com/user-attachments/assets/8d678c19-4ded-462b-bae4-1fa0840d9824">


