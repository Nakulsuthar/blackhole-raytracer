# blackhole-raytracer


A from scratch ray tracer simulating photon geodesics around a Schwarzschild (non rotating) black hole. Built as a summer project to visualise gravitational lensing, the photon sphere, and accretion disk emission using general relativity.

## Physics

Light paths are computed by integrating the Binet equation derived from the Schwarzschild metric:

$$\frac{d^2u}{d\phi^2} = -u + \frac{3}{2} r_s u^2$$

where $u = 1/r$ and $r_s = 2GM/c^2$ is the Schwarzschild radius.

## Roadmap

- [x] GR theory and derivation
- [ ] Single-ray integrator with RK4
- [ ] Camera model (pixel → impact parameter)
- [ ] Full-image render with gravitational lensing
- [ ] Accretion disk with redshift and Doppler beaming
- [ ] Polish and stretch goals (Kerr metric / GPU port)


## Project structure

    blackhole-raytracer/
    ├── main.py                # entry point — wires everything together
    ├── tracer/
    │   ├── integrator.py      # RK4 stepper + Binet ODE + trace_ray()
    │   ├── camera.py          # pixel → impact parameter mapping
    │   └── renderer.py        # ray result → pixel color (starfield, disk, redshift)
    ├── outputs/               # rendered images
    └── tests/
        └── test_integrator.py # validation checks