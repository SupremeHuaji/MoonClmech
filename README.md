# MoonClmech - Classical Mechanics Calculation Library

**MoonClmech** is a comprehensive classical mechanics calculation library written in **MoonBit**. It provides functions for particle motion, rigid body mechanics, Lagrangian and Hamiltonian mechanics, fluid dynamics, orbital mechanics, and more.

---

## Features

### Core Kinematics

* **Uniform Motion**: Displacement calculations for constant velocity
* **Accelerated Motion**: Displacement, velocity, and acceleration calculations with constant acceleration
* **Projectile Motion**: Horizontal and vertical displacement, maximum height, range, flight time, velocity components
* **Circular Motion**: Centripetal acceleration and force, angular velocity conversions, period and frequency
* **Relative Motion**: Velocity, displacement, and acceleration transformations between reference frames

---

### Dynamics and Forces

* **Newton's Laws**: Force, acceleration, and mass relationships
* **Momentum**: Linear momentum, momentum change, and impulse
* **Energy**: Kinetic energy, gravitational potential energy, elastic potential energy, mechanical energy
* **Work and Power**: Work calculations from force and displacement, power from force and velocity
* **Friction**: Static and kinetic friction, inclined plane forces and accelerations
* **Spring Forces**: Hooke's law, spring constant calculations

---

### Rigid Body Mechanics

* **Moments of Inertia**: Point mass, rod, ring, disk, sphere, cylinder, cone, rectangular plate, and more
* **Angular Motion**: Angular displacement, velocity, acceleration, and momentum
* **Torque**: Torque calculations from force and lever arm
* **Rotational Dynamics**: Rotational kinetic energy, rotational power, precession
* **Parallel Axis Theorem**: Moment of inertia for shifted axes
* **Rolling Motion**: Rolling velocity, kinetic energy, acceleration on inclined surfaces

---

### Simple Harmonic Motion and Oscillations

* **Simple Harmonic Motion**: Displacement, velocity, acceleration, angular frequency, period, frequency
* **Spring-Mass Systems**: Natural frequency and period calculations
* **Pendulums**: Simple and physical pendulum periods and frequencies
* **Damped Oscillations**: Damping coefficients, critical damping, damped angular frequency
* **Forced Oscillations**: Forced oscillation amplitude, resonance frequency
* **Coupled Oscillators**: Normal mode frequencies, normal coordinates, beat frequency

---

### Collisions and Conservation Laws

* **Elastic Collisions**: One-dimensional and two-dimensional elastic collision calculations
* **Inelastic Collisions**: Perfectly inelastic and partially inelastic collisions
* **Coefficient of Restitution**: Restitution coefficient calculation
* **Kinetic Energy Loss**: Energy loss calculations for inelastic collisions
* **Momentum Conservation**: Verification of momentum conservation in collisions
* **Center of Mass**: Center of mass position, velocity, and acceleration

---

### Orbital Mechanics and Central Forces

* **Kepler's Laws**: Orbital period from semi-major axis (Kepler's third law)
* **Orbital Velocities**: Circular orbital velocity, velocities at periapsis and apoapsis
* **Escape Velocity**: Escape velocity calculations
* **Cosmic Velocities**: First and second cosmic velocities
* **Orbital Energy**: Energy calculations for circular and elliptical orbits
* **Hohmann Transfer**: Delta-v calculations for Hohmann transfer orbits
* **Areal Velocity**: Areal velocity for orbital motion

---

### Lagrangian and Hamiltonian Mechanics

* **Lagrangian Formulation**: Lagrangian function from kinetic and potential energy
* **Generalized Forces**: Generalized force calculations
* **Hamiltonian Formulation**: Hamiltonian from kinetic and potential energy
* **Phase Space**: Phase space volume calculations

---

### Multi-body Systems

* **Center of Mass**: One-dimensional center of mass calculations
* **Total Momentum**: Total momentum for multiple objects
* **Total Kinetic Energy**: Total kinetic energy for multiple objects
* **Center of Mass Motion**: Center of mass velocity and acceleration

---

### Fluid Mechanics

* **Buoyancy**: Buoyant force (Archimedes' principle)
* **Hydrostatics**: Hydrostatic pressure at depth
* **Bernoulli's Principle**: Velocity and pressure calculations from Bernoulli's equation
* **Continuity Equation**: Velocity calculations from continuity equation
* **Viscosity**: Reynolds number, Stokes drag, viscous forces
* **Poiseuille Flow**: Flow rate in cylindrical pipes
* **Surface Tension**: Capillary rise, Laplace pressure
* **Drag Forces**: Linear and quadratic drag, terminal velocity

---

### Wave Mechanics

* **Wave Properties**: Wave displacement, velocity, frequency, wavelength relationships
* **Standing Waves**: Fundamental frequency, harmonic frequencies, node spacing
* **Doppler Effect**: Frequency shifts for moving sources and observers
* **Wave Intensity**: Wave intensity calculations
* **Wave Energy**: Energy and power in waves

---

### Extended Mechanics

* **Non-inertial Frames**: Centrifugal force, Coriolis force, inertial forces
* **Special Relativity**: Lorentz factor, relativistic mass, relativistic momentum, mass-energy equivalence
* **Continuum Mechanics**: Stress (normal and shear), strain, Young's modulus
* **Acoustics**: Sound intensity, sound intensity level, sound speed in ideal gases
* **Rocket Mechanics**: Rocket equation, thrust force, variable mass systems

---

### Physical Constants and Material Properties

* **Fundamental Constants**: Speed of light, gravitational constant, Planck constant
* **Planetary Constants**: Mass, radius, surface gravity for Earth, Moon, Sun, and other planets
* **Material Properties**: Density, viscosity, surface tension, Young's modulus for common materials
* **Astronomical Constants**: Astronomical unit, light-year, parsec, planetary distances

---

### Utility Functions

* **Vector Operations**: 3D vector creation, addition, subtraction, dot product, cross product, magnitude, normalization, projection, reflection
* **Unit Conversions**: 
  * Speed (m/s ↔ km/h)
  * Force (N ↔ lbf)
  * Energy (J ↔ eV, J ↔ cal, J ↔ BTU)
  * Power (W ↔ hp, W ↔ BTU/h)
  * Temperature (°C ↔ K ↔ °F, °F ↔ K)
  * Pressure (Pa ↔ atm, Pa ↔ bar, Pa ↔ psi)
* **Coordinate System Conversions**: 
  * 2D polar ↔ Cartesian coordinates
  * 3D spherical ↔ Cartesian coordinates
  * 3D cylindrical ↔ Cartesian coordinates
* **Dimensional Analysis**: Dimension checking for velocity, acceleration, force, energy, power, pressure
* **Physical Validation**: Range validation, non-relativistic velocity checks
* **Numerical Tools**: Maximum/minimum, interpolation, range mapping, range checking, averaging, RMS

---

## Usage

### Basic Kinematics

```moonbit
test "basic kinematics" {
  // Uniform motion: displacement = velocity * time
  let displacement = uniform_motion_displacement(10.0, 5.0)
  assert_eq(displacement, 50.0)

  // Accelerated motion: displacement with constant acceleration
  let s = uniform_acceleration_displacement(0.0, 2.0, 5.0)
  assert_eq(s, 25.0)

  // Projectile motion: maximum height and range
  let v0 = 20.0  // m/s
  let angle = 3.14159 / 4.0  // 45 degrees
  let g = 9.8  // m/s²
  let max_h = projectile_max_height(v0, angle, g)
  let range = projectile_range(v0, angle, g)
  assert_eq(max_h > 0.0, true)
  assert_eq(range > 0.0, true)
}
```

---

### Dynamics and Energy

```moonbit
test "dynamics and energy" {
  // Newton's second law: F = ma
  let force = force(5.0, 4.0)  // mass = 5 kg, acceleration = 4 m/s²
  assert_eq(force, 20.0)

  // Momentum and impulse
  let p = momentum(2.0, 10.0)  // mass = 2 kg, velocity = 10 m/s
  assert_eq(p, 20.0)
  let impulse_val = impulse(10.0, 5.0)  // force = 10 N, time = 5 s
  assert_eq(impulse_val, 50.0)

  // Energy calculations
  let ke = kinetic_energy(2.0, 10.0)  // mass = 2 kg, velocity = 10 m/s
  assert_eq(ke, 100.0)
  let pe = gravitational_potential_energy(5.0, 10.0, 9.8)  // mass = 5 kg, height = 10 m
  assert_eq(pe, 490.0)

  // Work and power
  let work_done = work(10.0, 5.0, 0.0)  // force = 10 N, displacement = 5 m, angle = 0°
  assert_eq(work_done, 50.0)
  let power_val = power(100.0, 5.0)  // work = 100 J, time = 5 s
  assert_eq(power_val, 20.0)
}
```

---

### Rigid Body Mechanics

```moonbit
test "rigid body mechanics" {
  // Moment of inertia for different shapes
  let i_disk = moment_of_inertia_disk(4.0, 2.0)  // mass = 4 kg, radius = 2 m
  assert_eq(i_disk, 8.0)
  let i_sphere = moment_of_inertia_sphere(5.0, 2.0)  // mass = 5 kg, radius = 2 m
  assert_eq(i_sphere, 8.0)

  // Angular motion
  let l = angular_momentum(10.0, 5.0)  // I = 10 kg·m², ω = 5 rad/s
  assert_eq(l, 50.0)
  let ke_rot = rotational_kinetic_energy(8.0, 3.0)  // I = 8 kg·m², ω = 3 rad/s
  assert_eq(ke_rot, 36.0)

  // Torque
  let tau = torque(10.0, 2.0, 3.14159 / 2.0)  // force = 10 N, lever = 2 m, angle = 90°
  assert_eq((tau - 20.0).abs() < 0.0001, true)
}
```

---

### Collisions

```moonbit
test "collisions" {
  // Elastic collision: equal masses
  let m1 = 1.0  // kg
  let v1i = 10.0  // m/s
  let m2 = 1.0  // kg
  let v2i = 0.0  // m/s
  
  let v1f = elastic_collision_v1(m1, v1i, m2, v2i)
  let v2f = elastic_collision_v2(m1, v1i, m2, v2i)
  assert_eq((v1f - 0.0).abs() < 0.0001, true)
  assert_eq((v2f - 10.0).abs() < 0.0001, true)

  // Verify momentum conservation
  let p_before = momentum(m1, v1i) + momentum(m2, v2i)
  let p_after = momentum(m1, v1f) + momentum(m2, v2f)
  assert_eq((p_before - p_after).abs() < 0.0001, true)

  // Inelastic collision
  let v_final = inelastic_collision_velocity(2.0, 10.0, 3.0, 0.0)  // m1=2 kg, v1=10 m/s, m2=3 kg, v2=0 m/s
  assert_eq(v_final, 4.0)
}
```

---

### Orbital Mechanics

```moonbit
test "orbital mechanics" {
  // Kepler's third law: orbital period
  let a = 1.496e11  // m (Earth-Sun distance)
  let m_sun = 1.989e30  // kg
  let g_const = 6.67430e-11
  let period = kepler_period(a, m_sun, g_const)
  assert_eq((period - 3.156e7).abs() < 0.1e7, true)  // Approximately 1 year

  // Circular orbital velocity
  let v_orbital = orbital_velocity_circular(a, m_sun, g_const)
  assert_eq((v_orbital - 29780.0).abs() < 100.0, true)

  // Escape velocity
  let r_earth = 6.371e6  // m
  let m_earth = 5.972e24  // kg
  let v_escape = escape_velocity(r_earth, m_earth, g_const)
  assert_eq((v_escape - 11186.0).abs() < 200.0, true)
}
```

---

### Simple Harmonic Motion

```moonbit
test "simple harmonic motion" {
  // Spring-mass system
  let k = 100.0  // N/m
  let m = 4.0  // kg
  let omega = shm_angular_frequency(k, m)
  assert_eq((omega - 5.0).abs() < 0.0001, true)

  // SHM displacement, velocity, acceleration
  let a = 5.0  // m (amplitude)
  let t = 0.0  // s
  let phi = 0.0  // initial phase
  let x = shm_displacement(a, omega, t, phi)
  let v = shm_velocity(a, omega, t, phi)
  assert_eq((x - 5.0).abs() < 0.0001, true)
  assert_eq((v - 0.0).abs() < 0.0001, true)

  // Pendulum
  let l = 1.0  // m
  let g = 9.8  // m/s²
  let period = pendulum_period(g, l)
  assert_eq((period - 2.0).abs() < 0.1, true)
}
```

---

### Fluid Mechanics

```moonbit
test "fluid mechanics" {
  // Buoyant force (Archimedes' principle)
  let rho_fluid = 1000.0  // kg/m³ (water)
  let volume = 0.1  // m³
  let f_buoyant = buoyant_force(rho_fluid, volume, 9.8)
  assert_eq((f_buoyant - 980.0).abs() < 0.0001, true)

  // Hydrostatic pressure
  let depth = 10.0  // m
  let p_surface = 101325.0  // Pa
  let pressure = hydrostatic_pressure(rho_fluid, depth, 9.8, p_surface)
  assert_eq(pressure, 199325.0)

  // Bernoulli's principle
  let p1 = 200000.0  // Pa
  let p2 = 150000.0  // Pa
  let v1 = 5.0  // m/s
  let v2 = bernoulli_velocity(p1, p2, rho_fluid, v1)
  assert_eq(v2 > 0.0, true)

  // Continuity equation
  let a1 = 0.1  // m²
  let a2 = 0.05  // m²
  let v2_cont = continuity_equation_velocity(a1, v1, a2)
  assert_eq(v2_cont, 10.0)
}
```

---

### Vector Operations

```moonbit
test "vector operations" {
  // Create 3D vectors
  let v1 = vec3(1.0, 2.0, 3.0)
  let v2 = vec3(4.0, 5.0, 6.0)

  // Vector operations
  let dot_product = vec3_dot(v1, v2)
  let magnitude = vec3_magnitude(v1)
  
  assert_eq(magnitude > 0.0, true)
  assert_eq(dot_product > 0.0, true)
  
  // Normalize vector
  let normalized = vec3_normalize(v1)
  let mag_norm = vec3_magnitude(normalized)
  assert_eq((mag_norm - 1.0).abs() < 0.0001, true)
}
```

---

### Unit Conversions

```moonbit
test "unit conversions" {
  // Speed conversions
  let speed_mps = 10.0  // m/s
  let speed_kmph = mps_to_kmph(speed_mps)
  assert_eq(speed_kmph, 36.0)
  let speed_back = kmph_to_mps(speed_kmph)
  assert_eq((speed_back - speed_mps).abs() < 0.0001, true)

  // Force conversions
  let force_n = 100.0  // N
  let force_lbf = n_to_lbf(force_n)
  let force_back = lbf_to_n(force_lbf)
  assert_eq((force_back - force_n).abs() < 1.0, true)

  // Energy conversions
  let energy_j = 1.602176634e-19  // J (1 eV)
  let energy_ev = joule_to_ev(energy_j)
  assert_eq((energy_ev - 1.0).abs() < 0.0001, true)
  
  let energy_cal = 1000.0  // cal
  let energy_j_from_cal = cal_to_joule(energy_cal)
  assert_eq((energy_j_from_cal - 4184.0).abs() < 0.1, true)
  
  let energy_btu = 1.0  // BTU
  let energy_j_from_btu = btu_to_joule(energy_btu)
  assert_eq((energy_j_from_btu - 1055.06).abs() < 0.1, true)

  // Power conversions
  let power_hp = 1.0  // hp
  let power_w = hp_to_watt(power_hp)
  assert_eq((power_w - 745.7).abs() < 0.1, true)
  
  let power_btu_h = 1000.0  // BTU/h
  let power_w_from_btu = btu_per_h_to_watt(power_btu_h)
  assert_eq((power_w_from_btu - 293.07).abs() < 0.1, true)

  // Temperature conversions
  let temp_c = 25.0  // Celsius
  let temp_k = celsius_to_kelvin(temp_c)
  assert_eq((temp_k - 298.15).abs() < 0.01, true)
  let temp_f = celsius_to_fahrenheit(temp_c)
  assert_eq((temp_f - 77.0).abs() < 0.1, true)
  let temp_k_from_f = fahrenheit_to_kelvin(temp_f)
  assert_eq((temp_k_from_f - temp_k).abs() < 0.1, true)

  // Pressure conversions
  let pressure_atm = 1.0  // atm
  let pressure_pa = atm_to_pascal(pressure_atm)
  assert_eq((pressure_pa - 101325.0).abs() < 0.1, true)
  
  let pressure_bar = 1.0  // bar
  let pressure_pa_from_bar = bar_to_pascal(pressure_bar)
  assert_eq((pressure_pa_from_bar - 100000.0).abs() < 0.1, true)
  
  let pressure_psi = 14.7  // psi
  let pressure_pa_from_psi = psi_to_pascal(pressure_psi)
  assert_eq((pressure_pa_from_psi - 101325.0).abs() < 100.0, true)
}
```

---

### Coordinate System Conversions

```moonbit
test "coordinate conversions" {
  // 2D Polar to Cartesian
  let r = 10.0
  let theta = 3.14159 / 4.0  // 45 degrees
  let (x, y) = polar_to_cartesian_2d(r, theta)
  assert_eq((x - 7.071).abs() < 0.1, true)
  assert_eq((y - 7.071).abs() < 0.1, true)
  
  // 2D Cartesian to Polar
  let (r_back, theta_back) = cartesian_to_polar_2d(x, y)
  assert_eq((r_back - r).abs() < 0.01, true)
  
  // 3D Spherical to Cartesian
  let r_sphere = 10.0
  let theta_sphere = 3.14159 / 3.0  // 60 degrees
  let phi_sphere = 3.14159 / 4.0    // 45 degrees
  let (x3, y3, z3) = spherical_to_cartesian_3d(r_sphere, theta_sphere, phi_sphere)
  assert_eq(x3 > 0.0, true)
  assert_eq(y3 > 0.0, true)
  assert_eq(z3 > 0.0, true)
  
  // 3D Cartesian to Spherical
  let (r_back3, theta_back3, phi_back3) = cartesian_to_spherical_3d(x3, y3, z3)
  assert_eq((r_back3 - r_sphere).abs() < 0.01, true)
  
  // 3D Cylindrical to Cartesian
  let rho = 10.0
  let phi_cyl = 3.14159 / 3.0  // 60 degrees
  let z_cyl = 5.0
  let (x_cyl, y_cyl, z_cyl_out) = cylindrical_to_cartesian_3d(rho, phi_cyl, z_cyl)
  assert_eq((z_cyl_out - z_cyl).abs() < 0.01, true)
  
  // 3D Cartesian to Cylindrical
  let (rho_back, phi_back, z_back) = cartesian_to_cylindrical_3d(x_cyl, y_cyl, z_cyl_out)
  assert_eq((rho_back - rho).abs() < 0.01, true)
  assert_eq((z_back - z_cyl).abs() < 0.01, true)
}
```

---

### Dimensional Analysis and Validation

```moonbit
test "dimensional analysis" {
  // Dimensional checks
  let velocity = dimensional_check_velocity(100.0, 10.0)  // 100 m / 10 s
  assert_eq(velocity, 10.0)
  
  let acceleration = dimensional_check_acceleration(50.0, 5.0)  // 50 m / (5 s)²
  assert_eq(acceleration, 2.0)
  
  let force = dimensional_check_force(10.0, 2.0)  // 10 kg * 2 m/s²
  assert_eq(force, 20.0)
  
  let energy = dimensional_check_energy(100.0, 5.0)  // 100 N * 5 m
  assert_eq(energy, 500.0)
  
  let power = dimensional_check_power(1000.0, 10.0)  // 1000 J / 10 s
  assert_eq(power, 100.0)
  
  let pressure = dimensional_check_pressure(1000.0, 2.0)  // 1000 N / 2 m²
  assert_eq(pressure, 500.0)
  
  // Physical range validation
  let is_valid = validate_physical_range(50.0, 0.0, 100.0)
  assert_eq(is_valid, true)
  
  let is_invalid = validate_physical_range(150.0, 0.0, 100.0)
  assert_eq(is_invalid, false)
  
  // Non-relativistic check
  let c = 299792458.0  // m/s (speed of light)
  let v_slow = 100.0  // m/s
  let is_non_rel = is_non_relativistic(v_slow, c, 0.1)  // threshold = 10% of c
  assert_eq(is_non_rel, true)
  
  let v_fast = 0.5 * c  // 50% of speed of light
  let is_rel = is_non_relativistic(v_fast, c, 0.1)
  assert_eq(is_rel, false)
}
```

---

## Parameter Ranges

### Valid Input Ranges

* **Velocity**: Any real number (typically 0 to ~10⁸ m/s for classical mechanics)
* **Acceleration**: Any real number (typically -100 to 100 m/s² for most applications)
* **Mass**: Positive values (typically 10⁻³⁰ to 10³⁰ kg)
* **Force**: Any real number (typically -10⁶ to 10⁶ N for most applications)
* **Time**: Non-negative values (typically 0 to 10¹⁸ s)
* **Distance/Length**: Non-negative values (typically 10⁻¹⁵ to 10²⁶ m)
* **Angle**: Any real number in radians (typically 0 to 2π for most applications)
* **Pressure**: Positive values (typically 10³ to 10⁹ Pa)
* **Density**: Positive values (typically 10⁻² to 10⁴ kg/m³)
* **Temperature**: Any real number in Kelvin (typically 1 to 10⁶ K)

---

### Typical Physical Values

* **Earth Surface Gravity**: 9.80665 m/s²
* **Speed of Light**: 2.998×10⁸ m/s
* **Gravitational Constant**: 6.674×10⁻¹¹ m³/(kg·s²)
* **Earth Mass**: 5.972×10²⁴ kg
* **Earth Radius**: 6.371×10⁶ m
* **Solar Mass**: 1.989×10³⁰ kg
* **Astronomical Unit**: 1.496×10¹¹ m
* **Standard Atmospheric Pressure**: 101325 Pa
* **Water Density**: ~1000 kg/m³

---

## Testing

The project includes a comprehensive test suite covering all major functionalities:

```bash
moon test
```

### Test Coverage

* Basic kinematics (uniform motion, accelerated motion, projectile motion)
* Dynamics (forces, momentum, energy, work, power)
* Rigid body mechanics (moments of inertia, angular motion, torque)
* Simple harmonic motion and oscillations
* Collisions (elastic, inelastic, conservation laws)
* Orbital mechanics (Kepler's laws, orbital velocities, escape velocity)
* Fluid mechanics (buoyancy, Bernoulli's principle, continuity equation)
* Wave mechanics (wave properties, Doppler effect, standing waves)
* Vector operations and unit conversions
* Coordinate system transformations (polar, spherical, cylindrical)
* Dimensional analysis and physical validation
* Boundary conditions and edge cases
* Conservation law verification

---

## Technical Details

### Standards and References

* **Classical Mechanics Textbooks**: Goldstein, Marion & Thornton, Taylor
* **NIST Physical Constants Database**: Fundamental physical constants
* **Standard Physics Formulas**: Well-established formulas from physics literature
* **SI Units**: All calculations use International System of Units

---

### Mathematical Methods

* **Algebraic Calculations**: Direct formula applications
* **Trigonometric Functions**: For angles, projectile motion, vector operations
* **Exponential Functions**: For damped oscillations, rocket equation
* **Numerical Integration**: For some extended calculations (e.g., CAPE in meteorology-inspired applications)
* **Iterative Methods**: For some complex calculations (where applicable)
* **Vector Mathematics**: 3D vector operations using standard vector algebra

---

## Notes

1. **Units**: All calculations use SI units (meters, kilograms, seconds, Newtons, Joules, etc.)
2. **Precision**: Floating-point calculations may have numerical precision limitations
3. **Boundary Conditions**: Values outside typical physical ranges may yield unexpected results
4. **Angles**: All angles must be in radians (use `@math.pi` or convert from degrees)
5. **Classical Mechanics**: This library focuses on classical (non-relativistic) mechanics. For relativistic effects, some functions are provided (e.g., `lorentz_factor`, `relativistic_mass`), but the primary focus is on classical physics
6. **Coordinate Systems**: Most functions assume standard Cartesian coordinate systems
7. **Conservation Laws**: Functions verify conservation of momentum and energy where applicable
8. **Gravity**: Default gravity value uses Earth's surface gravity (9.80665 m/s²), but can be specified for other planets

---

## Application Scenarios

* Physics education and teaching
* Engineering calculations (mechanics, dynamics, structural analysis)
* Simulation and modeling (projectile motion, orbital mechanics, collisions)
* Game development (physics engines, collision detection)
* Research in classical mechanics
* Scientific computing and data analysis
* Robotics (kinematics, dynamics, control systems)
* Aerospace engineering (orbital mechanics, rocket dynamics)
* Mechanical engineering (vibrations, rotating machinery, fluid systems)
* Educational physics demonstrations and visualizations

---

## Version Information

The current version implements approximately **300+ classical mechanics functions**, covering:

* Core kinematics and dynamics
* Rigid body mechanics (moments of inertia, angular motion)
* Simple harmonic motion and oscillations
* Collisions and conservation laws
* Orbital mechanics (Kepler's laws, orbital velocities)
* Fluid mechanics (buoyancy, Bernoulli's principle, viscosity)
* Wave mechanics (wave properties, Doppler effect)
* Lagrangian and Hamiltonian mechanics
* Extended features (relativity, continuum mechanics, acoustics)
* Physical constants and material properties
* Utility functions (vectors, unit conversions, coordinate transformations, dimensional analysis, numerical tools)

The library is actively developed and aims to provide comprehensive coverage of classical mechanics while being optimized for MoonBit.

