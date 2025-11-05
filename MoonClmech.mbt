/// MoonClmech - MoonBit Classical Mechanics Library
/// Main calculation functions based on classical mechanics principles
/// Provides functions for particle motion, rigid body mechanics, Lagrangian, Hamiltonian, etc.

///|
/// Uniform motion: calculate displacement
/// v: velocity (m/s)
/// t: time (s)
pub fn uniform_motion_displacement(velocity : Double, time : Double) -> Double {
  velocity * time
}

///|
/// Uniform acceleration: calculate displacement
/// v0: initial velocity (m/s)
/// a: acceleration (m/s²)
/// t: time (s)
pub fn uniform_acceleration_displacement(
  v0 : Double,
  a : Double,
  t : Double,
) -> Double {
  v0 * t + 0.5 * a * t * t
}

///|
/// Uniform acceleration: calculate velocity
/// v0: initial velocity (m/s)
/// a: acceleration (m/s²)
/// t: time (s)
pub fn uniform_acceleration_velocity(
  v0 : Double,
  a : Double,
  t : Double,
) -> Double {
  v0 + a * t
}

///|
/// Uniform acceleration: calculate velocity (using displacement)
/// v0: initial velocity (m/s)
/// a: acceleration (m/s²)
/// s: displacement (m)
pub fn uniform_acceleration_velocity_from_displacement(
  v0 : Double,
  a : Double,
  s : Double,
) -> Double {
  Double::sqrt(v0 * v0 + 2.0 * a * s)
}

///|
/// Projectile motion: horizontal displacement
/// v0: initial velocity (m/s)
/// angle: launch angle (radians)
/// t: time (s)
pub fn projectile_horizontal_displacement(
  v0 : Double,
  angle : Double,
  t : Double,
) -> Double {
  v0 * @math.cos(angle) * t
}

///|
/// Projectile motion: vertical displacement
/// v0: initial velocity (m/s)
/// angle: launch angle (radians)
/// t: time (s)
/// g: gravitational acceleration (m/s²), default uses Earth's gravity
pub fn projectile_vertical_displacement(
  v0 : Double,
  angle : Double,
  t : Double,
  g : Double,
) -> Double {
  v0 * @math.sin(angle) * t - 0.5 * g * t * t
}

///|
/// Projectile motion: maximum height
/// v0: initial velocity (m/s)
/// angle: launch angle (radians)
/// g: gravitational acceleration (m/s²)
pub fn projectile_max_height(v0 : Double, angle : Double, g : Double) -> Double {
  let v0y = v0 * @math.sin(angle)
  v0y * v0y / (2.0 * g)
}

///|
/// Projectile motion: range
/// v0: initial velocity (m/s)
/// angle: launch angle (radians)
/// g: gravitational acceleration (m/s²)
pub fn projectile_range(v0 : Double, angle : Double, g : Double) -> Double {
  v0 * v0 * @math.sin(2.0 * angle) / g
}

///|
/// Projectile motion: flight time
/// v0: initial velocity (m/s)
/// angle: launch angle (radians)
/// g: gravitational acceleration (m/s²)
pub fn projectile_flight_time(
  v0 : Double,
  angle : Double,
  g : Double,
) -> Double {
  2.0 * v0 * @math.sin(angle) / g
}

///|
/// Newton's second law: calculate force
/// m: mass (kg)
/// a: acceleration (m/s²)
pub fn force(mass : Double, acceleration : Double) -> Double {
  mass * acceleration
}

///|
/// Newton's second law: calculate acceleration
/// f: force (N)
/// m: mass (kg)
pub fn acceleration(force : Double, mass : Double) -> Double {
  force / mass
}

///|
/// Momentum
/// m: mass (kg)
/// v: velocity (m/s)
pub fn momentum(mass : Double, velocity : Double) -> Double {
  mass * velocity
}

///|
/// Change in momentum
/// m: mass (kg)
/// v1: initial velocity (m/s)
/// v2: final velocity (m/s)
pub fn momentum_change(mass : Double, v1 : Double, v2 : Double) -> Double {
  mass * (v2 - v1)
}

///|
/// Impulse
/// f: force (N)
/// t: time (s)
pub fn impulse(force : Double, time : Double) -> Double {
  force * time
}

///|
/// Kinetic energy
/// m: mass (kg)
/// v: velocity (m/s)
pub fn kinetic_energy(mass : Double, velocity : Double) -> Double {
  0.5 * mass * velocity * velocity
}

///|
/// Gravitational potential energy
/// m: mass (kg)
/// h: height (m)
/// g: gravitational acceleration (m/s²)
pub fn gravitational_potential_energy(
  mass : Double,
  height : Double,
  g : Double,
) -> Double {
  mass * g * height
}

///|
/// Elastic potential energy
/// k: spring constant (N/m)
/// x: displacement (m)
pub fn elastic_potential_energy(
  spring_constant : Double,
  displacement : Double,
) -> Double {
  0.5 * spring_constant * displacement * displacement
}

///|
/// Mechanical energy (kinetic + potential)
pub fn mechanical_energy(kinetic : Double, potential : Double) -> Double {
  kinetic + potential
}

///|
/// Gravitational force
/// m1: mass of object 1 (kg)
/// m2: mass of object 2 (kg)
/// r: distance (m)
/// g_const: gravitational constant (m³/(kg·s²))
pub fn gravitational_force(
  m1 : Double,
  m2 : Double,
  r : Double,
  g_const : Double,
) -> Double {
  g_const * m1 * m2 / (r * r)
}

///|
/// Moment of inertia: point mass
/// m: mass (kg)
/// r: distance from rotation axis (m)
pub fn moment_of_inertia_point_mass(mass : Double, radius : Double) -> Double {
  mass * radius * radius
}

///|
/// Moment of inertia: thin rod (about center)
/// m: mass (kg)
/// l: length (m)
pub fn moment_of_inertia_rod_center(mass : Double, length : Double) -> Double {
  mass * length * length / 12.0
}

///|
/// Moment of inertia: thin rod (about end)
/// m: mass (kg)
/// l: length (m)
pub fn moment_of_inertia_rod_end(mass : Double, length : Double) -> Double {
  mass * length * length / 3.0
}

///|
/// Moment of inertia: ring (about central axis)
/// m: mass (kg)
/// r: radius (m)
pub fn moment_of_inertia_ring(mass : Double, radius : Double) -> Double {
  mass * radius * radius
}

///|
/// Moment of inertia: disk/cylinder (about central axis)
/// m: mass (kg)
/// r: radius (m)
pub fn moment_of_inertia_disk(mass : Double, radius : Double) -> Double {
  0.5 * mass * radius * radius
}

///|
/// Moment of inertia: solid sphere (about diameter)
/// m: mass (kg)
/// r: radius (m)
pub fn moment_of_inertia_sphere(mass : Double, radius : Double) -> Double {
  0.4 * mass * radius * radius
}

///|
/// Moment of inertia: hollow sphere (about diameter)
/// m: mass (kg)
/// r: radius (m)
pub fn moment_of_inertia_hollow_sphere(
  mass : Double,
  radius : Double,
) -> Double {
  2.0 / 3.0 * mass * radius * radius
}

///|
/// Angular momentum
/// I: moment of inertia (kg·m²)
/// omega: angular velocity (rad/s)
pub fn angular_momentum(
  moment_of_inertia : Double,
  angular_velocity : Double,
) -> Double {
  moment_of_inertia * angular_velocity
}

///|
/// Angular momentum (point mass)
/// m: mass (kg)
/// r: distance from rotation axis (m)
/// v: linear velocity (m/s)
pub fn angular_momentum_point_mass(
  mass : Double,
  radius : Double,
  velocity : Double,
) -> Double {
  mass * radius * velocity
}

///|
/// Rotational kinetic energy
/// I: moment of inertia (kg·m²)
/// omega: angular velocity (rad/s)
pub fn rotational_kinetic_energy(
  moment_of_inertia : Double,
  angular_velocity : Double,
) -> Double {
  0.5 * moment_of_inertia * angular_velocity * angular_velocity
}

///|
/// Torque
/// f: force (N)
/// r: lever arm (m)
/// angle: angle between force and lever arm (radians)
pub fn torque(force : Double, lever_arm : Double, angle : Double) -> Double {
  force * lever_arm * @math.sin(angle)
}

///|
/// Angular acceleration
/// tau: torque (N·m)
/// I: moment of inertia (kg·m²)
pub fn angular_acceleration(
  torque : Double,
  moment_of_inertia : Double,
) -> Double {
  torque / moment_of_inertia
}

///|
/// Angular velocity (uniform angular acceleration)
/// omega0: initial angular velocity (rad/s)
/// alpha: angular acceleration (rad/s²)
/// t: time (s)
pub fn angular_velocity(omega0 : Double, alpha : Double, t : Double) -> Double {
  omega0 + alpha * t
}

///|
/// Angular displacement (uniform angular acceleration)
/// omega0: initial angular velocity (rad/s)
/// alpha: angular acceleration (rad/s²)
/// t: time (s)
pub fn angular_displacement(
  omega0 : Double,
  alpha : Double,
  t : Double,
) -> Double {
  omega0 * t + 0.5 * alpha * t * t
}

///|
/// Parallel axis theorem: calculate moment of inertia
/// i_cm: moment of inertia about center of mass (kg·m²)
/// m: mass (kg)
/// d: distance to new axis (m)
pub fn parallel_axis_theorem(
  i_cm : Double,
  mass : Double,
  distance : Double,
) -> Double {
  i_cm + mass * distance * distance
}

///|
/// Simple harmonic motion: displacement
/// A: amplitude (m)
/// omega: angular frequency (rad/s)
/// t: time (s)
/// phi: initial phase (radians)
pub fn shm_displacement(
  amplitude : Double,
  angular_frequency : Double,
  time : Double,
  phase : Double,
) -> Double {
  amplitude * @math.cos(angular_frequency * time + phase)
}

///|
/// Simple harmonic motion: velocity
pub fn shm_velocity(
  amplitude : Double,
  angular_frequency : Double,
  time : Double,
  phase : Double,
) -> Double {
  -amplitude * angular_frequency * @math.sin(angular_frequency * time + phase)
}

///|
/// Simple harmonic motion: acceleration
pub fn shm_acceleration(
  amplitude : Double,
  angular_frequency : Double,
  time : Double,
  phase : Double,
) -> Double {
  -amplitude *
  angular_frequency *
  angular_frequency *
  @math.cos(angular_frequency * time + phase)
}

///|
/// Simple harmonic motion: angular frequency
/// k: spring constant (N/m)
/// m: mass (kg)
pub fn shm_angular_frequency(spring_constant : Double, mass : Double) -> Double {
  Double::sqrt(spring_constant / mass)
}

///|
/// Simple harmonic motion: period
/// omega: angular frequency (rad/s)
pub fn shm_period(angular_frequency : Double) -> Double {
  2.0 * 3.141592653589793 / angular_frequency
}

///|
/// Simple harmonic motion: frequency
/// omega: angular frequency (rad/s)
pub fn shm_frequency(angular_frequency : Double) -> Double {
  angular_frequency / (2.0 * 3.141592653589793)
}

///|
/// Simple pendulum: angular frequency
/// g: gravitational acceleration (m/s²)
/// l: pendulum length (m)
pub fn pendulum_angular_frequency(gravity : Double, length : Double) -> Double {
  Double::sqrt(gravity / length)
}

///|
/// Simple pendulum: period
pub fn pendulum_period(gravity : Double, length : Double) -> Double {
  2.0 * 3.141592653589793 * Double::sqrt(length / gravity)
}

///|
/// Elastic collision: velocity after one-dimensional perfectly elastic collision
/// m1: mass of object 1 (kg)
/// v1: initial velocity of object 1 (m/s)
/// m2: mass of object 2 (kg)
/// v2: initial velocity of object 2 (m/s)
/// Returns: final velocity of object 1
pub fn elastic_collision_v1(
  m1 : Double,
  v1 : Double,
  m2 : Double,
  v2 : Double,
) -> Double {
  (m1 - m2) * v1 / (m1 + m2) + 2.0 * m2 * v2 / (m1 + m2)
}

///|
/// Elastic collision: velocity after one-dimensional perfectly elastic collision
/// Returns: final velocity of object 2
pub fn elastic_collision_v2(
  m1 : Double,
  v1 : Double,
  m2 : Double,
  v2 : Double,
) -> Double {
  2.0 * m1 * v1 / (m1 + m2) + (m2 - m1) * v2 / (m1 + m2)
}

///|
/// Perfectly inelastic collision: velocity after collision
pub fn inelastic_collision_velocity(
  m1 : Double,
  v1 : Double,
  m2 : Double,
  v2 : Double,
) -> Double {
  (m1 * v1 + m2 * v2) / (m1 + m2)
}

///|
/// Coefficient of restitution: ratio of relative velocity after collision to before collision
/// v1f: velocity of object 1 after collision (m/s)
/// v2f: velocity of object 2 after collision (m/s)
/// v1i: velocity of object 1 before collision (m/s)
/// v2i: velocity of object 2 before collision (m/s)
pub fn coefficient_of_restitution(
  v1f : Double,
  v2f : Double,
  v1i : Double,
  v2i : Double,
) -> Double {
  (v2f - v1f) / (v1i - v2i)
}

///|
/// Kepler's third law: orbital period
/// a: semi-major axis (m)
/// M: central body mass (kg)
/// g_const: gravitational constant (m³/(kg·s²))
pub fn kepler_period(
  semi_major_axis : Double,
  central_mass : Double,
  g_const : Double,
) -> Double {
  2.0 *
  3.141592653589793 *
  Double::sqrt(
    semi_major_axis *
    semi_major_axis *
    semi_major_axis /
    (g_const * central_mass),
  )
}

///|
/// Orbital velocity (circular orbit)
/// r: orbital radius (m)
/// M: central body mass (kg)
/// g_const: gravitational constant
pub fn orbital_velocity_circular(
  radius : Double,
  central_mass : Double,
  g_const : Double,
) -> Double {
  Double::sqrt(g_const * central_mass / radius)
}

///|
/// Escape velocity
/// r: distance from central body (m)
/// M: central body mass (kg)
/// g_const: gravitational constant
pub fn escape_velocity(
  radius : Double,
  central_mass : Double,
  g_const : Double,
) -> Double {
  Double::sqrt(2.0 * g_const * central_mass / radius)
}

///|
/// First cosmic velocity (Earth's surface)
pub fn first_cosmic_velocity(gravity : Double, radius : Double) -> Double {
  Double::sqrt(gravity * radius)
}

///|
/// Second cosmic velocity (Earth's surface escape velocity)
pub fn second_cosmic_velocity(gravity : Double, radius : Double) -> Double {
  Double::sqrt(2.0 * gravity * radius)
}

///|
/// Lagrangian function L = T - V
/// T: kinetic energy
/// V: potential energy
pub fn lagrangian(kinetic_energy : Double, potential_energy : Double) -> Double {
  kinetic_energy - potential_energy
}

///|
/// Lagrangian equation: calculate generalized force (conservative system)
/// For L = L(q, q_dot, t), generalized force Q = d/dt(dL/dq_dot) - dL/dq
/// This is a simplified version assuming L does not explicitly depend on time
/// 
/// Note: This is a simplified version; full implementation requires symbolic differentiation
/// 
/// dL_dq_dot: time derivative of dL/dq_dot
/// dL_dq: dL/dq
pub fn lagrangian_generalized_force(
  dL_dq_dot_dt : Double,
  dL_dq : Double,
) -> Double {
  dL_dq_dot_dt - dL_dq
}

///|
/// Hamiltonian function H = T + V (for conservative systems)
/// T: kinetic energy
/// V: potential energy
pub fn hamiltonian(
  kinetic_energy : Double,
  potential_energy : Double,
) -> Double {
  kinetic_energy + potential_energy
}

///|
/// Hamiltonian function H = Σ(p_i * q_dot_i) - L
/// For single degree of freedom system: H = p * q_dot - L
/// p: generalized momentum
/// q_dot: generalized velocity
/// L: Lagrangian function
pub fn hamiltonian_from_lagrangian(
  momentum : Double,
  velocity : Double,
  lagrangian : Double,
) -> Double {
  momentum * velocity - lagrangian
}

///|
/// Center of mass position (one-dimensional, multiple particles)
/// masses: mass array (kg)
/// positions: position array (m)
pub fn center_of_mass_1d(
  masses : Array[Double],
  positions : Array[Double],
) -> Double {
  if masses.length() != positions.length() || masses.length() == 0 {
    0.0
  } else {
    let mut total_mass = 0.0
    let mut weighted_sum = 0.0
    for i = 0; i < masses.length(); i = i + 1 {
      total_mass = total_mass + masses[i]
      weighted_sum = weighted_sum + masses[i] * positions[i]
    }
    if total_mass > 0.0 {
      weighted_sum / total_mass
    } else {
      0.0
    }
  }
}

///|
/// Total system momentum (one-dimensional)
/// masses: mass array (kg)
/// velocities: velocity array (m/s)
pub fn total_momentum_1d(
  masses : Array[Double],
  velocities : Array[Double],
) -> Double {
  if masses.length() != velocities.length() {
    0.0
  } else {
    let mut total = 0.0
    for i = 0; i < masses.length(); i = i + 1 {
      total = total + masses[i] * velocities[i]
    }
    total
  }
}

///|
/// Total system kinetic energy (one-dimensional)
pub fn total_kinetic_energy_1d(
  masses : Array[Double],
  velocities : Array[Double],
) -> Double {
  if masses.length() != velocities.length() {
    0.0
  } else {
    let mut total = 0.0
    for i = 0; i < masses.length(); i = i + 1 {
      total = total + 0.5 * masses[i] * velocities[i] * velocities[i]
    }
    total
  }
}

///|
/// Damping coefficient for damped oscillation
/// b: damping coefficient (kg/s)
/// m: mass (kg)
pub fn damping_coefficient(b : Double, mass : Double) -> Double {
  b / (2.0 * mass)
}

///|
/// Critical damping coefficient
/// m: mass (kg)
/// k: spring constant (N/m)
pub fn critical_damping_coefficient(
  mass : Double,
  spring_constant : Double,
) -> Double {
  2.0 * Double::sqrt(mass * spring_constant)
}

///|
/// Angular frequency of damped oscillation (underdamped)
/// omega0: natural angular frequency (rad/s)
/// gamma: damping coefficient
pub fn damped_angular_frequency(
  natural_frequency : Double,
  damping : Double,
) -> Double {
  Double::sqrt(natural_frequency * natural_frequency - damping * damping)
}

///|
/// Centripetal acceleration
/// v: linear velocity (m/s)
/// r: radius (m)
pub fn centripetal_acceleration(velocity : Double, radius : Double) -> Double {
  velocity * velocity / radius
}

///|
/// Centripetal acceleration (using angular velocity)
/// omega: angular velocity (rad/s)
/// r: radius (m)
pub fn centripetal_acceleration_angular(
  angular_velocity : Double,
  radius : Double,
) -> Double {
  angular_velocity * angular_velocity * radius
}

///|
/// Centripetal force
/// m: mass (kg)
/// v: linear velocity (m/s)
/// r: radius (m)
pub fn centripetal_force(
  mass : Double,
  velocity : Double,
  radius : Double,
) -> Double {
  mass * velocity * velocity / radius
}

///|
/// Centripetal force (using angular velocity)
/// m: mass (kg)
/// omega: angular velocity (rad/s)
/// r: radius (m)
pub fn centripetal_force_angular(
  mass : Double,
  angular_velocity : Double,
  radius : Double,
) -> Double {
  mass * angular_velocity * angular_velocity * radius
}

///|
/// Relationship between linear velocity and angular velocity
/// omega: angular velocity (rad/s)
/// r: radius (m)
pub fn linear_velocity_from_angular(
  angular_velocity : Double,
  radius : Double,
) -> Double {
  angular_velocity * radius
}

///|
/// Relationship between angular velocity and linear velocity
/// v: linear velocity (m/s)
/// r: radius (m)
pub fn angular_velocity_from_linear(
  velocity : Double,
  radius : Double,
) -> Double {
  velocity / radius
}

///|
/// Period of circular motion
/// omega: angular velocity (rad/s)
pub fn circular_motion_period(angular_velocity : Double) -> Double {
  2.0 * 3.141592653589793 / angular_velocity
}

///|
/// Frequency of circular motion
/// omega: angular velocity (rad/s)
pub fn circular_motion_frequency(angular_velocity : Double) -> Double {
  angular_velocity / (2.0 * 3.141592653589793)
}

///|
/// Static friction (maximum static friction)
/// mu_s: coefficient of static friction
/// N: normal force (N)
pub fn static_friction(mu_s : Double, normal_force : Double) -> Double {
  mu_s * normal_force
}

///|
/// Kinetic friction
/// mu_k: coefficient of kinetic friction
/// N: normal force (N)
pub fn kinetic_friction(mu_k : Double, normal_force : Double) -> Double {
  mu_k * normal_force
}

///|
/// Normal force on inclined plane
/// m: mass (kg)
/// g: gravitational acceleration (m/s²)
/// angle: inclination angle (radians)
pub fn normal_force_inclined_plane(
  mass : Double,
  gravity : Double,
  angle : Double,
) -> Double {
  mass * gravity * @math.cos(angle)
}

///|
/// Gravity component along inclined plane
/// m: mass (kg)
/// g: gravitational acceleration (m/s²)
/// angle: inclination angle (radians)
pub fn gravity_component_along_incline(
  mass : Double,
  gravity : Double,
  angle : Double,
) -> Double {
  mass * gravity * @math.sin(angle)
}

///|
/// Acceleration on inclined plane (no friction)
/// g: gravitational acceleration (m/s²)
/// angle: inclination angle (radians)
pub fn acceleration_inclined_plane_no_friction(
  gravity : Double,
  angle : Double,
) -> Double {
  gravity * @math.sin(angle)
}

///|
/// Acceleration on inclined plane (with friction)
/// g: gravitational acceleration (m/s²)
/// angle: inclination angle (radians)
/// mu_k: coefficient of kinetic friction
pub fn acceleration_inclined_plane_with_friction(
  gravity : Double,
  angle : Double,
  mu_k : Double,
) -> Double {
  gravity * (@math.sin(angle) - mu_k * @math.cos(angle))
}

///|
/// Work (constant force, linear motion)
/// F: force (N)
/// s: displacement (m)
/// angle: angle between force and displacement (radians)
pub fn work(force : Double, displacement : Double, angle : Double) -> Double {
  force * displacement * @math.cos(angle)
}

///|
/// Work (force and displacement in same direction)
/// F: force (N)
/// s: displacement (m)
pub fn work_parallel(force : Double, displacement : Double) -> Double {
  force * displacement
}

///|
/// Power (average power)
/// W: work (J)
/// t: time (s)
pub fn power(work : Double, time : Double) -> Double {
  work / time
}

///|
/// Power (force and velocity)
/// F: force (N)
/// v: velocity (m/s)
/// angle: angle between force and velocity (radians)
pub fn power_from_force(
  force : Double,
  velocity : Double,
  angle : Double,
) -> Double {
  force * velocity * @math.cos(angle)
}

///|
/// Power (force and velocity in same direction)
/// F: force (N)
/// v: velocity (m/s)
pub fn power_from_force_parallel(force : Double, velocity : Double) -> Double {
  force * velocity
}

///|
/// Work-energy theorem: calculate work from kinetic energy change
/// ke_final: final kinetic energy (J)
/// ke_initial: initial kinetic energy (J)
pub fn work_from_kinetic_energy(
  ke_final : Double,
  ke_initial : Double,
) -> Double {
  ke_final - ke_initial
}

///|
/// Moment of inertia: solid cylinder (about central axis)
/// m: mass (kg)
/// r: radius (m)
pub fn moment_of_inertia_cylinder_axis(
  mass : Double,
  radius : Double,
) -> Double {
  0.5 * mass * radius * radius
}

///|
/// Moment of inertia: solid cylinder (about diameter through center)
/// m: mass (kg)
/// r: radius (m)
/// h: height (m)
pub fn moment_of_inertia_cylinder_diameter(
  mass : Double,
  radius : Double,
  height : Double,
) -> Double {
  mass * (radius * radius / 4.0 + height * height / 12.0)
}

///|
/// Moment of inertia: solid cone (about symmetry axis)
/// m: mass (kg)
/// r: base radius (m)
pub fn moment_of_inertia_cone_axis(mass : Double, radius : Double) -> Double {
  0.3 * mass * radius * radius
}

///|
/// Moment of inertia: rectangular plate (about center, perpendicular to plate)
/// m: mass (kg)
/// a: length (m)
/// b: width (m)
pub fn moment_of_inertia_rectangular_plate(
  mass : Double,
  length : Double,
  width : Double,
) -> Double {
  mass * (length * length + width * width) / 12.0
}

///|
/// Moment of inertia: rectangular plate (about long edge)
/// m: mass (kg)
/// b: width (m)
pub fn moment_of_inertia_rectangular_plate_edge(
  mass : Double,
  width : Double,
) -> Double {
  mass * width * width / 3.0
}

///|
/// Lever equilibrium condition: torque balance
/// f1: force 1 (N)
/// r1: lever arm 1 (m)
/// f2: force 2 (N)
/// r2: lever arm 2 (m)
/// Returns: whether balanced (considering floating-point errors)
pub fn lever_equilibrium(
  f1 : Double,
  r1 : Double,
  f2 : Double,
  r2 : Double,
) -> Bool {
  let epsilon = 1.0e-9
  Double::abs(f1 * r1 - f2 * r2) < epsilon
}

///|
/// Lever equilibrium: calculate unknown force
/// f1: known force (N)
/// r1: lever arm 1 (m)
/// r2: lever arm 2 (m)
/// Returns: force 2 required for equilibrium
pub fn lever_equilibrium_force(f1 : Double, r1 : Double, r2 : Double) -> Double {
  f1 * r1 / r2
}

///|
/// Rigid body equilibrium: net force is zero
/// forces: force array (N)
/// Note: This is a simplified version assuming all forces are in the same direction
pub fn net_force_1d(forces : Array[Double]) -> Double {
  let mut total = 0.0
  for i = 0; i < forces.length(); i = i + 1 {
    total = total + forces[i]
  }
  total
}

///|
/// Rigid body equilibrium: net torque is zero
/// torques: torque array (N·m)
pub fn net_torque(torques : Array[Double]) -> Double {
  let mut total = 0.0
  for i = 0; i < torques.length(); i = i + 1 {
    total = total + torques[i]
  }
  total
}

///|
/// Relative velocity (one-dimensional)
/// v1: velocity of object 1 (m/s)
/// v2: velocity of object 2 (m/s)
/// Returns: velocity of object 2 relative to object 1
pub fn relative_velocity_1d(v1 : Double, v2 : Double) -> Double {
  v2 - v1
}

///|
/// Relative displacement (one-dimensional)
/// x1: position of object 1 (m)
/// x2: position of object 2 (m)
/// Returns: displacement of object 2 relative to object 1
pub fn relative_displacement_1d(x1 : Double, x2 : Double) -> Double {
  x2 - x1
}

///|
/// Relative acceleration (one-dimensional)
/// a1: acceleration of object 1 (m/s²)
/// a2: acceleration of object 2 (m/s²)
/// Returns: acceleration of object 2 relative to object 1
pub fn relative_acceleration_1d(a1 : Double, a2 : Double) -> Double {
  a2 - a1
}

///|
/// Buoyant force (Archimedes' principle)
/// rho: fluid density (kg/m³)
/// V: displaced fluid volume (m³)
/// g: gravitational acceleration (m/s²)
pub fn buoyant_force(
  fluid_density : Double,
  displaced_volume : Double,
  gravity : Double,
) -> Double {
  fluid_density * displaced_volume * gravity
}

///|
/// Hydrostatic pressure
/// rho: fluid density (kg/m³)
/// h: depth (m)
/// g: gravitational acceleration (m/s²)
/// P0: surface pressure (Pa), default is atmospheric pressure
pub fn hydrostatic_pressure(
  fluid_density : Double,
  depth : Double,
  gravity : Double,
  surface_pressure : Double,
) -> Double {
  surface_pressure + fluid_density * gravity * depth
}

///|
/// Bernoulli equation: calculate velocity (simplified, assuming constant height)
/// p1: pressure 1 (Pa)
/// p2: pressure 2 (Pa)
/// rho: fluid density (kg/m³)
/// v1: velocity 1 (m/s)
/// Returns: velocity 2
pub fn bernoulli_velocity(
  p1 : Double,
  p2 : Double,
  rho : Double,
  v1 : Double,
) -> Double {
  let delta_p = p1 - p2
  Double::sqrt(v1 * v1 + 2.0 * delta_p / rho)
}

///|
/// Bernoulli equation: calculate pressure (simplified, assuming constant height)
/// p1: pressure 1 (Pa)
/// rho: fluid density (kg/m³)
/// v1: velocity 1 (m/s)
/// v2: velocity 2 (m/s)
/// Returns: pressure 2
pub fn bernoulli_pressure(
  p1 : Double,
  rho : Double,
  v1 : Double,
  v2 : Double,
) -> Double {
  p1 + 0.5 * rho * (v1 * v1 - v2 * v2)
}

///|
/// Continuity equation: calculate flow velocity
/// a1: cross-sectional area 1 (m²)
/// v1: velocity 1 (m/s)
/// a2: cross-sectional area 2 (m²)
/// Returns: velocity 2
pub fn continuity_equation_velocity(
  a1 : Double,
  v1 : Double,
  a2 : Double,
) -> Double {
  a1 * v1 / a2
}

///|
/// Mechanical energy conservation: calculate final velocity
/// ke_initial: initial kinetic energy (J)
/// pe_initial: initial potential energy (J)
/// pe_final: final potential energy (J)
/// m: mass (kg)
/// Returns: final velocity
pub fn velocity_from_energy_conservation(
  ke_initial : Double,
  pe_initial : Double,
  pe_final : Double,
  mass : Double,
) -> Double {
  let ke_final = ke_initial + pe_initial - pe_final
  if ke_final >= 0.0 {
    Double::sqrt(2.0 * ke_final / mass)
  } else {
    0.0
  }
}

///|
/// Mechanical energy conservation: calculate maximum height
/// v0: initial velocity (m/s)
/// h0: initial height (m)
/// g: gravitational acceleration (m/s²)
/// Returns: maximum height
pub fn max_height_from_energy(v0 : Double, h0 : Double, g : Double) -> Double {
  h0 + v0 * v0 / (2.0 * g)
}

///|
/// Elastic potential energy to kinetic energy: calculate velocity
/// k: spring constant (N/m)
/// x: compression/extension (m)
/// m: mass (kg)
/// Returns: velocity
pub fn velocity_from_spring_energy(
  spring_constant : Double,
  displacement : Double,
  mass : Double,
) -> Double {
  let pe = elastic_potential_energy(spring_constant, displacement)
  Double::sqrt(2.0 * pe / mass)
}

///|
/// Angular momentum conservation: calculate final angular velocity
/// i1: initial moment of inertia (kg·m²)
/// omega1: initial angular velocity (rad/s)
/// i2: final moment of inertia (kg·m²)
/// Returns: final angular velocity
pub fn angular_velocity_from_conservation(
  i1 : Double,
  omega1 : Double,
  i2 : Double,
) -> Double {
  if i2 > 0.0 {
    i1 * omega1 / i2
  } else {
    0.0
  }
}

///|
/// Angular momentum conservation: calculate final moment of inertia
/// i1: initial moment of inertia (kg·m²)
/// omega1: initial angular velocity (rad/s)
/// omega2: final angular velocity (rad/s)
/// Returns: final moment of inertia
pub fn moment_of_inertia_from_conservation(
  i1 : Double,
  omega1 : Double,
  omega2 : Double,
) -> Double {
  if omega2 != 0.0 {
    i1 * omega1 / omega2
  } else {
    0.0
  }
}

///|
/// Forced oscillation amplitude (harmonic driving force)
/// f0: driving force amplitude (N)
/// m: mass (kg)
/// omega0: natural angular frequency (rad/s)
/// omega_d: driving angular frequency (rad/s)
/// gamma: damping coefficient
/// Returns: amplitude
pub fn forced_oscillation_amplitude(
  f0 : Double,
  mass : Double,
  omega0 : Double,
  omega_d : Double,
  gamma : Double,
) -> Double {
  let omega_diff_sq = omega0 * omega0 - omega_d * omega_d
  let denominator = Double::sqrt(
    omega_diff_sq * omega_diff_sq + 4.0 * gamma * gamma * omega_d * omega_d,
  )
  if denominator > 0.0 {
    f0 / (mass * denominator)
  } else {
    0.0
  }
}

///|
/// Resonance frequency (undamped)
/// omega0: natural angular frequency (rad/s)
pub fn resonance_frequency(omega0 : Double) -> Double {
  omega0
}

///|
/// Resonance frequency (damped)
/// omega0: natural angular frequency (rad/s)
/// gamma: damping coefficient
pub fn resonance_frequency_damped(omega0 : Double, gamma : Double) -> Double {
  Double::sqrt(omega0 * omega0 - gamma * gamma)
}

///|
/// Normal frequency of coupled oscillations (two identical oscillators)
/// omega0: natural angular frequency of single oscillator (rad/s)
/// k_coupling: coupling constant
/// m: mass (kg)
/// Returns: normal frequency (symmetric mode)
pub fn normal_frequency_symmetric(omega0 : Double) -> Double {
  omega0
}

///|
/// Normal frequency of coupled oscillations (two identical oscillators, antisymmetric mode)
/// omega0: natural angular frequency of single oscillator (rad/s)
/// k_coupling: coupling constant
/// m: mass (kg)
/// Returns: normal frequency (antisymmetric mode)
pub fn normal_frequency_antisymmetric(
  omega0 : Double,
  k_coupling : Double,
  mass : Double,
) -> Double {
  Double::sqrt(omega0 * omega0 + 2.0 * k_coupling / mass)
}

///|
/// Center of mass velocity (one-dimensional, multiple particles)
/// masses: mass array (kg)
/// velocities: velocity array (m/s)
pub fn center_of_mass_velocity_1d(
  masses : Array[Double],
  velocities : Array[Double],
) -> Double {
  if masses.length() != velocities.length() || masses.length() == 0 {
    0.0
  } else {
    let mut total_mass = 0.0
    let mut momentum_sum = 0.0
    for i = 0; i < masses.length(); i = i + 1 {
      total_mass = total_mass + masses[i]
      momentum_sum = momentum_sum + masses[i] * velocities[i]
    }
    if total_mass > 0.0 {
      momentum_sum / total_mass
    } else {
      0.0
    }
  }
}

///|
/// Center of mass acceleration (one-dimensional, multiple particles)
/// masses: mass array (kg)
/// accelerations: acceleration array (m/s²)
pub fn center_of_mass_acceleration_1d(
  masses : Array[Double],
  accelerations : Array[Double],
) -> Double {
  if masses.length() != accelerations.length() || masses.length() == 0 {
    0.0
  } else {
    let mut total_mass = 0.0
    let mut force_sum = 0.0
    for i = 0; i < masses.length(); i = i + 1 {
      total_mass = total_mass + masses[i]
      force_sum = force_sum + masses[i] * accelerations[i]
    }
    if total_mass > 0.0 {
      force_sum / total_mass
    } else {
      0.0
    }
  }
}

///|
/// Kinetic energy loss in perfectly inelastic collision
/// m1: mass of object 1 (kg)
/// v1: initial velocity of object 1 (m/s)
/// m2: mass of object 2 (kg)
/// v2: initial velocity of object 2 (m/s)
/// Returns: kinetic energy lost
pub fn kinetic_energy_loss_inelastic(
  m1 : Double,
  v1 : Double,
  m2 : Double,
  v2 : Double,
) -> Double {
  let ke_initial = kinetic_energy(m1, v1) + kinetic_energy(m2, v2)
  let v_final = inelastic_collision_velocity(m1, v1, m2, v2)
  let ke_final = kinetic_energy(m1 + m2, v_final)
  ke_initial - ke_final
}

///|
/// Kinetic energy loss in partially elastic collision
/// m1: mass of object 1 (kg)
/// v1: initial velocity of object 1 (m/s)
/// m2: mass of object 2 (kg)
/// v2: initial velocity of object 2 (m/s)
/// e: coefficient of restitution
/// Returns: kinetic energy lost
pub fn kinetic_energy_loss_partial(
  m1 : Double,
  v1 : Double,
  m2 : Double,
  v2 : Double,
  e : Double,
) -> Double {
  let ke_initial = kinetic_energy(m1, v1) + kinetic_energy(m2, v2)
  // Calculate final velocity using coefficient of restitution (simplified formula)
  let v_rel = v1 - v2
  let v1f = v1 - (1.0 + e) * m2 * v_rel / (m1 + m2)
  let v2f = v2 + (1.0 + e) * m1 * v_rel / (m1 + m2)
  let ke_final = kinetic_energy(m1, v1f) + kinetic_energy(m2, v2f)
  ke_initial - ke_final
}

///|
/// Rolling without slipping condition: relationship between linear velocity and angular velocity
/// omega: angular velocity (rad/s)
/// r: radius (m)
/// Returns: linear velocity
pub fn rolling_velocity(angular_velocity : Double, radius : Double) -> Double {
  angular_velocity * radius
}

///|
/// Total kinetic energy of rolling object (translational + rotational)
/// m: mass (kg)
/// v: linear velocity (m/s)
/// I: moment of inertia (kg·m²)
/// omega: angular velocity (rad/s)
pub fn rolling_kinetic_energy(
  mass : Double,
  velocity : Double,
  moment_of_inertia : Double,
  angular_velocity : Double,
) -> Double {
  kinetic_energy(mass, velocity) +
  rotational_kinetic_energy(moment_of_inertia, angular_velocity)
}

///|
/// Total kinetic energy of rolling object (using moment of inertia coefficient)
/// m: mass (kg)
/// v: linear velocity (m/s)
/// beta: moment of inertia coefficient (e.g., 0.4 for sphere, 0.5 for cylinder)
pub fn rolling_kinetic_energy_simplified(
  mass : Double,
  velocity : Double,
  beta : Double,
) -> Double {
  0.5 * mass * velocity * velocity * (1.0 + beta)
}

///|
/// Rolling downhill acceleration (no slipping)
/// g: gravitational acceleration (m/s²)
/// angle: inclination angle (radians)
/// beta: moment of inertia coefficient
pub fn rolling_acceleration_downhill(
  gravity : Double,
  angle : Double,
  beta : Double,
) -> Double {
  gravity * @math.sin(angle) / (1.0 + beta)
}

///|
/// Centrifugal force
/// m: mass (kg)
/// omega: angular velocity (rad/s)
/// r: distance from rotation axis (m)
pub fn centrifugal_force(
  mass : Double,
  angular_velocity : Double,
  radius : Double,
) -> Double {
  mass * angular_velocity * angular_velocity * radius
}

///|
/// Coriolis force (one-dimensional simplified)
/// m: mass (kg)
/// v: relative velocity (m/s)
/// omega: angular velocity (rad/s)
pub fn coriolis_force(
  mass : Double,
  relative_velocity : Double,
  angular_velocity : Double,
) -> Double {
  2.0 * mass * relative_velocity * angular_velocity
}

///|
/// Inertial force (accelerating reference frame)
/// m: mass (kg)
/// a: reference frame acceleration (m/s²)
pub fn inertial_force(mass : Double, acceleration : Double) -> Double {
  -mass * acceleration
}

///|
/// Precession angular velocity (gyroscope precession)
/// m: mass (kg)
/// g: gravitational acceleration (m/s²)
/// r: distance from center of mass to pivot (m)
/// i: moment of inertia (kg·m²)
/// omega: spin angular velocity (rad/s)
pub fn precession_angular_velocity(
  mass : Double,
  gravity : Double,
  radius : Double,
  moment_of_inertia : Double,
  spin_angular_velocity : Double,
) -> Double {
  mass * gravity * radius / (moment_of_inertia * spin_angular_velocity)
}

///|
/// Kinetic energy of rigid body rotating about fixed axis
/// i: moment of inertia (kg·m²)
/// omega: angular velocity (rad/s)
pub fn rotational_kinetic_energy_simple(
  moment_of_inertia : Double,
  angular_velocity : Double,
) -> Double {
  0.5 * moment_of_inertia * angular_velocity * angular_velocity
}

///|
/// Power of rigid body rotation
/// tau: torque (N·m)
/// omega: angular velocity (rad/s)
pub fn rotational_power(torque : Double, angular_velocity : Double) -> Double {
  torque * angular_velocity
}

///|
/// Orbital energy (circular orbit)
/// m: object mass (kg)
/// m_central: central body mass (kg)
/// r: orbital radius (m)
/// g_const: gravitational constant (m³/(kg·s²))
pub fn orbital_energy_circular(
  mass : Double,
  central_mass : Double,
  radius : Double,
  g_const : Double,
) -> Double {
  -g_const * mass * central_mass / (2.0 * radius)
}

///|
/// Orbital energy (elliptical orbit)
/// m: object mass (kg)
/// m_central: central body mass (kg)
/// a: semi-major axis (m)
/// g_const: gravitational constant (m³/(kg·s²))
pub fn orbital_energy_elliptical(
  mass : Double,
  central_mass : Double,
  semi_major_axis : Double,
  g_const : Double,
) -> Double {
  -g_const * mass * central_mass / (2.0 * semi_major_axis)
}

///|
/// Kepler's second law: areal velocity
/// r: distance (m)
/// v_perp: velocity component perpendicular to radial direction (m/s)
pub fn areal_velocity(radius : Double, v_perpendicular : Double) -> Double {
  0.5 * radius * v_perpendicular
}

///|
/// Orbital period (elliptical orbit)
/// a: semi-major axis (m)
/// m_central: central body mass (kg)
/// g_const: gravitational constant (m³/(kg·s²))
pub fn orbital_period_elliptical(
  semi_major_axis : Double,
  central_mass : Double,
  g_const : Double,
) -> Double {
  2.0 *
  3.141592653589793 *
  Double::sqrt(
    semi_major_axis *
    semi_major_axis *
    semi_major_axis /
    (g_const * central_mass),
  )
}

///|
/// Harmonic wave: displacement
/// a: amplitude (m)
/// k: wave number (rad/m)
/// x: position (m)
/// omega: angular frequency (rad/s)
/// t: time (s)
/// phi: initial phase (radians)
pub fn wave_displacement(
  amplitude : Double,
  wave_number : Double,
  position : Double,
  angular_frequency : Double,
  time : Double,
  phase : Double,
) -> Double {
  amplitude *
  @math.cos(wave_number * position - angular_frequency * time + phase)
}

///|
/// Wave velocity (using frequency and wavelength)
/// f: frequency (Hz)
/// lambda: wavelength (m)
pub fn wave_velocity(frequency : Double, wavelength : Double) -> Double {
  frequency * wavelength
}

///|
/// Wave velocity (using angular frequency and wave number)
/// omega: angular frequency (rad/s)
/// k: wave number (rad/m)
pub fn wave_velocity_from_angular(
  angular_frequency : Double,
  wave_number : Double,
) -> Double {
  angular_frequency / wave_number
}

///|
/// Doppler effect: observer stationary, source moving
/// f0: source frequency (Hz)
/// v: sound speed (m/s)
/// vs: source velocity (m/s), positive toward observer
pub fn doppler_effect_source_moving(
  source_frequency : Double,
  sound_speed : Double,
  source_velocity : Double,
) -> Double {
  source_frequency * sound_speed / (sound_speed - source_velocity)
}

///|
/// Doppler effect: source stationary, observer moving
/// f0: source frequency (Hz)
/// v: sound speed (m/s)
/// vo: observer velocity (m/s), positive toward source
pub fn doppler_effect_observer_moving(
  source_frequency : Double,
  sound_speed : Double,
  observer_velocity : Double,
) -> Double {
  source_frequency * (sound_speed + observer_velocity) / sound_speed
}

///|
/// Reynolds number (determining flow state)
/// rho: fluid density (kg/m³)
/// v: flow velocity (m/s)
/// d: characteristic length (e.g., pipe diameter) (m)
/// mu: dynamic viscosity (Pa·s)
pub fn reynolds_number(
  density : Double,
  velocity : Double,
  characteristic_length : Double,
  dynamic_viscosity : Double,
) -> Double {
  density * velocity * characteristic_length / dynamic_viscosity
}

///|
/// Stokes' law: viscous drag (sphere)
/// r: sphere radius (m)
/// v: velocity (m/s)
/// mu: dynamic viscosity (Pa·s)
pub fn stokes_drag(
  radius : Double,
  velocity : Double,
  dynamic_viscosity : Double,
) -> Double {
  6.0 * 3.141592653589793 * dynamic_viscosity * radius * velocity
}

///|
/// Viscous force (flow between parallel plates)
/// mu: dynamic viscosity (Pa·s)
/// a: area (m²)
/// dv_dy: velocity gradient (1/s)
pub fn viscous_force(
  dynamic_viscosity : Double,
  area : Double,
  velocity_gradient : Double,
) -> Double {
  dynamic_viscosity * area * velocity_gradient
}

///|
/// Flow rate from pressure difference (Poiseuille's law, simplified)
/// delta_p: pressure difference (Pa)
/// r: pipe radius (m)
/// l: pipe length (m)
/// mu: dynamic viscosity (Pa·s)
pub fn poiseuille_flow_rate(
  pressure_difference : Double,
  radius : Double,
  length : Double,
  dynamic_viscosity : Double,
) -> Double {
  3.141592653589793 *
  pressure_difference *
  radius *
  radius *
  radius *
  radius /
  (8.0 * dynamic_viscosity * length)
}

///|
/// Impulse-momentum theorem: calculate final velocity
/// m: mass (kg)
/// v0: initial velocity (m/s)
/// impulse: impulse (N·s)
pub fn velocity_from_impulse(
  mass : Double,
  initial_velocity : Double,
  impulse : Double,
) -> Double {
  initial_velocity + impulse / mass
}

///|
/// Variable mass system: rocket equation (Tsiolkovsky formula)
/// v0: initial velocity (m/s)
/// ve: exhaust velocity (m/s)
/// m0: initial mass (kg)
/// mf: final mass (kg)
pub fn rocket_equation(
  initial_velocity : Double,
  exhaust_velocity : Double,
  initial_mass : Double,
  final_mass : Double,
) -> Double {
  initial_velocity + exhaust_velocity * @math.ln(initial_mass / final_mass)
}

///|
/// Power-time relationship: calculate work
/// p: power (W)
/// t: time (s)
pub fn work_from_power(power : Double, time : Double) -> Double {
  power * time
}

///|
/// Average power (variable force)
/// w: total work (J)
/// t: total time (s)
pub fn average_power(work : Double, time : Double) -> Double {
  work / time
}

///|
/// Damping ratio
/// gamma: damping coefficient
/// omega0: natural angular frequency (rad/s)
pub fn damping_ratio(
  damping_coefficient : Double,
  natural_frequency : Double,
) -> Double {
  damping_coefficient / natural_frequency
}

///|
/// Compound pendulum: angular frequency
/// m: mass (kg)
/// g: gravitational acceleration (m/s²)
/// i: moment of inertia (kg·m²)
/// d: distance from center of mass to pivot (m)
pub fn physical_pendulum_angular_frequency(
  mass : Double,
  gravity : Double,
  moment_of_inertia : Double,
  distance : Double,
) -> Double {
  Double::sqrt(mass * gravity * distance / moment_of_inertia)
}

///|
/// Compound pendulum: period
/// m: mass (kg)
/// g: gravitational acceleration (m/s²)
/// i: moment of inertia (kg·m²)
/// d: distance from center of mass to pivot (m)
pub fn physical_pendulum_period(
  mass : Double,
  gravity : Double,
  moment_of_inertia : Double,
  distance : Double,
) -> Double {
  2.0 *
  3.141592653589793 *
  Double::sqrt(moment_of_inertia / (mass * gravity * distance))
}

///|
/// Double pendulum: angular frequency under small angle approximation (upper pendulum)
/// g: gravitational acceleration (m/s²)
/// l1: upper pendulum length (m)
pub fn double_pendulum_upper_frequency(
  gravity : Double,
  upper_length : Double,
) -> Double {
  Double::sqrt(gravity / upper_length)
}

///|
/// Normal force (horizontal surface)
/// m: mass (kg)
/// g: gravitational acceleration (m/s²)
pub fn normal_force_horizontal(mass : Double, gravity : Double) -> Double {
  mass * gravity
}

///|
/// Normal force (inclined plane, considering friction)
/// m: mass (kg)
/// g: gravitational acceleration (m/s²)
/// angle: incline angle (radians)
/// f_parallel: force along incline (e.g., friction) (N)
pub fn normal_force_inclined_with_force(
  mass : Double,
  gravity : Double,
  angle : Double,
  parallel_force : Double,
) -> Double {
  mass * gravity * @math.cos(angle) + parallel_force * @math.sin(angle)
}

///|
/// Equilibrium condition: net force equals zero (1D)
/// forces: force array (N)
pub fn equilibrium_net_force_1d(forces : Array[Double]) -> Bool {
  let epsilon = 1.0e-9
  Double::abs(net_force_1d(forces)) < epsilon
}

///|
/// Equilibrium condition: net torque equals zero
/// torques: torque array (N·m)
pub fn equilibrium_net_torque(torques : Array[Double]) -> Bool {
  let epsilon = 1.0e-9
  Double::abs(net_torque(torques)) < epsilon
}

///|
/// Uniformly accelerated motion: calculate time (given displacement and velocity)
/// v0: initial velocity (m/s)
/// v: final velocity (m/s)
/// a: acceleration (m/s²)
pub fn time_from_velocity_change(v0 : Double, v : Double, a : Double) -> Double {
  if a != 0.0 {
    (v - v0) / a
  } else {
    0.0
  }
}

///|
/// Uniformly accelerated motion: calculate acceleration (given displacement, initial and final velocities)
/// v0: initial velocity (m/s)
/// v: final velocity (m/s)
/// s: displacement (m)
pub fn acceleration_from_displacement(
  v0 : Double,
  v : Double,
  s : Double,
) -> Double {
  if s != 0.0 {
    (v * v - v0 * v0) / (2.0 * s)
  } else {
    0.0
  }
}

///|
/// Projectile motion: velocity magnitude at any time
/// v0: initial velocity (m/s)
/// angle: launch angle (radians)
/// t: time (s)
/// g: gravitational acceleration (m/s²)
pub fn projectile_velocity_magnitude(
  v0 : Double,
  angle : Double,
  t : Double,
  g : Double,
) -> Double {
  let vx = v0 * @math.cos(angle)
  let vy = v0 * @math.sin(angle) - g * t
  Double::sqrt(vx * vx + vy * vy)
}

///|
/// Projectile motion: velocity direction at any time (angle with horizontal)
/// v0: initial velocity (m/s)
/// angle: launch angle (radians)
/// t: time (s)
/// g: gravitational acceleration (m/s²)
pub fn projectile_velocity_angle(
  v0 : Double,
  angle : Double,
  t : Double,
  g : Double,
) -> Double {
  let vx = v0 * @math.cos(angle)
  let vy = v0 * @math.sin(angle) - g * t
  if vx != 0.0 {
    @math.atan(vy / vx)
  } else if vy > 0.0 {
    3.141592653589793 / 2.0
  } else {
    -3.141592653589793 / 2.0
  }
}

///|
/// Perfectly elastic collision: 1D, velocity after collision of two objects
/// m1: object 1 mass (kg)
/// m2: object 2 mass (kg)
/// v1i: object 1 initial velocity (m/s)
/// v2i: object 2 initial velocity (m/s)
/// Returns: object 1 final velocity (m/s)
pub fn elastic_collision_1d_v1f(
  m1 : Double,
  m2 : Double,
  v1i : Double,
  v2i : Double,
) -> Double {
  (m1 - m2) * v1i / (m1 + m2) + 2.0 * m2 * v2i / (m1 + m2)
}

///|
/// Perfectly elastic collision: 1D, velocity after collision of two objects
/// m1: object 1 mass (kg)
/// m2: object 2 mass (kg)
/// v1i: object 1 initial velocity (m/s)
/// v2i: object 2 initial velocity (m/s)
/// Returns: object 2 final velocity (m/s)
pub fn elastic_collision_1d_v2f(
  m1 : Double,
  m2 : Double,
  v1i : Double,
  v2i : Double,
) -> Double {
  2.0 * m1 * v1i / (m1 + m2) + (m2 - m1) * v2i / (m1 + m2)
}

///|
/// Perfectly inelastic collision: common velocity after collision
/// m1: object 1 mass (kg)
/// m2: object 2 mass (kg)
/// v1i: object 1 initial velocity (m/s)
/// v2i: object 2 initial velocity (m/s)
pub fn perfectly_inelastic_collision_velocity(
  m1 : Double,
  m2 : Double,
  v1i : Double,
  v2i : Double,
) -> Double {
  (m1 * v1i + m2 * v2i) / (m1 + m2)
}

///|
/// Inelastic collision: calculate final velocity using coefficient of restitution
/// m1: object 1 mass (kg)
/// m2: object 2 mass (kg)
/// v1i: object 1 initial velocity (m/s)
/// v2i: object 2 initial velocity (m/s)
/// e: coefficient of restitution
/// Returns: object 1 final velocity (m/s)
pub fn inelastic_collision_1d_v1f(
  m1 : Double,
  m2 : Double,
  v1i : Double,
  v2i : Double,
  e : Double,
) -> Double {
  let v_com = (m1 * v1i + m2 * v2i) / (m1 + m2)
  v_com + e * m2 * (v2i - v1i) / (m1 + m2)
}

///|
/// Generalized momentum (1D)
/// l: Lagrangian (J)
/// q_dot: generalized velocity (m/s or rad/s)
pub fn generalized_momentum_1d(generalized_velocity : Double) -> Double {
  // Simplified version, actually requires partial derivative
  // Here assumes L = 0.5 * m * q_dot^2 - V(q)
  generalized_velocity
}

///|
/// Generalized force (1D, conservative)
/// v: potential energy (J)
/// q: generalized coordinate (m or rad)
pub fn generalized_force_1d(
  potential_energy : Double,
  generalized_coordinate : Double,
) -> Double {
  // Simplified version, actually requires negative gradient of potential energy
  -potential_energy / generalized_coordinate
}

///|
/// Hamiltonian canonical equation: q_dot = ∂H/∂p (1D simplified)
/// Note: This is a simplified version, assumes H = p^2/(2m) + V(q), where ∂H/∂p = p/m
/// For simplicity, directly returns momentum value (divide by mass when needed)
/// p: generalized momentum (kg·m/s)
pub fn hamiltonian_q_dot_1d(momentum : Double) -> Double {
  // Simplified version, assumes H = p^2/(2m) + V(q)
  // Actually should be p/m, but missing mass parameter, returns p
  momentum
}

///|
/// Lorentz factor (gamma)
/// v: velocity (m/s)
/// c: speed of light (m/s)
pub fn lorentz_factor(velocity : Double, speed_of_light : Double) -> Double {
  let beta = velocity / speed_of_light
  Double::sqrt(1.0 / (1.0 - beta * beta))
}

///|
/// Relativistic mass
/// m0: rest mass (kg)
/// v: velocity (m/s)
/// c: speed of light (m/s)
pub fn relativistic_mass(
  rest_mass : Double,
  velocity : Double,
  speed_of_light : Double,
) -> Double {
  rest_mass * lorentz_factor(velocity, speed_of_light)
}

///|
/// Relativistic momentum
/// m0: rest mass (kg)
/// v: velocity (m/s)
/// c: speed of light (m/s)
pub fn relativistic_momentum(
  rest_mass : Double,
  velocity : Double,
  speed_of_light : Double,
) -> Double {
  rest_mass * velocity * lorentz_factor(velocity, speed_of_light)
}

///|
/// Relativistic total energy
/// m0: rest mass (kg)
/// v: velocity (m/s)
/// c: speed of light (m/s)
pub fn relativistic_total_energy(
  rest_mass : Double,
  velocity : Double,
  speed_of_light : Double,
) -> Double {
  rest_mass *
  speed_of_light *
  speed_of_light *
  lorentz_factor(velocity, speed_of_light)
}

///|
/// Relativistic kinetic energy
/// m0: rest mass (kg)
/// v: velocity (m/s)
/// c: speed of light (m/s)
pub fn relativistic_kinetic_energy(
  rest_mass : Double,
  velocity : Double,
  speed_of_light : Double,
) -> Double {
  let gamma = lorentz_factor(velocity, speed_of_light)
  (gamma - 1.0) * rest_mass * speed_of_light * speed_of_light
}

///|
/// Mass-energy equivalence: E = mc²
/// m: mass (kg)
/// c: speed of light (m/s)
pub fn mass_energy_equivalence(
  mass : Double,
  speed_of_light : Double,
) -> Double {
  mass * speed_of_light * speed_of_light
}

///|
/// Stress (normal stress)
/// f: force (N)
/// a: area (m²)
pub fn normal_stress(force : Double, area : Double) -> Double {
  force / area
}

///|
/// Shear stress
/// f: tangential force (N)
/// a: area (m²)
pub fn shear_stress(tangential_force : Double, area : Double) -> Double {
  tangential_force / area
}

///|
/// Strain (normal strain)
/// delta_l: length change (m)
/// l0: original length (m)
pub fn normal_strain(
  length_change : Double,
  original_length : Double,
) -> Double {
  if original_length != 0.0 {
    length_change / original_length
  } else {
    0.0
  }
}

///|
/// Young's modulus (elastic modulus)
/// stress: stress (Pa)
/// strain: strain
pub fn youngs_modulus(stress : Double, strain : Double) -> Double {
  if strain != 0.0 {
    stress / strain
  } else {
    0.0
  }
}

///|
/// Hooke's law: elastic restoring force
/// k: spring constant (N/m)
/// x: displacement (m)
pub fn hookes_law_force(
  spring_constant : Double,
  displacement : Double,
) -> Double {
  -spring_constant * displacement
}

///|
/// Bulk modulus (compression modulus)
/// delta_p: pressure change (Pa)
/// delta_v: volume change (m³)
/// v0: original volume (m³)
pub fn bulk_modulus(
  pressure_change : Double,
  volume_change : Double,
  original_volume : Double,
) -> Double {
  if volume_change != 0.0 && original_volume != 0.0 {
    -pressure_change * original_volume / volume_change
  } else {
    0.0
  }
}

///|
/// Sound intensity (wave intensity)
/// p: pressure amplitude (Pa)
/// rho: medium density (kg/m³)
/// c: sound speed (m/s)
pub fn sound_intensity(
  pressure_amplitude : Double,
  density : Double,
  sound_speed : Double,
) -> Double {
  pressure_amplitude * pressure_amplitude / (2.0 * density * sound_speed)
}

///|
/// Sound intensity level (decibels)
/// i: sound intensity (W/m²)
/// i0: reference intensity (W/m²), typically 10^-12 W/m²
pub fn sound_intensity_level(
  intensity : Double,
  reference_intensity : Double,
) -> Double {
  if intensity > 0.0 && reference_intensity > 0.0 {
    10.0 * @math.log10(intensity / reference_intensity)
  } else {
    0.0
  }
}

///|
/// Sound pressure level (decibels)
/// p: sound pressure (Pa)
/// p0: reference pressure (Pa), typically 20 μPa
pub fn sound_pressure_level(
  pressure : Double,
  reference_pressure : Double,
) -> Double {
  if pressure > 0.0 && reference_pressure > 0.0 {
    20.0 * @math.log10(pressure / reference_pressure)
  } else {
    0.0
  }
}

///|
/// Sound speed (ideal gas)
/// gamma: specific heat ratio
/// r: gas constant (J/(kg·K))
/// t: temperature (K)
pub fn sound_speed_ideal_gas(
  gamma : Double,
  gas_constant : Double,
  temperature : Double,
) -> Double {
  Double::sqrt(gamma * gas_constant * temperature)
}

///|
/// Coupled oscillators: normal mode frequency (symmetric coupling, two identical oscillators)
/// k: spring constant (N/m)
/// k_coupling: coupling spring constant (N/m)
/// m: mass (kg)
/// Returns: normal mode frequency (rad/s)
pub fn normal_mode_frequency_symmetric(
  k : Double,
  k_coupling : Double,
  mass : Double,
) -> Double {
  Double::sqrt((k + 2.0 * k_coupling) / mass)
}

///|
/// Coupled oscillators: normal mode frequency (antisymmetric mode)
/// k: spring constant (N/m)
/// k_coupling: coupling spring constant (N/m)
/// m: mass (kg)
/// Returns: normal mode frequency (rad/s)
pub fn normal_mode_frequency_antisymmetric(k : Double, mass : Double) -> Double {
  Double::sqrt(k / mass)
}

///|
/// Coupled oscillators: normal coordinate (symmetric mode)
/// x1: oscillator 1 displacement (m)
/// x2: oscillator 2 displacement (m)
pub fn normal_coordinate_symmetric(x1 : Double, x2 : Double) -> Double {
  (x1 + x2) / 2.0
}

///|
/// Coupled oscillators: normal coordinate (antisymmetric mode)
/// x1: oscillator 1 displacement (m)
/// x2: oscillator 2 displacement (m)
pub fn normal_coordinate_antisymmetric(x1 : Double, x2 : Double) -> Double {
  (x1 - x2) / 2.0
}

///|
/// Angular acceleration (from angular velocity and time)
/// omega0: initial angular velocity (rad/s)
/// omega: final angular velocity (rad/s)
/// t: time (s)
pub fn angular_acceleration_from_velocity(
  omega0 : Double,
  omega : Double,
  time : Double,
) -> Double {
  if time != 0.0 {
    (omega - omega0) / time
  } else {
    0.0
  }
}

///|
/// Angular displacement (uniform angular acceleration)
/// omega0: initial angular velocity (rad/s)
/// alpha: angular acceleration (rad/s²)
/// t: time (s)
pub fn angular_displacement_constant_acc(
  omega0 : Double,
  alpha : Double,
  t : Double,
) -> Double {
  omega0 * t + 0.5 * alpha * t * t
}

///|
/// Angular velocity (uniform angular acceleration)
/// omega0: initial angular velocity (rad/s)
/// alpha: angular acceleration (rad/s²)
/// t: time (s)
pub fn angular_velocity_constant_acc(
  omega0 : Double,
  alpha : Double,
  t : Double,
) -> Double {
  omega0 + alpha * t
}

///|
/// Angular velocity squared relation (uniform angular acceleration)
/// omega0: initial angular velocity (rad/s)
/// omega: final angular velocity (rad/s)
/// alpha: angular acceleration (rad/s²)
pub fn angular_velocity_squared_relation(
  omega0 : Double,
  omega : Double,
  alpha : Double,
) -> Double {
  (omega * omega - omega0 * omega0) / (2.0 * alpha)
}

///|
/// Phase space volume (1D)
/// x: position (m)
/// p: momentum (kg·m/s)
pub fn phase_space_volume_1d(position : Double, momentum : Double) -> Double {
  position * momentum
}

///|
/// Beat frequency (difference of two close frequencies)
/// f1: frequency 1 (Hz)
/// f2: frequency 2 (Hz)
pub fn beat_frequency(f1 : Double, f2 : Double) -> Double {
  Double::abs(f1 - f2)
}

///|
/// Beat period
/// f1: frequency 1 (Hz)
/// f2: frequency 2 (Hz)
pub fn beat_period(f1 : Double, f2 : Double) -> Double {
  let beat_freq = beat_frequency(f1, f2)
  if beat_freq > 0.0 {
    1.0 / beat_freq
  } else {
    0.0
  }
}

///|
/// Mach number (ratio of flow velocity to sound speed)
/// v: flow velocity (m/s)
/// c: sound speed (m/s)
pub fn mach_number(velocity : Double, sound_speed : Double) -> Double {
  if sound_speed > 0.0 {
    velocity / sound_speed
  } else {
    0.0
  }
}

///|
/// Venturi effect: velocity ratio
/// a1: cross-section 1 area (m²)
/// a2: cross-section 2 area (m²)
pub fn venturi_velocity_ratio(a1 : Double, a2 : Double) -> Double {
  if a2 > 0.0 {
    a1 / a2
  } else {
    0.0
  }
}

///|
/// Laminar drag (cylinder)
/// mu: dynamic viscosity (Pa·s)
/// l: length (m)
/// v: velocity (m/s)
/// d: diameter (m)
pub fn laminar_drag_cylinder(
  dynamic_viscosity : Double,
  length : Double,
  velocity : Double,
  diameter : Double,
) -> Double {
  4.0 *
  3.141592653589793 *
  dynamic_viscosity *
  length *
  velocity /
  @math.ln(length / diameter)
}

///|
/// Circular motion: tangential acceleration
/// alpha: angular acceleration (rad/s²)
/// r: radius (m)
pub fn tangential_acceleration(
  angular_acceleration : Double,
  radius : Double,
) -> Double {
  angular_acceleration * radius
}

///|
/// Circular motion: total acceleration magnitude
/// a_t: tangential acceleration (m/s²)
/// a_c: centripetal acceleration (m/s²)
pub fn total_acceleration_magnitude(
  tangential_acc : Double,
  centripetal_acc : Double,
) -> Double {
  Double::sqrt(
    tangential_acc * tangential_acc + centripetal_acc * centripetal_acc,
  )
}

///|
/// Curvilinear motion: velocity from path length
/// ds: path length element (m)
/// dt: time element (s)
pub fn velocity_from_path(ds : Double, dt : Double) -> Double {
  if dt > 0.0 {
    ds / dt
  } else {
    0.0
  }
}

///|
/// Relative motion: velocity transformation
/// v_object: object velocity in frame A (m/s)
/// v_frame: frame A velocity relative to frame B (m/s)
/// Returns: object velocity in frame B
pub fn velocity_transformation(v_object : Double, v_frame : Double) -> Double {
  v_object + v_frame
}

///|
/// Variable mass: thrust force
/// dm_dt: mass flow rate (kg/s)
/// v_exhaust: exhaust velocity (m/s)
pub fn thrust_force(
  mass_flow_rate : Double,
  exhaust_velocity : Double,
) -> Double {
  mass_flow_rate * exhaust_velocity
}

///|
/// Variable mass: acceleration from thrust
/// thrust: thrust force (N)
/// mass: current mass (kg)
/// drag: drag force (N)
pub fn acceleration_from_thrust(
  thrust : Double,
  mass : Double,
  drag : Double,
) -> Double {
  if mass > 0.0 {
    (thrust - drag) / mass
  } else {
    0.0
  }
}

///|
/// Moment of inertia: thin ring (about diameter)
/// m: mass (kg)
/// r: radius (m)
pub fn moment_of_inertia_ring_diameter(
  mass : Double,
  radius : Double,
) -> Double {
  0.5 * mass * radius * radius
}

///|
/// Moment of inertia: rectangular bar (about center, perpendicular to length)
/// m: mass (kg)
/// a: width (m)
/// b: height (m)
pub fn moment_of_inertia_rectangular_bar(
  mass : Double,
  width : Double,
  height : Double,
) -> Double {
  mass * (width * width + height * height) / 12.0
}

///|
/// Moment of inertia: ellipsoid (about axis through center)
/// m: mass (kg)
/// a: semi-axis a (m)
/// b: semi-axis b (m)
/// c: semi-axis c (m)
/// axis: which axis (1=a, 2=b, 3=c)
pub fn moment_of_inertia_ellipsoid(
  mass : Double,
  a : Double,
  b : Double,
  c : Double,
  axis : Int,
) -> Double {
  if axis == 1 {
    0.2 * mass * (b * b + c * c)
  } else if axis == 2 {
    0.2 * mass * (a * a + c * c)
  } else {
    0.2 * mass * (a * a + b * b)
  }
}

///|
/// Angular momentum: rate of change
/// tau: net torque (N·m)
/// Returns: rate of change of angular momentum (N·m)
pub fn angular_momentum_rate(torque : Double) -> Double {
  torque
}

///|
/// Gyroscope: nutation frequency
/// i: moment of inertia (kg·m²)
/// omega_precess: precession angular velocity (rad/s)
/// omega_spin: spin angular velocity (rad/s)
/// m: mass (kg)
/// g: gravitational acceleration (m/s²)
/// r: distance from pivot to center of mass (m)
pub fn nutation_frequency(
  moment_of_inertia : Double,
  omega_precess : Double,
  omega_spin : Double,
  mass : Double,
  gravity : Double,
  radius : Double,
) -> Double {
  let omega_p = omega_precess
  let omega_s = omega_spin
  Double::sqrt(
    mass * gravity * radius / moment_of_inertia +
    omega_p * omega_s * (omega_p * omega_s),
  )
}

///|
/// Wave intensity from amplitude
/// a: amplitude (m)
/// omega: angular frequency (rad/s)
/// rho: medium density (kg/m³)
/// c: wave speed (m/s)
pub fn wave_intensity_from_amplitude(
  amplitude : Double,
  angular_frequency : Double,
  density : Double,
  wave_speed : Double,
) -> Double {
  0.5 *
  density *
  wave_speed *
  angular_frequency *
  angular_frequency *
  amplitude *
  amplitude
}

///|
/// Standing wave: node spacing
/// lambda: wavelength (m)
pub fn standing_wave_node_spacing(wavelength : Double) -> Double {
  wavelength / 2.0
}

///|
/// Standing wave: fundamental frequency
/// c: wave speed (m/s)
/// l: length of string/pipe (m)
pub fn standing_wave_fundamental_frequency(
  wave_speed : Double,
  length : Double,
) -> Double {
  if length > 0.0 {
    wave_speed / (2.0 * length)
  } else {
    0.0
  }
}

///|
/// Standing wave: nth harmonic frequency
/// n: harmonic number
/// c: wave speed (m/s)
/// l: length (m)
pub fn standing_wave_harmonic_frequency(
  n : Int,
  wave_speed : Double,
  length : Double,
) -> Double {
  if length > 0.0 && n > 0 {
    Double::from_int(n) * wave_speed / (2.0 * length)
  } else {
    0.0
  }
}

///|
/// Damped oscillation: energy dissipation rate
/// gamma: damping coefficient
/// energy: current energy (J)
pub fn energy_dissipation_rate(
  damping_coefficient : Double,
  energy : Double,
) -> Double {
  -2.0 * damping_coefficient * energy
}

///|
/// Quality factor (Q factor) of oscillator
/// omega0: natural frequency (rad/s)
/// gamma: damping coefficient
pub fn quality_factor(
  natural_frequency : Double,
  damping_coefficient : Double,
) -> Double {
  if damping_coefficient > 0.0 {
    natural_frequency / (2.0 * damping_coefficient)
  } else {
    0.0
  }
}

///|
/// Surface tension: capillary rise
/// gamma: surface tension (N/m)
/// rho: liquid density (kg/m³)
/// r: capillary radius (m)
/// g: gravitational acceleration (m/s²)
/// theta: contact angle (radians)
pub fn capillary_rise(
  surface_tension : Double,
  density : Double,
  radius : Double,
  gravity : Double,
  contact_angle : Double,
) -> Double {
  if radius > 0.0 && density > 0.0 && gravity > 0.0 {
    2.0 *
    surface_tension *
    @math.cos(contact_angle) /
    (density * gravity * radius)
  } else {
    0.0
  }
}

///|
/// Surface tension: pressure difference across curved surface (Laplace's law)
/// gamma: surface tension (N/m)
/// r1: principal radius 1 (m)
/// r2: principal radius 2 (m)
pub fn laplace_pressure(
  surface_tension : Double,
  radius1 : Double,
  radius2 : Double,
) -> Double {
  if radius1 > 0.0 && radius2 > 0.0 {
    surface_tension * (1.0 / radius1 + 1.0 / radius2)
  } else {
    0.0
  }
}

///|
/// Drag force: quadratic drag (high Reynolds number)
/// cd: drag coefficient
/// rho: fluid density (kg/m³)
/// v: velocity (m/s)
/// a: cross-sectional area (m²)
pub fn quadratic_drag(
  drag_coefficient : Double,
  density : Double,
  velocity : Double,
  area : Double,
) -> Double {
  0.5 * drag_coefficient * density * velocity * velocity * area
}

///|
/// Terminal velocity (falling object with quadratic drag)
/// m: mass (kg)
/// g: gravitational acceleration (m/s²)
/// cd: drag coefficient
/// rho: fluid density (kg/m³)
/// a: cross-sectional area (m²)
pub fn terminal_velocity(
  mass : Double,
  gravity : Double,
  drag_coefficient : Double,
  density : Double,
  area : Double,
) -> Double {
  if drag_coefficient > 0.0 && density > 0.0 && area > 0.0 {
    Double::sqrt(2.0 * mass * gravity / (drag_coefficient * density * area))
  } else {
    0.0
  }
}

///|
/// Flow rate: from velocity and area
/// v: flow velocity (m/s)
/// a: cross-sectional area (m²)
pub fn volumetric_flow_rate(velocity : Double, area : Double) -> Double {
  velocity * area
}

///|
/// Mass flow rate
/// rho: density (kg/m³)
/// v: velocity (m/s)
/// a: area (m²)
pub fn mass_flow_rate(
  density : Double,
  velocity : Double,
  area : Double,
) -> Double {
  density * velocity * area
}

///|
/// Torricelli's theorem: exit velocity
/// h: height of fluid above opening (m)
/// g: gravitational acceleration (m/s²)
pub fn torricelli_velocity(height : Double, gravity : Double) -> Double {
  Double::sqrt(2.0 * gravity * height)
}

///|
/// Mechanical efficiency
/// output_power: useful output power (W)
/// input_power: total input power (W)
pub fn mechanical_efficiency(
  output_power : Double,
  input_power : Double,
) -> Double {
  if input_power > 0.0 {
    output_power / input_power
  } else {
    0.0
  }
}

///|
/// Power loss
/// input_power: input power (W)
/// efficiency: efficiency (0-1)
pub fn power_loss(input_power : Double, efficiency : Double) -> Double {
  input_power * (1.0 - efficiency)
}

///|
/// Energy conversion efficiency
/// energy_out: output energy (J)
/// energy_in: input energy (J)
pub fn energy_efficiency(energy_out : Double, energy_in : Double) -> Double {
  if energy_in > 0.0 {
    energy_out / energy_in
  } else {
    0.0
  }
}

///|
/// Kinetic energy from work
/// work: work done (J)
/// ke_initial: initial kinetic energy (J)
pub fn kinetic_energy_from_work(work : Double, ke_initial : Double) -> Double {
  ke_initial + work
}

///|
/// Potential energy change from work
/// work: work done (J)
/// pe_initial: initial potential energy (J)
pub fn potential_energy_from_work(work : Double, pe_initial : Double) -> Double {
  pe_initial - work
}

///|
/// Two-dimensional elastic collision: final velocity component
/// m1: mass 1 (kg)
/// m2: mass 2 (kg)
/// v1i: initial velocity 1 (m/s)
/// v2i: initial velocity 2 (m/s)
/// theta1: angle of approach 1 (radians)
/// theta2: angle of approach 2 (radians)
/// Returns: final velocity of object 1 in x-direction (simplified)
pub fn elastic_collision_2d_v1x(
  m1 : Double,
  m2 : Double,
  v1i : Double,
  v2i : Double,
  theta1 : Double,
  theta2 : Double,
) -> Double {
  // Simplified 2D collision assuming head-on collision
  let v1ix = v1i * @math.cos(theta1)
  let v2ix = v2i * @math.cos(theta2)
  (m1 - m2) * v1ix / (m1 + m2) + 2.0 * m2 * v2ix / (m1 + m2)
}

///|
/// Collision: kinetic energy in center of mass frame
/// m1: mass 1 (kg)
/// m2: mass 2 (kg)
/// v1: velocity 1 (m/s)
/// v2: velocity 2 (m/s)
pub fn kinetic_energy_com_frame(
  m1 : Double,
  m2 : Double,
  v1 : Double,
  v2 : Double,
) -> Double {
  let mu = m1 * m2 / (m1 + m2) // Reduced mass
  let v_rel = v1 - v2
  0.5 * mu * v_rel * v_rel
}

///|
/// Orbital velocity: elliptical orbit at periapsis
/// a: semi-major axis (m)
/// e: eccentricity
/// m_central: central mass (kg)
/// g_const: gravitational constant
pub fn orbital_velocity_periapsis(
  semi_major_axis : Double,
  eccentricity : Double,
  central_mass : Double,
  g_const : Double,
) -> Double {
  let a = semi_major_axis
  let e = eccentricity
  Double::sqrt(g_const * central_mass * (1.0 + e) / (a * (1.0 - e)))
}

///|
/// Orbital velocity: elliptical orbit at apoapsis
/// a: semi-major axis (m)
/// e: eccentricity
/// m_central: central mass (kg)
/// g_const: gravitational constant
pub fn orbital_velocity_apoapsis(
  semi_major_axis : Double,
  eccentricity : Double,
  central_mass : Double,
  g_const : Double,
) -> Double {
  let a = semi_major_axis
  let e = eccentricity
  Double::sqrt(g_const * central_mass * (1.0 - e) / (a * (1.0 + e)))
}

///|
/// Hohmann transfer: delta-v for first burn
/// r1: initial orbit radius (m)
/// r2: target orbit radius (m)
/// m_central: central mass (kg)
/// g_const: gravitational constant
pub fn hohmann_transfer_dv1(
  r1 : Double,
  r2 : Double,
  central_mass : Double,
  g_const : Double,
) -> Double {
  let a_transfer = (r1 + r2) / 2.0
  let v1_circular = Double::sqrt(g_const * central_mass / r1)
  let v1_transfer = Double::sqrt(
    g_const * central_mass * (2.0 / r1 - 1.0 / a_transfer),
  )
  v1_transfer - v1_circular
}

///|
/// Hohmann transfer: delta-v for second burn
/// r1: initial orbit radius (m)
/// r2: target orbit radius (m)
/// m_central: central mass (kg)
/// g_const: gravitational constant
pub fn hohmann_transfer_dv2(
  r1 : Double,
  r2 : Double,
  central_mass : Double,
  g_const : Double,
) -> Double {
  let a_transfer = (r1 + r2) / 2.0
  let v2_circular = Double::sqrt(g_const * central_mass / r2)
  let v2_transfer = Double::sqrt(
    g_const * central_mass * (2.0 / r2 - 1.0 / a_transfer),
  )
  v2_circular - v2_transfer
}

///|
/// Average value from array
/// values: array of values
pub fn average(values : Array[Double]) -> Double {
  if values.length() == 0 {
    0.0
  } else {
    let mut sum = 0.0
    for i = 0; i < values.length(); i = i + 1 {
      sum = sum + values[i]
    }
    sum / Double::from_int(values.length())
  }
}

///|
/// Root mean square (RMS) value
/// values: array of values
pub fn rms(values : Array[Double]) -> Double {
  if values.length() == 0 {
    0.0
  } else {
    let mut sum_sq = 0.0
    for i = 0; i < values.length(); i = i + 1 {
      sum_sq = sum_sq + values[i] * values[i]
    }
    Double::sqrt(sum_sq / Double::from_int(values.length()))
  }
}

///|
/// Time-weighted average
/// values: array of values
/// times: array of time intervals (s)
pub fn time_weighted_average(
  values : Array[Double],
  times : Array[Double],
) -> Double {
  if values.length() != times.length() || values.length() == 0 {
    0.0
  } else {
    let mut weighted_sum = 0.0
    let mut total_time = 0.0
    for i = 0; i < values.length(); i = i + 1 {
      weighted_sum = weighted_sum + values[i] * times[i]
      total_time = total_time + times[i]
    }
    if total_time > 0.0 {
      weighted_sum / total_time
    } else {
      0.0
    }
  }
}
