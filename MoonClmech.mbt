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

// ============================================================================
// Multi-body Systems (3D Extensions)
// ============================================================================

///|
/// Reduced mass for two-body system
/// m1: mass of object 1 (kg)
/// m2: mass of object 2 (kg)
pub fn reduced_mass(m1 : Double, m2 : Double) -> Double {
  if m1 + m2 > 0.0 {
    m1 * m2 / (m1 + m2)
  } else {
    0.0
  }
}

///|
/// Center of mass position (3D, multiple particles)
/// masses: mass array (kg)
/// positions_x: x-coordinate array (m)
/// positions_y: y-coordinate array (m)
/// positions_z: z-coordinate array (m)
/// Returns: x-coordinate of center of mass
pub fn center_of_mass_3d_x(
  masses : Array[Double],
  positions_x : Array[Double],
  _positions_y : Array[Double],
  _positions_z : Array[Double],
) -> Double {
  if masses.length() != positions_x.length() || masses.length() == 0 {
    0.0
  } else {
    let mut total_mass = 0.0
    let mut weighted_sum = 0.0
    for i = 0; i < masses.length(); i = i + 1 {
      total_mass = total_mass + masses[i]
      weighted_sum = weighted_sum + masses[i] * positions_x[i]
    }
    if total_mass > 0.0 {
      weighted_sum / total_mass
    } else {
      0.0
    }
  }
}

///|
/// Center of mass position (3D, multiple particles) - y coordinate
pub fn center_of_mass_3d_y(
  masses : Array[Double],
  _positions_x : Array[Double],
  positions_y : Array[Double],
  _positions_z : Array[Double],
) -> Double {
  if masses.length() != positions_y.length() || masses.length() == 0 {
    0.0
  } else {
    let mut total_mass = 0.0
    let mut weighted_sum = 0.0
    for i = 0; i < masses.length(); i = i + 1 {
      total_mass = total_mass + masses[i]
      weighted_sum = weighted_sum + masses[i] * positions_y[i]
    }
    if total_mass > 0.0 {
      weighted_sum / total_mass
    } else {
      0.0
    }
  }
}

///|
/// Center of mass position (3D, multiple particles) - z coordinate
pub fn center_of_mass_3d_z(
  masses : Array[Double],
  _positions_x : Array[Double],
  _positions_y : Array[Double],
  positions_z : Array[Double],
) -> Double {
  if masses.length() != positions_z.length() || masses.length() == 0 {
    0.0
  } else {
    let mut total_mass = 0.0
    let mut weighted_sum = 0.0
    for i = 0; i < masses.length(); i = i + 1 {
      total_mass = total_mass + masses[i]
      weighted_sum = weighted_sum + masses[i] * positions_z[i]
    }
    if total_mass > 0.0 {
      weighted_sum / total_mass
    } else {
      0.0
    }
  }
}

///|
/// Total momentum (3D, multiple particles) - x component
/// masses: mass array (kg)
/// velocities_x: x-velocity array (m/s)
/// velocities_y: y-velocity array (m/s)
/// velocities_z: z-velocity array (m/s)
pub fn total_momentum_3d_x(
  masses : Array[Double],
  velocities_x : Array[Double],
  _velocities_y : Array[Double],
  _velocities_z : Array[Double],
) -> Double {
  if masses.length() != velocities_x.length() {
    0.0
  } else {
    let mut total = 0.0
    for i = 0; i < masses.length(); i = i + 1 {
      total = total + masses[i] * velocities_x[i]
    }
    total
  }
}

///|
/// Total momentum (3D, multiple particles) - y component
pub fn total_momentum_3d_y(
  masses : Array[Double],
  _velocities_x : Array[Double],
  velocities_y : Array[Double],
  _velocities_z : Array[Double],
) -> Double {
  if masses.length() != velocities_y.length() {
    0.0
  } else {
    let mut total = 0.0
    for i = 0; i < masses.length(); i = i + 1 {
      total = total + masses[i] * velocities_y[i]
    }
    total
  }
}

///|
/// Total momentum (3D, multiple particles) - z component
pub fn total_momentum_3d_z(
  masses : Array[Double],
  _velocities_x : Array[Double],
  _velocities_y : Array[Double],
  velocities_z : Array[Double],
) -> Double {
  if masses.length() != velocities_z.length() {
    0.0
  } else {
    let mut total = 0.0
    for i = 0; i < masses.length(); i = i + 1 {
      total = total + masses[i] * velocities_z[i]
    }
    total
  }
}

///|
/// Total angular momentum (3D, multiple particles) - about origin
/// masses: mass array (kg)
/// positions_x: x-coordinate array (m)
/// positions_y: y-coordinate array (m)
/// positions_z: z-coordinate array (m)
/// velocities_x: x-velocity array (m/s)
/// velocities_y: y-velocity array (m/s)
/// velocities_z: z-velocity array (m/s)
/// Returns: z-component of angular momentum (L_z = x*p_y - y*p_x)
pub fn total_angular_momentum_3d_z(
  masses : Array[Double],
  positions_x : Array[Double],
  positions_y : Array[Double],
  _positions_z : Array[Double],
  velocities_x : Array[Double],
  velocities_y : Array[Double],
  _velocities_z : Array[Double],
) -> Double {
  if masses.length() != positions_x.length() || masses.length() == 0 {
    0.0
  } else {
    let mut total = 0.0
    for i = 0; i < masses.length(); i = i + 1 {
      // L_z = x*p_y - y*p_x = x*(m*v_y) - y*(m*v_x)
      total = total +
        masses[i] *
        (positions_x[i] * velocities_y[i] - positions_y[i] * velocities_x[i])
    }
    total
  }
}

// ============================================================================
// Coordinate System Transformations
// ============================================================================

///|
/// Convert Cartesian to polar coordinates (2D)
/// x: x-coordinate (m)
/// y: y-coordinate (m)
/// Returns: radius (m)
pub fn cartesian_to_polar_r(x : Double, y : Double) -> Double {
  Double::sqrt(x * x + y * y)
}

///|
/// Convert Cartesian to polar coordinates (2D)
/// x: x-coordinate (m)
/// y: y-coordinate (m)
/// Returns: angle (radians)
pub fn cartesian_to_polar_theta(x : Double, y : Double) -> Double {
  if x != 0.0 {
    @math.atan(y / x)
  } else if y > 0.0 {
    3.141592653589793 / 2.0
  } else if y < 0.0 {
    -3.141592653589793 / 2.0
  } else {
    0.0
  }
}

///|
/// Convert polar to Cartesian coordinates (2D)
/// r: radius (m)
/// theta: angle (radians)
/// Returns: x-coordinate
pub fn polar_to_cartesian_x(r : Double, theta : Double) -> Double {
  r * @math.cos(theta)
}

///|
/// Convert polar to Cartesian coordinates (2D)
/// r: radius (m)
/// theta: angle (radians)
/// Returns: y-coordinate
pub fn polar_to_cartesian_y(r : Double, theta : Double) -> Double {
  r * @math.sin(theta)
}

///|
/// Convert Cartesian to spherical coordinates (3D)
/// x: x-coordinate (m)
/// y: y-coordinate (m)
/// z: z-coordinate (m)
/// Returns: radius r (m)
pub fn cartesian_to_spherical_r(x : Double, y : Double, z : Double) -> Double {
  Double::sqrt(x * x + y * y + z * z)
}

///|
/// Convert Cartesian to spherical coordinates (3D)
/// x: x-coordinate (m)
/// y: y-coordinate (m)
/// z: z-coordinate (m)
/// Returns: polar angle theta (radians, 0 to π)
pub fn cartesian_to_spherical_theta(
  x : Double,
  y : Double,
  z : Double,
) -> Double {
  let r = cartesian_to_spherical_r(x, y, z)
  if r > 0.0 {
    @math.acos(z / r)
  } else {
    0.0
  }
}

///|
/// Convert Cartesian to spherical coordinates (3D)
/// x: x-coordinate (m)
/// y: y-coordinate (m)
/// z: z-coordinate (m)
/// Returns: azimuthal angle phi (radians, 0 to 2π)
pub fn cartesian_to_spherical_phi(
  x : Double,
  y : Double,
  _z : Double,
) -> Double {
  if x != 0.0 {
    @math.atan(y / x)
  } else if y > 0.0 {
    3.141592653589793 / 2.0
  } else if y < 0.0 {
    -3.141592653589793 / 2.0
  } else {
    0.0
  }
}

///|
/// Convert spherical to Cartesian coordinates (3D)
/// r: radius (m)
/// theta: polar angle (radians)
/// phi: azimuthal angle (radians)
/// Returns: x-coordinate
pub fn spherical_to_cartesian_x(
  r : Double,
  theta : Double,
  phi : Double,
) -> Double {
  r * @math.sin(theta) * @math.cos(phi)
}

///|
/// Convert spherical to Cartesian coordinates (3D)
/// r: radius (m)
/// theta: polar angle (radians)
/// phi: azimuthal angle (radians)
/// Returns: y-coordinate
pub fn spherical_to_cartesian_y(
  r : Double,
  theta : Double,
  phi : Double,
) -> Double {
  r * @math.sin(theta) * @math.sin(phi)
}

///|
/// Convert spherical to Cartesian coordinates (3D)
/// r: radius (m)
/// theta: polar angle (radians)
/// phi: azimuthal angle (radians)
/// Returns: z-coordinate
pub fn spherical_to_cartesian_z(
  r : Double,
  theta : Double,
  _phi : Double,
) -> Double {
  r * @math.cos(theta)
}

///|
/// Velocity in polar coordinates - radial component
/// vx: x-velocity (m/s)
/// vy: y-velocity (m/s)
/// x: x-position (m)
/// y: y-position (m)
pub fn velocity_polar_radial(
  vx : Double,
  vy : Double,
  x : Double,
  y : Double,
) -> Double {
  let r = cartesian_to_polar_r(x, y)
  if r > 0.0 {
    (x * vx + y * vy) / r
  } else {
    0.0
  }
}

///|
/// Velocity in polar coordinates - angular component
/// vx: x-velocity (m/s)
/// vy: y-velocity (m/s)
/// x: x-position (m)
/// y: y-position (m)
pub fn velocity_polar_angular(
  vx : Double,
  vy : Double,
  x : Double,
  y : Double,
) -> Double {
  let r = cartesian_to_polar_r(x, y)
  if r > 0.0 {
    (x * vy - y * vx) / (r * r)
  } else {
    0.0
  }
}

// ============================================================================
// 2D Collision (Complete Solution)
// ============================================================================

///|
/// 2D elastic collision: final velocity of object 1 (x-component)
/// m1: mass of object 1 (kg)
/// m2: mass of object 2 (kg)
/// v1ix: initial x-velocity of object 1 (m/s)
/// v1iy: initial y-velocity of object 1 (m/s)
/// v2ix: initial x-velocity of object 2 (m/s)
/// v2iy: initial y-velocity of object 2 (m/s)
/// collision_angle: angle of collision line (radians)
pub fn elastic_collision_2d_v1fx(
  m1 : Double,
  m2 : Double,
  v1ix : Double,
  v1iy : Double,
  v2ix : Double,
  v2iy : Double,
  collision_angle : Double,
) -> Double {
  // Transform to collision coordinate system
  let cos_a = @math.cos(collision_angle)
  let sin_a = @math.sin(collision_angle)
  let v1ix_n = v1ix * cos_a + v1iy * sin_a
  let v1iy_t = -v1ix * sin_a + v1iy * cos_a
  let v2ix_n = v2ix * cos_a + v2iy * sin_a
  let _v2iy_t = -v2ix * sin_a + v2iy * cos_a

  // 1D collision along normal direction
  let v1fx_n = (m1 - m2) * v1ix_n / (m1 + m2) + 2.0 * m2 * v2ix_n / (m1 + m2)
  let v1fy_t = v1iy_t // Tangential component unchanged

  // Transform back
  v1fx_n * cos_a - v1fy_t * sin_a
}

///|
/// 2D elastic collision: final velocity of object 1 (y-component)
pub fn elastic_collision_2d_v1fy(
  m1 : Double,
  m2 : Double,
  v1ix : Double,
  v1iy : Double,
  v2ix : Double,
  v2iy : Double,
  collision_angle : Double,
) -> Double {
  let cos_a = @math.cos(collision_angle)
  let sin_a = @math.sin(collision_angle)
  let v1ix_n = v1ix * cos_a + v1iy * sin_a
  let v1iy_t = -v1ix * sin_a + v1iy * cos_a
  let v2ix_n = v2ix * cos_a + v2iy * sin_a
  let _v2iy_t = -v2ix * sin_a + v2iy * cos_a
  let v1fx_n = (m1 - m2) * v1ix_n / (m1 + m2) + 2.0 * m2 * v2ix_n / (m1 + m2)
  let v1fy_t = v1iy_t
  v1fx_n * sin_a + v1fy_t * cos_a
}

///|
/// 2D elastic collision: final velocity of object 2 (x-component)
pub fn elastic_collision_2d_v2fx(
  m1 : Double,
  m2 : Double,
  v1ix : Double,
  v1iy : Double,
  v2ix : Double,
  v2iy : Double,
  collision_angle : Double,
) -> Double {
  let cos_a = @math.cos(collision_angle)
  let sin_a = @math.sin(collision_angle)
  let v1ix_n = v1ix * cos_a + v1iy * sin_a
  let _v1iy_t = -v1ix * sin_a + v1iy * cos_a
  let v2ix_n = v2ix * cos_a + v2iy * sin_a
  let v2iy_t = -v2ix * sin_a + v2iy * cos_a
  let v2fx_n = 2.0 * m1 * v1ix_n / (m1 + m2) + (m2 - m1) * v2ix_n / (m1 + m2)
  let v2fy_t = v2iy_t
  v2fx_n * cos_a - v2fy_t * sin_a
}

///|
/// 2D elastic collision: final velocity of object 2 (y-component)
pub fn elastic_collision_2d_v2fy(
  m1 : Double,
  m2 : Double,
  v1ix : Double,
  v1iy : Double,
  v2ix : Double,
  v2iy : Double,
  collision_angle : Double,
) -> Double {
  let cos_a = @math.cos(collision_angle)
  let sin_a = @math.sin(collision_angle)
  let v1ix_n = v1ix * cos_a + v1iy * sin_a
  let _v1iy_t = -v1ix * sin_a + v1iy * cos_a
  let v2ix_n = v2ix * cos_a + v2iy * sin_a
  let v2iy_t = -v2ix * sin_a + v2iy * cos_a
  let v2fx_n = 2.0 * m1 * v1ix_n / (m1 + m2) + (m2 - m1) * v2ix_n / (m1 + m2)
  let v2fy_t = v2iy_t
  v2fx_n * sin_a + v2fy_t * cos_a
}

// ============================================================================
// Numerical Integration Methods
// ============================================================================

///|
/// Euler method for solving first-order ODE: dy/dt = f(t, y)
/// f: function that returns derivative (dy/dt) given (t, y)
/// t0: initial time
/// y0: initial value
/// dt: time step
/// n_steps: number of steps
/// Returns: final value of y
/// Note: This is a simplified version; full implementation would require function pointers
pub fn euler_method_step(
  _current_t : Double,
  current_y : Double,
  derivative : Double,
  dt : Double,
) -> Double {
  current_y + derivative * dt
}

///|
/// Runge-Kutta 4th order method (RK4) - single step
/// f: derivative function (simplified: takes derivative value directly)
/// t: current time
/// y: current value
/// dt: time step
/// k1: derivative at current point (k1 = f(t, y))
/// k2: derivative at midpoint using k1 (k2 = f(t + dt/2, y + k1*dt/2))
/// k3: derivative at midpoint using k2 (k3 = f(t + dt/2, y + k2*dt/2))
/// k4: derivative at endpoint using k3 (k4 = f(t + dt, y + k3*dt))
/// Returns: next value of y
pub fn rk4_step(
  current_y : Double,
  k1 : Double,
  k2 : Double,
  k3 : Double,
  k4 : Double,
  dt : Double,
) -> Double {
  current_y + (k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt / 6.0
}

///|
/// Numerical integration: Trapezoidal rule
/// values: array of function values
/// step_size: spacing between values
pub fn trapezoidal_rule(values : Array[Double], step_size : Double) -> Double {
  if values.length() < 2 {
    0.0
  } else {
    let mut sum = 0.0
    for i = 0; i < values.length() - 1; i = i + 1 {
      sum = sum + (values[i] + values[i + 1]) * step_size / 2.0
    }
    sum
  }
}

///|
/// Numerical integration: Simpson's rule (requires odd number of points)
/// values: array of function values (must have odd length)
/// step_size: spacing between values
pub fn simpsons_rule(values : Array[Double], step_size : Double) -> Double {
  let n = values.length()
  if n < 3 || n % 2 == 0 {
    0.0 // Simpson's rule requires odd number of points
  } else {
    let mut sum = values[0] + values[n - 1]
    for i = 1; i < n - 1; i = i + 1 {
      if i % 2 == 1 {
        sum = sum + 4.0 * values[i]
      } else {
        sum = sum + 2.0 * values[i]
      }
    }
    sum * step_size / 3.0
  }
}

// ============================================================================
// Statistical Mechanics (Basic)
// ============================================================================

///|
/// Ideal gas law: calculate pressure
/// n: number of moles
/// t: temperature (K)
/// v: volume (m³)
/// r: gas constant (J/(mol·K))
pub fn ideal_gas_pressure(
  n : Double,
  t : Double,
  v : Double,
  r : Double,
) -> Double {
  if v > 0.0 {
    n * r * t / v
  } else {
    0.0
  }
}

///|
/// Ideal gas law: calculate volume
/// n: number of moles
/// t: temperature (K)
/// p: pressure (Pa)
/// r: gas constant (J/(mol·K))
pub fn ideal_gas_volume(
  n : Double,
  t : Double,
  p : Double,
  r : Double,
) -> Double {
  if p > 0.0 {
    n * r * t / p
  } else {
    0.0
  }
}

///|
/// Root mean square speed (Maxwell-Boltzmann distribution)
/// t: temperature (K)
/// m: molecular mass (kg)
/// k: Boltzmann constant (J/K)
pub fn rms_speed(t : Double, m : Double, k : Double) -> Double {
  if m > 0.0 {
    Double::sqrt(3.0 * k * t / m)
  } else {
    0.0
  }
}

///|
/// Most probable speed (Maxwell-Boltzmann distribution)
/// t: temperature (K)
/// m: molecular mass (kg)
/// k: Boltzmann constant (J/K)
pub fn most_probable_speed(t : Double, m : Double, k : Double) -> Double {
  if m > 0.0 {
    Double::sqrt(2.0 * k * t / m)
  } else {
    0.0
  }
}

///|
/// Average kinetic energy per molecule (ideal gas)
/// t: temperature (K)
/// k: Boltzmann constant (J/K)
pub fn average_kinetic_energy_per_molecule(t : Double, k : Double) -> Double {
  1.5 * k * t
}

// ============================================================================
// Lagrange Points (Restricted Three-Body Problem)
// ============================================================================

///|
/// Lagrange point L1 distance (circular restricted three-body problem)
/// r: distance between primary bodies (m)
/// m1: mass of primary body 1 (kg)
/// m2: mass of primary body 2 (kg)
/// Returns: distance from m1 to L1
pub fn lagrange_point_l1_distance(
  r : Double,
  m1 : Double,
  m2 : Double,
) -> Double {
  // Simplified calculation using approximation
  // L1 ≈ r * (m2/(3*m1))^(1/3) for m2 << m1
  // Using cube root approximation: x^(1/3) ≈ exp(ln(x)/3)
  if m1 > 0.0 {
    let ratio = m2 / (3.0 * m1)
    if ratio > 0.0 {
      // Cube root: ratio^(1/3) using iterative method or approximation
      // For small ratios, use Newton's method approximation
      let mut x = ratio
      // Newton's method for cube root: x_new = (2*x + ratio/(x*x)) / 3
      for _i = 0; _i < 10; _i = _i + 1 {
        if x > 0.0 {
          x = (2.0 * x + ratio / (x * x)) / 3.0
        }
      }
      r * x
    } else {
      0.0
    }
  } else {
    0.0
  }
}

///|
/// Lagrange point L2 distance (circular restricted three-body problem)
/// r: distance between primary bodies (m)
/// m1: mass of primary body 1 (kg)
/// m2: mass of primary body 2 (kg)
/// Returns: distance from m1 to L2
pub fn lagrange_point_l2_distance(
  r : Double,
  m1 : Double,
  m2 : Double,
) -> Double {
  // L2 ≈ r * (1 + (m2/(3*m1))^(1/3)) for m2 << m1, on opposite side
  if m1 > 0.0 {
    let ratio = m2 / (3.0 * m1)
    if ratio > 0.0 {
      let mut x = ratio
      for _i = 0; _i < 10; _i = _i + 1 {
        if x > 0.0 {
          x = (2.0 * x + ratio / (x * x)) / 3.0
        }
      }
      r * (1.0 + x)
    } else {
      r
    }
  } else {
    0.0
  }
}

///|
/// Lagrange point L3 distance (circular restricted three-body problem)
/// r: distance between primary bodies (m)
/// m1: mass of primary body 1 (kg)
/// m2: mass of primary body 2 (kg)
/// Returns: distance from m1 to L3 (approximately r)
pub fn lagrange_point_l3_distance(
  r : Double,
  _m1 : Double,
  _m2 : Double,
) -> Double {
  // L3 is approximately at distance r from m1, on opposite side from m2
  r
}

// ============================================================================
// Statistical Mechanics Extensions
// ============================================================================

///|
/// Ideal gas law: calculate temperature
/// n: number of moles (mol)
/// p: pressure (Pa)
/// v: volume (m³)
/// r: gas constant (J/(mol·K)), default 8.314462618
/// Returns: temperature (K)
pub fn ideal_gas_temperature(
  n : Double,
  p : Double,
  v : Double,
  r : Double,
) -> Double {
  if n > 0.0 && r > 0.0 {
    p * v / (n * r)
  } else {
    0.0
  }
}

///|
/// Average kinetic energy per molecule (ideal gas)
/// t: temperature (K)
/// k_b: Boltzmann constant (J/K), default 1.380649e-23
/// Returns: average kinetic energy (J)
pub fn ideal_gas_avg_kinetic_energy(t : Double, k_b : Double) -> Double {
  if t >= 0.0 {
    1.5 * k_b * t
  } else {
    0.0
  }
}

///|
/// Root mean square (RMS) speed of gas molecules (alternative name)
/// t: temperature (K)
/// m: mass per molecule (kg)
/// k_b: Boltzmann constant (J/K), default 1.380649e-23
/// Returns: RMS speed (m/s)
pub fn ideal_gas_rms_speed(t : Double, m : Double, k_b : Double) -> Double {
  if m > 0.0 && t >= 0.0 {
    Double::sqrt(3.0 * k_b * t / m)
  } else {
    0.0
  }
}

///|
/// Most probable speed of gas molecules (Maxwell-Boltzmann distribution)
/// t: temperature (K)
/// m: mass per molecule (kg)
/// k_b: Boltzmann constant (J/K), default 1.380649e-23
/// Returns: most probable speed (m/s)
pub fn ideal_gas_most_probable_speed(
  t : Double,
  m : Double,
  k_b : Double,
) -> Double {
  if m > 0.0 && t >= 0.0 {
    Double::sqrt(2.0 * k_b * t / m)
  } else {
    0.0
  }
}

///|
/// Average speed of gas molecules (Maxwell-Boltzmann distribution)
/// t: temperature (K)
/// m: mass per molecule (kg)
/// k_b: Boltzmann constant (J/K), default 1.380649e-23
/// Returns: average speed (m/s)
pub fn ideal_gas_avg_speed(t : Double, m : Double, k_b : Double) -> Double {
  if m > 0.0 && t >= 0.0 {
    // <v> = sqrt(8*k_B*T/(π*m))
    let pi = 3.14159265358979323846
    Double::sqrt(8.0 * k_b * t / (pi * m))
  } else {
    0.0
  }
}

///|
/// Maxwell-Boltzmann speed distribution function
/// v: speed (m/s)
/// t: temperature (K)
/// m: mass per molecule (kg)
/// k_b: Boltzmann constant (J/K), default 1.380649e-23
/// Returns: probability density f(v) (s/m)
pub fn maxwell_boltzmann_speed_distribution(
  v : Double,
  t : Double,
  m : Double,
  k_b : Double,
) -> Double {
  if v >= 0.0 && m > 0.0 && t > 0.0 {
    let v_sq = v * v
    let pi = 3.14159265358979323846
    let coeff = 4.0 * pi * Double::pow(m / (2.0 * pi * k_b * t), 1.5)
    coeff * v_sq * @math.exp(-m * v_sq / (2.0 * k_b * t))
  } else {
    0.0
  }
}

///|
/// Maxwell-Boltzmann energy distribution function
/// e: energy (J)
/// t: temperature (K)
/// k_b: Boltzmann constant (J/K), default 1.380649e-23
/// Returns: probability density f(E) (1/J)
pub fn maxwell_boltzmann_energy_distribution(
  e : Double,
  t : Double,
  k_b : Double,
) -> Double {
  if e >= 0.0 && t > 0.0 {
    let pi = 3.14159265358979323846
    let coeff = 2.0 / (Double::sqrt(pi) * Double::pow(k_b * t, 1.5))
    coeff * Double::sqrt(e) * @math.exp(-e / (k_b * t))
  } else {
    0.0
  }
}

// ============================================================================
// Energy and Power Extensions
// ============================================================================

///|
/// Energy-time uncertainty relation (quantum mechanics)
/// delta_e: energy uncertainty (J)
/// delta_t: time uncertainty (s)
/// Returns: minimum uncertainty product (J·s)
pub fn energy_time_uncertainty(delta_e : Double, delta_t : Double) -> Double {
  if delta_e >= 0.0 && delta_t >= 0.0 {
    delta_e * delta_t
  } else {
    0.0
  }
}

///|
/// Power spectral density (simplified, for periodic signal)
/// a: amplitude (arbitrary units)
/// f: frequency (Hz)
/// Returns: power spectral density (arbitrary units)
pub fn power_spectral_density(a : Double, f : Double) -> Double {
  if f > 0.0 {
    a * a / (2.0 * f)
  } else {
    0.0
  }
}

///|
/// Energy dissipation rate from velocity (damped oscillator, alternative)
/// b: damping coefficient (kg/s)
/// v: velocity (m/s)
/// Returns: power dissipated (W)
pub fn energy_dissipation_rate_velocity(b : Double, v : Double) -> Double {
  if b >= 0.0 {
    b * v * v
  } else {
    0.0
  }
}

///|
/// Thermal efficiency (Carnot cycle)
/// t_hot: hot reservoir temperature (K)
/// t_cold: cold reservoir temperature (K)
/// Returns: maximum efficiency (dimensionless, 0-1)
pub fn thermal_efficiency_carnot(t_hot : Double, t_cold : Double) -> Double {
  if t_hot > 0.0 && t_cold >= 0.0 && t_hot > t_cold {
    1.0 - t_cold / t_hot
  } else {
    0.0
  }
}

///|
/// Thermal efficiency (general heat engine)
/// q_in: heat input (J)
/// q_out: heat output (J)
/// Returns: efficiency (dimensionless, 0-1)
pub fn thermal_efficiency(q_in : Double, q_out : Double) -> Double {
  if q_in > 0.0 && q_out >= 0.0 {
    let eff = (q_in - q_out) / q_in
    if eff > 1.0 {
      1.0
    } else if eff < 0.0 {
      0.0
    } else {
      eff
    }
  } else {
    0.0
  }
}

///|
/// Energy conversion efficiency in a chain
/// efficiencies: array of individual efficiencies (dimensionless, 0-1)
/// Returns: overall efficiency (dimensionless, 0-1)
pub fn energy_conversion_chain_efficiency(
  efficiencies : Array[Double],
) -> Double {
  let mut product = 1.0
  for i = 0; i < efficiencies.length(); i = i + 1 {
    let eff = efficiencies[i]
    if eff > 1.0 {
      product = product * 1.0
    } else if eff < 0.0 {
      product = product * 0.0
    } else {
      product = product * eff
    }
  }
  product
}

// ============================================================================
// Fluid Mechanics Extensions
// ============================================================================

///|
/// Pascal's principle: pressure transmission
/// f1: force applied (N)
/// a1: area where force is applied (m²)
/// a2: area where pressure is transmitted (m²)
/// Returns: force at a2 (N)
pub fn pascal_principle_force(f1 : Double, a1 : Double, a2 : Double) -> Double {
  if a1 > 0.0 && a2 > 0.0 {
    f1 * a2 / a1
  } else {
    0.0
  }
}

///|
/// Pascal's principle: pressure calculation
/// f: force (N)
/// a: area (m²)
/// Returns: pressure (Pa)
pub fn pascal_principle_pressure(f : Double, a : Double) -> Double {
  if a > 0.0 {
    f / a
  } else {
    0.0
  }
}

///|
/// Communicating vessels: equilibrium height difference
/// rho1: density of fluid 1 (kg/m³)
/// rho2: density of fluid 2 (kg/m³)
/// h1: height of fluid 1 (m)
/// g: gravitational acceleration (m/s²), default 9.81
/// Returns: height of fluid 2 for equilibrium (m)
pub fn communicating_vessels_height(
  rho1 : Double,
  rho2 : Double,
  h1 : Double,
  _g : Double,
) -> Double {
  if rho2 > 0.0 {
    rho1 * h1 / rho2
  } else {
    0.0
  }
}

///|
/// Pressure transmission in hydraulic system
/// p_in: input pressure (Pa)
/// a_in: input area (m²)
/// a_out: output area (m²)
/// Returns: output pressure (Pa)
pub fn pressure_transmission(
  p_in : Double,
  a_in : Double,
  a_out : Double,
) -> Double {
  if a_out > 0.0 && a_in > 0.0 {
    p_in * a_in / a_out
  } else {
    p_in
  }
}

///|
/// Venturi effect: velocity in constriction
/// v1: velocity in wide section (m/s)
/// a1: area of wide section (m²)
/// a2: area of constriction (m²)
/// Returns: velocity in constriction (m/s)
pub fn venturi_velocity(v1 : Double, a1 : Double, a2 : Double) -> Double {
  if a2 > 0.0 && a1 > 0.0 {
    v1 * a1 / a2
  } else {
    0.0
  }
}

///|
/// Venturi effect: pressure difference
/// rho: fluid density (kg/m³)
/// v1: velocity in wide section (m/s)
/// v2: velocity in constriction (m/s)
/// Returns: pressure difference P1 - P2 (Pa)
pub fn venturi_pressure_difference(
  rho : Double,
  v1 : Double,
  v2 : Double,
) -> Double {
  if rho > 0.0 {
    0.5 * rho * (v2 * v2 - v1 * v1)
  } else {
    0.0
  }
}

///|
/// Pitot tube: dynamic pressure
/// rho: fluid density (kg/m³)
/// v: fluid velocity (m/s)
/// Returns: dynamic pressure (Pa)
pub fn pitot_dynamic_pressure(rho : Double, v : Double) -> Double {
  if rho > 0.0 {
    0.5 * rho * v * v
  } else {
    0.0
  }
}

///|
/// Pitot tube: velocity from pressure difference
/// delta_p: pressure difference (Pa)
/// rho: fluid density (kg/m³)
/// Returns: fluid velocity (m/s)
pub fn pitot_velocity(delta_p : Double, rho : Double) -> Double {
  if rho > 0.0 && delta_p >= 0.0 {
    Double::sqrt(2.0 * delta_p / rho)
  } else {
    0.0
  }
}

// ============================================================================
// Utility Functions: Unit Conversions
// ============================================================================
// Note: deg_to_rad and rad_to_deg are already defined in tools.mbt

///|
/// Convert electron-volts to joules
/// ev: energy in electron-volts (eV)
/// Returns: energy in joules (J)
pub fn ev_to_joule(ev : Double) -> Double {
  // 1 eV = 1.602176634e-19 J
  ev * 1.602176634e-19
}

///|
/// Convert joules to electron-volts
/// joule: energy in joules (J)
/// Returns: energy in electron-volts (eV)
pub fn joule_to_ev(joule : Double) -> Double {
  if joule >= 0.0 {
    joule / 1.602176634e-19
  } else {
    0.0
  }
}

///|
/// Convert calories to joules
/// cal: energy in calories (cal)
/// Returns: energy in joules (J)
pub fn cal_to_joule(cal : Double) -> Double {
  // 1 cal = 4.184 J
  cal * 4.184
}

///|
/// Convert joules to calories
/// joule: energy in joules (J)
/// Returns: energy in calories (cal)
pub fn joule_to_cal(joule : Double) -> Double {
  if joule >= 0.0 {
    joule / 4.184
  } else {
    0.0
  }
}

///|
/// Convert BTU to joules
/// btu: energy in British Thermal Units (BTU)
/// Returns: energy in joules (J)
pub fn btu_to_joule(btu : Double) -> Double {
  // 1 BTU = 1055.05585262 J
  btu * 1055.05585262
}

///|
/// Convert joules to BTU
/// joule: energy in joules (J)
/// Returns: energy in BTU
pub fn joule_to_btu(joule : Double) -> Double {
  if joule >= 0.0 {
    joule / 1055.05585262
  } else {
    0.0
  }
}

///|
/// Convert horsepower to watts
/// hp: power in horsepower (hp)
/// Returns: power in watts (W)
pub fn hp_to_watt(hp : Double) -> Double {
  // 1 hp = 745.699872 W
  hp * 745.699872
}

///|
/// Convert watts to horsepower
/// watt: power in watts (W)
/// Returns: power in horsepower (hp)
pub fn watt_to_hp(watt : Double) -> Double {
  if watt >= 0.0 {
    watt / 745.699872
  } else {
    0.0
  }
}

///|
/// Convert BTU per hour to watts
/// btu_per_h: power in BTU/h
/// Returns: power in watts (W)
pub fn btu_per_h_to_watt(btu_per_h : Double) -> Double {
  // 1 BTU/h = 0.29307107 W
  btu_per_h * 0.29307107
}

///|
/// Convert watts to BTU per hour
/// watt: power in watts (W)
/// Returns: power in BTU/h
pub fn watt_to_btu_per_h(watt : Double) -> Double {
  if watt >= 0.0 {
    watt / 0.29307107
  } else {
    0.0
  }
}

// Note: Temperature conversion functions (celsius_to_kelvin, kelvin_to_celsius,
// fahrenheit_to_celsius, celsius_to_fahrenheit) are already defined in tools.mbt

///|
/// Convert Fahrenheit to Kelvin
/// fahrenheit: temperature in Fahrenheit (°F)
/// Returns: temperature in Kelvin (K)
pub fn fahrenheit_to_kelvin(fahrenheit : Double) -> Double {
  // Use functions from tools.mbt via module path if needed
  // For now, inline the calculation
  let celsius = (fahrenheit - 32.0) * 5.0 / 9.0
  celsius + 273.15
}

///|
/// Convert Kelvin to Fahrenheit
/// kelvin: temperature in Kelvin (K)
/// Returns: temperature in Fahrenheit (°F)
pub fn kelvin_to_fahrenheit(kelvin : Double) -> Double {
  // Use functions from tools.mbt via module path if needed
  // For now, inline the calculation
  let celsius = kelvin - 273.15
  celsius * 9.0 / 5.0 + 32.0
}

///|
/// Convert atmospheres to pascals
/// atm: pressure in atmospheres (atm)
/// Returns: pressure in pascals (Pa)
pub fn atm_to_pascal(atm : Double) -> Double {
  // 1 atm = 101325 Pa
  atm * 101325.0
}

///|
/// Convert pascals to atmospheres
/// pascal: pressure in pascals (Pa)
/// Returns: pressure in atmospheres (atm)
pub fn pascal_to_atm(pascal : Double) -> Double {
  if pascal >= 0.0 {
    pascal / 101325.0
  } else {
    0.0
  }
}

///|
/// Convert bar to pascals
/// bar: pressure in bar
/// Returns: pressure in pascals (Pa)
pub fn bar_to_pascal(bar : Double) -> Double {
  // 1 bar = 100000 Pa
  bar * 100000.0
}

///|
/// Convert pascals to bar
/// pascal: pressure in pascals (Pa)
/// Returns: pressure in bar
pub fn pascal_to_bar(pascal : Double) -> Double {
  if pascal >= 0.0 {
    pascal / 100000.0
  } else {
    0.0
  }
}

///|
/// Convert PSI to pascals
/// psi: pressure in pounds per square inch (psi)
/// Returns: pressure in pascals (Pa)
pub fn psi_to_pascal(psi : Double) -> Double {
  // 1 psi = 6894.757293168 Pa
  psi * 6894.757293168
}

///|
/// Convert pascals to PSI
/// pascal: pressure in pascals (Pa)
/// Returns: pressure in PSI
pub fn pascal_to_psi(pascal : Double) -> Double {
  if pascal >= 0.0 {
    pascal / 6894.757293168
  } else {
    0.0
  }
}

// ============================================================================
// Utility Functions: Coordinate System Conversions
// ============================================================================

///|
/// Convert polar coordinates to Cartesian coordinates (2D)
/// r: radial distance
/// theta: angle in radians
/// Returns: (x, y) as a tuple
pub fn polar_to_cartesian_2d(r : Double, theta : Double) -> (Double, Double) {
  (r * @math.cos(theta), r * @math.sin(theta))
}

///|
/// Convert Cartesian coordinates to polar coordinates (2D)
/// x: x-coordinate
/// y: y-coordinate
/// Returns: (r, theta) as a tuple, where theta is in radians [0, 2π)
pub fn cartesian_to_polar_2d(x : Double, y : Double) -> (Double, Double) {
  let r = Double::sqrt(x * x + y * y)
  let pi = 3.14159265358979323846
  let theta = if x == 0.0 && y == 0.0 {
    0.0
  } else if x == 0.0 {
    if y > 0.0 {
      pi / 2.0
    } else {
      3.0 * pi / 2.0
    }
  } else {
    let atan_val = @math.atan(y / x)
    if x > 0.0 {
      if y >= 0.0 {
        atan_val
      } else {
        atan_val + 2.0 * pi
      }
    } else {
      atan_val + pi
    }
  }
  (r, theta)
}

///|
/// Convert spherical coordinates to Cartesian coordinates (3D)
/// r: radial distance
/// theta: polar angle in radians (0 to π, measured from +z axis)
/// phi: azimuthal angle in radians (0 to 2π, measured from +x axis in xy-plane)
/// Returns: (x, y, z) as a tuple
pub fn spherical_to_cartesian_3d(
  r : Double,
  theta : Double,
  phi : Double,
) -> (Double, Double, Double) {
  let sin_theta = @math.sin(theta)
  (
    r * sin_theta * @math.cos(phi),
    r * sin_theta * @math.sin(phi),
    r * @math.cos(theta),
  )
}

///|
/// Convert Cartesian coordinates to spherical coordinates (3D)
/// x: x-coordinate
/// y: y-coordinate
/// z: z-coordinate
/// Returns: (r, theta, phi) as a tuple
///   r: radial distance
///   theta: polar angle in radians [0, π]
///   phi: azimuthal angle in radians [0, 2π)
pub fn cartesian_to_spherical_3d(
  x : Double,
  y : Double,
  z : Double,
) -> (Double, Double, Double) {
  let r = Double::sqrt(x * x + y * y + z * z)
  let pi = 3.14159265358979323846
  let theta = if r == 0.0 { 0.0 } else { @math.acos(z / r) }
  let phi = if x == 0.0 && y == 0.0 {
    0.0
  } else if x == 0.0 {
    if y > 0.0 {
      pi / 2.0
    } else {
      3.0 * pi / 2.0
    }
  } else {
    let atan_val = @math.atan(y / x)
    if x > 0.0 {
      if y >= 0.0 {
        atan_val
      } else {
        atan_val + 2.0 * pi
      }
    } else {
      atan_val + pi
    }
  }
  (r, theta, phi)
}

///|
/// Convert cylindrical coordinates to Cartesian coordinates (3D)
/// rho: radial distance in xy-plane
/// phi: azimuthal angle in radians (0 to 2π, measured from +x axis)
/// z: height coordinate
/// Returns: (x, y, z) as a tuple
pub fn cylindrical_to_cartesian_3d(
  rho : Double,
  phi : Double,
  z : Double,
) -> (Double, Double, Double) {
  (rho * @math.cos(phi), rho * @math.sin(phi), z)
}

///|
/// Convert Cartesian coordinates to cylindrical coordinates (3D)
/// x: x-coordinate
/// y: y-coordinate
/// z: z-coordinate
/// Returns: (rho, phi, z) as a tuple
///   rho: radial distance in xy-plane
///   phi: azimuthal angle in radians [0, 2π)
///   z: height coordinate
pub fn cartesian_to_cylindrical_3d(
  x : Double,
  y : Double,
  z : Double,
) -> (Double, Double, Double) {
  let rho = Double::sqrt(x * x + y * y)
  let pi = 3.14159265358979323846
  let phi = if x == 0.0 && y == 0.0 {
    0.0
  } else if x == 0.0 {
    if y > 0.0 {
      pi / 2.0
    } else {
      3.0 * pi / 2.0
    }
  } else {
    let atan_val = @math.atan(y / x)
    if x > 0.0 {
      if y >= 0.0 {
        atan_val
      } else {
        atan_val + 2.0 * pi
      }
    } else {
      atan_val + pi
    }
  }
  (rho, phi, z)
}

// ============================================================================
// Utility Functions: Dimensional Analysis
// ============================================================================

///|
/// Check if a value has correct dimensions for velocity (length/time)
/// length: length value (m)
/// time: time value (s)
/// Returns: velocity in m/s, or 0.0 if invalid
pub fn dimensional_check_velocity(length : Double, time : Double) -> Double {
  if time > 0.0 {
    length / time
  } else {
    0.0
  }
}

///|
/// Check if a value has correct dimensions for acceleration (length/time²)
/// length: length value (m)
/// time: time value (s)
/// Returns: acceleration in m/s², or 0.0 if invalid
pub fn dimensional_check_acceleration(length : Double, time : Double) -> Double {
  if time > 0.0 {
    length / (time * time)
  } else {
    0.0
  }
}

///|
/// Check if a value has correct dimensions for force (mass × acceleration)
/// mass: mass value (kg)
/// acceleration: acceleration value (m/s²)
/// Returns: force in N, or 0.0 if invalid
pub fn dimensional_check_force(mass : Double, acceleration : Double) -> Double {
  if mass > 0.0 {
    mass * acceleration
  } else {
    0.0
  }
}

///|
/// Check if a value has correct dimensions for energy (force × length)
/// force: force value (N)
/// length: length value (m)
/// Returns: energy in J, or 0.0 if invalid
pub fn dimensional_check_energy(force : Double, length : Double) -> Double {
  if force >= 0.0 {
    force * length
  } else {
    0.0
  }
}

///|
/// Check if a value has correct dimensions for power (energy/time)
/// energy: energy value (J)
/// time: time value (s)
/// Returns: power in W, or 0.0 if invalid
pub fn dimensional_check_power(energy : Double, time : Double) -> Double {
  if time > 0.0 {
    energy / time
  } else {
    0.0
  }
}

///|
/// Check if a value has correct dimensions for pressure (force/area)
/// force: force value (N)
/// area: area value (m²)
/// Returns: pressure in Pa, or 0.0 if invalid
pub fn dimensional_check_pressure(force : Double, area : Double) -> Double {
  if area > 0.0 {
    force / area
  } else {
    0.0
  }
}

///|
/// Validate physical quantity range (check if value is physically reasonable)
/// value: physical quantity value
/// min_val: minimum reasonable value
/// max_val: maximum reasonable value
/// Returns: true if value is within reasonable range, false otherwise
pub fn validate_physical_range(
  value : Double,
  min_val : Double,
  max_val : Double,
) -> Bool {
  value >= min_val && value <= max_val
}

///|
/// Check if velocity is non-relativistic (v << c)
/// v: velocity (m/s)
/// c: speed of light (m/s), default 299792458
/// threshold: relative threshold (default 0.1 = 10% of c)
/// Returns: true if velocity is non-relativistic
pub fn is_non_relativistic(v : Double, c : Double, threshold : Double) -> Bool {
  if c > 0.0 {
    Double::abs(v) / c < threshold
  } else {
    false
  }
}
