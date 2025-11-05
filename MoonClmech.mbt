/// MoonClmech - MoonBit 经典力学库
/// 基于经典力学原理实现的主要计算功能
/// 提供粒子运动、刚体力学、拉格朗日、哈密顿等计算函数

// ============================================================================
// 1. 粒子运动学（Kinematics）
// ============================================================================

///|
/// 匀速直线运动：计算位移
/// v: 速度 (m/s)
/// t: 时间 (s)
pub fn uniform_motion_displacement(velocity : Double, time : Double) -> Double {
  velocity * time
}

///|
/// 匀加速直线运动：计算位移
/// v0: 初速度 (m/s)
/// a: 加速度 (m/s²)
/// t: 时间 (s)
pub fn uniform_acceleration_displacement(
  v0 : Double,
  a : Double,
  t : Double,
) -> Double {
  v0 * t + 0.5 * a * t * t
}

///|
/// 匀加速直线运动：计算速度
/// v0: 初速度 (m/s)
/// a: 加速度 (m/s²)
/// t: 时间 (s)
pub fn uniform_acceleration_velocity(
  v0 : Double,
  a : Double,
  t : Double,
) -> Double {
  v0 + a * t
}

///|
/// 匀加速直线运动：计算速度（使用位移）
/// v0: 初速度 (m/s)
/// a: 加速度 (m/s²)
/// s: 位移 (m)
pub fn uniform_acceleration_velocity_from_displacement(
  v0 : Double,
  a : Double,
  s : Double,
) -> Double {
  Double::sqrt(v0 * v0 + 2.0 * a * s)
}

///|
/// 抛体运动：水平位移
/// v0: 初速度 (m/s)
/// angle: 发射角度（弧度）
/// t: 时间 (s)
pub fn projectile_horizontal_displacement(
  v0 : Double,
  angle : Double,
  t : Double,
) -> Double {
  v0 * @math.cos(angle) * t
}

///|
/// 抛体运动：垂直位移
/// v0: 初速度 (m/s)
/// angle: 发射角度（弧度）
/// t: 时间 (s)
/// g: 重力加速度 (m/s²)，默认使用地球重力
pub fn projectile_vertical_displacement(
  v0 : Double,
  angle : Double,
  t : Double,
  g : Double,
) -> Double {
  v0 * @math.sin(angle) * t - 0.5 * g * t * t
}

///|
/// 抛体运动：最大高度
/// v0: 初速度 (m/s)
/// angle: 发射角度（弧度）
/// g: 重力加速度 (m/s²)
pub fn projectile_max_height(v0 : Double, angle : Double, g : Double) -> Double {
  let v0y = v0 * @math.sin(angle)
  v0y * v0y / (2.0 * g)
}

///|
/// 抛体运动：射程
/// v0: 初速度 (m/s)
/// angle: 发射角度（弧度）
/// g: 重力加速度 (m/s²)
pub fn projectile_range(v0 : Double, angle : Double, g : Double) -> Double {
  v0 * v0 * @math.sin(2.0 * angle) / g
}

///|
/// 抛体运动：飞行时间
/// v0: 初速度 (m/s)
/// angle: 发射角度（弧度）
/// g: 重力加速度 (m/s²)
pub fn projectile_flight_time(
  v0 : Double,
  angle : Double,
  g : Double,
) -> Double {
  2.0 * v0 * @math.sin(angle) / g
}

// ============================================================================
// 2. 动力学（Dynamics）
// ============================================================================

///|
/// 牛顿第二定律：计算力
/// m: 质量 (kg)
/// a: 加速度 (m/s²)
pub fn force(mass : Double, acceleration : Double) -> Double {
  mass * acceleration
}

///|
/// 牛顿第二定律：计算加速度
/// f: 力 (N)
/// m: 质量 (kg)
pub fn acceleration(force : Double, mass : Double) -> Double {
  force / mass
}

///|
/// 动量
/// m: 质量 (kg)
/// v: 速度 (m/s)
pub fn momentum(mass : Double, velocity : Double) -> Double {
  mass * velocity
}

///|
/// 动量变化
/// m: 质量 (kg)
/// v1: 初速度 (m/s)
/// v2: 末速度 (m/s)
pub fn momentum_change(mass : Double, v1 : Double, v2 : Double) -> Double {
  mass * (v2 - v1)
}

///|
/// 冲量
/// f: 力 (N)
/// t: 时间 (s)
pub fn impulse(force : Double, time : Double) -> Double {
  force * time
}

///|
/// 动能
/// m: 质量 (kg)
/// v: 速度 (m/s)
pub fn kinetic_energy(mass : Double, velocity : Double) -> Double {
  0.5 * mass * velocity * velocity
}

///|
/// 重力势能
/// m: 质量 (kg)
/// h: 高度 (m)
/// g: 重力加速度 (m/s²)
pub fn gravitational_potential_energy(
  mass : Double,
  height : Double,
  g : Double,
) -> Double {
  mass * g * height
}

///|
/// 弹性势能
/// k: 弹簧常数 (N/m)
/// x: 位移 (m)
pub fn elastic_potential_energy(
  spring_constant : Double,
  displacement : Double,
) -> Double {
  0.5 * spring_constant * displacement * displacement
}

///|
/// 机械能（动能 + 势能）
pub fn mechanical_energy(kinetic : Double, potential : Double) -> Double {
  kinetic + potential
}

///|
/// 万有引力
/// m1: 物体1质量 (kg)
/// m2: 物体2质量 (kg)
/// r: 距离 (m)
/// g_const: 万有引力常数 (m³/(kg·s²))
pub fn gravitational_force(
  m1 : Double,
  m2 : Double,
  r : Double,
  g_const : Double,
) -> Double {
  g_const * m1 * m2 / (r * r)
}

// ============================================================================
// 3. 刚体力学（Rigid Body Mechanics）
// ============================================================================

///|
/// 转动惯量：质点
/// m: 质量 (kg)
/// r: 距离转轴的距离 (m)
pub fn moment_of_inertia_point_mass(mass : Double, radius : Double) -> Double {
  mass * radius * radius
}

///|
/// 转动惯量：细杆（绕中心）
/// m: 质量 (kg)
/// l: 长度 (m)
pub fn moment_of_inertia_rod_center(mass : Double, length : Double) -> Double {
  mass * length * length / 12.0
}

///|
/// 转动惯量：细杆（绕一端）
/// m: 质量 (kg)
/// l: 长度 (m)
pub fn moment_of_inertia_rod_end(mass : Double, length : Double) -> Double {
  mass * length * length / 3.0
}

///|
/// 转动惯量：圆环（绕中心轴）
/// m: 质量 (kg)
/// r: 半径 (m)
pub fn moment_of_inertia_ring(mass : Double, radius : Double) -> Double {
  mass * radius * radius
}

///|
/// 转动惯量：圆盘/圆柱（绕中心轴）
/// m: 质量 (kg)
/// r: 半径 (m)
pub fn moment_of_inertia_disk(mass : Double, radius : Double) -> Double {
  0.5 * mass * radius * radius
}

///|
/// 转动惯量：实心球（绕直径）
/// m: 质量 (kg)
/// r: 半径 (m)
pub fn moment_of_inertia_sphere(mass : Double, radius : Double) -> Double {
  0.4 * mass * radius * radius
}

///|
/// 转动惯量：空心球（绕直径）
/// m: 质量 (kg)
/// r: 半径 (m)
pub fn moment_of_inertia_hollow_sphere(
  mass : Double,
  radius : Double,
) -> Double {
  2.0 / 3.0 * mass * radius * radius
}

///|
/// 角动量
/// I: 转动惯量 (kg·m²)
/// omega: 角速度 (rad/s)
pub fn angular_momentum(
  moment_of_inertia : Double,
  angular_velocity : Double,
) -> Double {
  moment_of_inertia * angular_velocity
}

///|
/// 角动量（质点）
/// m: 质量 (kg)
/// r: 距离转轴的距离 (m)
/// v: 线速度 (m/s)
pub fn angular_momentum_point_mass(
  mass : Double,
  radius : Double,
  velocity : Double,
) -> Double {
  mass * radius * velocity
}

///|
/// 转动动能
/// I: 转动惯量 (kg·m²)
/// omega: 角速度 (rad/s)
pub fn rotational_kinetic_energy(
  moment_of_inertia : Double,
  angular_velocity : Double,
) -> Double {
  0.5 * moment_of_inertia * angular_velocity * angular_velocity
}

///|
/// 力矩
/// f: 力 (N)
/// r: 力臂 (m)
/// angle: 力和力臂的夹角（弧度）
pub fn torque(force : Double, lever_arm : Double, angle : Double) -> Double {
  force * lever_arm * @math.sin(angle)
}

///|
/// 角加速度
/// tau: 力矩 (N·m)
/// I: 转动惯量 (kg·m²)
pub fn angular_acceleration(
  torque : Double,
  moment_of_inertia : Double,
) -> Double {
  torque / moment_of_inertia
}

///|
/// 角速度（匀角加速运动）
/// omega0: 初角速度 (rad/s)
/// alpha: 角加速度 (rad/s²)
/// t: 时间 (s)
pub fn angular_velocity(omega0 : Double, alpha : Double, t : Double) -> Double {
  omega0 + alpha * t
}

///|
/// 角位移（匀角加速运动）
/// omega0: 初角速度 (rad/s)
/// alpha: 角加速度 (rad/s²)
/// t: 时间 (s)
pub fn angular_displacement(
  omega0 : Double,
  alpha : Double,
  t : Double,
) -> Double {
  omega0 * t + 0.5 * alpha * t * t
}

///|
/// 平行轴定理：计算转动惯量
/// i_cm: 绕质心的转动惯量 (kg·m²)
/// m: 质量 (kg)
/// d: 到新轴的距离 (m)
pub fn parallel_axis_theorem(
  i_cm : Double,
  mass : Double,
  distance : Double,
) -> Double {
  i_cm + mass * distance * distance
}

// ============================================================================
// 4. 简谐振动（Simple Harmonic Motion）
// ============================================================================

///|
/// 简谐振动：位移
/// A: 振幅 (m)
/// omega: 角频率 (rad/s)
/// t: 时间 (s)
/// phi: 初相位（弧度）
pub fn shm_displacement(
  amplitude : Double,
  angular_frequency : Double,
  time : Double,
  phase : Double,
) -> Double {
  amplitude * @math.cos(angular_frequency * time + phase)
}

///|
/// 简谐振动：速度
pub fn shm_velocity(
  amplitude : Double,
  angular_frequency : Double,
  time : Double,
  phase : Double,
) -> Double {
  -amplitude * angular_frequency * @math.sin(angular_frequency * time + phase)
}

///|
/// 简谐振动：加速度
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
/// 简谐振动：角频率
/// k: 弹簧常数 (N/m)
/// m: 质量 (kg)
pub fn shm_angular_frequency(spring_constant : Double, mass : Double) -> Double {
  Double::sqrt(spring_constant / mass)
}

///|
/// 简谐振动：周期
/// omega: 角频率 (rad/s)
pub fn shm_period(angular_frequency : Double) -> Double {
  2.0 * 3.141592653589793 / angular_frequency
}

///|
/// 简谐振动：频率
/// omega: 角频率 (rad/s)
pub fn shm_frequency(angular_frequency : Double) -> Double {
  angular_frequency / (2.0 * 3.141592653589793)
}

///|
/// 单摆：角频率
/// g: 重力加速度 (m/s²)
/// l: 摆长 (m)
pub fn pendulum_angular_frequency(gravity : Double, length : Double) -> Double {
  Double::sqrt(gravity / length)
}

///|
/// 单摆：周期
pub fn pendulum_period(gravity : Double, length : Double) -> Double {
  2.0 * 3.141592653589793 * Double::sqrt(length / gravity)
}

// ============================================================================
// 5. 碰撞（Collisions）
// ============================================================================

///|
/// 弹性碰撞：一维完全弹性碰撞后的速度
/// m1: 物体1质量 (kg)
/// v1: 物体1初速度 (m/s)
/// m2: 物体2质量 (kg)
/// v2: 物体2初速度 (m/s)
/// 返回：物体1的末速度
pub fn elastic_collision_v1(
  m1 : Double,
  v1 : Double,
  m2 : Double,
  v2 : Double,
) -> Double {
  (m1 - m2) * v1 / (m1 + m2) + 2.0 * m2 * v2 / (m1 + m2)
}

///|
/// 弹性碰撞：一维完全弹性碰撞后的速度
/// 返回：物体2的末速度
pub fn elastic_collision_v2(
  m1 : Double,
  v1 : Double,
  m2 : Double,
  v2 : Double,
) -> Double {
  2.0 * m1 * v1 / (m1 + m2) + (m2 - m1) * v2 / (m1 + m2)
}

///|
/// 完全非弹性碰撞：碰撞后速度
pub fn inelastic_collision_velocity(
  m1 : Double,
  v1 : Double,
  m2 : Double,
  v2 : Double,
) -> Double {
  (m1 * v1 + m2 * v2) / (m1 + m2)
}

///|
/// 恢复系数：碰撞后相对速度与碰撞前相对速度的比值
/// v1f: 物体1碰撞后速度 (m/s)
/// v2f: 物体2碰撞后速度 (m/s)
/// v1i: 物体1碰撞前速度 (m/s)
/// v2i: 物体2碰撞前速度 (m/s)
pub fn coefficient_of_restitution(
  v1f : Double,
  v2f : Double,
  v1i : Double,
  v2i : Double,
) -> Double {
  (v2f - v1f) / (v1i - v2i)
}

// ============================================================================
// 6. 中心力场（Central Force Motion）
// ============================================================================

///|
/// 开普勒第三定律：轨道周期
/// a: 半长轴 (m)
/// M: 中心天体质量 (kg)
/// g_const: 万有引力常数 (m³/(kg·s²))
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
/// 轨道速度（圆形轨道）
/// r: 轨道半径 (m)
/// M: 中心天体质量 (kg)
/// g_const: 万有引力常数
pub fn orbital_velocity_circular(
  radius : Double,
  central_mass : Double,
  g_const : Double,
) -> Double {
  Double::sqrt(g_const * central_mass / radius)
}

///|
/// 逃逸速度
/// r: 距离中心天体的距离 (m)
/// M: 中心天体质量 (kg)
/// g_const: 万有引力常数
pub fn escape_velocity(
  radius : Double,
  central_mass : Double,
  g_const : Double,
) -> Double {
  Double::sqrt(2.0 * g_const * central_mass / radius)
}

///|
/// 第一宇宙速度（地球表面）
pub fn first_cosmic_velocity(gravity : Double, radius : Double) -> Double {
  Double::sqrt(gravity * radius)
}

///|
/// 第二宇宙速度（地球表面逃逸速度）
pub fn second_cosmic_velocity(gravity : Double, radius : Double) -> Double {
  Double::sqrt(2.0 * gravity * radius)
}

// ============================================================================
// 7. 拉格朗日力学（Lagrangian Mechanics）
// ============================================================================

///|
/// 拉格朗日函数 L = T - V
/// T: 动能
/// V: 势能
pub fn lagrangian(kinetic_energy : Double, potential_energy : Double) -> Double {
  kinetic_energy - potential_energy
}

///|
/// 拉格朗日方程：计算广义力（保守系统）
/// 对于 L = L(q, q_dot, t)，广义力 Q = d/dt(dL/dq_dot) - dL/dq
/// 这里提供简化版本，假设 L 不显式依赖于时间
/// 
/// 注意：这是简化版本，完整实现需要符号微分
/// 
/// dL_dq_dot: dL/dq_dot 对时间的导数
/// dL_dq: dL/dq
pub fn lagrangian_generalized_force(
  dL_dq_dot_dt : Double,
  dL_dq : Double,
) -> Double {
  dL_dq_dot_dt - dL_dq
}

// ============================================================================
// 8. 哈密顿力学（Hamiltonian Mechanics）
// ============================================================================

///|
/// 哈密顿函数 H = T + V（对于保守系统）
/// T: 动能
/// V: 势能
pub fn hamiltonian(
  kinetic_energy : Double,
  potential_energy : Double,
) -> Double {
  kinetic_energy + potential_energy
}

///|
/// 哈密顿函数 H = Σ(p_i * q_dot_i) - L
/// 对于单自由度系统：H = p * q_dot - L
/// p: 广义动量
/// q_dot: 广义速度
/// L: 拉格朗日函数
pub fn hamiltonian_from_lagrangian(
  momentum : Double,
  velocity : Double,
  lagrangian : Double,
) -> Double {
  momentum * velocity - lagrangian
}

// ============================================================================
// 9. 多体系统（Multi-body Systems）
// ============================================================================

///|
/// 质心位置（一维，多个质点）
/// masses: 质量数组 (kg)
/// positions: 位置数组 (m)
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
/// 系统总动量（一维）
/// masses: 质量数组 (kg)
/// velocities: 速度数组 (m/s)
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
/// 系统总动能（一维）
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

// ============================================================================
// 10. 阻尼振动（Damped Oscillations）
// ============================================================================

///|
/// 阻尼振动的阻尼系数
/// b: 阻尼系数 (kg/s)
/// m: 质量 (kg)
pub fn damping_coefficient(b : Double, mass : Double) -> Double {
  b / (2.0 * mass)
}

///|
/// 临界阻尼系数
/// m: 质量 (kg)
/// k: 弹簧常数 (N/m)
pub fn critical_damping_coefficient(
  mass : Double,
  spring_constant : Double,
) -> Double {
  2.0 * Double::sqrt(mass * spring_constant)
}

///|
/// 阻尼振动的角频率（欠阻尼）
/// omega0: 固有角频率 (rad/s)
/// gamma: 阻尼系数
pub fn damped_angular_frequency(
  natural_frequency : Double,
  damping : Double,
) -> Double {
  Double::sqrt(natural_frequency * natural_frequency - damping * damping)
}
