/// MoonClmech 测试用例
/// 全面测试经典力学库的各项功能

// ============================================================================
// 1. 粒子运动学测试（Kinematics）
// ============================================================================

///|
test "uniform_motion_displacement - basic calculation" {
  // 匀速直线运动：速度 10 m/s，时间 5 s，位移应为 50 m
  let result = uniform_motion_displacement(10.0, 5.0)
  assert_eq(result, 50.0)
}

///|
test "uniform_motion_displacement - zero velocity" {
  let result = uniform_motion_displacement(0.0, 5.0)
  assert_eq(result, 0.0)
}

///|
test "uniform_acceleration_displacement - basic calculation" {
  // 匀加速运动：初速度 0，加速度 2 m/s²，时间 5 s
  // 位移 = 0.5 * 2 * 5² = 25 m
  let result = uniform_acceleration_displacement(0.0, 2.0, 5.0)
  assert_eq(result, 25.0)
}

///|
test "uniform_acceleration_displacement - with initial velocity" {
  // 初速度 5 m/s，加速度 2 m/s²，时间 3 s
  // 位移 = 5*3 + 0.5*2*9 = 15 + 9 = 24 m
  let result = uniform_acceleration_displacement(5.0, 2.0, 3.0)
  assert_eq(result, 24.0)
}

///|
test "uniform_acceleration_velocity - basic calculation" {
  // 初速度 0，加速度 5 m/s²，时间 4 s
  // 速度 = 0 + 5*4 = 20 m/s
  let result = uniform_acceleration_velocity(0.0, 5.0, 4.0)
  assert_eq(result, 20.0)
}

///|
test "uniform_acceleration_velocity_from_displacement - basic calculation" {
  // 初速度 0，加速度 2 m/s²，位移 16 m
  // v² = 0 + 2*2*16 = 64, v = 8 m/s
  let result = uniform_acceleration_velocity_from_displacement(0.0, 2.0, 16.0)
  assert_eq((result - 8.0).abs() < 0.0001, true)
}

///|
test "projectile_horizontal_displacement - 45 degree angle" {
  // 初速度 10 m/s，角度 45°（π/4），时间 2 s
  // 水平位移 = 10 * cos(π/4) * 2 = 10 * 0.707 * 2 ≈ 14.14 m
  let angle = 3.141592653589793 / 4.0
  let result = projectile_horizontal_displacement(10.0, angle, 2.0)
  assert_eq((result - 14.142).abs() < 0.1, true)
}

///|
test "projectile_vertical_displacement - at peak" {
  // 初速度 20 m/s，角度 90°（垂直向上），时间 1 s，g = 9.8 m/s²
  // 垂直位移 = 20*sin(π/2)*1 - 0.5*9.8*1 = 20 - 4.9 = 15.1 m
  let angle = 3.141592653589793 / 2.0
  let g = 9.8
  let result = projectile_vertical_displacement(20.0, angle, 1.0, g)
  assert_eq((result - 15.1).abs() < 0.1, true)
}

///|
test "projectile_max_height - 45 degree angle" {
  // 初速度 20 m/s，角度 45°，g = 9.8 m/s²
  // v0y = 20*sin(45°) = 14.14 m/s
  // h_max = (14.14)² / (2*9.8) ≈ 10.2 m
  let angle = 3.141592653589793 / 4.0
  let g = 9.8
  let result = projectile_max_height(20.0, angle, g)
  assert_eq((result - 10.2).abs() < 0.5, true)
}

///|
test "projectile_range - 45 degree angle" {
  // 初速度 20 m/s，角度 45°，g = 9.8 m/s²
  // 射程 = 20² * sin(90°) / 9.8 = 400 / 9.8 ≈ 40.8 m
  let angle = 3.141592653589793 / 4.0
  let g = 9.8
  let result = projectile_range(20.0, angle, g)
  assert_eq((result - 40.8).abs() < 0.5, true)
}

///|
test "projectile_flight_time - vertical throw" {
  // 初速度 19.6 m/s，角度 90°，g = 9.8 m/s²
  // 飞行时间 = 2*19.6*sin(90°) / 9.8 = 39.2 / 9.8 = 4 s
  let angle = 3.141592653589793 / 2.0
  let g = 9.8
  let result = projectile_flight_time(19.6, angle, g)
  assert_eq((result - 4.0).abs() < 0.01, true)
}

// ============================================================================
// 2. 动力学测试（Dynamics）
// ============================================================================

///|
test "force - Newton's second law" {
  // 质量 5 kg，加速度 4 m/s²，力 = 20 N
  let result = force(5.0, 4.0)
  assert_eq(result, 20.0)
}

///|
test "acceleration - from force and mass" {
  // 力 30 N，质量 6 kg，加速度 = 5 m/s²
  let result = acceleration(30.0, 6.0)
  assert_eq(result, 5.0)
}

///|
test "momentum - basic calculation" {
  // 质量 2 kg，速度 10 m/s，动量 = 20 kg·m/s
  let result = momentum(2.0, 10.0)
  assert_eq(result, 20.0)
}

///|
test "momentum_change - velocity change" {
  // 质量 3 kg，初速度 5 m/s，末速度 15 m/s
  // 动量变化 = 3*(15-5) = 30 kg·m/s
  let result = momentum_change(3.0, 5.0, 15.0)
  assert_eq(result, 30.0)
}

///|
test "impulse - force and time" {
  // 力 10 N，时间 5 s，冲量 = 50 N·s
  let result = impulse(10.0, 5.0)
  assert_eq(result, 50.0)
}

///|
test "kinetic_energy - basic calculation" {
  // 质量 2 kg，速度 10 m/s，动能 = 0.5*2*100 = 100 J
  let result = kinetic_energy(2.0, 10.0)
  assert_eq(result, 100.0)
}

///|
test "gravitational_potential_energy - basic calculation" {
  // 质量 5 kg，高度 10 m，g = 9.8 m/s²
  // 势能 = 5*9.8*10 = 490 J
  let result = gravitational_potential_energy(5.0, 10.0, 9.8)
  assert_eq(result, 490.0)
}

///|
test "elastic_potential_energy - basic calculation" {
  // 弹簧常数 100 N/m，位移 0.5 m
  // 弹性势能 = 0.5*100*0.25 = 12.5 J
  let result = elastic_potential_energy(100.0, 0.5)
  assert_eq(result, 12.5)
}

///|
test "mechanical_energy - sum of kinetic and potential" {
  // 动能 100 J，势能 50 J，机械能 = 150 J
  let result = mechanical_energy(100.0, 50.0)
  assert_eq(result, 150.0)
}

///|
test "gravitational_force - Earth-Moon system" {
  // 地球质量 5.972e24 kg，月球质量 7.342e22 kg
  // 距离 3.844e8 m，G = 6.67430e-11
  // 万有引力 ≈ 1.98e20 N
  let m1 = 5.972e24
  let m2 = 7.342e22
  let r = 3.844e8
  let g_const = 6.67430e-11
  let result = gravitational_force(m1, m2, r, g_const)
  assert_eq((result - 1.98e20).abs() < 0.1e20, true)
}

// ============================================================================
// 3. 刚体力学测试（Rigid Body Mechanics）
// ============================================================================

///|
test "moment_of_inertia_point_mass - basic calculation" {
  // 质量 2 kg，距离 3 m，转动惯量 = 2*9 = 18 kg·m²
  let result = moment_of_inertia_point_mass(2.0, 3.0)
  assert_eq(result, 18.0)
}

///|
test "moment_of_inertia_rod_center - basic calculation" {
  // 质量 6 kg，长度 4 m，绕中心转动惯量 = 6*16/12 = 8 kg·m²
  let result = moment_of_inertia_rod_center(6.0, 4.0)
  assert_eq(result, 8.0)
}

///|
test "moment_of_inertia_rod_end - basic calculation" {
  // 质量 6 kg，长度 4 m，绕一端转动惯量 = 6*16/3 = 32 kg·m²
  let result = moment_of_inertia_rod_end(6.0, 4.0)
  assert_eq(result, 32.0)
}

///|
test "moment_of_inertia_ring - basic calculation" {
  // 质量 2 kg，半径 1 m，转动惯量 = 2*1 = 2 kg·m²
  let result = moment_of_inertia_ring(2.0, 1.0)
  assert_eq(result, 2.0)
}

///|
test "moment_of_inertia_disk - basic calculation" {
  // 质量 4 kg，半径 2 m，转动惯量 = 0.5*4*4 = 8 kg·m²
  let result = moment_of_inertia_disk(4.0, 2.0)
  assert_eq(result, 8.0)
}

///|
test "moment_of_inertia_sphere - basic calculation" {
  // 质量 5 kg，半径 2 m，转动惯量 = 0.4*5*4 = 8 kg·m²
  let result = moment_of_inertia_sphere(5.0, 2.0)
  assert_eq(result, 8.0)
}

///|
test "moment_of_inertia_hollow_sphere - basic calculation" {
  // 质量 3 kg，半径 2 m，转动惯量 = 2/3*3*4 = 8 kg·m²
  let result = moment_of_inertia_hollow_sphere(3.0, 2.0)
  assert_eq(result, 8.0)
}

///|
test "angular_momentum - basic calculation" {
  // 转动惯量 10 kg·m²，角速度 5 rad/s，角动量 = 50 kg·m²/s
  let result = angular_momentum(10.0, 5.0)
  assert_eq(result, 50.0)
}

///|
test "angular_momentum_point_mass - basic calculation" {
  // 质量 2 kg，距离 3 m，速度 4 m/s，角动量 = 2*3*4 = 24 kg·m²/s
  let result = angular_momentum_point_mass(2.0, 3.0, 4.0)
  assert_eq(result, 24.0)
}

///|
test "rotational_kinetic_energy - basic calculation" {
  // 转动惯量 8 kg·m²，角速度 3 rad/s
  // 转动动能 = 0.5*8*9 = 36 J
  let result = rotational_kinetic_energy(8.0, 3.0)
  assert_eq(result, 36.0)
}

///|
test "torque - perpendicular force" {
  // 力 10 N，力臂 2 m，角度 90°（π/2）
  // 力矩 = 10*2*sin(π/2) = 20 N·m
  let angle = 3.141592653589793 / 2.0
  let result = torque(10.0, 2.0, angle)
  assert_eq((result - 20.0).abs() < 0.0001, true)
}

///|
test "angular_acceleration - basic calculation" {
  // 力矩 30 N·m，转动惯量 6 kg·m²，角加速度 = 5 rad/s²
  let result = angular_acceleration(30.0, 6.0)
  assert_eq(result, 5.0)
}

///|
test "angular_velocity - uniform angular acceleration" {
  // 初角速度 2 rad/s，角加速度 3 rad/s²，时间 4 s
  // 角速度 = 2 + 3*4 = 14 rad/s
  let result = angular_velocity(2.0, 3.0, 4.0)
  assert_eq(result, 14.0)
}

///|
test "angular_displacement - uniform angular acceleration" {
  // 初角速度 1 rad/s，角加速度 2 rad/s²，时间 3 s
  // 角位移 = 1*3 + 0.5*2*9 = 3 + 9 = 12 rad
  let result = angular_displacement(1.0, 2.0, 3.0)
  assert_eq(result, 12.0)
}

///|
test "parallel_axis_theorem - basic calculation" {
  // 绕质心转动惯量 5 kg·m²，质量 2 kg，距离 3 m
  // 新轴转动惯量 = 5 + 2*9 = 23 kg·m²
  let result = parallel_axis_theorem(5.0, 2.0, 3.0)
  assert_eq(result, 23.0)
}

// ============================================================================
// 4. 简谐振动测试（Simple Harmonic Motion）
// ============================================================================

///|
test "shm_displacement - at t=0" {
  // 振幅 5 m，角频率 2 rad/s，时间 0 s，初相位 0
  // 位移 = 5*cos(0) = 5 m
  let result = shm_displacement(5.0, 2.0, 0.0, 0.0)
  assert_eq((result - 5.0).abs() < 0.0001, true)
}

///|
test "shm_velocity - at equilibrium" {
  // 振幅 5 m，角频率 2 rad/s，时间 0 s，初相位 0
  // 速度 = -5*2*sin(0) = 0 m/s
  let result = shm_velocity(5.0, 2.0, 0.0, 0.0)
  assert_eq((result - 0.0).abs() < 0.0001, true)
}

///|
test "shm_acceleration - at maximum displacement" {
  // 振幅 5 m，角频率 2 rad/s，时间 0 s，初相位 0
  // 加速度 = -5*4*cos(0) = -20 m/s²
  let result = shm_acceleration(5.0, 2.0, 0.0, 0.0)
  assert_eq((result - -20.0).abs() < 0.0001, true)
}

///|
test "shm_angular_frequency - spring-mass system" {
  // 弹簧常数 100 N/m，质量 4 kg
  // 角频率 = sqrt(100/4) = 5 rad/s
  let result = shm_angular_frequency(100.0, 4.0)
  assert_eq((result - 5.0).abs() < 0.0001, true)
}

///|
test "shm_period - from angular frequency" {
  // 角频率 2π rad/s，周期 = 2π/(2π) = 1 s
  let omega = 2.0 * 3.141592653589793
  let result = shm_period(omega)
  assert_eq((result - 1.0).abs() < 0.01, true)
}

///|
test "shm_frequency - from angular frequency" {
  // 角频率 2π rad/s，频率 = (2π)/(2π) = 1 Hz
  let omega = 2.0 * 3.141592653589793
  let result = shm_frequency(omega)
  assert_eq((result - 1.0).abs() < 0.01, true)
}

///|
test "pendulum_angular_frequency - basic calculation" {
  // 重力加速度 9.8 m/s²，摆长 1 m
  // 角频率 = sqrt(9.8/1) ≈ 3.13 rad/s
  let result = pendulum_angular_frequency(9.8, 1.0)
  assert_eq((result - 3.13).abs() < 0.1, true)
}

///|
test "pendulum_period - basic calculation" {
  // 重力加速度 9.8 m/s²，摆长 1 m
  // 周期 = 2π*sqrt(1/9.8) ≈ 2.0 s
  let result = pendulum_period(9.8, 1.0)
  assert_eq((result - 2.0).abs() < 0.1, true)
}

// ============================================================================
// 5. 碰撞测试（Collisions）
// ============================================================================

///|
test "elastic_collision_v1 - equal masses" {
  // 质量相等，m1=1 kg, v1=10 m/s, m2=1 kg, v2=0 m/s
  // 物体1末速度 = (1-1)*10/(1+1) + 2*1*0/(1+1) = 0 m/s
  let result = elastic_collision_v1(1.0, 10.0, 1.0, 0.0)
  assert_eq((result - 0.0).abs() < 0.0001, true)
}

///|
test "elastic_collision_v2 - equal masses" {
  // 质量相等，m1=1 kg, v1=10 m/s, m2=1 kg, v2=0 m/s
  // 物体2末速度 = 2*1*10/(1+1) + (1-1)*0/(1+1) = 10 m/s
  let result = elastic_collision_v2(1.0, 10.0, 1.0, 0.0)
  assert_eq((result - 10.0).abs() < 0.0001, true)
}

///|
test "elastic_collision_v1 - heavy object hits light object" {
  // m1=3 kg, v1=5 m/s, m2=1 kg, v2=0 m/s
  // 物体1末速度 = (3-1)*5/4 + 0 = 2.5 m/s
  let result = elastic_collision_v1(3.0, 5.0, 1.0, 0.0)
  assert_eq((result - 2.5).abs() < 0.0001, true)
}

///|
test "inelastic_collision_velocity - basic calculation" {
  // m1=2 kg, v1=10 m/s, m2=3 kg, v2=0 m/s
  // 碰撞后速度 = (2*10 + 3*0)/(2+3) = 20/5 = 4 m/s
  let result = inelastic_collision_velocity(2.0, 10.0, 3.0, 0.0)
  assert_eq(result, 4.0)
}

///|
test "coefficient_of_restitution - perfectly elastic" {
  // 完全弹性碰撞：v1f=0, v2f=10, v1i=10, v2i=0
  // 恢复系数 = (10-0)/(10-0) = 1.0
  let result = coefficient_of_restitution(0.0, 10.0, 10.0, 0.0)
  assert_eq(result, 1.0)
}

///|
test "coefficient_of_restitution - perfectly inelastic" {
  // 完全非弹性碰撞：v1f=5, v2f=5, v1i=10, v2i=0
  // 恢复系数 = (5-5)/(10-0) = 0.0
  let result = coefficient_of_restitution(5.0, 5.0, 10.0, 0.0)
  assert_eq(result, 0.0)
}

// ============================================================================
// 6. 中心力场测试（Central Force Motion）
// ============================================================================

///|
test "kepler_period - Earth orbit" {
  // 地球轨道半长轴 1.496e11 m，太阳质量 1.989e30 kg，G = 6.67430e-11
  // 周期 ≈ 3.156e7 s (约 1 年)
  let a = 1.496e11
  let m = 1.989e30
  let g_const = 6.67430e-11
  let result = kepler_period(a, m, g_const)
  assert_eq((result - 3.156e7).abs() < 0.1e7, true)
}

///|
test "orbital_velocity_circular - Earth orbit" {
  // 轨道半径 1.496e11 m，太阳质量 1.989e30 kg，G = 6.67430e-11
  // 轨道速度 ≈ 29780 m/s
  let r = 1.496e11
  let m = 1.989e30
  let g_const = 6.67430e-11
  let result = orbital_velocity_circular(r, m, g_const)
  assert_eq((result - 29780.0).abs() < 100.0, true)
}

///|
test "escape_velocity - Earth surface" {
  // 地球半径 6.371e6 m，地球质量 5.972e24 kg，G = 6.67430e-11
  // 逃逸速度 ≈ 11186 m/s
  let r = 6.371e6
  let m = 5.972e24
  let g_const = 6.67430e-11
  let result = escape_velocity(r, m, g_const)
  assert_eq((result - 11186.0).abs() < 100.0, true)
}

///|
test "first_cosmic_velocity - Earth surface" {
  // 重力加速度 9.8 m/s²，地球半径 6.371e6 m
  // 第一宇宙速度 = sqrt(9.8*6.371e6) ≈ 7900 m/s
  let g = 9.8
  let r = 6.371e6
  let result = first_cosmic_velocity(g, r)
  assert_eq((result - 7900.0).abs() < 100.0, true)
}

///|
test "second_cosmic_velocity - Earth surface" {
  // 重力加速度 9.8 m/s²，地球半径 6.371e6 m
  // 第二宇宙速度 = sqrt(2*9.8*6.371e6) ≈ 11186 m/s
  let g = 9.8
  let r = 6.371e6
  let result = second_cosmic_velocity(g, r)
  assert_eq((result - 11186.0).abs() < 200.0, true)
}

// ============================================================================
// 7. 拉格朗日力学测试（Lagrangian Mechanics）
// ============================================================================

///|
test "lagrangian - basic calculation" {
  // 动能 100 J，势能 50 J，拉格朗日函数 = 100 - 50 = 50 J
  let result = lagrangian(100.0, 50.0)
  assert_eq(result, 50.0)
}

///|
test "lagrangian_generalized_force - basic calculation" {
  // dL/dq_dot 对时间的导数 = 10，dL/dq = 5
  // 广义力 = 10 - 5 = 5
  let result = lagrangian_generalized_force(10.0, 5.0)
  assert_eq(result, 5.0)
}

// ============================================================================
// 8. 哈密顿力学测试（Hamiltonian Mechanics）
// ============================================================================

///|
test "hamiltonian - conservative system" {
  // 动能 100 J，势能 50 J，哈密顿函数 = 100 + 50 = 150 J
  let result = hamiltonian(100.0, 50.0)
  assert_eq(result, 150.0)
}

///|
test "hamiltonian_from_lagrangian - basic calculation" {
  // 动量 10 kg·m/s，速度 5 m/s，拉格朗日函数 50 J
  // 哈密顿函数 = 10*5 - 50 = 50 - 50 = 0 J
  let result = hamiltonian_from_lagrangian(10.0, 5.0, 50.0)
  assert_eq(result, 0.0)
}

///|
test "hamiltonian_from_lagrangian - non-zero result" {
  // 动量 20 kg·m/s，速度 10 m/s，拉格朗日函数 50 J
  // 哈密顿函数 = 20*10 - 50 = 200 - 50 = 150 J
  let result = hamiltonian_from_lagrangian(20.0, 10.0, 50.0)
  assert_eq(result, 150.0)
}

// ============================================================================
// 9. 多体系统测试（Multi-body Systems）
// ============================================================================

///|
test "center_of_mass_1d - two equal masses" {
  // 两个质量相等的物体：m=[1,1], x=[0,2]
  // 质心 = (1*0 + 1*2)/(1+1) = 1.0 m
  let masses = [1.0, 1.0]
  let positions = [0.0, 2.0]
  let result = center_of_mass_1d(masses, positions)
  assert_eq(result, 1.0)
}

///|
test "center_of_mass_1d - unequal masses" {
  // 两个质量不等的物体：m=[1,3], x=[0,2]
  // 质心 = (1*0 + 3*2)/(1+3) = 6/4 = 1.5 m
  let masses = [1.0, 3.0]
  let positions = [0.0, 2.0]
  let result = center_of_mass_1d(masses, positions)
  assert_eq(result, 1.5)
}

///|
test "center_of_mass_1d - three masses" {
  // 三个物体：m=[1,2,3], x=[0,2,4]
  // 质心 = (1*0 + 2*2 + 3*4)/(1+2+3) = 16/6 ≈ 2.67 m
  let masses = [1.0, 2.0, 3.0]
  let positions = [0.0, 2.0, 4.0]
  let result = center_of_mass_1d(masses, positions)
  assert_eq((result - 2.667).abs() < 0.01, true)
}

///|
test "center_of_mass_1d - empty arrays" {
  // 空数组应返回 0
  let masses = Array::new()
  let positions = Array::new()
  let result = center_of_mass_1d(masses, positions)
  assert_eq(result, 0.0)
}

///|
test "total_momentum_1d - two objects" {
  // 两个物体：m=[2,3], v=[5,10]
  // 总动量 = 2*5 + 3*10 = 10 + 30 = 40 kg·m/s
  let masses = [2.0, 3.0]
  let velocities = [5.0, 10.0]
  let result = total_momentum_1d(masses, velocities)
  assert_eq(result, 40.0)
}

///|
test "total_momentum_1d - opposite velocities" {
  // 两个物体：m=[1,1], v=[10,-10]
  // 总动量 = 1*10 + 1*(-10) = 0 kg·m/s
  let masses = [1.0, 1.0]
  let velocities = [10.0, -10.0]
  let result = total_momentum_1d(masses, velocities)
  assert_eq(result, 0.0)
}

///|
test "total_kinetic_energy_1d - two objects" {
  // 两个物体：m=[2,3], v=[5,10]
  // 总动能 = 0.5*2*25 + 0.5*3*100 = 25 + 150 = 175 J
  let masses = [2.0, 3.0]
  let velocities = [5.0, 10.0]
  let result = total_kinetic_energy_1d(masses, velocities)
  assert_eq(result, 175.0)
}

///|
test "total_kinetic_energy_1d - single object" {
  // 单个物体：m=[5], v=[10]
  // 总动能 = 0.5*5*100 = 250 J
  let masses = [5.0]
  let velocities = [10.0]
  let result = total_kinetic_energy_1d(masses, velocities)
  assert_eq(result, 250.0)
}

// ============================================================================
// 10. 阻尼振动测试（Damped Oscillations）
// ============================================================================

///|
test "damping_coefficient - basic calculation" {
  // 阻尼系数 b = 4 kg/s，质量 m = 2 kg
  // 阻尼系数 gamma = 4/(2*2) = 1.0 s⁻¹
  let result = damping_coefficient(4.0, 2.0)
  assert_eq(result, 1.0)
}

///|
test "critical_damping_coefficient - basic calculation" {
  // 质量 4 kg，弹簧常数 100 N/m
  // 临界阻尼系数 = 2*sqrt(4*100) = 2*20 = 40 kg/s
  let result = critical_damping_coefficient(4.0, 100.0)
  assert_eq(result, 40.0)
}

///|
test "damped_angular_frequency - underdamped" {
  // 固有角频率 5 rad/s，阻尼系数 3 s⁻¹
  // 阻尼角频率 = sqrt(25 - 9) = sqrt(16) = 4 rad/s
  let result = damped_angular_frequency(5.0, 3.0)
  assert_eq((result - 4.0).abs() < 0.0001, true)
}

///|
test "damped_angular_frequency - critical damping" {
  // 固有角频率 5 rad/s，阻尼系数 5 s⁻¹（临界阻尼）
  // 阻尼角频率 = sqrt(25 - 25) = 0 rad/s
  let result = damped_angular_frequency(5.0, 5.0)
  assert_eq((result - 0.0).abs() < 0.0001, true)
}

// ============================================================================
// 11. 综合测试和边界条件
// ============================================================================

///|
test "energy conservation - projectile motion" {
  // 测试抛体运动中机械能守恒
  // 初速度 20 m/s，角度 45°，质量 1 kg，g = 9.8 m/s²
  let v0 = 20.0
  let angle = 3.141592653589793 / 4.0
  let m = 1.0
  let g = 9.8
  let h0 = 0.0

  // 初始动能
  let ke0 = kinetic_energy(m, v0)
  // 初始势能
  let pe0 = gravitational_potential_energy(m, h0, g)
  // 初始机械能
  let me0 = mechanical_energy(ke0, pe0)

  // 最大高度处的动能和势能
  let h_max = projectile_max_height(v0, angle, g)
  let v_max = v0 * @math.cos(angle) // 水平速度不变
  let ke_max = kinetic_energy(m, v_max)
  let pe_max = gravitational_potential_energy(m, h_max, g)
  let me_max = mechanical_energy(ke_max, pe_max)

  // 机械能应该守恒（考虑浮点误差）
  assert_eq((me0 - me_max).abs() < 0.1, true)
}

///|
test "momentum conservation - elastic collision" {
  // 测试弹性碰撞中动量守恒
  // m1=2 kg, v1=10 m/s, m2=3 kg, v2=0 m/s
  let m1 = 2.0
  let v1 = 10.0
  let m2 = 3.0
  let v2 = 0.0

  // 碰撞前总动量
  let p_before = momentum(m1, v1) + momentum(m2, v2)

  // 碰撞后速度
  let v1f = elastic_collision_v1(m1, v1, m2, v2)
  let v2f = elastic_collision_v2(m1, v1, m2, v2)

  // 碰撞后总动量
  let p_after = momentum(m1, v1f) + momentum(m2, v2f)

  // 动量应该守恒
  assert_eq((p_before - p_after).abs() < 0.0001, true)
}

///|
test "kinetic energy conservation - elastic collision" {
  // 测试完全弹性碰撞中动能守恒
  // m1=1 kg, v1=10 m/s, m2=1 kg, v2=0 m/s
  let m1 = 1.0
  let v1 = 10.0
  let m2 = 1.0
  let v2 = 0.0

  // 碰撞前总动能
  let ke_before = kinetic_energy(m1, v1) + kinetic_energy(m2, v2)

  // 碰撞后速度
  let v1f = elastic_collision_v1(m1, v1, m2, v2)
  let v2f = elastic_collision_v2(m1, v1, m2, v2)

  // 碰撞后总动能
  let ke_after = kinetic_energy(m1, v1f) + kinetic_energy(m2, v2f)

  // 动能应该守恒（完全弹性碰撞）
  assert_eq((ke_before - ke_after).abs() < 0.0001, true)
}

///|
test "zero cases - various functions" {
  // 测试各种零值情况
  assert_eq(uniform_motion_displacement(0.0, 10.0), 0.0)
  assert_eq(uniform_acceleration_displacement(0.0, 0.0, 5.0), 0.0)
  assert_eq(momentum(0.0, 10.0), 0.0)
  assert_eq(kinetic_energy(5.0, 0.0), 0.0)
  assert_eq(gravitational_potential_energy(5.0, 0.0, 9.8), 0.0)
  assert_eq(elastic_potential_energy(100.0, 0.0), 0.0)
}

///|
test "relationship checks - angular motion" {
  // 测试角运动量之间的关系
  // 力矩 20 N·m，转动惯量 5 kg·m²
  let tau = 20.0
  let i = 5.0

  // 角加速度
  let alpha = angular_acceleration(tau, i)
  assert_eq(alpha, 4.0)

  // 时间 2 s 后的角速度（假设初角速度为 0）
  let omega = angular_velocity(0.0, alpha, 2.0)
  assert_eq(omega, 8.0)

  // 角动量
  let l = angular_momentum(i, omega)
  assert_eq(l, 40.0)

  // 转动动能
  let ke_rot = rotational_kinetic_energy(i, omega)
  assert_eq(ke_rot, 160.0)
}
