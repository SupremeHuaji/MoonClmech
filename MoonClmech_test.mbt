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

// ============================================================================
// 12. 圆周运动测试（Circular Motion）
// ============================================================================

///|
test "centripetal_acceleration - basic calculation" {
  // 线速度 10 m/s，半径 5 m
  // 向心加速度 = 10²/5 = 100/5 = 20 m/s²
  let result = centripetal_acceleration(10.0, 5.0)
  assert_eq(result, 20.0)
}

///|
test "centripetal_acceleration_angular - basic calculation" {
  // 角速度 2 rad/s，半径 5 m
  // 向心加速度 = 2²*5 = 4*5 = 20 m/s²
  let result = centripetal_acceleration_angular(2.0, 5.0)
  assert_eq(result, 20.0)
}

///|
test "centripetal_force - basic calculation" {
  // 质量 2 kg，线速度 10 m/s，半径 5 m
  // 向心力 = 2*10²/5 = 2*20 = 40 N
  let result = centripetal_force(2.0, 10.0, 5.0)
  assert_eq(result, 40.0)
}

///|
test "centripetal_force_angular - basic calculation" {
  // 质量 2 kg，角速度 2 rad/s，半径 5 m
  // 向心力 = 2*2²*5 = 2*20 = 40 N
  let result = centripetal_force_angular(2.0, 2.0, 5.0)
  assert_eq(result, 40.0)
}

///|
test "linear_velocity_from_angular - basic calculation" {
  // 角速度 5 rad/s，半径 2 m
  // 线速度 = 5*2 = 10 m/s
  let result = linear_velocity_from_angular(5.0, 2.0)
  assert_eq(result, 10.0)
}

///|
test "angular_velocity_from_linear - basic calculation" {
  // 线速度 10 m/s，半径 2 m
  // 角速度 = 10/2 = 5 rad/s
  let result = angular_velocity_from_linear(10.0, 2.0)
  assert_eq(result, 5.0)
}

///|
test "circular_motion_period - basic calculation" {
  // 角速度 2π rad/s，周期 = 2π/(2π) = 1 s
  let omega = 2.0 * 3.141592653589793
  let result = circular_motion_period(omega)
  assert_eq((result - 1.0).abs() < 0.01, true)
}

///|
test "circular_motion_frequency - basic calculation" {
  // 角速度 2π rad/s，频率 = (2π)/(2π) = 1 Hz
  let omega = 2.0 * 3.141592653589793
  let result = circular_motion_frequency(omega)
  assert_eq((result - 1.0).abs() < 0.01, true)
}

// ============================================================================
// 13. 摩擦力测试（Friction）
// ============================================================================

///|
test "static_friction - basic calculation" {
  // 静摩擦系数 0.5，正压力 100 N
  // 最大静摩擦力 = 0.5*100 = 50 N
  let result = static_friction(0.5, 100.0)
  assert_eq(result, 50.0)
}

///|
test "kinetic_friction - basic calculation" {
  // 动摩擦系数 0.3，正压力 100 N
  // 动摩擦力 = 0.3*100 = 30 N
  let result = kinetic_friction(0.3, 100.0)
  assert_eq(result, 30.0)
}

///|
test "normal_force_inclined_plane - basic calculation" {
  // 质量 10 kg，重力加速度 9.8 m/s²，角度 30°
  // 正压力 = 10*9.8*cos(30°) ≈ 84.87 N
  let angle = 3.141592653589793 / 6.0
  let result = normal_force_inclined_plane(10.0, 9.8, angle)
  assert_eq((result - 84.87).abs() < 0.1, true)
}

///|
test "gravity_component_along_incline - basic calculation" {
  // 质量 10 kg，重力加速度 9.8 m/s²，角度 30°
  // 重力分量 = 10*9.8*sin(30°) = 49 N
  let angle = 3.141592653589793 / 6.0
  let result = gravity_component_along_incline(10.0, 9.8, angle)
  assert_eq((result - 49.0).abs() < 0.1, true)
}

///|
test "acceleration_inclined_plane_no_friction - basic calculation" {
  // 重力加速度 9.8 m/s²，角度 30°
  // 加速度 = 9.8*sin(30°) = 4.9 m/s²
  let angle = 3.141592653589793 / 6.0
  let result = acceleration_inclined_plane_no_friction(9.8, angle)
  assert_eq((result - 4.9).abs() < 0.1, true)
}

///|
test "acceleration_inclined_plane_with_friction - basic calculation" {
  // 重力加速度 9.8 m/s²，角度 30°，动摩擦系数 0.2
  // 加速度 = 9.8*(sin(30°) - 0.2*cos(30°)) ≈ 3.20 m/s²
  let angle = 3.141592653589793 / 6.0
  let result = acceleration_inclined_plane_with_friction(9.8, angle, 0.2)
  assert_eq((result - 3.20).abs() < 0.1, true)
}

// ============================================================================
// 14. 功和功率测试（Work and Power）
// ============================================================================

///|
test "work - parallel force" {
  // 力 10 N，位移 5 m，角度 0°（同向）
  // 功 = 10*5*cos(0) = 50 J
  let result = work(10.0, 5.0, 0.0)
  assert_eq((result - 50.0).abs() < 0.0001, true)
}

///|
test "work - perpendicular force" {
  // 力 10 N，位移 5 m，角度 90°（垂直）
  // 功 = 10*5*cos(90°) = 0 J
  let angle = 3.141592653589793 / 2.0
  let result = work(10.0, 5.0, angle)
  assert_eq((result - 0.0).abs() < 0.0001, true)
}

///|
test "work_parallel - basic calculation" {
  // 力 10 N，位移 5 m，功 = 50 J
  let result = work_parallel(10.0, 5.0)
  assert_eq(result, 50.0)
}

///|
test "power - basic calculation" {
  // 功 100 J，时间 5 s，功率 = 20 W
  let result = power(100.0, 5.0)
  assert_eq(result, 20.0)
}

///|
test "power_from_force - parallel" {
  // 力 10 N，速度 5 m/s，角度 0°
  // 功率 = 10*5*cos(0) = 50 W
  let result = power_from_force(10.0, 5.0, 0.0)
  assert_eq((result - 50.0).abs() < 0.0001, true)
}

///|
test "power_from_force_parallel - basic calculation" {
  // 力 10 N，速度 5 m/s，功率 = 50 W
  let result = power_from_force_parallel(10.0, 5.0)
  assert_eq(result, 50.0)
}

///|
test "work_from_kinetic_energy - basic calculation" {
  // 末动能 100 J，初动能 50 J
  // 功 = 100 - 50 = 50 J
  let result = work_from_kinetic_energy(100.0, 50.0)
  assert_eq(result, 50.0)
}

// ============================================================================
// 15. 更多转动惯量测试（Additional Moments of Inertia）
// ============================================================================

///|
test "moment_of_inertia_cylinder_axis - basic calculation" {
  // 质量 4 kg，半径 2 m，转动惯量 = 0.5*4*4 = 8 kg·m²
  let result = moment_of_inertia_cylinder_axis(4.0, 2.0)
  assert_eq(result, 8.0)
}

///|
test "moment_of_inertia_cylinder_diameter - basic calculation" {
  // 质量 6 kg，半径 1 m，高度 2 m
  // 转动惯量 = 6*(1/4 + 4/12) = 6*(0.25 + 0.333) = 3.5 kg·m²
  let result = moment_of_inertia_cylinder_diameter(6.0, 1.0, 2.0)
  assert_eq((result - 3.5).abs() < 0.1, true)
}

///|
test "moment_of_inertia_cone_axis - basic calculation" {
  // 质量 3 kg，半径 2 m，转动惯量 = 0.3*3*4 = 3.6 kg·m²
  let result = moment_of_inertia_cone_axis(3.0, 2.0)
  assert_eq((result - 3.6).abs() < 0.1, true)
}

///|
test "moment_of_inertia_rectangular_plate - basic calculation" {
  // 质量 12 kg，长度 4 m，宽度 3 m
  // 转动惯量 = 12*(16+9)/12 = 12*25/12 = 25 kg·m²
  let result = moment_of_inertia_rectangular_plate(12.0, 4.0, 3.0)
  assert_eq(result, 25.0)
}

///|
test "moment_of_inertia_rectangular_plate_edge - basic calculation" {
  // 质量 6 kg，宽度 2 m，转动惯量 = 6*4/3 = 8 kg·m²
  let result = moment_of_inertia_rectangular_plate_edge(6.0, 2.0)
  assert_eq(result, 8.0)
}

// ============================================================================
// 16. 刚体平衡测试（Rigid Body Equilibrium）
// ============================================================================

///|
test "lever_equilibrium - balanced" {
  // 力1=10 N，力臂1=2 m，力2=20 N，力臂2=1 m
  // 力矩平衡：10*2 = 20*1 = 20 N·m
  let result = lever_equilibrium(10.0, 2.0, 20.0, 1.0)
  assert_eq(result, true)
}

///|
test "lever_equilibrium_force - basic calculation" {
  // 力1=10 N，力臂1=2 m，力臂2=1 m
  // 平衡所需的力2 = 10*2/1 = 20 N
  let result = lever_equilibrium_force(10.0, 2.0, 1.0)
  assert_eq(result, 20.0)
}

///|
test "net_force_1d - basic calculation" {
  // 力数组 [10, -5, 3]
  // 合力 = 10 - 5 + 3 = 8 N
  let forces = [10.0, -5.0, 3.0]
  let result = net_force_1d(forces)
  assert_eq(result, 8.0)
}

///|
test "net_torque - basic calculation" {
  // 力矩数组 [10, -5, 3]
  // 合力矩 = 10 - 5 + 3 = 8 N·m
  let torques = [10.0, -5.0, 3.0]
  let result = net_torque(torques)
  assert_eq(result, 8.0)
}

// ============================================================================
// 17. 相对运动测试（Relative Motion）
// ============================================================================

///|
test "relative_velocity_1d - basic calculation" {
  // 物体1速度 5 m/s，物体2速度 10 m/s
  // 物体2相对于物体1的速度 = 10 - 5 = 5 m/s
  let result = relative_velocity_1d(5.0, 10.0)
  assert_eq(result, 5.0)
}

///|
test "relative_displacement_1d - basic calculation" {
  // 物体1位置 0 m，物体2位置 5 m
  // 物体2相对于物体1的位移 = 5 - 0 = 5 m
  let result = relative_displacement_1d(0.0, 5.0)
  assert_eq(result, 5.0)
}

///|
test "relative_acceleration_1d - basic calculation" {
  // 物体1加速度 2 m/s²，物体2加速度 5 m/s²
  // 物体2相对于物体1的加速度 = 5 - 2 = 3 m/s²
  let result = relative_acceleration_1d(2.0, 5.0)
  assert_eq(result, 3.0)
}

// ============================================================================
// 18. 流体力学基础测试（Fluid Mechanics Basics）
// ============================================================================

///|
test "buoyant_force - basic calculation" {
  // 流体密度 1000 kg/m³，排开体积 0.1 m³，重力加速度 9.8 m/s²
  // 浮力 = 1000*0.1*9.8 = 980 N
  let result = buoyant_force(1000.0, 0.1, 9.8)
  assert_eq((result - 980.0).abs() < 0.0001, true)
}

///|
test "hydrostatic_pressure - basic calculation" {
  // 流体密度 1000 kg/m³，深度 10 m，重力加速度 9.8 m/s²，表面压力 101325 Pa
  // 静压力 = 101325 + 1000*9.8*10 = 101325 + 98000 = 199325 Pa
  let result = hydrostatic_pressure(1000.0, 10.0, 9.8, 101325.0)
  assert_eq(result, 199325.0)
}

///|
test "bernoulli_velocity - basic calculation" {
  // 压力1=200000 Pa，压力2=150000 Pa，密度=1000 kg/m³，速度1=5 m/s
  // 速度2 = sqrt(25 + 2*50000/1000) = sqrt(125) ≈ 11.18 m/s
  let result = bernoulli_velocity(200000.0, 150000.0, 1000.0, 5.0)
  assert_eq((result - 11.18).abs() < 0.1, true)
}

///|
test "bernoulli_pressure - basic calculation" {
  // 压力1=200000 Pa，密度=1000 kg/m³，速度1=10 m/s，速度2=5 m/s
  // 压力2 = 200000 + 0.5*1000*(100-25) = 200000 + 37500 = 237500 Pa
  let result = bernoulli_pressure(200000.0, 1000.0, 10.0, 5.0)
  assert_eq(result, 237500.0)
}

///|
test "continuity_equation_velocity - basic calculation" {
  // 面积1=0.1 m²，速度1=5 m/s，面积2=0.05 m²
  // 速度2 = 0.1*5/0.05 = 10 m/s
  let result = continuity_equation_velocity(0.1, 5.0, 0.05)
  assert_eq(result, 10.0)
}

// ============================================================================
// 19. 能量守恒应用测试（Energy Conservation Applications）
// ============================================================================

///|
test "velocity_from_energy_conservation - basic calculation" {
  // 初动能 100 J，初势能 50 J，末势能 100 J，质量 2 kg
  // 末动能 = 100 + 50 - 100 = 50 J
  // 末速度 = sqrt(2*50/2) = sqrt(50) ≈ 7.07 m/s
  let result = velocity_from_energy_conservation(100.0, 50.0, 100.0, 2.0)
  assert_eq((result - 7.07).abs() < 0.1, true)
}

///|
test "max_height_from_energy - basic calculation" {
  // 初速度 20 m/s，初始高度 0 m，重力加速度 9.8 m/s²
  // 最大高度 = 0 + 400/(2*9.8) ≈ 20.41 m
  let result = max_height_from_energy(20.0, 0.0, 9.8)
  assert_eq((result - 20.41).abs() < 0.1, true)
}

///|
test "velocity_from_spring_energy - basic calculation" {
  // 弹簧常数 100 N/m，压缩量 0.5 m，质量 2 kg
  // 弹性势能 = 0.5*100*0.25 = 12.5 J
  // 速度 = sqrt(2*12.5/2) = sqrt(12.5) ≈ 3.54 m/s
  let result = velocity_from_spring_energy(100.0, 0.5, 2.0)
  assert_eq((result - 3.54).abs() < 0.1, true)
}

// ============================================================================
// 20. 角动量守恒测试（Angular Momentum Conservation）
// ============================================================================

///|
test "angular_velocity_from_conservation - basic calculation" {
  // 初转动惯量 10 kg·m²，初角速度 5 rad/s，末转动惯量 5 kg·m²
  // 末角速度 = 10*5/5 = 10 rad/s
  let result = angular_velocity_from_conservation(10.0, 5.0, 5.0)
  assert_eq(result, 10.0)
}

///|
test "moment_of_inertia_from_conservation - basic calculation" {
  // 初转动惯量 10 kg·m²，初角速度 5 rad/s，末角速度 10 rad/s
  // 末转动惯量 = 10*5/10 = 5 kg·m²
  let result = moment_of_inertia_from_conservation(10.0, 5.0, 10.0)
  assert_eq(result, 5.0)
}

// ============================================================================
// 21. 受迫振动和共振测试（Forced Oscillations and Resonance）
// ============================================================================

///|
test "forced_oscillation_amplitude - basic calculation" {
  // 驱动力幅值 10 N，质量 2 kg，固有角频率 5 rad/s，驱动角频率 4 rad/s，阻尼系数 1 s⁻¹
  // 振幅计算（简化验证）
  let result = forced_oscillation_amplitude(10.0, 2.0, 5.0, 4.0, 1.0)
  assert_eq(result > 0.0, true)
}

///|
test "resonance_frequency - basic calculation" {
  // 固有角频率 5 rad/s，共振频率 = 5 rad/s
  let result = resonance_frequency(5.0)
  assert_eq(result, 5.0)
}

///|
test "resonance_frequency_damped - basic calculation" {
  // 固有角频率 5 rad/s，阻尼系数 3 s⁻¹
  // 共振频率 = sqrt(25 - 9) = sqrt(16) = 4 rad/s
  let result = resonance_frequency_damped(5.0, 3.0)
  assert_eq((result - 4.0).abs() < 0.0001, true)
}

///|
test "normal_frequency_antisymmetric - basic calculation" {
  // 固有角频率 5 rad/s，耦合常数 10，质量 2 kg
  // 反对称模式频率 = sqrt(25 + 2*10/2) = sqrt(35) ≈ 5.92 rad/s
  let result = normal_frequency_antisymmetric(5.0, 10.0, 2.0)
  assert_eq((result - 5.92).abs() < 0.1, true)
}

// ============================================================================
// 23. 质心运动定理测试（Center of Mass Motion）
// ============================================================================

///|
test "center_of_mass_velocity_1d - basic calculation" {
  // 质量数组 [2, 3]，速度数组 [5, 10]
  // 质心速度 = (2*5 + 3*10)/(2+3) = 40/5 = 8 m/s
  let masses = [2.0, 3.0]
  let velocities = [5.0, 10.0]
  let result = center_of_mass_velocity_1d(masses, velocities)
  assert_eq(result, 8.0)
}

///|
test "center_of_mass_acceleration_1d - basic calculation" {
  // 质量数组 [2, 3]，加速度数组 [4, 6]
  // 质心加速度 = (2*4 + 3*6)/(2+3) = 26/5 = 5.2 m/s²
  let masses = [2.0, 3.0]
  let accelerations = [4.0, 6.0]
  let result = center_of_mass_acceleration_1d(masses, accelerations)
  assert_eq((result - 5.2).abs() < 0.01, true)
}

// ============================================================================
// 24. 碰撞的动能损失测试（Kinetic Energy Loss in Collisions）
// ============================================================================

///|
test "kinetic_energy_loss_inelastic - basic calculation" {
  // m1=2 kg, v1=10 m/s, m2=3 kg, v2=0 m/s
  // 初动能 = 0.5*2*100 + 0 = 100 J
  // 末速度 = (2*10 + 3*0)/5 = 4 m/s
  // 末动能 = 0.5*5*16 = 40 J
  // 动能损失 = 100 - 40 = 60 J
  let result = kinetic_energy_loss_inelastic(2.0, 10.0, 3.0, 0.0)
  assert_eq((result - 60.0).abs() < 0.1, true)
}

///|
test "kinetic_energy_loss_partial - basic calculation" {
  // m1=1 kg, v1=10 m/s, m2=1 kg, v2=0 m/s, 恢复系数 0.5
  // 动能损失应该大于0（验证非弹性碰撞）
  let result = kinetic_energy_loss_partial(1.0, 10.0, 1.0, 0.0, 0.5)
  assert_eq(result > 0.0, true)
}

// ============================================================================
// 25. 滚动运动测试（Rolling Motion）
// ============================================================================

///|
test "rolling_velocity - basic calculation" {
  // 角速度 5 rad/s，半径 2 m
  // 线速度 = 5*2 = 10 m/s
  let result = rolling_velocity(5.0, 2.0)
  assert_eq(result, 10.0)
}

///|
test "rolling_kinetic_energy - basic calculation" {
  // 质量 2 kg，线速度 5 m/s，转动惯量 4 kg·m²，角速度 2.5 rad/s
  // 平动动能 = 0.5*2*25 = 25 J
  // 转动动能 = 0.5*4*6.25 = 12.5 J
  // 总动能 = 37.5 J
  let result = rolling_kinetic_energy(2.0, 5.0, 4.0, 2.5)
  assert_eq(result, 37.5)
}

///|
test "rolling_kinetic_energy_simplified - basic calculation" {
  // 质量 2 kg，线速度 5 m/s，转动惯量系数 0.5
  // 总动能 = 0.5*2*25*(1+0.5) = 25*1.5 = 37.5 J
  let result = rolling_kinetic_energy_simplified(2.0, 5.0, 0.5)
  assert_eq(result, 37.5)
}

///|
test "rolling_acceleration_downhill - basic calculation" {
  // 重力加速度 9.8 m/s²，角度 30°，转动惯量系数 0.5
  // 加速度 = 9.8*sin(30°)/(1+0.5) = 4.9/1.5 ≈ 3.27 m/s²
  let angle = 3.141592653589793 / 6.0
  let result = rolling_acceleration_downhill(9.8, angle, 0.5)
  assert_eq((result - 3.27).abs() < 0.1, true)
}

// ============================================================================
// 26. 非惯性参考系测试（Non-inertial Reference Frames）
// ============================================================================

///|
test "centrifugal_force - basic calculation" {
  // 质量 2 kg，角速度 3 rad/s，半径 2 m
  // 离心力 = 2*9*2 = 36 N
  let result = centrifugal_force(2.0, 3.0, 2.0)
  assert_eq(result, 36.0)
}

///|
test "coriolis_force - basic calculation" {
  // 质量 2 kg，相对速度 5 m/s，角速度 3 rad/s
  // 科里奥利力 = 2*2*5*3 = 60 N
  let result = coriolis_force(2.0, 5.0, 3.0)
  assert_eq(result, 60.0)
}

///|
test "inertial_force - basic calculation" {
  // 质量 5 kg，参考系加速度 4 m/s²
  // 惯性力 = -5*4 = -20 N
  let result = inertial_force(5.0, 4.0)
  assert_eq(result, -20.0)
}

// ============================================================================
// 27. 刚体动力学扩展测试（Extended Rigid Body Dynamics）
// ============================================================================

///|
test "precession_angular_velocity - basic calculation" {
  // 质量 1 kg，重力加速度 9.8 m/s²，距离 0.1 m，转动惯量 0.1 kg·m²，自转角速度 10 rad/s
  // 进动角速度 = 1*9.8*0.1/(0.1*10) = 0.98 rad/s
  let result = precession_angular_velocity(1.0, 9.8, 0.1, 0.1, 10.0)
  assert_eq((result - 0.98).abs() < 0.01, true)
}

///|
test "rotational_power - basic calculation" {
  // 力矩 20 N·m，角速度 5 rad/s
  // 功率 = 20*5 = 100 W
  let result = rotational_power(20.0, 5.0)
  assert_eq(result, 100.0)
}

// ============================================================================
// 28. 轨道力学扩展测试（Extended Orbital Mechanics）
// ============================================================================

///|
test "orbital_energy_circular - basic calculation" {
  // 质量 1000 kg，中心质量 5.972e24 kg，半径 6.371e6 m，G = 6.67430e-11
  // 轨道能量 = -6.67430e-11*1000*5.972e24/(2*6.371e6) ≈ -3.13e10 J
  let result = orbital_energy_circular(1000.0, 5.972e24, 6.371e6, 6.67430e-11)
  assert_eq(result < 0.0, true)
}

///|
test "orbital_energy_elliptical - basic calculation" {
  // 质量 1000 kg，中心质量 5.972e24 kg，半长轴 7.0e6 m，G = 6.67430e-11
  // 轨道能量 < 0（椭圆轨道）
  let result = orbital_energy_elliptical(1000.0, 5.972e24, 7.0e6, 6.67430e-11)
  assert_eq(result < 0.0, true)
}

///|
test "areal_velocity - basic calculation" {
  // 距离 1.496e11 m，垂直速度 29780 m/s
  // 面积速度 = 0.5*1.496e11*29780 ≈ 2.23e15 m²/s
  let result = areal_velocity(1.496e11, 29780.0)
  assert_eq((result - 2.23e15).abs() < 0.1e15, true)
}

// ============================================================================
// 29. 波动测试（Waves）
// ============================================================================

///|
test "wave_displacement - basic calculation" {
  // 振幅 5 m，波数 2 rad/m，位置 0 m，角频率 4 rad/s，时间 0 s，初相位 0
  // 位移 = 5*cos(0) = 5 m
  let result = wave_displacement(5.0, 2.0, 0.0, 4.0, 0.0, 0.0)
  assert_eq((result - 5.0).abs() < 0.0001, true)
}

///|
test "wave_velocity - basic calculation" {
  // 频率 2 Hz，波长 5 m，波速 = 2*5 = 10 m/s
  let result = wave_velocity(2.0, 5.0)
  assert_eq(result, 10.0)
}

///|
test "wave_velocity_from_angular - basic calculation" {
  // 角频率 10 rad/s，波数 2 rad/m，波速 = 10/2 = 5 m/s
  let result = wave_velocity_from_angular(10.0, 2.0)
  assert_eq(result, 5.0)
}

///|
test "doppler_effect_source_moving - source approaching" {
  // 声源频率 1000 Hz，声速 340 m/s，声源速度 20 m/s（朝向观察者）
  // 观察频率 = 1000*340/(340-20) = 340000/320 ≈ 1062.5 Hz
  let result = doppler_effect_source_moving(1000.0, 340.0, 20.0)
  assert_eq((result - 1062.5).abs() < 1.0, true)
}

///|
test "doppler_effect_observer_moving - observer approaching" {
  // 声源频率 1000 Hz，声速 340 m/s，观察者速度 20 m/s（朝向声源）
  // 观察频率 = 1000*(340+20)/340 = 360000/340 ≈ 1058.8 Hz
  let result = doppler_effect_observer_moving(1000.0, 340.0, 20.0)
  assert_eq((result - 1058.8).abs() < 1.0, true)
}

// ============================================================================
// 30. 流体力学扩展测试（Extended Fluid Mechanics）
// ============================================================================

///|
test "reynolds_number - basic calculation" {
  // 密度 1000 kg/m³，流速 2 m/s，特征长度 0.1 m，动力粘度 0.001 Pa·s
  // 雷诺数 = 1000*2*0.1/0.001 = 200000
  let result = reynolds_number(1000.0, 2.0, 0.1, 0.001)
  assert_eq(result, 200000.0)
}

///|
test "stokes_drag - basic calculation" {
  // 半径 0.01 m，速度 0.1 m/s，动力粘度 0.001 Pa·s
  // 斯托克斯阻力 = 6*π*0.001*0.01*0.1 ≈ 1.88e-5 N
  let result = stokes_drag(0.01, 0.1, 0.001)
  assert_eq(result > 0.0, true)
}

///|
test "viscous_force - basic calculation" {
  // 动力粘度 0.001 Pa·s，面积 0.1 m²，速度梯度 10 s⁻¹
  // 粘性力 = 0.001*0.1*10 = 0.001 N
  let result = viscous_force(0.001, 0.1, 10.0)
  assert_eq(result, 0.001)
}

///|
test "poiseuille_flow_rate - basic calculation" {
  // 压力差 1000 Pa，半径 0.01 m，长度 1 m，动力粘度 0.001 Pa·s
  // 流量 = π*1000*0.01^4/(8*0.001*1) ≈ 3.93e-6 m³/s
  let result = poiseuille_flow_rate(1000.0, 0.01, 1.0, 0.001)
  assert_eq(result > 0.0, true)
}

// ============================================================================
// 31. 能量和动量扩展测试（Extended Energy and Momentum）
// ============================================================================

///|
test "velocity_from_impulse - basic calculation" {
  // 质量 2 kg，初速度 5 m/s，冲量 10 N·s
  // 末速度 = 5 + 10/2 = 10 m/s
  let result = velocity_from_impulse(2.0, 5.0, 10.0)
  assert_eq(result, 10.0)
}

///|
test "rocket_equation - basic calculation" {
  // 初始速度 0 m/s，排气速度 2000 m/s，初始质量 1000 kg，最终质量 500 kg
  // 最终速度 = 0 + 2000*ln(1000/500) = 2000*ln(2) ≈ 1386.3 m/s
  let result = rocket_equation(0.0, 2000.0, 1000.0, 500.0)
  assert_eq((result - 1386.3).abs() < 10.0, true)
}

///|
test "work_from_power - basic calculation" {
  // 功率 100 W，时间 5 s，功 = 100*5 = 500 J
  let result = work_from_power(100.0, 5.0)
  assert_eq(result, 500.0)
}

// ============================================================================
// 32. 振动系统扩展测试（Extended Vibration Systems）
// ============================================================================

///|
test "damping_ratio - basic calculation" {
  // 阻尼系数 2 s⁻¹，固有角频率 5 rad/s
  // 阻尼比 = 2/5 = 0.4
  let result = damping_ratio(2.0, 5.0)
  assert_eq(result, 0.4)
}

///|
test "physical_pendulum_angular_frequency - basic calculation" {
  // 质量 1 kg，重力加速度 9.8 m/s²，转动惯量 0.1 kg·m²，距离 0.1 m
  // 角频率 = sqrt(1*9.8*0.1/0.1) = sqrt(9.8) ≈ 3.13 rad/s
  let result = physical_pendulum_angular_frequency(1.0, 9.8, 0.1, 0.1)
  assert_eq((result - 3.13).abs() < 0.1, true)
}

///|
test "physical_pendulum_period - basic calculation" {
  // 质量 1 kg，重力加速度 9.8 m/s²，转动惯量 0.1 kg·m²，距离 0.1 m
  // 周期 = 2π*sqrt(0.1/(1*9.8*0.1)) = 2π*sqrt(0.102) ≈ 2.0 s
  let result = physical_pendulum_period(1.0, 9.8, 0.1, 0.1)
  assert_eq((result - 2.0).abs() < 0.1, true)
}

///|
test "double_pendulum_upper_frequency - basic calculation" {
  // 重力加速度 9.8 m/s²，上摆长度 1 m
  // 角频率 = sqrt(9.8/1) ≈ 3.13 rad/s
  let result = double_pendulum_upper_frequency(9.8, 1.0)
  assert_eq((result - 3.13).abs() < 0.1, true)
}

// ============================================================================
// 33. 刚体静力学测试（Rigid Body Statics）
// ============================================================================

///|
test "normal_force_horizontal - basic calculation" {
  // 质量 10 kg，重力加速度 9.8 m/s²
  // 支撑力 = 10*9.8 = 98 N
  let result = normal_force_horizontal(10.0, 9.8)
  assert_eq(result, 98.0)
}

///|
test "equilibrium_net_force_1d - balanced" {
  // 力数组 [10, -5, -5]
  // 合力 = 0，应该平衡
  let forces = [10.0, -5.0, -5.0]
  let result = equilibrium_net_force_1d(forces)
  assert_eq(result, true)
}

///|
test "equilibrium_net_torque - balanced" {
  // 力矩数组 [10, -5, -5]
  // 合力矩 = 0，应该平衡
  let torques = [10.0, -5.0, -5.0]
  let result = equilibrium_net_torque(torques)
  assert_eq(result, true)
}

// ============================================================================
// 34. 运动学扩展测试（Extended Kinematics）
// ============================================================================

///|
test "time_from_velocity_change - basic calculation" {
  // 初速度 0 m/s，末速度 20 m/s，加速度 5 m/s²
  // 时间 = (20-0)/5 = 4 s
  let result = time_from_velocity_change(0.0, 20.0, 5.0)
  assert_eq(result, 4.0)
}

///|
test "acceleration_from_displacement - basic calculation" {
  // 初速度 0 m/s，末速度 10 m/s，位移 25 m
  // 加速度 = (100-0)/(2*25) = 2 m/s²
  let result = acceleration_from_displacement(0.0, 10.0, 25.0)
  assert_eq(result, 2.0)
}

///|
test "projectile_velocity_magnitude - basic calculation" {
  // 初速度 20 m/s，角度 45°，时间 1 s，重力加速度 9.8 m/s²
  // 速度大小应该小于初速度
  let angle = 3.141592653589793 / 4.0
  let result = projectile_velocity_magnitude(20.0, angle, 1.0, 9.8)
  assert_eq(result < 20.0, true)
}

///|
test "projectile_velocity_angle - basic calculation" {
  // 初速度 20 m/s，角度 45°，时间 1 s，重力加速度 9.8 m/s²
  // 速度角度应该小于45°
  let angle = 3.141592653589793 / 4.0
  let result = projectile_velocity_angle(20.0, angle, 1.0, 9.8)
  assert_eq(result < angle, true)
}

// ============================================================================
// 35. 碰撞力学测试（Collision Mechanics）
// ============================================================================

///|
test "elastic_collision_1d_v1f - equal masses" {
  // 质量相等，m1=1 kg, m2=1 kg, v1i=10 m/s, v2i=0 m/s
  // 物体1末速度 = (1-1)*10/2 + 2*1*0/2 = 0 m/s
  let result = elastic_collision_1d_v1f(1.0, 1.0, 10.0, 0.0)
  assert_eq((result - 0.0).abs() < 0.0001, true)
}

///|
test "elastic_collision_1d_v2f - equal masses" {
  // 质量相等，m1=1 kg, m2=1 kg, v1i=10 m/s, v2i=0 m/s
  // 物体2末速度 = 2*1*10/2 + (1-1)*0/2 = 10 m/s
  let result = elastic_collision_1d_v2f(1.0, 1.0, 10.0, 0.0)
  assert_eq((result - 10.0).abs() < 0.0001, true)
}

///|
test "perfectly_inelastic_collision_velocity - basic calculation" {
  // m1=2 kg, v1i=10 m/s, m2=3 kg, v2i=0 m/s
  // 碰撞后速度 = (2*10 + 3*0)/5 = 4 m/s
  let result = perfectly_inelastic_collision_velocity(2.0, 3.0, 10.0, 0.0)
  assert_eq(result, 4.0)
}

///|
test "inelastic_collision_1d_v1f - basic calculation" {
  // m1=2 kg, m2=3 kg, v1i=10 m/s, v2i=0 m/s, 恢复系数 0.5
  // 验证结果合理性
  let result = inelastic_collision_1d_v1f(2.0, 3.0, 10.0, 0.0, 0.5)
  assert_eq(result >= 0.0 && result <= 10.0, true)
}

///|
test "generalized_force_1d - basic calculation" {
  // 势能 50 J，广义坐标 5 m
  // 广义力 = -50/5 = -10 N（简化版本）
  let result = generalized_force_1d(50.0, 5.0)
  assert_eq(result, -10.0)
}

// ============================================================================
// 37. 哈密顿力学扩展测试（Extended Hamiltonian Mechanics）
// ============================================================================

///|
test "hamiltonian_q_dot_1d - basic calculation" {
  // 动量 10 kg·m/s
  // q_dot = 10（简化版本，实际应为 p/m）
  let result = hamiltonian_q_dot_1d(10.0)
  assert_eq(result, 10.0)
}

// ============================================================================
// 38. 狭义相对论力学测试（Special Relativity）
// ============================================================================

///|
test "lorentz_factor - low velocity" {
  // 速度 3e6 m/s，光速 3e8 m/s
  // 洛伦兹因子 ≈ 1.00005（低速近似）
  let result = lorentz_factor(3.0e6, 3.0e8)
  assert_eq((result - 1.0).abs() < 0.01, true)
}

///|
test "relativistic_mass - low velocity" {
  // 静质量 1 kg，速度 3e6 m/s，光速 3e8 m/s
  // 相对论质量 ≈ 1.00005 kg（低速近似）
  let result = relativistic_mass(1.0, 3.0e6, 3.0e8)
  assert_eq((result - 1.0).abs() < 0.01, true)
}

///|
test "relativistic_momentum - basic calculation" {
  // 静质量 1 kg，速度 3e6 m/s，光速 3e8 m/s
  // 相对论动量 ≈ 3e6 kg·m/s（低速近似）
  let result = relativistic_momentum(1.0, 3.0e6, 3.0e8)
  assert_eq((result - 3.0e6).abs() < 1000.0, true)
}

///|
test "mass_energy_equivalence - basic calculation" {
  // 质量 1 kg，光速 3e8 m/s
  // 能量 = 1*(3e8)² = 9e16 J
  let result = mass_energy_equivalence(1.0, 3.0e8)
  assert_eq(result, 9.0e16)
}

// ============================================================================
// 39. 连续介质力学测试（Continuum Mechanics）
// ============================================================================

///|
test "normal_stress - basic calculation" {
  // 力 1000 N，面积 0.01 m²
  // 正应力 = 1000/0.01 = 100000 Pa
  let result = normal_stress(1000.0, 0.01)
  assert_eq(result, 100000.0)
}

///|
test "shear_stress - basic calculation" {
  // 切向力 500 N，面积 0.01 m²
  // 剪应力 = 500/0.01 = 50000 Pa
  let result = shear_stress(500.0, 0.01)
  assert_eq(result, 50000.0)
}

///|
test "normal_strain - basic calculation" {
  // 长度变化 0.01 m，原始长度 1 m
  // 正应变 = 0.01/1 = 0.01
  let result = normal_strain(0.01, 1.0)
  assert_eq(result, 0.01)
}

///|
test "youngs_modulus - basic calculation" {
  // 应力 100000 Pa，应变 0.001
  // 杨氏模量 = 100000/0.001 = 1e8 Pa
  let result = youngs_modulus(100000.0, 0.001)
  assert_eq(result, 1.0e8)
}

///|
test "hookes_law_force - basic calculation" {
  // 弹簧常数 100 N/m，位移 0.5 m
  // 弹性恢复力 = -100*0.5 = -50 N
  let result = hookes_law_force(100.0, 0.5)
  assert_eq(result, -50.0)
}

// ============================================================================
// 40. 声学测试（Acoustics）
// ============================================================================

///|
test "sound_intensity - basic calculation" {
  // 声压振幅 1 Pa，密度 1.2 kg/m³，声速 340 m/s
  // 声强 = 1²/(2*1.2*340) ≈ 0.001225 W/m²
  let result = sound_intensity(1.0, 1.2, 340.0)
  assert_eq((result - 0.001225).abs() < 0.0001, true)
}

///|
test "sound_intensity_level - basic calculation" {
  // 声强 0.001 W/m²，参考声强 1e-12 W/m²
  // 声强级 = 10*log10(0.001/1e-12) = 10*log10(1e9) = 90 dB
  let result = sound_intensity_level(0.001, 1.0e-12)
  assert_eq((result - 90.0).abs() < 1.0, true)
}

///|
test "sound_speed_ideal_gas - basic calculation" {
  // 比热容比 1.4，气体常数 287 J/(kg·K)，温度 300 K
  // 声速 = sqrt(1.4*287*300) ≈ 347 m/s
  let result = sound_speed_ideal_gas(1.4, 287.0, 300.0)
  assert_eq((result - 347.0).abs() < 5.0, true)
}

// ============================================================================
// 41. 耦合振荡器测试（Coupled Oscillators）
// ============================================================================

///|
test "normal_mode_frequency_symmetric - basic calculation" {
  // 弹簧常数 100 N/m，耦合弹簧常数 50 N/m，质量 4 kg
  // 对称模式频率 = sqrt((100+2*50)/4) = sqrt(50) ≈ 7.07 rad/s
  let result = normal_mode_frequency_symmetric(100.0, 50.0, 4.0)
  assert_eq((result - 7.07).abs() < 0.1, true)
}

///|
test "normal_mode_frequency_antisymmetric - basic calculation" {
  // 弹簧常数 100 N/m，耦合弹簧常数 50 N/m，质量 4 kg
  // 反对称模式频率 = sqrt(100/4) = 5 rad/s
  let result = normal_mode_frequency_antisymmetric(100.0, 4.0)
  assert_eq(result, 5.0)
}

///|
test "normal_coordinate_symmetric - basic calculation" {
  // 振子1位移 5 m，振子2位移 3 m
  // 对称坐标 = (5+3)/2 = 4 m
  let result = normal_coordinate_symmetric(5.0, 3.0)
  assert_eq(result, 4.0)
}

///|
test "normal_coordinate_antisymmetric - basic calculation" {
  // 振子1位移 5 m，振子2位移 3 m
  // 反对称坐标 = (5-3)/2 = 1 m
  let result = normal_coordinate_antisymmetric(5.0, 3.0)
  assert_eq(result, 1.0)
}

// ============================================================================
// 42. 刚体运动学扩展测试（Extended Rigid Body Kinematics）
// ============================================================================

///|
test "angular_acceleration_from_velocity - basic calculation" {
  // 初始角速度 2 rad/s，末角速度 10 rad/s，时间 4 s
  // 角加速度 = (10-2)/4 = 2 rad/s²
  let result = angular_acceleration_from_velocity(2.0, 10.0, 4.0)
  assert_eq(result, 2.0)
}

///|
test "angular_displacement_constant_acc - basic calculation" {
  // 初始角速度 1 rad/s，角加速度 2 rad/s²，时间 3 s
  // 角位移 = 1*3 + 0.5*2*9 = 3 + 9 = 12 rad
  let result = angular_displacement_constant_acc(1.0, 2.0, 3.0)
  assert_eq(result, 12.0)
}

///|
test "angular_velocity_constant_acc - basic calculation" {
  // 初始角速度 2 rad/s，角加速度 3 rad/s²，时间 4 s
  // 角速度 = 2 + 3*4 = 14 rad/s
  let result = angular_velocity_constant_acc(2.0, 3.0, 4.0)
  assert_eq(result, 14.0)
}

// ============================================================================
// 43. 相空间和守恒量测试（Phase Space and Conserved Quantities）
// ============================================================================

///|
test "phase_space_volume_1d - basic calculation" {
  // 位置 5 m，动量 10 kg·m/s
  // 相空间体积 = 5*10 = 50
  let result = phase_space_volume_1d(5.0, 10.0)
  assert_eq(result, 50.0)
}

// ============================================================================
// 44. 振动耦合系统测试（Vibrational Coupling Systems）
// ============================================================================

///|
test "beat_frequency - basic calculation" {
  // 频率1 100 Hz，频率2 102 Hz
  // 拍频 = |100-102| = 2 Hz
  let result = beat_frequency(100.0, 102.0)
  assert_eq(result, 2.0)
}

///|
test "beat_period - basic calculation" {
  // 频率1 100 Hz，频率2 102 Hz
  // 拍周期 = 1/2 = 0.5 s
  let result = beat_period(100.0, 102.0)
  assert_eq(result, 0.5)
}

// ============================================================================
// 45. 流体动力学扩展测试（Extended Fluid Dynamics）
// ============================================================================

///|
test "mach_number - basic calculation" {
  // 流速 680 m/s，声速 340 m/s
  // 马赫数 = 680/340 = 2.0
  let result = mach_number(680.0, 340.0)
  assert_eq(result, 2.0)
}

///|
test "venturi_velocity_ratio - basic calculation" {
  // 截面1面积 0.1 m²，截面2面积 0.05 m²
  // 速度比 = 0.1/0.05 = 2.0
  let result = venturi_velocity_ratio(0.1, 0.05)
  assert_eq(result, 2.0)
}

// ============================================================================
// 46. 单位转换测试：能量（Unit Conversions: Energy）
// ============================================================================

///|
test "ev_to_joule - basic calculation" {
  // 1 eV = 1.602176634e-19 J
  // 1000 eV = 1.602176634e-16 J
  let result = ev_to_joule(1000.0)
  assert_eq((result - 1.602176634e-16).abs() < 1.0e-25, true)
}

///|
test "joule_to_ev - basic calculation" {
  // 1.602176634e-19 J = 1 eV
  let result = joule_to_ev(1.602176634e-19)
  assert_eq((result - 1.0).abs() < 0.0001, true)
}

///|
test "cal_to_joule - basic calculation" {
  // 1 cal = 4.184 J
  // 100 cal = 418.4 J
  let result = cal_to_joule(100.0)
  assert_eq((result - 418.4).abs() < 0.0001, true)
}

///|
test "joule_to_cal - basic calculation" {
  // 4.184 J = 1 cal
  let result = joule_to_cal(4.184)
  assert_eq((result - 1.0).abs() < 0.0001, true)
}

///|
test "btu_to_joule - basic calculation" {
  // 1 BTU = 1055.05585262 J
  // 2 BTU = 2110.11170524 J
  let result = btu_to_joule(2.0)
  assert_eq((result - 2110.11170524).abs() < 0.01, true)
}

///|
test "joule_to_btu - basic calculation" {
  // 1055.05585262 J = 1 BTU
  let result = joule_to_btu(1055.05585262)
  assert_eq((result - 1.0).abs() < 0.0001, true)
}

// ============================================================================
// 47. 单位转换测试：功率（Unit Conversions: Power）
// ============================================================================

///|
test "hp_to_watt - basic calculation" {
  // 1 hp = 745.699872 W
  // 2 hp = 1491.399744 W
  let result = hp_to_watt(2.0)
  assert_eq((result - 1491.399744).abs() < 0.01, true)
}

///|
test "watt_to_hp - basic calculation" {
  // 745.699872 W = 1 hp
  let result = watt_to_hp(745.699872)
  assert_eq((result - 1.0).abs() < 0.0001, true)
}

///|
test "btu_per_h_to_watt - basic calculation" {
  // 1 BTU/h = 0.29307107 W
  // 100 BTU/h = 29.307107 W
  let result = btu_per_h_to_watt(100.0)
  assert_eq((result - 29.307107).abs() < 0.0001, true)
}

///|
test "watt_to_btu_per_h - basic calculation" {
  // 0.29307107 W = 1 BTU/h
  let result = watt_to_btu_per_h(0.29307107)
  assert_eq((result - 1.0).abs() < 0.0001, true)
}

// ============================================================================
// 48. 单位转换测试：温度（Unit Conversions: Temperature）
// ============================================================================

///|
test "fahrenheit_to_kelvin - basic calculation" {
  // 32°F = 273.15 K (水的冰点)
  let result = fahrenheit_to_kelvin(32.0)
  assert_eq((result - 273.15).abs() < 0.01, true)
}

///|
test "kelvin_to_fahrenheit - basic calculation" {
  // 273.15 K = 32°F
  let result = kelvin_to_fahrenheit(273.15)
  assert_eq((result - 32.0).abs() < 0.01, true)
}

///|
test "fahrenheit_to_kelvin - boiling point" {
  // 212°F = 373.15 K (水的沸点)
  let result = fahrenheit_to_kelvin(212.0)
  assert_eq((result - 373.15).abs() < 0.01, true)
}

// ============================================================================
// 49. 单位转换测试：压力（Unit Conversions: Pressure）
// ============================================================================

///|
test "atm_to_pascal - basic calculation" {
  // 1 atm = 101325 Pa
  // 2 atm = 202650 Pa
  let result = atm_to_pascal(2.0)
  assert_eq(result, 202650.0)
}

///|
test "pascal_to_atm - basic calculation" {
  // 101325 Pa = 1 atm
  let result = pascal_to_atm(101325.0)
  assert_eq((result - 1.0).abs() < 0.0001, true)
}

///|
test "bar_to_pascal - basic calculation" {
  // 1 bar = 100000 Pa
  // 2 bar = 200000 Pa
  let result = bar_to_pascal(2.0)
  assert_eq(result, 200000.0)
}

///|
test "pascal_to_bar - basic calculation" {
  // 100000 Pa = 1 bar
  let result = pascal_to_bar(100000.0)
  assert_eq((result - 1.0).abs() < 0.0001, true)
}

///|
test "psi_to_pascal - basic calculation" {
  // 1 psi = 6894.757293168 Pa
  // 10 psi = 68947.57293168 Pa
  let result = psi_to_pascal(10.0)
  assert_eq((result - 68947.57293168).abs() < 0.01, true)
}

///|
test "pascal_to_psi - basic calculation" {
  // 6894.757293168 Pa = 1 psi
  let result = pascal_to_psi(6894.757293168)
  assert_eq((result - 1.0).abs() < 0.0001, true)
}

// ============================================================================
// 50. 坐标系转换测试：2D极坐标（Coordinate Conversions: 2D Polar）
// ============================================================================

///|
test "polar_to_cartesian_2d - basic calculation" {
  // r=5, theta=0 (沿x轴)
  // x=5*cos(0)=5, y=5*sin(0)=0
  let (x, y) = polar_to_cartesian_2d(5.0, 0.0)
  assert_eq((x - 5.0).abs() < 0.0001, true)
  assert_eq((y - 0.0).abs() < 0.0001, true)
}

///|
test "polar_to_cartesian_2d - 45 degrees" {
  // r=10, theta=π/4
  // x=10*cos(π/4)≈7.071, y=10*sin(π/4)≈7.071
  let pi = 3.141592653589793
  let (x, y) = polar_to_cartesian_2d(10.0, pi / 4.0)
  assert_eq((x - 7.071).abs() < 0.1, true)
  assert_eq((y - 7.071).abs() < 0.1, true)
}

///|
test "cartesian_to_polar_2d - basic calculation" {
  // x=3, y=4
  // r=5, theta≈0.9273 rad (约53.13°)
  let (r, theta) = cartesian_to_polar_2d(3.0, 4.0)
  assert_eq((r - 5.0).abs() < 0.0001, true)
  assert_eq((theta - 0.9273).abs() < 0.01, true)
}

///|
test "cartesian_to_polar_2d - round trip" {
  // 测试往返转换：极坐标 -> 笛卡尔 -> 极坐标
  let r0 = 10.0
  let theta0 = 3.141592653589793 / 3.0 // 60°
  let (x, y) = polar_to_cartesian_2d(r0, theta0)
  let (r1, theta1) = cartesian_to_polar_2d(x, y)
  assert_eq((r1 - r0).abs() < 0.0001, true)
  // theta可能相差2π的倍数，需要归一化比较
  let pi = 3.141592653589793
  let theta_diff = Double::abs(theta1 - theta0)
  assert_eq(theta_diff < 0.01 || (theta_diff - 2.0 * pi).abs() < 0.01, true)
}

// ============================================================================
// 51. 坐标系转换测试：3D球坐标（Coordinate Conversions: 3D Spherical）
// ============================================================================

///|
test "spherical_to_cartesian_3d - basic calculation" {
  // r=5, theta=0 (沿z轴), phi=0
  // x=0, y=0, z=5
  let (x, y, z) = spherical_to_cartesian_3d(5.0, 0.0, 0.0)
  assert_eq((x - 0.0).abs() < 0.0001, true)
  assert_eq((y - 0.0).abs() < 0.0001, true)
  assert_eq((z - 5.0).abs() < 0.0001, true)
}

///|
test "spherical_to_cartesian_3d - xy plane" {
  // r=10, theta=π/2 (xy平面), phi=π/4
  // x=10*sin(π/2)*cos(π/4)≈7.071, y=10*sin(π/2)*sin(π/4)≈7.071, z=0
  let pi = 3.141592653589793
  let (x, y, z) = spherical_to_cartesian_3d(10.0, pi / 2.0, pi / 4.0)
  assert_eq((x - 7.071).abs() < 0.1, true)
  assert_eq((y - 7.071).abs() < 0.1, true)
  assert_eq((z - 0.0).abs() < 0.0001, true)
}

///|
test "cartesian_to_spherical_3d - basic calculation" {
  // x=0, y=0, z=5
  // r=5, theta=0, phi=0
  let (r, theta, phi) = cartesian_to_spherical_3d(0.0, 0.0, 5.0)
  assert_eq((r - 5.0).abs() < 0.0001, true)
  assert_eq((theta - 0.0).abs() < 0.0001, true)
  assert_eq((phi - 0.0).abs() < 0.0001, true)
}

///|
test "cartesian_to_spherical_3d - round trip" {
  // 测试往返转换：球坐标 -> 笛卡尔 -> 球坐标
  let r0 = 10.0
  let theta0 = 3.141592653589793 / 3.0 // 60°
  let phi0 = 3.141592653589793 / 4.0 // 45°
  let (x, y, z) = spherical_to_cartesian_3d(r0, theta0, phi0)
  let (r1, theta1, phi1) = cartesian_to_spherical_3d(x, y, z)
  assert_eq((r1 - r0).abs() < 0.0001, true)
  assert_eq((theta1 - theta0).abs() < 0.01, true)
  // phi可能相差2π的倍数
  let pi = 3.141592653589793
  let phi_diff = Double::abs(phi1 - phi0)
  assert_eq(phi_diff < 0.01 || (phi_diff - 2.0 * pi).abs() < 0.01, true)
}

// ============================================================================
// 52. 坐标系转换测试：3D柱坐标（Coordinate Conversions: 3D Cylindrical）
// ============================================================================

///|
test "cylindrical_to_cartesian_3d - basic calculation" {
  // rho=5, phi=0, z=3
  // x=5*cos(0)=5, y=5*sin(0)=0, z=3
  let (x, y, z) = cylindrical_to_cartesian_3d(5.0, 0.0, 3.0)
  assert_eq((x - 5.0).abs() < 0.0001, true)
  assert_eq((y - 0.0).abs() < 0.0001, true)
  assert_eq((z - 3.0).abs() < 0.0001, true)
}

///|
test "cylindrical_to_cartesian_3d - 45 degrees" {
  // rho=10, phi=π/4, z=2
  // x=10*cos(π/4)≈7.071, y=10*sin(π/4)≈7.071, z=2
  let pi = 3.141592653589793
  let (x, y, z) = cylindrical_to_cartesian_3d(10.0, pi / 4.0, 2.0)
  assert_eq((x - 7.071).abs() < 0.1, true)
  assert_eq((y - 7.071).abs() < 0.1, true)
  assert_eq((z - 2.0).abs() < 0.0001, true)
}

///|
test "cartesian_to_cylindrical_3d - basic calculation" {
  // x=3, y=4, z=5
  // rho=5, phi≈0.9273, z=5
  let (rho, phi, z) = cartesian_to_cylindrical_3d(3.0, 4.0, 5.0)
  assert_eq((rho - 5.0).abs() < 0.0001, true)
  assert_eq((phi - 0.9273).abs() < 0.01, true)
  assert_eq((z - 5.0).abs() < 0.0001, true)
}

///|
test "cartesian_to_cylindrical_3d - round trip" {
  // 测试往返转换：柱坐标 -> 笛卡尔 -> 柱坐标
  let rho0 = 10.0
  let phi0 = 3.141592653589793 / 3.0 // 60°
  let z0 = 5.0
  let (x, y, z) = cylindrical_to_cartesian_3d(rho0, phi0, z0)
  let (rho1, phi1, z1) = cartesian_to_cylindrical_3d(x, y, z)
  assert_eq((rho1 - rho0).abs() < 0.0001, true)
  assert_eq((z1 - z0).abs() < 0.0001, true)
  // phi可能相差2π的倍数
  let pi = 3.141592653589793
  let phi_diff = Double::abs(phi1 - phi0)
  assert_eq(phi_diff < 0.01 || (phi_diff - 2.0 * pi).abs() < 0.01, true)
}

// ============================================================================
// 53. 量纲分析测试（Dimensional Analysis）
// ============================================================================

///|
test "dimensional_check_velocity - basic calculation" {
  // 长度 100 m，时间 10 s
  // 速度 = 100/10 = 10 m/s
  let result = dimensional_check_velocity(100.0, 10.0)
  assert_eq(result, 10.0)
}

///|
test "dimensional_check_velocity - zero time" {
  // 时间=0，应返回0
  let result = dimensional_check_velocity(100.0, 0.0)
  assert_eq(result, 0.0)
}

///|
test "dimensional_check_acceleration - basic calculation" {
  // 长度 50 m，时间 5 s
  // 加速度 = 50/(5²) = 2 m/s²
  let result = dimensional_check_acceleration(50.0, 5.0)
  assert_eq(result, 2.0)
}

///|
test "dimensional_check_force - basic calculation" {
  // 质量 5 kg，加速度 4 m/s²
  // 力 = 5*4 = 20 N
  let result = dimensional_check_force(5.0, 4.0)
  assert_eq(result, 20.0)
}

///|
test "dimensional_check_force - zero mass" {
  // 质量=0，应返回0
  let result = dimensional_check_force(0.0, 4.0)
  assert_eq(result, 0.0)
}

///|
test "dimensional_check_energy - basic calculation" {
  // 力 10 N，长度 5 m
  // 能量 = 10*5 = 50 J
  let result = dimensional_check_energy(10.0, 5.0)
  assert_eq(result, 50.0)
}

///|
test "dimensional_check_power - basic calculation" {
  // 能量 100 J，时间 5 s
  // 功率 = 100/5 = 20 W
  let result = dimensional_check_power(100.0, 5.0)
  assert_eq(result, 20.0)
}

///|
test "dimensional_check_power - zero time" {
  // 时间=0，应返回0
  let result = dimensional_check_power(100.0, 0.0)
  assert_eq(result, 0.0)
}

///|
test "dimensional_check_pressure - basic calculation" {
  // 力 1000 N，面积 0.01 m²
  // 压力 = 1000/0.01 = 100000 Pa
  let result = dimensional_check_pressure(1000.0, 0.01)
  assert_eq(result, 100000.0)
}

///|
test "dimensional_check_pressure - zero area" {
  // 面积=0，应返回0
  let result = dimensional_check_pressure(1000.0, 0.0)
  assert_eq(result, 0.0)
}

// ============================================================================
// 54. 物理量验证测试（Physical Quantity Validation）
// ============================================================================

///|
test "validate_physical_range - within range" {
  // 值 50，范围 [0, 100]
  let result = validate_physical_range(50.0, 0.0, 100.0)
  assert_eq(result, true)
}

///|
test "validate_physical_range - at minimum" {
  // 值 0，范围 [0, 100]
  let result = validate_physical_range(0.0, 0.0, 100.0)
  assert_eq(result, true)
}

///|
test "validate_physical_range - at maximum" {
  // 值 100，范围 [0, 100]
  let result = validate_physical_range(100.0, 0.0, 100.0)
  assert_eq(result, true)
}

///|
test "validate_physical_range - below minimum" {
  // 值 -10，范围 [0, 100]
  let result = validate_physical_range(-10.0, 0.0, 100.0)
  assert_eq(result, false)
}

///|
test "validate_physical_range - above maximum" {
  // 值 150，范围 [0, 100]
  let result = validate_physical_range(150.0, 0.0, 100.0)
  assert_eq(result, false)
}

///|
test "is_non_relativistic - low velocity" {
  // 速度 3e6 m/s，光速 3e8 m/s，阈值 0.1
  // v/c = 0.01 < 0.1，非相对论性
  let result = is_non_relativistic(3.0e6, 3.0e8, 0.1)
  assert_eq(result, true)
}

///|
test "is_non_relativistic - high velocity" {
  // 速度 1e8 m/s，光速 3e8 m/s，阈值 0.1
  // v/c = 0.333 > 0.1，相对论性
  let result = is_non_relativistic(1.0e8, 3.0e8, 0.1)
  assert_eq(result, false)
}

///|
test "is_non_relativistic - zero speed of light" {
  // 光速=0，应返回false
  let result = is_non_relativistic(100.0, 0.0, 0.1)
  assert_eq(result, false)
}

///|
test "is_non_relativistic - negative velocity" {
  // 负速度，绝对值应被使用
  // -3e6 m/s，光速 3e8 m/s，阈值 0.1
  // |v|/c = 0.01 < 0.1，非相对论性
  let result = is_non_relativistic(-3.0e6, 3.0e8, 0.1)
  assert_eq(result, true)
}
