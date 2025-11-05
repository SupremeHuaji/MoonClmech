/// 经典力学物理常数模块
/// 提供经典力学计算中常用的物理常量

///|
/// 基础物理常量
/// 光速 (m/s)
pub let speed_of_light : Double = 299792458.0

///|
/// 万有引力常数 G (m³/(kg·s²))
pub let gravitational_constant : Double = 6.67430e-11

///|
/// 普朗克常数 h (J·s)
pub let planck_constant : Double = 6.62607015e-34

///|
/// 约化普朗克常数 ħ = h/(2π) (J·s)
pub let reduced_planck_constant : Double = 1.054571817e-34

///|
/// 地球相关常数
/// 地球表面重力加速度 g (m/s²)
pub let earth_gravity : Double = 9.80665

///|
/// 地球平均半径 (m)
pub let earth_radius : Double = 6371000.0

///|
/// 地球质量 (kg)
pub let earth_mass : Double = 5.9722e24

///|
/// 地球平均角速度 (rad/s)
pub let earth_angular_velocity : Double = 7.292115e-5

///|
/// 太阳质量 (kg)
pub let sun_mass : Double = 1.9891e30

///|
/// 月球质量 (kg)
pub let moon_mass : Double = 7.342e22

///|
/// 常用别名 - 重力加速度
pub let g : Double = earth_gravity

///|
/// 常用别名 - 万有引力常数
pub const G : Double = 6.67430e-11

///|
/// 常用别名 - 光速
pub let c : Double = speed_of_light
