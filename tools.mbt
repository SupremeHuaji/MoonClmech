/// 工具函数模块
/// 提供向量运算、单位转换、数值计算等辅助函数

// ============================================================================
// 向量类型定义
// ============================================================================

///|
/// 三维向量结构体
pub struct Vec3 {
  x : Double
  y : Double
  z : Double
}

///|
/// 创建三维向量
pub fn vec3(x : Double, y : Double, z : Double) -> Vec3 {
  Vec3::{ x, y, z }
}

///|
/// 零向量
pub fn vec3_zero() -> Vec3 {
  Vec3::{ x: 0.0, y: 0.0, z: 0.0 }
}

///|
/// 单位向量（x方向）
pub fn vec3_unit_x() -> Vec3 {
  Vec3::{ x: 1.0, y: 0.0, z: 0.0 }
}

///|
/// 单位向量（y方向）
pub fn vec3_unit_y() -> Vec3 {
  Vec3::{ x: 0.0, y: 1.0, z: 0.0 }
}

///|
/// 单位向量（z方向）
pub fn vec3_unit_z() -> Vec3 {
  Vec3::{ x: 0.0, y: 0.0, z: 1.0 }
}

// ============================================================================
// 向量运算
// ============================================================================

///|
/// 向量加法
pub fn vec3_add(a : Vec3, b : Vec3) -> Vec3 {
  Vec3::{ x: a.x + b.x, y: a.y + b.y, z: a.z + b.z }
}

///|
/// 向量减法
pub fn vec3_sub(a : Vec3, b : Vec3) -> Vec3 {
  Vec3::{ x: a.x - b.x, y: a.y - b.y, z: a.z - b.z }
}

///|
/// 向量数乘
pub fn vec3_scale(v : Vec3, s : Double) -> Vec3 {
  Vec3::{ x: v.x * s, y: v.y * s, z: v.z * s }
}

///|
/// 向量点积（内积）
pub fn vec3_dot(a : Vec3, b : Vec3) -> Double {
  a.x * b.x + a.y * b.y + a.z * b.z
}

///|
/// 向量叉积（外积）
pub fn vec3_cross(a : Vec3, b : Vec3) -> Vec3 {
  Vec3::{
    x: a.y * b.z - a.z * b.y,
    y: a.z * b.x - a.x * b.z,
    z: a.x * b.y - a.y * b.x,
  }
}

///|
/// 向量模长（长度）
pub fn vec3_magnitude(v : Vec3) -> Double {
  Double::sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
}

///|
/// 向量模长的平方（避免开方运算）
pub fn vec3_magnitude_sq(v : Vec3) -> Double {
  v.x * v.x + v.y * v.y + v.z * v.z
}

///|
/// 向量归一化（单位化）
pub fn vec3_normalize(v : Vec3) -> Vec3 {
  let mag = vec3_magnitude(v)
  if mag > 0.0 {
    Vec3::{ x: v.x / mag, y: v.y / mag, z: v.z / mag }
  } else {
    vec3_zero()
  }
}

///|
/// 向量夹角（弧度）
pub fn vec3_angle(a : Vec3, b : Vec3) -> Double {
  let dot = vec3_dot(a, b)
  let mag_a = vec3_magnitude(a)
  let mag_b = vec3_magnitude(b)
  if mag_a > 0.0 && mag_b > 0.0 {
    let cos_angle = dot / (mag_a * mag_b)
    // 限制在 [-1, 1] 范围内以避免数值误差
    let clamped = if cos_angle > 1.0 {
      1.0
    } else if cos_angle < -1.0 {
      -1.0
    } else {
      cos_angle
    }
    @math.acos(clamped)
  } else {
    0.0
  }
}

// ============================================================================
// 单位转换
// ============================================================================

///|
/// 角度转弧度
pub fn deg_to_rad(degrees : Double) -> Double {
  degrees * 3.141592653589793 / 180.0
}

///|
/// 弧度转角度
pub fn rad_to_deg(radians : Double) -> Double {
  radians * 180.0 / 3.141592653589793
}

///|
/// 米转千米
pub fn m_to_km(meters : Double) -> Double {
  meters / 1000.0
}

///|
/// 千米转米
pub fn km_to_m(kilometers : Double) -> Double {
  kilometers * 1000.0
}

///|
/// 千克转克
pub fn kg_to_g(kilograms : Double) -> Double {
  kilograms * 1000.0
}

///|
/// 克转千克
pub fn g_to_kg(grams : Double) -> Double {
  grams / 1000.0
}

///|
/// 秒转小时
pub fn s_to_h(seconds : Double) -> Double {
  seconds / 3600.0
}

///|
/// 小时转秒
pub fn h_to_s(hours : Double) -> Double {
  hours * 3600.0
}

// ============================================================================
// 数值计算工具
// ============================================================================

///|
/// 符号函数
pub fn sign(x : Double) -> Int {
  if x > 0.0 {
    1
  } else if x < 0.0 {
    -1
  } else {
    0
  }
}

///|
/// 将值限制在 [min, max] 范围内
pub fn clamp(value : Double, min_val : Double, max_val : Double) -> Double {
  if value < min_val {
    min_val
  } else if value > max_val {
    max_val
  } else {
    value
  }
}

///|
/// 线性插值
pub fn lerp(a : Double, b : Double, t : Double) -> Double {
  a + (b - a) * t
}

///|
/// 判断两个浮点数是否近似相等（考虑浮点误差）
pub fn approx_eq(a : Double, b : Double, epsilon : Double) -> Bool {
  Double::abs(a - b) < epsilon
}

///|
/// 默认精度下的近似相等判断（epsilon = 1e-9）
pub fn approx_equal(a : Double, b : Double) -> Bool {
  approx_eq(a, b, 1.0e-9)
}
