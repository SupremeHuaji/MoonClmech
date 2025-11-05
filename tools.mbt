/// Utility functions module
/// Provides auxiliary functions for vector operations, unit conversions, numerical calculations, etc.

///|
/// 3D vector structure
pub struct Vec3 {
  x : Double
  y : Double
  z : Double
}

///|
/// Create a 3D vector
pub fn vec3(x : Double, y : Double, z : Double) -> Vec3 {
  Vec3::{ x, y, z }
}

///|
/// Zero vector
pub fn vec3_zero() -> Vec3 {
  Vec3::{ x: 0.0, y: 0.0, z: 0.0 }
}

///|
/// Unit vector (x direction)
pub fn vec3_unit_x() -> Vec3 {
  Vec3::{ x: 1.0, y: 0.0, z: 0.0 }
}

///|
/// Unit vector (y direction)
pub fn vec3_unit_y() -> Vec3 {
  Vec3::{ x: 0.0, y: 1.0, z: 0.0 }
}

///|
/// Unit vector (z direction)
pub fn vec3_unit_z() -> Vec3 {
  Vec3::{ x: 0.0, y: 0.0, z: 1.0 }
}

///|
/// Vector addition
pub fn vec3_add(a : Vec3, b : Vec3) -> Vec3 {
  Vec3::{ x: a.x + b.x, y: a.y + b.y, z: a.z + b.z }
}

///|
/// Vector subtraction
pub fn vec3_sub(a : Vec3, b : Vec3) -> Vec3 {
  Vec3::{ x: a.x - b.x, y: a.y - b.y, z: a.z - b.z }
}

///|
/// Vector scalar multiplication
pub fn vec3_scale(v : Vec3, s : Double) -> Vec3 {
  Vec3::{ x: v.x * s, y: v.y * s, z: v.z * s }
}

///|
/// Vector dot product (inner product)
pub fn vec3_dot(a : Vec3, b : Vec3) -> Double {
  a.x * b.x + a.y * b.y + a.z * b.z
}

///|
/// Vector cross product (outer product)
pub fn vec3_cross(a : Vec3, b : Vec3) -> Vec3 {
  Vec3::{
    x: a.y * b.z - a.z * b.y,
    y: a.z * b.x - a.x * b.z,
    z: a.x * b.y - a.y * b.x,
  }
}

///|
/// Vector magnitude (length)
pub fn vec3_magnitude(v : Vec3) -> Double {
  Double::sqrt(v.x * v.x + v.y * v.y + v.z * v.z)
}

///|
/// Square of vector magnitude (avoids square root operation)
pub fn vec3_magnitude_sq(v : Vec3) -> Double {
  v.x * v.x + v.y * v.y + v.z * v.z
}

///|
/// Vector normalization (unit vector)
pub fn vec3_normalize(v : Vec3) -> Vec3 {
  let mag = vec3_magnitude(v)
  if mag > 0.0 {
    Vec3::{ x: v.x / mag, y: v.y / mag, z: v.z / mag }
  } else {
    vec3_zero()
  }
}

///|
/// Angle between vectors (radians)
pub fn vec3_angle(a : Vec3, b : Vec3) -> Double {
  let dot = vec3_dot(a, b)
  let mag_a = vec3_magnitude(a)
  let mag_b = vec3_magnitude(b)
  if mag_a > 0.0 && mag_b > 0.0 {
    let cos_angle = dot / (mag_a * mag_b)
    // Clamp to [-1, 1] range to avoid numerical errors
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

///|
/// Convert degrees to radians
pub fn deg_to_rad(degrees : Double) -> Double {
  degrees * 3.141592653589793 / 180.0
}

///|
/// Convert radians to degrees
pub fn rad_to_deg(radians : Double) -> Double {
  radians * 180.0 / 3.141592653589793
}

///|
/// Convert meters to kilometers
pub fn m_to_km(meters : Double) -> Double {
  meters / 1000.0
}

///|
/// Convert kilometers to meters
pub fn km_to_m(kilometers : Double) -> Double {
  kilometers * 1000.0
}

///|
/// Convert kilograms to grams
pub fn kg_to_g(kilograms : Double) -> Double {
  kilograms * 1000.0
}

///|
/// Convert grams to kilograms
pub fn g_to_kg(grams : Double) -> Double {
  grams / 1000.0
}

///|
/// Convert seconds to hours
pub fn s_to_h(seconds : Double) -> Double {
  seconds / 3600.0
}

///|
/// Convert hours to seconds
pub fn h_to_s(hours : Double) -> Double {
  hours * 3600.0
}

///|
/// Sign function
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
/// Clamp value to [min, max] range
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
/// Linear interpolation
pub fn lerp(a : Double, b : Double, t : Double) -> Double {
  a + (b - a) * t
}

///|
/// Check if two floating-point numbers are approximately equal (considering floating-point errors)
pub fn approx_eq(a : Double, b : Double, epsilon : Double) -> Bool {
  Double::abs(a - b) < epsilon
}

///|
/// Approximate equality check with default precision (epsilon = 1e-9)
pub fn approx_equal(a : Double, b : Double) -> Bool {
  approx_eq(a, b, 1.0e-9)
}

///|
/// Vector projection: project vector a onto vector b
pub fn vec3_project(a : Vec3, b : Vec3) -> Vec3 {
  let mag_b_sq = vec3_magnitude_sq(b)
  if mag_b_sq > 0.0 {
    let scale = vec3_dot(a, b) / mag_b_sq
    vec3_scale(b, scale)
  } else {
    vec3_zero()
  }
}

///|
/// Vector rejection: component of a perpendicular to b
pub fn vec3_reject(a : Vec3, b : Vec3) -> Vec3 {
  vec3_sub(a, vec3_project(a, b))
}

///|
/// Vector reflection: reflect vector a about normal n
pub fn vec3_reflect(a : Vec3, normal : Vec3) -> Vec3 {
  let n = vec3_normalize(normal)
  let dot = vec3_dot(a, n)
  vec3_sub(a, vec3_scale(n, 2.0 * dot))
}

///|
/// Distance between two points (3D)
pub fn vec3_distance(a : Vec3, b : Vec3) -> Double {
  vec3_magnitude(vec3_sub(b, a))
}

///|
/// Square of distance between two points (avoids square root)
pub fn vec3_distance_sq(a : Vec3, b : Vec3) -> Double {
  vec3_magnitude_sq(vec3_sub(b, a))
}

///|
/// Linear interpolation between two vectors
pub fn vec3_lerp(a : Vec3, b : Vec3, t : Double) -> Vec3 {
  Vec3::{ x: lerp(a.x, b.x, t), y: lerp(a.y, b.y, t), z: lerp(a.z, b.z, t) }
}

///|
/// Vector component along direction
pub fn vec3_component(v : Vec3, direction : Vec3) -> Double {
  let dir_norm = vec3_normalize(direction)
  vec3_dot(v, dir_norm)
}

///|
/// Convert meters per second to kilometers per hour
pub fn mps_to_kmph(mps : Double) -> Double {
  mps * 3.6
}

///|
/// Convert kilometers per hour to meters per second
pub fn kmph_to_mps(kmph : Double) -> Double {
  kmph / 3.6
}

///|
/// Convert Newtons to pounds-force
pub fn n_to_lbf(newtons : Double) -> Double {
  newtons / 4.44822
}

///|
/// Convert pounds-force to Newtons
pub fn lbf_to_n(pounds : Double) -> Double {
  pounds * 4.44822
}

///|
/// Convert joules to calories
pub fn j_to_cal(joules : Double) -> Double {
  joules / 4.184
}

///|
/// Convert calories to joules
pub fn cal_to_j(calories : Double) -> Double {
  calories * 4.184
}

///|
/// Convert watts to horsepower
pub fn w_to_hp(watts : Double) -> Double {
  watts / 745.7
}

///|
/// Convert horsepower to watts
pub fn hp_to_w(horsepower : Double) -> Double {
  horsepower * 745.7
}

///|
/// Convert Celsius to Kelvin
pub fn celsius_to_kelvin(celsius : Double) -> Double {
  celsius + 273.15
}

///|
/// Convert Kelvin to Celsius
pub fn kelvin_to_celsius(kelvin : Double) -> Double {
  kelvin - 273.15
}

///|
/// Convert Fahrenheit to Celsius
pub fn fahrenheit_to_celsius(fahrenheit : Double) -> Double {
  (fahrenheit - 32.0) * 5.0 / 9.0
}

///|
/// Convert Celsius to Fahrenheit
pub fn celsius_to_fahrenheit(celsius : Double) -> Double {
  celsius * 9.0 / 5.0 + 32.0
}

///|
/// Maximum value in array
pub fn max_value(values : Array[Double]) -> Double {
  if values.length() == 0 {
    0.0
  } else {
    let mut max_val = values[0]
    for i = 1; i < values.length(); i = i + 1 {
      if values[i] > max_val {
        max_val = values[i]
      }
    }
    max_val
  }
}

///|
/// Minimum value in array
pub fn min_value(values : Array[Double]) -> Double {
  if values.length() == 0 {
    0.0
  } else {
    let mut min_val = values[0]
    for i = 1; i < values.length(); i = i + 1 {
      if values[i] < min_val {
        min_val = values[i]
      }
    }
    min_val
  }
}

///|
/// Absolute maximum value
pub fn abs_max(a : Double, b : Double) -> Double {
  if Double::abs(a) > Double::abs(b) {
    Double::abs(a)
  } else {
    Double::abs(b)
  }
}

///|
/// Absolute minimum value
pub fn abs_min(a : Double, b : Double) -> Double {
  if Double::abs(a) < Double::abs(b) {
    Double::abs(a)
  } else {
    Double::abs(b)
  }
}

///|
/// Smooth step interpolation (smooth transition between 0 and 1)
pub fn smoothstep(edge0 : Double, edge1 : Double, x : Double) -> Double {
  let t = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0)
  t * t * (3.0 - 2.0 * t)
}

///|
/// Check if value is within range [min, max]
pub fn in_range(value : Double, min_val : Double, max_val : Double) -> Bool {
  value >= min_val && value <= max_val
}

///|
/// Map value from one range to another
pub fn map_range(
  value : Double,
  from_min : Double,
  from_max : Double,
  to_min : Double,
  to_max : Double,
) -> Double {
  if from_max != from_min {
    let normalized = (value - from_min) / (from_max - from_min)
    to_min + normalized * (to_max - to_min)
  } else {
    to_min
  }
}
