//! Unit quaternion type.

use vector3d::Vector3d;

/// A unit quaternion.
#[allow(missing_docs)]
#[derive(Copy, Clone)]
pub struct UnitQuaternion {
    // vector part
    pub x: f64,
    pub y: f64,
    pub z: f64,
    // scalar part
    pub w: f64,
}

impl UnitQuaternion {
    /// The vector part.
    pub fn vector(self) -> Vector3d<f64> {
        Vector3d::new(self.x, self.y, self.z)
    }
}

impl std::ops::Mul for UnitQuaternion {
    type Output = Self;
    fn mul(self, rhs: Self) -> Self {
        let a = self.w;
        let b = self.x;
        let c = self.y;
        let d = self.z;
        let e = rhs.w;
        let f = rhs.x;
        let g = rhs.y;
        let h = rhs.z;
        Self {
            w: a * e - b * f - c * g - d * h,
            x: a * f + b * e + c * h - d * g,
            y: a * g - b * h + c * e + d * f,
            z: a * h + b * g - c * f + d * e,
        }
    }
}
