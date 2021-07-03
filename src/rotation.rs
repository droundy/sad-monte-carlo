//! 3D rotation type.

use crate::unit_quaternion::UnitQuaternion;
use vector3d::Vector3d;

/// 3D rotation.
#[derive(Copy, Clone)]
pub struct Rotation {
    /// The underlying quaternion.
    pub q: UnitQuaternion,
}

impl<T> std::ops::Mul<Vector3d<T>> for Rotation
where
    T: Copy
        + std::ops::Mul<f64, Output = T>
        + std::ops::Add<Output = T>
        + std::ops::Sub<Output = T>, // TODO: maybe change this to Length
    f64: std::ops::Mul<T, Output = T> + std::ops::Mul<f64, Output = f64>,
{
    type Output = Vector3d<T>;
    fn mul(self, v: Vector3d<T>) -> Vector3d<T> {
        // https://gamedev.stackexchange.com/a/50545
        let u = self.q.vector();
        let s = self.q.w;
        u * u.dot(v) * 2.0 + v * (s * s - u.norm2()) + u.cross(v) * s * 2.0
    }
}
