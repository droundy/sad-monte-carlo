//! Systems are things that have energy and can be changed into
//! different configurations.

pub mod ising;

/// A physical system, which has some energy, and to which we can make
/// some changes.

pub trait System {
    /// Returns the energy of the system.
    fn energy(&self) -> f64;
}
