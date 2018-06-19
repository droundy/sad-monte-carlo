//! Systems are things that have energy and can be changed into
//! different configurations.

pub mod units;
pub mod ising;

/// An energy
pub type Energy = units::Energy<f64>;

/// A distance
pub type Length = units::Distance<f64>;

/// A force
pub type Force = units::Force<f64>;

/// A physical system, which has some energy, and to which we can make
/// some changes.
pub trait System {
    /// Returns the energy of the system.
    fn energy(&self) -> Energy;
    /// The number of possible energies?
    fn num_energies(&self) -> Energy;
}

/// A system that can be moved.
pub trait MovableSystem : System {
    /// Moves an atom, and returns the change in energy of the system.
    fn move_once(&mut self, mean_distance: Length) -> Energy;
    /// Revert the move
    fn undo_move(&mut self);
}

/// A system that can gain or lose atoms?
pub trait GrandSystem : System {
    /// Adds an atom, and returns the change in energy of the system.
    fn add_atom(&mut self) -> Energy;
    /// Removes an atom, and returns the change in energy of the system.
    fn remove_atom(&mut self) -> Energy;
    /// Revert the add/remove
    fn undo_grand(&mut self);
}
