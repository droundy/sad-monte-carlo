//! Systems are things that have energy and can be changed into
//! different configurations.

use super::rng::MyRng;

pub mod units;
pub mod ising;
pub mod square;
pub mod optsquare;

/// A unitless number.  This is equivalent to a f64, but makes clear
/// that it is going to be interpreted as a dimensionless quantity.
pub type Unitless = units::Unitless<f64>;

/// An energy
pub type Energy = units::Energy<f64>;

/// A distance
pub type Length = units::Distance<f64>;

/// An area
pub type Area = units::Area<f64>;

/// A volume
pub type Volume = units::Volume<f64>;

/// A force
pub type Force = units::Force<f64>;

/// A physical system, which has some energy, and to which we can make
/// some changes.
pub trait System : ::serde::Serialize + ::serde::de::DeserializeOwned {
    /// Returns the energy of the system, and is fast.  This should
    /// just access a cached variable.
    fn energy(&self) -> Energy;
    /// Returns the energy of the system the hard way.  This is slow,
    /// and should only be used as a test that the changes of energy
    /// are being tracked properly.
    fn compute_energy(&self) -> Energy;
    /// The "native" energy step size, assuming one exists.  The
    /// default implementation says the energy is continuous, i.e. it
    /// returns a `None` value.
    fn delta_energy(&self) -> Option<Energy> {
        None
    }
    /// The greatest possible energy, if such a thing exists and can
    /// be estimated.
    fn greatest_possible_energy(&self) -> Option<Energy> {
        None
    }
    /// The lowest possible energy, if such a thing exists and can be
    /// estimated.
    fn lowest_possible_energy(&self) -> Option<Energy> {
        None
    }
}


/// A system that can be reverted after a change is made.
pub trait UndoSystem : System {
    /// Revert the change, whatever it may have been.
    fn undo(&mut self);
}

/// A system that can be moved.
pub trait MovableSystem : UndoSystem {
    /// Moves an atom, and returns the change in energy of the system.
    /// If, however, the move fails (i.e. has infinite energy), `None`
    /// is returned, and no change is made to the system.
    fn move_once(&mut self, &mut MyRng, mean_distance: Length) -> Option<Energy>;
}

/// A system that can gain or lose atoms?
pub trait GrandSystem : UndoSystem {
    /// Adds an atom, and returns the change in energy of the system.
    /// If, however, the move fails (i.e. has infinite energy), `None`
    /// is returned, and no change is made to the system.
    fn add_atom(&mut self, &mut MyRng) -> Option<Energy>;
    /// Removes an atom, and returns the change in energy of the system.
    fn remove_atom(&mut self, &mut MyRng) -> Energy;
}
