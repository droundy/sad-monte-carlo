//! Systems are things that have energy and can be changed into
//! different configurations.

use super::rng::MyRng;
use auto_args::AutoArgs;

pub mod cell;
pub mod optcell;
pub mod units;

pub mod ising;
pub mod lattice_gas;
pub mod lj;
pub mod water;
pub mod optsquare;
pub mod square;
pub mod wca;

pub mod any;

pub mod erfinv;
pub mod fake;

pub use crate::mc::binning::Interned;

/// A unitless number.  This is equivalent to a f64, but makes clear
/// that it is going to be interpreted as a dimensionless quantity.
pub type Unitless = units::Unitless<f64>;

/// An energy
pub type Energy = units::Energy<f64>;

/// Per energy
pub type PerEnergy = units::PerEnergy<f64>;

/// An energy squared
pub type EnergySquared = units::EnergySquared<f64>;

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
pub trait System {
    /// Returns the energy of the system, and is fast.  This should
    /// just access a cached variable.
    fn energy(&self) -> Energy;
    /// Completely randomizes the system
    fn randomize(&mut self, _: &mut MyRng) -> Energy;
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
    /// Update any cached info after resume.  This is to make it
    /// possible for a system to avoid saving information that is easy
    /// and safe to recompute when restarting.
    fn update_caches(&mut self) {}
    /// Verify as well as we can that the energy is currently correct.
    fn verify_energy(&self) {}
    /// Collect some data for the current state of the system if we
    /// want to do so.
    fn data_to_collect(&self, _iter: u64) -> Vec<(Interned, f64)> {
        Vec::new()
    }
    /// How many moves at the minimum could change all the coordinates of the
    /// system
    fn min_moves_to_randomize(&self) -> u64;
    /// Dimensions of the configuration space
    fn dimensionality(&self) -> u64;
}

/// A system that can have a change confirmed.
pub trait ConfirmSystem: System {
    /// Confirm the change, whatever it may have been.
    fn confirm(&mut self);
    /// Print something descriptive about this system.
    fn describe(&self) -> String {
        "".to_string()
    }
}

/// A system that can be moved.
pub trait MovableSystem: ConfirmSystem {
    /// Considers moving an atom, and returns the resulting energy of
    /// the system.  If, however, the move is impossible (i.e. has
    /// infinite energy), `None` is returned, and no change is made to
    /// the system.  The atom is not actually moved until the change
    /// is confirmed.
    fn plan_move(&mut self, _: &mut MyRng, mean_distance: Length) -> Option<Energy>;
    /// A maximum reasonable value for mean_distance, i.e. the size of the configuration space.
    fn max_size(&self) -> Length;
}

/// A system that can gain or lose atoms?
pub trait GrandSystem: MovableSystem {
    /// Considers adding an atom, and returns the resulting energy of
    /// the system.  If the add is impossible (i.e. has infinite
    /// energy), `None` is returned, and no change is made to the
    /// system.  The atom is not actually added until the change is
    /// confirmed.
    fn plan_add(&mut self, _: &mut MyRng) -> Option<Energy>;
    /// Considers removing an atom, and returns the change in energy
    /// of the system.  The atom is not actually removed until the
    /// change is confirmed.
    fn plan_remove(&mut self, _: &mut MyRng) -> Energy;
    /// The number of atoms
    fn num_atoms(&self) -> usize;
}

/// A system that can gain or lose atoms?
pub trait GrandReplicaSystem: GrandSystem {
    /// Consider swapping an atom from self to other.  Returns which atom was contemplated to be eremoved,
    /// and the final energies of the two systems.  The first energy is for the system which
    /// lost an atom, and the second is for the system that gained an atom.
    fn plan_swap_atom(&self, other: &Self, _: &mut MyRng) -> Option<(usize, Energy, Energy)>;
    /// Swap an atom from self to the other system.
    fn swap_atom(&mut self, other: &mut Self, which: usize);
}
