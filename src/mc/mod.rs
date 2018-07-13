//! These are different Monte Carlo algorithms.

pub mod sad;

use ::system::*;

/// A Monte Carlo algorithm.
pub trait MonteCarlo {
    /// Make one random move, collecting appropriate statistics.
    /// Return the energy of the system after the move.
    fn move_once(&mut self) -> Energy;

    /// The number of moves that have been made.
    fn num_moves(&self) -> u64;

    /// The number of rejected moves.
    fn num_rejected_moves(&self) -> u64;
}
