//! These are different Monte Carlo algorithms.

pub mod sad;

use ::system::*;
use clapme::ClapMe;

#[derive(ClapMe)]
enum Params<MP, SP> {
    ResumeFrom(::std::path::PathBuf),
    _Params {
        _sys: SP,
        _mc: MP,
    },
}

/// A Monte Carlo algorithm.
pub trait MonteCarlo: Sized {
    /// A type defining a new Monte Carlo.
    type Params: ClapMe;
    /// A type defining the corresponding system
    type System: System;
    /// Create this MonteCarlo from its parameters
    fn from_params(params: Self::Params, system: Self::System) -> Self;

    /// Create a new simulation from command-line flags.
    fn from_args() -> Self {
        match <Params<Self::Params, <Self::System as System>::Params>>::from_args() {
            Params::_Params { _sys, _mc } => {
                let s = Self::System::from_params(_sys);
                Self::from_params(_mc, s)
            },
            Params::ResumeFrom(p) => {
                panic!("unable to open file {:?}", p)
            },
        }
    }


    /// Make one random move, collecting appropriate statistics.
    /// Return the energy of the system after the move.
    fn move_once(&mut self) -> Energy;

    /// The number of moves that have been made.
    fn num_moves(&self) -> u64;

    /// The number of rejected moves.
    fn num_rejected_moves(&self) -> u64;
}
