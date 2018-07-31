//! These are different Monte Carlo algorithms.

pub mod sad;

use ::system::*;
use clapme::ClapMe;

use serde_yaml;
use serde;

#[derive(ClapMe)]
enum Params<MP, SP> {
    _ResumeFrom(::std::path::PathBuf),
    _Params {
        _sys: SP,
        _mc: MP,
        save_as: Option<::std::path::PathBuf>,
    },
}

/// A Monte Carlo algorithm.
pub trait MonteCarlo: Sized + serde::Serialize {
    /// A type defining a new Monte Carlo.
    type Params: ClapMe;
    /// A type defining the corresponding system
    type System: System;
    /// Create this MonteCarlo from its parameters
    fn from_params(params: Self::Params, system: Self::System, save_as: ::std::path::PathBuf) -> Self;

    /// Create a new simulation from command-line flags.
    fn from_args() -> Self {
        match <Params<Self::Params, <Self::System as System>::Params>>::from_args() {
            Params::_Params { _sys, _mc, save_as } => {
                let s = Self::System::from_params(_sys);
                let save_as = save_as.unwrap_or(::std::path::PathBuf::from("resume.yaml"));
                Self::from_params(_mc, s, save_as)
            },
            Params::ResumeFrom(p) => {
                panic!("unable to open file {:?}", p)
            },
        }
    }

    /// Create a simulation checkpoint.
    fn checkpoint(&self) {
        let f = ::std::fs::File::create(self.save_as())
            .expect(&format!("error creating file {:?}", self.save_as()));
        serde_yaml::to_writer(&f, self).expect("error writing checkpoint?!")
    }


    /// Make one random move, collecting appropriate statistics.
    /// Return the energy of the system after the move.
    fn move_once(&mut self) -> Energy;

    /// The number of moves that have been made.
    fn num_moves(&self) -> u64;

    /// The number of rejected moves.
    fn num_rejected_moves(&self) -> u64;

    /// The path to save the file at.
    fn save_as(&self) -> ::std::path::PathBuf;
}
