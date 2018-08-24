//! These are different Monte Carlo algorithms.

#[macro_use]
pub mod plugin;
pub mod energy;

use ::system::*;
use clapme::ClapMe;

use serde_yaml;
use serde;

use atomicfile::AtomicFile;

#[derive(ClapMe)]
enum Params<MP, SP> {
    ResumeFrom(::std::path::PathBuf),
    _Params {
        _sys: SP,
        _mc: MP,
        save_as: Option<::std::path::PathBuf>,
    },
}

/// A Monte Carlo algorithm.
pub trait MonteCarlo: Sized + serde::Serialize + ::serde::de::DeserializeOwned {
    /// A type defining a new Monte Carlo.
    type Params: ClapMe;
    /// A type defining the corresponding system
    type System: System;
    /// Create this MonteCarlo from its parameters
    fn from_params(params: Self::Params, system: Self::System, save_as: ::std::path::PathBuf) -> Self;

    /// This method is called when the program is resumed, so certain
    /// flags may be updated, e.g. the maximum number of iterations to
    /// extend a simulation.
    fn update_from_params(&mut self, _params: Self::Params) {}

    /// Create a new simulation from command-line flags.
    fn from_args<S: ClapMe + Into<Self::System>>() -> Self {
        match <Params<Self::Params, S>>::from_args() {
            Params::_Params { _sys, _mc, save_as } => {
                if let Some(ref save_as) = save_as {
                    if let Ok(f) = ::std::fs::File::open(save_as) {
                        if let Ok(mut s) = serde_yaml::from_reader::<_,Self>(&f) {
                            println!("Resuming from file {:?}", save_as);
                            s.update_from_params(_mc);
                            return s;
                        }
                    }
                }
                let save_as = save_as.unwrap_or(::std::path::PathBuf::from("resume.yaml"));
                Self::from_params(_mc, _sys.into(), save_as)
            },
            Params::ResumeFrom(p) => {
                let f = ::std::fs::File::open(&p)
                    .expect(&format!("error reading file {:?}", &p));
                serde_yaml::from_reader(&f).expect("error reading checkpoint?!")
            },
        }
    }

    /// Create a simulation checkpoint.
    fn checkpoint(&self) {
        let f = AtomicFile::create(self.save_as())
            .expect(&format!("error creating file {:?}", self.save_as()));
        serde_yaml::to_writer(&f, self).expect("error writing checkpoint?!")
    }


    /// Make one random move, collecting appropriate statistics.
    fn move_once(&mut self);

    /// Return the system!
    fn system(&self) -> &Self::System;

    /// The number of moves that have been made.
    fn num_moves(&self) -> u64;

    /// The number of accepted moves.
    fn num_accepted_moves(&self) -> u64;

    /// The path to save the file at.
    fn save_as(&self) -> ::std::path::PathBuf;
}
