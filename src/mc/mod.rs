//! These are different Monte Carlo algorithms.

pub mod binning;
pub mod energy;
pub mod energy_binning;
pub mod energy_number;
pub mod energy_transposed;
pub mod grand_transposed;
pub mod number;
pub mod plugin;

use crate::system::*;
use auto_args::AutoArgs;

use serde;

use crate::atomicfile::AtomicFile;

#[derive(AutoArgs)]
enum Params<MP, SP> {
    ResumeFrom(::std::path::PathBuf),
    _Params {
        _sys: SP,
        _mc: MP,
        save_as: Option<::std::path::PathBuf>,
        /// The maximum number of threads to use (specify 0 for using all cores)
        num_threads: Option<usize>,
    },
}

const VERSION: &str = git_version::git_describe!("--always", "--dirty");

/// A Monte Carlo algorithm.
pub trait MonteCarlo: Sized + serde::Serialize + ::serde::de::DeserializeOwned {
    /// A type defining a new Monte Carlo.
    type Params: AutoArgs;
    /// A type defining the corresponding system
    type System: System;
    /// Create this MonteCarlo from its parameters
    fn from_params(
        params: Self::Params,
        system: Self::System,
        save_as: ::std::path::PathBuf,
    ) -> Self;

    /// This method is called when the program is resumed, so certain
    /// flags may be updated, e.g. the maximum number of iterations to
    /// extend a simulation.
    fn update_from_params(&mut self, _params: Self::Params) {}

    /// Create a new simulation from command-line flags.
    fn from_args<S: AutoArgs + Into<Self::System>>() -> Self {
        println!("git version: {}", VERSION);
        match <Params<Self::Params, S>>::from_args() {
            Params::_Params {
                _sys,
                _mc,
                save_as,
                num_threads,
            } => {
                if let Some(num_threads) = num_threads {
                    rayon::ThreadPoolBuilder::new()
                        .num_threads(num_threads)
                        .build_global()
                        .unwrap()
                }
                if let Some(ref save_as) = save_as {
                    if let Ok(f) = ::std::fs::File::open(save_as) {
                        let mut s = match save_as.extension().and_then(|x| x.to_str()) {
                            Some("yaml") => serde_yaml::from_reader::<_, Self>(&f)
                                .expect("error parsing save-as file"),
                            Some("json") => serde_json::from_reader::<_, Self>(&f)
                                .expect("error parsing save-as file"),
                            Some("cbor") => serde_cbor::from_reader::<Self, _>(&f)
                                .expect("error parsing save-as file"),
                            _ => panic!("I don't know how to read file {:?}", f),
                        };
                        println!("Resuming from file {:?}", save_as);
                        s.update_from_params(_mc);
                        s.system_mut().update_caches();
                        return s;
                    } else {
                        return Self::from_params(_mc, _sys.into(), save_as.clone());
                    }
                }
                let save_as = save_as.unwrap_or(::std::path::PathBuf::from("resume.yaml"));
                Self::from_params(_mc, _sys.into(), save_as)
            }
            Params::ResumeFrom(p) => {
                let f = ::std::fs::File::open(&p).expect(&format!("error reading file {:?}", &p));
                match p.extension().and_then(|x| x.to_str()) {
                    Some("yaml") => {
                        serde_yaml::from_reader(&f).expect("error reading checkpoint?!")
                    }
                    Some("json") => {
                        serde_json::from_reader(&f).expect("error reading checkpoint?!")
                    }
                    Some("cbor") => {
                        serde_cbor::from_reader(&f).expect("error reading checkpoint?!")
                    }
                    _ => panic!("I don't know how to read file {:?}", f),
                }
            }
        }
    }

    /// Create a simulation checkpoint.
    fn checkpoint(&self) {
        let f = AtomicFile::create(self.save_as())
            .expect(&format!("error creating file {:?}", self.save_as()));
        match self.save_as().extension().and_then(|x| x.to_str()) {
            Some("yaml") => serde_yaml::to_writer(&f, self).expect("error writing checkpoint?!"),
            Some("json") => serde_json::to_writer(&f, self).expect("error writing checkpoint?!"),
            Some("cbor") => serde_cbor::to_writer(&f, self).expect("error writing checkpoint?!"),
            _ => panic!("I don't know how to create file {:?}", self.save_as()),
        }
    }

    /// Make one random move, collecting appropriate statistics.
    fn move_once(&mut self);

    /// Return the system!
    fn system(&self) -> &Self::System;

    /// Return the system, but mutable!
    fn system_mut(&mut self) -> &mut Self::System;

    /// The number of moves that have been made.
    fn num_moves(&self) -> u64;

    /// The number of accepted moves.
    fn num_accepted_moves(&self) -> u64;

    /// The path to save the file at.
    fn save_as(&self) -> ::std::path::PathBuf;
}
