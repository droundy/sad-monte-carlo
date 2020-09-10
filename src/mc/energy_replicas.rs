//! An implementation of SAD (Statistical Association, Dynamical version).
//!
//! This version using a [Binning]() approach so that we can use
//! something nicer than a histogram.  The idea is also to use this
//! for a 2D histogram, which will be even fancier.

#![allow(non_snake_case)]

use super::*;
use crate::system::*;
use rayon::prelude::*;

use rand::{Rng, SeedableRng};
use std::default::Default;

use std::collections::HashSet;

/// The parameters needed to configure a simulation.
#[derive(Debug, AutoArgs, Clone)]
pub struct MCParams {
    /// The lowest temperature of interest
    min_T: Energy,
    /// The seed for the random number generator.
    pub seed: Option<u64>,
    /// report input
    pub _movies: plugin::MovieParams,
    /// report input
    pub _save: plugin::SaveParams,
}

impl Default for MCParams {
    fn default() -> Self {
        MCParams {
            min_T: 0.2 * units::EPSILON,
            seed: None,
            _save: plugin::SaveParams::default(),
            _movies: plugin::MovieParams::default(),
        }
    }
}

/// A simple simulation, with a maximum energy, a single system, and an energy
/// divider.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Replica<S> {
    /// A unique integer indicating which replica it is
    pub id: usize,
    /// The maximum allowed energy
    pub max_energy: Energy,
    /// The cutoff energy
    pub cutoff_energy: Energy,
    /// The number of rejected moves
    pub rejected_count: u64,
    /// The counts above
    pub above_count: u64,
    /// The counts below
    pub below_count: u64,
    /// The total energy above (for computing the mean)
    pub above_total: Energy,
    /// The total energy below (for computing the mean)
    pub below_total: Energy,
    /// The system itself
    pub system: S,
    /// The set of regions visited by this system
    pub system_visited: HashSet<usize>,
    /// The number of independent systems seen
    pub independent_count: u64,
    /// The random number generator.
    pub rng: crate::rng::MyRng,

    /// The current translation scale
    pub translation_scale: Length,
}

impl<S: MovableSystem> Replica<S> {
    fn new(
        id: usize,
        max_energy: Energy,
        cutoff_energy: Energy,
        system: S,
        mut system_visited: HashSet<usize>,
        rng: crate::rng::MyRng,
    ) -> Self {
        system_visited.insert(id);
        Replica {
            id,
            max_energy,
            cutoff_energy,
            above_count: 0,
            below_count: 0,
            rejected_count: 0,
            independent_count: 1,
            above_total: Energy::new(0.0),
            below_total: Energy::new(0.0),
            system,
            system_visited,
            rng,

            translation_scale: Length::new(0.05),
        }
    }
    fn new_system(&mut self) {
        if !self.system_visited.contains(&self.id) {
            self.independent_count += 1;
            self.system_visited.insert(self.id);
        }
    }
    fn run_once(&mut self) {
        if let Some(e) = self.system.plan_move(&mut self.rng, Length::new(0.05)) {
            if e < self.max_energy {
                self.system.confirm();
            } else {
                self.rejected_count += 1;
            }
        } else {
            self.rejected_count += 1;
        }
        let e = self.system.energy();
        if e > self.cutoff_energy {
            self.above_count += 1;
            self.above_total += e;
        } else {
            self.below_count += 1;
            self.below_total += e;
        }
    }
    fn occasional_update(&mut self) {
        let moves = self.above_count + self.below_count;
        let rejection_rate = self.rejected_count as f64 / moves as f64;
        if (rejection_rate - 0.5).abs() > 0.1 {
            let old_scale = self.translation_scale;
            self.translation_scale *= 0.5 / rejection_rate;
            println!(
                "    ({:.5}) Updating translation scale from {:.2} -> {:.2} [rejection rate {:.2}",
                self.max_energy.pretty(),
                old_scale.pretty(),
                self.translation_scale.pretty(),
                crate::prettyfloat::PrettyFloat(rejection_rate),
            );
        }
    }
}

/// A simulation with many replicas
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct MC<S> {
    /// The minimum temperature of interest
    min_T: Energy,
    /// The random number generator.
    pub rng: crate::rng::MyRng,
    /// Where to save the resume file.
    pub save_as: ::std::path::PathBuf,
    /// The number of MC moves we have done
    pub moves: u64,

    /// The relative sizes of the bins
    pub replicas: Vec<Replica<S>>,
    /// How frequently to save...
    #[serde(skip, default)]
    save: plugin::Save,
    /// Movie state
    #[serde(skip, default)]
    movie: plugin::Movie,
}

const MAX_INIT: usize = 10; // Looking for 10 initial quantiles requires maybe a dozen megabytes of memory.

impl<
        S: Clone
            + ConfirmSystem
            + MovableSystem
            + Sync
            + Send
            + serde::Serialize
            + serde::de::DeserializeOwned,
    > MC<S>
{
    fn from_params(params: MCParams, mut system: S, save_as: ::std::path::PathBuf) -> Self {
        let mut rng = crate::rng::MyRng::seed_from_u64(params.seed.unwrap_or(0));
        let min_T = params.min_T;

        // First let's brute-force to find where a number of the quantiles are
        // We look for at most MAX_INIT quantiles.
        let mut energies = Vec::with_capacity(1 << (2 * MAX_INIT));
        for _ in 0..1 << (2 * MAX_INIT) {
            energies.push(system.randomize(&mut rng));
        }
        energies.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let replicas = vec![Replica::new(
            0,
            energies[1 << MAX_INIT],
            energies[1 << (MAX_INIT + 1)],
            system,
            HashSet::new(),
            rng.clone(),
        )];
        rng.jump();
        std::mem::drop(energies); // Just to save on RAM...
        MC {
            min_T,
            replicas,
            moves: 0,

            rng,
            save_as: save_as,
            save: plugin::Save::from(params._save),
            movie: plugin::Movie::from(params._movies),
        }
    }

    /// Create a new simulation from command-line flags.
    pub fn from_args<A: AutoArgs + Into<S>>() -> Self {
        println!("git version: {}", VERSION);
        match <Params<MCParams, A>>::from_args() {
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
                        let s = match save_as.extension().and_then(|x| x.to_str()) {
                            Some("yaml") => serde_yaml::from_reader::<_, Self>(&f)
                                .expect("error parsing save-as file"),
                            Some("json") => serde_json::from_reader::<_, Self>(&f)
                                .expect("error parsing save-as file"),
                            Some("cbor") => serde_cbor::from_reader::<Self, _>(&f)
                                .expect("error parsing save-as file"),
                            _ => panic!("I don't know how to read file {:?}", f),
                        };
                        println!("Resuming from file {:?}", save_as);
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
    pub fn checkpoint(&mut self) {
        for r in self.replicas.iter_mut() {
            r.occasional_update();
        }
        let f = AtomicFile::create(&self.save_as)
            .expect(&format!("error creating file {:?}", self.save_as));
        match self.save_as.extension().and_then(|x| x.to_str()) {
            Some("yaml") => serde_yaml::to_writer(&f, self).expect("error writing checkpoint?!"),
            Some("json") => serde_json::to_writer(&f, self).expect("error writing checkpoint?!"),
            Some("cbor") => serde_cbor::to_writer(&f, self).expect("error writing checkpoint?!"),
            _ => panic!("I don't know how to create file {:?}", self.save_as),
        }
        println!("Saved to file {:?}", self.save_as);
    }

    /// Run a simulation
    pub fn run_once(&mut self) {
        // First run a few steps of the simulation for each replica
        self.replicas.par_iter_mut().for_each(|r| {
            r.run_once();
        });

        // Now let us try swapping if we can.
        let iterator = if self.rng.gen::<bool>() {
            // Try swapping odd onces
            self.replicas.chunks_exact_mut(2)
        } else {
            // Try swapping even ones.
            self.replicas[1..].chunks_exact_mut(2)
        };
        iterator.for_each(|chunk| {
            // for chunk in iterator {
            if let [r0, r1] = chunk {
                if r0.system.energy() < r1.max_energy && r1.system.energy() < r0.max_energy {
                    std::mem::swap(&mut r0.system, &mut r1.system);
                    std::mem::swap(&mut r0.system_visited, &mut r1.system_visited);
                    r0.new_system();
                    r1.new_system();
                }
            }
        });

        self.moves += 1;
        let movie_time = self.movie.shall_i_save(self.moves);
        if movie_time {
            self.movie.save_frame(&self.save_as, self.moves, &self);
        }
        if self.save.shall_i_save(self.moves) || movie_time {
            self.checkpoint();
        }
    }
}
