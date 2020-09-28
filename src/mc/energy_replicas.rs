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

use std::collections::HashMap;

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
    /// report input
    pub _report: plugin::ReportParams,
}

impl Default for MCParams {
    fn default() -> Self {
        MCParams {
            min_T: 0.2 * units::EPSILON,
            seed: None,
            _save: plugin::SaveParams::default(),
            _movies: plugin::MovieParams::default(),
            _report: plugin::ReportParams::default(),
        }
    }
}

/// An estimator for a median
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct MedianEstimator {
    energies: Vec<Energy>,
    step: usize,
    skip: usize,
}
const ESTIMATOR_SIZE: usize = 32;
impl MedianEstimator {
    fn new(e: Energy) -> Self {
        let mut energies = Vec::with_capacity(ESTIMATOR_SIZE);
        energies.push(e);
        MedianEstimator {
            energies,
            step: 0,
            skip: 1,
        }
    }
    fn reset(&mut self, e: Energy) {
        self.energies.truncate(0);
        self.energies.push(e);
        self.step = 0;
        self.skip = 1;
    }
    fn add_energy(&mut self, e: Energy) {
        if self.step == 0 {
            if self.energies.len() < ESTIMATOR_SIZE {
                self.energies.push(e);
            } else {
                // Keep the middle 1/2 of our energies.
                self.energies.sort_by(|a, b| a.partial_cmp(b).unwrap());
                {
                    // Move the third quarter of the energies back into the
                    // first quarter, to the middle half of energies will be
                    // in the first half of the vec.
                    let (left, right) = self.energies.split_at_mut(ESTIMATOR_SIZE / 2);
                    left[0..ESTIMATOR_SIZE / 4].copy_from_slice(&right[0..ESTIMATOR_SIZE / 4]);
                }
                self.energies.truncate(ESTIMATOR_SIZE / 2);
                self.skip *= 2;
            }
        }
        self.step = (self.step + 1) % self.skip;
    }
    fn median(&mut self) -> Energy {
        self.energies.sort_by(|a, b| a.partial_cmp(b).unwrap());
        if self.energies.len() & 1 == 0 {
            self.energies[self.energies.len() / 2]
        } else {
            0.5 * (self.energies[self.energies.len() / 2]
                + self.energies[self.energies.len() / 2 + 1])
        }
    }
}

/// A simple simulation, with a maximum energy, a single system, and an energy
/// divider.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Replica<S> {
    /// The maximum allowed energy
    pub max_energy: Energy,
    /// The cutoff energy
    pub cutoff_energy: Energy,
    /// The number of rejected moves since we adjusted the translation_scale.
    pub rejected_count: u64,
    /// The number of accepted moves since we adjusted the translation_scale
    pub accepted_count: u64,
    /// The counts above
    pub above_count: u64,
    /// The counts below
    pub below_count: u64,
    /// The total energy above (for computing the mean)
    pub above_total: Energy,
    /// The total energy below (for computing the mean)
    pub below_total: Energy,
    /// Any extra data we might want to collect, total and count
    pub above_extra: HashMap<Interned, (f64, u64)>,
    /// The system itself
    pub system: S,
    /// The lowest `max_energy` that the system has visited
    pub system_lowest_max_energy: Option<Energy>,
    /// The number of completely and provably indepdendent systems that have
    /// us
    pub unique_visitors: u64,
    /// The random number generator.
    pub rng: crate::rng::MyRng,

    /// The current translation scale
    pub translation_scale: Length,
}

impl<S: MovableSystem> Replica<S> {
    fn new(max_energy: Energy, cutoff_energy: Energy, system: S, rng: crate::rng::MyRng) -> Self {
        Replica {
            max_energy,
            cutoff_energy,
            above_count: 0,
            below_count: 0,
            rejected_count: 0,
            accepted_count: 0,
            above_total: Energy::new(0.0),
            below_total: Energy::new(0.0),
            above_extra: HashMap::new(),
            system,
            system_lowest_max_energy: None,
            unique_visitors: 0,
            rng,

            translation_scale: Length::new(0.05),
        }
    }
    fn run_once(&mut self, moves: u64) {
        if self.max_energy.value_unsafe.is_finite() {
            if let Some(e) = self.system.plan_move(&mut self.rng, self.translation_scale) {
                if e < self.max_energy {
                    self.system.confirm();
                    self.accepted_count += 1;
                } else {
                    self.rejected_count += 1;
                }
            } else {
                self.rejected_count += 1;
            }
        } else {
            // If there is no upper bound, we can just entirely randomize the
            // system.
            self.system.randomize(&mut self.rng);
            self.system_lowest_max_energy = Some(self.max_energy);
        }
        let e = self.system.energy();
        if e > self.cutoff_energy {
            self.above_count += 1;
            self.above_total += e;
            for (k, d) in self.system.data_to_collect(moves).into_iter() {
                if let Some(p) = self.above_extra.get_mut(&k) {
                    p.0 += d;
                    p.1 += 1;
                } else {
                    self.above_extra.insert(k, (d, 1));
                }
            }
        } else {
            self.below_count += 1;
            self.below_total += e;
        }
    }
    fn occasional_update(&mut self) {
        if self.rejected_count > 128
            && self.accepted_count > 128
            && self.max_energy.value_unsafe.is_finite()
        {
            let acceptance_ratio = self.accepted_count as f64 / self.rejected_count as f64;
            let max_acceptance_ratio = self.system.min_moves_to_randomize() as f64;
            // We don't mess with the translation scale unless the acceptance
            // is more than a factor of two away from our goal range of between
            // 50% and 1/min_moves_to_randomize.  This gives us a pretty wide
            // target range, and ensures that the translation scale does not
            // too frequently.
            if acceptance_ratio < 0.5 || acceptance_ratio > 2.0 * max_acceptance_ratio {
                let old_translation_scale = self.translation_scale;
                let mut adjustment = if acceptance_ratio < 0.5 {
                    acceptance_ratio / max_acceptance_ratio.sqrt()
                } else {
                    acceptance_ratio * max_acceptance_ratio.sqrt()
                };
                if adjustment > 2.0 {
                    adjustment = 2.0;
                } else if adjustment < 0.5 {
                    adjustment = 0.5;
                }
                self.translation_scale *= adjustment;
                println!(
                    "      ({:.5}) Updating translation scale from {:.2} -> {:.2} [{:.1}]",
                    self.max_energy.pretty(),
                    old_translation_scale.pretty(),
                    self.translation_scale.pretty(),
                    crate::prettyfloat::PrettyFloat(
                        self.accepted_count as f64 / self.rejected_count as f64
                    ),
                );
                // Reset the counts for the next time around only if we have
                // made a change.
                self.accepted_count = 0;
                self.rejected_count = 0;
            }
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

    /// An estimator for the median energy below the lowest bin
    pub median: MedianEstimator,

    /// The relative sizes of the bins
    pub replicas: Vec<Replica<S>>,
    /// How frequently to save...
    save: plugin::Save,
    /// Movie state
    movie: plugin::Movie,
    /// Movie state
    report: plugin::Report,
}

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

        const MAX_INIT: usize = 1 << 15; // Number of energies to use in computing the first couple of quantiles.

        // First let's brute-force to find where a number of the quantiles are
        // We look for at most MAX_INIT quantiles.
        let mut energies = Vec::with_capacity(MAX_INIT);
        for _ in 0..MAX_INIT {
            energies.push(system.randomize(&mut rng));
        }
        let mut high_system = system.clone();
        high_system.randomize(&mut rng);
        while system.energy() > energies[energies.len() / 2] {
            system.randomize(&mut rng);
        }
        energies.sort_by(|a, b| a.partial_cmp(b).unwrap());
        for (i, e) in energies.iter().enumerate() {
            println!("  {:3}: {}", i, e.pretty());
        }
        let replicas = vec![
            Replica::new(
                Energy::new(std::f64::INFINITY),
                energies[energies.len() / 2],
                high_system,
                rng.clone(),
            ),
            Replica::new(
                energies[energies.len() / 2],
                energies[energies.len() / 4],
                system,
                rng.clone(),
            ),
        ];
        rng.jump();
        MC {
            min_T,
            replicas,
            moves: 0,
            median: MedianEstimator::new(energies[energies.len() / 4]),

            rng,
            save_as: save_as,
            save: plugin::Save::from(params._save),
            movie: plugin::Movie::from(params._movies),
            report: plugin::Report::from(params._report),
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
                        for r in s.replicas.iter_mut() {
                            r.system.update_caches();
                        }
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
        self.report.print(self.moves);
        println!("    [{} replicas]", self.replicas.len());

        for r in self.replicas.iter() {
            let percent = r.above_count as f64 / (r.above_count as f64 + r.below_count as f64);
            println!(
                "       < {:9.5} [{:.2}%] {:.2} unique",
                r.cutoff_energy.pretty(),
                crate::prettyfloat::PrettyFloat(100.0 * percent),
                crate::prettyfloat::PrettyFloat(r.unique_visitors as f64),
            );
        }
        if let Some(r) = self.replicas.last() {
            let mean_below = r.below_total / r.below_count as f64;
            println!("        mean_below: {:9.5}", mean_below.pretty());
        };

        let f = AtomicFile::create(&self.save_as)
            .expect(&format!("error creating file {:?}", self.save_as));
        match self.save_as.extension().and_then(|x| x.to_str()) {
            Some("yaml") => serde_yaml::to_writer(&f, self).expect("error writing checkpoint?!"),
            Some("json") => serde_json::to_writer(&f, self).expect("error writing checkpoint?!"),
            Some("cbor") => serde_cbor::to_writer(&f, self).expect("error writing checkpoint?!"),
            _ => panic!("I don't know how to create file {:?}", self.save_as),
        }
    }

    /// Run a simulation
    pub fn run_once(&mut self) {
        let moves = self.moves;
        let steps = self.replicas[0].system.min_moves_to_randomize();
        // First run a few steps of the simulation for each replica
        self.replicas.par_iter_mut().for_each(|r| {
            for i in 0..steps {
                r.run_once(moves + i);
            }
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
                assert_eq!(r0.cutoff_energy, r1.max_energy);
                // We will swap them if both systems can go into the lower bin,
                // the higher energy system is not a clone of some other
                // We want the clones to all end up at the top where they can
                // be annihilated, leaving us with independent replicas.
                if r0.system_lowest_max_energy.is_some() && r0.system.energy() < r1.max_energy {
                    std::mem::swap(&mut r0.system, &mut r1.system);
                    std::mem::swap(
                        &mut r0.system_lowest_max_energy,
                        &mut r1.system_lowest_max_energy,
                    );
                    if r1.system_lowest_max_energy.unwrap() > r1.max_energy {
                        r1.unique_visitors += 1;
                        r1.system_lowest_max_energy = Some(r1.max_energy);
                    }
                }
            }
        });
        const INDEPENDENT_SYSTEMS_BEFORE_NEW_BIN: u64 = 8;
        if let Some((cutoff, energy)) = self
            .replicas
            .last()
            .and_then(|r| Some((r.cutoff_energy, r.system.energy())))
        {
            if energy < cutoff {
                self.median.add_energy(energy);
            }
        }
        let new_replica = if let Some(r) = self.replicas.last() {
            if r.unique_visitors >= INDEPENDENT_SYSTEMS_BEFORE_NEW_BIN {
                let mean_below = r.below_total / r.below_count as f64;
                if mean_below + self.min_T < r.cutoff_energy && r.system.energy() < r.cutoff_energy
                {
                    let median_below = self.median.median();
                    self.median.reset(r.system.energy());
                    let mut newr = Replica::new(
                        r.cutoff_energy,
                        median_below,
                        r.system.clone(),
                        self.rng.clone(),
                    );
                    newr.translation_scale = r.translation_scale;
                    println!(
                        "      New bin: {:5} with mean {:5} and max {:.5}",
                        median_below.pretty(),
                        mean_below.pretty(),
                        r.cutoff_energy.pretty()
                    );
                    self.rng.jump();
                    Some(newr)
                } else {
                    None
                }
            } else {
                None
            }
        } else {
            unreachable!()
        };
        if let Some(r) = new_replica {
            self.replicas.push(r);
        }

        self.replicas.par_iter_mut().for_each(|r| r.occasional_update());

        for _ in 0..self.replicas.len() as u64 * steps {
            self.moves += 1;

            let movie_time = self.movie.shall_i_save(self.moves);
            if movie_time {
                self.movie.save_frame(&self.save_as, self.moves, &self);
            }
            if self.report.am_all_done(self.moves) {
                self.checkpoint();
                println!("All done!");
                ::std::process::exit(0);
            }
            if self.save.shall_i_save(self.moves) || movie_time {
                self.checkpoint();
            }
        }
    }
}
