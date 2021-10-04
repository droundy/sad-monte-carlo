//! An implementation of a new algorithm, with no decent name.

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
    /// The number of independent systems required before creating a new bin (default 64).
    pub independent_systems_before_new_bin: Option<u64>,
    /// report input
    pub _movies: plugin::MovieParams,
    /// report input
    pub _save: plugin::SaveParams,
    /// report input
    pub _report: plugin::ReportParams,
    /// Use the mean as an approximation of the median
    pub mean_for_median: bool,
    /// Always obey detailed balance as far as we can
    pub always_balance: bool,
}

impl Default for MCParams {
    fn default() -> Self {
        MCParams {
            min_T: 0.2 * units::EPSILON,
            seed: None,
            independent_systems_before_new_bin: None,
            _save: plugin::SaveParams::default(),
            _movies: plugin::MovieParams::default(),
            _report: plugin::ReportParams::default(),
            mean_for_median: false,
            always_balance: false,
        }
    }
}

/// An estimator for a median
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct MedianEstimator {
    energies: Vec<Energy>,
}
const ESTIMATOR_SIZE: usize = 4096; // 1024;
impl MedianEstimator {
    fn new(e: Energy) -> Self {
        let mut energies = Vec::with_capacity(ESTIMATOR_SIZE);
        energies.push(e);
        MedianEstimator { energies }
    }
    fn reset(&mut self, e: Energy) {
        self.energies.truncate(0);
        self.energies.push(e);
    }
    fn add_energy(&mut self, e: Energy, rng: &mut crate::rng::MyRng) {
        if self.energies.len() < ESTIMATOR_SIZE {
            self.energies.push(e);
        } else {
            if rng.gen::<f64>() < 1.0 / (self.energies.len() as f64 + 1.0) {
                let i = rng.gen_range(0, self.energies.len());
                self.energies[i] = e;
            }
        }
    }
    fn median(&mut self) -> Energy {
        self.energies.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let middle = self.energies.len() / 2;
        let e_middle = self.energies[middle];
        // We attempt to pick a median that is between a pair of our energies,
        // so that if we have a discrete set of possible energies our median
        // will not be one of them, which would be a bad place to put a bin
        // divider.
        if let Some(next) = self.energies[middle + 1..]
            .iter()
            .cloned()
            .filter(|&e| e != e_middle)
            .next()
        {
            0.5 * (e_middle + next)
        } else if let Some(prev) = self.energies[..middle]
            .iter()
            .rev()
            .cloned()
            .filter(|&e| e != e_middle)
            .next()
        {
            0.5 * (e_middle + prev)
        } else {
            println!("not happy to pick our only ever energy as the median... :(");
            e_middle
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
    /// The total count from systems that have already reached the lowest energy
    /// bin.
    #[serde(default)]
    pub upwelling_count: u64,
    /// The total energy above (for computing the mean)
    pub above_total: Energy,
    /// The total energy below (for computing the mean)
    pub below_total: Energy,
    /// The total energy squared above (for computing the standard deviation)
    pub above_total_squared: EnergySquared,
    /// The total energy squared below (for computing the standard deviation)
    pub below_total_squared: EnergySquared,
    /// Any extra data we might want to collect, total and count
    pub above_extra: HashMap<Interned, (f64, u64)>,
    /// The lowest `max_energy` that the system has visited, and the system itself
    pub system_with_lowest_max_energy: Option<(S, Energy)>,
    /// The number of completely and provably indepdendent systems that have
    /// us
    pub unique_visitors: u64,
    /// Are we collecting data for now?
    pub collecting_data: bool,
    /// The random number generator.
    pub rng: crate::rng::MyRng,

    /// The current translation scale
    pub translation_scale: Length,
}

impl<S: MovableSystem + Clone> Replica<S> {
    fn new(max_energy: Energy, cutoff_energy: Energy, system: S, rng: crate::rng::MyRng) -> Self {
        let mut r = Replica::empty(max_energy, cutoff_energy, rng);
        r.translation_scale = system.max_size();
        r.system_with_lowest_max_energy = Some((system, max_energy));
        r.unique_visitors = 1;
        r
    }
    fn empty(max_energy: Energy, cutoff_energy: Energy, rng: crate::rng::MyRng) -> Self {
        Replica {
            max_energy,
            cutoff_energy,
            above_count: 0,
            below_count: 0,
            upwelling_count: 0,
            rejected_count: 0,
            accepted_count: 0,
            above_total: Energy::new(0.0),
            below_total: Energy::new(0.0),
            above_total_squared: EnergySquared::new(0.0),
            below_total_squared: EnergySquared::new(0.0),
            above_extra: HashMap::new(),

            translation_scale: Length::new(1.0),
            system_with_lowest_max_energy: None,
            unique_visitors: 0,
            collecting_data: true,
            rng,
        }
    }
    fn decimate(&mut self) {
        self.upwelling_count = 0;
        if self.above_count > 1 {
            self.above_total /= self.above_count as f64;
            self.above_total_squared /= self.above_count as f64;
            self.above_count = 1;
        }
        if self.below_count > 1 {
            self.below_total /= self.below_count as f64;
            self.below_total_squared /= self.below_count as f64;
            self.below_count = 1;
        }
        for (_, e) in self.above_extra.iter_mut() {
            if e.1 > 1 {
                e.0 /= e.1 as f64;
                e.1 = 1;
            }
        }
        if let Some((_, me)) = &mut self.system_with_lowest_max_energy {
            // Restart tracking unique visitors.
            *me = self.max_energy;
        }
        self.accepted_count = 1;
        self.rejected_count = 1;
        self.unique_visitors = 1;
    }
    fn above_fraction(&self) -> f64 {
        self.above_count as f64 / (self.above_count as f64 + self.below_count as f64)
    }
    fn energy(&self) -> Option<Energy> {
        self.system_with_lowest_max_energy
            .as_ref()
            .map(|(s, _)| s.energy())
    }
    fn run_once(&mut self, moves: u64, very_lowest_max_energy: Energy) {
        if let Some((system, lowest_max_energy)) = &mut self.system_with_lowest_max_energy {
            if self.max_energy.value_unsafe.is_finite() {
                if let Some(e) = system.plan_move(&mut self.rng, self.translation_scale) {
                    if e < self.max_energy {
                        system.confirm();
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
                system.randomize(&mut self.rng);
                *lowest_max_energy = self.max_energy;
            }
            let e = system.energy();
            if self.collecting_data {
                if e > self.cutoff_energy {
                    self.above_count += 1;
                    self.above_total += e;
                    self.above_total_squared += e * e;
                    for (k, d) in system.data_to_collect(moves).into_iter() {
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
                    self.below_total_squared += e * e;
                }
                if *lowest_max_energy == very_lowest_max_energy {
                    self.upwelling_count += 1;
                }
            }
        }
    }
    fn occasional_update(&mut self) {
        if let Some((system, _)) = &mut self.system_with_lowest_max_energy {
            if self.rejected_count > 128
                && self.accepted_count > 128
                && self.max_energy.value_unsafe.is_finite()
            {
                let acceptance_ratio = self.accepted_count as f64 / self.rejected_count as f64;
                let max_acceptance_ratio = system.min_moves_to_randomize() as f64;
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
    fn split(&mut self) -> Replica<S> {
        print!("Splitting ");
        self.printme();
        // let mean = self.above_total / self.above_count as f64;
        let middle = 0.5 * (self.max_energy + self.cutoff_energy);
        let mut next = Replica::empty(middle, self.cutoff_energy, self.rng.clone());
        next.rng.jump();
        next.translation_scale = self.translation_scale;
        *self = Replica {
            max_energy: self.max_energy,
            cutoff_energy: middle,
            above_count: 0,
            below_count: 0,
            upwelling_count: 0,
            rejected_count: 0,
            accepted_count: 0,
            above_total: Energy::new(0.0),
            below_total: Energy::new(0.0),
            above_total_squared: EnergySquared::new(0.0),
            below_total_squared: EnergySquared::new(0.0),
            above_extra: HashMap::new(),

            translation_scale: self.translation_scale,
            system_with_lowest_max_energy: std::mem::replace(
                &mut self.system_with_lowest_max_energy,
                None,
            ),
            unique_visitors: 1,
            collecting_data: true,
            rng: self.rng.clone(),
        };
        self.printme();
        next.printme();
        next
    }
    fn split_clone(&mut self) -> Replica<S> {
        print!("Cloning ");
        self.printme();
        let middle = 0.5 * (self.max_energy + self.cutoff_energy);
        let mut next = self.clone();
        next.cutoff_energy = middle;
        next.max_energy = self.cutoff_energy;
        next.rng.jump();
        // Ensure that neither clone can count as unique after this
        if let Some((_, lme)) = &mut self.system_with_lowest_max_energy {
            *lme = Energy::new(f64::NEG_INFINITY);
        }
        if let Some((_, lme)) = &mut next.system_with_lowest_max_energy {
            *lme = Energy::new(f64::NEG_INFINITY);
        }
        self.printme();
        next.printme();
        next
    }
    fn printme(&self) {
        print!(
            "   {} {:2.1}% > {:9.5} {:.2} unique, 𝚫E = {:.2} ({:.3}% up)",
            if self.energy().is_none() {
                ":("
            } else if self.system_with_lowest_max_energy.is_some()
                && self
                    .system_with_lowest_max_energy
                    .as_ref()
                    .unwrap()
                    .1
                    .value_unsafe
                    .is_infinite()
            {
                ":/"
            } else {
                "  "
            },
            crate::prettyfloat::PrettyFloat(100.0 * self.above_fraction()),
            self.cutoff_energy.pretty(),
            crate::prettyfloat::PrettyFloat(self.unique_visitors as f64),
            (self.max_energy - self.cutoff_energy).pretty(),
            crate::prettyfloat::PrettyFloat(
                100.0 * self.upwelling_count as f64
                    / (self.above_count as f64 + self.below_count as f64)
            ),
        );
        if let Some((sys, _)) = &self.system_with_lowest_max_energy {
            sys.print_debug();
        }
        println!();
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
    /// The number of independent systems required before we create a new bin
    #[serde(default)]
    pub independent_systems_before_new_bin: u64,

    /// An estimator for the median energy below the lowest bin
    pub median: MedianEstimator,

    /// Use the mean as an approximation of the median
    pub mean_for_median: bool,
    /// Always obey detailed balance as far as we can
    pub always_balance: bool,

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
            print!("Random: ");
            system.print_debug();
            println!(" E: {:.3}", system.energy().pretty());
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
            mean_for_median: params.mean_for_median,
            always_balance: params.always_balance,
            independent_systems_before_new_bin: params
                .independent_systems_before_new_bin
                .unwrap_or(64),

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
                            if let Some((system, _)) = &mut r.system_with_lowest_max_energy {
                                system.update_caches();
                            }
                        }
                        s.report.update_from(_mc._report);
                        s.save.update_from(_mc._save);
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
        self.report
            .print(self.moves, self.replicas.last().unwrap().unique_visitors);
        println!(
            "    [{} walkers in {} zones]",
            self.replicas
                .iter()
                .filter(|r| r.system_with_lowest_max_energy.is_some())
                .count(),
            self.replicas.len()
        );

        let num_replicas = self.replicas.len();
        let mut need_dots = false;
        let mut print_one_more = false;
        for (which, r) in self.replicas.iter().enumerate() {
            let percent = r.above_count as f64 / (r.above_count as f64 + r.below_count as f64);
            let am_crazy = percent > 0.75 || percent < 0.25 || r.energy().is_none();
            if which < 5 || which + 15 >= num_replicas || am_crazy || print_one_more {
                need_dots = true;
                print_one_more = am_crazy; // print one more bin after a crazy bin.
                r.printme();
            } else if need_dots {
                println!("      ...");
                need_dots = false;
            }
        }
        if let Some(r) = self.replicas.last() {
            let mean_below = r.below_total / r.below_count as f64;
            println!(
                "        mean_below: {:9.5}, T_below {:.2}",
                mean_below.pretty(),
                (r.cutoff_energy - mean_below).pretty()
            );
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
        // The unwrap below is safe because the unbounded high energy bin will always be occupied.
        let steps = self.replicas[0]
            .system_with_lowest_max_energy
            .as_ref()
            .unwrap()
            .0
            .min_moves_to_randomize();

        let these_moves = std::sync::atomic::AtomicU64::new(0);

        let any_missing = self.replicas.iter().any(|r| r.energy().is_none());
        let lowest_max_energy = self.replicas.last().unwrap().max_energy;
        let always_balance = self.always_balance;
        // Run a few steps of the simulation for each replica
        self.replicas.par_iter_mut().for_each(|r| {
            if any_missing
                && !always_balance
                && r.max_energy != lowest_max_energy
                && r.max_energy.value_unsafe.is_finite()
            {
                // If we have any bubbles, and we are not looking at the highest or lowest
                // zones, we should stop collecting data for a while to avoid violations
                // of detailed balance while we let the bubbles rise.
                r.collecting_data = always_balance;
            }
            if r.system_with_lowest_max_energy.is_some() {
                these_moves.fetch_add(steps, std::sync::atomic::Ordering::Relaxed);
                if r.max_energy.value_unsafe.is_finite() {
                    for i in 0..steps {
                        r.run_once(moves + i, lowest_max_energy);
                    }
                } else {
                    // For the unbounded high energy bin, we only need to fully randomize once.
                    // This is important, because fully randomizing is way more expensive
                    // than a single move (but even more effective).
                    r.run_once(moves, lowest_max_energy);
                }
            }
        });

        // Now let us try swapping if we can.
        if any_missing && !self.always_balance {
            // Bubble up as rapidly as posssible
            for i in (1..self.replicas.len()).rev() {
                if let [r0, r1] = &mut self.replicas[i - 1..i + 1] {
                    if r1.energy().is_none() {
                        if let Some(e) = r0.energy() {
                            if e < r1.max_energy {
                                if r0.max_energy.value_unsafe.is_finite() {
                                    std::mem::swap(
                                        &mut r0.system_with_lowest_max_energy,
                                        &mut r1.system_with_lowest_max_energy,
                                    );
                                } else {
                                    // If the higher-energy bin is unbounded, then leave the existing system
                                    // there, since it'll get randomized in a moment anyhow.
                                    r1.system_with_lowest_max_energy =
                                        r0.system_with_lowest_max_energy.clone();
                                }
                            }
                        }
                    }
                }
            }
        }
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
                // We will swap them if both systems can go into the lower bin.
                if (r0.energy().is_some() && r0.energy().unwrap() < r1.max_energy)
                    // If the higher energy is a bubble, move it down with 50% probability,
                    // so bubbles will freely diffuse.
                    || (always_balance && r0.energy().is_none() && r0.rng.gen())
                {
                    if r0.max_energy.value_unsafe.is_finite() {
                        std::mem::swap(
                            &mut r0.system_with_lowest_max_energy,
                            &mut r1.system_with_lowest_max_energy,
                        );
                        // We start collecting data after the zone has done an ordinary swap.
                        // This is a heuristic to avoid the bias that
                        // comes from systems dropping down to lower energies and not
                        // returning because we have created new low-energy zones.
                        r0.collecting_data = true;
                        r1.collecting_data = true;
                    } else {
                        // If the higher-energy bin is unbounded, then leave the existing system
                        // there, since it'll get randomized in a moment anyhow.
                        r1.system_with_lowest_max_energy = r0.system_with_lowest_max_energy.clone();
                    }
                    if let Some((_, me)) = &mut r1.system_with_lowest_max_energy {
                        if *me > r1.max_energy {
                            // if r1.collecting_data && !any_missing {
                            //     // Only increment the number of unique visitors if we're
                            //     // currently actually collecting statistics.  This may
                            //     // undercount the unique visitors, but that beats overcounting.
                            //     r1.unique_visitors += 1;
                            // }

                            // Always increment the number of unique visitors.  Otherwise we can never
                            // know how many there have been, and any visitors that pass by while
                            // we aren't collecting data will probably come back again later.
                            r1.unique_visitors += 1;
                            *me = r1.max_energy;
                        }
                    }
                }
            }
        });
        if let Some((cutoff, energy)) = self.replicas.last().and_then(|r| {
            r.system_with_lowest_max_energy
                .as_ref()
                .map(|(system, _)| (r.cutoff_energy, system.energy()))
        }) {
            if energy < cutoff {
                self.median.add_energy(energy, &mut self.rng);
            }
        }
        let independent_systems_before_new_bin = self.independent_systems_before_new_bin;
        let new_replica = if let Some(r) = self.replicas.last() {
            // We will create a new bin if we have had
            // INDEPENDENT_SYSTEMS_BEFORE_NEW_BIN unique visitors at the lowest
            // energy.  Also, we need the mean energy below our lowest cutoff to
            // be lower than min_T below that cutoff (so we aren't at lower
            // microcanonical temperature than min_T), the system energy must be
            // below the median so that we have a to put into our new bin.
            if let Some((system, _)) = &r.system_with_lowest_max_energy {
                if r.unique_visitors >= independent_systems_before_new_bin
                // && self.replicas.iter().all(|rr| rr.energy().is_some())
                {
                    let mean_below = r.below_total / r.below_count as f64;
                    if mean_below + self.min_T < r.cutoff_energy
                        && system.energy() < r.cutoff_energy
                    {
                        let median_below = if self.mean_for_median {
                            mean_below
                        } else {
                            self.median.median()
                        };
                        self.median.reset(median_below);
                        let mut newr = if always_balance {
                            let mut newr = r.clone();
                            newr.max_energy = r.cutoff_energy;
                            newr.cutoff_energy = median_below;
                            // Ensure that neither clone can count as unique after this
                            // if let Some((_, lme)) = &mut r.system_with_lowest_max_energy {
                            //     *lme = Energy::new(f64::NEG_INFINITY);
                            // }
                            if let Some((_, lme)) = &mut newr.system_with_lowest_max_energy {
                                *lme = Energy::new(f64::NEG_INFINITY);
                            }
                            newr
                        } else {
                            Replica::empty(r.cutoff_energy, median_below, self.rng.clone())
                        };
                        newr.translation_scale =
                            r.translation_scale * 0.5f64.powf(1.0 / system.dimensionality() as f64);
                        println!(
                            "      New bin: {:5} with mean {:5} and max {:.5} with {} unique",
                            median_below.pretty(),
                            mean_below.pretty(),
                            r.cutoff_energy.pretty(),
                            r.unique_visitors,
                        );
                        // println!("       below count: {}", r.below_count);
                        // println!("       above count: {}", r.above_count);
                        // println!("       max energy: {}", r.max_energy.pretty());
                        newr.rng.jump();
                        Some(newr)
                    } else {
                        None
                    }
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
            // We now choose to throw away data, to account for the fact that
            // we're now spending some time violating detailed balance.
            if !self.always_balance {
                for r in self.replicas.iter_mut() {
                    r.decimate();
                }
            }
        }

        self.replicas
            .par_iter_mut()
            .for_each(|r| r.occasional_update());

        if false && self.replicas.iter().all(|rr| rr.energy().is_some()) {
            // Split bins only when we've got the "decks" cleared in terms of bubbles.
            // This helps ensure we don't go overboard.
            let need_splitting = self
                .replicas
                .iter()
                .enumerate()
                .filter(|(_, r)| {
                    r.above_fraction() > 0.75
                        && 2.0 / (r.unique_visitors as f64).sqrt() < r.above_fraction() - 0.75
                })
                .map(|(i, _)| i)
                .collect::<Vec<_>>();
            for i in need_splitting.iter().rev().cloned() {
                let next = if always_balance {
                    self.replicas[i].split_clone()
                } else {
                    self.replicas[i].split()
                };
                self.replicas.insert(i + 1, next);
            }
        }

        let mut moves = self.moves;
        let these_moves = these_moves.load(std::sync::atomic::Ordering::Relaxed);
        self.moves += these_moves;
        for _ in 0..these_moves {
            moves += 1;

            let movie_time = self.movie.shall_i_save(moves);
            if movie_time {
                self.movie.save_frame(&self.save_as, moves, &self);
            }
            if self
                .report
                .am_all_done(moves, self.replicas.last().unwrap().unique_visitors)
            {
                self.checkpoint();
                println!("All done!");
                ::std::process::exit(0);
            }
            if self.save.shall_i_save(moves) || movie_time {
                self.checkpoint();
            }
        }
    }
}
