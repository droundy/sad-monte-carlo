//! An implementation of a new algorithm, with no decent name.

#![allow(non_snake_case)]

use super::*;
use crate::system::*;
use dimensioned::Dimensionless;
use rayon::prelude::*;

use rand::{Rng, SeedableRng};
use std::default::Default;

/// The parameters needed to configure a simulation.
#[derive(Debug, AutoArgs, Clone)]
pub struct MCParams {
    /// A temperature of interest
    T: Vec<f64>,
    /// The seed for the random number generator.
    pub seed: Option<u64>,
    /// report input
    pub _movies: plugin::MovieParams,
    /// saving
    pub _save: plugin::SaveParams,
    /// report input
    pub _report: plugin::ReportParams,
}

impl Default for MCParams {
    fn default() -> Self {
        MCParams {
            T: vec![0.2],
            seed: None,
            _save: plugin::SaveParams::default(),
            _movies: plugin::MovieParams::default(),
            _report: plugin::ReportParams::default(),
        }
    }
}

/// A simple simulation, with a maximum energy, a single system, and an energy
/// divider.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Replica<S> {
    /// The temperature of this replica
    pub T: Energy,
    /// The number of rejected moves since we adjusted the translation_scale.
    pub rejected_count: u64,
    /// The number of accepted moves since we adjusted the translation_scale
    pub accepted_count: u64,
    /// The lowest `max_energy` that the system has visited, and the system itself
    pub system: S,
    /// The random number generator.
    pub rng: crate::rng::MyRng,

    /// The total energy over all moves
    pub total_energy: Energy,
    /// The total energy squared over all moves
    pub total_energy_squared: EnergySquared,

    /// The current translation scale
    pub translation_scale: Length,
}

impl<S: MovableSystem> Replica<S> {
    fn new(T: Energy, system: S, rng: crate::rng::MyRng) -> Self {
        Replica {
            T,
            system,
            rejected_count: 0,
            accepted_count: 0,

            total_energy: Energy::new(0.0),
            total_energy_squared: EnergySquared::new(0.0),

            translation_scale: Length::new(1.0),
            rng,
        }
    }
    fn energy(&self) -> Energy {
        self.system.energy()
    }
    fn run_once(&mut self) {
        if let Some(e) = self.system.plan_move(&mut self.rng, self.translation_scale) {
            let beta_delta_e = *((e - self.energy())/self.T).value();
            if beta_delta_e < 0.0 {
                // always allow energy to drop
                self.system.confirm();
                self.accepted_count += 1;//Probably not right?
            } else if self.rng.gen::<f64>() < (-beta_delta_e).exp() {
                self.system.confirm();
                self.accepted_count += 1;//Probably not right?
            } else {
                self.rejected_count += 1;
            }
            let e = self.energy();
            // collect data
            self.total_energy += e;
            self.total_energy_squared += e * e;
        }
    }
    fn printme(&self) {
        // print!(
        //     "   {} {:2.1}% > {:9.5} {:.2} unique, 𝚫E = {:.2} ({:.3}% up)",
        //     if self.energy().is_none() { ":(" } else { "  " },
        //     crate::prettyfloat::PrettyFloat(100.0 * self.above_fraction()),
        //     self.cutoff_energy.pretty(),
        //     crate::prettyfloat::PrettyFloat(self.unique_visitors as f64),
        //     (self.max_energy - self.cutoff_energy).pretty(),
        //     crate::prettyfloat::PrettyFloat(
        //         100.0 * self.upwelling_count as f64
        //             / (self.above_count as f64 + self.below_count as f64)
        //     ),
        // );
        // if let Some((sys, _)) = &self.system_with_lowest_max_energy {
        //     sys.print_debug();
        // }
        println!();
    }
}

/// A simulation with many replicas
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct MC<S> {
    /// The temperatures of interest
    T: Vec<Energy>,
    /// The random number generator.
    pub rng: crate::rng::MyRng,
    /// Where to save the resume file.
    pub save_as: ::std::path::PathBuf,
    /// The number of MC moves we have done
    pub moves: u64,

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
        let T:Vec<Energy> = params.T.iter().cloned().map(|e| Energy::new(e)).collect();

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
        let mut replicas = Vec::<Replica<S>>::with_capacity(T.len());
        for t in T{
            replicas.push(Replica::new(t, system, rng.clone()));
        }
            
        rng.jump();
        MC {
            T,
            replicas,
            moves: 0,

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
        self.report.print(self.moves, 0);
        println!(
            "    [{} walkers in {} zones]",
            self.replicas.len(),
            self.replicas.len()
        );

        let num_replicas = self.replicas.len();
        let mut need_dots = false;
        for (which, r) in self.replicas.iter().enumerate() {
            if which < 5 || which + 15 >= num_replicas {
                need_dots = true;
                r.printme();
            } else if need_dots {
                println!("      ...");
                need_dots = false;
            }
        }

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
        // The unwrap below is safe because the unbounded high energy bin will always be occupied.
        let steps = self.replicas[0].system.min_moves_to_randomize();

        let these_moves = std::sync::atomic::AtomicU64::new(0);

        // Run a few steps of the simulation for each replica
        self.replicas.par_iter_mut().for_each(|r| {
            for _ in 0..steps {
                r.run_once();
            }
        });

        // Now let us try swapping if we can.
        let iterator = if self.rng.gen::<bool>() {
            // Try swapping odd ones
            self.replicas.chunks_exact_mut(2)
        } else {
            // Try swapping even ones.
            self.replicas[1..].chunks_exact_mut(2)
        };
        iterator.for_each(|chunk| {
            // for chunk in iterator {
            if let [r0, r1] = chunk {
                // We will swap them if both systems can go into the lower bin.
                let de_db = *((r0.energy() - r1.energy())*(1./r0.T - 1./r1.T)).value();
                if de_db >= 0. {
                    std::mem::swap(&mut r0.system, &mut r1.system);
                } else if r1.rng.gen::<f64>() < de_db.exp(){
                    std::mem::swap(&mut r0.system, &mut r1.system)
                }//else don't swap
            }
        });

        let mut moves = self.moves;
        let these_moves = these_moves.load(std::sync::atomic::Ordering::Relaxed);
        self.moves += these_moves;
        for _ in 0..these_moves {
            moves += 1;

            let movie_time = self.movie.shall_i_save(moves);
            if movie_time {
                self.movie.save_frame(&self.save_as, moves, &self);
            }
            if self.report.am_all_done(moves, 0) {
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
