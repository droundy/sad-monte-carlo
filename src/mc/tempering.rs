//! An implementation of a new algorithm, with no decent name.

#![allow(non_snake_case)]

use super::*;
use crate::system::*;
use dimensioned::{Dimensionless}; //TODO: What was the unused import si::E 
use rayon::prelude::*;

use rand::{Rng, SeedableRng};
use std::default::Default;

/// The parameters needed to configure a simulation.
#[derive(Debug, AutoArgs, Clone)]
pub struct MCParams {
    /// The temperatures of interest
    T: Vec<f64>,
    /// The number of Canonical steps to take between swaps
    canonical_steps: u64,
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
            T: vec![0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024],
            canonical_steps: 1,
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
    /// The number of rejected moves
    pub rejected_count: u64,
    /// The number of accepted moves
    pub accepted_count: u64,
    /// The number of rejected swaps involving this replica
    pub rejected_swap_count: u64,
    /// The number of accepted swaps involving this replica
    pub accepted_swap_count: u64,
    /// The number of ignored moves involving this replica
    pub ignored_count: u64,
    /// The system to analyze
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

            rejected_swap_count: 0,
            accepted_swap_count: 0,

            ignored_count: 0,

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
            if beta_delta_e < 0.0 || self.rng.gen::<f64>() < (-beta_delta_e).exp() {
                self.system.confirm();
                self.accepted_count += 1;
            } else {
                self.rejected_count += 1;
            }
            let e = self.energy();
            // collect data
            self.total_energy += e;
            self.total_energy_squared += e * e;
            if e >= Energy::new(0.){
                self.ignored_count += 1
            }
        }
    }
    fn printme(&self) {
        println!(
            "E = {:3.1},   T = {:2.1}",
            self.energy().pretty(),
            self.T.pretty()
        );
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
    /// The number of Canonical steps to take between swaps
    pub canonical_steps: u64,
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
    fn from_params(params: MCParams, system: S, save_as: ::std::path::PathBuf) -> Self {
        let mut rng = crate::rng::MyRng::seed_from_u64(params.seed.unwrap_or(0));
        println!("T args{:?}", params.T);
        let T:Vec<Energy> = params.T.iter().cloned().map(|e| Energy::new(e)).collect();
        

        let mut replicas = Vec::<Replica<S>>::with_capacity(T.len());
        for t in T.iter().copied() {
            replicas.push(Replica::new(t, system.clone(), rng.clone()));
            println!("Creating new Replica with temperature  {:3}", t);
        } 
        rng.jump();
        MC {
            T,
            replicas,
            moves: 0,
            canonical_steps: params.canonical_steps,

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
            "    [{} walkers]",
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
        let steps = self.replicas[0].system.min_moves_to_randomize()*self.canonical_steps;

        let these_moves = std::sync::atomic::AtomicU64::new(0);

        // Run a few steps of the simulation for each replica
        self.replicas.par_iter_mut().for_each(|r| {
            these_moves.fetch_add(steps, std::sync::atomic::Ordering::Relaxed);
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
                // We will swap them if both systems accept the swap.
                let de_db = *((r0.energy() - r1.energy())*(1./r0.T - 1./r1.T)).value();
                if de_db >= 0. || r1.rng.gen::<f64>() < de_db.exp() {
                    r0.accepted_swap_count += 1;
                    r1.accepted_swap_count += 1;
                    std::mem::swap(&mut r0.system, &mut r1.system);
                }else{//else don't swap
                    r0.rejected_swap_count += 1;
                    r1.rejected_swap_count += 1;
                }
                let e0=r0.energy();
                let e1=r1.energy();
                r0.total_energy += e0;
                r1.total_energy += e1;

                r0.total_energy_squared += e0 * e0;
                r1.total_energy_squared += e1 * e1;
                if e0 >= Energy::new(0.){
                    r0.ignored_count += 1
                }
                if e1 >= Energy::new(0.){
                    r1.ignored_count += 1
                }
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
