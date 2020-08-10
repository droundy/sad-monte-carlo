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

use energy_transposed::{EnergyMC, EnergyMCParams};

#[derive(Serialize, Deserialize, Debug, Clone, Copy, AutoArgs)]
enum ReplicaPlan {
    /// Time between replica swaps
    SwapTime(u64),
    /// Do not ever swap replicas
    NoSwap,
}
impl ReplicaPlan {
    fn it_is_time(self, moves: u64) -> bool {
        match self {
            ReplicaPlan::SwapTime(period) => moves % period == 0,
            ReplicaPlan::NoSwap => false,
        }
    }
}
/// The parameters needed to configure a simulation.
#[derive(Debug, Clone, AutoArgs)]
pub struct GrandMCParams {
    /// The normal parameters
    pub _energy: EnergyMCParams,
    /// The maximum number of atoms
    max_N: usize,
    /// Time between replicas
    _replica_plan: Option<ReplicaPlan>,
}

impl Default for GrandMCParams {
    fn default() -> Self {
        GrandMCParams {
            _energy: EnergyMCParams::default(),
            max_N: 5,
            _replica_plan: None,
        }
    }
}

/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct GrandMC<S> {
    /// The energy mcs
    pub mc: Vec<EnergyMC<S>>,
    /// The number of moves that have been made.
    pub moves: u64,
    /// Time between replicas
    replica_plan: ReplicaPlan,
    /// The number of moves that have been accepted.
    pub accepted_swaps: u64,
    /// Where to save the resume file.
    pub save_as: ::std::path::PathBuf,
    /// The random number generator.
    pub random: crate::rng::MyRng,

    #[serde(skip, default)]
    next_output: std::cell::Cell<u64>,
    /// This is when and where the simulation started.
    #[serde(skip, default)]
    start: std::cell::Cell<Option<(std::time::Instant, u64)>>,
    /// How frequently to save...
    #[serde(default)]
    save_time_seconds: Option<f64>,
}

impl<S: Clone + ConfirmSystem + GrandReplicaSystem + Sync + Send> GrandMC<S> {
    /// read from params
    pub fn from_params(_mc: GrandMCParams, _sys: S, save_as: std::path::PathBuf) -> Self {
        let mut random = crate::rng::MyRng::seed_from_u64(_mc._energy.seed.unwrap_or(0));
        let mut mc = Vec::new();
        // First create the output directory if it does not yet exist.
        std::fs::create_dir_all(save_as.file_stem().unwrap()).ok();
        let mut energy_params = _mc._energy.clone();
        energy_params._report.quiet = true;
        for n in 1.._mc.max_N {
            let mut s = EnergyMC::from_params(
                energy_params.clone(),
                _sys.clone(),
                format!(
                    "{}/N-{}.cbor",
                    save_as.file_stem().unwrap().to_string_lossy(),
                    n
                )
                .into(),
            );
            for _ in 0..n {
                while s.system.plan_add(&mut random).is_none() {}
                s.system.confirm();
            }
            mc.push(s);
            println!("Finished initializing system with {} atoms...", n);
        }
        GrandMC {
            mc,
            moves: 0,
            replica_plan: _mc._replica_plan.unwrap_or(ReplicaPlan::SwapTime(1)),
            accepted_swaps: 0,
            save_as: save_as.clone(),
            random,

            next_output: std::cell::Cell::new(0),
            start: std::cell::Cell::new(Some((std::time::Instant::now(), 0))),
            save_time_seconds: Some(60.0 * 60.0), // 1 hour
        }
    }
    /// Create a new simulation from command-line flags.
    pub fn from_args<A: AutoArgs + Into<S>>() -> Self {
        println!("git version: {}", VERSION);
        match <Params<GrandMCParams, A>>::from_args() {
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
                        for system in s.mc.iter_mut() {
                            system.update_from_params(_mc._energy.clone());
                            system.system_mut().update_caches();
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
    pub fn checkpoint(&self) {
        let f = AtomicFile::create(&self.save_as)
            .expect(&format!("error creating file {:?}", self.save_as));
        match self.save_as.extension().and_then(|x| x.to_str()) {
            Some("yaml") => serde_yaml::to_writer(&f, self).expect("error writing checkpoint?!"),
            Some("json") => serde_json::to_writer(&f, self).expect("error writing checkpoint?!"),
            Some("cbor") => serde_cbor::to_writer(&f, self).expect("error writing checkpoint?!"),
            _ => panic!("I don't know how to create file {:?}", self.save_as),
        }
    }
    /// Save, if it is time.
    fn save(&self) {
        if self.next_output.get() <= self.moves {
            self.checkpoint();
            if let Some(period) = self.save_time_seconds {
                match self.start.get() {
                    Some((start_time, start_iter)) => {
                        let moves = self.moves;
                        let runtime = start_time.elapsed();
                        let time_per_move = duration_to_secs(runtime) / (moves - start_iter) as f64;
                        let moves_per_period = 1 + (period / time_per_move) as u64;
                        self.next_output.set(moves + moves_per_period);
                    }
                    None => {
                        self.start
                            .set(Some((std::time::Instant::now(), self.moves)));
                        self.next_output.set(self.moves + (1 << 20));
                    }
                }
            } else {
                self.next_output.set(self.next_output.get() * 2)
            }
        }
    }
    /// Run a simulation
    pub fn run_once(&mut self) {
        self.mc.par_iter_mut().for_each(|mc| {
            for _ in 0..mc.system.num_atoms() {
                mc.move_once();
            }
        });

        if self.replica_plan.it_is_time(self.moves) {
            for (i, mc) in self.mc.iter().enumerate() {
                assert_eq!(i + 1, mc.system.num_atoms());
            }
            let iterator = if self.random.gen::<f64>() > 0.5 {
                // Try swapping odd onces
                self.mc.chunks_exact_mut(2)
            } else {
                // Try swapping even ones.
                self.mc[1..].chunks_exact_mut(2)
            };
            self.accepted_swaps += iterator.map(|chunk| {
            // for chunk in iterator {
                if let [mc0, mc1] = chunk {
                    assert_eq!(mc0.system.num_atoms() + 1, mc1.system.num_atoms());
                    if let Some((which, e1, e0)) =
                        mc1.system.plan_swap_atom(&mc0.system, &mut mc0.rng)
                    {
                        let old_lnw =
                            mc1.e_to_lnw(mc1.system.energy()) + mc0.e_to_lnw(mc0.system.energy());
                        let new_lnw = mc1.e_to_lnw(e0) + mc0.e_to_lnw(e1);
                        if (old_lnw - new_lnw).exp() > mc0.rng.gen::<f64>() {
                            // println!(
                            //     "I am swapping an atom from {} to {}!",
                            //     mc1.system.num_atoms(),
                            //     mc0.system.num_atoms()
                            // );
                            mc1.system.swap_atom(&mut mc0.system, which);
                            std::mem::swap(&mut mc0.system, &mut mc1.system);
                            // self.accepted_swaps += 1;
                            return 1; // one swap accepted!
                        }
                    }
                }
                0
            }).sum::<u64>();
        }
        self.moves += 1;
        self.save();
    }
}

fn duration_to_secs(t: std::time::Duration) -> f64 {
    t.as_secs() as f64 + t.subsec_nanos() as f64 * 1e-9
}
