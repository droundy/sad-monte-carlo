//! An implementation of SAD (Statistical Association, Dynamical version).
//!
//! This version using a [Binning]() approach so that we can use
//! something nicer than a histogram.  The idea is also to use this
//! for a 2D histogram, which will be even fancier.

#![allow(non_snake_case)]

use super::*;
use crate::system::*;

use super::plugin::Plugin;
use rand::{Rng, SeedableRng};
use std::default::Default;

use dimensioned::Dimensionless;

/// Parameters to configure the moves.
#[derive(Serialize, Deserialize, Debug, AutoArgs)]
pub enum MoveParams {
    /// The rms distance of moves
    TranslationScale(Length),
    /// Adjust translation scale to reach acceptance rate.
    AcceptanceRate(f64),
}

/// The parameters needed to configure a simulation.
#[derive(Debug, AutoArgs)]
pub struct EnergyMCParams {
    /// The rejection fraction.
    pub f: f64,
    /// The lowest temperature of interest
    min_T: Energy,
    /// After gamma drops to this value, begin a "production" run
    min_gamma: Option<f64>,
    /// The seed for the random number generator.
    pub seed: Option<u64>,
    /// The highest energy to allow.
    max_allowed_energy: Option<Energy>,
    _moves: MoveParams,
    _report: plugin::ReportParams,
    _movies: plugin::MovieParams,
    _save: plugin::SaveParams,
}

impl Default for EnergyMCParams {
    fn default() -> Self {
        EnergyMCParams {
            min_gamma: None,
            min_T: 0.2 * units::EPSILON,
            f: 0.5,
            seed: None,
            max_allowed_energy: None,
            _moves: MoveParams::TranslationScale(0.05 * units::SIGMA),
            _report: plugin::ReportParams::default(),
            _save: plugin::SaveParams::default(),
            _movies: plugin::MovieParams::default(),
        }
    }
}

/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct EnergyMC<S> {
    /// The rejection fraction.
    pub f: f64,
    /// The minimum temperature of interest
    min_T: Energy,
    /// The system we are simulating.
    pub system: S,
    /// The current gamma value for each energy separator
    pub gamma: f64,
    /// The minimum gamma before we start a production run
    pub min_gamma: Option<f64>,
    /// The number of moves that have been made.
    pub moves: u64,
    /// The number of moves that have been accepted.
    pub accepted_moves: u64,
    /// The highest energy to allow.
    max_allowed_energy: Option<Energy>,
    /// The move plan
    pub move_plan: MoveParams,
    /// The current translation scale
    pub translation_scale: Length,
    /// The "recent" acceptance rate.
    pub acceptance_rate: f64,
    /// The random number generator.
    pub rng: crate::rng::MyRng,
    /// Where to save the resume file.
    pub save_as: ::std::path::PathBuf,
    report: plugin::Report,
    save: plugin::Save,
    movies: plugin::Movie,
    manager: plugin::PluginManager,

    /// Min energy for bins
    pub min_energy: Energy,
    /// Max energy
    pub max_energy: Energy,
    /// The relative sizes of the bins
    pub rel_bins: Vec<f64>,
    /// The normalization of the bin sizes
    pub bin_norm: f64,
    /// The log weights for the bins
    pub lnw: Vec<f64>,
    /// The histogram counts
    pub histogram: Vec<u64>,
    /// Whether we have visited here since decreasing gamma
    pub have_seen: Vec<bool>,
    /// The total energy in each bin
    pub total_energy: Vec<Energy>,
}

fn min_of(stuff: &[f64]) -> f64 {
    stuff.iter().cloned().fold(0. / 0., f64::min)
}
fn max_of(stuff: &[f64]) -> f64 {
    stuff.iter().cloned().fold(0. / 0., f64::max)
}

impl<S: MovableSystem + ConfirmSystem> EnergyMC<S> {
    /// This decides whether to reject the move based on the actual
    /// method in use.
    fn reject_move(&mut self, e1: Energy, e2: Energy) -> bool {
        let lnw1 = self.lnw[self.e_to_idx(e1)];
        let lnw2 = self.lnw[self.e_to_idx(e2)];
        lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp()
    }
    fn e_to_idx(&self, energy: Energy) -> usize {
        let mut e = self.max_energy;
        let dedw = (self.max_energy - self.min_energy)/self.bin_norm;
        for (i, w) in self.rel_bins.iter().cloned().enumerate() {
            if energy > e {
                return i;
            }
            e -= dedw*w;
        }
        if energy > e {
            return self.rel_bins.len();
        }
        self.rel_bins.len()+1
    }
    /// This updates the lnw based on the actual method in use.
    fn update_weights(&mut self, energy: Energy) {
        let i = self.e_to_idx(energy);
        self.histogram[i] += 1;
        self.total_energy[i] += energy;
        if i == 0 {
            // We are in the widthless arbitrarily high-energy bin.
            let de = (self.max_energy - self.min_energy)*self.gamma;
            self.max_energy += de;
        } else if i == self.rel_bins.len()+1 {
            // We are in the low-energy bin.
            let de = (self.max_energy - self.min_energy)*self.gamma;
            self.min_energy -= de;

            // Now decide whether we want a new bin here.
            let min_de: Energy = -self.min_T*self.f.ln();
            assert!(min_de > Energy::new(0.));
            let min_w = *(self.bin_norm*min_de/(self.max_energy-self.min_energy)).value();
            if min_w < 0. {
                println!(" {} vs {}", self.max_energy, self.min_energy);
                println!("this gives min_w of {}", min_w);
                assert!(false);
            }
            let de = self.min_energy - self.total_energy[i]/self.histogram[i] as f64;
            if de > self.min_T && min_of(&self.rel_bins) > min_w {
                self.rel_bins.push(min_w);
                self.bin_norm += min_w;
                println!("new norm is {} from {}", self.bin_norm, min_w);
                self.histogram.push(0);
                self.total_energy.push(Energy::new(0.));
                self.have_seen.push(false);
                let lnw_last = self.lnw[i-1];
                self.lnw[i] = lnw_last + self.f.ln(); // make the last bin a regular bin
                self.lnw.push(lnw_last + self.f.ln() + (self.f/(1. - self.f)).ln());
                for h in self.have_seen.iter_mut() {
                    *h = false;
                }
                self.gamma = 0.25; // reset gamma since we have just discovered something potentially important.
                println!("opened up a new bin: current energy {}", energy.pretty());
                Logger.log(self, &self.system);
            }
        } else {
            let w = self.rel_bins[i-1];
            let dw = self.gamma*w;
            self.rel_bins[i-1] -= dw;
            self.bin_norm -= dw;
            let de = (self.max_energy - self.min_energy)*self.gamma/ self.rel_bins.len() as f64;
            self.max_energy -= de;
            self.min_energy += de;
        }

        if !self.have_seen[i] {
            self.have_seen[i] = true;
            if self.have_seen.iter().all(|h| *h) {
                // We have visited all the energies, so we can reduce gamma.
                for h in self.have_seen.iter_mut() {
                    *h = false;
                }
                self.gamma *= 0.5;
                // Let's re-normalize the bin widths occasionally, so we won't risk
                // underflow due to the bins getting unreasonably small.
                for w in self.rel_bins.iter_mut() {
                    *w /= self.bin_norm;
                }
                // The following should be about 1, but due to roundoff error it is
                // not, and we don't want any (unavoidable) errors in our binning.
                self.bin_norm = self.rel_bins.iter().cloned().sum::<f64>();
                if self.gamma < self.rel_bins.len() as f64 / self.moves as f64 {
                    self.gamma = self.rel_bins.len() as f64 / self.moves as f64;
                } else {
                    println!("found all energies: gamma -> {}", self.gamma)
                }
            }
        }
    }
}

impl<S: MovableSystem> MonteCarlo for EnergyMC<S> {
    type Params = EnergyMCParams;
    type System = S;
    fn from_params(params: EnergyMCParams, mut system: S, save_as: ::std::path::PathBuf) -> Self {
        let mut rng = crate::rng::MyRng::seed_from_u64(params.seed.unwrap_or(0));
        // Let's spend a little effort getting an energy that is
        // within our range of interest.  We are only aiming downward,
        // because it is unusual that we have trouble getting an
        // energy that is high enough.
        if let Some(maxe) = params.max_allowed_energy {
            for _ in 0..1e8 as u64 {
                if let Some(newe) = system.plan_move(&mut rng, 0.05 * units::SIGMA) {
                    if newe < system.energy() {
                        system.confirm();
                    }
                    if system.energy() < maxe {
                        break;
                    }
                }
            }
        }
        let f = params.f;
        let min_T = params.min_T;
        EnergyMC {
            min_gamma: params.min_gamma,
            moves: 0,
            min_T,
            f,
            accepted_moves: 0,
            acceptance_rate: 0.5, // arbitrary starting guess.
            max_energy: system.energy(),
            min_energy: system.energy() - 2.0*params.min_T,
            max_allowed_energy: params.max_allowed_energy,

            rel_bins: vec![1.0,1.0],
            bin_norm: 2.0,
            gamma: 0.25,
            histogram: vec![0; 4],
            have_seen: vec![false; 4],
            lnw: vec![
                0.,
                (1. - f).ln(),
                (1. - f).ln() + f.ln(),
                (1. - f).ln() + f.ln() + (f / (1. - f)).ln(),
            ],
            total_energy: vec![Energy::new(0.0); 4],

            translation_scale: match params._moves {
                MoveParams::TranslationScale(x) => x,
                _ => 0.05 * units::SIGMA,
            },
            move_plan: params._moves,
            system: system,

            rng,
            save_as: save_as,
            report: plugin::Report::from(params._report),
            save: plugin::Save::from(params._save),
            movies: plugin::Movie::from(params._movies),
            manager: plugin::PluginManager::new(),
        }
    }
    fn update_from_params(&mut self, params: Self::Params) {
        self.report.update_from(params._report);
        self.save.update_from(params._save);
    }

    fn move_once(&mut self) {
        self.moves += 1;
        // self.system.collect_data();
        if self.moves % 100000000 == 0 {
            self.system.verify_energy();
        }
        let e1 = self.system.energy();
        let recent_scale = (1.0 / self.moves as f64).sqrt();
        self.acceptance_rate *= 1. - recent_scale;
        if let Some(e2) = self.system.plan_move(&mut self.rng, self.translation_scale) {
            let mut out_of_bounds = false;
            if let Some(maxe) = self.max_allowed_energy {
                out_of_bounds = e2 > maxe && e2 > e1;
            }
            if !out_of_bounds {
                if !self.reject_move(e1, e2) {
                    self.accepted_moves += 1;
                    self.acceptance_rate += recent_scale;
                    self.system.confirm();
                }
            }
        }
        let energy = self.system.energy();

        self.update_weights(energy);

        let plugins = [
            &self.report as &dyn Plugin<Self>,
            &Logger,
            &self.movies,
            &self.save,
        ];
        self.manager.run(self, &self.system, &plugins);
    }
    fn system(&self) -> &Self::System {
        &self.system
    }
    fn system_mut(&mut self) -> &mut Self::System {
        &mut self.system
    }
    fn num_moves(&self) -> u64 {
        self.moves
    }
    fn num_accepted_moves(&self) -> u64 {
        self.accepted_moves
    }
    fn save_as(&self) -> ::std::path::PathBuf {
        self.save_as.clone()
    }
}

#[derive(Serialize, Deserialize, Debug)]
struct Logger;
impl<S: MovableSystem> Plugin<EnergyMC<S>> for Logger {
    fn log(&self, mc: &EnergyMC<S>, sys: &S) {
        // print!("    ");
        let print_am_here = |i| {
            if i == mc.e_to_idx(sys.energy()) {
                print!(" >");
            } else {
                print!("  ");
            }
        };
        let mut etop = mc.max_energy;
        let mut ebot;
        let mut i = 0;
        print_am_here(i);
        println!("  {:8.5} -> infty   : {}", etop.pretty(), mc.histogram[i]);
        let dedw = (mc.max_energy - mc.min_energy)/mc.bin_norm;
        while i < mc.rel_bins.len() {
            ebot = etop - dedw*mc.rel_bins[i];
            print_am_here(i+1);
            println!("  {:8.5} -> {:8.5}: {}", ebot.pretty(), etop.pretty(), mc.histogram[i+1]);
            i += 1;
            etop = ebot;
        }
        print_am_here(i+1);
        println!("  -infty   -> {:8.5}: {}", etop.pretty(), mc.histogram[i+1]);
        println!(
            " [gamma = {:.2}]",
            crate::prettyfloat::PrettyFloat(mc.gamma)
        );
        println!("        norm: {} vs. {}", mc.bin_norm, mc.rel_bins.iter().cloned().sum::<f64>());
    }
}
