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
    /// The lowest energy to allow.
    min_allowed_energy: Option<Energy>,
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
            min_allowed_energy: None,
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
    pub gamma: Vec<f64>,
    /// The minimum gamma before we start a production run
    pub min_gamma: Option<f64>,
    /// The number of moves that have been made.
    pub moves: u64,
    /// The number of moves that have been accepted.
    pub accepted_moves: u64,
    /// The lowest energy to allow.
    min_allowed_energy: Option<Energy>,
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

    /// The energies dividing bins
    pub energies: Vec<Energy>,
    /// The log weights for the bins
    pub lnw: Vec<f64>,
    /// The histogram counts
    pub histogram: Vec<u64>,
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
        for (i, e) in self.energies.iter().cloned().enumerate() {
            if energy > e {
                return i;
            }
        }
        self.energies.len()
    }
    /// This updates the lnw based on the actual method in use.
    fn update_weights(&mut self, energy: Energy) {
        let i = self.e_to_idx(energy);
        let old_hist_here = self.histogram[i];
        self.histogram[i] += 1;
        self.total_energy[i] += energy;
        self.histogram.len();
        let de_left = if i == 0 {
            self.energies[i] - self.energies[i + 1]
        } else if i == 1 {
            self.energies[i - 1] - self.energies[i]
        } else if i == self.energies.len() {
            self.energies[i - 2] - self.energies[i - 1]
        } else if self.energies[i - 1] - self.energies[i]
            < self.energies[i - 2] - self.energies[i - 1]
        {
            self.energies[i - 1] - self.energies[i]
        } else {
            self.energies[i - 2] - self.energies[i - 1]
        };
        assert!(de_left > Energy::new(0.));
        let de_right = if i == self.energies.len() {
            self.energies[i - 2] - self.energies[i - 1]
        } else if i == self.energies.len() - 1 {
            self.energies[i - 1] - self.energies[i]
        } else if i == 0 {
            self.energies[i] - self.energies[i + 1]
        } else if self.energies[i - 1] - self.energies[i] < self.energies[i] - self.energies[i + 1]
        {
            self.energies[i - 1] - self.energies[i]
        } else {
            self.energies[i] - self.energies[i + 1]
        };
        assert!(de_right > Energy::new(0.));
        // println!("currently {}", energy);
        if i < self.energies.len() {
            // println!("adding {:.3} to {}", (de_right * self.gamma[i]).pretty(), i);
            self.energies[i] += de_right * self.gamma[i];
        }
        if i > 0 {
            // println!("subtracting {:.3} from {}", (de_left * self.gamma[i - 1]).pretty(), i-1);
            self.energies[i - 1] -= de_left * self.gamma[i - 1];
        }
        // We only are willing to add a new bin *after* we know that every bin has been
        // re-rexplored at least once since the last time we added a bin.  Otberwise, we
        // can accumulate far too many bins.
        if i == self.histogram.len() - 1 && self.gamma[i-1] < 0.25 {
            // Now decide whether to add a new bin
            let elow = self.total_energy[i] / self.histogram[i] as f64;
            if self.energies[i - 1] - elow > self.min_T {
                // We have room for another bin (we think)
                let T = self.energies[i - 1] - elow;
                let enew = self.energies[i - 1] + T*(1.-self.f).ln();
                self.energies.push(enew);
                self.gamma.push(0.25);
                self.lnw.pop(); // get rid of the last "big" bin
                let lnw_last = self.lnw.last().copied().unwrap();
                if self.lnw.len() > 1 {
                    self.lnw.push(lnw_last + self.f.ln()); // insert a regular bin
                    self.lnw
                        .push(lnw_last + self.f.ln() + (1. / self.f - 1.).ln());
                // insert a new "big" bin.
                } else {
                    // We had just two bins to start with, so we need to treat the first one specially
                    // to maintain the large "negative temperature" bin.
                    self.lnw.push(lnw_last + (1. - self.f).ln());
                    self.lnw.push(lnw_last + (1. / self.f - 1.).ln());
                }
                self.histogram.pop();
                self.total_energy.pop();
                for _ in 0..2 {
                    self.total_energy.push(Energy::new(0.));
                    self.histogram.push(0);
                }
                assert_eq!(self.lnw.len(), self.energies.len() + 1);
                assert_eq!(self.histogram.len(), self.lnw.len());
                assert_eq!(self.total_energy.len(), self.lnw.len());
            }
        }
        if let Some(min_gamma) = self.min_gamma {
            if max_of(&self.gamma) < min_gamma {
                // We are in a production run.
                return;
            }
        }
        for i in 0..self.energies.len() - 1 {
            if self.energies[i + 1] >= self.energies[i] {
                println!(
                    "craziness at i={} with {} and {}",
                    i,
                    self.energies[i + 1],
                    self.energies[i]
                );
            }
            assert!(self.energies[i + 1] < self.energies[i]);
        }
        if old_hist_here == 0 && self.histogram.len() > 2 && !self.histogram.iter().any(|&x| x == 0)
        {
            // We have explored all the bins!
            for g in self.gamma.iter_mut() {
                *g *= 0.5;
            }
            for h in self.histogram.iter_mut() {
                *h = 0;
            }
            for e in self.total_energy.iter_mut() {
                *e = Energy::new(0.);
            }
            println!(
                "[{:14}] WL:  We have found all {} bins, {:.3} <= gamma <= {:.3}",
                crate::prettyfloat::PrettyFloat(self.moves as f64),
                self.histogram.len(),
                crate::prettyfloat::PrettyFloat(min_of(&self.gamma)),
                crate::prettyfloat::PrettyFloat(max_of(&self.gamma))
            );
        }
    }
}

impl<S: MovableSystem> MonteCarlo for EnergyMC<S> {
    type Params = EnergyMCParams;
    type System = S;
    fn from_params(params: EnergyMCParams, mut system: S, save_as: ::std::path::PathBuf) -> Self {
        // center zero energy in a bin!
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
            min_allowed_energy: params.min_allowed_energy,
            max_allowed_energy: params.max_allowed_energy,

            energies: vec![
                system.energy(),
                system.energy() - min_T,
                system.energy() - 2. * min_T,
            ],
            gamma: vec![0.25; 3],
            histogram: vec![0; 4],
            lnw: vec![0., (1.-f).ln(), (1.-f).ln()+ f.ln(), (1.-f).ln()+ f.ln() + (1. / f - 1.).ln()],
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
            if let Some(mine) = self.min_allowed_energy {
                out_of_bounds = out_of_bounds || (e2 < mine && e2 < e1)
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
        for i in 0..mc.energies.len() {
            if sys.energy() > mc.energies[i] && (i == 0 || sys.energy() < mc.energies[i-1]) {
                print!(">>>");
            } else {
                print!("   ");
            }
            println!(" {:8.3}: {}", mc.energies[i].pretty(), mc.histogram[i]);
        }
        println!(
            "    {:8.3}: {}  [currently {:.3} after {:.2} moves]",
            "",
            mc.histogram[mc.histogram.len() - 1],
            sys.energy().pretty(),
            crate::prettyfloat::PrettyFloat(mc.moves as f64)
        );
        println!(
            " [biggest gamma = {:.2}]",
            crate::prettyfloat::PrettyFloat(max_of(&mc.gamma))
        );
    }
}
