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
    /// The current gamma value
    pub gamma: f64,
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
        if self.energies.len() < 2 {
            // No way to shift things just yet, as we do not have a length scale.
        } else if i == 0 {
            // shift the right side if we are on far left
            let de = self.energies[i + 1] - self.energies[i];
            self.energies[i] += de * self.gamma;
        } else if i == self.histogram.len() - 1 {
            // shift the left side if we are on the far right
            let de = self.energies[i - 1] - self.energies[i - 2];
            self.energies[i - 1] += de * self.gamma;
            // Now decide whether to add a new bin
            let elow = self.total_energy[i] / self.histogram[i] as f64;
            if self.energies[i - 1] - elow > self.min_T {
                // We have room for another bin (we think)
                self.energies.push(elow); // FIXME
                self.histogram.pop();
                self.histogram.push(0);
                self.histogram.push(0);
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
                self.total_energy[self.histogram.len() - 1] = Energy::new(0.);
                self.total_energy.push(Energy::new(0.));
            }
        } else if i < self.energies.len() {
            // shift both sides of the bin
            let de_plus = self.energies[i + 1] - self.energies[i];
            let de_minus = self.energies[i - 1] - self.energies[i - 2];
            let de_self = self.energies[i] - self.energies[i - 1];
            // Left side
            let de = if de_self < de_minus {
                de_self
            } else {
                de_minus
            };
            self.energies[i - 1] += de * self.gamma;
            // Right side
            let de = if de_self < de_plus { de_self } else { de_plus };
            self.energies[i] += de * self.gamma;
        }
        if let Some(min_gamma) = self.min_gamma {
            if self.gamma < min_gamma {
                // We are in a production run.
                return;
            }
        }
        if old_hist_here == 0 && !self.histogram.iter().any(|&x| x == 0) {
            // We have explored all the bins!
            self.gamma *= 0.5;
            for h in self.histogram.iter_mut() {
                *h = 0;
            }
            for e in self.total_energy.iter_mut() {
                *e = Energy::new(0.);
            }
            println!("    WL:  We have found everything, gamma = {}", self.gamma);
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
        EnergyMC {
            gamma: 0.5,
            min_gamma: params.min_gamma,
            moves: 0,
            min_T: params.min_T,
            f: params.f,
            accepted_moves: 0,
            acceptance_rate: 0.5, // arbitrary starting guess.
            min_allowed_energy: params.min_allowed_energy,
            max_allowed_energy: params.max_allowed_energy,

            energies: vec![system.energy()],
            histogram: vec![0, 0],
            lnw: vec![0., 0.5f64.ln()],
            total_energy: vec![Energy::new(0.0), Energy::new(0.0)],

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
    fn log(&self, mc: &EnergyMC<S>, _sys: &S) {
        print!("    ");
        println!(
            " [gamma = {:.2}]",
            crate::prettyfloat::PrettyFloat(mc.gamma)
        );
    }
}
