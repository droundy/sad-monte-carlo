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

use dimensioned::{Abs, Dimensionless};

/// Parameters to configure the moves.
#[derive(Serialize, Deserialize, Debug, AutoArgs, Clone)]
pub enum MoveParams {
    /// The rms distance of moves
    TranslationScale(Length),
    /// Adjust translation scale to reach acceptance rate.
    AcceptanceRate(f64),
}

/// The parameters needed to configure a simulation.
#[derive(Debug, AutoArgs, Clone)]
pub struct EnergyMCParams {
    /// The rejection fraction.
    pub f: f64,
    /// The lowest temperature of interest
    min_T: Energy,
    /// After gamma drops to this value, begin a "production" run
    min_gamma: Option<f64>,
    /// The seed for the random number generator.
    pub seed: Option<u64>,
    _moves: MoveParams,
    /// report input
    pub _report: plugin::ReportParams,
    /// report input
    pub _movies: plugin::MovieParams,
    /// report input
    pub _save: plugin::SaveParams,
}

impl Default for EnergyMCParams {
    fn default() -> Self {
        EnergyMCParams {
            min_gamma: None,
            min_T: 0.2 * units::EPSILON,
            f: 0.5,
            seed: None,
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
    /// The histogram counts since added a bin
    pub histogram_since_adding_bin: Vec<u64>,
    /// Whether we have visited here since decreasing gamma
    pub have_seen: Vec<bool>,
    /// The total energy in each bin
    pub total_energy: Vec<Energy>,
    /// The total energy in our average of the min_energy
    pub min_energy_total: Energy,
    /// The count of the min_energy average
    pub min_energy_count: u64,
    /// The recent mean of the acceptance probability.  This is a
    /// logarithmically waverage.
    pub recent_acceptance_probability: f64,
}

const RECENT_ACCEPTANCE: f64 = 1.0 / 16.0;

impl<S: MovableSystem + ConfirmSystem> EnergyMC<S> {
    /// This decides whether to reject the move based on the actual
    /// method in use.
    fn reject_move(&mut self, e1: Energy, e2: Energy) -> bool {
        let lnw1 = self.lnw[self.e_to_idx(e1)];
        let lnw2 = self.lnw[self.e_to_idx(e2)];
        if self.rel_bins[self.rel_bins.len() - 1] < 1e-10 * self.bin_norm {
            println!(
                "  We are in crazy land {} -> {}:\n  {:10} -> {:10}\n  {:10} -> {:10}\n  {:10} -> {:10}",
                crate::prettyfloat::PrettyFloat(self.rel_bins[self.rel_bins.len() - 1] / self.bin_norm),
                crate::prettyfloat::PrettyFloat((lnw1 - lnw2).exp()),
                e1.pretty(),
                e2.pretty(),
                lnw1,
                lnw2,
                self.e_to_idx(e1),
                self.e_to_idx(e2)
            );
            assert!(false);
        }
        if lnw2 > lnw1 {
            let prob = (lnw1 - lnw2).exp();
            self.recent_acceptance_probability = self.recent_acceptance_probability
                * (1.0 - RECENT_ACCEPTANCE)
                + prob * RECENT_ACCEPTANCE;
            self.rng.gen::<f64>() > prob
        } else {
            false
        }
    }
    fn e_to_idx(&self, energy: Energy) -> usize {
        let mut e = self.min_energy;
        let dedw = (self.max_energy - self.min_energy) / self.bin_norm;
        for (i, w) in self.rel_bins.iter().cloned().enumerate().rev() {
            if energy < e {
                return i + 2;
            }
            e += dedw * w;
        }
        if energy < self.max_energy {
            return 1;
        }
        0
    }
    fn e_to_idx_and_boundaries(&self, energy: Energy) -> (usize, Option<Energy>, Option<Energy>) {
        let mut e = self.min_energy;
        let mut lower_bound = None;
        let dedw = (self.max_energy - self.min_energy) / self.bin_norm;
        for (i, w) in self.rel_bins.iter().cloned().enumerate().rev() {
            if energy < e {
                return (i + 2, lower_bound, Some(e));
            }
            lower_bound = Some(e);
            e += dedw * w;
        }
        if energy < self.max_energy {
            return (1, lower_bound, Some(self.max_energy));
        }
        (0, Some(self.max_energy), None)
    }
    /// Find out the lnw for a given energy
    pub fn e_to_lnw(&self, energy: Energy) -> f64 {
        self.lnw[self.e_to_idx(energy)]
    }
    /// This updates the energy boundaries, not the weights.
    fn update_weights(&mut self, energy: Energy) {
        let (i, elo, ehi) = self.e_to_idx_and_boundaries(energy);
        self.histogram[i] += 1;
        self.histogram_since_adding_bin[i] += 1;
        self.total_energy[i] += energy;
        let mean_here = self.total_energy[i] / self.histogram[i] as f64;
        // If the mean is not within the bin at all, we just reset the histogram
        // and total energy for this bin.
        if let Some(elo) = elo {
            if mean_here < elo {
                self.histogram[i] = 1;
                self.total_energy[i] = energy;
            }
        }
        if let Some(ehi) = ehi {
            if mean_here > ehi {
                self.histogram[i] = 1;
                self.total_energy[i] = energy;
            }
        }
        if i == 0 {
            // We are in the widthless arbitrarily high-energy bin.
            let w = self.rel_bins[0];
            let old_breadth = self.max_energy - self.min_energy;
            let wid = w * old_breadth / self.bin_norm;
            let de = self.gamma * wid;
            self.max_energy += de;
            // We expand just the highest energy bounded bin.
            self.bin_norm *= *((self.max_energy - self.min_energy) / old_breadth).value();
            self.rel_bins[0] =
                *((wid + de) / (self.max_energy - self.min_energy)).value() * self.bin_norm;
        } else if i == self.rel_bins.len() + 1 {
            // We are in the low-energy bin.
            self.min_energy_count += 1;
            self.min_energy_total += self.min_energy;
            let min_energy_mean = self.min_energy_total / self.min_energy_count as f64;
            if (min_energy_mean - self.min_energy).abs() > self.min_T {
                // The minimum energy has shifted a fair amount...
                self.histogram_since_adding_bin
                    .iter_mut()
                    .map(|x| *x = 0)
                    .count(); // reset histogram
                self.min_energy_count = 1;
                self.min_energy_total = self.min_energy;
            }

            let w = self.rel_bins[i - 2];
            let old_breadth = self.max_energy - self.min_energy;
            let wid = w * old_breadth / self.bin_norm;
            let de = self.gamma * wid;
            self.min_energy -= de;
            // We expand just the highest energy bounded bin.
            self.bin_norm *= *((self.max_energy - self.min_energy) / old_breadth).value();
            self.rel_bins[i - 2] =
                *((wid + de) / (self.max_energy - self.min_energy)).value() * self.bin_norm;

            // Now decide whether we want a new bin here.
            let min_de: Energy = -self.min_T * self.f.ln();
            assert!(min_de > Energy::new(0.));
            let min_w = *(self.bin_norm * min_de / (self.max_energy - self.min_energy)).value();
            if min_w < 0. {
                println!(
                    " maximum energy {} vs minimum energy {}",
                    self.max_energy, self.min_energy
                );
                println!("this gives min_de of {}", min_de);
                println!(
                    "this gives bin_norm of {} compared with {}",
                    self.bin_norm,
                    self.rel_bins.iter().cloned().sum::<f64>()
                );
                println!("this gives min_w of {}", min_w);
                assert!(false);
            }
            if self.recent_acceptance_probability < self.gamma {
                println!("\nWe are in crazy land already! Recent acceptance probability = {:.3}, while gamma = {:.3}, E = {}",
                    crate::prettyfloat::PrettyFloat(self.recent_acceptance_probability),
                    crate::prettyfloat::PrettyFloat(self.gamma),
                    energy.pretty(),
                );
                // assert!(false);
            }

            let de = min_energy_mean - mean_here;
            if de > self.min_T
                && self.recent_acceptance_probability > self.gamma.sqrt()
                && self
                    .histogram_since_adding_bin
                    .iter()
                    .cloned()
                    .min()
                    .unwrap()
                    > 0
                && self.rel_bins[self.rel_bins.len() - 1] > min_w
            {
                // We add a new low energy bin under the following circumstances:
                //
                // 1. The mean energy in the unbounded low-energy bin must be at least min_T
                //    below the upper value of that bin.  This means that the temperature in this
                //    bin must be greater than the minimum temperature we are aiming for.
                //
                // 2. We must have visited all bins since the last time we added a new bin.
                //    This is important because otherwise we could end up adding way too many
                //    bins if we start out at some insanely high energy.
                //
                // 3. The smallest energy bin must not be must be smaller than min_w, which is itself
                //    proportional to min_T.

                // Let us actually add several bins each time we realize we need more.
                // This should allow us to more rapidly explore lower energies.
                let my_relw = -*(self.bin_norm * de / (self.max_energy - self.min_energy)).value()
                    * self.f.ln();
                if true {
                    // FIXME Set the above to false for less verbosity, when I finish debugging.
                    println!("\nAdding new bin with width {}", de.pretty());
                    self.log_me();
                    println!("");
                }
                // remove the unbound bin
                self.lnw.pop();
                // for _ in 0..10 {
                self.rel_bins.push(my_relw);
                self.bin_norm += my_relw;
                self.histogram.push(0);
                self.total_energy.push(Energy::new(0.));
                self.have_seen.push(false);
                // add a regular bin
                self.lnw.push(self.lnw.last().unwrap() + self.f.ln());
                // }
                // add the extra-big last bin.
                self.lnw
                    .push(self.lnw.last().unwrap() + (self.f / (1. - self.f)).ln());
                for h in self.have_seen.iter_mut() {
                    *h = false;
                }
                self.histogram_since_adding_bin = vec![0; self.histogram.len()];
                let newgamma = 4.0 * self.gamma;
                if newgamma <= GAMMA_INIT {
                    self.gamma = newgamma; // make gamma a bit bigger to let us explore this new region of energy.
                }
            }
        } else if i == self.rel_bins.len() {
            let w = self.rel_bins[i - 1];
            let wid = self.rel_bins[i - 1] * (self.max_energy - self.min_energy) / self.bin_norm;
            let dw = 2.0 * self.gamma * w;
            self.rel_bins[i - 1] -= dw;
            assert!(self.rel_bins[i - 1] > 0.0);
            self.bin_norm -= dw;
            self.min_energy += self.gamma * wid;
        } else if i == 1 {
            let w = self.rel_bins[i - 1];
            let wid = self.rel_bins[i - 1] * (self.max_energy - self.min_energy) / self.bin_norm;
            let dw = 2.0 * self.gamma * w;
            self.rel_bins[i - 1] -= dw;
            assert!(self.rel_bins[i - 1] > 0.0);
            self.bin_norm -= dw;
            self.max_energy -= self.gamma * wid;
        } else {
            let w = self.rel_bins[i - 1];
            let dw = self.gamma * w;
            self.rel_bins[i - 1] -= dw;
            assert!(self.rel_bins[i - 1] > 0.0);
            self.bin_norm -= dw;
        }
        if self.bin_norm <= 1e-8 {
            self.bin_norm = self.rel_bins.iter().cloned().sum::<f64>();
            // Let's re-normalize the bin widths, since we are getting a bit small.
            for w in self.rel_bins.iter_mut() {
                *w /= self.bin_norm;
            }
            // The following should be about 1, but due to roundoff error it is
            // not, and we don't want any (unavoidable) errors in our binning.
            self.bin_norm = self.rel_bins.iter().cloned().sum::<f64>();
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
                    // We are in the 1/t stage.
                    self.gamma = self.rel_bins.len() as f64 / self.moves as f64;
                }
            }
        }
    }
}

const MAX_INIT: usize = 10; // Looking for 10 initial quantiles requires maybe a dozen megabytes of memory.
const GAMMA_INIT: f64 = 0.25;

impl<S: MovableSystem + serde::Serialize + serde::de::DeserializeOwned> MonteCarlo for EnergyMC<S> {
    type Params = EnergyMCParams;
    type System = S;
    fn from_params(params: EnergyMCParams, mut system: S, save_as: ::std::path::PathBuf) -> Self {
        let mut rng = crate::rng::MyRng::seed_from_u64(params.seed.unwrap_or(0));
        let f = params.f;
        let min_T = params.min_T;

        // First let's brute-force to find where a number of the quantiles are
        // We look for at most MAX_INIT quantiles.
        let mut energies = Vec::with_capacity(1 << (2 * MAX_INIT));
        for _ in 0..1 << (2 * MAX_INIT) {
            energies.push(system.randomize(&mut rng));
        }
        energies.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut energy_boundaries = Vec::with_capacity(MAX_INIT);
        energy_boundaries.push(energies[1 << MAX_INIT]);
        for i in (0..MAX_INIT).rev() {
            let here = energies[1 << i];
            if *energy_boundaries.last().unwrap() > here + params.min_T
                || energy_boundaries.len() < 3
            {
                energy_boundaries.push(here);
            } else {
                break;
            }
        }
        std::mem::drop(energies); // Just to save on RAM...
        let max_energy = energy_boundaries[0];
        let min_energy = energy_boundaries[energy_boundaries.len() - 1];
        let mut lnw = Vec::with_capacity(energy_boundaries.len() + 2);
        lnw.push((1.0 - f).ln());
        for _ in 0..energy_boundaries.len() - 1 {
            lnw.push(lnw[lnw.len() - 1] + f.ln());
        }
        lnw.push(lnw[lnw.len() - 1] + (f / (1. - f)).ln());

        let mut rel_bins = Vec::with_capacity(lnw.len());
        for i in 0..energy_boundaries.len() - 1 {
            rel_bins.push(
                *((energy_boundaries[i] - energy_boundaries[i + 1]) / (max_energy - min_energy))
                    .value(),
            );
        }

        assert_eq!(rel_bins.len() + 2, lnw.len());
        EnergyMC {
            min_gamma: params.min_gamma,
            moves: 0,
            min_T,
            f,
            accepted_moves: 0,
            acceptance_rate: 0.5, // arbitrary starting guess.
            max_energy,
            min_energy,

            bin_norm: rel_bins.iter().cloned().sum::<f64>(),
            rel_bins,
            gamma: GAMMA_INIT,
            histogram: vec![0; lnw.len()],
            histogram_since_adding_bin: vec![0; lnw.len()],
            have_seen: vec![false; lnw.len()],
            total_energy: vec![Energy::new(0.0); lnw.len()],
            lnw,

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
            min_energy_count: 1,
            min_energy_total: min_energy,
            recent_acceptance_probability: 1.0,
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
            if !self.reject_move(e1, e2) {
                self.accepted_moves += 1;
                self.acceptance_rate += recent_scale;
                self.system.confirm();
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

impl<S: MovableSystem> EnergyMC<S> {
    fn log_me(&self) {
        println!(
            "        E_hi  = {:10.5} (mean {:8.3} from {:.3}): {:.3}",
            self.max_energy.pretty(),
            (self.total_energy[0] / self.histogram[0] as f64).pretty(),
            crate::prettyfloat::PrettyFloat(self.histogram[0] as f64),
            crate::prettyfloat::PrettyFloat(self.histogram_since_adding_bin[0] as f64),
        );
        let i_current = self.e_to_idx(self.system.energy());
        println!(
            "        E_now = {:10.5} (mean {:8.3} from {:.3}): {:.3} i = {}",
            self.system.energy().pretty(),
            (self.total_energy[i_current] / self.histogram[i_current] as f64).pretty(),
            crate::prettyfloat::PrettyFloat(self.histogram[i_current] as f64),
            crate::prettyfloat::PrettyFloat(self.histogram_since_adding_bin[i_current] as f64),
            i_current,
        );
        let de_po = (self.max_energy - self.min_energy) * self.rel_bins[self.rel_bins.len() - 1]
            / self.bin_norm;
        println!(
            "  lo Delta E  = {:10.5} (mean {:8.3} from {:.3}): {:.3}",
            de_po.pretty(),
            (self.total_energy[self.histogram.len() - 2]
                / self.histogram[self.histogram.len() - 2] as f64)
                .pretty(),
            crate::prettyfloat::PrettyFloat(self.histogram[self.histogram.len() - 2] as f64),
            crate::prettyfloat::PrettyFloat(
                self.histogram_since_adding_bin[self.histogram.len() - 2] as f64
            ),
        );
        println!(
            "        E_lo  = {:10.5} (mean {:8.3} from {:.3}): {:.3}",
            self.min_energy.pretty(),
            (self.total_energy[self.histogram.len() - 1]
                / self.histogram[self.histogram.len() - 1] as f64)
                .pretty(),
            crate::prettyfloat::PrettyFloat(self.histogram[self.histogram.len() - 1] as f64),
            crate::prettyfloat::PrettyFloat(
                self.histogram_since_adding_bin[self.histogram.len() - 1] as f64
            ),
        );
        println!(
            "        {} [n_bins = {}, gamma = {:.2}, P_recent = {:.2}]",
            self.system.describe(),
            self.rel_bins.len() + 2,
            crate::prettyfloat::PrettyFloat(self.gamma),
            crate::prettyfloat::PrettyFloat(self.recent_acceptance_probability),
        );
    }
}

#[derive(Serialize, Deserialize, Debug)]
struct Logger;
impl<S: MovableSystem + serde::Serialize + serde::de::DeserializeOwned> Plugin<EnergyMC<S>>
    for Logger
{
    fn log(&self, mc: &EnergyMC<S>, _sys: &S) {
        mc.log_me();
    }
}
