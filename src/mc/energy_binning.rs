//! An implementation of SAD (Statistical Association, Dynamical version).
//!
//! This version using a [Binning]() approach so that we can use
//! something nicer than a histogram.  The idea is also to use this
//! for a 2D histogram, which will be even fancier.

#![allow(non_snake_case)]

use super::*;
use crate::system::*;

use super::plugin::Plugin;
use dimensioned::Dimensionless;
use rand::{Rng, SeedableRng};
use std::default::Default;

use binning::Binning;

/// Parameters to configure a particular MC.
#[derive(Debug, AutoArgs)]
#[allow(non_camel_case_types)]
pub enum MethodParams {
    /// Sad
    Sad {
        /// Use the SAD algorithm, with the specified minimum temperature of interest.
        min_T: Energy,
    },
    /// Samc
    Samc {
        /// Use the SAMC algorithm, with the specified t0 parameter
        t0: f64,
    },
    /// Use the Wang-Landau algorithm
    WL {
        /// After gamma drops to this value, begin a "production" run
        min_gamma: Option<f64>,
    },
    /// Use the 1/t-Wang-Landau algorithm
    inv_t_WL,
}

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
    /// The actual method.
    pub _method: MethodParams,
    /// The seed for the random number generator.
    pub seed: Option<u64>,
    /// The lowest energy to allow.
    min_allowed_energy: Option<Energy>,
    /// The highest energy to allow.
    max_allowed_energy: Option<Energy>,
    /// The width of an optional high-resolution histogram
    pub high_resolution_de: Option<Energy>,
    _binning: binning::BinningParams,
    _moves: MoveParams,
    _report: plugin::ReportParams,
    _movies: plugin::MovieParams,
    _save: plugin::SaveParams,
}

impl Default for EnergyMCParams {
    fn default() -> Self {
        EnergyMCParams {
            _method: MethodParams::Sad {
                min_T: 0.2 * units::EPSILON,
            },
            seed: None,
            min_allowed_energy: None,
            max_allowed_energy: None,
            high_resolution_de: None,
            _moves: MoveParams::TranslationScale(0.05 * units::SIGMA),
            _report: plugin::ReportParams::default(),
            _save: plugin::SaveParams::default(),
            _movies: plugin::MovieParams::default(),
            _binning: binning::BinningParams::default(),
        }
    }
}

/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct EnergyMC<S> {
    /// The system we are simulating.
    pub system: S,
    /// The method we use
    method: Method,
    /// The number of moves that have been made.
    pub moves: u64,
    /// The last move where we discovered a new energy.
    pub time_L: u64,
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

    /// The parameters describing the bins
    pub bins: binning::Bins,
    /// The width of an optional high-resolution histogram
    pub high_resolution: Option<binning::histogram::Bins>,
}

#[derive(Serialize, Deserialize, Debug)]
enum Method {
    /// Sad
    Sad {
        num_states: usize,
        min_T: Energy,
        too_lo: Energy,
        too_hi: Energy,
        tL: u64,
        tF: f64,
        latest_parameter: f64,
    },
    /// Samc
    Samc { t0: f64 },
    /// Wang Landau
    WL {
        gamma: f64,
        inv_t: bool,
        min_gamma: Option<f64>,
    },
}

impl Method {
    fn new(p: MethodParams, E: Energy) -> Self {
        match p {
            MethodParams::Sad { min_T } => Method::Sad {
                num_states: 0,
                min_T,
                too_lo: E,
                too_hi: E,
                tL: 0,
                tF: 0.,
                latest_parameter: 0.,
            },
            MethodParams::Samc { t0 } => Method::Samc { t0 },
            MethodParams::WL { min_gamma } => Method::WL {
                gamma: 1.0,
                inv_t: false,
                min_gamma,
            },
            MethodParams::inv_t_WL => Method::WL {
                gamma: 1.0,
                inv_t: true,
                min_gamma: None,
            },
        }
    }
    fn report<S: MovableSystem>(&self, mc: &EnergyMC<S>) {
        print!("    ");
        match self {
            Method::Sad { too_lo, too_hi, .. } => {
                print!(
                    "SAD: {:.5} ({:.3}) -> {:.5} ({:.3})",
                    too_lo.pretty(),
                    mc.bins.get_count(*too_lo).pretty(),
                    too_hi.pretty(),
                    mc.bins.get_count(*too_hi).pretty(),
                );
            }
            Method::Samc { .. } => {}
            Method::WL { .. } => {
                let hist = "hist".into();
                if mc.min_allowed_energy.is_some()
                    && mc.bins.get_count(mc.min_allowed_energy.unwrap()) == PerEnergy::new(0.)
                {
                    println!(
                        "    WL:  We only got down to {} > {}",
                        mc.bins.min_energy().pretty(),
                        mc.min_allowed_energy.unwrap().pretty()
                    );
                }
                if mc.max_allowed_energy.is_some()
                    && mc.bins.get_count(mc.max_allowed_energy.unwrap()) == PerEnergy::new(0.)
                {
                    println!(
                        "    WL:  We only got up to {} < {}",
                        mc.bins.max_energy().pretty(),
                        mc.max_allowed_energy.unwrap().pretty()
                    );
                }
                print!(
                    "WL:  We have reached flatness {:.2} min: E={}!",
                    (mc.bins.min_count_extra(hist) / mc.bins.mean_count_extra(hist)).pretty(),
                    mc.bins.min_count_extra_energy(hist).pretty()
                );
            }
        }
        println!(
            " [gamma = {:.2}]",
            crate::prettyfloat::PrettyFloat(mc.gamma())
        );
    }
    // fn entropy(&self, bins: &impl Binning, energies: &[Energy]) -> Vec<f64> {
    //     let mut entropy: Vec<_> = energies.iter().map(|&e| bins.get_lnw(e)).collect();
    //     if let Method::Sad {
    //         min_T,
    //         too_lo,
    //         too_hi,
    //         ..
    //     } = *self
    //     {
    //         let mut meanhist = PerEnergy::new(0.0);
    //         let mut meancount = 0.0;
    //         let too_lo_entropy = bins.get_lnw(too_lo);
    //         let too_hi_entropy = bins.get_lnw(too_hi);
    //         for e in energies.iter().cloned() {
    //             if e >= too_lo && e <= too_hi {
    //                 meanhist += bins.get_count(e);
    //                 meancount += 1.;
    //             }
    //         }
    //         meanhist /= meancount;
    //         for (i, s) in entropy.iter_mut().enumerate() {
    //             let e = energies[i];
    //             if bins.get_count(e) == PerEnergy::new(0.) {
    //                 *s = 0.0;
    //             } else if e < too_lo {
    //                 *s = (bins.get_count(e) / meanhist).ln()
    //                     + too_lo_entropy
    //                     + *((e - too_lo) / min_T).value();
    //             } else if e > too_hi {
    //                 *s = (bins.get_count(e) / meanhist).ln() + too_hi_entropy;
    //             } else {
    //                 *s = bins.get_lnw(e);
    //             }
    //         }
    //     } else if let Method::WL { gamma, .. } = self {
    //         if *gamma == 0.0 {
    //             // We are in a production run, so we should modify the
    //             // lnw based on the histogram.
    //             let hist = "hist".into();
    //             if bins.min_count_extra(hist) > PerEnergy::new(0.) {
    //                 // We have explored everything at least once, so
    //                 // we can at least take a log of everything...
    //                 for (i, s) in entropy.iter_mut().enumerate() {
    //                     let e = energies[i];
    //                     *s += (bins.count_extra(hist, e) * units::EPSILON).ln();
    //                 }
    //             }
    //         }
    //     }
    //     entropy
    // }
}

impl<S: MovableSystem + ConfirmSystem> EnergyMC<S> {
    /// This decides whether to reject the move based on the actual
    /// method in use.
    fn reject_move(&mut self, e1: Energy, e2: Energy) -> bool {
        let lnw1 = self.bins.get_lnw(e1);
        let lnw2 = self.bins.get_lnw(e2);
        match self.method {
            Method::Sad {
                too_lo,
                too_hi,
                min_T,
                ..
            } => {
                let lnw1 = if e1 < too_lo {
                    self.bins.get_lnw(too_lo) + *((e1 - too_lo) / min_T).value()
                } else if e1 > too_hi {
                    self.bins.get_lnw(too_hi)
                } else {
                    lnw1
                };
                let lnw2 = if e2 < too_lo {
                    self.bins.get_lnw(too_lo) + *((e2 - too_lo) / min_T).value()
                } else if e2 > too_hi {
                    self.bins.get_lnw(too_hi)
                } else {
                    lnw2
                };
                let rejected = lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp();
                if !rejected
                    && self.bins.get_count(e2) == PerEnergy::new(0.)
                    && e2 < too_hi
                    && e2 > too_lo
                {
                    // Here we do changes that need only happen when
                    // we encounter an energy in our important range
                    // that we have never seen before.
                    match &mut self.method {
                        Method::Sad { ref mut tL, .. } => {
                            *tL = self.moves;
                        }
                        _ => unreachable!(),
                    }
                }
                rejected
            }
            Method::Samc { .. } => lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp(),
            Method::WL { .. } => lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp(),
        }
    }
    /// This updates the lnw based on the actual method in use.
    fn update_weights(&mut self, energy: Energy) {
        let gamma = self.gamma(); // compute gamma out front...
        let old_highest_hist = self.bins.max_count();
        let old_hist_here = self.bins.get_count(energy);
        self.bins.increment_count(energy, gamma);
        if let Some(bins) = &mut self.high_resolution {
            bins.increment_count(energy, 0.);
        }
        let mut gamma_changed = false;
        let mut switch_to_samc: Option<f64> = None;
        match self.method {
            Method::Sad {
                min_T,
                ref mut too_lo,
                ref mut too_hi,
                ref mut tL,
                ref mut tF,
                ref mut num_states,
                ref mut latest_parameter,
                ..
            } => {
                let hist_here = self.bins.get_count(energy);
                if old_hist_here == PerEnergy::new(0.) {
                    let tfound = "t_found".into();
                    self.bins
                        .accumulate_extra(tfound, energy, self.moves as f64);
                }
                if hist_here > old_highest_hist {
                    if energy > *too_hi {
                        let lnw_too_hi = self.bins.get_lnw(*too_hi);
                        {
                            let too_hi = *too_hi;
                            self.bins.set_lnw(|e, count| {
                                if e > too_hi && count > PerEnergy::new(0.) {
                                    Some(lnw_too_hi)
                                } else {
                                    None
                                }
                            });
                        }
                        *latest_parameter = *((energy - *too_lo) / min_T).value();
                        *tL = self.moves;
                        // FIXME We probably should round the energy to one
                        // of the bins.  Add that to the Binning trait?
                        *too_hi = energy;
                        *num_states = self.bins.count_states(|e, count| {
                            e >= *too_lo && e <= *too_hi && count >= PerEnergy::new(0.)
                        })
                    } else if energy < *too_lo {
                        // We should set the lnw below energy too_low
                        // to something reasonable.
                        let lnw_too_lo = self.bins.get_lnw(*too_lo);
                        {
                            let too_lo = *too_lo;
                            self.bins.set_lnw(move |e, count| {
                                if e < too_lo && count > PerEnergy::new(0.) {
                                    Some(lnw_too_lo + *((e - too_lo) / min_T).value())
                                } else {
                                    None
                                }
                            });
                        }
                        *latest_parameter = *((*too_hi - energy) / min_T).value();
                        *tL = self.moves;
                        // FIXME We probably should round the energy to one
                        // of the bins.  Add that to the Binning trait?
                        *too_lo = energy;
                        *num_states = self.bins.count_states(|e, count| {
                            e >= *too_lo && e <= *too_hi && count >= PerEnergy::new(0.)
                        })
                    }
                }
                if *tL == self.moves {
                    gamma_changed = true;
                    // Set tF to the latest discovery time in the
                    // range of energies that we actually care about.
                    let old_tF = *tF;
                    let tfound = "t_found".into();
                    *tF = self.bins.max_total_extra(tfound);
                    if old_tF == *tF {
                        // We didn't change gamma after all!
                        gamma_changed = false;
                    } else {
                        // We just discovered a new important energy.
                        // Let's take this as an opportunity to revise our
                        // translation scale, and also to log the news.
                        if let MoveParams::AcceptanceRate(r) = self.move_plan {
                            let s = self.acceptance_rate / r;
                            let s = if s < 0.8 {
                                0.8
                            } else if s > 1.2 {
                                1.2
                            } else {
                                s
                            };
                            self.translation_scale *= s;
                            if !self.report.quiet {
                                println!(
                                    "        new translation scale: {:.3}",
                                    self.translation_scale
                                );
                                println!(
                                    "        acceptance rate {:.1}% [long-term: {:.1}%]",
                                    100.0 * self.acceptance_rate,
                                    100.0 * self.accepted_moves as f64 / self.moves as f64
                                );
                            }
                        }
                        if !self.report.quiet {
                            println!(
                                "    sad: [{}]  {}:  {:.7} < {:.7} ... {:.7} < {:.7}",
                                *tF,
                                self.bins.num_states(),
                                self.bins.min_energy().pretty(),
                                too_lo.pretty(),
                                too_hi.pretty(),
                                self.bins.max_energy().pretty()
                            );
                        }
                    }
                }
            }
            Method::Samc { .. } => {}
            Method::WL {
                ref mut gamma,
                inv_t,
                min_gamma,
            } => {
                let hist = "hist".into();
                let old_lowest_hist = self.bins.min_count_extra(hist);
                self.bins.accumulate_extra(hist, energy, 0.0);
                if let Some(min_gamma) = min_gamma {
                    if *gamma < min_gamma {
                        // We are in a production run.
                        return;
                    }
                }
                let lowest_hist = self.bins.min_count_extra(hist);

                if lowest_hist > old_lowest_hist
                    && (self.min_allowed_energy.is_none()
                        || self.bins.get_count(self.min_allowed_energy.unwrap())
                            > PerEnergy::new(0.))
                    && (self.max_allowed_energy.is_none()
                        || self.bins.get_count(self.max_allowed_energy.unwrap())
                            > PerEnergy::new(0.))
                {
                    if (inv_t && lowest_hist > PerEnergy::new(0.))
                        || lowest_hist >= 0.8 * self.bins.mean_count_extra(hist)
                    {
                        gamma_changed = true;
                        *gamma *= 0.5;
                        let flatness = lowest_hist / self.bins.mean_count_extra(hist);
                        self.bins.zero_out_extra(hist);
                        if let Some(min_gamma) = min_gamma {
                            if *gamma < min_gamma {
                                // Switch to a "production" run.
                                *gamma = 0.0;
                            }
                        }
                        println!(
                            "    WL:  We have reached flatness {:.2}! gamma = {}",
                            flatness.pretty(),
                            *gamma
                        );
                    }
                    if inv_t && *gamma < (self.bins.num_states() as f64) / (self.moves as f64) {
                        println!("    1/t-WL:  Switching to 1/t!");
                        switch_to_samc = Some(self.bins.num_states() as f64);
                    }
                }
            }
        }
        if let Some(t0) = switch_to_samc {
            self.method = Method::Samc { t0 };
        }
        if gamma_changed {
            // self.movies.new_gamma(self.moves, gamma);
            // self.movies.new_gamma(self.moves, self.gamma());
        }
    }
}

impl<S: MovableSystem + ConfirmSystem> EnergyMC<S> {
    fn gamma(&self) -> f64 {
        match self.method {
            Method::Sad {
                tF,
                latest_parameter,
                num_states,
                ..
            } => {
                let num_states = num_states as f64;
                if latest_parameter * tF as f64 * num_states == 0.0 {
                    0.0
                } else {
                    let t = self.moves as f64;
                    (latest_parameter + t / tF) / (latest_parameter + t / num_states * (t / tF))
                }
            }
            Method::Samc { t0 } => {
                let t = self.moves as f64;
                if t > t0 {
                    t0 / t
                } else {
                    1.0
                }
            }
            Method::WL { gamma, .. } => gamma,
        }
    }
}

impl<S: MovableSystem + serde::Serialize + serde::de::DeserializeOwned> MonteCarlo for EnergyMC<S> {
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
            method: Method::new(params._method, system.energy()),
            moves: 0,
            time_L: 0,
            accepted_moves: 0,
            acceptance_rate: 0.5, // arbitrary starting guess.
            min_allowed_energy: params.min_allowed_energy,
            max_allowed_energy: params.max_allowed_energy,

            bins: binning::Bins::from_params(&system, params._binning),
            high_resolution: params
                .high_resolution_de
                .map(|de| binning::histogram::Bins::new(system.energy(), de)),

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
        self.bins
            .accumulate_extra("energy".into(), e1, e1.value_unsafe);
        for (k, d) in self.system.data_to_collect(self.moves).into_iter() {
            self.bins.accumulate_extra(k, e1, d);
        }
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
impl<S: MovableSystem + serde::Serialize + serde::de::DeserializeOwned> Plugin<EnergyMC<S>> for Logger {
    fn log(&self, mc: &EnergyMC<S>, _sys: &S) {
        mc.method.report(&mc);
    }
}
