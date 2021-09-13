//! An implementation of SAD (Statistical Association, Dynamical version).

#![allow(non_snake_case)]

use super::*;
use crate::system::*;

use super::plugin::Plugin;
use crate::prettyfloat::PrettyFloat;
use dimensioned::Dimensionless;
use rand::{Rng, SeedableRng};
use std::default::Default;

/// Which experimental version of SAD are we doing?
#[derive(Serialize, Deserialize, Debug, AutoArgs, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SadVersion {
    /// Sad
    Sad,
}

impl SadVersion {
    fn compute_gamma(&self, latest_parameter: f64, t: f64, tF: f64, num_states: f64) -> f64 {
        if latest_parameter * tF * num_states == 0.0 {
            0.0
        } else {
            match *self {
                SadVersion::Sad => {
                    let g = (latest_parameter + t / tF)
                        / (latest_parameter + t / num_states * (t / tF));
                    if g.is_nan() {
                        println!("problem with {} and {}", tF, num_states);
                        assert!(!g.is_nan());
                    }
                    g
                }
            }
        }
    }
}

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
    Inv_t_WL,
    /// Use the 1/t-Wang-Landau algorithm, nicer spelling
    inv_t_wl,
    /// A canonical simulation
    _Canonical {
        /// Use a canonical simulation with this temperature
        T: Energy,
    },
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
    /// The energy binsize.
    energy_bin: Option<Energy>,
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
            _method: MethodParams::Sad {
                min_T: 0.2 * units::EPSILON,
            },
            seed: None,
            min_allowed_energy: None,
            max_allowed_energy: None,
            energy_bin: None,
            _moves: MoveParams::TranslationScale(0.05 * units::SIGMA),
            _report: plugin::ReportParams::default(),
            _movies: plugin::MovieParams::default(),
            _save: plugin::SaveParams::default(),
        }
    }
}

/// This defines a "state".  In this case it is just an energy, but it
/// should make this code easier to transition to a grand ensemble.
#[derive(Serialize, Deserialize, Debug, Clone, Copy, PartialEq)]
pub struct State {
    /// The energy
    pub E: Energy,
}
impl ::std::fmt::Display for State {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "{}", PrettyFloat(*(self.E / units::EPSILON).value()))
    }
}
impl State {
    /// Find the state of a system
    pub fn new<S: System>(sys: &S) -> State {
        State { E: sys.energy() }
    }
}
/// A set of counts for a variable
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct BinCounts {
    /// The total of the thing
    total: Vec<f64>,
    /// The count in each bin
    count: Vec<u64>,
}

/// Where we store the info about the energy grid
#[derive(Serialize, Deserialize, Debug)]
pub struct Bins {
    /// The lowest allowed energy in any bin.
    pub min: Energy,
    /// The energy bin size.
    pub width: Energy,
    /// The number of times we have been at each energy.
    pub histogram: Vec<u64>,
    /// The iteration when we found each energy.
    pub t_found: Vec<u64>,
    /// The ln weight for each energy bin.
    pub lnw: Vec<Unitless>,
    /// The total of all energes found in each bin
    pub energy_total: Vec<Energy>,
    /// The total square of all energes found in each bin
    pub energy_squared_total: Vec<EnergySquared>,
    /// Extra data we might want to collect occasionally
    pub extra: std::collections::HashMap<Interned, BinCounts>,
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
    movies: plugin::Movie,
    save: plugin::Save,
    manager: plugin::PluginManager,

    // The following were formerly part of Bins.  I joined them all
    // together in EnergyMC.
    /// The parameters describing the bins
    pub bins: Bins,

    /// Whether we have seen this since the last visit to maxentropy.
    have_visited_since_maxentropy: Vec<bool>,
    /// How many round trips have we seen at this energy.
    round_trips: Vec<u64>,
    /// The maximum entropy we have seen.
    max_S: Unitless,
    /// The index with the maximum entropy.
    max_S_index: usize,
}

#[derive(Serialize, Deserialize, Debug)]
enum Method {
    /// Sad
    Sad {
        min_T: Energy,
        too_lo: Energy,
        too_hi: Energy,
        tL: u64,
        tF: u64,
        num_states: u64,
        highest_hist: u64,
        version: SadVersion,
        latest_parameter: f64,
    },
    /// Samc
    Samc { t0: f64 },
    /// Wang Landau
    WL {
        gamma: f64,
        lowest_hist: u64,
        highest_hist: u64,
        total_hist: u64,
        num_states: f64,
        hist: Vec<u64>,
        min_energy: Energy,
        inv_t: bool,
        min_gamma: Option<f64>,
    },
    /// Canonical
    Canonical {
        temperature: Energy,
    }
}

impl Method {
    fn new(
        p: MethodParams,
        E: Energy,
        dE: Energy,
        min_allowed_energy: Option<Energy>,
        max_allowed_energy: Option<Energy>,
    ) -> Self {
        match p {
            MethodParams::Sad { min_T } => Method::Sad {
                min_T,
                too_lo: E,
                too_hi: E,
                tL: 0,
                tF: 0,
                num_states: 1,
                highest_hist: 1,
                version: SadVersion::Sad,
                latest_parameter: 0.,
            },
            MethodParams::Samc { t0 } => Method::Samc { t0 },
            MethodParams::WL { min_gamma } => Method::WL {
                gamma: 1.0,
                lowest_hist: if min_allowed_energy.is_some() && max_allowed_energy.is_some() {
                    0
                } else {
                    1
                },
                highest_hist: 1,
                total_hist: 0,
                num_states: if let (Some(mine), Some(maxe)) =
                    (min_allowed_energy, max_allowed_energy)
                {
                    *((maxe - mine) / dE).value()
                } else {
                    1.0
                },
                hist: Vec::new(),
                min_energy: E,
                inv_t: false,
                min_gamma,
            },
            MethodParams::Inv_t_WL | MethodParams::inv_t_wl => Method::WL {
                gamma: 1.0,
                lowest_hist: if min_allowed_energy.is_some() && max_allowed_energy.is_some() {
                    0
                } else {
                    1
                },
                highest_hist: 1,
                total_hist: 0,
                num_states: if let (Some(mine), Some(maxe)) =
                    (min_allowed_energy, max_allowed_energy)
                {
                    *((maxe - mine) / dE).value()
                } else {
                    1.0
                },
                hist: Vec::new(),
                min_energy: E,
                inv_t: true,
                min_gamma: None,
            },
            MethodParams::_Canonical { T } => Method::Canonical {
                temperature: T,
            },
        }
    }
    // fn entropy(&self, bins: &Bins) -> Vec<f64> {
    //     let mut entropy = bins.lnw.clone();
    //     if let Method::Sad {
    //         min_T,
    //         too_lo,
    //         too_hi,
    //         num_states,
    //         ..
    //     } = *self
    //     {
    //         let mut meanhist = 0.0;
    //         let too_lo_entropy = bins.lnw[bins.state_to_index(State { E: too_lo })];
    //         let too_hi_entropy = bins.lnw[bins.state_to_index(State { E: too_hi })];
    //         for (i, h) in bins.histogram.iter().cloned().enumerate() {
    //             if bins.index_to_state(i).E >= too_lo && bins.index_to_state(i).E <= too_hi {
    //                 meanhist += h as f64;
    //             }
    //         }
    //         // Round too_lo to the center of the energy bin.
    //         let too_lo = bins
    //             .index_to_state(bins.state_to_index(State { E: too_lo }))
    //             .E;
    //         meanhist /= num_states as f64;
    //         for (i, s) in entropy.iter_mut().enumerate() {
    //             let e = bins.index_to_state(i).E;
    //             if bins.histogram[i] == 0 {
    //                 *s *= 0.0;
    //             } else if e < too_lo {
    //                 *s = (bins.histogram[i] as f64 / meanhist).ln()
    //                     + too_lo_entropy
    //                     + (e - too_lo) / min_T;
    //             } else if e > too_hi {
    //                 *s = (bins.histogram[i] as f64 / meanhist).ln() + too_hi_entropy;
    //             }
    //         }
    //     } else if let Method::WL { gamma, hist, .. } = self {
    //         if *gamma == 0.0 {
    //             // We are in a production run, so we should modify the
    //             // lnw based on the histogram.
    //             return entropy
    //                 .iter()
    //                 .map(|x| *x.value())
    //                 .zip(hist.iter().cloned())
    //                 .map(|(lnw, h)| lnw + if h > 0 { (h as f64).ln() } else { 0.0 })
    //                 .collect();
    //         }
    //     }
    //     entropy.iter().map(|x| *(x.value())).collect()
    // }
}

impl Bins {
    fn index_to_state(&self, i: usize) -> State {
        State {
            E: self.min + (i as f64 + 0.5) * self.width,
        }
    }
    fn state_to_index(&self, s: State) -> usize {
        *((s.E - self.min) / self.width).value() as usize
    }
    fn accumulate_extra(&mut self, k: Interned, idx: usize, value: f64) {
        if let Some(data) = self.extra.get_mut(&k) {
            data.count[idx] += 1;
            data.total[idx] += value;
        } else {
            let values = BinCounts {
                count: vec![0; self.lnw.len()],
                total: vec![0.0; self.lnw.len()],
            };
            self.extra.insert(k, values);
            self.accumulate_extra(k, idx, value); // sloppy recursion...
        }
    }
}

impl<S: System> EnergyMC<S> {
    /// Find the index corresponding to a given energy.  This should
    /// panic if the energy is less than `min`.
    pub fn state_to_index(&self, s: State) -> usize {
        self.bins.state_to_index(s)
    }
    /// Find the energy corresponding to a given index.
    pub fn index_to_state(&self, i: usize) -> State {
        self.bins.index_to_state(i)
    }
    /// Make room in our arrays for a new energy value
    pub fn prepare_for_state(&mut self, e: State) {
        let e = e.E;
        assert!(self.bins.width > Energy::new(0.0));
        while e < self.bins.min {
            // this is a little wasteful, but seems the easiest way to
            // ensure we end up with enough room.
            self.bins.histogram.insert(0, 0);
            self.bins.t_found.insert(0, 0);
            self.bins.lnw.insert(0, Unitless::new(0.0));
            self.bins.energy_total.insert(0, Energy::new(0.0));
            self.bins
                .energy_squared_total
                .insert(0, EnergySquared::new(0.0));
            for v in self.bins.extra.iter_mut() {
                v.1.count.insert(0, 0);
                v.1.total.insert(0, 0.0);
            }
            self.have_visited_since_maxentropy.insert(0, true);
            self.round_trips.insert(0, 1);
            self.bins.min -= self.bins.width;
        }
        while e >= self.bins.min + self.bins.width * (self.bins.lnw.len() as f64) {
            self.bins.lnw.push(Unitless::new(0.0));
            self.bins.histogram.push(0);
            self.bins.t_found.push(0);
            for v in self.bins.extra.iter_mut() {
                v.1.count.push(0);
                v.1.total.push(0.0);
            }
            self.bins.energy_total.push(Energy::new(0.0));
            self.bins.energy_squared_total.push(EnergySquared::new(0.0));
            self.have_visited_since_maxentropy.push(true);
            self.round_trips.push(1);
        }
    }
}

impl<S: System> EnergyMC<S> {
    /// This decides whether to reject the move based on the actual
    /// method in use.
    fn reject_move(&mut self, e1: State, e2: State) -> bool {
        let i1 = self.state_to_index(e1);
        let i2 = self.state_to_index(e2);
        match self.method {
            Method::Sad {
                too_lo,
                too_hi,
                min_T,
                ..
            } => {
                let lnw = &self.bins.lnw;
                let lnw1 = if e1.E < too_lo {
                    lnw[self.state_to_index(State { E: too_lo })] + (e1.E - too_lo) / min_T
                } else if e1.E > too_hi {
                    lnw[self.state_to_index(State { E: too_hi })]
                } else {
                    lnw[i1]
                };
                let lnw2 = if e2.E < too_lo {
                    lnw[self.state_to_index(State { E: too_lo })] + (e2.E - too_lo) / min_T
                } else if e2.E > too_hi {
                    lnw[self.state_to_index(State { E: too_hi })]
                } else {
                    lnw[i2]
                };
                let rejected = lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp();
                if !rejected && self.bins.histogram[i2] == 0 && e2.E < too_hi && e2.E > too_lo {
                    // Here we do changes that need only happen when
                    // we encounter an energy in our important range
                    // that we have never seen before.
                    match &mut self.method {
                        Method::Sad {
                            ref mut num_states,
                            ref mut tL,
                            ref mut latest_parameter,
                            ..
                        } => {
                            *latest_parameter = *((too_hi - too_lo) / min_T).value();
                            *num_states += 1;
                            *tL = self.moves;
                        }
                        _ => unreachable!(),
                    }
                }
                rejected
            }
            Method::Samc { .. } => {
                let lnw1 = self.bins.lnw[i1].value();
                let lnw2 = self.bins.lnw[i2].value();
                lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp()
            }
            Method::WL {
                ref mut num_states,
                lowest_hist,
                ..
            } => {
                let lnw1 = self.bins.lnw[i1].value();
                let lnw2 = self.bins.lnw[i2].value();
                let rejected = lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp();
                if !rejected && self.bins.histogram[i2] == 0 && lowest_hist > 0 {
                    *num_states += 1.0;
                }
                rejected
            }
            Method::Canonical { temperature } => {
                if e1.E >= e2.E {
                    false
                } else {
                    self.rng.gen::<f64>() > ((e1.E - e2.E)/temperature).exp()
                }
            }
        }
    }
    /// This updates the lnw based on the actual method in use.
    fn update_weights(&mut self, energy: State) {
        let i = self.state_to_index(energy);
        let gamma = self.gamma(); // compute gamma out front...
        let old_lnw = self.bins.lnw[i];
        self.bins.lnw[i] += gamma;
        // let mut gamma_changed = false;
        let mut switch_to_samc: Option<f64> = None;
        match self.method {
            Method::Canonical {..} => (), // Nothing to update!
            Method::Sad {
                min_T,
                ref mut too_lo,
                ref mut too_hi,
                ref mut num_states,
                ref mut tL,
                ref mut tF,
                ref mut highest_hist,
                ref mut latest_parameter,
                ..
            } => {
                let histogram = &self.bins.histogram;
                if *too_lo > *too_hi || energy.E < *too_lo || energy.E > *too_hi {
                    // Ooops, we didn't want to add gamma after all...
                    self.bins.lnw[i] = old_lnw;
                }

                if histogram[i] > *highest_hist {
                    *highest_hist = histogram[i];
                    if energy.E > *too_hi {
                        let ihi = self.bins.state_to_index(State { E: *too_hi });
                        for j in 0..histogram.len() {
                            let ej = self.bins.index_to_state(j).E;
                            let lnw = &mut self.bins.lnw;
                            if ej > *too_hi && ej <= energy.E {
                                if histogram[j] != 0 {
                                    lnw[j] = lnw[ihi];
                                    *num_states += 1;
                                } else {
                                    lnw[j] = Unitless::new(0.0);
                                }
                            }
                        }
                        *latest_parameter = *((energy.E - *too_lo) / min_T).value();
                        *tL = self.moves;
                        // The following rounds the energy to one of the bins.
                        let bin_e = self.bins.index_to_state(self.bins.state_to_index(energy)).E;
                        *too_hi = bin_e;
                    } else if energy.E < *too_lo {
                        let ilo = self.bins.state_to_index(State { E: *too_lo });
                        for j in 0..histogram.len() {
                            let ej = self.bins.index_to_state(j).E;
                            let lnw = &mut self.bins.lnw;
                            if ej < *too_lo && ej >= energy.E {
                                if histogram[j] != 0 {
                                    lnw[j] = lnw[ilo] + (ej - *too_lo) / min_T;
                                    if lnw[j] < Unitless::new(0.) {
                                        lnw[j] = Unitless::new(0.);
                                    }
                                    *num_states += 1;
                                } else {
                                    lnw[j] = Unitless::new(0.0);
                                }
                            }
                        }
                        *latest_parameter = *((*too_hi - energy.E) / min_T).value();
                        *tL = self.moves;
                        // The following rounds the energy to one of the bins.
                        let bin_e = self.bins.index_to_state(self.bins.state_to_index(energy)).E;
                        *too_lo = bin_e;
                    }
                }
                if *tL == self.moves {
                    // gamma_changed = true;
                    // Set tF to the latest discovery time in the
                    // range of energies that we actually care about.
                    let ilo = self.bins.state_to_index(State { E: *too_lo });
                    let ihi = self.bins.state_to_index(State { E: *too_hi });
                    let old_tF = *tF;
                    *tF = *self.bins.t_found[ilo..ihi + 1].iter().max().unwrap();
                    if old_tF == *tF {
                        // We didn't change gamma after all!
                        // gamma_changed = false;
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
                                num_states,
                                self.bins.min.pretty(),
                                too_lo.pretty(),
                                too_hi.pretty(),
                                (self.bins.min + self.bins.width * (histogram.len() - 1) as f64)
                                    .pretty()
                            );
                        }
                    }
                }
            }
            Method::Samc { .. } => {}
            Method::WL {
                ref mut gamma,
                ref mut lowest_hist,
                ref mut highest_hist,
                ref mut total_hist,
                ref mut hist,
                ref mut min_energy,
                num_states,
                inv_t,
                min_gamma,
            } => {
                if let Some(min_gamma) = min_gamma {
                    if *gamma < min_gamma {
                        // We are in a production run.
                        hist[i] += 1;
                        return;
                    }
                }
                if hist.len() != self.bins.lnw.len() {
                    // Oops, we found a new energy, so let's regroup.
                    if hist.len() == 0 || (*gamma != 1.0 && *lowest_hist > 0) {
                        println!(
                            "    WL: Starting fresh with {} energies",
                            self.bins.lnw.len()
                        );
                        // gamma_changed = true;
                        *gamma = 1.0;
                        *lowest_hist = 0;
                        *highest_hist = 0;
                        *total_hist = 0;
                        *hist = vec![0; self.bins.lnw.len()];
                        *min_energy = self.bins.min;
                    } else {
                        if *min_energy > self.bins.min {
                            println!("Found a new energy minimum: {:.5}", self.bins.min.pretty());
                        }
                        // We have to adjust our hist Vec, but can
                        // keep our counts! We just pretend we already
                        // knew there were more energies to be
                        // found...
                        while *min_energy > self.bins.min {
                            *min_energy -= self.bins.width;
                            hist.insert(0, 0);
                        }
                        while hist.len() < self.bins.lnw.len() {
                            hist.push(0);
                        }
                        *lowest_hist = hist.iter().cloned().min().unwrap();
                    }
                }
                hist[i] += 1;
                if hist[i] > *highest_hist {
                    *highest_hist = hist[i];
                }
                *total_hist += 1;
                let histogram = &self.bins.histogram;
                let max_energy = *min_energy + (hist.len() as f64) * self.bins.width;
                if hist[i] == *lowest_hist + 1
                    && hist.len() > 1
                    && (self.min_allowed_energy.is_none()
                        || self.min_allowed_energy.unwrap() >= *min_energy)
                    && (self.max_allowed_energy.is_none()
                        || self.max_allowed_energy.unwrap() <= max_energy)
                    && hist
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| histogram[*i] != 0)
                        .map(|(_, &h)| h)
                        .min()
                        == Some(*lowest_hist + 1)
                {
                    *lowest_hist = hist[i];
                    if (inv_t && *lowest_hist > 0)
                        || *lowest_hist as f64 >= 0.8 * *total_hist as f64 / num_states
                    {
                        // gamma_changed = true;
                        *gamma *= 0.5;
                        if *gamma > 1e-16 {
                            println!(
                                "    WL: ({}) We have reached flatness {:.2} with min {}!",
                                PrettyFloat(*gamma),
                                PrettyFloat(*lowest_hist as f64 * num_states / *total_hist as f64),
                                *lowest_hist
                            );
                            report_wl_flatness(
                                *lowest_hist,
                                *highest_hist,
                                *total_hist,
                                num_states,
                                hist,
                                &self.bins,
                                *min_energy,
                            );
                            println!();
                        }
                        for h in hist.iter_mut() {
                            *h = 0;
                        }
                        *total_hist = 0;
                        *lowest_hist = 0;
                        *highest_hist = 0;
                        if let Some(min_gamma) = min_gamma {
                            if *gamma < min_gamma {
                                // Switch to a "production" run.
                                *gamma = 0.0;
                                println!("      WL: BEGINNING PRODUCTION MODE!");
                            }
                        }
                    }
                    if inv_t && *gamma < (num_states as f64) / (self.moves as f64) {
                        println!("    1/t-WL:  Switching to 1/t!");
                        switch_to_samc = Some(num_states as f64);
                    }
                }
            }
        }
        if let Some(t0) = switch_to_samc {
            self.method = Method::Samc { t0 };
        }
        // if gamma_changed {
        //     self.movies.new_gamma(self.moves, gamma);
        //     self.movies.new_gamma(self.moves, self.gamma());
        // }
    }

    /// Estimate the temperature for a given energy
    pub fn temperature(&self, energy: State) -> Energy {
        let i = self.state_to_index(energy);
        let energy = energy.E;
        let lnwi = self.bins.lnw[i];
        let mut Tlo = Energy::new(0.);
        for j in 0..i {
            let lnwj = self.bins.lnw[j];
            if lnwj > Unitless::new(0.) {
                let Tj = (energy - self.index_to_state(j).E) / (lnwi - lnwj);
                if Tj > Tlo {
                    Tlo = Tj;
                }
            }
        }
        let mut Thi = Energy::new(1e300);
        for j in i + 1..self.bins.lnw.len() {
            let lnwj = self.bins.lnw[j];
            if lnwj > Unitless::new(0.) {
                let Tj = (energy - self.index_to_state(j).E) / (lnwi - lnwj);
                if Tj < Thi {
                    Thi = Tj;
                }
            }
        }
        if Thi > Tlo && Tlo > Energy::new(0.) {
            return 0.5 * (Thi + Tlo);
        }
        if Tlo > Energy::new(0.) {
            return Tlo;
        }
        return Thi;
    }
}

impl<S: System> EnergyMC<S> {
    fn gamma(&self) -> f64 {
        match self.method {
            Method::Canonical {..} => 0.0,
            Method::Sad {
                num_states,
                tF,
                version,
                latest_parameter,
                ..
            } => version.compute_gamma(
                latest_parameter,
                self.moves as f64,
                tF as f64,
                num_states as f64,
            ),
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
        let ewidth = params
            .energy_bin
            .unwrap_or(system.delta_energy().unwrap_or(Energy::new(1.0)));
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
        let emin = ((system.energy() / ewidth).value().round() - 0.5) * ewidth;
        EnergyMC {
            method: Method::new(
                params._method,
                system.energy(),
                ewidth,
                params.min_allowed_energy,
                params.max_allowed_energy,
            ),
            moves: 0,
            time_L: 0,
            accepted_moves: 0,
            acceptance_rate: 0.5, // arbitrary starting guess.
            min_allowed_energy: params.min_allowed_energy,
            max_allowed_energy: params.max_allowed_energy,

            bins: Bins {
                histogram: vec![1],
                t_found: vec![0],
                lnw: vec![Unitless::new(0.0)],
                energy_total: vec![system.energy()],
                energy_squared_total: vec![system.energy() * system.energy()],
                min: emin,
                width: ewidth,
                extra: std::collections::HashMap::new(),
            },

            have_visited_since_maxentropy: vec![false],
            round_trips: vec![1],
            max_S: Unitless::new(0.),
            max_S_index: 0,

            translation_scale: match params._moves {
                MoveParams::TranslationScale(x) => x,
                _ => 0.05 * units::SIGMA,
            },
            move_plan: params._moves,
            system: system,

            rng,
            save_as: save_as,
            report: plugin::Report::from(params._report),
            movies: plugin::Movie::from(params._movies),
            save: plugin::Save::from(params._save),
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
        if self.moves % (self.bins.histogram.len() as u64 * self.bins.histogram.len() as u64 * 1000)
            == 0
        {
            self.system.verify_energy();
        }
        let e1 = State::new(&self.system);
        let recent_scale = (1.0 / self.moves as f64).sqrt();
        self.acceptance_rate *= 1. - recent_scale;
        if let Some(e2) = self.system.plan_move(&mut self.rng, self.translation_scale) {
            let mut out_of_bounds = false;
            if let Some(maxe) = self.max_allowed_energy {
                out_of_bounds = e2 > maxe && e2 > e1.E;
            }
            if let Some(mine) = self.min_allowed_energy {
                out_of_bounds = out_of_bounds || (e2 < mine && e2 < e1.E)
            }
            if !out_of_bounds {
                let e2 = State { E: e2 };
                self.prepare_for_state(e2);

                if !self.reject_move(e1, e2) {
                    self.accepted_moves += 1;
                    self.acceptance_rate += recent_scale;
                    self.system.confirm();
                }
            }
        }
        let energy = State::new(&self.system);
        let i = self.state_to_index(energy);

        // track the time we found each energy.
        if self.bins.histogram[i] == 0 {
            self.bins.t_found[i] = self.moves;
        }
        self.bins.histogram[i] += 1;
        self.bins.energy_total[i] += energy.E;
        self.bins.energy_squared_total[i] += energy.E * energy.E;
        for (k, d) in self.system.data_to_collect(self.moves).into_iter() {
            self.bins.accumulate_extra(k, i, d);
        }

        self.update_weights(energy);

        if self.bins.lnw[i] > self.max_S {
            self.max_S = self.bins.lnw[i];
            self.max_S_index = i;
            for x in self.have_visited_since_maxentropy.iter_mut() {
                *x = true;
            }
        } else if i == self.max_S_index {
            if self.state_to_index(e1) != i {
                for x in self.have_visited_since_maxentropy.iter_mut() {
                    *x = false;
                }
            }
        } else if !self.have_visited_since_maxentropy[i] {
            self.have_visited_since_maxentropy[i] = true;
            self.round_trips[i] += 1;
        }

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

fn report_wl_flatness(
    lowest_hist: u64,
    highest_hist: u64,
    total_hist: u64,
    num_states: f64,
    hist: &[u64],
    bins: &Bins,
    min_energy: Energy,
) {
    if total_hist > 0 && hist.len() > 1 {
        let mut lowest = 111111111;
        let mut highest = 111111111;
        for (i, &h) in hist.iter().enumerate() {
            if h == lowest_hist && bins.histogram[i] != 0 {
                lowest = i;
            }
            if h == highest_hist && bins.histogram[i] != 0 {
                highest = i;
            }
        }
        let lowest_energy = if lowest == 111111111 {
            bins.index_to_state(0).E - bins.width
        } else {
            bins.index_to_state(lowest).E
        };
        let highest_energy = if highest == 111111111 {
            bins.index_to_state(0).E - bins.width
        } else {
            bins.index_to_state(highest).E
        };
        if lowest_hist == 0 {
            let zeros = hist.iter().cloned().filter(|&x| x == 0).count();
            print!(
                "        WL:  {:.3}/{:.3} bins have been unexplored, down to {:.10}",
                PrettyFloat(zeros as f64),
                PrettyFloat(num_states),
                min_energy.pretty(),
            );
        } else {
            print!("        WL:  flatness {:.1} with min {:.2} at {:.7} and max {:.2} at {:.7} (with total {:.2})!",
                     PrettyFloat(lowest_hist as f64*num_states as f64
                                 / total_hist as f64),
                     PrettyFloat(lowest_hist as f64), lowest_energy.pretty(),
                     PrettyFloat(highest_hist as f64), highest_energy.pretty(),
                     PrettyFloat(total_hist as f64));
        }
    }
}

impl Method {
    fn report<S: MovableSystem>(&self, mc: &EnergyMC<S>) {
        print!("    ");
        match self {
            Method::Canonical {..} => (), // Nothing special to report!
            Method::Sad { too_lo, too_hi, .. } => {
                let too_lo_count = mc.bins.histogram[mc.bins.state_to_index(State { E: *too_lo })];
                let too_hi_count = mc.bins.histogram[mc.bins.state_to_index(State { E: *too_hi })];
                print!(
                    "SAD: {:.5} ({:.3}) -> {:.5} ({:.3})",
                    too_lo.pretty(),
                    too_lo_count,
                    too_hi.pretty(),
                    too_hi_count,
                );
            }
            Method::Samc { .. } => {}
            Method::WL {
                lowest_hist,
                highest_hist,
                total_hist,
                num_states,
                ref hist,
                min_energy,
                ..
            } => {
                report_wl_flatness(
                    *lowest_hist,
                    *highest_hist,
                    *total_hist,
                    *num_states,
                    hist,
                    &mc.bins,
                    *min_energy,
                );
            }
        }
        println!(
            " [gamma = {:.2}]",
            crate::prettyfloat::PrettyFloat(mc.gamma())
        );
    }
}

#[derive(Serialize, Deserialize, Debug)]
struct Logger;
impl<S: MovableSystem + serde::Serialize + serde::de::DeserializeOwned> Plugin<EnergyMC<S>>
    for Logger
{
    fn log(&self, mc: &EnergyMC<S>, _sys: &S) {
        mc.method.report(&mc);
    }
}
