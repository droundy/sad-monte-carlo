//! An implementation of SAD (Statistical Association, Dynamical version).

#![allow(non_snake_case)]

use ::system::*;
use super::*;

use super::plugin::Plugin;
use dimensioned::Dimensionless;
use rand::Rng;
use std::default::Default;
use std::cell::{RefCell,Cell};
use ::prettyfloat::PrettyFloat;

fn min(a: f64, b: f64) -> f64 { if a < b { a } else { b } }
fn max(a: f64, b: f64) -> f64 { if a > b { a } else { b } }

/// Which experimental version of SAD are we doing?
#[derive(Serialize, Deserialize, Debug, ClapMe, Clone, Copy)]
pub enum SadVersion {
    /// Sad
    Sad,
    /// New sad with 1/t behavior
    SadOverT,
    /// New sad with 1/t^2 behavior
    SadOverT2
}

impl SadVersion {
    fn compute_gamma(&self, latest_parameter: f64, t: f64, tL: f64, num_states: f64) -> f64 {
        if latest_parameter == 0.0 {
            0.0
        } else {
            match *self {
                SadVersion::Sad =>
                    (latest_parameter + t/tL)/(latest_parameter + t/num_states*(t/tL)),
                SadVersion::SadOverT => latest_parameter/(latest_parameter + t),
                SadVersion::SadOverT2 =>
                    (latest_parameter+t*num_states)/(latest_parameter + t*t),
            }
        }
    }
    fn compute_latest_parameter(&self, dE_over_T: f64, Ns: f64, tL: f64,
                                old_dE_over_T: f64, old_Ns: f64, old_tL: f64,
                                previous_parameter: f64) -> f64 {
        match *self {
            SadVersion::Sad => dE_over_T,
            SadVersion::SadOverT => {
                let old_Stot = old_dE_over_T*(1./3.)*old_Ns;
                let new_Stot = dE_over_T*(1./3.)*Ns;
                let f = old_Stot*min(1., (tL/old_tL).ln());
                max(previous_parameter, new_Stot - f)
            }
            SadVersion::SadOverT2 => {
                (dE_over_T*Ns*0.5 - old_dE_over_T*old_Ns*0.5)*tL + previous_parameter
            }
        }
    }
}

/// Parameters to configure a particular MC.
#[derive(Debug, ClapMe)]
pub enum MethodParams {
    /// Sad
    Sad {
        /// The minimum temperature we are interested in.
        min_T: Energy,
        /// Which experimental version of SAD?
        version: SadVersion,
    },
    /// Samc
    Samc {
        /// The t0 parameter, determining how long to leave gamma=1.
        t0: f64,
    },
    /// Wang-Landau
    WL,
}

/// Parameters to configure the moves.
#[derive(Serialize, Deserialize, Debug, ClapMe)]
pub enum MoveParams {
    /// The rms distance of moves
    TranslationScale(Length),
    /// Adjust translation scale to reach acceptance rate.
    AcceptanceRate(f64),
}

/// The parameters needed to configure a simulation.
#[derive(Debug, ClapMe)]
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
    _movies: MoviesParams,
    _save: plugin::SaveParams,
}

impl Default for EnergyMCParams {
    fn default() -> Self {
        EnergyMCParams {
            _method: MethodParams::Sad {
                min_T: 0.2*units::EPSILON,
                version: SadVersion::Sad,
            },
            seed: None,
            min_allowed_energy: None,
            max_allowed_energy: None,
            energy_bin: None,
            _moves: MoveParams::TranslationScale(0.05*units::SIGMA),
            _report: plugin::ReportParams::default(),
            _movies: MoviesParams::default(),
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
        write!(f, "{}", PrettyFloat(*(self.E/units::EPSILON).value()))
    }
}
impl State {
    /// Find the state of a system
    pub fn new<S: System>(sys: &S) -> State {
        State { E: sys.energy() }
    }
}

/// This defines the energy bins.
#[derive(Serialize, Deserialize, Debug)]
pub struct Bins {
    /// The lowest allowed energy in any bin.
    pub min: Energy,
    /// The energy bin size.
    pub width: Energy,
    /// The number of times we have been at each energy.
    pub histogram: Vec<u64>,
    /// The ln weight for each energy bin.
    pub lnw: Vec<Unitless>,

    /// Whether we have seen this since the last visit to maxentropy.
    have_visited_since_maxentropy: Vec<bool>,
    /// How many round trips have we seen at this energy.
    round_trips: Vec<u64>,
    /// The maximum entropy we have seen.
    max_S: Unitless,
    /// The index with the maximum entropy.
    max_S_index: usize,
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
    /// The energy bins.
    pub bins: Bins,
    /// The move plan
    pub move_plan: MoveParams,
    /// The current translation scale
    pub translation_scale: Length,
    /// The "recent" acceptance rate.
    pub acceptance_rate: f64,
    /// The random number generator.
    pub rng: ::rng::MyRng,
    /// Where to save the resume file.
    pub save_as: ::std::path::PathBuf,
    report: plugin::Report,
    movies: Movies,
    save: plugin::Save,
    manager: plugin::PluginManager,
}

#[derive(Serialize, Deserialize, Debug)]
enum Method {
    /// Sad
    Sad {
        min_T: Energy,
        too_lo: Energy,
        too_hi: Energy,
        tL: u64,
        num_states: u64,
        highest_hist: u64,
        version: SadVersion,
        latest_parameter: f64,
    },
    /// Samc
    Samc {
        t0: f64,
    },
    /// Wang Landau
    WL {
        gamma: f64,
        lowest_hist: u64,
        highest_hist: u64,
        total_hist: u64,
        num_states: f64,
        hist: Vec<u64>,
        min_energy: Energy,
    }
}

impl Method {
    fn new(p: MethodParams, E: Energy, dE: Energy,
           min_allowed_energy: Option<Energy>,
           max_allowed_energy: Option<Energy>) -> Self {
        match p {
            MethodParams::Sad { min_T, version } =>
                Method::Sad {
                    min_T,
                    too_lo: E,
                    too_hi: E,
                    tL: 0,
                    num_states: 1,
                    highest_hist: 1,
                    version: version,
                    latest_parameter: 0.,
                },
            MethodParams::Samc { t0 } => Method::Samc { t0 },
            MethodParams::WL => Method::WL {
                gamma: 1.0,
                lowest_hist: if min_allowed_energy.is_some() && max_allowed_energy.is_some() { 0 } else { 1 },
                highest_hist: 1,
                total_hist: 0,
                num_states: if let (Some(mine), Some(maxe)) = (min_allowed_energy, max_allowed_energy) {
                    *((maxe - mine)/dE).value()
                } else {
                    1.0
                },
                hist: Vec::new(),
                min_energy: E,
            },
        }
    }
}

impl Bins {
    /// Find the index corresponding to a given energy.  This should
    /// panic if the energy is less than `min`.
    pub fn state_to_index(&self, s: State) -> usize {
        *((s.E - self.min)/self.width).value() as usize
    }
    /// Find the energy corresponding to a given index.
    pub fn index_to_state(&self, i: usize) -> State {
        State { E: self.min + (i as f64 + 0.5)*self.width }
    }
    /// Make room in our arrays for a new energy value
    pub fn prepare_for_state(&mut self, e: State) {
        let e = e.E;
        assert!(self.width > Energy::new(0.0));
        while e < self.min {
            // this is a little wasteful, but seems the easiest way to
            // ensure we end up with enough room.
            self.histogram.insert(0, 0);
            self.lnw.insert(0, Unitless::new(0.0));
            self.have_visited_since_maxentropy.insert(0, true);
            self.round_trips.insert(0, 1);
            self.min -= self.width;
        }
        while e >= self.min + self.width*(self.lnw.len() as f64) {
            self.lnw.push(Unitless::new(0.0));
            self.histogram.push(0);
            self.have_visited_since_maxentropy.push(true);
            self.round_trips.push(1);
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

    /// This decides whether to reject the move based on the actual
    /// method in use.
    fn reject_move(&mut self, e1: State, e2: State) -> bool {
        let i1 = self.state_to_index(e1);
        let i2 = self.state_to_index(e2);
        match self.method {
            Method::Sad { too_lo, too_hi, min_T, version,  .. } => {
                let lnw = &self.bins.lnw;
                let lnw1 = if e1.E < too_lo {
                    lnw[self.state_to_index(State { E: too_lo })] + (e1.E - too_lo)/min_T
                } else if e1.E > too_hi {
                    lnw[self.state_to_index(State { E: too_hi })]
                } else {
                    lnw[i1]
                };
                let lnw2 = if e2.E < too_lo {
                    lnw[self.state_to_index(State { E: too_lo })] + (e2.E - too_lo)/min_T
                } else if e2.E > too_hi {
                    lnw[self.state_to_index(State { E: too_hi })]
                } else {
                    lnw[i2]
                };
                let rejected = lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp();
                if !rejected && self.bins.histogram[i2] == 0 && e2.E < too_hi && e2.E > too_lo {
                    // Here we do changes that need only happen when
                    // we encounter an energy we have never seen before.
                    match &mut self.method {
                        Method::Sad { ref mut num_states, ref mut tL,
                                      ref mut latest_parameter, .. } => {
                            *latest_parameter = version.compute_latest_parameter(
                                *((too_hi - too_lo)/min_T).value(),
                                (*num_states + 1) as f64,
                                self.moves as f64,
                                *((too_hi - too_lo)/min_T).value(),
                                *num_states as f64,
                                *tL as f64,
                                *latest_parameter);
                            *num_states += 1;
                            *tL = self.moves;
                        }
                        _ => unreachable!()
                    }
                }
                rejected
            }
            Method::Samc { .. } => {
                let lnw1 = self.bins.lnw[i1].value();
                let lnw2 = self.bins.lnw[i2].value();
                lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp()
            }
            Method::WL { ref mut num_states, lowest_hist, .. } => {
                let lnw1 = self.bins.lnw[i1].value();
                let lnw2 = self.bins.lnw[i2].value();
                let rejected = lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp();
                if !rejected && self.bins.histogram[i2] == 0 && lowest_hist > 0 {
                    *num_states += 1.0;
                }
                rejected
            }
        }
    }
    /// This updates the lnw based on the actual method in use.
    fn update_weights(&mut self, energy: State) {
        let i = self.state_to_index(energy);
        let gamma = self.gamma(); // compute gamma out front...
        let old_lnw = self.bins.lnw[i];
        self.bins.lnw[i] += gamma;
        let mut gamma_changed = false;
        match self.method {
            Method::Sad { min_T, ref mut too_lo, ref mut too_hi, ref mut num_states,
                          ref mut tL, ref mut highest_hist,
                          version, ref mut latest_parameter, .. } => {
                let histogram = &self.bins.histogram;
                if *too_lo > *too_hi || energy.E < *too_lo || energy.E > *too_hi {
                    // Ooops, we didn't want to add gamma after all...
                    self.bins.lnw[i] = old_lnw;
                }

                if histogram[i] > *highest_hist {
                    *highest_hist = histogram[i];
                    let old_num_states = *num_states;
                    if energy.E > *too_hi {
                        let ihi = self.bins.state_to_index(State { E: *too_hi });
                        for j in 0 .. histogram.len() {
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
                        *latest_parameter = version.compute_latest_parameter(
                            *((energy.E - *too_lo)/min_T).value(),
                            *num_states as f64,
                            self.moves as f64,
                            *((*too_hi - *too_lo)/min_T).value(),
                            old_num_states as f64,
                            *tL as f64,
                            *latest_parameter);
                        *tL = self.moves;
                        *too_hi = energy.E;
                    } else if energy.E < *too_lo {
                        let ilo = self.bins.state_to_index(State { E: *too_lo });
                        for j in 0 .. histogram.len() {
                            let ej = self.bins.index_to_state(j).E;
                            let lnw = &mut self.bins.lnw;
                            if ej < *too_lo && ej >= energy.E {
                                if histogram[j] != 0 {
                                    lnw[j] = lnw[ilo] + (ej-*too_lo)/min_T;
                                    if lnw[j] < Unitless::new(0.) {
                                        lnw[j] = Unitless::new(0.);
                                    }
                                    *num_states += 1;
                                } else {
                                    lnw[j] = Unitless::new(0.0);
                                }
                            }
                        }
                        *latest_parameter = version.compute_latest_parameter(
                            *((*too_hi - energy.E)/min_T).value(),
                            *num_states as f64,
                            self.moves as f64,
                            *((*too_hi - *too_lo)/min_T).value(),
                            old_num_states as f64,
                            *tL as f64,
                            *latest_parameter);
                        *tL = self.moves;
                        *too_lo = energy.E;
                    }
                }
                if *tL == self.moves {
                    gamma_changed = true;
                    // We just discovered a new important energy.
                    // Let's take this as an opportunity to revise our
                    // translation scale, and also to log the news.
                    if let MoveParams::AcceptanceRate(r) = self.move_plan {
                        let s = self.acceptance_rate/r;
                        let s = if s < 0.8 { 0.8 } else if s > 1.2 { 1.2 } else { s };
                        self.translation_scale *= s;
                        if !self.report.quiet {
                            println!("        new translation scale: {:.3}",
                                     self.translation_scale);
                            println!("        acceptance rate {:.1}% [long-term: {:.1}%]",
                                     100.0*self.acceptance_rate,
                                     100.0*self.accepted_moves as f64
                                     /self.moves as f64);
                        }
                    }
                    if !self.report.quiet {
                        println!("    sad: [{}]  {}:  {} < {} ... {} < {}",
                                 self.moves, num_states,
                                 self.bins.min.value_unsafe,
                                 too_lo.value_unsafe,
                                 too_hi.value_unsafe,
                                 (self.bins.min
                                  + self.bins.width*(histogram.len()-1) as f64).value_unsafe);
                    }
                }
            }
            Method::Samc { .. } => {}
            Method::WL { ref mut gamma,
                         ref mut lowest_hist, ref mut highest_hist, ref mut total_hist,
                         ref mut hist, ref mut min_energy, num_states } => {
                if hist.len() != self.bins.lnw.len() {
                    // Oops, we found a new energy, so let's regroup.
                    if hist.len() == 0 || (*gamma != 1.0 && *lowest_hist > 0) {
                        println!("    WL: Starting fresh with {} energies",
                                 self.bins.lnw.len());
                        gamma_changed = true;
                        *gamma = 1.0;
                        *lowest_hist = 0;
                        *highest_hist = 0;
                        *total_hist = 0;
                        *hist = vec![0; self.bins.lnw.len()];
                        *min_energy = self.bins.min;
                    } else {
                        // We have to adjust our hist Vec, but can
                        // keep our counts! We just pretend we already
                        // knew there were more energies to be
                        // found...
                        while *min_energy > self.bins.min {
                            *min_energy -= self.bins.width;
                            hist.insert(0,0);
                        }
                        while hist.len() < self.bins.lnw.len() {
                            hist.push(0);
                        }
                        for (i,h) in hist.iter().cloned().enumerate() {
                            if self.bins.histogram[i] > 0 && h < *lowest_hist {
                                *lowest_hist = h;
                            }
                        }
                    }
                }
                hist[i] += 1;
                if hist[i] > *highest_hist {
                    *highest_hist = hist[i];
                }
                *total_hist += 1;
                let histogram = &self.bins.histogram;
                let max_energy = *min_energy + (hist.len() as f64)*self.bins.width;
                if hist[i] == *lowest_hist + 1
                    && hist.len() > 1
                    && (self.min_allowed_energy.is_none() || self.min_allowed_energy.unwrap() >= *min_energy)
                    && (self.max_allowed_energy.is_none() || self.max_allowed_energy.unwrap() <= max_energy)
                    && hist.iter().enumerate()
                          .filter(|(i,_)| histogram[*i] != 0)
                          .map(|(_,&h)|h).min() == Some(*lowest_hist+1)
                {
                    *lowest_hist = hist[i];
                    if *lowest_hist as f64 >= 0.8**total_hist as f64 / num_states {
                        gamma_changed = true;
                        *gamma *= 0.5;
                        if *gamma > 1e-16 {
                            println!("    WL:  We have reached flatness {:.2} with min {}!",
                                     PrettyFloat(*lowest_hist as f64*num_states
                                                 / *total_hist as f64),
                                     *lowest_hist);
                            report_wl_flatness(*lowest_hist, *highest_hist, *total_hist,
                                               num_states, hist, &self.bins);
                        }
                        for h in hist.iter_mut() {
                            *h = 0;
                        }
                        *total_hist = 0;
                        *lowest_hist = 0;
                        *highest_hist = 0;
                    }
                }
            }
        }
        if gamma_changed {
            self.movies.new_gamma(self.moves, gamma);
            self.movies.new_gamma(self.moves, self.gamma());
        }
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
                let Tj = (energy - self.index_to_state(j).E)/(lnwi - lnwj);
                if Tj > Tlo { Tlo = Tj; }
            }
        }
        let mut Thi = Energy::new(1e300);
        for j in i+1 .. self.bins.lnw.len() {
            let lnwj = self.bins.lnw[j];
            if lnwj > Unitless::new(0.) {
                let Tj = (energy - self.index_to_state(j).E)/(lnwi - lnwj);
                if Tj < Thi { Thi = Tj; }
            }
        }
        if Thi > Tlo && Tlo > Energy::new(0.) {
            return 0.5*(Thi + Tlo);
        }
        if Tlo > Energy::new(0.) {
            return Tlo;
        }
        return Thi;
    }
}

impl<S> EnergyMC<S> {
    fn gamma(&self) -> f64 {
        match self.method {
            Method::Sad { num_states, tL, version, latest_parameter, .. } => {
                version.compute_gamma(latest_parameter, self.moves as f64,
                                      tL as f64, num_states as f64)
            }
            Method::Samc { t0 } => {
                let t = self.moves as f64;
                if t > t0 { t0/t } else { 1.0 }
            }
            Method::WL { gamma, .. } => {
                gamma
            }
        }
    }
}

impl<S: MovableSystem> MonteCarlo for EnergyMC<S> {
    type Params = EnergyMCParams;
    type System = S;
    fn from_params(params: EnergyMCParams, system: S, save_as: ::std::path::PathBuf) -> Self {
        let ewidth = params.energy_bin.unwrap_or(system.delta_energy().unwrap_or(Energy::new(1.0)));
        EnergyMC {
            method: Method::new(params._method, system.energy(), ewidth,
                                params.min_allowed_energy, params.max_allowed_energy),
            moves: 0,
            time_L: 0,
            accepted_moves: 0,
            acceptance_rate: 0.5, // arbitrary starting guess.
            min_allowed_energy: params.min_allowed_energy,
            max_allowed_energy: params.max_allowed_energy,
            bins: Bins {
                histogram: vec![1],
                lnw: vec![Unitless::new(0.0)],
                min: system.energy() - ewidth*0.5, // center initial energy in a bin!
                width: ewidth,
                have_visited_since_maxentropy: vec![false],
                round_trips: vec![1],
                max_S: Unitless::new(0.),
                max_S_index: 0,
            },
            translation_scale: match params._moves {
                MoveParams::TranslationScale(x) => x,
                _ => 0.05*units::SIGMA,
            },
            move_plan: params._moves,
            system: system,

            rng: ::rng::MyRng::from_u64(params.seed.unwrap_or(0)),
            save_as: save_as,
            report: plugin::Report::from(params._report),
            movies: Movies::from(params._movies),
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
        let e1 = State::new(&self.system);
        let recent_scale = (1.0/self.moves as f64).sqrt();
        self.acceptance_rate *= 1. - recent_scale;
        if let Some(e2) = self.system.plan_move(&mut self.rng, self.translation_scale) {
            let mut out_of_bounds = false;
            if let Some(maxe) = self.max_allowed_energy {
                out_of_bounds = e2 > maxe;
            }
            if let Some(mine) = self.min_allowed_energy {
                if e2 < mine {
                    out_of_bounds = true;
                }
            }
            if !out_of_bounds {
                let e2 = State { E: e2 };
                self.bins.prepare_for_state(e2);

                if !self.reject_move(e1,e2) {
                    self.accepted_moves += 1;
                    self.acceptance_rate += recent_scale;
                    self.system.confirm();
                }
            }
        }
        let energy = State::new(&self.system);
        let i = self.state_to_index(energy);

        self.bins.histogram[i] += 1;
        self.update_weights(energy);

        if self.bins.lnw[i] > self.bins.max_S {
            self.bins.max_S = self.bins.lnw[i];
            self.bins.max_S_index = i;
            for x in self.bins.have_visited_since_maxentropy.iter_mut() {
                *x = true;
            }
        } else if i == self.bins.max_S_index {
            if self.state_to_index(e1) != i {
                for x in self.bins.have_visited_since_maxentropy.iter_mut() {
                    *x = false;
                }
            }
        } else if !self.bins.have_visited_since_maxentropy[i] {
            self.bins.have_visited_since_maxentropy[i] = true;
            self.bins.round_trips[i] += 1;
        }

        let plugins = [&self.report as &Plugin<Self>,
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


/// Do we want movies? Where?
#[derive(ClapMe, Debug)]
pub struct MoviesParams {
    // How often (logarithmically) do we want a movie frame? If this
    // is 2.0, it means we want a frame every time the number of
    // iterations doubles.
    /// 2.0 means a frame every time iterations double.
    pub movie_time: Option<f64>,
}

impl Default for MoviesParams {
    fn default() -> Self {
        MoviesParams {
            movie_time: None,
        }
    }
}

/// A plugin that saves movie data.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Movies {
    movie_time: Option<f64>,
    which_frame: Cell<i32>,
    period: Cell<plugin::TimeToRun>,
    entropy: RefCell<Vec<Vec<f64>>>,
    histogram: RefCell<Vec<Vec<u64>>>,
    time: RefCell<Vec<u64>>,
    energy: RefCell<Vec<Energy>>,
    #[serde(default)]
    gamma: RefCell<Vec<f64>>,
    #[serde(default)]
    gamma_time: RefCell<Vec<u64>>,
}

impl From<MoviesParams> for Movies {
    fn from(params: MoviesParams) -> Self {
        Movies {
            movie_time: params.movie_time,
            which_frame: Cell::new(0),
            period: Cell::new(if params.movie_time.is_some() {
                plugin::TimeToRun::TotalMoves(1)
            } else {
                plugin::TimeToRun::Never
            }),
            entropy: RefCell::new(Vec::new()),
            histogram: RefCell::new(Vec::new()),
            time: RefCell::new(Vec::new()),
            energy: RefCell::new(Vec::new()),
            gamma: RefCell::new(Vec::new()),
            gamma_time: RefCell::new(Vec::new()),
        }
    }
}
impl Movies {
    fn new_gamma(&self, t: u64, gamma: f64) {
        self.gamma.borrow_mut().push(gamma);
        self.gamma_time.borrow_mut().push(t);
    }
}
impl<S: MovableSystem> Plugin<EnergyMC<S>> for Movies {
    fn run(&self, mc: &EnergyMC<S>, _sys: &S) -> plugin::Action {
        if let Some(movie_time) = self.movie_time {
            let moves = mc.num_moves();
            if plugin::TimeToRun::TotalMoves(moves) == self.period.get() {
                println!("Saving movie...");
                self.new_gamma(moves, mc.gamma());
                // First, let's create the arrays for the time and
                // energy indices.
                self.time.borrow_mut().push(moves);
                let new_energy: Vec<_> =
                    (0 .. mc.bins.lnw.len()).map(|i| mc.index_to_state(i).E).collect();
                let old_energy = self.energy.replace(new_energy.clone());

                let histogram = &mc.bins.histogram;
                let lnw = &mc.bins.lnw;
                let entropy: Vec<_> = (0..new_energy.len()).map(|i| {
                    if let Method::Sad { too_lo, too_hi, highest_hist, .. } = mc.method {
                        if histogram[i] == 0 {
                            return 0.0;
                        }
                        let e = mc.index_to_state(i).E;
                        if e < too_lo {
                            return lnw[mc.state_to_index(State { E: too_lo })].value()
                                + (histogram[i] as f64/highest_hist as f64).ln();
                        }
                        if e > too_hi {
                            return lnw[mc.state_to_index(State { E: too_hi })].value()
                                + (histogram[i] as f64/highest_hist as f64).ln();
                        }
                    }
                    *lnw[i].value()
                }).collect();
                if new_energy == old_energy {
                    // We can just add a row.
                    let mut S = self.entropy.borrow_mut();
                    S.push(entropy);
                    let mut hist_movie = self.histogram.borrow_mut();
                    hist_movie.push(histogram.clone());
                } else {
                    let energy = self.energy.borrow().clone();
                    let mut left_zeros = 0;
                    for e in energy.iter() {
                        if !old_energy.contains(e) {
                            left_zeros += 1;
                        } else {
                            break;
                        }
                    }
                    let mut right_zeros = 0;
                    for e in energy.iter().rev() {
                        if !old_energy.contains(e) {
                            right_zeros += 1;
                        } else {
                            break;
                        }
                    }
                    let mut S = self.entropy.borrow_mut();
                    for v in S.iter_mut() {
                        for _ in 0..right_zeros {
                            v.push(0.);
                        }
                        for _ in 0..left_zeros {
                            v.insert(0,0.);
                        }
                    }
                    let mut hist_movie = self.histogram.borrow_mut();
                    for v in hist_movie.iter_mut() {
                        for _ in 0..right_zeros {
                            v.push(0);
                        }
                        for _ in 0..left_zeros {
                            v.insert(0,0);
                        }
                    }
                    S.push(entropy);
                    hist_movie.push(histogram.clone());
                }

                // Now decide when we need the next frame to be.
                let mut which_frame = self.which_frame.get() + 1;
                let mut next_time = (movie_time.powi(which_frame) + 0.5) as u64;
                while next_time <= moves {
                    which_frame += 1;
                    next_time = (movie_time.powi(which_frame) + 0.5) as u64;
                }
                self.which_frame.set(which_frame);
                self.period.set(plugin::TimeToRun::TotalMoves(next_time));
                return plugin::Action::Save;
            }
        }
        plugin::Action::None
    }
    fn run_period(&self) -> plugin::TimeToRun { self.period.get() }

    /// This isn't really a movies thing, but there isn't a great
    /// reason to create yet another plugin for an EnergyMC-specific
    /// log message.
    fn log(&self, mc: &EnergyMC<S>, sys: &S) {
        let mut one_trip: Option<State> = None;
        let mut ten_trips: Option<State> = None;
        let mut hundred_trips: Option<State> = None;
        let mut thousand_trips: Option<State> = None;
        for (i, &trips) in mc.bins.round_trips.iter().enumerate() {
            if one_trip.is_none() && trips >= 1 {
                one_trip = Some(mc.index_to_state(i));
            }
            if ten_trips.is_none() && trips >= 10 {
                ten_trips = Some(mc.index_to_state(i));
            }
            if hundred_trips.is_none() && trips >= 100 {
                hundred_trips = Some(mc.index_to_state(i));
            }
            if thousand_trips.is_none() && trips >= 1000 {
                thousand_trips = Some(mc.index_to_state(i));
            }
        }
        let thousand_T = thousand_trips
            .map(|e| format!(" ({:.1})",
                             PrettyFloat(*(mc.temperature(e)/units::EPSILON).value())))
            .unwrap_or("".to_string());
        let hundred_T = hundred_trips
            .map(|e| format!(" ({:.1})",
                             PrettyFloat(*(mc.temperature(e)/units::EPSILON).value())))
            .unwrap_or("".to_string());
        let ten_T = ten_trips
            .map(|e| format!(" ({:.1})",
                             PrettyFloat(*(mc.temperature(e)/units::EPSILON).value())))
            .unwrap_or("".to_string());
        let thousand_trips = thousand_trips.map(|e| format!("{}", e))
            .unwrap_or("-".to_string());
        let ten_trips = ten_trips.map(|e| format!("{}", e))
            .unwrap_or("-".to_string());
        let one_trip = one_trip.map(|e| format!("{}", e))
            .unwrap_or("-".to_string());
        let hundred_trips = hundred_trips.map(|e| format!("{}", e))
            .unwrap_or("-".to_string());
        if !mc.report.quiet {
            println!("   {} * {}{} * {}{} * {}{} | {} currently {}",
                     one_trip,
                     ten_trips, ten_T,
                     hundred_trips, hundred_T,
                     thousand_trips, thousand_T,
                     mc.index_to_state(mc.bins.max_S_index).E/units::EPSILON,
                     sys.energy()/units::EPSILON,
            );
            if let Method::WL { lowest_hist, highest_hist, total_hist, num_states,
                                ref hist, .. } = mc.method {
                report_wl_flatness(lowest_hist, highest_hist, total_hist, num_states,
                                   hist, &mc.bins);
            }
        }
    }
}

fn report_wl_flatness(lowest_hist: u64, highest_hist: u64, total_hist: u64,
                      num_states: f64, hist: &[u64], bins: &Bins) {
    if total_hist > 0 && hist.len() > 1 {
        let mut lowest = 111111111;
        let mut highest = 111111111;
        for (i,&h) in hist.iter().enumerate() {
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
        println!("        WL:  flatness {:.1} with min {:.2} at {} and max {:.2} at {} (with total {})!",
                 PrettyFloat(lowest_hist as f64*num_states as f64
                             / total_hist as f64),
                 PrettyFloat(lowest_hist as f64), lowest_energy,
                 PrettyFloat(highest_hist as f64), highest_energy, total_hist);
    }
}
