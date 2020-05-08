//! An implementation of SAD (Statistical Association, Dynamical version).

#![allow(non_snake_case)]

use crate::system::*;
use super::*;

use super::plugin::Plugin;
use dimensioned::Dimensionless;
use rand::{Rng, SeedableRng};
use std::default::Default;
use crate::prettyfloat::PrettyFloat;

/// Parameters to configure a particular MC.
#[derive(Debug, AutoArgs)]
pub enum MethodParams {
    /// Samc
    Samc {
        /// Use the SAMC algorithm, with the specified t0 parameter
        t0: f64,
    },
    /// Use the Wang-Landau algorithm
    WL,
}

/// The parameters needed to configure a simulation.
#[derive(Debug, AutoArgs)]
pub struct EnergyNumberMCParams {
    /// The actual method.
    pub _method: MethodParams,
    /// The seed for the random number generator.
    pub seed: Option<u64>,
    /// The energy binsize.
    energy_bin: Option<Energy>,
    _report: plugin::ReportParams,
    _movies: plugin::MovieParams,
    _save: plugin::SaveParams,
}

impl Default for EnergyNumberMCParams {
    fn default() -> Self {
        EnergyNumberMCParams {
            _method: MethodParams::WL,
            seed: None,
            energy_bin: None,
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
    /// The number of atoms
    pub N: usize,
}
impl ::std::fmt::Display for State {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "{}:{}", self.N, PrettyFloat(*(self.E/units::EPSILON).value()))
    }
}
impl State {
    /// Find the state of a system
    pub fn new<S: GrandSystem>(sys: &S) -> State {
        State { E: sys.energy(), N: sys.num_atoms() }
    }
}

/// This defines the energy bins.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Bins {
    /// The lowest allowed energy in any bin.
    pub min: Energy,
    /// The energy bin size.
    pub width: Energy,
    /// The number of times we have been at each state.
    pub histogram: Vec<u64>,
    /// The ln weight for each state.
    pub lnw: Vec<Unitless>,

    /// Whether we have seen this since the last visit to maxentropy.
    have_visited_since_maxentropy: Vec<bool>,
    /// How many round trips have we seen at this energy.
    round_trips: Vec<u64>,
    /// The maximum entropy we have seen.
    max_S: Unitless,
    /// The index with the maximum entropy.
    max_S_index: usize,
    /// The maximum number of atoms
    max_N: usize,
    /// The number of energy bins
    num_E: usize,
    /// The number of bins occupied
    num_states: usize,
    /// Time of last new discovery
    t_last: u64,
}

/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct EnergyNumberMC<S> {
    /// The system we are simulating.
    pub system: S,
    /// The method we use
    method: Method,
    /// The number of moves that have been made.
    pub moves: u64,
    /// The number of moves that have been accepted.
    pub accepted_moves: u64,
    /// The energy bins.
    pub bins: Bins,
    /// The random number generator.
    pub rng: crate::rng::MyRng,
    /// Where to save the resume file.
    pub save_as: ::std::path::PathBuf,
    report: plugin::Report,
    movies: plugin::Movie,
    save: plugin::Save,
    manager: plugin::PluginManager,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
enum Method {
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
        bins: Bins,
    }
}

impl Method {
    fn new<S: GrandSystem>(p: MethodParams, sys: &S, ewidth: Energy) -> Self {
        match p {
            MethodParams::Samc { t0 } => Method::Samc { t0 },
            MethodParams::WL => Method::WL {
                gamma: 1.0,
                lowest_hist: 1,
                highest_hist: 1,
                total_hist: 0,
                bins: Bins {
                    histogram: vec![1],
                    lnw: vec![Unitless::new(0.0)],
                    min: sys.energy() - ewidth*0.5, // center initial energy in a bin!
                    width: ewidth,
                    have_visited_since_maxentropy: vec![false],
                    round_trips: vec![1],
                    max_S: Unitless::new(0.),
                    max_S_index: 0,
                    max_N: 0,
                    num_E: 1,
                    num_states: 1,
                    t_last: 1,
                }
            },
        }
    }
    fn report<S: GrandSystem>(&self, mc: &EnergyNumberMC<S>) {
        print!("    ");
        // match self {
        //     Method::Samc { .. } => {
        //     }
        // }
        println!(" [gamma = {:.2}]", crate::prettyfloat::PrettyFloat(mc.gamma()));
    }
}

impl Bins {
    /// Find the index corresponding to a given energy.  This should
    /// panic if the energy is less than `min`.
    pub fn state_to_index(&self, s: State) -> usize {
        (*((s.E - self.min)/self.width).value() as usize)*(self.max_N+1) + s.N
    }
    /// Find the energy corresponding to a given index.
    pub fn index_to_state(&self, i: usize) -> State {
        State {
            E: self.min + ((i/(self.max_N+1)) as f64 + 0.5)*self.width,
            N: i % (self.max_N+1),
        }
    }
    fn emax(&self) -> Energy {
        self.min + self.width*(self.num_E as f64)
    }
    fn nbins(&self) -> usize {
        self.num_E*(self.max_N+1)
    }
    /// Make room in our arrays for a new energy value.  Returns true
    /// if we had to make changes.
    pub fn prepare_for_state(&mut self, e: State) -> bool {
        assert!(self.width > Energy::new(0.0));
        if e.N > self.max_N || e.E < self.min || e.E >= self.emax() {
            // this is a little wasteful, but seems the easiest way to
            // ensure we end up with enough room.
            let mut newbins = self.clone();
            newbins.max_N = if e.N > self.max_N { e.N } else { self.max_N };
            while newbins.min > e.E {
                newbins.min -= self.width;
                newbins.num_E += 1;
            }
            while newbins.emax() <= e.E {
                newbins.num_E += 1;
            }
            newbins.lnw = vec![Unitless::new(0.); newbins.nbins()];
            newbins.histogram = vec![0; newbins.nbins()];
            newbins.have_visited_since_maxentropy = vec![false; newbins.nbins()];
            newbins.round_trips = vec![0; newbins.nbins()];
            for i in 0..self.lnw.len() {
                let s = self.index_to_state(i);
                let j = newbins.state_to_index(s);
                newbins.lnw[j] = self.lnw[i];
                newbins.histogram[j] = self.histogram[i];
                newbins.have_visited_since_maxentropy[j] = self.have_visited_since_maxentropy[i];
                newbins.round_trips[j] = self.round_trips[i];
            }
            *self = newbins;
            true
        } else {
            false
        }
    }
}

impl<S: GrandSystem> EnergyNumberMC<S> {
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
            Method::Samc { .. } => {
                let lnw1 = self.bins.lnw[i1].value();
                let lnw2 = self.bins.lnw[i2].value();
                lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp()
            }
            Method::WL { ref mut bins, ref mut gamma, ref mut lowest_hist, ref mut highest_hist, ref mut total_hist, .. } => {
                if bins.prepare_for_state(e2) {
                    // Oops, we found a new energy, so let's regroup.
                    if *gamma != 1.0 || bins.histogram.len() == 0 {
                        // self.movies.new_gamma(self.moves, *gamma);
                        // self.movies.new_gamma(self.moves, 1.0);
                        println!("    WL: Starting fresh with {} energies",
                                 self.bins.lnw.len());
                        *gamma = 1.0;
                        *lowest_hist = 0;
                        *highest_hist = 0;
                        *total_hist = 0;
                    } else {
                        // We can keep our counts!
                    }
                }
                let lnw1 = self.bins.lnw[i1].value();
                let lnw2 = self.bins.lnw[i2].value();
                let rejected = lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp();
                rejected
            }
        }
    }
    /// This updates the lnw based on the actual method in use.
    fn update_weights(&mut self, energy: State) {
        let i = self.state_to_index(energy);
        let gamma = self.gamma(); // compute gamma out front...
        self.bins.lnw[i] += gamma;
        let mut gamma_changed = false;
        match self.method {
            Method::Samc { .. } => {}
            Method::WL { ref mut gamma,
                         ref mut lowest_hist, ref mut highest_hist, ref mut total_hist, ref mut bins } => {
                let num_states = self.bins.num_states as f64;
                bins.histogram[i] += 1;
                if bins.histogram[i] > *highest_hist {
                    *highest_hist = bins.histogram[i];
                }
                *total_hist += 1;
                let histogram = &self.bins.histogram;
                if bins.histogram[i] == *lowest_hist + 1
                    && bins.histogram.iter().enumerate()
                          .filter(|(i,_)| histogram[*i] != 0)
                          .map(|(_,&h)|h).min() == Some(*lowest_hist+1)
                {
                    *lowest_hist = bins.histogram[i];
                    if bins.histogram.len() > 1 && *lowest_hist as f64 >= 0.8**total_hist as f64 / num_states {
                        gamma_changed = true;
                        *gamma *= 0.5;
                        if *gamma > 1e-16 {
                            println!("    WL:  We have reached flatness {:.2} with min {}!",
                                     PrettyFloat(*lowest_hist as f64*num_states
                                                 / *total_hist as f64),
                                     *lowest_hist);
                            println!("         gamma => {}", PrettyFloat(*gamma));
                        }
                        bins.histogram.iter_mut().map(|x| *x = 0).count();
                        *total_hist = 0;
                        *lowest_hist = 0;
                    }
                }
            }
        }
        if gamma_changed {
            // self.movies.new_gamma(self.moves, gamma);
            // self.movies.new_gamma(self.moves, self.gamma());
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

impl<S> EnergyNumberMC<S> {
    fn gamma(&self) -> f64 {
        match self.method {
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

impl<S: GrandSystem> MonteCarlo for EnergyNumberMC<S> {
    type Params = EnergyNumberMCParams;
    type System = S;
    fn from_params(params: EnergyNumberMCParams, system: S, save_as: ::std::path::PathBuf) -> Self {
        let ewidth = params.energy_bin.unwrap_or(system.delta_energy().unwrap_or(Energy::new(1.0)));
        EnergyNumberMC {
            method: Method::new(params._method, &system, ewidth),
            moves: 0,
            accepted_moves: 0,
            bins: Bins {
                histogram: vec![1],
                lnw: vec![Unitless::new(0.0)],
                min: system.energy() - ewidth*0.5, // center initial energy in a bin!
                width: ewidth,
                have_visited_since_maxentropy: vec![false],
                round_trips: vec![1],
                max_S: Unitless::new(0.),
                max_S_index: 0,
                max_N: 0,
                num_E: 1,
                num_states: 1,
                t_last: 1,
            },
            system: system,

            rng: crate::rng::MyRng::seed_from_u64(params.seed.unwrap_or(0)),
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
        let e1 = State::new(&self.system);

        let option_e2 = if self.rng.gen::<u64>() & 1 == 0 {
            // add
            self.system.plan_add(&mut self.rng).map(|e2| State { E: e2, N: e1.N+1 })
        } else {
            // remove
            if e1.N > 0 {
                Some(State { E: self.system.plan_remove(&mut self.rng), N: e1.N-1})
            } else {
                None
            }
        };
        if let Some(e2) = option_e2 {
            self.bins.prepare_for_state(e2);
            if !self.reject_move(e1,e2) {
                self.accepted_moves += 1;
                self.system.confirm();
            }
        }
        let energy = State::new(&self.system);
        let i = self.state_to_index(energy);

        if self.bins.histogram[i] == 0 {
            self.bins.num_states += 1;
            self.bins.t_last = self.moves;
        }
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

        let plugins = [&self.report as &dyn Plugin<Self>,
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
impl<S: GrandSystem> Plugin<EnergyNumberMC<S>> for Logger {
    fn log(&self, mc: &EnergyNumberMC<S>, _sys: &S) {
        mc.method.report(&mc);
    }
}
