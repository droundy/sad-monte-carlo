//! An implementation of SAD (Statistical Association, Dynamical version).

#![allow(non_snake_case)]

use ::system::*;
use super::*;

use super::plugin::Plugin;
use dimensioned::Dimensionless;
use rand::Rng;
use std::default::Default;
use std::cell::{RefCell,Cell};

/// Parameters to configure a particular MC.
#[derive(Debug, ClapMe)]
pub enum MethodParams {
    /// Sad
    Sad {
        /// The minimum temperature we are interested in.
        min_T: Energy,
    },
    /// Samc
    Samc {
        /// The t0 parameter, determining how long to leave gamma=1.
        t0: u64,
    },
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
    _moves: MoveParams,
    _report: plugin::ReportParams,
    _movies: MoviesParams,
    _save: plugin::SaveParams,
}

impl Default for EnergyMCParams {
    fn default() -> Self {
        EnergyMCParams {
            _method: MethodParams::Sad { min_T: 0.2*units::EPSILON },
            seed: None,
            _moves: MoveParams::TranslationScale(0.05*units::SIGMA),
            _report: plugin::ReportParams::default(),
            _movies: MoviesParams::default(),
            _save: plugin::SaveParams::default(),
        }
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
    },
    /// Samc
    Samc {
        t0: u64,
    },
}

impl Method {
    fn new(p: MethodParams, E: Energy) -> Self {
        match p {
            MethodParams::Sad { min_T } =>
                Method::Sad {
                    min_T,
                    too_lo: E,
                    too_hi: E,
                    tL: 0,
                    num_states: 1,
                    highest_hist: 1,
                },
            MethodParams::Samc { t0 } => Method::Samc { t0 },
        }
    }
}

impl Bins {
    /// Find the index corresponding to a given energy.  This should
    /// panic if the energy is less than `min_energy_bin`.
    pub fn energy_to_index(&self, e: Energy) -> usize {
        *((e - self.min)/self.width).value() as usize
    }
    /// Find the energy corresponding to a given index.
    pub fn index_to_energy(&self, i: usize) -> Energy {
        self.min + (i as f64)*self.width
    }
    /// Make room in our arrays for a new energy value
    pub fn prepare_for_energy(&mut self, e: Energy) {
        assert!(self.width > Energy::new(0.0));
        while e < self.min {
            // this is a little wasteful, but seems the easiest way to
            // ensure we end up with enough room.
            self.histogram.insert(0, 0);
            self.lnw.insert(0, Unitless::new(0.0));
            self.min -= self.width;
        }
        while e >= self.min + self.width*(self.lnw.len() as f64) {
            self.lnw.push(Unitless::new(0.0));
            self.histogram.push(0);
        }
    }
}

impl<S: System> EnergyMC<S> {
    /// Find the index corresponding to a given energy.  This should
    /// panic if the energy is less than `min_energy_bin`.
    pub fn energy_to_index(&self, e: Energy) -> usize {
        self.bins.energy_to_index(e)
    }
    /// Find the energy corresponding to a given index.
    pub fn index_to_energy(&self, i: usize) -> Energy {
        self.bins.index_to_energy(i)
    }

    /// This decides whether to reject the move based on the actual
    /// method in use.
    fn reject_move(&mut self, e1: Energy, e2: Energy) -> bool {
        let i1 = self.energy_to_index(e1);
        let i2 = self.energy_to_index(e2);
        match self.method {
            Method::Sad { too_lo, too_hi, min_T,  .. } => {
                let lnw = &self.bins.lnw;
                let lnw1 = if e1 < too_lo {
                    lnw[self.energy_to_index(too_lo)] + (e1 - too_lo)/min_T
                } else if e1 > too_hi {
                    lnw[self.energy_to_index(too_hi)]
                } else {
                    lnw[i1]
                };
                let lnw2 = if e2 < too_lo {
                    lnw[self.energy_to_index(too_lo)] + (e2 - too_lo)/min_T
                } else if e2 > too_hi {
                    lnw[self.energy_to_index(too_hi)]
                } else {
                    lnw[i2]
                };
                let rejected = lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp();
                if !rejected && self.bins.histogram[i2] == 0 && e2 < too_hi && e2 > too_lo {
                    // Here we do changes that need only happen when
                    // we encounter an energy we have never seen before.
                    match &mut self.method {
                        Method::Sad { ref mut num_states, ref mut tL, .. } => {
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
        }
    }
    /// This updates the lnw based on the actual method in use.
    fn update_weights(&mut self, energy: Energy) {
        let i = self.energy_to_index(energy);
        match self.method {
            Method::Sad { min_T, ref mut too_lo, ref mut too_hi, ref mut num_states,
                          ref mut tL, ref mut highest_hist, .. } => {
                let histogram = &self.bins.histogram;
                if *too_lo < *too_hi && energy >= *too_lo && energy <= *too_hi {
                    let t = self.moves as f64;
                    let tL = *tL as f64;
                    let dE = *too_hi - *too_lo;
                    let n = *num_states as f64;

                    let gamma = dE/(3.0*min_T*t)*(n*n + t*(t/tL - 1.0) + n*t)/(n*n + t*(t/tL - 1.0) + t);
                    self.bins.lnw[i] += gamma;
                }

                if histogram[i] > *highest_hist {
                    *highest_hist = histogram[i];
                    if energy > *too_hi {
                        *tL = self.moves;
                        let hmax = histogram[i] as f64;
                        let ihi = self.bins.energy_to_index(*too_hi);
                        for j in 0 .. histogram.len() {
                            let ej = self.bins.index_to_energy(j);
                            let lnw = &mut self.bins.lnw;
                            if ej > *too_hi && ej <= energy {
                                if histogram[j] != 0 {
                                    lnw[j] = lnw[ihi] + (histogram[j] as f64/hmax).ln();
                                    *num_states += 1;
                                } else {
                                    lnw[j] = Unitless::new(0.0);
                                }
                            }
                        }
                        *too_hi = energy;
                    } else if energy < *too_lo {
                        *tL = self.moves;
                        let hmax = Unitless::new(histogram[i] as f64);
                        let ilo = self.bins.energy_to_index(*too_lo);
                        for j in 0 .. histogram.len() {
                            let ej = self.bins.index_to_energy(j);
                            let lnw = &mut self.bins.lnw;
                            if ej < *too_lo && ej >= energy {
                                if histogram[j] != 0 {
                                    lnw[j] = lnw[ilo] + log(histogram[j] as f64/hmax) + (ej-*too_lo)/min_T;
                                    *num_states += 1;
                                } else {
                                    lnw[j] = Unitless::new(0.0);
                                }
                            }
                        }
                        *too_lo = energy;
                    }
                }
                if *tL == self.moves {
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
            Method::Samc { t0 } => {
                let t = self.moves;
                self.bins.lnw[i] += if t > t0 { t0 as f64/t as f64 } else { 1.0 };
            }
        }
    }
}

impl<S: MovableSystem> MonteCarlo for EnergyMC<S> {
    type Params = EnergyMCParams;
    type System = S;
    fn from_params(params: EnergyMCParams, system: S, save_as: ::std::path::PathBuf) -> Self {
        EnergyMC {
            method: Method::new(params._method, system.energy()),
            moves: 0,
            time_L: 0,
            accepted_moves: 0,
            acceptance_rate: 0.5, // arbitrary starting guess.
            bins: Bins {
                histogram: vec![1],
                lnw: vec![Unitless::new(0.0)],
                min: system.energy(),
                width: system.delta_energy().unwrap_or(Energy::new(1.0)),
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

    fn move_once(&mut self) {
        self.moves += 1;
        let e1 = self.system.energy();
        let recent_scale = (1.0/self.moves as f64).sqrt();
        self.acceptance_rate *= 1. - recent_scale;
        if let Some(_) = self.system.move_once(&mut self.rng, self.translation_scale) {
            let e2 = self.system.energy();
            self.bins.prepare_for_energy(e2);

            if self.reject_move(e1,e2) {
                self.system.undo();
            } else {
                self.accepted_moves += 1;
                self.acceptance_rate += recent_scale;
            }
        }
        let energy = self.system.energy();
        let i = self.energy_to_index(energy);

        self.bins.histogram[i] += 1;
        self.update_weights(e1);

        let plugins = [&self.report as &Plugin<Self>,
                       &self.movies,
                       &self.save,
        ];
        self.manager.run(self, &self.system, &plugins);
    }
    fn system(&self) -> &Self::System {
        &self.system
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

fn log(x: Unitless) -> Unitless {
    Unitless::new(x.ln())
}



/// Do we want movies? Where?
#[derive(ClapMe, Debug)]
pub struct MoviesParams {
    /// How often (logarithmically) do we want a movie frame? If this
    /// is 2.0, it means we want a frame every time the number of
    /// iterations doubles.
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
    next_frame: Cell<u64>,
    period: Cell<plugin::TimeToRun>,
    entropy: RefCell<Vec<Vec<f64>>>,
    time: RefCell<Vec<u64>>,
    energy: RefCell<Vec<Energy>>,
}

impl From<MoviesParams> for Movies {
    fn from(params: MoviesParams) -> Self {
        Movies {
            movie_time: params.movie_time,
            next_frame: Cell::new(1),
            period: Cell::new(if params.movie_time.is_some() {
                plugin::TimeToRun::Period(1)
            } else {
                plugin::TimeToRun::Never
            }),
            entropy: RefCell::new(Vec::new()),
            time: RefCell::new(Vec::new()),
            energy: RefCell::new(Vec::new()),
        }
    }
}
impl<S: MovableSystem> Plugin<EnergyMC<S>> for Movies {
    fn run(&self, mc: &EnergyMC<S>, _sys: &S) -> plugin::Action {
        if let Some(movie_time) = self.movie_time {
            let next_frame = self.next_frame.get();
            if mc.num_moves() == next_frame {
                // First, let's create the arrays for the time and
                // energy indices.
                self.time.borrow_mut().push(next_frame);
                let new_energy: Vec<_> =
                    (0 .. mc.bins.lnw.len()).map(|i| mc.index_to_energy(i)).collect();
                let old_energy = self.energy.replace(new_energy.clone());

                let histogram = &mc.bins.histogram;
                let lnw = &mc.bins.lnw;
                let entropy: Vec<_> = (0..new_energy.len()).map(|i| {
                    if let Method::Sad { too_lo, too_hi, highest_hist, .. } = mc.method {
                        if histogram[i] == 0 {
                            return 0.0;
                        }
                        let e = mc.index_to_energy(i);
                        if e < too_lo {
                            return lnw[mc.energy_to_index(too_lo)].value()
                                + (histogram[i] as f64/highest_hist as f64).ln();
                        }
                        if e > too_hi {
                            return lnw[mc.energy_to_index(too_hi)].value()
                                + (histogram[i] as f64/highest_hist as f64).ln();
                        }
                    }
                    *lnw[i].value()
                }).collect();
                if new_energy == old_energy {
                    // We can just add a row.
                    let mut S = self.entropy.borrow_mut();
                    S.push(entropy);
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

                    S.push(entropy);
                }

                // Now decide when we need the next frame to be.
                let following_frame = (next_frame as f64*movie_time) as u64;
                self.next_frame.set(following_frame);
                self.period.set(plugin::TimeToRun::Period(following_frame-next_frame));
            }
        }
        plugin::Action::None
    }
    fn run_period(&self) -> plugin::TimeToRun { self.period.get() }
}
