//! An implementation of SAD (Statistical Association, Dynamical version).

#![allow(non_snake_case)]

use ::system::*;
use super::*;

use super::plugin::Plugin;
use dimensioned::Dimensionless;
use rand::Rng;
use std::default::Default;
use ndarray::{Array2, Axis};
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

/// The parameters needed to configure a simulation.
#[derive(Debug, ClapMe)]
pub struct EnergyMCParams {
    /// The actual method.
    pub _method: MethodParams,
    /// The seed for the random number generator.
    pub seed: Option<u64>,
    _report: plugin::ReportParams,
    _movies: MoviesParams,
    _save: plugin::SaveParams,
}

impl Default for EnergyMCParams {
    fn default() -> Self {
        EnergyMCParams {
            _method: MethodParams::Sad { min_T: 0.2*units::EPSILON },
            seed: None,
            _report: plugin::ReportParams::default(),
            _movies: MoviesParams::default(),
            _save: plugin::SaveParams::default(),
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
    /// The number of moves that have been rejected.
    pub rejected_moves: u64,
    /// The number of times we have been at each energy.
    pub histogram: Vec<u64>,
    /// The ln weight for each energy bin.
    pub lnw: Vec<Unitless>,
    /// The lowest allowed energy in any bin.
    pub min_energy_bin: Energy,
    /// The energy bin size.
    pub energy_bin: Energy,
    /// The max-entropy energy.
    pub max_entropy_energy: Energy,
    /// The max-entropy energy.
    pub max_S: Unitless,


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
        n_found: u64,
        higheset_hist: u64,
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
                    n_found: 1,
                    higheset_hist: 1,
                },
            MethodParams::Samc { t0 } => Method::Samc { t0 },
        }
    }
}

impl<S: System> EnergyMC<S> {
    /// Find the index corresponding to a given energy.  This should
    /// panic if the energy is less than `min_energy_bin`.
    pub fn energy_to_index(&self, e: Energy) -> usize {
        *((e - self.min_energy_bin)/self.energy_bin).value() as usize
    }
    /// Find the energy corresponding to a given index.
    pub fn index_to_energy(&self, i: usize) -> Energy {
        self.min_energy_bin + (i as f64)*self.energy_bin
    }
    /// Make room in our arrays for a new energy value
    pub fn prepare_for_energy(&mut self, e: Energy) {
        assert!(self.energy_bin > Energy::new(0.0));
        while e < self.min_energy_bin {
            // this is a little wasteful, but seems the easiest way to
            // ensure we end up with enough room.
            self.histogram.insert(0, 0);
            self.lnw.insert(0, Unitless::new(0.0));
            self.min_energy_bin -= self.energy_bin;
        }
        while e >= self.min_energy_bin + self.energy_bin*(self.lnw.len() as f64) {
            self.lnw.push(Unitless::new(0.0));
            self.histogram.push(0);
        }
    }

    /// This decides whether to reject the move based on the actual
    /// method in use.
    fn reject_move(&mut self, e1: Energy, e2: Energy) -> bool {
        let i1 = self.energy_to_index(e1);
        let i2 = self.energy_to_index(e2);
        match self.method {
            Method::Sad { too_lo, too_hi, min_T,  .. } => {
                let lnw1 = if e1 < too_lo {
                    self.lnw[self.energy_to_index(too_lo)] + (e1 - too_lo)/min_T
                } else if e1 > too_hi {
                    self.lnw[self.energy_to_index(too_hi)]
                } else {
                    self.lnw[i1]
                };
                let lnw2 = if e2 < too_lo {
                    self.lnw[self.energy_to_index(too_lo)] + (e2 - too_lo)/min_T
                } else if e2 > too_hi {
                    self.lnw[self.energy_to_index(too_hi)]
                } else {
                    self.lnw[i2]
                };
                let rejected = lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp();
                if !rejected && self.histogram[i2] == 0 {
                    // Here we do changes that need only happen when
                    // we encounter an energy we have never seen before.
                    match &mut self.method {
                        Method::Sad { ref mut n_found, ref mut tL, .. } => {
                            *n_found += 1;
                            println!("    sad: [{}]  {}:  {} < {} ... {} < {} < {}",
                                     self.moves, n_found,
                                     self.min_energy_bin.value_unsafe,
                                     too_lo.value_unsafe,
                                     self.max_entropy_energy.value_unsafe,
                                     too_hi.value_unsafe,
                                     (self.min_energy_bin
                                      + self.energy_bin*(self.lnw.len()-1) as f64).value_unsafe);
                            *tL = self.moves;
                        }
                        _ => unreachable!()
                    }
                }
                rejected
            }
            Method::Samc { .. } => {
                let lnw1 = self.lnw[i1].value();
                let lnw2 = self.lnw[i2].value();
                lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp()
            }
        }
    }
    /// This updates the lnw based on the actual method in use.
    fn update_weights(&mut self, energy: Energy) {
        let i = self.energy_to_index(energy);
        match self.method {
            Method::Sad { too_lo, too_hi, min_T, n_found, tL, higheset_hist, .. } => {
                if too_lo < too_hi && energy >= too_lo && energy <= too_hi {
                    let t = self.moves as f64;
                    let tL = tL as f64;
                    let dE = too_hi - too_lo;
                    let n = n_found as f64;

                    let gamma = dE/(3.0*min_T*t)*(n*n + t*(t/tL - 1.0) + n*t)/(n*n + t*(t/tL - 1.0) + t);
                    self.lnw[i] += gamma;
                }

                let mut want_print = false;
                if self.histogram[i] > higheset_hist {
                    match &mut self.method {
                        Method::Sad { ref mut higheset_hist, .. } => {
                            *higheset_hist = self.histogram[i];
                        }
                        _ => unreachable!()
                    }
                    if energy > too_hi {
                        want_print = true;
                        let hmax = self.histogram[i] as f64;
                        let ihi = self.energy_to_index(too_hi);
                        for j in 0 .. self.histogram.len() {
                            let ej = self.index_to_energy(j);
                            if ej > too_hi && ej <= energy {
                                if self.histogram[j] != 0 {
                                    self.lnw[j] = self.lnw[ihi] + (self.histogram[j] as f64/hmax).ln();
                                } else {
                                    self.lnw[j] = Unitless::new(0.0);
                                }
                            }
                        }
                        match &mut self.method {
                            Method::Sad { ref mut too_hi, .. } => {
                                *too_hi = energy;
                            }
                            _ => unreachable!()
                        }
                    } else if energy < too_lo {
                        want_print = true;
                        let hmax = Unitless::new(self.histogram[i] as f64);
                        let ilo = self.energy_to_index(too_lo);
                        for j in 0 .. self.histogram.len() {
                            let ej = self.index_to_energy(j);
                            if ej < too_lo && ej >= energy {
                                if self.histogram[j] != 0 {
                                    self.lnw[j] = self.lnw[ilo] + log(self.histogram[j] as f64/hmax) + (ej-too_lo)/min_T;
                                } else {
                                    self.lnw[j] = Unitless::new(0.0);
                                }
                            }
                        }
                        match &mut self.method {
                            Method::Sad { ref mut too_lo, .. } => {
                                *too_lo = energy;
                            }
                            _ => unreachable!()
                        }
                    }
                }
                if want_print {
                    println!("    sad: [{}]  {}:  {} < {} ... {} < {} < {}",
                             self.moves, n_found,
                             self.min_energy_bin.value_unsafe,
                             too_lo.value_unsafe,
                             self.max_entropy_energy.value_unsafe,
                             too_hi.value_unsafe,
                             (self.min_energy_bin
                              + self.energy_bin*(self.lnw.len()-1) as f64).value_unsafe);
                }
            }
            Method::Samc { t0 } => {
                let t = self.moves;
                self.lnw[i] += if t > t0 { t0 as f64/t as f64 } else { 1.0 };
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
            rejected_moves: 0,
            histogram: vec![1],
            lnw: vec![Unitless::new(0.0)],
            min_energy_bin: system.energy(),
            max_entropy_energy: system.energy(),
            max_S: Unitless::new(0.0),
            energy_bin: system.delta_energy().unwrap_or(Energy::new(1.0)),
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
        if let Some(_) = self.system.move_once(&mut self.rng, Length::new(0.1)) {
            let e2 = self.system.energy();
            self.prepare_for_energy(e2);

            if self.reject_move(e1,e2) {
                self.system.undo();
                self.rejected_moves += 1;
            }
        } else {
            // The system itself rejected the move.
            self.rejected_moves += 1;
        }
        let energy = self.system.energy();
        let i = self.energy_to_index(energy);

        self.histogram[i] += 1;
        self.update_weights(e1);

        if self.lnw[i] > self.max_S {
            self.max_S = self.lnw[i];
            self.max_entropy_energy = energy;
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
    fn num_moves(&self) -> u64 {
        self.moves
    }
    fn num_rejected_moves(&self) -> u64 {
        self.rejected_moves
    }
    fn save_as(&self) -> ::std::path::PathBuf {
        self.save_as.clone()
    }
}

fn log(x: Unitless) -> Unitless {
    Unitless::new(x.ln())
}
fn exp(x: Unitless) -> Unitless {
    Unitless::new(x.exp())
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
    period: Cell<Option<u64>>,
    entropy: RefCell<Array2<f64>>,
    time: RefCell<Vec<u64>>,
    energy: RefCell<Vec<Energy>>,
}

impl From<MoviesParams> for Movies {
    fn from(params: MoviesParams) -> Self {
        Movies {
            movie_time: params.movie_time,
            next_frame: Cell::new(1),
            period: Cell::new(params.movie_time.map(|_| 1)),
            entropy: RefCell::new(Array2::zeros((0,0))),
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
                    (0 .. mc.histogram.len()).map(|i| mc.index_to_energy(i)).collect();
                let old_energy = self.energy.replace(new_energy.clone());

                let entropy = Array2::from_shape_fn((1,new_energy.len()), |(_,i)| {
                    if let Method::Sad { too_lo, too_hi, higheset_hist, .. } = mc.method {
                        if mc.histogram[i] == 0 {
                            return 0.0;
                        }
                        let e = mc.index_to_energy(i);
                        if e < too_lo {
                            return mc.lnw[mc.energy_to_index(too_lo)].value()
                                + (mc.histogram[i] as f64/higheset_hist as f64).ln();
                        }
                        if e > too_hi {
                            return mc.lnw[mc.energy_to_index(too_hi)].value()
                                + (mc.histogram[i] as f64/higheset_hist as f64).ln();
                        }
                    }
                    *mc.lnw[i].value()
                });
                if new_energy == old_energy {
                    // We can just add a row.
                    let mut S = self.entropy.borrow_mut();
                    *S = stack!(Axis(0), *S, entropy);
                } else {
                    let old_S = self.entropy.replace(
                        Array2::zeros((self.time.borrow().len(),
                                       new_energy.len())));
                    // We need to change the shape of the array.
                }

                // Now decide when we need the next frame to be.
                let following_frame = (next_frame as f64*movie_time) as u64;
                self.next_frame.set(following_frame);
                self.period.set(Some(following_frame-next_frame));
            } else {
                println!("wrong time at {} vs {}", mc.num_moves(), next_frame);
            }
        }
        plugin::Action::None
    }
    fn run_period(&self) -> Option<u64> { self.period.get() }
}
