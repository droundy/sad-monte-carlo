//! An implementation of broad histogram methods that have fixed
//! temperature and varying number of particles.

#![allow(non_snake_case)]

use super::*;
use crate::system::*;

use super::plugin::Plugin;
use crate::prettyfloat::PrettyFloat;
use dimensioned::Dimensionless;
use rand::{Rng, SeedableRng};
use std::cell::{Cell, RefCell};
use std::default::Default;

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
    /// Use the (experimental) SAD algorithm
    Sad {
        /// The maximum chemical potential that we care about.
        mu_max: Energy,
    },
}

/// Parameters to configure the moves.
#[derive(Serialize, Deserialize, Debug, AutoArgs, Clone, Copy)]
pub enum MoveParams {
    /// This means you chose to be explicit about translation scale etc.
    _Explicit {
        /// The rms distance of translations
        translation_scale: Length,
        /// The fraction of attempted changes that are adds or removes
        addremove_probability: f64,
    },
    /// Adjust translation scale and add/remove probability to reach acceptance rate.
    AcceptanceRate(f64),
}

/// The parameters needed to configure a simulation.
#[derive(Debug, AutoArgs)]
pub struct NumberMCParams {
    /// The actual method.
    pub _method: MethodParams,
    /// The seed for the random number generator.
    pub seed: Option<u64>,
    /// The temperature.
    T: Energy,
    /// The maximum allowed number of atoms
    max_N: Option<usize>,
    _moves: MoveParams,
    _report: plugin::ReportParams,
    movie: MoviesParams,
    _save: plugin::SaveParams,
}

impl Default for NumberMCParams {
    fn default() -> Self {
        NumberMCParams {
            _method: MethodParams::WL,
            T: Energy::new(1.0),
            seed: None,
            max_N: None,
            _moves: MoveParams::_Explicit {
                translation_scale: 0.05 * units::SIGMA,
                addremove_probability: 0.05,
            },
            _report: plugin::ReportParams::default(),
            movie: MoviesParams::default(),
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
        if self.E == Energy::new(0.) {
            write!(f, "{}", self.N)
        } else {
            write!(
                f,
                "{}:{}",
                self.N,
                PrettyFloat(*(self.E / units::EPSILON).value())
            )
        }
    }
}
impl State {
    /// Find the state of a system
    pub fn new<S: GrandSystem>(sys: &S) -> State {
        State {
            E: sys.energy(),
            N: sys.num_atoms(),
        }
    }
}

/// This defines the energy bins.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Bins {
    /// The number of times we have been at each state.
    pub histogram: Vec<u64>,
    /// The ln weight for each state.
    pub lnw: Vec<Unitless>,
    /// The total energy for each state.
    pub total_energy: Vec<Energy>,
    /// The total energy squared for each state.
    pub total_energy_squared: Vec<EnergySquared>,
    /// The current translation scale for each number of atoms
    pub translation_scale: Vec<Length>,
    /// The number of translation moves tried at each number of atoms
    pub num_translation_attempts: Vec<u64>,
    /// The number of translation moves accepted at each number of atoms
    pub num_translation_accepted: Vec<u64>,
    /// The iteration when we found each number of atoms.
    pub t_found: Vec<u64>,
    /// The number of adds/removes tried at each number of atoms
    pub num_addremove_attempts: u64,
    /// The number of adds/removes accepted at each number of atoms
    pub num_addremove_accepted: u64,

    /// Whether we have seen this since the last visit to maxlnw.
    have_visited_since_maxlnw: Vec<bool>,
    /// How many round trips have we seen at this energy.
    round_trips: Vec<u64>,
    /// The maximum lnw we have seen.
    maxlnw: Unitless,
    /// The index with the maximum lnw.
    maxlnw_index: usize,
    /// The number of bins occupied
    num_states: usize,
    /// Time of last new discovery
    t_last: u64,
}

/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct NumberMC<S> {
    /// The system we are simulating.
    pub system: S,
    /// The temperature.
    T: Energy,
    /// The method we use
    method: Method,
    /// The number of moves that have been made.
    pub moves: u64,
    /// The number of moves that have been accepted.
    pub accepted_moves: u64,
    /// The energy bins.
    pub bins: Bins,
    /// The maximum allowed number of atoms
    pub max_N: Option<usize>,
    /// The move plan
    pub move_plan: MoveParams,
    /// The current add/remove probability
    pub addremove_probability: f64,
    /// The random number generator.
    pub rng: crate::rng::MyRng,
    /// Where to save the resume file.
    pub save_as: ::std::path::PathBuf,
    report: plugin::Report,
    movies: Movies,
    save: plugin::Save,
    manager: plugin::PluginManager,
}

#[derive(Serialize, Deserialize, Debug, Clone)]
enum Method {
    /// Samc
    Samc { t0: f64 },
    /// Samc
    Sad {
        mu_max: Energy,
        too_many: usize,
        tF: u64,
        highest_hist: u64,
        latest_parameter: f64,
    },
    /// Wang Landau
    WL {
        gamma: f64,
        lowest_hist: u64,
        highest_hist: u64,
        total_hist: u64,
        bins: Bins,
    },
}

impl Method {
    fn new(p: MethodParams) -> Self {
        match p {
            MethodParams::Samc { t0 } => Method::Samc { t0 },
            MethodParams::Sad { mu_max } => Method::Sad {
                mu_max,
                too_many: 1,
                tF: 0,
                highest_hist: 1,
                latest_parameter: 0.,
            },
            MethodParams::WL => Method::WL {
                gamma: 1.0,
                lowest_hist: 1,
                highest_hist: 1,
                total_hist: 0,
                bins: Bins {
                    histogram: vec![1],
                    lnw: vec![Unitless::new(0.0)],
                    total_energy: vec![Energy::new(0.0)],
                    total_energy_squared: vec![EnergySquared::new(0.0)],
                    translation_scale: vec![0.05 * units::SIGMA],
                    num_translation_attempts: vec![0],
                    num_translation_accepted: vec![0],
                    t_found: vec![0],
                    num_addremove_attempts: 0,
                    num_addremove_accepted: 0,
                    have_visited_since_maxlnw: vec![false],
                    round_trips: vec![1],
                    maxlnw: Unitless::new(0.),
                    maxlnw_index: 0,
                    num_states: 1,
                    t_last: 1,
                },
            },
        }
    }
}

impl Bins {
    /// Find the index corresponding to a given energy.  This should
    /// panic if the energy is less than `min`.
    pub fn state_to_index(&self, s: State) -> usize {
        s.N
    }
    /// Find the energy corresponding to a given index.
    pub fn index_to_state(&self, i: usize) -> State {
        State {
            E: Energy::new(0.),
            N: i,
        }
    }
    /// Make room in our arrays for a new energy value.  Returns true
    /// if we had to make changes.
    pub fn prepare_for_state(&mut self, e: State) -> bool {
        let mut made_change = false;
        while e.N + 1 > self.lnw.len() {
            self.lnw.push(Unitless::new(0.));
            self.total_energy.push(Energy::new(0.));
            self.total_energy_squared.push(EnergySquared::new(0.));
            self.histogram.push(0);
            self.num_translation_accepted.push(0);
            self.num_translation_attempts.push(0);
            self.round_trips.push(1);
            self.t_found.push(0);
            self.have_visited_since_maxlnw.push(false);
            let last = self.translation_scale.pop().unwrap();
            self.translation_scale.push(last);
            self.translation_scale.push(last);
            made_change = true;
        }
        made_change
    }
}

impl<S: GrandSystem> NumberMC<S> {
    /// Find the index corresponding to a given state.
    pub fn state_to_index(&self, s: State) -> usize {
        self.bins.state_to_index(s)
    }
    /// Find the state corresponding to a given index.
    pub fn index_to_state(&self, i: usize) -> State {
        self.bins.index_to_state(i)
    }

    /// This decides whether to reject the move based on the actual
    /// method in use.
    fn reject_move(&mut self, e1: State, e2: State) -> bool {
        let i1 = self.state_to_index(e1);
        let i2 = self.state_to_index(e2);
        let lnw1 = *(self.bins.lnw[i1] + e1.E / self.T).value();
        let lnw2 = *(self.bins.lnw[i2] + e2.E / self.T).value();
        match self.method {
            Method::Samc { .. } => lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp(),
            Method::Sad {
                too_many, mu_max, ..
            } => {
                let lnw = &self.bins.lnw;
                let lnw1 = if e1.N > too_many {
                    let max_slope = *(mu_max / self.T).value();
                    *(lnw[too_many] + (e1.N - too_many) as f64 / max_slope).value()
                } else {
                    lnw1
                };
                let lnw2 = if e2.N > too_many {
                    let max_slope = *(mu_max / self.T).value();
                    *(lnw[too_many] + (e2.N - too_many) as f64 / max_slope).value()
                } else {
                    lnw2
                };
                lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp()
            }
            Method::WL {
                ref mut bins,
                ref mut gamma,
                ref mut lowest_hist,
                ref mut highest_hist,
                ref mut total_hist,
                ..
            } => {
                if bins.prepare_for_state(e2) {
                    // Oops, we found a new energy, so let's regroup.
                    if *gamma != 1.0 || bins.histogram.len() == 0 {
                        self.movies.new_gamma(self.moves, *gamma);
                        self.movies.new_gamma(self.moves, 1.0);
                        println!(
                            "    WL: Starting fresh with {} numbers",
                            self.bins.lnw.len()
                        );
                        *gamma = 1.0;
                        *lowest_hist = 0;
                        *highest_hist = 0;
                        *total_hist = 0;
                    } else {
                        // We can keep our counts!
                    }
                }
                let rejected = lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp();
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
        let n = i; // number of atoms is also the index.
        match self.method {
            Method::Samc { .. } => {}
            Method::Sad {
                mu_max,
                ref mut too_many,
                ref mut tF,
                ref mut highest_hist,
                ref mut latest_parameter,
                ..
            } => {
                let histogram = &self.bins.histogram;
                if n > *too_many {
                    // Ooops, we didn't want to add gamma after all...
                    self.bins.lnw[i] = old_lnw;
                }

                if histogram[i] > *highest_hist {
                    *highest_hist = histogram[i];
                    if n > *too_many {
                        let ihi = *too_many;
                        for j in ihi + 1..histogram.len() {
                            let lnw = &mut self.bins.lnw;
                            lnw[j] = if histogram[j] != 0 {
                                lnw[ihi] + (j - ihi) as f64 * mu_max / self.T
                            } else {
                                Unitless::new(0.0)
                            };
                        }
                        *latest_parameter = *(*too_many as f64 * mu_max / self.T).value();
                        *too_many = n;

                        gamma_changed = true;
                        // Set tF to the latest discovery time in the
                        // range of energies that we actually care about.
                        *tF = *self.bins.t_found[0..n + 1].iter().max().unwrap();
                        // We just discovered a new important number of atoms.
                        // Should we take this as an opportunity to
                        // revise our translation scale? We should definitely
                        // log the news.
                        if !self.report.quiet {
                            println!(
                                "    sad: [{}]:  0 < {} < {}",
                                self.moves,
                                *too_many,
                                self.bins.histogram.len() + 1
                            );
                        }
                    }
                }
            }
            Method::WL {
                ref mut gamma,
                ref mut lowest_hist,
                ref mut highest_hist,
                ref mut total_hist,
                ref mut bins,
            } => {
                let num_states = self.bins.num_states as f64;
                bins.histogram[i] += 1;
                if bins.histogram[i] > *highest_hist {
                    *highest_hist = bins.histogram[i];
                }
                *total_hist += 1;
                let histogram = &self.bins.histogram;
                if bins.histogram[i] == *lowest_hist + 1
                    && bins
                        .histogram
                        .iter()
                        .enumerate()
                        .filter(|(i, _)| histogram[*i] != 0)
                        .map(|(_, &h)| h)
                        .min()
                        == Some(*lowest_hist + 1)
                {
                    *lowest_hist = bins.histogram[i];
                    if bins.histogram.len() > 1
                        && *lowest_hist as f64 >= 0.8 * *total_hist as f64 / num_states
                    {
                        gamma_changed = true;
                        *gamma *= 0.5;
                        if *gamma > 1e-16 {
                            println!(
                                "    WL:  We have reached flatness {:.2} with min {}!",
                                PrettyFloat(*lowest_hist as f64 * num_states / *total_hist as f64),
                                *lowest_hist
                            );
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
            self.movies.new_gamma(self.moves, gamma);
            self.movies.new_gamma(self.moves, self.gamma());
        }
    }
}

impl<S> NumberMC<S> {
    fn gamma(&self) -> f64 {
        match self.method {
            Method::Samc { t0 } => {
                let t = self.moves as f64;
                if t > t0 {
                    t0 / t
                } else {
                    1.0
                }
            }
            Method::Sad {
                latest_parameter,
                too_many,
                tF,
                ..
            } => {
                let t = self.moves as f64;
                let tF = tF as f64;
                let too_many = too_many as f64;
                if latest_parameter == 0.0 {
                    0.0
                } else {
                    (latest_parameter + t / tF) / (latest_parameter + t / too_many * (t / tF))
                }
            }
            Method::WL { gamma, .. } => gamma,
        }
    }
}

impl<S: GrandSystem + serde::Serialize + serde::de::DeserializeOwned> MonteCarlo for NumberMC<S> {
    type Params = NumberMCParams;
    type System = S;
    fn from_params(params: NumberMCParams, system: S, save_as: ::std::path::PathBuf) -> Self {
        NumberMC {
            method: Method::new(params._method),
            moves: 0,
            T: params.T,
            accepted_moves: 0,
            max_N: params.max_N,
            bins: Bins {
                histogram: vec![1],
                lnw: vec![Unitless::new(0.0)],
                total_energy: vec![Energy::new(0.0)],
                total_energy_squared: vec![EnergySquared::new(0.0)],
                translation_scale: vec![match params._moves {
                    MoveParams::_Explicit {
                        translation_scale, ..
                    } => translation_scale,
                    _ => 0.05 * units::SIGMA,
                }],
                num_translation_attempts: vec![0],
                num_translation_accepted: vec![0],
                t_found: vec![0],
                num_addremove_attempts: 0,
                num_addremove_accepted: 0,
                have_visited_since_maxlnw: vec![false],
                round_trips: vec![1],
                maxlnw: Unitless::new(0.),
                maxlnw_index: 0,
                num_states: 1,
                t_last: 1,
            },
            move_plan: params._moves,
            addremove_probability: {
                if let MoveParams::_Explicit {
                    addremove_probability,
                    ..
                } = params._moves
                {
                    addremove_probability
                } else {
                    // This is the case where we will dynamically
                    // adjust this probability.  Start out very often
                    // adding or removing, since we begin with zero
                    // atoms.  Also because adding and removing is the
                    // most "random" way to change the system, so
                    // modulo rejection we should prefer to
                    // add/remove.
                    0.5
                }
            },
            system: system,

            rng: crate::rng::MyRng::seed_from_u64(params.seed.unwrap_or(0)),
            save_as: save_as,
            report: plugin::Report::from(params._report),
            movies: Movies::from(params.movie),
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

        if self.rng.gen::<f64>() > self.addremove_probability {
            self.bins.num_translation_attempts[e1.N] += 1;
            if let Some(e2) = self
                .system
                .plan_move(&mut self.rng, self.bins.translation_scale[e1.N])
            {
                let e2 = State { E: e2, N: e1.N };
                self.bins.prepare_for_state(e2);
                if !self.reject_move(e1, e2) {
                    self.accepted_moves += 1;
                    self.bins.num_translation_accepted[e1.N] += 1;
                    self.system.confirm();
                }
            }
        } else {
            self.bins.num_addremove_attempts += 1;
            let option_e2 = if self.rng.gen::<u64>() & 1 == 0 {
                // add
                if let Some(max_N) = self.max_N {
                    if e1.N == max_N {
                        None
                    } else {
                        self.system
                            .plan_add(&mut self.rng)
                            .map(|e2| State { E: e2, N: e1.N + 1 })
                    }
                } else {
                    self.system
                        .plan_add(&mut self.rng)
                        .map(|e2| State { E: e2, N: e1.N + 1 })
                }
            } else {
                // remove
                if e1.N > 0 {
                    Some(State {
                        E: self.system.plan_remove(&mut self.rng),
                        N: e1.N - 1,
                    })
                } else {
                    None
                }
            };
            if let Some(e2) = option_e2 {
                self.bins.prepare_for_state(e2);
                if !self.reject_move(e1, e2) {
                    self.accepted_moves += 1;
                    self.bins.num_addremove_accepted += 1;
                    self.system.confirm();
                }
            }
        }
        let energy = State::new(&self.system);
        let i = self.state_to_index(energy);

        if self.bins.histogram[i] == 0 {
            self.bins.num_states += 1;
            self.bins.t_last = self.moves;
            self.bins.t_found[i] = self.moves;
        }
        if self.moves % self.bins.t_last == 0 {
            if let MoveParams::AcceptanceRate(r) = self.move_plan {
                let old_addremove_prob = self.addremove_probability;
                let overall_acceptance_rate = self.accepted_moves as f64 / self.moves as f64;
                let addremove_acceptance_rate = self.bins.num_addremove_accepted as f64
                    / self.bins.num_addremove_attempts as f64;
                let translation_goal = if addremove_acceptance_rate > r {
                    // addremove has too high an acceptance rate, so
                    // we want translations to have too low an
                    // acceptance rate.
                    r * r
                } else {
                    // addremove has too low an acceptance rate, so we
                    // want translations to have a higher acceptance
                    // rate.
                    r.sqrt()
                };
                let translation_acceptance_rate =
                    (self.accepted_moves - self.bins.num_addremove_accepted) as f64
                        / (self.moves - self.bins.num_addremove_attempts) as f64;
                // Fine-tune the translation scales.
                for i in 2..self.bins.lnw.len() {
                    if self.bins.num_translation_attempts[i] == 0 {
                        continue;
                    }
                    let old_scale = self.bins.translation_scale[i];
                    let my_rate = self.bins.num_translation_accepted[i] as f64
                        / self.bins.num_translation_attempts[i] as f64;
                    let s = (my_rate / translation_goal).sqrt();
                    let s = if s < 0.8 {
                        0.8
                    } else if s > 1.2 {
                        1.2
                    } else {
                        s
                    };
                    self.bins.translation_scale[i] *= s;
                    if self.bins.translation_scale[i] > self.system.max_size() {
                        self.bins.translation_scale[i] = self.system.max_size();
                    }
                    if !self.report.quiet && self.bins.translation_scale[i] != old_scale {
                        println!(
                            "        new translation scale for N={}: {:.3}",
                            i, self.bins.translation_scale[i]
                        );
                        println!("          acceptance rate {:.1}%", 100.0 * my_rate);
                    }
                }
                // Now adjust the addremove probability to try to get
                // the net acceptance rate to match our goal "r".
                let new_p = (r - translation_acceptance_rate)
                    / (addremove_acceptance_rate - translation_acceptance_rate);
                if !new_p.is_nan() {
                    if new_p > 0.99 {
                        self.addremove_probability = 0.99;
                    } else if new_p < 0.01 {
                        self.addremove_probability = 0.01;
                    } else {
                        self.addremove_probability = new_p;
                    }
                }
                if !self.report.quiet
                    && (self.addremove_probability - old_addremove_prob).abs() > 0.01
                {
                    println!(
                        "        new addremove_probability: {:.2}",
                        self.addremove_probability
                    );
                    println!(
                        "        acceptance rate {:.1}% [add/remove: {:.1}%]",
                        100.0 * overall_acceptance_rate,
                        100.0 * addremove_acceptance_rate
                    );
                }
            }
        }
        self.bins.histogram[i] += 1;
        self.bins.total_energy[i] += energy.E;
        self.bins.total_energy_squared[i] += energy.E * energy.E;
        self.update_weights(energy);

        if self.bins.lnw[i] > self.bins.maxlnw {
            self.bins.maxlnw = self.bins.lnw[i];
            self.bins.maxlnw_index = i;
            for x in self.bins.have_visited_since_maxlnw.iter_mut() {
                *x = true;
            }
        } else if i == self.bins.maxlnw_index {
            if self.state_to_index(e1) != i {
                for x in self.bins.have_visited_since_maxlnw.iter_mut() {
                    *x = false;
                }
            }
        } else if !self.bins.have_visited_since_maxlnw[i] {
            self.bins.have_visited_since_maxlnw[i] = true;
            self.bins.round_trips[i] += 1;
        }

        let plugins = [&self.report as &dyn Plugin<Self>, &self.movies, &self.save];
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
#[derive(AutoArgs, Debug)]
pub struct MoviesParams {
    // How often (logarithmically) do we want a movie frame? If this
    // is 2.0, it means we want a frame every time the number of
    // iterations doubles.
    /// 2.0 means a frame every time iterations double.
    pub time: Option<f64>,
}

impl Default for MoviesParams {
    fn default() -> Self {
        MoviesParams { time: None }
    }
}

/// A plugin that saves movie data.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Movies {
    lnw: RefCell<Vec<Vec<f64>>>,
    histogram: RefCell<Vec<Vec<u64>>>,

    movie_time: Option<f64>,
    which_frame: Cell<i32>,
    period: Cell<plugin::TimeToRun>,
    time: RefCell<Vec<u64>>,
    gamma: RefCell<Vec<f64>>,
    gamma_time: RefCell<Vec<u64>>,
}

impl From<MoviesParams> for Movies {
    fn from(params: MoviesParams) -> Self {
        Movies {
            movie_time: params.time,
            lnw: RefCell::new(Vec::new()),
            histogram: RefCell::new(Vec::new()),
            which_frame: Cell::new(0),
            period: Cell::new(if params.time.is_some() {
                plugin::TimeToRun::TotalMoves(1)
            } else {
                plugin::TimeToRun::Never
            }),
            time: RefCell::new(Vec::new()),
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
impl<S: GrandSystem + serde::Serialize + serde::de::DeserializeOwned> Plugin<NumberMC<S>>
    for Movies
{
    fn run(&self, mc: &NumberMC<S>, sys: &S) -> plugin::Action {
        if let Some(movie_time) = self.movie_time {
            let moves = mc.num_moves();
            if plugin::TimeToRun::TotalMoves(moves) == self.period.get() {
                println!("Saving movie...");
                assert_eq!(sys.energy(), sys.compute_energy());
                assert!(sys.energy() <= Energy::new(0.0));
                self.new_gamma(moves, mc.gamma());
                // First, let's create the arrays for the time and
                // energy indices.
                self.time.borrow_mut().push(moves);

                let histogram = mc.bins.histogram.clone();
                let lnw: Vec<_> = mc.bins.lnw.iter().map(|v| *v.value()).collect();
                for val in self.lnw.borrow_mut().iter_mut() {
                    while val.len() < lnw.len() {
                        val.push(0.);
                    }
                }
                for h in self.histogram.borrow_mut().iter_mut() {
                    while h.len() < lnw.len() {
                        h.push(0);
                    }
                }

                self.lnw.borrow_mut().push(lnw);
                let mut hist_movie = self.histogram.borrow_mut();
                hist_movie.push(histogram);

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
    fn run_period(&self) -> plugin::TimeToRun {
        self.period.get()
    }

    /// This isn't really a movies thing, but there isn't a great
    /// reason to create yet another plugin for an NumberMC-specific
    /// log message.
    fn log(&self, mc: &NumberMC<S>, sys: &S) {
        let mut one_trip: Option<State> = None;
        let mut ten_trips: Option<State> = None;
        let mut hundred_trips: Option<State> = None;
        let mut thousand_trips: Option<State> = None;
        for (i, &trips) in mc.bins.round_trips.iter().enumerate().rev() {
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
        let thousand_trips = thousand_trips
            .map(|e| format!("{}", e))
            .unwrap_or("-".to_string());
        let ten_trips = ten_trips
            .map(|e| format!("{}", e))
            .unwrap_or("-".to_string());
        let one_trip = one_trip
            .map(|e| format!("{}", e))
            .unwrap_or("-".to_string());
        let hundred_trips = hundred_trips
            .map(|e| format!("{}", e))
            .unwrap_or("-".to_string());
        if !mc.report.quiet {
            println!(
                "   {} * {} * {} * {} | {} currently {} {{{} {:.2}}}",
                one_trip,
                ten_trips,
                hundred_trips,
                thousand_trips,
                mc.index_to_state(mc.bins.maxlnw_index),
                State::new(sys),
                mc.bins.num_states,
                PrettyFloat(mc.moves as f64 / mc.bins.t_last as f64),
            );
            if let Method::WL {
                lowest_hist,
                highest_hist,
                total_hist,
                ref bins,
                ..
            } = mc.method
            {
                if total_hist > 0 {
                    let mut lowest = 123456789;
                    let mut highest = 123456789;
                    for (i, &h) in bins.histogram.iter().enumerate() {
                        if h == lowest_hist && mc.bins.histogram[i] != 0 {
                            lowest = i;
                        }
                        if h == highest_hist && mc.bins.histogram[i] != 0 {
                            highest = i;
                        }
                    }
                    println!(
                        "        WL:  flatness {:.1} with min {:.2} at {} and max {:.2} at {}!",
                        PrettyFloat(
                            lowest_hist as f64 * mc.bins.num_states as f64 / total_hist as f64
                        ),
                        PrettyFloat(lowest_hist as f64),
                        mc.index_to_state(lowest),
                        PrettyFloat(highest_hist as f64),
                        mc.index_to_state(highest)
                    );
                }
            }
        }
    }
}
