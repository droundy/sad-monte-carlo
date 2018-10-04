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

/// Parameters to configure a particular MC.
#[derive(Debug, ClapMe)]
pub enum MethodParams {
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
pub struct NumberMCParams {
    /// The actual method.
    pub _method: MethodParams,
    /// The seed for the random number generator.
    pub seed: Option<u64>,
    /// The temperature.
    T: Energy,
    _moves: MoveParams,
    _report: plugin::ReportParams,
    _movies: MoviesParams,
    _save: plugin::SaveParams,
}

impl Default for NumberMCParams {
    fn default() -> Self {
        NumberMCParams {
            _method: MethodParams::Samc { t0: 1e3 },
            seed: None,
            T: 1.0*units::EPSILON,
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
    /// The number of atoms
    pub N: usize,
}
impl ::std::fmt::Display for State {
    fn fmt(&self, f: &mut ::std::fmt::Formatter) -> ::std::fmt::Result {
        write!(f, "{}", self.N)
    }
}
impl State {
    /// Find the state of a system
    pub fn new<S: GrandSystem>(sys: &S) -> State {
        State { N: sys.num_atoms() }
    }
}

/// This defines the energy bins.
#[derive(Serialize, Deserialize, Debug)]
pub struct Bins {
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
pub struct NumberMC<S> {
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
    /// The temperature.
    pub T: Energy,
    report: plugin::Report,
    movies: Movies,
    save: plugin::Save,
    manager: plugin::PluginManager,
}

#[derive(Serialize, Deserialize, Debug)]
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
        num_states: f64,
        hist: Vec<u64>,
    }
}

impl Method {
    fn new(p: MethodParams) -> Self {
        match p {
            MethodParams::Samc { t0 } => Method::Samc { t0 },
            MethodParams::WL => Method::WL {
                gamma: 1.0,
                lowest_hist: 1,
                highest_hist: 1,
                total_hist: 0,
                num_states: 1.0,
                hist: Vec::new(),
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
        State { N: i }
    }
    /// Make room in our arrays for a new energy value
    pub fn prepare_for_state(&mut self, e: State) {
        let n = e.N;
        while n >= self.lnw.len() {
            self.lnw.push(Unitless::new(0.0));
            self.histogram.push(0);
            self.have_visited_since_maxentropy.push(true);
            self.round_trips.push(1);
        }
    }
}

impl<S: System> NumberMC<S> {
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
    fn reject_move(&mut self, n1: usize, e1: Energy, n2: usize, e2: Energy) -> bool {
        let i1 = self.state_to_index(State { N: n1 });
        let i2 = self.state_to_index(State { N: n2 });
        let lnw1 = *self.bins.lnw[i1].value() + e1/self.T;
        let lnw2 = *self.bins.lnw[i2].value() + e2/self.T;
        match self.method {
            Method::Samc { .. } => {
                lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp()
            }
            Method::WL { ref mut num_states, .. } => {
                let rejected = lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp();
                if !rejected && self.bins.histogram[i2] == 0 {
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
        self.bins.lnw[i] += gamma;
        let mut gamma_changed = false;
        match self.method {
            Method::Samc { .. } => {}
            Method::WL { ref mut gamma,
                         ref mut lowest_hist, ref mut highest_hist, ref mut total_hist,
                         ref mut hist, num_states } => {
                if hist.len() != self.bins.lnw.len() {
                    // Oops, we found a new number, so let's regroup.
                    if *gamma != 1.0 || hist.len() == 0 {
                        println!("    WL: Starting fresh with {} energies",
                                 self.bins.lnw.len());
                        gamma_changed = true;
                        *gamma = 1.0;
                        *lowest_hist = 0;
                        *highest_hist = 0;
                        *total_hist = 0;
                        *hist = vec![0; self.bins.lnw.len()];
                    } else {
                        // We have to adjust our hist Vec, but can
                        // keep our counts! We just pretend we already
                        // knew there were more energies to be
                        // found...
                        while hist.len() < self.bins.lnw.len() {
                            hist.push(0);
                        }
                        *lowest_hist = 0;
                    }
                }
                hist[i] += 1;
                if hist[i] > *highest_hist {
                    *highest_hist = hist[i];
                }
                *total_hist += 1;
                let histogram = &self.bins.histogram;
                if hist[i] == *lowest_hist + 1
                    && hist.len() > 1
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
                            println!("         gamma => {}", PrettyFloat(*gamma));
                        }
                        hist.iter_mut().map(|x| *x = 0).count();
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
                if t > t0 { t0/t } else { 1.0 }
            }
            Method::WL { gamma, .. } => {
                gamma
            }
        }
    }
}

impl<S: GrandSystem> MonteCarlo for NumberMC<S> {
    type Params = NumberMCParams;
    type System = S;
    fn from_params(params: NumberMCParams, system: S, save_as: ::std::path::PathBuf) -> Self {
        NumberMC {
            method: Method::new(params._method),
            moves: 0,
            time_L: 0,
            accepted_moves: 0,
            acceptance_rate: 0.5, // arbitrary starting guess.
            T: params.T,
            bins: Bins {
                histogram: vec![1],
                lnw: vec![Unitless::new(0.0)],
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
        let e1 = self.system.energy();
        let n1 = self.system.num_atoms();
        let i1 = self.state_to_index(State { N: n1 });
        let recent_scale = (1.0/self.moves as f64).sqrt();
        self.acceptance_rate *= 1. - recent_scale;
        if let Some(e2) = self.system.plan_move(&mut self.rng, self.translation_scale) {
            if !self.reject_move(n1, e1, n1, e2) {
                self.accepted_moves += 1;
                self.acceptance_rate += recent_scale;
                self.system.confirm();
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
            if i1 != i {
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
impl<S: GrandSystem> Plugin<NumberMC<S>> for Movies {
    fn run(&self, mc: &NumberMC<S>, _sys: &S) -> plugin::Action {
        if let Some(movie_time) = self.movie_time {
            let moves = mc.num_moves();
            if plugin::TimeToRun::TotalMoves(moves) == self.period.get() {
                println!("Saving movie...");
                self.new_gamma(moves, mc.gamma());
                // First, let's create the arrays for the time and
                // energy indices.
                self.time.borrow_mut().push(moves);

                let histogram = mc.bins.histogram.clone();
                let entropy: Vec<_> = mc.bins.lnw.iter().map(|v| *v.value()).collect();
                for entr in self.entropy.borrow_mut().iter_mut() {
                    while entr.len() < entropy.len() {
                        entr.push(0.);
                    }
                }
                for h in self.histogram.borrow_mut().iter_mut() {
                    while h.len() < entropy.len() {
                        h.push(0);
                    }
                }
                let mut S = self.entropy.borrow_mut();
                S.push(entropy);
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
    fn run_period(&self) -> plugin::TimeToRun { self.period.get() }

    /// This isn't really a movies thing, but there isn't a great
    /// reason to create yet another plugin for an NumberMC-specific
    /// log message.
    fn log(&self, mc: &NumberMC<S>, sys: &S) {
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
        let thousand_trips = thousand_trips.map(|e| format!("{}", e))
            .unwrap_or("-".to_string());
        let ten_trips = ten_trips.map(|e| format!("{}", e))
            .unwrap_or("-".to_string());
        let one_trip = one_trip.map(|e| format!("{}", e))
            .unwrap_or("-".to_string());
        let hundred_trips = hundred_trips.map(|e| format!("{}", e))
            .unwrap_or("-".to_string());
        if !mc.report.quiet {
            println!("   {} * {} * {} * {} | {} currently {}",
                     one_trip, ten_trips, hundred_trips, thousand_trips,
                     mc.index_to_state(mc.bins.max_S_index).N,
                     sys.energy()/units::EPSILON,
            );
            if let Method::WL { lowest_hist, highest_hist, total_hist, num_states,
                                ref hist, .. } = mc.method {
                let mut lowest = 111111111;
                let mut highest = 111111111;
                for (i,&h) in hist.iter().enumerate() {
                    if h == lowest_hist && mc.bins.histogram[i] != 0 {
                        lowest = i;
                    }
                    if h == highest_hist && mc.bins.histogram[i] != 0 {
                        highest = i;
                    }
                }
                println!("        WL:  flatness {:.1} with min {:.2} at {} and max {:.2} at {}!",
                         PrettyFloat(lowest_hist as f64*num_states as f64
                                     / total_hist as f64),
                         PrettyFloat(lowest_hist as f64),
                         mc.index_to_state(lowest),
                         PrettyFloat(highest_hist as f64),
                         mc.index_to_state(highest));
            }
        }
    }
}
