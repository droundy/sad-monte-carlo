//! A plugin architecture to enable reusing of interfaces and
//! implementation for different Monte Carlo algorithms.

use super::*;

use std::cell::Cell;
use std::default::Default;
use std::time;
use prettyfloat::PrettyFloat;

/// A `Plugin` is an object that can be used to configure a MonteCarlo
/// simulation.  The plugin will be called regularly, and will have a
/// chance to save data (e.g. collect statistics) and/or terminate the
/// simulation.
pub trait Plugin<MC: MonteCarlo> {
    /// Run and do something.  If the simulation needs to be
    /// terminated, `None` is returned.  If you want to modify
    /// information, you will have to use interior mutability, because
    /// I can't figure out any practical way to borrow `self` mutably
    /// while still giving read access to the `MC`.
    fn run(&self, _mc: &MC, _sys: &MC::System) -> Action { Action::None }
    /// How often we need the plugin to run.  A `None` value means
    /// that this plugin never needs to run.  Note that it is expected
    /// that this period may change any time the plugin is called, so
    /// this should be a cheap call as it may happen frequently.  Also
    /// note that this is an upper, not a lower bound.
    fn run_period(&self) -> TimeToRun { TimeToRun::Never }
    /// We might be about to die, so please do any cleanup or saving.
    /// Note that the plugin state is stored on each checkpoint.  This
    /// is called in response to `Action::Save` and `Action::Exit`.
    fn save(&self, _mc: &MC, _sys: &MC::System) {}
    /// Log to stdout any interesting data we think our user might
    /// care about.  This is called in response to `Action::Save`,
    /// `Action::Log` and `Action::Exit`.
    fn log(&self, _mc: &MC, _sys: &MC::System) {}
}

/// A time when we want to be run.
#[derive(Serialize, Deserialize, Copy, Clone, Debug, PartialOrd, Ord, PartialEq, Eq)]
pub enum TimeToRun {
    /// Don't stop on our behalf!
    Never,
    /// After this many moves in total.
    TotalMoves(u64),
    /// This often.
    Period(u64),
}

/// An action that should be taken based on this plugin's decision.
#[derive(Copy, Clone, Debug, PartialOrd, Ord, PartialEq, Eq)]
pub enum Action {
    /// Nothing special need be done.
    None,
    /// Log interesting information.
    Log,
    /// Save things.
    Save,
    /// Exit the program.
    Exit,
}
impl Action {
    /// Do both of two actions.
    pub fn and(self, other: Action) -> Action {
        ::std::cmp::max(self, other)
    }
}

/// A helper to enable Monte Carlo implementations to easily run their
/// plugins without duplicating code.
#[derive(Serialize, Deserialize, Debug)]
pub struct PluginManager {
    #[serde(skip, default)]
    period: Cell<u64>,
    #[serde(skip, default)]
    moves: Cell<u64>,
}

impl PluginManager {
    /// Create a plugin manager.
    pub fn new() -> PluginManager {
        PluginManager { period: Cell::new(1), moves: Cell::new(0) }
    }
    /// Run all the plugins, if needed.  This should always be called
    /// with the same set of plugins.  If you want different sets of
    /// plugins, use different managers.
    pub fn run<MC: MonteCarlo>(&self, mc: &MC, sys: &MC::System,
                               plugins: &[&Plugin<MC>]) {
        let moves = self.moves.get() + 1;
        self.moves.set(moves);
        if moves >= self.period.get() {
            self.moves.set(0);
            let mut todo = plugin::Action::None;
            for p in plugins.iter() {
                todo = todo.and(p.run(mc, sys));
            }
            if todo >= plugin::Action::Log {
                for p in plugins.iter() {
                    p.log(mc, sys);
                }
            }
            if todo >= plugin::Action::Save {
                mc.checkpoint();
                for p in plugins.iter() {
                    p.save(mc, sys);
                }
            }
            if todo >= plugin::Action::Exit {
                ::std::process::exit(0);
            }
            // run plugins every trillion iterations minimum
            let mut new_period = 1u64 << 40;
            for p in plugins.iter() {
                match p.run_period() {
                    TimeToRun::Never => (),
                    TimeToRun::TotalMoves(moves) => {
                        if moves > mc.num_moves() && moves - mc.num_moves() < new_period {
                            new_period = moves - mc.num_moves();
                        }
                    }
                    TimeToRun::Period(period) => {
                        if period < new_period {
                            new_period = period;
                        }
                    }
                }
            }
            self.period.set(new_period);
        }
    }
}

/// A plugin that terminates the simulation after a fixed number of iterations.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Report {
    max_iter: TimeToRun,
    /// This is when and where the simulation started.
    #[serde(skip, default)]
    start: Cell<Option<(time::Instant, u64)>>,
    /// The user has requested that nothing be printed!
    pub quiet: bool,
}

/// The parameters to define the report information as well as stop
/// time (which is part of the report).
#[derive(ClapMe, Debug)]
pub struct ReportParams {
    /// The maximum number of iterations to run.
    pub max_iter: Option<u64>,
    /// Do not make reports!
    pub quiet: bool,
}

impl Default for ReportParams {
    fn default() -> Self {
        ReportParams {
            max_iter: None,
            quiet: true,
        }
    }
}

impl From<ReportParams> for Report {
    fn from(params: ReportParams) -> Self {
        Report {
            max_iter: if let Some(mi) = params.max_iter {
                TimeToRun::TotalMoves(mi)
            } else {
                TimeToRun::Never
            },
            start: Cell::new(Some((time::Instant::now(), 0))),
            quiet: params.quiet,
        }
    }
}
impl Report {
    /// Allows a resuming simulation to get updated report parameters
    /// from the flags.
    pub fn update_from(&mut self, params: ReportParams) {
        let other = Self::from(params);
        self.max_iter = other.max_iter;
        self.quiet = other.quiet;
    }
}
impl<MC: MonteCarlo> Plugin<MC> for Report {
    fn run(&self, mc: &MC, _sys: &MC::System) -> Action {
        if let TimeToRun::TotalMoves(maxiter) = self.max_iter {
            if mc.num_moves() >= maxiter {
                return Action::Exit;
            }
        }
        Action::None
    }
    fn run_period(&self) -> TimeToRun { self.max_iter }
    fn log(&self, mc: &MC, _sys: &MC::System) {
        if self.quiet { return; }
        match self.start.get() {
            Some((start_time, start_iter)) => {
                let moves = mc.num_moves();
                let runtime = start_time.elapsed();
                let time_per_move = duration_to_secs(runtime)/(moves - start_iter) as f64;
                if let TimeToRun::TotalMoves(max) = self.max_iter {
                    let frac_complete = moves as f64/max as f64;
                    let moves_left = if max >= moves { max - moves } else { 0 };
                    let time_left = (time_per_move*moves_left as f64) as u64;
                    println!("[{}] {}% complete after {} ({} left, {:.1}us per move)",
                             PrettyFloat(moves as f64),
                             (100.*frac_complete) as isize,
                             format_duration(runtime.as_secs()),
                             format_duration(time_left),
                             PrettyFloat(time_per_move*1e6));
                } else {
                    println!("[{}] after {} ({:.1}us per move)",
                             PrettyFloat(moves as f64),
                             format_duration(runtime.as_secs()),
                             PrettyFloat(time_per_move*1e6));
                }
            }
            None => {
                self.start.set(Some((time::Instant::now(), mc.num_moves())));
            }
        }
    }
    fn save(&self, mc: &MC, _sys: &MC::System) {
        if self.quiet { return; }
        let accepted = mc.num_accepted_moves();
        let moves = mc.num_moves();
        println!("        Accepted {:.2}/{:.2} = {:.0}% of the moves",
                 PrettyFloat(accepted as f64), PrettyFloat(moves as f64),
                 100.0*accepted as f64/moves as f64);
    }
}


/// A plugin that schedules when to save
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Save {
    #[serde(skip, default)]
    next_output: Cell<u64>,
    /// This is when and where the simulation started.
    #[serde(skip, default)]
    start: Cell<Option<(time::Instant, u64)>>,
    /// How frequently to save...
    #[serde(default)]
    save_time_seconds: Option<f64>,
}

/// The parameter to define the save schedule
#[derive(ClapMe, Debug)]
pub struct SaveParams {
    /// Maximum time between saves in hours
    pub save_time: Option<f64>,
}

impl Default for SaveParams {
    fn default() -> Self {
        SaveParams { save_time: None, }
    }
}
impl From<SaveParams> for Save {
    fn from(params: SaveParams) -> Self {
        Save {
            next_output: Cell::new(1),
            start: Cell::new(Some((time::Instant::now(), 0))),
            save_time_seconds: params.save_time.map(|h| 60.*60.*h),
        }
    }
}
impl Save {
    /// Allows a resuming simulation to get updated save parameters
    /// from the flags.
    pub fn update_from(&mut self, params: SaveParams) {
        self.save_time_seconds = params.save_time.map(|h| 60.*60.*h);
    }
}
impl<MC: MonteCarlo> Plugin<MC> for Save {
    fn run(&self, mc: &MC, _sys: &MC::System) -> Action {
        if mc.num_moves() >= self.next_output.get() {
            Action::Save
        } else {
            Action::None
        }
    }
    fn run_period(&self) -> TimeToRun {
        TimeToRun::TotalMoves(self.next_output.get())
    }
    fn save(&self, mc: &MC, _sys: &MC::System) {
        if let Some(period) = self.save_time_seconds {
            match self.start.get() {
                Some((start_time, start_iter)) => {
                    let moves = mc.num_moves();
                    let runtime = start_time.elapsed();
                    let time_per_move =
                        duration_to_secs(runtime)/(moves - start_iter) as f64;
                    let moves_per_period = 1 + (period/time_per_move) as u64;
                    self.next_output.set(moves + moves_per_period);
                }
                None => {
                    self.start.set(Some((time::Instant::now(), mc.num_moves())));
                    self.next_output.set(mc.num_moves() + (1<<20));
                }
            }
        } else {
            self.next_output.set(self.next_output.get()*2)
        }
    }
}

fn format_duration(secs: u64) -> String {
    let mins = secs / 60;
    let hours = mins / 60;
    let mins = mins % 60;
    if hours > 50 {
        format!("{} hours", hours)
    } else if mins < 1 {
        format!("{} seconds", secs)
    } else if mins == 1 {
        format!("1 minute {} seconds", secs % 60)
    } else if hours < 1 {
        format!("{} minutes", mins)
    } else if hours < 2 {
        format!("1 hour, {} minutes", mins)
    } else {
        format!("{} hours, {} minutes", hours, mins)
    }
}
fn duration_to_secs(t: time::Duration) -> f64 {
    t.as_secs() as f64 + t.subsec_nanos() as f64*1e-9
}
