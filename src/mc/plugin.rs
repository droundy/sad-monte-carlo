//! A plugin architecture to enable reusing of interfaces and
//! implementation for different Monte Carlo algorithms.

use super::*;
use std::cell::Cell;

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
    fn run(&self, mc: &MC, sys: &MC::System) -> Action { Action::None }
    /// How often we need the plugin to run.  A `None` value means
    /// that this plugin never needs to run.  Note that it is expected
    /// that this period may change any time the plugin is called, so
    /// this should be a cheap call as it may happen frequently.  Also
    /// note that this is an upper, not a lower bound.
    fn run_period(&self) -> Option<u64> { None }
    /// We might be about to die, so please do any cleanup or saving.
    /// Note that the plugin state is stored on each checkpoint.  This
    /// is called in response to `Action::Save` and `Action::Exit`.
    fn save(&self, _mc: &MC, _sys: &MC::System) {}
    /// Log to stdout any interesting data we think our user might
    /// care about.  This is called in response to `Action::Save`,
    /// `Action::Log` and `Action::Exit`.
    fn log(&self, _mc: &MC, _sys: &MC::System) {}
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
    period: Cell<u64>,
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
                if let Some(period) = p.run_period() {
                    if period < new_period {
                        new_period = period;
                    }
                }
            }
            self.period.set(new_period);
        }
    }
}

/// A plugin that terminates the simulation after a fixed number of iterations.
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct MaxIter {
    max_iter: Option<u64>,
}

/// The parameter to define the maximum number of iterations.
#[derive(ClapMe, Debug)]
pub struct MaxIterParams {
    max_iter: Option<u64>,
}

impl From<MaxIterParams> for MaxIter {
    fn from(params: MaxIterParams) -> Self {
        MaxIter { max_iter: params.max_iter }
    }
}
impl<MC: MonteCarlo> Plugin<MC> for MaxIter {
    fn run(&self, mc: &MC, _sys: &MC::System) -> Action {
        if let Some(maxiter) = self.max_iter {
            if mc.num_moves() >= maxiter {
                return Action::Exit;
            }
        }
        Action::None
    }
    fn run_period(&self) -> Option<u64> { self.max_iter }
}

/// A plugin that terminates the simulation after a fixed number of iterations.
#[derive(Serialize, Deserialize, Debug, Clone, Copy)]
pub struct FinalReport;

/// The parameter to define the maximum number of iterations.
#[derive(ClapMe, Debug)]
pub struct FinalReportParams;

impl From<FinalReportParams> for FinalReport {
    fn from(params: FinalReportParams) -> Self {
        FinalReport
    }
}
impl<MC: MonteCarlo> Plugin<MC> for FinalReport {
    fn save(&self, mc: &MC, _sys: &MC::System) {
        let rejects = mc.num_rejected_moves();
        let moves = mc.num_moves();
        println!("Rejected {}/{} = {}% of the moves",
                 rejects, moves, 100.0*rejects as f64/moves as f64);
    }
}
