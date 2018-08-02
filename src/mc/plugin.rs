//! A plugin architecture to enable reusing of interfaces and
//! implementation for different Monte Carlo algorithms.

use super::*;

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
