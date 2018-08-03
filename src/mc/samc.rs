//! An implementation of SAD (Statistical Association, Dynamical version).

use ::system::*;
use super::*;

use super::plugin::Plugin;
use dimensioned::Dimensionless;
use rand::Rng;
use std::default::Default;

/// The parameters needed to configure a SAD simulation.
#[allow(non_snake_case)]
#[derive(Debug, ClapMe)]
pub struct SamcParams {
    /// The t0 parameter.
    pub t0: u64,
    /// The seed for the random number generator.
    pub seed: Option<u64>,
    _maxiter: plugin::MaxIterParams,
    _final_report: plugin::FinalReportParams,
}

impl Default for SamcParams {
    fn default() -> Self {
        SamcParams {
            t0: 1000,
            seed: None,
            _maxiter: plugin::MaxIterParams::default(),
            _final_report: plugin::FinalReportParams::default(),
        }
    }
}

#[allow(non_snake_case)]
/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct Samc<S> {
    /// The system we are simulating.
    pub system: S,
    /// The minimum temperature we are interested in.
    pub t0: u64,
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
    maxiter: plugin::MaxIter,
    final_report: plugin::FinalReport,
    manager: plugin::PluginManager,
}

impl<S: System> Samc<S> {
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
}

impl<S: MovableSystem> MonteCarlo for Samc<S> {
    type Params = SamcParams;
    type System = S;
    fn from_params(params: SamcParams, system: S, save_as: ::std::path::PathBuf) -> Self {
        Samc {
            t0: params.t0,
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
            maxiter: plugin::MaxIter::from(params._maxiter),
            final_report: plugin::FinalReport::from(params._final_report),
            manager: plugin::PluginManager::new(),
        }
    }

    #[allow(non_snake_case)]
    fn move_once(&mut self) {
        self.moves += 1;
        let e1 = self.system.energy();
        if let Some(_) = self.system.move_once(&mut self.rng, Length::new(0.1)) {
            let e2 = self.system.energy();
            self.prepare_for_energy(e2);
            let i1 = self.energy_to_index(e1);
            let i2 = self.energy_to_index(e2);

            let lnw1 = self.lnw[i1].value();
            let lnw2 = self.lnw[i2].value();
            if lnw2 > lnw1 && self.rng.gen::<f64>() > (lnw1 - lnw2).exp() {
                // we reject this move!
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
        let t = self.moves;
        let gamma = if t > self.t0 { self.t0 as f64/t as f64 } else { 1.0 };
        self.lnw[i] += gamma;

        if self.lnw[i] > self.max_S {
            self.max_S = self.lnw[i];
            self.max_entropy_energy = energy;
        }
        let plugins = [&self.maxiter as &Plugin<Self>,
                       &self.final_report,
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
