//! An implementation of SAD (Statistical Association, Dynamical version).

use ::system::*;
use super::*;

use dimensioned::Dimensionless;

/// The parameters needed to configure a SAD simulation.
#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug, ClapMe)]
pub struct SadParams {
    /// The minimum temperature we are interested in.
    pub min_T: Energy,
    /// The seed for the random number generator.
    pub seed: u64,
}

#[allow(non_snake_case)]
/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct Sad<S> {
    /// The system we are simulating.
    pub system: S,
    /// The minimum temperature we are interested in.
    pub min_T: Energy,
    /// The number of moves that have been made.
    pub moves: u64,
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
    /// The too-low energy.
    pub too_lo: Energy,
    /// The too-high energy.
    pub too_hi: Energy,
    /// The max-entropy energy.
    pub max_entropy_energy: Energy,
    /// The max-entropy energy.
    pub max_S: Unitless,
    /// The max-entropy energy.
    pub min_important_energy: Energy,
    /// Number of states found.
    pub n_found: u64,

    /// The random number generator.
    pub rng: rng::Rng,
}

impl<S: System> Sad<S> {
    /// Find the index corresponding to a given energy.  This should
    /// panic if the energy is less than `min_energy_bin`.
    pub fn energy_to_index(&self, e: Energy) -> usize {
        *((e - self.min_energy_bin)/self.energy_bin).value() as usize
    }
    /// Find the energy corresponding to a given index.
    pub fn index_to_energy(&self, i: usize) -> Energy {
        self.min_energy_bin + (i as f64)*self.energy_bin
    }
}

impl<S: MovableSystem> MonteCarlo for Sad<S> {
    type Params = SadParams;
    type System = S;
    fn from_params(params: SadParams, system: S) -> Self {
        Sad {
            min_T: params.min_T,
            moves: 0,
            rejected_moves: 0,
            histogram: Vec::new(),
            lnw: Vec::new(),
            n_found: 0,
            min_energy_bin: system.energy(),
            too_lo: system.energy(),
            too_hi: system.energy(),
            max_entropy_energy: system.energy(),
            min_important_energy: system.energy(),
            max_S: Unitless::new(0.0),
            energy_bin: system.delta_energy().unwrap_or(Energy::new(1.0)),
            system: system,

            rng: rng::Rng::from_u64(params.seed),
        }
    }

    fn move_once(&mut self) -> Energy {
        self.system.energy()
    }
    fn num_moves(&self) -> u64 {
        self.moves
    }
    fn num_rejected_moves(&self) -> u64 {
        self.rejected_moves
    }
}
