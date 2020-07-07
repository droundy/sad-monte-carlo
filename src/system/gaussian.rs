//! A test model with fake energy.

use super::*;

use rand::prelude::*;

/// The parameters needed to configure an fake model.
///
/// These parameters are normally set via command-line arguments.
#[derive(Serialize, Deserialize, Debug, AutoArgs)]
#[allow(non_snake_case)]
pub struct Parameters {
    /// the mean energy
    pub mean_energy: Energy,
}

#[allow(non_snake_case)]
/// An Fake model.
#[derive(Serialize, Deserialize, Debug)]
pub struct Gaussian {
    /// The state of the system
    pub position: Vec<f64>,
    /// The function itself
    pub parameters: Parameters,
    /// The last change we made (and might want to undo).
    possible_change: Vec<f64>,
}

impl From<Parameters> for Gaussian {
    fn from(parameters: Parameters) -> Gaussian {
        Gaussian {
            position: Vec::new(),
            possible_change: Vec::new(),
            parameters: parameters,
        }
    }
}

impl Gaussian {
    fn find_energy(&self, position: &[f64]) -> Energy {
        position
            .iter()
            .map(|&x| Energy::new(statrs::function::erf::erf_inv(x)))
            .sum::<Energy>()
    }
}

impl System for Gaussian {
    type CollectedData = ();
    fn energy(&self) -> Energy {
        self.find_energy(&self.position)
    }
    fn compute_energy(&self) -> Energy {
        self.energy()
    }
}

impl ConfirmSystem for Gaussian {
    fn confirm(&mut self) {
        self.position = self.possible_change.clone();
    }
}

impl MovableSystem for Gaussian {
    fn plan_move(&mut self, rng: &mut MyRng, d: Length) -> Option<Energy> {
        let i = rng.gen_range(0, self.position.len());
        self.possible_change = self.position.clone();
        let v: f64 = rng.sample(rand_distr::StandardNormal);
        self.possible_change[i] += v * d.value_unsafe;
        Some(self.find_energy(&self.possible_change))
    }
}
