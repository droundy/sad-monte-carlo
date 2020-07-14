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
pub struct ErfInv {
    /// The state of the system
    pub position: Vec<f64>,
    /// The function itself
    pub parameters: Parameters,
    /// The last change we made (and might want to undo).
    possible_change: Vec<f64>,
}

impl From<Parameters> for ErfInv {
    fn from(parameters: Parameters) -> ErfInv {
        ErfInv {
            position: Vec::new(),
            possible_change: Vec::new(),
            parameters: parameters,
        }
    }
}

impl ErfInv {
    fn find_energy(&self, position: &[f64]) -> Energy {
        Energy::new(
            position
                .iter()
                .map(|&x| self .parameters.mean_energy.value_unsafe+statrs::function::erf::erf_inv(x))
                .sum::<f64>(),
        )
    }
}

impl System for ErfInv {
    type CollectedData = ();
    fn energy(&self) -> Energy {
        self.find_energy(&self.position)
    }
    fn compute_energy(&self) -> Energy {
        self.energy()
    }
}

impl ConfirmSystem for ErfInv {
    fn confirm(&mut self) {
        self.position = self.possible_change.clone();
    }
}

impl MovableSystem for ErfInv {
    fn plan_move(&mut self, rng: &mut MyRng, d: Length) -> Option<Energy> {
        let i = rng.gen_range(0, self.position.len());
        self.possible_change = self.position.clone();
        let v: f64 = rng.sample(rand_distr::StandardNormal);
        self.possible_change[i] += v * d.value_unsafe;
        Some(self.find_energy(&self.possible_change))
    }
}

impl GrandSystem for ErfInv {
    fn plan_add(&mut self, rng: &mut MyRng) -> Option<Energy> {
        self.possible_change = self.position.clone();
        self.possible_change.push( rng.gen_range(0., 1.));
        Some(self.find_energy(&self.possible_change))
    }
    fn plan_remove(&mut self, rng: &mut MyRng) -> Energy {
        let i = rng.gen_range(0, self.position.len());
        self.possible_change = self.position.clone();
        self.possible_change.swap_remove(i);
        self.find_energy(&self.possible_change)
    }
    fn num_atoms(&self) -> usize {
        self.position.len()
    }
}
