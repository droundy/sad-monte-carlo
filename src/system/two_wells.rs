//! A test model with fake energy.

use super::*;

use rand::prelude::*;
/// The parameters needed to configure an fake model with N dimensions.
///
/// These parameters are normally set via command-line arguments.
#[derive(Serialize, Deserialize, Debug, AutoArgs, Clone)]
#[allow(non_snake_case)]
pub struct Parameters {
    /// the mean energy
    pub mean_energy: Energy,
    /// the number of dimensions
    pub N: usize,
}

#[allow(non_snake_case)]
/// An Fake model.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TwoWells {
    /// The state of the system
    pub position: Vec<f64>,
    /// The function itself
    pub parameters: Parameters,
    /// The last change we made (and might want to undo).
    possible_change: Vec<f64>,
}

impl From<Parameters> for TwoWells {
    fn from(parameters: Parameters) -> TwoWells {
        TwoWells {
            position: vec![0.5; parameters.N],
            possible_change: Vec::new(),
            parameters,
        }
    }
}

impl TwoWells {
    fn find_energy(&self, position: &[f64]) -> Energy {
        Energy::new(
            position
                .iter()
                .map(|&x| self .parameters.mean_energy.value_unsafe+statrs::function::erf::erf_inv(x))
                .sum::<f64>(),
        )
    }
}

impl System for TwoWells {
    fn energy(&self) -> Energy {
        self.find_energy(&self.position)
    }
    fn compute_energy(&self) -> Energy {
        self.energy()
    }
    fn randomize(&mut self, rng: &mut MyRng) -> Energy {
        for x in self.position.iter_mut() {
            *x = rng.gen_range(-1.0, 1.0);
        }
        self.energy()
    }
    fn min_moves_to_randomize(&self) -> u64 {
        self.dimensionality() // FIXME /3
    }
    fn dimensionality(&self) -> u64 {
        self.position.len() as u64
    }
}

impl ConfirmSystem for TwoWells {
    fn confirm(&mut self) {
        self.position = self.possible_change.clone();
    }
}

impl MovableSystem for TwoWells {
    fn plan_move(&mut self, rng: &mut MyRng, d: Length) -> Option<Energy> {
        // FIXME move three coordinates at once.
        let i = rng.gen_range(0, self.position.len());
        self.possible_change = self.position.clone();
        let v: f64 = rng.sample(rand_distr::StandardNormal);
        self.possible_change[i] += v * d.value_unsafe;
        if self.possible_change[i] >= 1.0 || self.possible_change[i] <= -1.0 {
            return None;
        }
        Some(self.find_energy(&self.possible_change))
    }
    fn max_size(&self) -> Length {
        Length::new(2.0)
    }
}