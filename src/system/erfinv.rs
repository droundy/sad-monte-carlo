//! A test model with fake energy.

use super::*;

use rand::prelude::*;

/// The parameters needed to configure an fake model.
///
/// These parameters are normally set via command-line arguments.
#[derive(Serialize, Deserialize, Debug, AutoArgs, Clone)]
#[allow(non_snake_case)]
pub struct Parameters {
    /// the mean energy
    pub mean_energy: Energy,
}
/// The parameters needed to configure an fake model with N dimensions.
///
/// These parameters are normally set via command-line arguments.
#[derive(Serialize, Deserialize, Debug, AutoArgs, Clone)]
#[allow(non_snake_case)]
pub struct ParametersN {
    /// the fundamental parameters
    pub _parameters: Parameters,
    /// the number of dimensions
    pub N: usize,
}

#[allow(non_snake_case)]
/// An Fake model.
#[derive(Serialize, Deserialize, Debug, Clone)]
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

impl From<ParametersN> for ErfInv {
    fn from(parameters: ParametersN) -> ErfInv {
        ErfInv {
            position: vec![0.5; parameters.N],
            possible_change: Vec::new(),
            parameters: parameters._parameters,
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
        self.position.len() as u64
    }
    fn dimensionality(&self) -> u64 {
        self.min_moves_to_randomize()*3
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
        if self.possible_change[i] >= 1.0 || self.possible_change[i] <= -1.0 {
            return None;
        }
        Some(self.find_energy(&self.possible_change))
    }
    fn max_size(&self) -> Length {
        Length::new(0.5)
    }
}

impl GrandSystem for ErfInv {
    fn plan_add(&mut self, rng: &mut MyRng) -> Option<Energy> {
        self.possible_change = self.position.clone();
        self.possible_change.push( rng.gen_range(-1., 1.));
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

impl GrandReplicaSystem for ErfInv {
    fn plan_swap_atom(&self, other: &Self, rng: &mut MyRng) -> Option<(usize, Energy, Energy)> {
        let which = rng.gen_range(0, self.num_atoms());
        let mut p = self.position.clone();
        let moved = p.swap_remove(which);
        let e_self = self.find_energy(&p);
        p = other.position.clone();
        p.push(moved);
        let e_other = other.find_energy(&p);
        Some((which, e_self, e_other))
    }
    fn swap_atom(&mut self, other: &mut Self, which: usize) {
        let r = self.position.swap_remove(which);
        other.position.push(r);
    }
}