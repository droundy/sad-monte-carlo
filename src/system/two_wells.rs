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
    /// Height of the first well
    pub h_1: Energy,
    /// Height of the second well
    pub h_2: Energy,
    /// Radius of the first well
    pub r_1: Length,
    /// Radius of the first well
    pub r_2: Length,
    /// Center of the Second Well
    pub center_2: Vec<f64>,
}

#[allow(non_snake_case)]
/// A Fake model.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct TwoWells {
    /// The state of the system
    pub position: Vec<Length>,
    /// The function itself
    pub parameters: Parameters,
    /// The last change we made (and might want to undo).
    possible_change: Vec<Length>,
}

impl From<Parameters> for TwoWells {
    fn from(parameters: Parameters) -> TwoWells {
        TwoWells {
            position: vec![Length::new(0.5); parameters.N],
            possible_change: Vec::new(),
            parameters,
        }
    }
}

impl TwoWells {
    fn find_energy(&self, position: &[Length]) -> Energy {
        self.parameters.h_1
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
            *x = Length::new(rng.gen_range(-1.0, 1.0));
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
        self.possible_change[i] += v * d;
        if self.possible_change[i] >= Length::new(1.0) || self.possible_change[i] <= Length::new(-1.0) {
            return None;
        }
        Some(self.find_energy(&self.possible_change))
    }
    fn max_size(&self) -> Length {
        Length::new(2.0)
    }
}

#[test]
fn a_test() {
    println!("I am testing");
    assert!(true);
}