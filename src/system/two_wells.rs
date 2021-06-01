//! A test model with fake energy.

use super::*;

use rand::prelude::*;
/// The parameters needed to configure a fake model with N dimensions.
///
/// These parameters are normally set via command-line arguments.
#[derive(Serialize, Deserialize, Debug, AutoArgs, Clone)] //AutoArgs in incompatable with Vec<Length>
#[allow(non_snake_case)]
pub struct Parameters {
    /// the number of dimensions
    pub N: usize,
    /// Ratio of depths of the wells (h_2/h_1)
    pub h2_to_h1: f64,
    /// Radius of the second well
    pub r2: f64,
}

#[allow(non_snake_case)]
/// A Fake model.
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
        //r_1 is assumed to be 1 and center_2 is such that the two wells are touching

        let mut d_1_squared = 0.;
        let mut d_2_squared = 0.;

        let center_2 = (1.+self.parameters.r2) / ((self.parameters.N as f64).sqrt());

        for i in 0..(self.parameters.N){
            d_1_squared += (position[i])*(position[i]);
            d_2_squared += (position[i] - center_2)*(position[i] - center_2);
        }

        let e_1 = Energy::new(1.)*(d_1_squared - 1.);
        let e_2 = Energy::new(self.parameters.h2_to_h1)*(d_2_squared/ (self.parameters.r2*self.parameters.r2) - 1.);

        if e_1 < e_2{
            e_1
        }
        else{
            e_2
        }
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

#[test]
fn a_test() {
    println!("I am testing");
    assert!(true);
}