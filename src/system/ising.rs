//! A square well fluid.

use super::*;

use dimensioned::Dimensionless;
use dimensioned::{Cbrt,Abs};
use vector3d::Vector3d;
use rand::prelude::*;
use rand::distributions::Uniform;
use std::f64::consts::PI;
use std::default::Default;

/// The parameters needed to configure a square well system.
#[derive(Serialize, Deserialize, Debug, ClapMe)]
#[allow(non_snake_case)]
pub struct IsingParams {
    /// Width of the square grid
    N: usize,
}

#[allow(non_snake_case)]
/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct Ising {
    /// The energy of the system
    E: Energy,
    /// The dimensions of the box.
    pub N: usize,
    /// The spins themselves
    S: Vec<i8>,
    /// The last change we made (and might want to undo).
    last_changed: Option<(usize, Energy)>,
}

fn modulus(i: isize, sz: usize) -> usize {
    let sz = sz as isize;
    (((i % sz) + sz) % sz) as usize
}

impl From<IsingParams> for Ising {
    fn from(params: IsingParams) -> Ising {
        let mut ising = Ising {
            E: Energy::new(0.),
            N: params.N,
            S: vec![1; params.N*params.N],
            last_changed: None,
        };
        ising.E = ising.compute_energy();
        ising
    }
}

impl System for Ising {
    fn energy(&self) -> Energy {
        self.E
    }
    fn compute_energy(&self) -> Energy {
        let mut e: Energy = units::EPSILON*0.0;
        for i1 in 0 .. self.N {
            for j1 in 0 .. self.N {
                let j2 = (j1 + 1) % self.N;
                let mut neighbor_tot = self.S[i1 + j2*self.N];
                // let j2 = modulus(j - 1, self.N);
                // neighbor_tot += self.S[i1 + j2*self.N];
                let i2 = (i1 + 1) % self.N;
                neighbor_tot += self.S[i2 + j1*self.N];
                // let i2 = modulus(i - 1, self.N);
                // neighbor_tot += self.S[i1 + j2*self.N];
                e += neighbor_tot as f64 * Energy::new(self.S[i1+j1*self.N] as f64);
            }
        }
        e
    }
    fn delta_energy(&self) -> Option<Energy> {
        Some(2.*units::EPSILON)
    }
}

impl UndoSystem for Ising {
    fn undo(&mut self) {
        if let Some((i, e)) = self.last_changed {
            self.S[i] *= -1;
            self.E = e;
        }
    }
}

impl MovableSystem for Ising {
    fn move_once(&mut self, rng: &mut MyRng, _: Length) -> Option<Energy> {
        Some(self.E)
    }
}
