//! A trait and implementations for storing binned information.

#![allow(non_snake_case)]

use crate::system::{units,Energy,PerEnergy};

use dimensioned::Dimensionless;
use std::default::Default;

trait Binning : Default {
    fn increment_counts(&mut self, e: Energy, gamma: f64);
    fn get_lnw(&self, e: Energy) -> f64;
    fn get_counts(&self, e: Energy) -> PerEnergy;
}

#[cfg(test)]
fn test_binning<B: Binning>() {
    let eps = units::EPSILON;
    let mut b = B::default();
    assert_eq!(b.get_counts(eps), 0.0/eps);
    b.increment_counts(eps, 1.0);
    assert!(b.get_counts(eps) > 0.0/eps);
}

/// Where we store the info about the energy grid
#[derive(Serialize, Deserialize, Debug)]
pub struct Bins {
    /// The lowest allowed energy in any bin.
    pub min: Energy,
    /// The energy bin size.
    pub width: Energy,
    /// The number of times we have been at each energy.
    pub histogram: Vec<u64>,
    // /// The iteration when we found each energy.
    // pub t_found: Vec<u64>,
    /// The ln weight for each energy bin.
    pub lnw: Vec<f64>,
}

impl Default for Bins {
    fn default() -> Self {
        Bins {
            min: units::EPSILON*0.0,
            width: units::EPSILON,
            histogram: Vec::new(),
            lnw: Vec::new(),
        }
    }
}

impl Bins {
    fn index_to_energy(&self, i: usize) -> Energy {
        self.min + (i as f64 + 0.5)*self.width
    }
    fn energy_to_index(&self, e: Energy) -> usize {
        let i = *((e - self.min)/self.width).value() as usize;
        i
    }
}

impl Binning for Bins {
    fn increment_counts(&mut self, e: Energy, gamma: f64) {
    }
    fn get_lnw(&self, e: Energy) -> f64 {
        1.0
    }
    fn get_counts(&self, e: Energy) -> PerEnergy {
        1.0/units::EPSILON
    }
}
