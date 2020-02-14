//! A trait and implementations for storing binned information.

#![allow(non_snake_case)]

use crate::system::{units,Energy,PerEnergy};

use dimensioned::Dimensionless;
use std::default::Default;

/// The `Binning` trait defines a type that can hold histogram-like information
pub trait Binning : Default {
    /// Increment the counts stored for a measurement of energy `e`
    /// with update factor `gamma`.
    fn increment_counts(&mut self, e: Energy, gamma: f64);
    /// Find the current estimate of the entropy, `lnw`.
    ///
    /// Note that any energy value is permissible, and should not
    /// result in a panic.
    fn get_lnw(&self, e: Energy) -> f64;
    /// Find the current estimate of the number of counts per energy.
    ///
    /// This is not just an integer in order to abstract out the bin
    /// size.
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

/// A standard energy histogram
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

#[test]
fn test_bins() {
    test_binning::<Bins>();
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
        *((e - self.min)/self.width).value() as usize
    }
}

impl Binning for Bins {
    fn increment_counts(&mut self, e: Energy, gamma: f64) {
        assert!(self.width > Energy::new(0.0));
        assert_eq!(self.lnw.len(), self.histogram.len());
        while e < self.min {
            // this is a little wasteful, but seems the easiest way to
            // ensure we end up with enough room.
            self.histogram.insert(0, 0);
            self.lnw.insert(0, 0.0);
            self.min -= self.width;
        }
        while e >= self.min + self.width*(self.lnw.len() as f64) {
            self.lnw.push(0.0);
            self.histogram.push(0);
        }
        let idx = self.energy_to_index(e);
        self.histogram[idx] += 1;
        self.lnw[idx] += gamma;
    }
    fn get_lnw(&self, e: Energy) -> f64 {
        self.lnw[self.energy_to_index(e)]
    }
    fn get_counts(&self, e: Energy) -> PerEnergy {
        let idx = self.energy_to_index(e);
        if idx < self.histogram.len() {
            self.histogram[idx] as f64/self.width
        } else {
            0.0/units::EPSILON
        }
    }
}
