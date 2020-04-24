//! A histogram binning solution

#![allow(non_snake_case)]

use crate::system::{units,Energy,PerEnergy};
use crate::mc::binning::{ Binning, Interned };

use dimensioned::Dimensionless;
use std::default::Default;

/// A set of counts for a variable
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct BinCounts {
    /// The total of the thing
    total: Vec<f64>,
    /// The least value of the thing
    min_total: f64,
    /// The greatest value of the thing
    max_total: f64,
    /// The energy corresponding to the greatest value of the thing
    e_max_total: Energy,
    /// The count in each bin
    count: Vec<u64>,
    /// The lowest count
    min_count: u64,
    /// The greatest count
    max_count: u64,
    /// The energy that has the greatest count
    e_max_count: Energy,
    /// The total count
    total_count: u64,
}

impl BinCounts {
    fn new(sz: usize) -> Self {
        BinCounts {
            total: vec![0.; sz],
            total_count: 0,
            min_total: 0.,
            max_total: 0.,
            e_max_total: Energy::new(-std::f64::INFINITY),
            count: vec![0; sz],
            min_count: 0,
            max_count: 0,
            e_max_count: Energy::new(-std::f64::INFINITY),
        }
    }
    fn insert_zero(&mut self) {
        self.total.insert(0, 0.);
        self.count.insert(0, 0);
        self.min_total = 0.;
        self.min_count = 0;
    }
    fn push_zero(&mut self) {
        self.total.push(0.);
        self.count.push(0);
        self.min_total = 0.;
        self.min_count = 0;
    }
    fn increment_count(&mut self, e: Energy, idx: usize, value: f64) {
        self.total_count += 1;
        let old_count = self.count[idx];
        self.count[idx] += 1;
        let old_total = self.total[idx];
        self.total[idx] += value;
        if self.total[idx] > self.max_total
            && max_of(&self.total) == self.total[idx]
        {
            self.max_total = self.total[idx];
            self.e_max_total = e;
        }
        if old_total == self.min_total {
            self.min_total = min_of(&self.total);
        }
        if old_count == self.min_count {
            self.min_count = self.count.iter().cloned().min().unwrap();
        }
        if self.count[idx] > self.max_count {
            self.max_count = self.count[idx];
            self.e_max_count = e;
        }
    }

    fn get_total(&self, idx: usize) -> f64 {
        if idx < self.total.len() {
            self.total[idx]
        } else {
            0.0
        }
    }
    fn get_count(&self, idx: usize) -> u64 {
        if idx < self.count.len() {
            self.count[idx]
        } else {
            0
        }
    }
}

/// A standard energy histogram
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct Bins {
    /// The lowest allowed energy in any bin.
    pub min: Energy,
    min_e: Energy,
    max_e: Energy,
    /// The energy bin size.
    pub width: Energy,
    /// The ln weight for each energy bin.
    pub lnw: BinCounts,
    /// The extra totals
    pub extra: std::collections::HashMap<Interned, BinCounts>,
}

#[test]
fn test_histogram() {
    crate::mc::binning::test_binning::<Bins>();
}

impl Default for Bins {
    fn default() -> Self {
        Bins {
            min: Energy::new(std::f64::INFINITY),
            width: units::EPSILON,
            lnw: BinCounts::new(0),
            extra: std::collections::HashMap::new(),
            min_e: Energy::new(std::f64::INFINITY),
            max_e: Energy::new(-std::f64::INFINITY),
        }
    }
}

impl Bins {
    fn index_to_energy(&self, i: usize) -> Energy {
        self.min + (i as f64 + 0.5)*self.width
    }
    fn energy_to_index(&self, e: Energy) -> usize {
        if e < self.min {
            std::usize::MAX // this is a bogus value for negative-index cases
        } else {
            *((e - self.min)/self.width).value() as usize
        }
    }
    fn prep_for_e(&mut self, e: Energy) {
        assert!(self.width > Energy::new(0.0));
        while e < self.min {
            // this is a little wasteful, but seems the easiest way to
            // ensure we end up with enough room.
            self.lnw.insert_zero();
            for d in self.extra.values_mut() {
                d.insert_zero();
            }
            self.min -= self.width;
        }
        while e >= self.min + self.width*(self.lnw.count.len() as f64) {
            self.lnw.push_zero();
            for d in self.extra.values_mut() {
                d.push_zero();
            }
        }
    }
}

impl Binning for Bins {
    fn new(e: Energy, width: Energy) -> Self {
        let min = ((e/width).value().round() - 0.5)*width;
        Bins {
            min,
            width,
            lnw: BinCounts::new(0),
            extra: std::collections::HashMap::new(),
            min_e: e,
            max_e: e,
        }
    }
    fn increment_count(&mut self, e: Energy, gamma: f64) {
        if e > self.max_e {
            self.max_e = e;
        }
        if e < self.min_e {
            self.min_e = e;
        }
        self.prep_for_e(e);
        let idx = self.energy_to_index(e);
        self.lnw.increment_count(e, idx, gamma);
    }
    fn set_lnw<F: Fn(Energy, PerEnergy) -> Option<f64>>(&mut self, f: F) {
        for i in 0..self.lnw.count.len() {
            let e = self.index_to_energy(i);
            if let Some(lnw) = f(e, self.get_count(e)) {
                self.lnw.total[i] = lnw;
                self.lnw.count[i] = 0;
            }
        }
    }
    fn count_states<F: Fn(Energy, PerEnergy) -> bool>(&self, f:F) -> usize {
        let mut total = 0;
        for i in 0..self.lnw.count.len() {
            let e = self.index_to_energy(i);
            if f(e, self.get_count(e)) {
                total += 1;
            }
        }
        total
    }
    fn num_states(&self) -> usize {
        self.lnw.count.len()
    }

    fn get_lnw(&self, e: Energy) -> f64 {
        self.lnw.get_total(self.energy_to_index(e))
    }
    fn get_count(&self, e: Energy) -> PerEnergy {
        let idx = self.energy_to_index(e);
        self.lnw.get_count(idx) as f64/self.width
    }
    fn max_lnw(&self) -> f64 {
        self.lnw.max_total
    }
    fn min_lnw(&self) -> f64 {
        self.lnw.min_total
    }
    fn max_count(&self) -> PerEnergy {
        self.lnw.max_count as f64/self.width
    }
    fn min_count(&self) -> PerEnergy {
        self.lnw.min_count as f64/self.width
    }

    fn accumulate_extra(&mut self, name: Interned, e: Energy, value: f64) {
        self.prep_for_e(e);
        let idx = self.energy_to_index(e);
        assert!(idx < self.lnw.total.len());
        if let Some(data) = self.extra.get_mut(&name) {
            data.total_count += 1;
            data.count[idx] += 1;
            if data.count[idx] == data.min_count + 1 {
                data.min_count = data.count.iter().cloned().min().unwrap();
            }
            if data.count[idx] > data.max_count {
                data.max_count = data.count[idx];
            }
            let old_total = data.total[idx];
            data.total[idx] = old_total + value;
            if data.total[idx] > data.max_total {
                data.max_total = data.total[idx];
            }
            if old_total == data.min_total {
                data.min_total = min_of(&data.total);
            }
        } else {
            let values = BinCounts::new(self.lnw.total.len());
            self.extra.insert(name, values);
            self.accumulate_extra(name, e, value); // sloppy recursion...
        }
    }
    fn zero_out_extra(&mut self, name: Interned) {
        if let Some(data) = self.extra.get_mut(&name) {
            for v in data.count.iter_mut() {
                *v = 0;
            }
            for v in data.total.iter_mut() {
                *v = 0.;
            }
            data.min_total = 0.;
            data.max_total = -std::f64::INFINITY;
            data.min_count = 0;
            data.max_count = 0;
            data.total_count = 0;
        }
    }
    fn mean_extra(&self, name: Interned, e: Energy) -> f64 {
        if let Some(data) = self.extra.get(&name) {
            let idx = self.energy_to_index(e);
            if data.count[idx] > 0 {
                data.get_total(idx) / data.count[idx] as f64
            } else {
                0.0
            }
        } else {
            0.0
        }
    }
    fn total_extra(&self, name: Interned, e: Energy) -> f64 {
        if let Some(data) = self.extra.get(&name) {
            let idx = self.energy_to_index(e);
            data.get_total(idx)
        } else {
            0.0
        }
    }
    fn max_total_extra(&self, name: Interned) -> f64 {
        if let Some(data) = self.extra.get(&name) {
            data.max_total
        } else {
            0.0
        }
    }
    fn min_total_extra(&self, name: Interned) -> f64 {
        if let Some(data) = self.extra.get(&name) {
            data.min_total
        } else {
            0.0
        }
    }
    fn count_extra(&self, name: Interned, e: Energy) -> PerEnergy {
        if let Some(data) = self.extra.get(&name) {
            let idx = self.energy_to_index(e);
            data.get_count(idx) as f64 / self.width
        } else {
            PerEnergy::new(0.0)
        }
    }
    fn total_count(&self) -> u64 {
        self.lnw.total_count
    }
    fn total_count_extra(&self, name: Interned) -> u64 {
        if let Some(data) = self.extra.get(&name) {
            data.total_count
        } else {
            0
        }
    }
    fn mean_count_extra(&self, name: Interned) -> PerEnergy {
        if let Some(data) = self.extra.get(&name) {
            data.total_count as f64/(self.width*self.num_states() as f64)
        } else {
            PerEnergy::new(0.)
        }
    }
    fn min_count_extra(&self, name: Interned) -> PerEnergy {
        if let Some(data) = self.extra.get(&name) {
            data.min_count as f64/self.width
        } else {
            PerEnergy::new(0.)
        }
    }

    fn max_energy(&self) -> Energy {
        self.max_e
    }
    fn min_energy(&self) -> Energy {
        self.min_e
    }
}

fn max_of(stuff: &[f64]) -> f64 {
    stuff.iter().cloned().fold(0./0., f64::max)
}

fn min_of(stuff: &[f64]) -> f64 {
    stuff.iter().cloned().fold(0./0., f64::min)
}
