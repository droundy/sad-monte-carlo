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
    count: Vec<f64>,
    /// The lowest count
    min_count: f64,
    /// The greatest count
    max_count: f64,
    /// The energy that has the greatest count
    e_max_count: Energy,
    /// The total count
    total_count: u64,
}

#[derive(Debug,Clone,Copy)]
enum Idx {
    None,
    One(usize, f64),
    Two(usize, f64),
}

impl Iterator for Idx {
    type Item = (usize, f64);
    fn next(&mut self) -> Option<(usize, f64)> {
        match *self {
            Idx::None => None,
            Idx::One(i,o) => {
                *self = Idx::None;
                Some((i, 1.-o))
            }
            Idx::Two(i,o) => {
                *self = Idx::One(i,o);
                Some((i+1, o))
            }
        }
    }
}

impl BinCounts {
    fn new(sz: usize) -> Self {
        BinCounts {
            total: vec![0.; sz],
            total_count: 0,
            min_total: 0.,
            max_total: 0.,
            e_max_total: Energy::new(-std::f64::INFINITY),
            count: vec![0.; sz],
            min_count: 0.,
            max_count: 0.,
            e_max_count: Energy::new(-std::f64::INFINITY),
        }
    }
    fn insert_zero(&mut self) {
        self.total.insert(0, 0.);
        self.count.insert(0, 0.);
        self.min_total = 0.;
        self.min_count = 0.;
    }
    fn push_zero(&mut self) {
        self.total.push(0.);
        self.count.push(0.);
        self.min_total = 0.;
        self.min_count = 0.;
    }
    fn interpret_float_index(&self, fidx: f64) -> Idx {
        if fidx < -1.0 {
            Idx::None
        } else if fidx < 0. {
            Idx::One(0, -fidx)
        } else if fidx < self.total.len() as f64 - 1. {
            let i = fidx as usize;
            Idx::Two(i, fidx - i as f64)
        } else if fidx < self.total.len() as f64 {
            Idx::One(self.total.len() - 1, fidx - (self.total.len() as f64 - 1.))
        } else {
            Idx::None
        }
    }
    fn increment_count(&mut self, elo: Energy, ehi: Energy, fidx: f64, value: f64) {
        let idx = fidx as usize;
        let offset = fidx - idx as f64;

        self.total_count += 1;
        let old_count = self.count[idx];
        let old_plus_count = self.count[idx+1];
        self.count[idx] += 1. - offset;
        self.count[idx+1] += offset;

        let old_total = self.total[idx];
        let old_plus_total = self.total[idx+1];
        self.total[idx] += value*(1.-offset);
        self.total[idx+1] += value*offset;
        if self.total[idx] > self.max_total
            && max_of(&self.total) == self.total[idx]
        {
            self.max_total = self.total[idx];
            self.e_max_total = elo;
        }
        if self.total[idx+1] > self.max_total
            && max_of(&self.total) == self.total[idx+1]
        {
            self.max_total = self.total[idx+1];
            self.e_max_total = ehi;
        }
        if old_total == self.min_total || old_plus_total == self.min_total {
            self.min_total = min_of(&self.total);
        }
        if old_count == self.min_count || old_plus_count == self.min_count {
            self.min_count = min_of(&self.count);
        }
        if self.count[idx] > self.max_count {
            self.max_count = self.count[idx];
            self.e_max_count = elo;
        }
        if self.count[idx+1] > self.max_count {
            self.max_count = self.count[idx+1];
            self.e_max_count = ehi;
        }
    }

    fn get_total(&self, fidx: f64) -> f64 {
        let mut total = 0.;
        for (i,f) in self.interpret_float_index(fidx) {
            total += self.total[i]*f;
        }
        total
    }
    fn get_count(&self, fidx: f64) -> f64 {
        let mut total = 0.;
        for (i,f) in self.interpret_float_index(fidx) {
            total += self.count[i]*f;
        }
        total
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
fn test_linear() {
    crate::mc::binning::test_binning::<Bins>();
    println!("passed test_binning");
    let mut bins = Bins::default();
    bins.increment_count(Energy::new(1.9), 1.0);
    assert!(bins.get_lnw(Energy::new(1.95)) > bins.get_lnw(Energy::new(1.85)));
    assert!(bins.get_lnw(Energy::new(1.0)) > 0.0);
    assert!(bins.get_lnw(Energy::new(0.51)) > 0.0);
    assert!(bins.get_lnw(Energy::new(2.49)) > 0.0);
    assert_eq!(bins.get_lnw(Energy::new(3.01)), 0.0);
    for e in -100..100 {
        let e = Energy::new(e as f64*0.125);
        println!("{}: {}", e, bins.get_lnw(e));
    }
    assert_eq!(bins.get_lnw(Energy::new(-0.01)), 0.0);
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
    fn energy_to_index(&self, e: Energy) -> f64 {
        *((e - self.min)/self.width).value()
    }
    fn prep_for_e(&mut self, e: Energy) {
        assert!(self.width > Energy::new(0.0));
        if self.lnw.count.len() == 0 {
            self.min = (e/self.width).value().floor()*self.width;
        }
        while e < self.min {
            // this is a little wasteful, but seems the easiest way to
            // ensure we end up with enough room.
            self.lnw.insert_zero();
            for d in self.extra.values_mut() {
                d.insert_zero();
            }
            self.min -= self.width;
        }
        while e >= self.min + self.width*(self.lnw.count.len() as f64 - 1.0) {
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
        let int_idx = idx as usize;
        let rescaled_gamma = *(gamma*Energy::new(1.)/self.width).value();
        self.lnw.increment_count(self.index_to_energy(int_idx),
                                 self.index_to_energy(int_idx+1), idx, rescaled_gamma);
    }
    fn set_lnw<F: Fn(Energy, PerEnergy) -> Option<f64>>(&mut self, f: F) {
        for i in 0..self.lnw.count.len() {
            let e = self.index_to_energy(i);
            if let Some(lnw) = f(e, self.get_count(e)) {
                self.lnw.total[i] = lnw;
                self.lnw.count[i] = 0.;
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
        self.lnw.get_count(idx)/self.width
    }
    fn max_lnw(&self) -> f64 {
        self.lnw.max_total
    }
    fn min_lnw(&self) -> f64 {
        self.lnw.min_total
    }
    fn max_count(&self) -> PerEnergy {
        self.lnw.max_count/self.width
    }
    fn min_count(&self) -> PerEnergy {
        self.lnw.min_count/self.width
    }

    fn accumulate_extra(&mut self, name: Interned, e: Energy, value: f64) {
        self.prep_for_e(e);
        let idx = self.energy_to_index(e);
        assert!(idx < self.lnw.total.len() as f64);
        if let Some(data) = self.extra.get_mut(&name) {
            let int_idx = idx as usize;
            let offset = idx - int_idx as f64;

            data.total_count += 1;

            let old_count = data.count[int_idx];
            let old_plus_count = data.count[int_idx+1];
            data.count[int_idx] += 1.0 - offset;
            data.count[int_idx+1] += offset;
            if old_count == data.min_count || old_plus_count == data.min_count {
                data.min_count = min_of(&data.count);
            }
            if data.count[int_idx] > data.max_count {
                data.max_count = data.count[int_idx];
            }
            if data.count[int_idx+1] > data.max_count {
                data.max_count = data.count[int_idx+1];
            }

            let old_total = data.total[int_idx];
            let old_plus_total = data.total[int_idx+1];
            data.total[int_idx] = old_total + value*(1.-offset);
            data.total[int_idx+1] = old_plus_total + value*offset;
            if data.total[int_idx] > data.max_total {
                data.max_total = data.total[int_idx];
            }
            if data.total[int_idx+1] > data.max_total {
                data.max_total = data.total[int_idx+1];
            }
            if old_total == data.min_total || old_plus_total == data.min_total {
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
                *v = 0.;
            }
            for v in data.total.iter_mut() {
                *v = 0.;
            }
            data.min_total = 0.;
            data.max_total = -std::f64::INFINITY;
            data.min_count = 0.;
            data.max_count = 0.;
            data.total_count = 0;
        }
    }
    fn mean_extra(&self, name: Interned, e: Energy) -> f64 {
        if let Some(data) = self.extra.get(&name) {
            let idx = self.energy_to_index(e);
            if data.get_count(idx) > 0. {
                data.get_total(idx) / data.get_count(idx)
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
