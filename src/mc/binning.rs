//! A trait and implementations for storing binned information.

#![allow(non_snake_case)]

use crate::system::{units,Energy,PerEnergy};

use dimensioned::Dimensionless;
use std::default::Default;
use internment::Intern;
use auto_args::AutoArgs;

/// Parameters to decide how to do binning.
#[derive(Debug, AutoArgs)]
#[allow(non_camel_case_types)]
pub enum BinningParams {
    /// Use an ordinary histogram.
    Histogram {
        /// A histogram with this bin width
        bin: Energy,
    },
}

impl Default for BinningParams {
    fn default() -> Self {
        BinningParams::Histogram { bin: Energy::new(1.0) }
    }
}

/// Any Binning thing
#[derive(Serialize, Deserialize, Debug, Clone)]
pub enum Bins {
    /// Just a plain histogram
    Histogram( Histogram ),
}

impl Default for Bins {
    fn default() -> Self {
        Bins::Histogram( Histogram::default() )
    }
}

impl Bins {
    /// Construct a Bins from an energy and a set of command-line flags
    pub fn from_params<S: crate::system::System>(sys: &S, params: BinningParams) -> Self {
        match params {
            BinningParams::Histogram { bin } => {
                Bins::Histogram(Histogram::new(sys.energy(), bin))
            }
        }
    }
}

impl Binning for Bins {
// impl Bins {
    fn num_states(&self) -> usize {
        match self {
            Bins::Histogram(b) => b.num_states(),
        }
    }
    fn count_states<F: Fn(Energy, PerEnergy) -> bool>(&self, f: F) -> usize {
        match self {
            Bins::Histogram(b) => b.count_states(f),
        }
    }

    fn increment_count(&mut self, e: Energy, gamma: f64) {
        match self {
            Bins::Histogram(b) => b.increment_count(e, gamma),
        }
    }

    fn get_lnw(&self, e: Energy) -> f64 {
        match self {
            Bins::Histogram(b) => b.get_lnw(e),
        }
    }
    fn get_count(&self, e: Energy) -> PerEnergy {
        match self {
            Bins::Histogram(b) => b.get_count(e),
        }
    }

    fn set_lnw<F: Fn(Energy, PerEnergy) -> Option<f64>>(&mut self, f: F) {
        match self {
            Bins::Histogram(b) => b.set_lnw(f),
        }
    }
    fn max_lnw(&self) -> f64 {
        match self {
            Bins::Histogram(b) => b.max_lnw(),
        }
    }
    fn min_lnw(&self) -> f64 {
        match self {
            Bins::Histogram(b) => b.min_lnw(),
        }
    }

    fn total_count(&self) -> u64 {
        match self {
            Bins::Histogram(b) => b.total_count(),
        }
    }
    fn max_count(&self) -> PerEnergy {
        match self {
            Bins::Histogram(b) => b.max_count(),
        }
    }
    fn min_count(&self) -> PerEnergy {
        match self {
            Bins::Histogram(b) => b.min_count(),
        }
    }

    fn accumulate_extra(&mut self, name: Intern<String>, e: Energy, value: f64) {
        match self {
            Bins::Histogram(b) => b.accumulate_extra(name, e, value),
        }
    }
    fn zero_out_extra(&mut self, name: Intern<String>) {
        match self {
            Bins::Histogram(b) => b.zero_out_extra(name),
        }
    }
    fn total_extra(&self, name: Intern<String>, e: Energy) -> f64 {
        match self {
            Bins::Histogram(b) => b.total_extra(name, e),
        }
    }
    fn mean_count_extra(&self, extra: Intern<String>) -> PerEnergy {
        match self {
            Bins::Histogram(b) => b.mean_count_extra(extra),
        }
    }

    fn mean_extra(&self, name: Intern<String>, e: Energy) -> f64 {
        match self {
            Bins::Histogram(b) => b.mean_extra(name, e),
        }
    }

    fn min_energy(&self) -> Energy {
        match self {
            Bins::Histogram(b) => b.min_energy(),
        }
    }
    fn max_energy(&self) -> Energy {
        match self {
            Bins::Histogram(b) => b.max_energy(),
        }
    }

    fn min_total_extra(&self, name: Intern<String>) -> f64 {
        match self {
            Bins::Histogram(b) => b.min_total_extra(name),
        }
    }
    fn max_total_extra(&self, name: Intern<String>) -> f64 {
        match self {
            Bins::Histogram(b) => b.max_total_extra(name),
        }
    }
    fn min_count_extra(&self, name: Intern<String>) -> PerEnergy {
        match self {
            Bins::Histogram(b) => b.min_count_extra(name),
        }
    }
    fn total_count_extra(&self, name: Intern<String>) -> u64 {
        match self {
            Bins::Histogram(b) => b.total_count_extra(name),
        }
    }
    fn count_extra(&self, name: Intern<String>, e: Energy) -> PerEnergy {
        match self {
            Bins::Histogram(b) => b.count_extra(name, e),
        }
    }
    fn new(e: Energy, de: Energy) -> Self {
        Bins::Histogram(Histogram::new(e,de))
    }
}

#[test]
fn test_bins() {
    test_binning::<Bins>();
}

/// The `Binning` trait defines a type that can hold histogram-like information
pub trait Binning : Default + serde::Serialize + serde::de::DeserializeOwned+ std::fmt::Debug + Clone {
    /// The number of times we must call increment_count in order to
    /// uniformly shift the lnw.
    fn num_states(&self) -> usize;
    /// The number of times we must call increment_count in order to
    /// uniformly shift the lnw.
    fn count_states<F: Fn(Energy, PerEnergy) -> bool>(&self, f: F) -> usize;

    /// Increment the counts stored for a measurement of energy `e`
    /// with update factor `gamma`.
    fn increment_count(&mut self, e: Energy, gamma: f64);

    /// Find the current estimate of the entropy, `lnw`.
    ///
    /// Note that any energy value is permissible, and should not
    /// result in a panic.
    fn get_lnw(&self, e: Energy) -> f64;
    /// Find the current estimate of the number of counts per energy.
    ///
    /// This is not just an integer in order to abstract out the bin
    /// size.
    fn get_count(&self, e: Energy) -> PerEnergy;

    /// Set the lnw
    fn set_lnw<F: Fn(Energy, PerEnergy) -> Option<f64>>(&mut self, f: F);
    /// Find the maximum of the entropy, `lnw`.
    fn max_lnw(&self) -> f64;
    /// Find the minimum of the entropy, `lnw`.
    fn min_lnw(&self) -> f64;

    /// Total counts in the histogram.
    fn total_count(&self) -> u64;
    /// Find the maximum of the histogram.
    fn max_count(&self) -> PerEnergy;
    /// Find the minimum of the histogram.
    fn min_count(&self) -> PerEnergy;

    /// Accumulate some other datum that we might want to keep track
    /// of.  This is for things like pressure for instance, and we do
    /// not assume that they are incremented with each count, so they
    /// need their own count.  Note that this does require storing a
    /// HashMap of extra data, but that seems like the only really
    /// flexible way of doing this.
    fn accumulate_extra(&mut self, name: Intern<String>, e: Energy, value: f64);
    /// Reset the accumulated extra data.
    fn zero_out_extra(&mut self, name: Intern<String>);
    /// Find the total accumulated extra value that we accumulated at
    /// this energy.
    fn total_extra(&self, name: Intern<String>, e: Energy) -> f64;
    /// Find the average extra value that we accumulated at this
    /// energy.
    fn mean_extra(&self, name: Intern<String>, e: Energy) -> f64;
    /// Find the maximum extra value.
    fn max_total_extra(&self, name: Intern<String>) -> f64;
    /// Find the minimum extra value.
    fn min_total_extra(&self, name: Intern<String>) -> f64;
    /// Find the minimum extra count.
    fn min_count_extra(&self, name: Intern<String>) -> PerEnergy;
    /// Find the number of times we have counted this thing
    fn count_extra(&self, name: Intern<String>, e: Energy) -> PerEnergy;
    /// Find the number of times we have counted this thing
    fn total_count_extra(&self, name: Intern<String>) -> u64;
    /// Find the mean value of count_extra
    fn mean_count_extra(&self, name: Intern<String>) -> PerEnergy;
    /// Creates a Binning with the desired energy and energy width.
    ///
    /// It may be a no-op for a binning scheme that does not require a
    /// width.
    fn new(e: Energy, de: Energy) -> Self;
    /// The highest energy that has been seen
    fn max_energy(&self) -> Energy;
    /// The lowest energy that has been seen
    fn min_energy(&self) -> Energy;
}

#[cfg(test)]
fn test_binning<B: Binning>() {
    let eps = units::EPSILON;
    let mut b = B::default();
    assert_eq!(b.get_count(eps), 0.0/eps);
    b.increment_count(eps, 1.0);
    assert!(b.get_count(eps) > 0.0/eps);
    let mydat = Intern::new("datum".to_string());
    b.accumulate_extra(mydat, eps, 7.0);
    assert_eq!(b.mean_extra(mydat, eps), 7.0);
    assert_eq!(b.total_extra(mydat, eps), 7.0);
    assert_eq!(b.max_total_extra(mydat), 7.0);
}

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
pub struct Histogram {
    /// The lowest allowed energy in any bin.
    pub min: Energy,
    min_e: Energy,
    max_e: Energy,
    /// The energy bin size.
    pub width: Energy,
    /// The ln weight for each energy bin.
    pub lnw: BinCounts,
    /// The extra totals
    pub extra: std::collections::HashMap<Intern<String>, BinCounts>,
}

#[test]
fn test_histogram() {
    test_binning::<Histogram>();
}

impl Default for Histogram {
    fn default() -> Self {
        Histogram {
            min: Energy::new(std::f64::INFINITY),
            width: units::EPSILON,
            lnw: BinCounts::new(0),
            extra: std::collections::HashMap::new(),
            min_e: Energy::new(std::f64::INFINITY),
            max_e: Energy::new(-std::f64::INFINITY),
        }
    }
}

impl Histogram {
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

impl Binning for Histogram {
    fn new(e: Energy, width: Energy) -> Self {
        let min = ((e/width).value().round() - 0.5)*width;
        Histogram {
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

    fn accumulate_extra(&mut self, name: Intern<String>, e: Energy, value: f64) {
        self.prep_for_e(e);
        let idx = self.energy_to_index(e);
        if let Some(data) = self.extra.get_mut(&name) {
            data.count[idx] += 1;
            data.total[idx] += value;
        } else {
            let values = BinCounts::new(self.lnw.total.len());
            self.extra.insert(name, values);
            self.accumulate_extra(name, e, value); // sloppy recursion...
        }
    }
    fn zero_out_extra(&mut self, name: Intern<String>) {
        if let Some(data) = self.extra.get_mut(&name) {
            for v in data.count.iter_mut() {
                *v = 0;
            }
            for v in data.total.iter_mut() {
                *v = 0.;
            }
            data.min_total = std::f64::INFINITY;
            data.max_total = -std::f64::INFINITY;
            data.min_count = 0;
            data.max_count = 0;
            data.total_count = 0;
        }
    }
    fn mean_extra(&self, name: Intern<String>, e: Energy) -> f64 {
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
    fn total_extra(&self, name: Intern<String>, e: Energy) -> f64 {
        if let Some(data) = self.extra.get(&name) {
            let idx = self.energy_to_index(e);
            data.get_total(idx)
        } else {
            0.0
        }
    }
    fn max_total_extra(&self, name: Intern<String>) -> f64 {
        if let Some(data) = self.extra.get(&name) {
            data.max_total
        } else {
            0.0
        }
    }
    fn min_total_extra(&self, name: Intern<String>) -> f64 {
        if let Some(data) = self.extra.get(&name) {
            data.min_total
        } else {
            0.0
        }
    }
    fn count_extra(&self, name: Intern<String>, e: Energy) -> PerEnergy {
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
    fn total_count_extra(&self, name: Intern<String>) -> u64 {
        if let Some(data) = self.extra.get(&name) {
            data.total_count
        } else {
            0
        }
    }
    fn mean_count_extra(&self, name: Intern<String>) -> PerEnergy {
        if let Some(data) = self.extra.get(&name) {
            data.total_count as f64/(self.width*self.num_states() as f64)
        } else {
            PerEnergy::new(0.)
        }
    }
    fn min_count_extra(&self, name: Intern<String>) -> PerEnergy {
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
