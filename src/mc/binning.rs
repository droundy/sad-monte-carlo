//! A trait and implementations for storing binned information.

#![allow(non_snake_case)]

use crate::system::{Energy,PerEnergy};

use std::default::Default;
use auto_args::AutoArgs;

pub mod histogram;
pub mod linear;

/// A constant string that we want to use as a parameter name.
#[derive(Hash, Ord, PartialOrd, Eq, PartialEq, Serialize, Deserialize, Clone, Copy)]
pub struct Interned(internment::Intern<std::borrow::Cow<'static,str>>);

impl From<&'static str> for Interned {
    fn from(x: &'static str) -> Self {
        Interned(internment::Intern::new(x.into()))
    }
}
impl From<String> for Interned {
    fn from(x: String) -> Self {
        Interned(internment::Intern::new(x.into()))
    }
}

impl std::fmt::Display for Interned {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        self.0.fmt(f)
    }
}

impl std::fmt::Debug for Interned {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> Result<(), std::fmt::Error> {
        let mystr: &str = self.0.as_ref();
        mystr.fmt(f)
    }
}

#[test]
fn test_interned() {
    let x: Interned = "hello world".into();
    assert_eq!(format!("{}", x), "hello world");
    assert_eq!(format!("{:?}", x), "\"hello world\"");
}

/// Parameters to decide how to do binning.
#[derive(Debug, AutoArgs)]
#[allow(non_camel_case_types)]
pub enum BinningParams {
    /// Use an ordinary histogram.
    Histogram {
        /// A histogram with this bin width
        bin: Energy,
    },
    /// Use a linear interpolation
    Linear {
        /// A regularly spaced linear interpolation with this bin width
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
    Histogram( histogram::Bins ),
    /// A regularly spaced linear interpolation
    Linear( linear::Bins ),
}

impl Default for Bins {
    fn default() -> Self {
        Bins::Histogram( histogram::Bins::default() )
    }
}

impl Bins {
    /// Construct a Bins from an energy and a set of command-line flags
    pub fn from_params<S: crate::system::System>(sys: &S, params: BinningParams) -> Self {
        match params {
            BinningParams::Histogram { bin } => {
                Bins::Histogram(histogram::Bins::new(sys.energy(), bin))
            }
            BinningParams::Linear { bin } => {
                Bins::Linear(linear::Bins::new(sys.energy(), bin))
            }
        }
    }
}

impl Binning for Bins {
// impl Bins {
    fn num_states(&self) -> usize {
        match self {
            Bins::Histogram(b) => b.num_states(),
            Bins::Linear(b) => b.num_states(),
        }
    }
    fn count_states<F: Fn(Energy, PerEnergy) -> bool>(&self, f: F) -> usize {
        match self {
            Bins::Histogram(b) => b.count_states(f),
            Bins::Linear(b) => b.count_states(f),
        }
    }

    fn increment_count(&mut self, e: Energy, gamma: f64) {
        match self {
            Bins::Histogram(b) => b.increment_count(e, gamma),
            Bins::Linear(b) => b.increment_count(e, gamma),
        }
    }

    fn get_lnw(&self, e: Energy) -> f64 {
        match self {
            Bins::Histogram(b) => b.get_lnw(e),
            Bins::Linear(b) => b.get_lnw(e),
        }
    }
    fn get_count(&self, e: Energy) -> PerEnergy {
        match self {
            Bins::Histogram(b) => b.get_count(e),
            Bins::Linear(b) => b.get_count(e),
        }
    }

    fn set_lnw<F: Fn(Energy, PerEnergy) -> Option<f64>>(&mut self, f: F) {
        match self {
            Bins::Histogram(b) => b.set_lnw(f),
            Bins::Linear(b) => b.set_lnw(f),
        }
    }
    fn max_lnw(&self) -> f64 {
        match self {
            Bins::Histogram(b) => b.max_lnw(),
            Bins::Linear(b) => b.max_lnw(),
        }
    }
    fn min_lnw(&self) -> f64 {
        match self {
            Bins::Histogram(b) => b.min_lnw(),
            Bins::Linear(b) => b.min_lnw(),
        }
    }

    fn total_count(&self) -> u64 {
        match self {
            Bins::Histogram(b) => b.total_count(),
            Bins::Linear(b) => b.total_count(),
        }
    }
    fn max_count(&self) -> PerEnergy {
        match self {
            Bins::Histogram(b) => b.max_count(),
            Bins::Linear(b) => b.max_count(),
        }
    }
    fn min_count(&self) -> PerEnergy {
        match self {
            Bins::Histogram(b) => b.min_count(),
            Bins::Linear(b) => b.min_count(),
        }
    }

    fn accumulate_extra(&mut self, name: Interned, e: Energy, value: f64) {
        match self {
            Bins::Histogram(b) => b.accumulate_extra(name, e, value),
            Bins::Linear(b) => b.accumulate_extra(name, e, value),
        }
    }
    fn zero_out_extra(&mut self, name: Interned) {
        match self {
            Bins::Histogram(b) => b.zero_out_extra(name),
            Bins::Linear(b) => b.zero_out_extra(name),
        }
    }
    fn total_extra(&self, name: Interned, e: Energy) -> f64 {
        match self {
            Bins::Histogram(b) => b.total_extra(name, e),
            Bins::Linear(b) => b.total_extra(name, e),
        }
    }
    fn mean_count_extra(&self, extra: Interned) -> PerEnergy {
        match self {
            Bins::Histogram(b) => b.mean_count_extra(extra),
            Bins::Linear(b) => b.mean_count_extra(extra),
        }
    }

    fn mean_extra(&self, name: Interned, e: Energy) -> f64 {
        match self {
            Bins::Histogram(b) => b.mean_extra(name, e),
            Bins::Linear(b) => b.mean_extra(name, e),
        }
    }

    fn min_energy(&self) -> Energy {
        match self {
            Bins::Histogram(b) => b.min_energy(),
            Bins::Linear(b) => b.min_energy(),
        }
    }
    fn max_energy(&self) -> Energy {
        match self {
            Bins::Histogram(b) => b.max_energy(),
            Bins::Linear(b) => b.max_energy(),
        }
    }

    fn min_total_extra(&self, name: Interned) -> f64 {
        match self {
            Bins::Histogram(b) => b.min_total_extra(name),
            Bins::Linear(b) => b.min_total_extra(name),
        }
    }
    fn max_total_extra(&self, name: Interned) -> f64 {
        match self {
            Bins::Histogram(b) => b.max_total_extra(name),
            Bins::Linear(b) => b.max_total_extra(name),
        }
    }
    fn min_count_extra(&self, name: Interned) -> PerEnergy {
        match self {
            Bins::Histogram(b) => b.min_count_extra(name),
            Bins::Linear(b) => b.min_count_extra(name),
        }
    }
    fn min_count_extra_energy(&self, name: Interned) -> Energy {
        match self {
            Bins::Histogram(b) => b.min_count_extra_energy(name),
            Bins::Linear(b) => b.min_count_extra_energy(name),
        }
    }
    fn total_count_extra(&self, name: Interned) -> u64 {
        match self {
            Bins::Histogram(b) => b.total_count_extra(name),
            Bins::Linear(b) => b.total_count_extra(name),
        }
    }
    fn count_extra(&self, name: Interned, e: Energy) -> PerEnergy {
        match self {
            Bins::Histogram(b) => b.count_extra(name, e),
            Bins::Linear(b) => b.count_extra(name, e),
        }
    }
    fn new(e: Energy, de: Energy) -> Self {
        Bins::Histogram(histogram::Bins::new(e,de))
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
    fn accumulate_extra(&mut self, name: Interned, e: Energy, value: f64);
    /// Reset the accumulated extra data.
    fn zero_out_extra(&mut self, name: Interned);
    /// Find the total accumulated extra value that we accumulated at
    /// this energy.
    fn total_extra(&self, name: Interned, e: Energy) -> f64;
    /// Find the average extra value that we accumulated at this
    /// energy.
    fn mean_extra(&self, name: Interned, e: Energy) -> f64;
    /// Find the maximum extra value.
    fn max_total_extra(&self, name: Interned) -> f64;
    /// Find the minimum extra value.
    fn min_total_extra(&self, name: Interned) -> f64;
    /// Find the minimum extra count.
    fn min_count_extra(&self, name: Interned) -> PerEnergy;
    /// Find the energy with minimum count
    fn min_count_extra_energy(&self, name: Interned) -> Energy;
    /// Find the number of times we have counted this thing
    fn count_extra(&self, name: Interned, e: Energy) -> PerEnergy;
    /// Find the number of times we have counted this thing
    fn total_count_extra(&self, name: Interned) -> u64;
    /// Find the mean value of count_extra
    fn mean_count_extra(&self, name: Interned) -> PerEnergy;
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

/// Test that a binning works as expected
#[cfg(test)]
pub fn test_binning<B: Binning>() {
    use crate::system::units;
    println!("running a binning test!");
    let eps = units::EPSILON;
    let mut b = B::default();
    println!("made default");
    assert_eq!(b.get_count(eps), 0.0/eps);
    println!("tested get_count");
    b.increment_count(eps, 1.0);
    println!("incremented count!");
    assert!(b.get_count(eps) > 0.0/eps);
    println!("ran get_count!");
    let mydat = "datum".into();
    println!("interned");
    b.accumulate_extra(mydat, eps, 7.0);
    println!("accumulated extra!");
    assert_eq!(b.mean_extra(mydat, eps), 7.0);
    assert_eq!(b.total_extra(mydat, eps), 7.0);
    assert_eq!(b.max_total_extra(mydat), 7.0);
}
