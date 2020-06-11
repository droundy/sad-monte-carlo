//! This crate is for sad monte carlo.

#![cfg_attr(feature = "strict", deny(warnings))]
#![deny(missing_docs)]

#[macro_use]
extern crate dimensioned;
#[macro_use]
extern crate serde_derive;

pub mod atomicfile;
pub mod mc;
pub mod prettyfloat;
pub mod rng;
pub mod system;
