//! This crate is for sad monte carlo.

#![cfg_attr(feature = "strict", deny(warnings))]
#![deny(missing_docs)]

#[macro_use] extern crate dimensioned;
#[macro_use] extern crate serde_derive;
#[macro_use] extern crate clapme;

pub mod system;
pub mod mc;
pub mod rng;
pub mod atomicfile;
pub mod prettyfloat;
