//! This crate is for sad monte carlo.

#![cfg_attr(feature = "strict", deny(warnings))]
#![deny(missing_docs)]

#[macro_use] extern crate dimensioned;
#[macro_use] extern crate serde_derive;
extern crate serde;
extern crate serde_yaml;
#[macro_use] extern crate clapme;

extern crate vector3d;
extern crate rand;
extern crate rand_core;
extern crate tempfile;

pub mod system;
pub mod mc;
pub mod rng;
pub mod atomicfile;
pub mod prettyfloat;
