//! This crate is for sad monte carlo.

#![cfg_attr(feature = "strict", deny(warnings))]
#![cfg_attr(feature = "strict", deny(missing_docs))]

#[macro_use] extern crate dimensioned;
#[macro_use] extern crate serde_derive;
extern crate serde;
extern crate serde_yaml;

extern crate vector3d;

pub mod system;
