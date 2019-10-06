extern crate sadmc;

use sadmc::system::ising::{Ising, IsingParams};
use sadmc::system::System;

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy::EnergyMC;

fn main() {
    let mut mc = EnergyMC::<Ising, <Ising as System>::CollectedData>::from_args::<IsingParams>();
    loop {
        mc.move_once();
    }
}
