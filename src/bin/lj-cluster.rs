extern crate sadmc;

use sadmc::system::lj::*;
use sadmc::system::System;

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy::EnergyMC;

fn main() {
    let mut mc = EnergyMC::<Lj, <Lj as System>::CollectedData>::from_args::<LjParams>();
    loop {
        mc.move_once();
    }
}
