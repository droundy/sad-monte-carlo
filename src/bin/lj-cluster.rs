use sadmc::system::lj::*;
use sadmc::system::System;

use sadmc::mc::energy::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<Lj, <Lj as System>::CollectedData>::from_args::<LjParams>();
    loop {
        mc.move_once();
    }
}
