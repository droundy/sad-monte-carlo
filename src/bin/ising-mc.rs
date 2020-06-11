use sadmc::system::ising::{Ising, IsingParams};
use sadmc::system::System;

use sadmc::mc::energy::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<Ising, <Ising as System>::CollectedData>::from_args::<IsingParams>();
    loop {
        mc.move_once();
    }
}
