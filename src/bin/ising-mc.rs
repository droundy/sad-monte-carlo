extern crate sadmc;

use sadmc::system::ising::{Ising, IsingParams};

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy::EnergyMC;

fn main() {
    let mut mc = EnergyMC::<Ising>::from_args::<IsingParams>();
    loop {
        mc.move_once();
    }
}
