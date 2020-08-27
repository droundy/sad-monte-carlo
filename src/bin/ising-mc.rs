use sadmc::system::ising::{Ising, IsingParams};

use sadmc::mc::energy::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<Ising>::from_args::<IsingParams>();
    loop {
        mc.move_once();
    }
}
