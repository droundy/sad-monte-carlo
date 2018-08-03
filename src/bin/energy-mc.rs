extern crate sadmc;

use sadmc::system::square::*;

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy::EnergyMC;

fn main() {
    let mut mc = EnergyMC::<SquareWell>::from_args::<SquareWellNParams>();
    loop {
        mc.move_once();
    }
}
