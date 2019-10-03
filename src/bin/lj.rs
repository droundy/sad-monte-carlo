extern crate sadmc;

use sadmc::system::lj::*;

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy::EnergyMC;

fn main() {
    let mut mc = EnergyMC::<Lj>::from_args::<LjParams>();
    loop {
        mc.move_once();
    }
}
