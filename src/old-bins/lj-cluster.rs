use sadmc::system::lj::*;

use sadmc::mc::energy::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<Lj>::from_args::<LjParams>();
    loop {
        mc.move_once();
    }
}
