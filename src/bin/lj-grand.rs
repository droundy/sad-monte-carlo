use sadmc::system::lj::*;

use sadmc::mc::energy_number::EnergyNumberMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyNumberMC::<Lj>::from_args::<GrandLjParams>();
    loop {
        mc.move_once();
    }
}
