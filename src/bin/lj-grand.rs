use sadmc::system::lj::*;

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy_number::EnergyNumberMC;

fn main() {
    let mut mc = EnergyNumberMC::<Lj>::from_args::<GrandLjParams>();
    loop {
        mc.move_once();
    }
}
