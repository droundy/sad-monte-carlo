use sadmc::system::optsquare::*;

use sadmc::mc::energy_number::EnergyNumberMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyNumberMC::<SquareWell>::from_args::<SquareWellParams>();
    loop {
        mc.move_once();
    }
}
