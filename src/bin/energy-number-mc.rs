use sadmc::system::optsquare::*;

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy_number::EnergyNumberMC;

fn main() {
    let mut mc = EnergyNumberMC::<SquareWell>::from_args::<SquareWellParams>();
    loop {
        mc.move_once();
    }
}
