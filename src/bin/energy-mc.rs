use sadmc::system::optsquare::*;

use sadmc::mc::energy::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<SquareWell>::from_args::<
        SquareWellNParams,
    >();
    loop {
        mc.move_once();
    }
}
