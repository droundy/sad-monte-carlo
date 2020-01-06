use sadmc::system::optsquare::*;
use sadmc::system::System;

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy::EnergyMC;

fn main() {
    let mut mc = EnergyMC::<SquareWell, <SquareWell as System>::CollectedData>::from_args::<SquareWellNParams>();
    loop {
        mc.move_once();
    }
}
