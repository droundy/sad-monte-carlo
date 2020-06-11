use sadmc::system::optsquare::*;
use sadmc::system::System;

use sadmc::mc::energy::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<SquareWell, <SquareWell as System>::CollectedData>::from_args::<
        SquareWellNParams,
    >();
    loop {
        mc.move_once();
    }
}
