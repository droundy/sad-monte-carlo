use sadmc::system::wca::*;
use sadmc::system::System;

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy::EnergyMC;

fn main() {
    let mut mc = EnergyMC::<Wca, <Wca as System>::CollectedData>::from_args::<WcaNParams>();
    loop {
        mc.move_once();
    }
}
