use sadmc::system::wca::*;
use sadmc::system::System;

use sadmc::mc::energy::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<Wca, <Wca as System>::CollectedData>::from_args::<WcaNParams>();
    loop {
        mc.move_once();
    }
}
