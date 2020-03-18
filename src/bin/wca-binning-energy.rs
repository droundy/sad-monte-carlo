use sadmc::system::wca::{ Wca, WcaNParams };
use sadmc::system::System;

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy_binning::EnergyMC;

fn main() {
    let mut mc = EnergyMC::<Wca>::from_args::<WcaNParams>();
    loop {
        mc.move_once();
    }
}
