extern crate sadmc;

use sadmc::system::optsquare::*;

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy::EnergyMC;

fn main() {
    let mut mc = EnergyMC::<Wca>::from_args::<WcaNParams>();
    loop {
        mc.move_once();
    }
}
