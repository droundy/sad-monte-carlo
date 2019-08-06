extern crate sadmc;

use sadmc::system::optsquare::*;

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy_number::EnergyNumberMC;

fn main() {
    let mut mc = EnergyNumberMC::<Wca>::from_args::<WcaParams>();
    loop {
        mc.move_once();
    }
}
