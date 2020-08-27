use sadmc::system::wca::*;

use sadmc::mc::energy::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<Wca>::from_args::<WcaNParams>();
    loop {
        mc.move_once();
    }
}
