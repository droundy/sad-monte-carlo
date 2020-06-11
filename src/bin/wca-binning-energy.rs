use sadmc::system::wca::{Wca, WcaNParams};

use sadmc::mc::energy_binning::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<Wca>::from_args::<WcaNParams>();
    loop {
        mc.move_once();
    }
}
