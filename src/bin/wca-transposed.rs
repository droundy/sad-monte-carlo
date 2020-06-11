use sadmc::system::wca::{ Wca, WcaNParams };

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy_transposed::EnergyMC;

fn main() {
    let mut mc = EnergyMC::<Wca>::from_args::<WcaNParams>();
    loop {
        mc.move_once();
    }
}
