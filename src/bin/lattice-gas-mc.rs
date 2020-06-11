use sadmc::system::lattice_gas::*;

use sadmc::mc::energy_number_no_translation::EnergyNumberMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyNumberMC::<LatticeGas>::from_args::<LatticeGasParams>();
    loop {
        mc.move_once();
    }
}
