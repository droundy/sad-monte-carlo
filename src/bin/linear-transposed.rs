use sadmc::system::linear_energy::{Linear, LinearParams};

use sadmc::mc::energy_transposed::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<Linear>::from_args::<LinearParams>();
    loop {
        mc.move_once();
    }
}
