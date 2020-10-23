use sadmc::system::any::{AnyParams, Any};

use sadmc::mc::energy_transposed::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<Any>::from_args::<AnyParams>();
    loop {
        mc.move_once();
    }
}
