use sadmc::system::any::{Any, AnyParams};

use sadmc::mc::energy::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<Any>::from_args::<AnyParams>();
    loop {
        mc.move_once();
    }
}
