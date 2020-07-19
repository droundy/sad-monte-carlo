use sadmc::system::lj::{Lj, LjParams};

use sadmc::mc::energy_transposed::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<Lj>::from_args::<LjParams>();
    loop {
        mc.move_once();
    }
}
