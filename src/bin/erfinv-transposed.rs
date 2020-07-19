use sadmc::system::erfinv::{Parameters, ErfInv};

use sadmc::mc::energy_transposed::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<ErfInv>::from_args::<Parameters>();
    loop {
        mc.move_once();
    }
}
