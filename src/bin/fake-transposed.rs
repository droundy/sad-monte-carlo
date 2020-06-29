use sadmc::system::fake::{FakeParams, Fake};

use sadmc::mc::energy_transposed::EnergyMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = EnergyMC::<Fake>::from_args::<FakeParams>();
    loop {
        mc.move_once();
    }
}
