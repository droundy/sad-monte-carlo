use sadmc::system::erfinv::{Parameters, ErfInv};

use sadmc::mc::grand_transposed::GrandMC;

fn main() {
    let mut mc = GrandMC::<ErfInv>::from_args::<Parameters>();
    loop {
        mc.run_once();
    }
}
