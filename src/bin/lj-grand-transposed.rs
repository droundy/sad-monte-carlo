use sadmc::system::lj::{GrandLjParams, Lj};

use sadmc::mc::grand_transposed::GrandMC;

fn main() {
    let mut mc = GrandMC::<Lj>::from_args::<GrandLjParams>();
    loop {
        mc.run_once();
    }
}
