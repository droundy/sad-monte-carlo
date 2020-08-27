use sadmc::system::any::{AnyGrandParams, AnyGrand};

use sadmc::mc::grand_transposed::GrandMC;

fn main() {
    let mut mc = GrandMC::<AnyGrand>::from_args::<AnyGrandParams>();
    loop {
        mc.run_once();
    }
}
