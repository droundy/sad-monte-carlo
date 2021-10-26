use sadmc::system::any::{Any, AnyParams};

use sadmc::mc::tempering::MC;

fn main() {
    let mut mc = MC::<Any>::from_args::<AnyParams>();
    loop {
        mc.run_once();
    }
}
