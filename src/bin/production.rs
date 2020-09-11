use sadmc::system::any::{AnyParams, Any};

use sadmc::mc::energy_production::MC;

fn main() {
    let mut mc = MC::<Any>::from_args::<AnyParams>();
    loop {
        mc.run_once();
    }
}
