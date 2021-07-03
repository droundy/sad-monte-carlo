use sadmc::system::any::{Any, AnyParams};

use sadmc::mc::energy_replicas::MC;

fn main() {
    let mut mc = MC::<Any>::from_args::<AnyParams>();
    loop {
        mc.run_once();
    }
}
