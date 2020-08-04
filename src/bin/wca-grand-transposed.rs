use sadmc::system::wca::{WcaParams, Wca};

use sadmc::mc::grand_transposed::GrandMC;

fn main() {
    let mut mc = GrandMC::<Wca>::from_args::<WcaParams>();
    loop {
        mc.run_once();
    }
}
