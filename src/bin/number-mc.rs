extern crate sadmc;

use sadmc::system::optsquare::*;

use sadmc::mc::MonteCarlo;
use sadmc::mc::number::NumberMC;

fn main() {
    let mut mc = NumberMC::<Wca>::from_args::<WcaParams>();
    loop {
        mc.move_once();
    }
}
