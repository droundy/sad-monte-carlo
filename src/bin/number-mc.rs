extern crate sadmc;

use sadmc::system::optsquare::*;

use sadmc::mc::MonteCarlo;
use sadmc::mc::number::NumberMC;

fn main() {
    let mut mc = NumberMC::<SquareWell>::from_args::<SquareWellParams>();
    loop {
        mc.move_once();
    }
}
