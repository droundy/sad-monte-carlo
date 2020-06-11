extern crate sadmc;

use sadmc::system::optsquare::*;

use sadmc::mc::number::NumberMC;
use sadmc::mc::MonteCarlo;

fn main() {
    let mut mc = NumberMC::<SquareWell>::from_args::<SquareWellParams>();
    loop {
        mc.move_once();
    }
}
