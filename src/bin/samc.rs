extern crate sadmc;

use sadmc::system::square::*;

use sadmc::mc::MonteCarlo;
use sadmc::mc::samc::Samc;

fn main() {
    let mut mc = Samc::<SquareWell>::from_args::<SquareWellNParams>();
    loop {
        mc.move_once();
    }
}
