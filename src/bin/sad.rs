extern crate sadmc;

use sadmc::system::square::*;

use sadmc::mc::MonteCarlo;
use sadmc::mc::sad::Sad;

fn main() {
    Sad::<SquareWell>::from_args();
    println!("hello world? hello");
}
