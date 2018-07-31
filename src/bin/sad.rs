extern crate sadmc;

use sadmc::system::square::*;

use sadmc::mc::MonteCarlo;
use sadmc::mc::sad::Sad;

fn main() {
    let mut sad = Sad::<SquareWell>::from_args();
    for _ in 0..1000 {
        sad.move_once();
    }
    println!("Rejected {}/{} = {}% of the moves",
             sad.rejected_moves, sad.moves,
             100.0*sad.rejected_moves as f64/sad.moves as f64);
    sad.checkpoint();
}
