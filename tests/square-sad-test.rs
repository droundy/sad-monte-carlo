extern crate sadmc;

use std::default::Default;

use sadmc::system::{optsquare, square, System};

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy::{EnergyMC, EnergyMCParams};

#[test]
fn energies_agree() {
    let mut mcnew = EnergyMC::from_params(EnergyMCParams::default(),
                                          optsquare::SquareWell::from(optsquare::SquareWellNParams::default()),
                                          ::std::path::PathBuf::from("/dev/null"));
    let mut mcold = EnergyMC::from_params(EnergyMCParams::default(),
                                          square::SquareWell::from(square::SquareWellNParams::default()),
                                          ::std::path::PathBuf::from("/dev/null"));
    assert_eq!(mcnew.system.energy(), mcold.system.energy());
    for i in 0..1000000 {
        mcnew.move_once();
        mcold.move_once();
        println!("iteration {}", i);
        println!("new:");
        for r in mcnew.system.cell.positions.iter() {
            println!("  {}", r);
        }
        println!("old:");
        for r in mcold.system.positions.iter() {
            println!("  {} ({})", r, i);
            assert!(mcnew.system.cell.positions.contains(r));
        }
        assert_eq!(mcnew.system.energy(), mcold.system.energy());
    }
}
