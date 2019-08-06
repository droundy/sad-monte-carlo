extern crate sadmc;
extern crate vector3d;

use std::default::Default;

use vector3d::Vector3d;
use sadmc::system::units::SIGMA;
use sadmc::system::{optsquare, square, System};

use sadmc::mc::MonteCarlo;
use sadmc::mc::energy::{EnergyMC, EnergyMCParams};

#[test]
fn energies_agree() {
    let cell_width = 6.0*SIGMA;
    let box_diag = Vector3d::new(cell_width, cell_width, cell_width);
    let mut sparam = square::WcaNParams::default();
    let mut oparam = optsquare::WcaNParams::default();
    oparam._dim = optsquare::CellDimensionsGivenNumber::CellWidth(box_diag);
    sparam._dim = square::CellDimensionsGivenNumber::CellWidth(box_diag);
    let tempd = tempfile::TempDir::new().unwrap();
    let yaml_file = tempd.path().join("test.yaml");
    let mut mcnew = EnergyMC::from_params(EnergyMCParams::default(),
                                          optsquare::Wca::from(oparam), yaml_file.clone());
    let mut mcold = EnergyMC::from_params(EnergyMCParams::default(),
                                          square::Wca::from(sparam), yaml_file.clone());
    assert_eq!(mcnew.system.energy(), mcold.system.energy());
    for i in 0..10000 {
        mcnew.move_once();
        mcold.move_once();
        println!("iteration {}", i);
        println!("new vs old");
        for i in 0..mcnew.system.cell.positions.len() {
            assert_eq!(mcnew.system.cell.positions[i], mcold.system.positions[i]);
            // if mcnew.system.cell.positions[i] != mcold.system.positions[i] || true {
            //     println!("{:>width$} vs {:width$}", mcnew.system.cell.positions[i], mcold.system.positions[i],
            //              width=60);
            // }
        }
        for r in mcold.system.positions.iter() {
            assert!(mcnew.system.cell.positions.contains(r));
        }
        assert_eq!(mcnew.system.energy(), mcold.system.energy());
    }
}
