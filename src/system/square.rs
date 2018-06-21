//! A square well fluid.

use super::*;

use dimensioned::Dimensionless;
use vector3d::Vector3d;

/// A square well fluid.
#[derive(Serialize, Deserialize, Debug)]
pub struct SquareWell {
    /// The dimensionless well width.
    well_width: Unitless,
    /// The atom positions
    positions: Vec<Vector3d<Length>>,
}

impl SquareWell {
    /// Create a new square well fluid.
    pub fn new(well_width: Unitless) -> SquareWell {
        SquareWell {
            well_width: well_width,
            positions: Vec::new(),
        }
    }
    fn max_interaction(&self) -> u64 {
        max_balls_within(self.well_width*units::SIGMA)
    }
}

impl System for SquareWell {
    fn energy(&self) -> Energy {
        unimplemented!()
    }
    fn delta_energy(&self) -> Option<Energy> {
        Some(units::EPSILON)
    }
    fn greatest_possible_energy(&self) -> Option<Energy> {
        Some(0.0*units::EPSILON)
    }
    fn lowest_possible_energy(&self) -> Option<Energy> {
        Some(-(self.positions.len() as f64)*(self.max_interaction() as f64)*units::EPSILON)
    }
}

fn max_balls_within(mut distance: Length) -> u64 {
    distance += 1e-10*units::SIGMA; // add a tiny, but necessary margin of error
    let a = 2.0_f64.sqrt()*units::SIGMA; // fcc lattice constant
    let c = (*(distance/a).value()).ceil() as i64 + 1; // number of cubic fcc cells to go out from center
    let mut num: i64 = -1; // number of balls within a given radius; don't count the center ball
    let d2 = distance*distance;

    for n in -c .. c+1 {
        for m in -c .. c+1 {
            for l in -c .. c+1 {
                let x0 = ((m+l) as f64)*a;
                let y0 = ((n+l) as f64)*a;
                let z0 = ((m+n) as f64)*a;
                if x0*x0 + y0*y0 + z0*z0 <= d2 {
                    num += 1;
                }
                if (x0 + 0.5*a)*(x0 + 0.5*a) + (y0 + 0.5*a)*(y0 + 0.5*a) + z0*z0 <= d2 {
                    num += 1;
                }
                if (x0 + 0.5*a)*(x0 + 0.5*a) + y0*y0 + (z0 + 0.5*a)*(z0 + 0.5*a) <= d2 {
                    num += 1;
                }
                if x0*x0 + (y0 + 0.5*a)*(y0 + 0.5*a) + (z0 + 0.5*a)*(z0 + 0.5*a) <= d2 {
                    num += 1;
                }
            }
        }
    }
    num as u64
}
